/****************************************************************************
 *                                                                          *
 *   Lennard-Jones Simulator                                                *
 *   Copyright (C) 2020-2024 Ivo Filot <i.a.w.filot@tue.nl>                 *
 *                                                                          *
 *   This program is free software: you can redistribute it and/or modify   *
 *   it under the terms of the GNU Lesser General Public License as         *
 *   published by the Free Software Foundation, either version 3 of the     *
 *   License, or (at your option) any later version.                        *
 *                                                                          *
 *   This program is distributed in the hope that it will be useful,        *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *   GNU General Public License for more details.                           *
 *                                                                          *
 *   You should have received a copy of the GNU General Public license      *
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>. *
 *                                                                          *
 ****************************************************************************/

#include "lennardjonessimulation.h"

#include <array>
#include <cmath>
#include <cstring>

/**
 * @brief Construct a Lennard-Jones simulation.
 * @param _params Shared pointer to the parameter set.
 */
LennardJonesSimulation::LennardJonesSimulation(const std::shared_ptr<LennardJonesParameters>& _params) :
    params(_params)
{
    this->ttime = std::chrono::system_clock::now();

    double l = this->params->get_param<double>("cell_length");
    this->dims = glm::dvec3(l, l, l);

    const size_t n = (size_t)this->params->get_param<int>("nr_particles");

    // Allocate SoA buffers
    for (int b = 0; b < 2; ++b) {
        pos_buf[b].x.resize(n);
        pos_buf[b].y.resize(n);
        pos_buf[b].z.resize(n);
    }

    vx.resize(n); vy.resize(n); vz.resize(n);
    fx.resize(n); fy.resize(n); fz.resize(n);

    dijx.resize(n); dijy.resize(n); dijz.resize(n);

    // Initialize into write buffer (1), then publish it as render buffer (0/1 doesn’t matter yet)
    render_idx.store(0, std::memory_order_release);
    sim_write_idx = 0; // write initial into 0
    initialize();

    // neighbor list + forces use "current positions"
    build_neighbor_list(sim_write_idx);
    calculate_forces(sim_write_idx);

    publish_state(sim_write_idx); // publish initialized buffer
}

/**
 * @brief Perform one velocity-Verlet integration step.
 * @param step Current simulation step index.
 * @param dt   Time step size.
 */
void LennardJonesSimulation::integrate(unsigned int step, double dt) {
    const int read_idx  = render_idx.load(std::memory_order_acquire);
    const int write_idx = 1 - read_idx;

    // -------------------------
    // OpenMP: keep multiple loops inside one parallel region
    // -------------------------
    update_velocities_half_dt(dt);

    // thermostat uses ekin; scale velocities in parallel inside one region
    apply_berendsen_thermostat(dt);

    update_positions(dt, read_idx, write_idx);

    // The following parts are not OpenMP-parallel in your current design (due to data races),
    // but are still the main compute hot spot; we keep them as-is but SoA + manual math helps.
    calculate_forces(write_idx);
    apply_boundary_conditions(write_idx);

    if (check_update_neighbor_list()) {
        build_neighbor_list(write_idx);
    }

    update_velocities_half_dt(dt);

    etot.store(ekin.load(std::memory_order_relaxed) + epot.load(std::memory_order_relaxed),
               std::memory_order_relaxed);

    publish_state(write_idx);

    if (step % 1000 == 0) {
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - ttime;

        QString logMessage = QString("Step %1  Ekin = %2  Epot = %3  Etot = %4  Time = %5 s")
                                .arg(step, 4, 10, QChar('0'))
                                .arg(this->ekin.load(std::memory_order_relaxed), 0, 'f', 6)
                                .arg(this->epot.load(std::memory_order_relaxed), 0, 'f', 6)
                                .arg(this->etot.load(std::memory_order_relaxed), 0, 'f', 6)
                                .arg(elapsed_seconds.count(), 0, 'f', 4);

        qDebug() << logMessage;
        this->ttime = std::chrono::system_clock::now();
    }
}

/**
 * @brief Write current simulation state to a movie file.
 * @param moviefile Path to output file.
 * @param create    If true, file is truncated; otherwise appended.
 */
void LennardJonesSimulation::write_to_movie_file(const std::string& moviefile, bool create) {
    std::ofstream outfile;
    uint32_t n = (uint32_t)this->params->get_param<int>("nr_particles");

    const int ridx = render_idx.load(std::memory_order_acquire);
    const auto& px = pos_buf[ridx].x;
    const auto& py = pos_buf[ridx].y;
    const auto& pz = pos_buf[ridx].z;

    if(create) {
        outfile.open(moviefile, std::ios_base::trunc | std::ios::binary);
        char bufferN[sizeof(uint32_t)];
        std::memcpy(bufferN, &n, sizeof(uint32_t));
        outfile.write(bufferN, sizeof(uint32_t));
    } else {
        outfile.open(moviefile, std::ios_base::app | std::ios::binary);
    }

    // copy unit cell dimensions
    char buffer[sizeof(double) * 7];
    std::memcpy(&buffer[0], &this->dims[0], sizeof(double) * 3);
    outfile.write(buffer, sizeof(double) * 3);

    for(unsigned int i=0; i<n; i++) {
        double pos3[3] = { px[i], py[i], pz[i] };
        double vel3[3] = { vx[i], vy[i], vz[i] };
        double speed = std::sqrt(vel3[0]*vel3[0] + vel3[1]*vel3[1] + vel3[2]*vel3[2]);

        std::memcpy(&buffer[0], &pos3[0], sizeof(double) * 3);
        std::memcpy(&buffer[3 * sizeof(double)], &vel3[0], sizeof(double) * 3);
        std::memcpy(&buffer[6 * sizeof(double)], &speed, sizeof(double));

        outfile.write(buffer, sizeof(double) * 7);
    }

    outfile.close();
}

std::vector<glm::dvec3> LennardJonesSimulation::get_positions_copy() const {
    const int idx = render_idx.load(std::memory_order_acquire);
    const auto& px = pos_buf[idx].x;
    const auto& py = pos_buf[idx].y;
    const auto& pz = pos_buf[idx].z;

    std::vector<glm::dvec3> out(px.size());
    #pragma omp parallel for
    for (int i = 0; i < (int)px.size(); ++i) {
        out[i] = glm::dvec3(px[i], py[i], pz[i]);
    }
    return out;
}

std::vector<glm::dvec3> LennardJonesSimulation::get_velocities_copy() const {
    std::vector<glm::dvec3> out(vx.size());
    #pragma omp parallel for
    for (int i = 0; i < (int)vx.size(); ++i) {
        out[i] = glm::dvec3(vx[i], vy[i], vz[i]);
    }
    return out;
}

/**
 * @brief Initialize particle positions and velocities.
 */
void LennardJonesSimulation::initialize() {
    const size_t n = (size_t)this->params->get_param<int>("nr_particles");
    const double volume = this->dims.x * this->dims.y * this->dims.z;
    const double dl = std::pow(volume / (double)n, 1.0 / 3.0);

    size_t nx = (size_t) std::ceil(this->dims.x / dl);
    size_t ny = (size_t) std::ceil(this->dims.y / dl);
    size_t nz = (size_t) std::ceil(this->dims.z / dl);

    const double dx = dims.x/(double)nx;
    const double dy = dims.y/(double)ny;
    const double dz = dims.z/(double)nz;

    auto& px = pos_buf[sim_write_idx].x;
    auto& py = pos_buf[sim_write_idx].y;
    auto& pz = pos_buf[sim_write_idx].z;

    size_t count = 0;
    for(size_t i = 0; i<nx; i++) {
        const double x = (i+0.5)*dx;
        for(size_t j = 0; j<ny; j++) {
            const double y = (j+0.5)*dy;
            for(size_t k = 0; k<nz; k++) {
                const double z = (k+0.5)*dz;

                if (count >= n) break;

                px[count] = x;
                py[count] = y;
                pz[count] = z;
                count++;
            }
        }
    }

    // velocities init
    const double kT = this->params->get_param<double>("kT");
    const double m  = this->params->get_param<double>("mass");
    const double sqrtktm = std::sqrt(kT/m);

    double sumx = 0.0, sumy = 0.0, sumz = 0.0;
    for (size_t i=0; i<n; i++) {
        vx[i] = sqrtktm * get_gauss();
        vy[i] = sqrtktm * get_gauss();
        vz[i] = sqrtktm * get_gauss();
        sumx += vx[i]; sumy += vy[i]; sumz += vz[i];
    }

    const double invn = 1.0 / (double)n;
    sumx *= invn; sumy *= invn; sumz *= invn;

    for (size_t i=0; i<n; i++) {
        vx[i] -= sumx;
        vy[i] -= sumy;
        vz[i] -= sumz;
    }

    // reset dij
    std::fill(dijx.begin(), dijx.end(), 0.0);
    std::fill(dijy.begin(), dijy.end(), 0.0);
    std::fill(dijz.begin(), dijz.end(), 0.0);
}

double LennardJonesSimulation::get_gauss() const {
    static std::default_random_engine generator;
    static std::normal_distribution<double> distribution(0.0, 1.0);
    return distribution(generator);
}

/**
 * @brief Compute forces and potential energy.
 * @param cur_idx Index of position buffer to use.
 */
void LennardJonesSimulation::calculate_forces(int cur_idx) {
    const size_t n = (size_t)this->params->get_param<int>("nr_particles");

    const double rcut = this->params->get_param<double>("rcut");
    const double rcutsq = rcut * rcut;
    const double sigma = this->params->get_param<double>("sigma");
    const double sigmasq = sigma * sigma;

    const double epsilon = this->params->get_param<double>("epsilon");
    const double prefctr = 24.0 * epsilon;

    // cutoff shift
    double sr2 = sigmasq / rcutsq;
    double sr6 = sr2 * sr2 * sr2;
    double sr12 = sr6 * sr6;
    const double epot_cutoff = sr12 - sr6;

    // zero forces
    std::fill(fx.begin(), fx.end(), 0.0);
    std::fill(fy.begin(), fy.end(), 0.0);
    std::fill(fz.begin(), fz.end(), 0.0);

    const auto& px = pos_buf[cur_idx].x;
    const auto& py = pos_buf[cur_idx].y;
    const auto& pz = pos_buf[cur_idx].z;

    double epotsum = 0.0;

    for (int i = 0; i < (int)n; i++) {
        for (auto it = neighbor_list.begin() + neighbor_list[i]; *it != (unsigned)i; ++it) {
            const unsigned j = *it;

            // minimum image
            double dx = px[i] - px[j];
            double dy = py[i] - py[j];
            double dz = pz[i] - pz[j];

            if (dx >= dims.x * 0.5) dx -= dims.x;
            else if (dx < -dims.x * 0.5) dx += dims.x;

            if (dy >= dims.y * 0.5) dy -= dims.y;
            else if (dy < -dims.y * 0.5) dy += dims.y;

            if (dz >= dims.z * 0.5) dz -= dims.z;
            else if (dz < -dims.z * 0.5) dz += dims.z;

            const double d2 = dx*dx + dy*dy + dz*dz;
            if (d2 >= rcutsq) continue;

            const double invd2 = 1.0 / d2;
            const double s2 = sigmasq * invd2;
            const double s6 = s2 * s2 * s2;
            const double s12 = s6 * s6;

            epotsum += (s12 - s6 - epot_cutoff);

            const double fr = (2.0 * s12 - s6) * invd2;
            fx[i] += fr * dx;
            fy[i] += fr * dy;
            fz[i] += fr * dz;
        }
    }

    for (size_t i = 0; i < n; i++) {
        fx[i] *= prefctr;
        fy[i] *= prefctr;
        fz[i] *= prefctr;
    }

    double ep = epotsum * 4.0 * epsilon;
    ep *= 0.5; // double counting
    epot.store(ep, std::memory_order_relaxed);
}

/**
 * @brief Update velocities by half a timestep.
 * @param dt Time step size.
 */
void LennardJonesSimulation::update_velocities_half_dt(double dt) {
    const size_t n = (size_t)this->params->get_param<int>("nr_particles");
    const double m = this->params->get_param<double>("mass");
    const double factor = 0.5 * dt / m;

    double ekinsum = 0.0;

    // single parallel region to reduce overhead
    #pragma omp parallel
    {
        double local_sum = 0.0;

        #pragma omp for schedule(static)
        for (int i = 0; i < (int)n; i++) {
            vx[i] += factor * fx[i];
            vy[i] += factor * fy[i];
            vz[i] += factor * fz[i];

            local_sum += (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
        }

        #pragma omp atomic
        ekinsum += local_sum;
    }

    ekin.store(0.5 * m * ekinsum, std::memory_order_relaxed);
}

/**
 * @brief Apply Berendsen thermostat scaling.
 * @param dt Time step size.
 */
void LennardJonesSimulation::apply_berendsen_thermostat(double dt) {
    const size_t n = (size_t)this->params->get_param<int>("nr_particles");
    const double tau = this->params->get_param<double>("tau");
    const double kT = this->params->get_param<double>("kT");

    const double ek = ekin.load(std::memory_order_relaxed);
    const double ekin0 = 1.5 * (double)(n - 1) * kT;
    const double lambda = std::sqrt(1.0 + dt / tau * (ekin0 / ek - 1.0));

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < (int)n; i++) {
        vx[i] *= lambda;
        vy[i] *= lambda;
        vz[i] *= lambda;
    }
}

/**
 * @brief Update positions: read from read_idx, write to write_idx (true double buffer)
 */
void LennardJonesSimulation::update_positions(double dt, int read_idx, int write_idx) {
    const size_t n = (size_t)this->params->get_param<int>("nr_particles");

    const auto& prx = pos_buf[read_idx].x;
    const auto& pry = pos_buf[read_idx].y;
    const auto& prz = pos_buf[read_idx].z;

    auto& pwx = pos_buf[write_idx].x;
    auto& pwy = pos_buf[write_idx].y;
    auto& pwz = pos_buf[write_idx].z;

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < (int)n; i++) {
        const double incx = vx[i] * dt;
        const double incy = vy[i] * dt;
        const double incz = vz[i] * dt;

        pwx[i] = prx[i] + incx;
        pwy[i] = pry[i] + incy;
        pwz[i] = prz[i] + incz;

        dijx[i] += incx;
        dijy[i] += incy;
        dijz[i] += incz;
    }
}

/**
 * @brief Apply periodic boundary conditions.
 * @param cur_idx Index of position buffer.
 */
void LennardJonesSimulation::apply_boundary_conditions(int cur_idx) {
    const size_t n = (size_t)this->params->get_param<int>("nr_particles");

    auto& px = pos_buf[cur_idx].x;
    auto& py = pos_buf[cur_idx].y;
    auto& pz = pos_buf[cur_idx].z;

    const double invLx = 1.0 / dims.x;
    const double invLy = 1.0 / dims.y;
    const double invLz = 1.0 / dims.z;

    for (size_t i = 0; i < n; i++) {
        px[i] -= dims.x * std::floor(px[i] * invLx);
        py[i] -= dims.y * std::floor(py[i] * invLy);
        pz[i] -= dims.z * std::floor(pz[i] * invLz);
    }
}

/**
 * @brief Build Verlet neighbor list.
 * @param cur_idx Index of position buffer.
 */
void LennardJonesSimulation::build_neighbor_list(int cur_idx) {
    const double cutoff = this->params->get_param<double>("rcut");
    const double shell  = this->params->get_param<double>("shell");
    const double csc    = cutoff + shell;
    const double cutsq  = csc * csc;

    unsigned int nx = (unsigned int)(this->dims.x / csc);
    unsigned int ny = (unsigned int)(this->dims.y / csc);
    unsigned int nz = (unsigned int)(this->dims.z / csc);

    nx = std::max(1u, nx);
    ny = std::max(1u, ny);
    nz = std::max(1u, nz);

    const double dx = dims.x / (double)nx;
    const double dy = dims.y / (double)ny;
    const double dz = dims.z / (double)nz;

    const size_t npos = pos_buf[cur_idx].x.size();

    const auto& px = pos_buf[cur_idx].x;
    const auto& py = pos_buf[cur_idx].y;
    const auto& pz = pos_buf[cur_idx].z;

    std::vector<std::vector<unsigned int>> subcells(nx * ny * nz);

    for (unsigned int i = 0; i < (unsigned)npos; i++) {
        unsigned int ix = (unsigned int)(px[i] / dx);
        unsigned int iy = (unsigned int)(py[i] / dy);
        unsigned int iz = (unsigned int)(pz[i] / dz);

        ix = std::min(ix, nx - 1);
        iy = std::min(iy, ny - 1);
        iz = std::min(iz, nz - 1);

        unsigned int cellid = iz * nx * ny + iy * nx + ix;
        subcells[cellid].push_back(i);
    }

    std::vector<std::vector<unsigned int>> nnlist(npos);

    #pragma omp parallel for
    for (int i = 0; i < (int)npos; i++) {
        dijx[i] = 0.0;
        dijy[i] = 0.0;
        dijz[i] = 0.0;

        unsigned int ix = (unsigned int)(px[i] / dx);
        unsigned int iy = (unsigned int)(py[i] / dy);
        unsigned int iz = (unsigned int)(pz[i] / dz);

        ix = std::min(ix, nx - 1);
        iy = std::min(iy, ny - 1);
        iz = std::min(iz, nz - 1);

        unsigned int lx = ix == 0 ? (nx - 1) : (ix - 1);
        unsigned int ly = iy == 0 ? (ny - 1) : (iy - 1);
        unsigned int lz = iz == 0 ? (nz - 1) : (iz - 1);

        unsigned int hx = ix == (nx - 1) ? 0 : (ix + 1);
        unsigned int hy = iy == (ny - 1) ? 0 : (iy + 1);
        unsigned int hz = iz == (nz - 1) ? 0 : (iz + 1);

        std::array<unsigned int, 3> xx = {lx, ix, hx};
        std::array<unsigned int, 3> yy = {ly, iy, hy};
        std::array<unsigned int, 3> zz = {lz, iz, hz};

        std::array<unsigned int, 27> cellids;
        unsigned int count = 0;
        for (unsigned int cx : xx) {
            for (unsigned int cy : yy) {
                for (unsigned int cz : zz) {
                    cellids[count++] = cz * nx * ny + cy * nx + cx;
                }
            }
        }

        for (unsigned int cid : cellids) {
            for (unsigned int pid : subcells[cid]) {
                if (pid == (unsigned)i) continue;

                double dxp = px[i] - px[pid];
                double dyp = py[i] - py[pid];
                double dzp = pz[i] - pz[pid];

                if (dxp >= dims.x * 0.5) dxp -= dims.x;
                else if (dxp < -dims.x * 0.5) dxp += dims.x;

                if (dyp >= dims.y * 0.5) dyp -= dims.y;
                else if (dyp < -dims.y * 0.5) dyp += dims.y;

                if (dzp >= dims.z * 0.5) dzp -= dims.z;
                else if (dzp < -dims.z * 0.5) dzp += dims.z;

                const double dist2 = dxp*dxp + dyp*dyp + dzp*dzp;
                if (dist2 <= cutsq) {
                    nnlist[i].push_back(pid);
                }
            }
        }
    }

    neighbor_list.clear();
    neighbor_list.resize(npos);

    for (unsigned int i = 0; i < (unsigned)npos; i++) {
        neighbor_list[i] = (unsigned int)neighbor_list.size();
        neighbor_list.insert(neighbor_list.end(), nnlist[i].begin(), nnlist[i].end());
        neighbor_list.push_back(i); // sentinel
    }
}

/**
 * @brief Check whether neighbor list rebuild is required.
 * @return True if rebuild is needed.
 */
bool LennardJonesSimulation::check_update_neighbor_list() const {
    const double shell = this->params->get_param<double>("shell");
    const double shellsq = shell * shell;
    const double thresh = shellsq * 0.25;

    const size_t n = dijx.size();
    for (size_t i = 0; i < n; i++) {
        const double d2 = dijx[i]*dijx[i] + dijy[i]*dijy[i] + dijz[i]*dijz[i];
        if (d2 > thresh) return true;
    }
    return false;
}

/**
 * @brief Publish simulation buffer as render buffer.
 * @param new_render_idx Buffer index to publish.
 */
void LennardJonesSimulation::publish_state(int new_render_idx) {
    render_idx.store(new_render_idx, std::memory_order_release);
}