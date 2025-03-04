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

/**
 * @brief      default constructor
 *
 * @param[in]  _params  set of parameters
 */
LennardJonesSimulation::LennardJonesSimulation(const std::shared_ptr<LennardJonesParameters>& _params) :
    params(_params)
{
    this->ttime = std::chrono::system_clock::now();

    double l = this->params->get_param<double>("cell_length");
    this->dims = glm::dvec3(l,l,l);

    this->initialize();
    this->build_neighbor_list();
    this->calculate_forces();
}

/**
 * @brief      Perform integration by timestep dt
 *
 * @param[in]  step  The step number
 * @param[in]  dt    timestep
 */
void LennardJonesSimulation::integrate(unsigned int step, double dt) {
    this->update_velocities_half_dt(dt);
    this->apply_berendsen_thermostat(dt);
    this->update_positions(dt);
    this->calculate_forces();
    this->apply_boundary_conditions();

    // check if neighbor list needs to be rebuild
    if(this->check_update_neighbor_list()) {
        this->build_neighbor_list();
    }

    this->update_velocities_half_dt(dt);
    this->etot = this->ekin + this->epot;

    if (step % 1000 == 0) {
        // Calculate elapsed time
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - ttime;
    
        // Format string using QString::arg()
        QString logMessage = QString("Step %1  Ekin = %2  Epot = %3  Etot = %4  Time = %5 s")
                                .arg(step, 4, 10, QChar('0'))  // Step, padded to 4 digits
                                .arg(this->ekin, 0, 'f', 6)   // Kinetic energy, 6 decimal places
                                .arg(this->epot, 0, 'f', 6)   // Potential energy, 6 decimal places
                                .arg(this->etot, 0, 'f', 6)   // Total energy, 6 decimal places
                                .arg(elapsed_seconds.count(), 0, 'f', 4); // Time, 4 decimal places
    
        qDebug() << logMessage;
    
        // Reset timer
        this->ttime = std::chrono::system_clock::now();
    }
}

/**
 * @brief      Write to movie file for current state
 *
 * @param[in]  moviefile  Path to the movie file
 * @param[in]  create     Whether to truncate (true) or to append (false)
 */
void LennardJonesSimulation::write_to_movie_file(const std::string& moviefile, bool create) {
    std::ofstream outfile;
    uint32_t nr_particles = this->params->get_param<int>("nr_particles");

    if(create) {
        outfile.open(moviefile, std::ios_base::trunc | std::ios::binary);
        char buffer[sizeof(uint32_t)];
        std::memcpy(buffer, &nr_particles, sizeof(uint32_t));
        outfile.write(buffer, sizeof(uint32_t));
    } else {
        outfile.open(moviefile, std::ios_base::app | std::ios::binary);
    }

    // copy unit cell dimensions
    char buffer[sizeof(double) * 7];
    std::memcpy(&buffer[0], &this->dims[0], sizeof(double) * 3);
    outfile.write(buffer, sizeof(double) * 3);

    for(unsigned int i=0; i<nr_particles; i++) {
        std::memcpy(&buffer[0], &this->positions[i][0], sizeof(double) * 3);
        std::memcpy(&buffer[3 * sizeof(double)], &this->velocities[i][0], sizeof(double) * 3);
        double velocity = glm::length(this->velocities[i]);
        std::memcpy(&buffer[6 * sizeof(double)], &velocity, sizeof(double));

        outfile.write(buffer, sizeof(double) * 7);
    }

    outfile.close();
}

/**
 * @brief      Initialize unit cell
 */
void LennardJonesSimulation::initialize() {
    const size_t nr_particles = this->params->get_param<int>("nr_particles");
    const double volume = this->dims.x * this->dims.y * this->dims.z;
    const double dl = std::pow(volume / (double)nr_particles, 1.0 / 3.0);

    // calculate dimensions of rectangular lattice to position particles on
    size_t nx = (size_t) std::ceil(this->dims.x / dl);
    size_t ny = (size_t) std::ceil(this->dims.y / dl);
    size_t nz = (size_t) std::ceil(this->dims.z / dl);
    double dx = dims.x/(double)nx;
    double dy = dims.y/(double)ny;
    double dz = dims.z/(double)nz;

    // expand vectors to match number of particles
    this->positions.resize(nr_particles);
    this->velocities.resize(nr_particles);
    this->forces.resize(nr_particles);

    // position all particles on a rectangular lattice
    size_t count = 0;
    for(size_t i = 0; i<nx; i++) {
        const double x = (i+0.5)*dx;
        for(size_t j = 0; j<ny; j++) {
            const double y = (j+0.5)*dy;
            for(size_t k = 0; k<nz; k++) {
                const double z = (k+0.5)*dz;

                if (count >= nr_particles) {
                    break;
                }

                this->positions[count] = glm::dvec3(x,y,z);
                count++;
            }
        }
    }

    // obtain data from params object
    double kT = this->params->get_param<double>("kT");
    double m = this->params->get_param<double>("mass");

    double sqrtktm = std::sqrt(kT/m);

    glm::dvec3 sum;

    // randomly set velocities
    for (size_t i=0; i<nr_particles; i++) {
        this->velocities[i] = sqrtktm * glm::dvec3(this->get_gauss(), this->get_gauss(), this->get_gauss());
        sum += this->velocities[i];
    }

    sum /= (double)nr_particles;

    // subtract velocity average to remove nett momentum
    for (size_t i=0; i<nr_particles; i++) {
        this->velocities[i] -= sum;
    }

    this->dij_list.resize(nr_particles);
}

/**
 * @brief      Get number from normal distribution between -1 and 1
 *
 * @return     Random number
 */
double LennardJonesSimulation::get_gauss() const {
    // construct RNG
    static std::default_random_engine generator;
    static std::normal_distribution<double> distribution(0.0, 1.0); // avg 0, sigma 1

    return distribution(generator);
}

/**
 * @brief      Calculates the forces.
 */
void LennardJonesSimulation::calculate_forces() {
    const size_t nr_particles = this->params->get_param<int>("nr_particles");
    const double rcut = this->params->get_param<double>("rcut");
    const double rcutsq = rcut * rcut;
    const double sigma = this->params->get_param<double>("sigma");
    const double sigmasq = sigma * sigma;

    // set potential energy to zero
    this->epot = 0.0;

    double sr2 = sigmasq / rcutsq;
    double sr6 = sr2 * sr2 * sr2;
    double sr12 = sr6 * sr6;
    const double epot_cutoff = sr12 - sr6;

    // set forces to zero
    memset(&this->forces[0], 0, nr_particles * sizeof(double) * 3);

    double epotsum = 0.0;

    // TODO: this function could potentially benefit from multi-hreading
    // though this function suffers from a race condition
    for(int i=0; i<nr_particles; i++) {
        for(auto it = this->neighbor_list.begin() + this->neighbor_list[i]; *it != i; ++it) {
            glm::dvec3 rij = this->positions[i] - this->positions[*it];

            for(unsigned int k=0; k<3; k++) {
                if(rij[k] >= this->dims[k] * 0.5) {
                    rij[k] -= this->dims[k];
                    continue;
                }

                if(rij[k] < -this->dims[k] * 0.5) {
                    rij[k] += this->dims[k];
                }
            }

            double distsq = glm::length2(rij);

            if(distsq >= rcutsq) {
                continue;
            }

            sr2 = sigmasq / distsq;
            sr6 = sr2 * sr2 * sr2;
            sr12 = sr6 * sr6;

            epotsum += (sr12 - sr6 - epot_cutoff);
            double fr = (2.0 * sr12 - sr6) / distsq;
            this->forces[i] += fr * rij;
        }
    }

    const double epsilon = this->params->get_param<double>("epsilon");
    const double prefctr = 24.0 * epsilon;

    for(unsigned int i=0; i<nr_particles; i++) {
        this->forces[i] *= prefctr;
    }

    this->epot = epotsum * 4.0 * epsilon;

    // multiply by 0.5 because we have double counting
    this->epot *= 0.5;
}

/**
 * @brief      Update velocities by half a timestep
 *
 * @param[in]  dt    Timestep (is divided by half INSIDE the function)
 */
void LennardJonesSimulation::update_velocities_half_dt(double dt) {
    const double m = this->params->get_param<double>("mass");
    const size_t nr_particles = this->params->get_param<int>("nr_particles");

    this->ekin = 0.0;
    const double factor = 0.5 / m * dt;

    double ekinsum = 0.0;

    #pragma omp parallel for reduction(+: ekinsum)
    for(int i=0; i<nr_particles; i++) {
        this->velocities[i] += factor * this->forces[i];
        ekinsum += glm::length2(this->velocities[i]);
    }

    this->ekin = ekinsum * 0.5 * m;
}

/**
 * @brief      Update positions by timestep dt
 *
 * @param[in]  dt    Timestep
 */
void LennardJonesSimulation::update_positions(double dt) {
    const size_t nr_particles = this->params->get_param<int>("nr_particles");

    #pragma omp parallel for schedule(static)
    for(size_t i=0; i<nr_particles; i++) {
        const glm::dvec3 inc = this->velocities[i] * dt;
        this->positions[i] += inc;
        this->dij_list[i] += inc;
    }
}

/**
 * @brief      Applies periodic boundary conditions to the positions
 */
void LennardJonesSimulation::apply_boundary_conditions() {
    const size_t nr_particles = this->params->get_param<int>("nr_particles");

    const glm::dvec3 invL = 1.0 / this->dims;

    for(unsigned int i=0; i<nr_particles; i++) {
        this->positions[i].x -= this->dims.x * std::floor(this->positions[i].x * invL.x);
        this->positions[i].y -= this->dims.y * std::floor(this->positions[i].y * invL.y);
        this->positions[i].z -= this->dims.z * std::floor(this->positions[i].z * invL.z);
    }
}

/**
 * @brief      Perform velocity scaling using Berendsen thermostat
 *
 * @param[in]  dt    Timestep
 */
void LennardJonesSimulation::apply_berendsen_thermostat(double dt) {
    const size_t nr_particles = this->params->get_param<int>("nr_particles");
    const double tau = this->params->get_param<double>("tau");
    const double kT = this->params->get_param<double>("kT");

    double ekin0 = 1.5*((double)(nr_particles-1))*kT;
    const double lambda = std::sqrt(1.0 + dt / tau * (ekin0 / this->ekin - 1.0));

    #pragma omp parallel for schedule(static)
    for(int i=0; i<nr_particles; i++) {
        this->velocities[i] *= lambda;
    }
}

/**
 * @brief      Build neighbor list
 */
void LennardJonesSimulation::build_neighbor_list() {
    const double cutoff = this->params->get_param<double>("rcut");
    const double shell = this->params->get_param<double>("shell");
    const double csc = cutoff + shell;
    const double cutsq = csc * csc;

    unsigned int nx = this->dims.x / csc;
    unsigned int ny = this->dims.y / csc;
    unsigned int nz = this->dims.z / csc;

    // specify subcell sizes
    double dx = dims.x / (double)nx;
    double dy = dims.y / (double)ny;
    double dz = dims.z / (double)nz;

    // place all particles in subcells
    std::vector<std::vector<unsigned int> > subcells(nx * ny * nz);

    for(unsigned int i=0; i<this->positions.size(); i++) {
        unsigned int ix = positions[i].x / dx;
        unsigned int iy = positions[i].y / dy;
        unsigned int iz = positions[i].z / dz;

        unsigned int cellid = iz * nx * ny + iy * nx + ix;
        subcells[cellid].push_back(i);
    }

    // collect neighbors for each particle
    std::vector<std::vector<unsigned int>> nnlist(this->positions.size());
    #pragma omp parallel for
    for(int i=0; i<this->positions.size(); i++) {
        this->dij_list[i] = glm::dvec3(0.0, 0.0, 0.0);

        unsigned int ix = positions[i].x / dx;
        unsigned int iy = positions[i].y / dy;
        unsigned int iz = positions[i].z / dz;

        // build list of neighboring cells
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
        for(unsigned int cx : xx) {
            for(unsigned int cy : yy) {
                for(unsigned int cz : zz) {
                    cellids[count] = cz * nx * ny + cy * nx + cx;
                    count++;
                }
            }
        }

        for(unsigned int cid : cellids) {
            for(unsigned int pid : subcells[cid]) {
                if(pid == i) {
                    continue;
                }

                glm::dvec3 rij = this->positions[i] - this->positions[pid];

                for(unsigned int k=0; k<3; k++) {
                    if(rij[k] >= this->dims[k] * 0.5) {
                        rij[k] -= this->dims[k];
                        continue;
                    }

                    if(rij[k] < -this->dims[k] * 0.5) {
                        rij[k] += this->dims[k];
                    }
                }

                const double dist2 = glm::length2(rij);

                if(dist2 <= cutsq) {
                    nnlist[i].push_back(pid);
                }
            }
        }
    }

    // build continuous list
    this->neighbor_list.resize(this->positions.size());
    for(unsigned int i=0; i<this->positions.size(); i++) {
        this->neighbor_list[i] = this->neighbor_list.size();
        this->neighbor_list.insert(this->neighbor_list.end(), nnlist[i].begin(), nnlist[i].end());
        this->neighbor_list.push_back(i);
    }
}

/**
 * @brief      Check whether neighbor list needs to be updated
 *
 * @return     Whether neighbor list needs to be updated
 */
bool LennardJonesSimulation::check_update_neighbor_list() const {
    const double shell = this->params->get_param<double>("shell");
    const double shellsq = shell * shell;

    for(const auto& dij : this->dij_list) {
        if(glm::length2(dij) > shellsq * 0.25) {
            return true;
        }
    }

    return false;
}
