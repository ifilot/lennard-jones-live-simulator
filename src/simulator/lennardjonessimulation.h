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

/****************************************************************************
 *                                                                          *
 *   Lennard-Jones Simulator                                                *
 *                                                                          *
 ****************************************************************************/

#ifndef LENNARDJONESSIMULATION_H
#define LENNARDJONESSIMULATION_H

#include <vector>
#include <memory>
#include <random>
#include <iostream>
#include <fstream>
#include <chrono>
#include <atomic>
#include <algorithm>

#include <QDebug>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>

#include "lennardjonesparameters.h"

/**
 * @class LennardJonesSimulation
 * @brief Molecular dynamics simulator for particles interacting via the Lennard-Jones potential.
 *
 * This class implements a velocity-Verlet integrator combined with:
 * - Periodic boundary conditions
 * - Verlet neighbor lists
 * - Berendsen thermostat
 *
 * ## Performance-oriented design
 * - Structure-of-Arrays (SoA) layout for particle state (x[], y[], z[]) for cache efficiency
 * - Double-buffered positions to decouple simulation and rendering threads
 * - Lock-free publishing of simulation state via atomic buffer index
 *
 * ## Thread model
 * - Simulation thread writes to a "write buffer"
 * - Rendering/UI threads read from a stable "render buffer"
 * - Buffer swap is atomic → no blocking or data races
 */
class LennardJonesSimulation {
private:

    /// Simulation box dimensions (periodic unit cell)
    glm::dvec3 dims;

    /**
     * @struct PosBuffer
     * @brief Structure-of-arrays container for particle positions.
     *
     * This layout improves memory locality and vectorization compared to AoS (glm::dvec3).
     */
    struct PosBuffer {
        std::vector<double> x, y, z;
    };

    /// Double buffer for particle positions (read/write separation)
    PosBuffer pos_buf[2];

    /// Index of buffer currently visible to renderer (atomic, lock-free)
    std::atomic<int> render_idx{0};

    /// Index of buffer simulation thread writes into
    int sim_write_idx = 1;

    // --- Particle dynamic state (SoA) ---
    std::vector<double> vx, vy, vz;   ///< Particle velocities
    std::vector<double> fx, fy, fz;   ///< Particle forces

    /// Accumulated displacements since last neighbor-list build
    std::vector<double> dijx, dijy, dijz;

    std::shared_ptr<LennardJonesParameters> params;
    std::vector<unsigned int> neighbor_list; ///< Verlet neighbor list (compressed format)

    // --- Energies (atomic for thread-safe readout) ---
    std::atomic<double> epot{0.0};
    std::atomic<double> ekin{0.0};
    std::atomic<double> etot{0.0};

    /// Timer used for performance logging
    std::chrono::time_point<std::chrono::system_clock> ttime;

public:

    /**
     * @brief Construct a Lennard-Jones simulation.
     * @param _params Shared pointer to the parameter set.
     */
    LennardJonesSimulation(const std::shared_ptr<LennardJonesParameters>& _params);

    /**
     * @brief Perform one integration step using the velocity-Verlet scheme.
     * @param step Simulation step index.
     * @param dt Time step size.
     */
    void integrate(unsigned int step, double dt);

    /// Access simulation parameters
    inline const auto& get_params() const { return this->params; }

    /// Access box dimensions
    inline const auto& get_dims() const { return this->dims; }

    // ---------- Rendering access (SoA) ----------

    /// Get current render buffer index (thread-safe)
    inline int get_render_index() const {
        return render_idx.load(std::memory_order_acquire);
    }

    /// Pointer accessors for direct GPU upload (no copies)
    inline const double* get_render_x() const { return pos_buf[get_render_index()].x.data(); }
    inline const double* get_render_y() const { return pos_buf[get_render_index()].y.data(); }
    inline const double* get_render_z() const { return pos_buf[get_render_index()].z.data(); }

    /// Number of particles
    inline size_t get_particle_count() const {
        return pos_buf[get_render_index()].x.size();
    }

    /// Energy getters
    inline double get_etot() const { return etot.load(std::memory_order_relaxed); }
    inline double get_ekin() const { return ekin.load(std::memory_order_relaxed); }
    inline double get_epot() const { return epot.load(std::memory_order_relaxed); }

    /// Compatibility functions (AoS conversion for GUI)
    std::vector<glm::dvec3> get_positions_copy() const;
    std::vector<glm::dvec3> get_velocities_copy() const;

    /// Write trajectory snapshot to movie file
    void write_to_movie_file(const std::string& moviefile, bool create = false);

private:

    /**
     * @brief Initialize particle positions and velocities.
     */
    void initialize();

    /**
     * @brief Generate a Gaussian-distributed random number.
     * @return Normally distributed random value.
     */
    double get_gauss() const;

     /**
     * @brief Compute forces and potential energy.
     * @param cur_idx Index of position buffer to use.
     */
    void calculate_forces(int cur_idx);

    /**
     * @brief Update velocities by half a timestep.
     * @param dt Time step size.
     */
    void update_velocities_half_dt(double dt);

    /**
     * @brief Update particle positions.
     * @param dt        Time step size.
     * @param read_idx  Index of buffer to read positions from.
     * @param write_idx Index of buffer to write updated positions to.
     */
    void update_positions(double dt, int read_idx, int write_idx);

    /**
     * @brief Apply periodic boundary conditions.
     * @param cur_idx Index of position buffer.
     */
    void apply_boundary_conditions(int cur_idx);

    /**
     * @brief Apply Berendsen thermostat scaling.
     * @param dt Time step size.
     */
    void apply_berendsen_thermostat(double dt);

    /**
     * @brief Build Verlet neighbor list.
     * @param cur_idx Index of position buffer.
     */
    void build_neighbor_list(int cur_idx);

    /**
     * @brief Check whether neighbor list rebuild is required.
     * @return True if rebuild is needed.
     */
    bool check_update_neighbor_list() const;

    /**
     * @brief Publish simulation buffer as render buffer.
     * @param new_render_idx Buffer index to publish.
     */
    void publish_state(int new_render_idx);
};

#endif // LENNARDJONESSIMULATION_H
