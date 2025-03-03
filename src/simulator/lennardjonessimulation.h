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

#ifndef LENNARDJONESSIMULATION_H
#define LENNARDJONESSIMULATION_H

#include <vector>
#include <memory>
#include <random>
#include <iostream>
#include <fstream>
#include <chrono>

#include <QDebug>

#define GLM_ENABLE_EXPERIMENTAL

#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>
#include <glm/gtx/string_cast.hpp>

#include "lennardjonesparameters.h"

/**
 * @brief      Unit cell containing all particles
 */
class LennardJonesSimulation {
private:
    glm::dvec3 dims;                                        // dimensions of the unit cell
    std::vector<glm::dvec3> positions;                      // positions of the particles in the unit cell
    std::vector<glm::dvec3> velocities;                     // velocities of the particles in the unit cell
    std::vector<glm::dvec3> forces;                         // forces on the particles in the unit cell
    std::shared_ptr<LennardJonesParameters> params;         // pointer to the parameters object
    std::vector<unsigned int> neighbor_list;                // neighbor list
    std::vector<glm::dvec3> dij_list;                       // progessive distance (used for updating neighbor list)

    double epot = 0.0;  // potential energy
    double ekin = 0.0;  // kinetic energy
    double etot = 0.0;  // total energy

    std::chrono::time_point<std::chrono::system_clock> ttime;

public:
    /**
     * @brief      default constructor
     *
     * @param[in]  _params  set of parameters
     */
    LennardJonesSimulation(const std::shared_ptr<LennardJonesParameters>& _params);

    /**
     * @brief      Perform integration by timestep dt
     *
     * @param[in]  step  The step number
     * @param[in]  dt    timestep
     */
    void integrate(unsigned int step, double dt);

    /**
     * @brief Get Lennard Jones parameters
     * @return Parameter object
     */
    inline const auto& get_params() const {
        return this->params;
    }

    /**
     * @brief      get the dimensions vector
     *
     * @return     dimensions vector
     */
    inline const auto& get_dims() const {
        return this->dims;
    }

    /**
     * @brief      Gets the positions.
     *
     * @return     The positions.
     */
    inline const auto& get_positions() const {
        return this->positions;
    }

    /**
     * @brief      Gets the velocities.
     *
     * @return     The velocities.
     */
    inline const auto& get_velocities() const {
        return this->velocities;
    }

    /**
     * @brief      Gets the total energy
     *
     * @return     The total energy
     */
    inline double get_etot() const {
        return this->etot;
    }

    /**
     * @brief      Gets the kinetic energy
     *
     * @return     The kinetic energy
     */
    inline double get_ekin() const {
        return this->ekin;
    }

    /**
     * @brief      Get the potential energy
     *
     * @return     The potential energy
     */
    inline double get_epot() const {
        return this->epot;
    }

    /**
     * @brief      Write to movie file for current state
     *
     * @param[in]  moviefile  Path to the movie file
     * @param[in]  create     Whether to truncate (true) or to append (false)
     */
    void write_to_movie_file(const std::string& moviefile, bool create = false);

private:
    /**
     * @brief      Initialize unit cell
     */
    void initialize();

    /**
     * @brief      Get number from normal distribution between -1 and 1
     *
     * @return     Random number
     */
    double get_gauss() const;

    /**
     * @brief      Calculates the forces.
     */
    void calculate_forces();

    /**
     * @brief      Update velocities by half a timestep
     *
     * @param[in]  dt    Timestep (is divided by half INSIDE the function)
     */
    void update_velocities_half_dt(double dt);

    /**
     * @brief      Update positions by timestep dt
     *
     * @param[in]  dt    Timestep
     */
    void update_positions(double dt);

    /**
     * @brief      Applies periodic boundary conditions to the positions
     */
    void apply_boundary_conditions();

    /**
     * @brief      Perform velocity scaling using Berendsen thermostat
     *
     * @param[in]  dt    Timestep
     */
    void apply_berendsen_thermostat(double dt);

    /**
     * @brief      Build neighbor list
     */
    void build_neighbor_list();

    /**
     * @brief      Check whether neighbor list needs to be updated
     *
     * @return     Whether neighbor list needs to be updated
     */
    bool check_update_neighbor_list() const;
};

#endif // LENNARDJONESSIMULATION_H
