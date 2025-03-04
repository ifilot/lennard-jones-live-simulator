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

#pragma once

#include <QDebug>
#include <QVector3D>
#include <QMatrix4x4>
#include <QGenericMatrix>
#include <vector>
#include <QString>
#include <glm/glm.hpp>

#include "matrixmath.h"

/**
 * @brief      This class describes a chemical structure.
 */
class Structure {

private:
    std::vector<glm::dvec3> positions;
    std::vector<glm::dvec3> velocities;
    MatrixUnitcell unitcell;

public:
    /**
     * @brief      Constructs a new instance.
     */
    Structure(const MatrixUnitcell& unitcell);

    /**
     * @brief      Constructs a new instance.
     */
    Structure(unsigned int elnr);

    /**
     * Set positions
     */
    void set_particle_positions(const std::vector<glm::dvec3>& _positions);

    /**
     * Set positions
     */
     void set_particle_velocities(const std::vector<glm::dvec3>& _velocities);

    /**
     * Get positions
     */
    const auto& get_positions() const {
        return this->positions;
    }

    /**
     * Get velocities
     */
    const auto& get_velocities() const {
        return this->velocities;
    }

    /**
     * @brief      Gets the unitcell.
     *
     * @return     The unitcell.
     */
    inline const auto& get_unitcell() const {
        return this->unitcell;
    }

    //********************************************
    // [END BLOCK] DATA GETTERS AND SETTERS
    //********************************************

    /**
     * @brief      Update data based on contents;
     */
    void update();

    /**
     * @brief      Get the centering vector
     *
     * @return     Vector that puts unitcell at the origin
     */
    QVector3D get_center_vector() const;

private:
    /**
     * @brief      Gets the unitcell matrix.
     *
     * @param[in]  matrix  The matrix
     *
     * @return     The casted matrix.
     */
    QMatrix3x3 get_matrix3x3(const MatrixUnitcell& matrix) const;
};
