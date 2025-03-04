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

#include "structure.h"

/**
 * @brief      Constructs a new instance.
 */
Structure::Structure(const MatrixUnitcell& _unitcell) :
unitcell(_unitcell) {

}

/**
 * Set positions
 */
void Structure::set_particle_positions(const std::vector<glm::dvec3>& _positions) {
    this->positions = _positions;
}

/**
 * Set velocities
 */
void Structure::set_particle_velocities(const std::vector<glm::dvec3>& _velocities) {
    this->velocities = _velocities;
}

/**
 * @brief      Get the centering vector
 *
 * @return     Vector that puts unitcell at the origin
 */
 QVector3D Structure::get_center_vector() const {
    auto ctr = this->unitcell.transpose() * VectorPosition::Ones() * 0.5;
    return QVector3D(-ctr(0), -ctr(1), -ctr(2));
}

/**
 * @brief      Gets the unitcell matrix.
 *
 * @param[in]  matrix  The matrix
 *
 * @return     The casted matrix.
 */
QMatrix3x3 Structure::get_matrix3x3(const MatrixUnitcell& matrix) const {
    std::vector<float> values(9, 0.0);
    for(unsigned int i=0; i<3; i++) {
        for(unsigned int j=0; j<3; j++) {
            values[i*3+j] = (float)matrix(i,j);
        }
    }
    return QMatrix3x3(&values[0]);
}