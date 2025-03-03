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
 * @brief      Constructs a new instance.
 */
Structure::Structure(unsigned int elnr) {
    this->unitcell = MatrixUnitcell::Identity() * 2.5f;
    this->atoms.emplace_back(elnr, 0.0, 0.0, 0.0);
    this->center();
}

/**
 * Set positions
 */
void Structure::set_particle_positions(const std::vector<glm::dvec3>& _positions) {
    this->atoms.resize(_positions.size());
    for(unsigned int i=0; i<_positions.size(); i++) {
        atoms[i].x = _positions[i][0];
        atoms[i].y = _positions[i][1];
        atoms[i].z = _positions[i][2];
        atoms[i].atnr = 1;
    }
}

/**
 * Set velocities
 */
void Structure::set_particle_velocities(const std::vector<glm::dvec3>& _velocities) {
    this->atoms.resize(_velocities.size());
    for(unsigned int i=0; i<_velocities.size(); i++) {
        atoms[i].vx = _velocities[i][0];
        atoms[i].vy = _velocities[i][1];
        atoms[i].vz = _velocities[i][2];
    }
}

/**
 * @brief      Add an atom to the structure including forces
 *
 * @param[in]  atnr  Atom number
 * @param[in]  x     x coordinate
 * @param[in]  y     y coordinate
 * @param[in]  z     z coordinate
 */
void Structure::add_atom(unsigned int atnr, double x, double y, double z) {
    this->atoms.emplace_back(atnr, x, y, z);
}

/**
 * @brief      Center the structure at the origin
 */
void Structure::center() {
    double sumx = 0.0;
    double sumy = 0.0;
    double sumz = 0.0;

    #pragma omp parallel for reduction(+: sumx)
    for(unsigned int i=0; i<this->atoms.size(); i++) {
        sumx += this->atoms[i].x;
    }

    #pragma omp parallel for reduction(+: sumy)
    for(unsigned int i=0; i<this->atoms.size(); i++) {
        sumy += this->atoms[i].y;
    }

    #pragma omp parallel for reduction(+: sumz)
    for(unsigned int i=0; i<this->atoms.size(); i++) {
        sumz += this->atoms[i].z;
    }

    sumx /= (float)this->atoms.size();
    sumy /= (float)this->atoms.size();
    sumz /= (float)this->atoms.size();

    auto cv = get_center_vector();

    #pragma omp parallel for
    for(unsigned int i=0; i<this->atoms.size(); i++) {
        this->atoms[i].x -= sumx + cv[0];
        this->atoms[i].y -= sumy + cv[1];
        this->atoms[i].z -= sumz + cv[2];
    }
}

/**
 * @brief      Get the largest distance from the origin
 *
 * @return     The largest distance.
 */
QVector3D Structure::get_largest_distance() const {
    unsigned int idx = 0;
    float dist;

    #pragma omp parallel for
    for(unsigned int i=0; i<this->atoms.size(); i++) {
        float vdist = this->atoms[i].get_pos_qtvec().lengthSquared();

        if(vdist > dist) {
            #pragma omp critical
            {
                dist = vdist;
                idx = i;
            }
        }
    }

    return this->atoms[idx].get_pos_qtvec();
}

/**
 * @brief      Gets the elements in this structure as a string
 *
 * @return     String holding comma seperated list of elements
 */
std::string Structure::get_elements_string() const {
    std::string result;

    for(const auto& item : this->element_types) {
        result += QString("%1 (%2); ").arg(QString(item.first.c_str())).arg(item.second).toStdString();
    }

    // remove last two characters
    result.pop_back();
    result.pop_back();

    return result;
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
 * @brief      Update data based on contents;
 */
void Structure::update() {
    this->count_elements();
}

/**
 * @brief      Count the number of elements
 */
void Structure::count_elements() {
    this->element_types.clear();

    for(const auto& atom : this->atoms) {
        std::string atomname = AtomSettings::get().get_name_from_elnr(atom.atnr);
        auto got = this->element_types.find(atomname);
        if(got != this->element_types.end()) {
            got->second++;
        } else {
            this->element_types.emplace(atomname, 1);
        }
    }
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