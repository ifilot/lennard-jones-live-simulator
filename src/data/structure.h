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

#include "atom_settings.h"
#include "matrixmath.h"
#include "atom.h"

/**
 * @brief      This class describes a chemical structure.
 */
class Structure {

private:
    std::vector<Atom> atoms;            // atoms in the structure

    double energy = 0.0;                // energy of the structure (if known, zero otherwise)

    MatrixUnitcell unitcell;            // matrix describing the unit cell
    std::vector<double> radii;          // radii of the atoms

    std::unordered_map<std::string, unsigned int> element_types;    // elements present in the structure

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
     * @brief      Get all atoms from the structure
     *
     * @return     The atoms.
     */
    inline const auto& get_atoms() const {
        return this->atoms;
    }

    /**
     * @brief      Get specific atom
     *
     * @param[in]  idx   The index
     *
     * @return     The atom.
     */
    inline const Atom& get_atom(unsigned int idx) const {
        return this->atoms[idx];
    }

    /**
     * @brief      Gets the unitcell.
     *
     * @return     The unitcell.
     */
    inline const auto& get_unitcell() const {
        return this->unitcell;
    }

    /**
     * @brief      Gets the atomic radius.
     *
     * @param[in]  idx   The index
     *
     * @return     The radius.
     */
    inline double get_radius(unsigned int idx) const {
        return this->radii[idx];
    }

    /**
     * @brief      Add an atom to the structure
     *
     * @param[in]  atnr  Atom number
     * @param[in]  x     x coordinate
     * @param[in]  y     y coordinate
     * @param[in]  z     z coordinate
     */
    void add_atom(unsigned int atnr, double x, double y, double z);

    /**
     * @brief      Gets the total number of atoms.
     *
     * @return     The number of atoms
     */
    inline size_t get_nr_atoms() const {
        return this->atoms.size();
    }

    ~Structure() {
        qDebug() << "Deleting structure ("
                 << QString("0x%1").arg((size_t)this, 0, 16)
                 << "; " << this->atoms.size() << " atoms ).";
    }

    //********************************************
    // [END BLOCK] DATA GETTERS AND SETTERS
    //********************************************

    /**
     * @brief      Center the structure at the origin
     */
    void center();

    /**
     * @brief      Get the largest distance from the origin
     *
     * @return     The largest distance.
     */
    QVector3D get_largest_distance() const;

    /**
     * @brief      Gets the elements in this structure as a string
     *
     * @return     String holding comma seperated list of elements
     */
    std::string get_elements_string() const;

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
     * @brief      Count the number of elements
     */
    void count_elements();

    /**
     * @brief      Construct the bonds
     */
    void construct_bonds();

    /**
     * @brief      Gets the unitcell matrix.
     *
     * @param[in]  matrix  The matrix
     *
     * @return     The casted matrix.
     */
    QMatrix3x3 get_matrix3x3(const MatrixUnitcell& matrix) const;
};
