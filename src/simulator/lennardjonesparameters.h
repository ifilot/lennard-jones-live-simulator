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

#ifndef LENNARDJONESPARAMETERS_H
#define LENNARDJONESPARAMETERS_H

#include <QMap>
#include <QString>
#include <QVariant>
#include <stdexcept>

/**
 * @brief      Class for Lennard-Jones parameters using Qt.
 */
class LennardJonesParameters {
private:
    QMap<QString, QVariant> params; // map containing parameters and their values

public:
    /**
     * @brief      Constructs the object.
     */
    LennardJonesParameters();

    /**
     * @brief      Get a parameter value as a specific type.
     *
     * @param[in]  name  The name of the parameter
     *
     * @tparam     T     The type to convert the parameter value to
     *
     * @return     The parameter value converted to type T.
     */
    template<typename T>
    T get_param(const QString& name) const {
        if (params.contains(name)) {
            return params.value(name).value<T>();
        } else {
            throw std::runtime_error("Could not find parameter: " + name.toStdString());
        }
    }

    /**
     * @brief      Sets the parameter.
     *
     * @param[in]  name  Name of the parameter
     * @param[in]  value The parameter value
     */
    inline void set_param(const QString& name, const QVariant& value) {
        params.insert(name, value);
    }
};

#endif // LENNARDJONESPARAMETERS_H