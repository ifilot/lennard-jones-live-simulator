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

#ifndef THREADINTEGRATE_H
#define THREADINTEGRATE_H

#include <QThread>
#include <QObject>
#include <memory>
#include <QTimer>
#include <chrono>

#include "lennardjonessimulation.h"

class ThreadIntegrate : public QThread {

    Q_OBJECT

private:
    std::shared_ptr<LennardJonesSimulation> ljsim;
    bool keeprunning = true;
    bool graph_update_wait = false;
    unsigned int iterator;
    unsigned int local_iterator;

public:
    /**
     * @brief Default constructor
     * @param Pointer to Lennard Jones Simulation
     */
    ThreadIntegrate(const std::shared_ptr<LennardJonesSimulation> _ljsim) :
        ljsim(_ljsim) {}

    /**
     * @brief Perform integration
     */
    void run();

    inline void set_simulation(const std::shared_ptr<LennardJonesSimulation> _ljsim) {
        this->iterator = 0;
        this->local_iterator = 0;
        this->ljsim = _ljsim;
    }

public slots:
    /**
     * @brief stop
     */
    void stop();

    /**
     * @brief Pause/unpause the simulation
     */
    void toggle_pause();

    /**
     * @brief Wait for response call from GraphWidget to continue the simulation
     */
    void slot_simulation_unlock();

signals:
    /**
     * @brief Signal that a new integration step has been made
     */
    void signal_integration_step();

    /**
     * @brief Add new kinetic energy data point
     * @param ekin
     */
    void signal_ekin(double t, double ekin);

    /**
     * @brief Add new potential energy data point
     * @param ekin
     */
    void signal_epot(double t, double epot);

    /**
     * @brief Add new total energy data point
     * @param ekin
     */
    void signal_etot(double t, double etot);

    /**
     * @brief Signal for new set of particle velocities
     * @param velocities
     */
    void signal_velocities();
};

#endif // THREADINTEGRATE_H
