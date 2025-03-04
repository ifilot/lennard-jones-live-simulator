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

#include "threadintegrate.h"

/**
 * @brief Perform integration
 */
void ThreadIntegrate::run() {
    qDebug() << "Starting simulation";
    auto dt = this->ljsim->get_params()->get_param<double>("stepsize");
    this->iterator = 0;
    this->local_iterator = 0;
    const double fps = 60.0;
    const double ips = 2000;

    auto ttime = std::chrono::system_clock::now();

    while(true) {
        if(this->keeprunning && !this->graph_update_wait) {
            if(this->local_iterator < ips / fps) {
                this->ljsim->integrate(this->iterator, dt);
                this->local_iterator++;
                iterator++;
            }

            auto end = std::chrono::system_clock::now();
            std::chrono::duration<double, std::milli> elapsed_ms = end - ttime;

            if(elapsed_ms.count() > (1000.0 / fps)) {
                ttime = std::chrono::system_clock::now();
                this->local_iterator = 0;
                emit(signal_integration_step());
            }

            if(iterator % 1000 == 0) {
                emit(signal_ekin(iterator * dt, this->ljsim->get_ekin()));
                emit(signal_epot(iterator * dt, this->ljsim->get_epot()));
                emit(signal_etot(iterator * dt, this->ljsim->get_etot()));
                emit(signal_velocities());
                this->graph_update_wait = true;
            }
        } else {
            QThread::sleep(1); // prevents thread from crashing
        }
    }

    exec();
}

/**
 * @brief stop
 */
void ThreadIntegrate::stop() {
    this->keeprunning = false;
}

/**
 * @brief Pause/unpause the simulation
 */
void ThreadIntegrate::toggle_pause() {
    this->keeprunning = !this->keeprunning;
    qDebug() << this->keeprunning;
}

/**
 * @brief Receive slot to continue the simulation
 */
void ThreadIntegrate::slot_simulation_unlock() {
    this->graph_update_wait = false;
}
