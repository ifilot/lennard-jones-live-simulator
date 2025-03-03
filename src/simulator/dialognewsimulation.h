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

#ifndef DIALOGNEWSIMULATION_H
#define DIALOGNEWSIMULATION_H

#include <QDialog>
#include <QObject>
#include <QGridLayout>
#include <QLabel>
#include <QDoubleSpinBox>
#include <QVBoxLayout>
#include <QPushButton>

#include "lennardjonesparameters.h"

class DialogNewSimulation : public QDialog {

    Q_OBJECT

private:
    QVector<QDoubleSpinBox*> spinboxes;
    std::shared_ptr<LennardJonesParameters> params;

    const QStringList variables = {
        "cell_length",
        "nr_particles",
        "kT",
        "mass",
        "sigma",
        "rcut",
        "epsilon",
        "tau",
        "shell",
        "stepsize"
    };

    const QVector<double> values = {
        14.938,
        1000,
        1.0,
        1.0,
        1.0,
        2.5,
        1.0,
        0.5,
        0.4,
        0.001
    };

public:
    DialogNewSimulation(const std::shared_ptr<LennardJonesParameters>& _params, QWidget* parent = nullptr);

private slots:
    void slot_confirm();

    void slot_cancel();
};

#endif // DIALOGNEWSIMULATION_H
