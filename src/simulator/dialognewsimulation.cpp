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

#include "dialognewsimulation.h"

DialogNewSimulation::DialogNewSimulation(const std::shared_ptr<LennardJonesParameters>& _params, QWidget *parent) :
    QDialog(parent),
    params(_params) {

    QVBoxLayout* layout = new QVBoxLayout();
    this->setLayout(layout);


    QGridLayout* gridlayout = new QGridLayout();
    QWidget* gridcontainer = new QWidget();
    gridcontainer->setLayout(gridlayout);
    layout->addWidget(gridcontainer);

    this->setLayout(layout);

    for(int i=0; i<this->variables.size(); i++) {
        QLabel* label = new QLabel(variables[i]);
        gridlayout->addWidget(label, i, 0);
        this->spinboxes.push_back(new QDoubleSpinBox());
        if(i > 1) {
            this->spinboxes.back()->setMinimum(0);
            this->spinboxes.back()->setMaximum(5);
            this->spinboxes.back()->setDecimals(4);
        }
        this->spinboxes.back()->setValue(this->values[i]);
        gridlayout->addWidget(this->spinboxes.back(), i, 1);
    }

    // custom settings
    this->spinboxes[1]->setMinimum(1000);
    this->spinboxes[1]->setMaximum(10000);

    QHBoxLayout* buttonlayout = new QHBoxLayout();
    QWidget* buttoncontainer = new QWidget();
    layout->addWidget(buttoncontainer);
    buttoncontainer->setLayout(buttonlayout);

    QPushButton* button_ok = new QPushButton("OK");
    buttonlayout->addWidget(button_ok);

    QPushButton* button_cancel = new QPushButton("Cancel");
    buttonlayout->addWidget(button_cancel);

    this->setWindowTitle(tr("Build new simulation"));

    connect(button_ok, SIGNAL(released()), this, SLOT(slot_confirm()));
    connect(button_cancel, SIGNAL(released()), this, SLOT(slot_cancel()));
}


void DialogNewSimulation::slot_confirm() {
    for(int i=0; i<this->variables.size(); i++) {
        this->params->set_param(this->variables[i], tr("%1").arg(this->spinboxes[i]->value()));
    }
    this->done(QDialog::Accepted);
}

void DialogNewSimulation::slot_cancel() {
    this->done(QDialog::Rejected);
}
