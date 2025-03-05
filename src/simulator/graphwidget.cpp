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

#include "graphwidget.h"

GraphWidget::GraphWidget(QWidget *parent) : QWidget(parent) {
    QVBoxLayout* layout = new QVBoxLayout();
    this->setLayout(layout);
    this->setMinimumWidth(400);

    // data kinetic energy
    this->chartview_kinetic_energy = new QChartView();
    this->chartview_kinetic_energy->setFixedHeight(200);
    layout->addWidget(this->chartview_kinetic_energy);
    this->chart_kinetic_energy = new QChart();
    this->chartview_kinetic_energy->setChart(this->chart_kinetic_energy);

    this->chartview_potential_energy = new QChartView();
    this->chartview_potential_energy->setFixedHeight(200);
    layout->addWidget(this->chartview_potential_energy);
    this->chart_potential_energy = new QChart();
    this->chartview_potential_energy->setChart(this->chart_potential_energy);

    this->chartview_total_energy = new QChartView();
    this->chartview_total_energy->setFixedHeight(200);
    layout->addWidget(this->chartview_total_energy);
    this->chart_total_energy = new QChart();
    this->chartview_total_energy->setChart(this->chart_total_energy);

    this->chartview_particle_speed = new QChartView();
    layout->addWidget(this->chartview_particle_speed);
    this->chart_particle_speed = new QChart();
    this->chartview_particle_speed->setChart(this->chart_particle_speed);

    this->build_empty_graphs();
}

void GraphWidget::set_params(double temperature, double m, int numparts) {
    this->series_maxwell_boltzmann->clear();
    this->series_maxwell_boltzmann->setPointsVisible();
    this->series_maxwell_boltzmann->setName("Maxwell-Boltzmann Speed Distribution");
    double maxnumparts = 0;
    double binsize = MAXSPEED / (double)this->NUMBINS;
    for(unsigned int i=0; i<this->NUMBINS; i++) {
        double v = binsize * i;
        double f = 4.0 * M_PI * std::pow(m / (2.0 * M_PI * temperature), 3./2.) * v * v * std::exp(-m * v * v / (2.0 * temperature));
        double np = f * numparts * binsize;
        maxnumparts = std::max(np, maxnumparts);

        this->series_maxwell_boltzmann->append(v/MAXSPEED, np);
    }

    this->chart_particle_speed->axes(Qt::Vertical).first()->setRange(0, maxnumparts * 1.1);
    qobject_cast<QValueAxis*>(this->chart_particle_speed->axes(Qt::Vertical).first())->applyNiceNumbers();
    chart_particle_speed->update();
}

void GraphWidget::build_empty_graphs() {
    // data kinetic energy
    this->build_extended_graph(&this->chart_kinetic_energy,
                               &this->series_data_kinetic_energy,
                               &this->series_data_average_kinetic_energy,
                               &this->series_data_stdlower_kinetic_energy,
                               &this->series_data_stdupper_kinetic_energy,
                               &this->series_area_kinetic_energy,
                               "Kinetic energy");

    // data potential energy
    this->build_extended_graph(&this->chart_potential_energy,
                               &this->series_data_potential_energy,
                               &this->series_data_average_potential_energy,
                               &this->series_data_stdlower_potential_energy,
                               &this->series_data_stdupper_potential_energy,
                               &this->series_area_potential_energy,
                               "Potential energy");

    // data total energy
    this->build_extended_graph(&this->chart_total_energy,
                               &this->series_data_total_energy,
                               &this->series_data_average_total_energy,
                               &this->series_data_stdlower_total_energy,
                               &this->series_data_stdupper_total_energy,
                               &this->series_area_total_energy,
                               "Total energy");

    // data particle speed
    this->set_bar_particle_speed = new QBarSet("Particle speed");
    this->series_bar_particle_speed = new QBarSeries();
    this->series_maxwell_boltzmann = new QLineSeries();
    this->series_bar_particle_speed->append(this->set_bar_particle_speed);
    this->chart_particle_speed->addSeries(this->series_bar_particle_speed);
    this->chart_particle_speed->addSeries(this->series_maxwell_boltzmann);
    this->chart_particle_speed->createDefaultAxes();

    this->chart_particle_speed->axes(Qt::Horizontal).first()->setRange(0,NUMBINS);
    this->chart_particle_speed->axes(Qt::Horizontal).first()->hide();
    this->chart_particle_speed->setTitle("Histogram: particle speed");

    // specify axes
    QValueAxis *axisX = new QValueAxis();
    axisX->setRange(0,MAXSPEED);
    this->series_maxwell_boltzmann->attachAxis(axisX);
    this->chart_particle_speed->addAxis(axisX, Qt::AlignBottom);
}

void GraphWidget::reset_graphs() {
    this->kinetic_energy.clear();
    this->chart_kinetic_energy->removeAllSeries();

    this->potential_energy.clear();
    this->chart_potential_energy->removeAllSeries();

    this->total_energy.clear();
    this->chart_total_energy->removeAllSeries();

    this->chart_particle_speed->removeAllSeries();

    this->build_empty_graphs();
    this->no_speed_axis = true;
}

void GraphWidget::add_data_item_kinetic_energy(double t, double ekin) {
    this->series_data_kinetic_energy->append(t, ekin);
    this->chart_kinetic_energy->axes(Qt::Horizontal).first()->setRange(0, t+1);

    this->kinetic_energy.push_back(ekin);

    double emin = *std::min_element(this->kinetic_energy.begin(), this->kinetic_energy.end());
    double emax = *std::max_element(this->kinetic_energy.begin(), this->kinetic_energy.end());
    double avg = std::accumulate(this->kinetic_energy.begin(), this->kinetic_energy.end(), 0) / (double)this->kinetic_energy.size();
    double variance = this->calculate_variance(this->kinetic_energy);

    this->series_data_average_kinetic_energy->replace(0, 0, avg);
    this->series_data_average_kinetic_energy->replace(1, t, avg);
    this->series_data_stdlower_kinetic_energy->replace(0, 0, avg - std::sqrt(variance));
    this->series_data_stdlower_kinetic_energy->replace(1, t, avg - std::sqrt(variance));
    this->series_data_stdupper_kinetic_energy->replace(0, 0, avg + std::sqrt(variance));
    this->series_data_stdupper_kinetic_energy->replace(1, t, avg + std::sqrt(variance));
    
    this->chart_kinetic_energy->axes(Qt::Vertical).first()->setRange(emin, emax);
    qobject_cast<QValueAxis*>(this->chart_kinetic_energy->axes(Qt::Vertical).first())->applyNiceNumbers();
}

void GraphWidget::add_data_item_potential_energy(double t, double epot) {
    this->series_data_potential_energy->append(t, epot);
    this->chart_potential_energy->axes(Qt::Horizontal).first()->setRange(0, t+1);

    this->potential_energy.push_back(epot);

    double emin = *std::min_element(this->potential_energy.begin(), this->potential_energy.end());
    double emax = *std::max_element(this->potential_energy.begin(), this->potential_energy.end());
    double avg = std::accumulate(this->potential_energy.begin(), this->potential_energy.end(), 0) / (double)this->potential_energy.size();
    double variance = this->calculate_variance(this->potential_energy);

    this->series_data_average_potential_energy->replace(0, 0, avg);
    this->series_data_average_potential_energy->replace(1, t, avg);
    this->series_data_stdlower_potential_energy->replace(0, 0, avg - std::sqrt(variance));
    this->series_data_stdlower_potential_energy->replace(1, t, avg - std::sqrt(variance));
    this->series_data_stdupper_potential_energy->replace(0, 0, avg + std::sqrt(variance));
    this->series_data_stdupper_potential_energy->replace(1, t, avg + std::sqrt(variance));

    this->chart_potential_energy->axes(Qt::Vertical).first()->setRange(emin, emax);
    qobject_cast<QValueAxis*>(this->chart_potential_energy->axes(Qt::Vertical).first())->applyNiceNumbers();
}

void GraphWidget::add_data_item_total_energy(double t, double epot) {
    this->series_data_total_energy->append(t, epot);
    this->chart_total_energy->axes(Qt::Horizontal).first()->setRange(0, t+1);

    this->total_energy.push_back(epot);

    double emin = *std::min_element(this->total_energy.begin(), this->total_energy.end());
    double emax = *std::max_element(this->total_energy.begin(), this->total_energy.end());
    double avg = std::accumulate(this->total_energy.begin(), this->total_energy.end(), 0) / (double)this->total_energy.size();
    double variance = this->calculate_variance(this->total_energy);

    this->series_data_average_total_energy->replace(0, 0, avg);
    this->series_data_average_total_energy->replace(1, t, avg);
    this->series_data_stdlower_total_energy->replace(0, 0, avg - std::sqrt(variance));
    this->series_data_stdlower_total_energy->replace(1, t, avg - std::sqrt(variance));
    this->series_data_stdupper_total_energy->replace(0, 0, avg + std::sqrt(variance));
    this->series_data_stdupper_total_energy->replace(1, t, avg + std::sqrt(variance));

    this->chart_total_energy->axes(Qt::Vertical).first()->setRange(emin, emax);
    qobject_cast<QValueAxis*>(this->chart_total_energy->axes(Qt::Vertical).first())->applyNiceNumbers();
}

void GraphWidget::set_particle_speed(const std::vector<glm::dvec3>& velocities) {
    if(this->no_speed_axis) {
        for(unsigned int i=0; i<this->NUMBINS; i++) {
            this->set_bar_particle_speed->append(0.0);
        }

        QValueAxis *axisX = new QValueAxis();
        axisX->setRange(0, this->NUMBINS);
        chart_particle_speed->addAxis(axisX, Qt::AlignBottom);
        this->series_bar_particle_speed->attachAxis(axisX);
        axisX->hide();
        this->series_bar_particle_speed->setBarWidth(1.0);

        this->no_speed_axis = false;
    }

    // create histogram
    for(unsigned int i=0; i<NUMBINS; i++) {
        this->set_bar_particle_speed->replace(i, 0);
    }

    for(unsigned int i=0; i<velocities.size(); i++) {
        double speed = glm::length(velocities[i]);
        int bin = std::min((int)(speed / MAXSPEED * NUMBINS), this->set_bar_particle_speed->count()-1);
        this->set_bar_particle_speed->replace(bin, this->set_bar_particle_speed->at(bin)+1);
    }

    emit(signal_simulation_continue());
}

double GraphWidget::calculate_variance(const std::vector<double>& values) {
    double avg = std::accumulate(values.begin(), values.end(), 0) / (double)this->kinetic_energy.size();
    double variance = 0.0;
    for(auto val : values) {
        variance += (val - avg) * (val - avg);
    }
    variance /= (double)(values.size() - 1);

    if(!qIsFinite(variance)) {
        return 0.0;
    }

    return variance;
}

void GraphWidget::build_extended_graph(QChart** chart,
                                       QLineSeries** base_value,
                                       QLineSeries** average_value,
                                       QLineSeries** std_lower,
                                       QLineSeries** std_upper,
                                       QAreaSeries** area,
                                       const QString& title) {
    *base_value = new QLineSeries();
    *average_value = new QLineSeries();
    (*average_value)->append(0,0);
    (*average_value)->append(0,0);
    *std_lower = new QLineSeries();
    (*std_lower)->append(0,0);
    (*std_lower)->append(0,0);
    (*std_upper) = new QLineSeries();
    (*std_upper)->append(0,0);
    (*std_upper)->append(0,0);

    *area = new QAreaSeries(*std_lower, *std_upper);
    QPen pen(0x222222);
    pen.setStyle(Qt::DotLine);
    pen.setColor(QColor(0xFF,0x7F,0x0e,80));
    (*area)->setPen(pen);
    QLinearGradient gradient(QPointF(0, 0), QPointF(0, 1));
    gradient.setColorAt(0.0, QColor(0xFF,0x7F,0x0e,50));
    gradient.setColorAt(1.0, QColor(0xFF,0x7F,0x0e,50));
    gradient.setCoordinateMode(QGradient::ObjectBoundingMode);
    (*area)->setBrush(gradient);

    QPen pen2(QColor(0,0,0));
    pen2.setColor(QColor(0xFF,0x7F,0x0e,100));
    pen2.setStyle(Qt::DashLine);
    (*average_value)->setPen(pen2);
    (*base_value)->setPointsVisible();

    (*chart)->addSeries(*area);
    (*chart)->addSeries(*average_value);
    (*chart)->addSeries(*base_value);
    (*chart)->setTitle(title);
    (*chart)->legend()->hide();
    (*chart)->createDefaultAxes();
    (*chart)->axes(Qt::Horizontal).first()->setTitleText("Time [a.u.]");
    (*chart)->axes(Qt::Vertical).first()->setTitleText("Energy [a.u.]");
}
