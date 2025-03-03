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

#ifndef GRAPHWIDGET_H
#define GRAPHWIDGET_H

#include <QObject>
#include <QWidget>
#include <QtCharts>
#include <QVBoxLayout>
#include <QBarSet>

#include <numeric>

#include <glm/glm.hpp>

class GraphWidget : public QWidget
{
    Q_OBJECT

private:
    QVBoxLayout* layout;

    QChartView* chartview_kinetic_energy;
    QChart* chart_kinetic_energy;
    QLineSeries* series_data_kinetic_energy;
    QLineSeries* series_data_average_kinetic_energy;
    QLineSeries* series_data_stdlower_kinetic_energy;
    QLineSeries* series_data_stdupper_kinetic_energy;
    QAreaSeries* series_area_kinetic_energy;
    std::vector<double> kinetic_energy;

    QChartView* chartview_potential_energy;
    QChart* chart_potential_energy;
    QLineSeries* series_data_potential_energy;
    QLineSeries* series_data_average_potential_energy;
    QLineSeries* series_data_stdlower_potential_energy;
    QLineSeries* series_data_stdupper_potential_energy;
    QAreaSeries* series_area_potential_energy;
    std::vector<double> potential_energy;

    QChartView* chartview_total_energy;
    QChart* chart_total_energy;
    QLineSeries* series_data_total_energy;
    QLineSeries* series_data_average_total_energy;
    QLineSeries* series_data_stdlower_total_energy;
    QLineSeries* series_data_stdupper_total_energy;
    QAreaSeries* series_area_total_energy;
    std::vector<double> total_energy;

    QChartView* chartview_particle_speed;
    QChart* chart_particle_speed;
    QBarSeries* series_bar_particle_speed;
    QBarSet* set_bar_particle_speed;

    static const unsigned int NUMBINS = 100;
    static constexpr double MAXSPEED = 10.0;
    QLineSeries* series_maxwell_boltzmann;

    bool no_speed_axis = true;

public:
    explicit GraphWidget(QWidget *parent = nullptr);

    void set_params(double temperature, double m, int numparts);

private:
    void build_empty_graphs();

    double calculate_variance(const std::vector<double>& values);

    void build_extended_graph(QChart** chart,
                              QLineSeries** base_value,
                              QLineSeries** average_value,
                              QLineSeries** std_lower,
                              QLineSeries** std_upper,
                              QAreaSeries** area,
                              const QString& title);
signals:
    void signal_simulation_continue();

public slots:
    void reset_graphs();

    void add_data_item_kinetic_energy(double t, double ekin);

    void add_data_item_potential_energy(double t, double ekin);

    void add_data_item_total_energy(double t, double etot);

    void set_particle_speed(const std::vector<glm::dvec3>& velocities);
};

#endif // GRAPHWIDGET_H
