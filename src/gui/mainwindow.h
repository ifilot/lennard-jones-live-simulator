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

#include <QMainWindow>
#include <QMenuBar>
#include <QMenu>
#include <QMessageBox>
#include <QFileDialog>
#include <QStatusBar>
#include <QString>
#include <QtWidgets/QApplication>
#include <QFileInfo>
#include <QMimeData>
#include <QTimer>
#include <QStringList>
#include <QLabel>
#include <memory>

#include "logwindow.h"
#include "anaglyph_widget.h"
#include "../simulator/lennardjonessimulation.h"
#include "../simulator/threadintegrate.h"
#include "../simulator/graphwidget.h"
#include "../simulator/dialognewsimulation.h"

#include "../config.h"

/**
 * @brief      Class for main window.
 */
class MainWindow : public QMainWindow {
    Q_OBJECT

private:
    AnaglyphWidget* anaglyph_widget;
    QLabel* statusbar_projection_icon;
    QTimer* statusbar_timer;

    // Lennard-Jones simulation
    std::shared_ptr<LennardJonesSimulation> ljsim;
    ThreadIntegrate* thread = nullptr;
    GraphWidget *graph_widget = nullptr;

    // control interface
    QWidget *control_interface;
    QPushButton *button_toggle_axis;
    QPushButton *button_rotate_model;
    QPushButton *button_simulation_pause;
    QPushButton *button_relative_coloring;

    // storage for log messages
    std::shared_ptr<QStringList> log_messages;

    // window for log messages
    std::unique_ptr<LogWindow> log_window;

    QTimer *rotation_timer;

    bool flag_rotate = false;

public:
    /**
     * @brief      Constructs the object.
     */
    MainWindow(const std::shared_ptr<QStringList> _log_messages,
               QWidget *parent = nullptr);
    ~MainWindow() override;

protected:
    void moveEvent(QMoveEvent*) Q_DECL_OVERRIDE;

private slots:
    // /**
    //  * @brief      Open a new object file
    //  */
    // void open();

    // /**
    //  * @brief      Open a new object file
    //  */
    // void save();

    /**
     * @brief      Close the application
     */
    void exit();

    /**
     * @brief      Display about menu
     */
    void about();

    /**
     * @brief      Set stereo projection
     */
    void set_stereo(QString fragment_shader);

    /**
     * @brief      Show a message on the statusbar
     *
     * @param[in]  message  The message
     */
    void show_message_statusbar(const QString& message);

    /**
     * @brief      Clear statusbar when timer runs out
     */
    void statusbar_timeout();

    /**
     * @brief Show an about window
     */
    void slot_debug_log();

    /**
     * @brief Construct a new simulation
     */
     void new_simulation();

    // viewport interfacing
    void toggle_world_axes();

    /*
     * Enable / disable scene rotation
     */
    void toggle_rotation();

    /**
     * @brief Perform integration step update call
     */
    void integration_step();

    /**
     * @brief Transmit velocities from simulation to Graph Widget
     */
    void slot_transmit_velocities();

    /**
    * @brief Handles rotation timer
    */
    void rotation_timer_trigger();

private:
    void stop_thread();

    /**
     * @brief Build the Lennard-Jones simulation
     */
    void build_simulation(const std::shared_ptr<LennardJonesParameters>& params);

     /**
      * @brief Build the control interface
      */
    void build_control_interface();

    /**
     * Create the dropdown menu and associate all actions
     */
    void create_dropdown_menu();

    /**
      * @brief Sets up the main interface layout for the application.
      *
      * This function initializes and arranges the widgets within the main window.
      * It creates a central container with a horizontal layout, which consists of:
      * - A left vertical container holding the AnaglyphWidget (maximized horizontally)
      *   and the control interface (minimized as much as possible).
      * - A GraphWidget on the right.
      */
    void set_widgets();

    /**
    * Setup a simulation
    */
    void setup_simulation(const std::shared_ptr<LennardJonesParameters>& params);
};
