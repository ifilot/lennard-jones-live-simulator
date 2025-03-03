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

#include "mainwindow.h"

/**
 * @brief      Class for main window.
 */
MainWindow::MainWindow(const std::shared_ptr<QStringList> _log_messages,
                       QWidget *parent)
    : QMainWindow(parent),
    log_messages(_log_messages) {

    qDebug() << "Constructing Main Window";

    // log window
    this->log_window = std::make_unique<LogWindow>(this->log_messages);

    // create and activate the dropdown menu
    this->create_dropdown_menu();

    // create all widgets
    this->set_widgets();
    
    this->statusbar_timer = new QTimer(this);
    connect(this->statusbar_timer, SIGNAL(timeout()), this, SLOT(statusbar_timeout()));

    // add projection icon to status bar
    this->statusbar_projection_icon = new QLabel();
    this->statusbar_projection_icon->setFixedSize(16, 16);
    this->statusbar_projection_icon->setPixmap(QPixmap(":/assets/icon/two_dimensional_32.png").scaled(16, 16));
    statusBar()->addPermanentWidget(this->statusbar_projection_icon);

    // display status message
    statusBar()->showMessage(QString(PROGRAM_NAME) + " " + QString(PROGRAM_VERSION));
    this->statusbar_timer->start(1000);

    // set icon
    setWindowIcon(QIcon(":/assets/icon/atom_architect_256.ico"));

    // set Window properties
    this->setWindowTitle(QString(PROGRAM_NAME) + " " + QString(PROGRAM_VERSION));
    this->resize(1280,640);
}

/**
 * @brief      Close the application
 */
void MainWindow::exit() {
    QMessageBox msgBox;
    msgBox.setText("Exit program.");
    msgBox.setStyleSheet("QLabel{min-width: 350px; font-weight: normal;}");
    msgBox.setInformativeText("Are you sure you want to quit? Your progress will be <b>unsaved</b>.");
    msgBox.setStandardButtons(QMessageBox::Ok | QMessageBox::Cancel);
    msgBox.setDefaultButton(QMessageBox::Cancel);
    msgBox.setWindowIcon(QIcon(":/assets/icon/atom_architect_256.ico"));
    int ret = msgBox.exec();

    switch (ret) {
      case QMessageBox::Ok:
          QApplication::quit();
          return;
      case QMessageBox::Cancel:
          return;
      default:
          // should never be reached
          return;
    }
}

/**
 * @brief      Display about menu
 */
void MainWindow::about() {
    QMessageBox message_box;
    message_box.setStyleSheet("QLabel{min-width: 250px; font-weight: normal;}");
    message_box.setText(PROGRAM_NAME
                        " version "
                        PROGRAM_VERSION
                        ".\n\nAuthor:\nIvo Filot <i.a.w.filot@tue.nl>\n\n"
                        PROGRAM_NAME " is licensed under the GPLv3 license.\n\n"
                        PROGRAM_NAME " is dynamically linked to Qt, which is licensed under LGPLv3.\n");
    message_box.setIcon(QMessageBox::Information);
    message_box.setWindowTitle("About " PROGRAM_NAME);
    message_box.setWindowIcon(QIcon(":/assets/icon/atom_architect_256.ico"));
    message_box.setIconPixmap(QPixmap(":/assets/icon/atom_architect_256.ico"));
    message_box.exec();
}

/**
 * @brief      Set stereo projection
 */
void MainWindow::set_stereo(QString stereo_name) {
    this->anaglyph_widget->set_stereo(stereo_name);

    // set icon
    QString icon_name;
    if (stereo_name.startsWith("stereo")) {
        icon_name = ":/assets/icon/" + stereo_name.replace("stereo_", "") + "_32.png";
    } else {
        icon_name = ":/assets/icon/two_dimensional_32.png";
    }

    QPixmap pixmap(icon_name);
    this->statusbar_projection_icon->setPixmap(pixmap.scaled(16, 16));
}

/**
 * @brief      Handle windows move event
 *
 * Updates anaglyph window
 *
 * @param      event  MoveEvent
 */
void MainWindow::moveEvent(QMoveEvent* event) {
    this->anaglyph_widget->window_move_event();

    // this is just to silence a warning message, we are not using event here, but this is the callback form
    (void)event;
}

/**
 * @brief      Clear statusbar when timer runs out
 */
void MainWindow::statusbar_timeout() {
    statusBar()->showMessage("");
}

/**
 * @brief      Show a message on the statusbar
 *
 * @param[in]  message  The message
 */
void MainWindow::show_message_statusbar(const QString& message) {
    statusBar()->showMessage(message);
    this->statusbar_timer->start(1000);
}

/**
 * @brief Show an about window
 */
void MainWindow::slot_debug_log() {
    this->log_window->show();
}

void MainWindow::toggle_world_axes() {
    // this->anaglyph_widget->toggle_world_axes();

    // if(this->anaglyph_widget->get_flag_world_axes()) {
    //     this->button_toggle_axis->setIcon(QIcon(":/assets/icon/axes_32.png"));
    // } else {
    //     this->button_toggle_axis->setIcon(QIcon(":/assets/icon/axes_gray_32.png"));
    // }
}

void MainWindow::toggle_rotation() {
    // this->anaglyph_widget->toggle_rotation_z();

    // if(this->anaglyph_widget->get_flag_rotation()) {
    //     this->button_rotate_model->setIcon(QIcon(":/assets/icon/rotation_32.png"));
    // } else {
    //     this->button_rotate_model->setIcon(QIcon(":/assets/icon/rotation_gray_32.png"));
    // }
}

/**
 * @brief Perform integration step update call
 */
 void MainWindow::integration_step() {
    // this->anaglyph_widget->set_particle_positions(this->ljsim->get_positions());
    // this->anaglyph_widget->set_particle_velocities(this->ljsim->get_velocities());
    // this->anaglyph_widget->update();
}

/**
 * @brief Transmit velocities from simulation to Graph Widget
 */
void MainWindow::slot_transmit_velocities() {
    this->graph_widget->set_particle_speed(this->ljsim->get_velocities());
}

/**
 * @brief Build the Lennard-Jones simulation
 */
 void MainWindow::build_simulation(const std::shared_ptr<LennardJonesParameters>& params) {
    qDebug() << "Build new simulation";
    this->ljsim = std::make_shared<LennardJonesSimulation>(params);
    auto positions = this->ljsim->get_positions();
    auto dimensions = this->ljsim->get_dims();
    // this->anaglyph_widget->set_unitcell_scale(dimensions[0]);
    // this->anaglyph_widget->set_particle_positions(positions);
    // this->anaglyph_widget->set_particle_velocities(this->ljsim->get_velocities());

    if(this->graph_widget != nullptr) {
        qDebug() << "Resetting graphs";
        this->graph_widget->reset_graphs();
    }
    this->graph_widget->set_params(params->get_param<double>("kT"), params->get_param<double>("mass"), this->ljsim->get_velocities().size());
}

/**
 * @brief Build the control interface
 */
void MainWindow::build_control_interface() {
    // Initialization and timer for sequence
    this->control_interface = new QWidget();

    // horizontal buttons
    QBoxLayout *horizontalBox = new QHBoxLayout;
    this->control_interface->setLayout(horizontalBox);

    // add toggle axes button to horizontal row
    this->button_toggle_axis = new QPushButton;
    this->button_toggle_axis->setIcon(QIcon(":/assets/icon/axes_32.png"));
    this->button_toggle_axis->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
    this->button_toggle_axis->setToolTip(QString("Toggle coordination axes in viewport."));
    horizontalBox->addWidget(this->button_toggle_axis);
    connect(this->button_toggle_axis, SIGNAL(clicked()), this, SLOT(toggle_world_axes()));

    // add toggle axes button to horizontal row
    this->button_rotate_model = new QPushButton;
    this->button_rotate_model->setIcon(QIcon(":/assets/icon/rotation_gray_32.png"));
    this->button_rotate_model->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
    this->button_rotate_model->setToolTip(QString("Toggle rotation of model around z-axis."));
    horizontalBox->addWidget(this->button_rotate_model);
    connect(this->button_rotate_model, SIGNAL(clicked()), this, SLOT(toggle_rotation()));

    // add button to pause the simulation
    this->button_simulation_pause = new QPushButton("Pause/Unpause");
    horizontalBox->addWidget(this->button_simulation_pause);

    // add button to pause the simulation
    this->button_relative_coloring = new QPushButton("Change coloring");
    horizontalBox->addWidget(this->button_relative_coloring);
    connect(this->button_relative_coloring, SIGNAL(released()), this->anaglyph_widget, SLOT(toggle_relative_coloring()));

    QFrame *frame = new QFrame;
    frame->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Minimum);
    horizontalBox->addWidget(frame);
}

/**
 * Create the dropdown menu and associate all actions
 */
void MainWindow::create_dropdown_menu() {
    // dropdown menu bar
    QMenuBar *menuBar = new QMenuBar;

    // add drop-down menus
    QMenu *menu_file = menuBar->addMenu(tr("&File"));
    QMenu *menu_view = menuBar->addMenu(tr("&View"));
    QMenu *menu_help = menuBar->addMenu(tr("&Help"));

    // actions for file menu
    QAction *action_new = new QAction(menu_file);
    QAction *action_quit = new QAction(menu_file);

    // actions for projection menu
    QMenu *menu_camera = new QMenu(tr("Camera"));

    // camera alignment
    QMenu *menu_camera_align = new QMenu(tr("Align"));
    QAction *action_camera_default = new QAction(menu_camera_align);
    QAction *action_camera_top = new QAction(menu_camera_align);
    QAction *action_camera_bottom = new QAction(menu_camera_align);
    QAction *action_camera_left = new QAction(menu_camera_align);
    QAction *action_camera_right = new QAction(menu_camera_align);
    QAction *action_camera_front = new QAction(menu_camera_align);
    QAction *action_camera_back = new QAction(menu_camera_align);

    // camera modes
    QMenu *menu_camera_mode = new QMenu(tr("Mode"));
    QAction *action_camera_perspective = new QAction(menu_camera_mode);
    QAction *action_camera_orthographic = new QAction(menu_camera_mode);

    // camera projections
    QMenu *menu_projection = new QMenu(tr("Projection"));
    QAction *action_projection_two_dimensional = new QAction(menu_projection);
    QAction *action_projection_anaglyph_red_cyan = new QAction(menu_projection);

    // interlaced projections
    QMenu *menu_projection_interlaced = new QMenu(tr("Interlaced"));
    QAction *action_projection_interlaced_rows_lr = new QAction(menu_projection_interlaced);
    QAction *action_projection_interlaced_rows_rl = new QAction(menu_projection_interlaced);
    QAction *action_projection_interlaced_columns_lr = new QAction(menu_projection_interlaced);
    QAction *action_projection_interlaced_columns_rl = new QAction(menu_projection_interlaced);
    QAction *action_projection_interlaced_checkerboard_lr = new QAction(menu_projection_interlaced);
    QAction *action_projection_interlaced_checkerboard_rl = new QAction(menu_projection_interlaced);

    // actions for help menu
    QAction *action_about = new QAction(menu_help);

    // debug log
    QAction *action_debug_log = new QAction(menu_help);
    action_debug_log->setText(tr("Debug Log"));
    action_debug_log ->setShortcut(Qt::Key_F2);
    menu_help->addAction(action_debug_log);
    connect(action_debug_log, &QAction::triggered, this, &MainWindow::slot_debug_log);

    // create actions for file menu
    action_new->setText(tr("New"));
    action_new->setShortcuts(QKeySequence::New);
    action_new->setShortcut(Qt::CTRL | Qt::Key_N);
    action_quit->setText(tr("Quit"));
    action_quit->setShortcuts(QKeySequence::Quit);
    action_quit->setShortcut(Qt::CTRL | Qt::Key_Q);

    // create actions for projection menu
    action_camera_default->setText(tr("Default"));
    action_camera_default->setData(QVariant((int)CameraAlignment::DEFAULT));
    action_camera_default->setShortcut(Qt::Key_0);
    action_camera_top->setText(tr("Top"));
    action_camera_top->setData(QVariant((int)CameraAlignment::TOP));
    action_camera_top->setShortcut(Qt::Key_7);
    action_camera_bottom->setText(tr("Bottom"));
    action_camera_bottom->setData(QVariant((int)CameraAlignment::BOTTOM));
    action_camera_bottom->setShortcut(Qt::CTRL | Qt::Key_7);
    action_camera_left->setText(tr("Left"));
    action_camera_left->setData(QVariant((int)CameraAlignment::LEFT));
    action_camera_left->setShortcut(Qt::Key_3);
    action_camera_right->setText(tr("Right"));
    action_camera_right->setData(QVariant((int)CameraAlignment::RIGHT));
    action_camera_right->setShortcut(Qt::CTRL | Qt::Key_3);
    action_camera_front->setText(tr("Front"));
    action_camera_front->setData(QVariant((int)CameraAlignment::FRONT));
    action_camera_front->setShortcut(Qt::Key_1);
    action_camera_back->setText(tr("Back"));
    action_camera_back->setData(QVariant((int)CameraAlignment::BACK));
    action_camera_back->setShortcut(Qt::CTRL | Qt::Key_1);
    action_camera_perspective->setText(tr("Perspective"));
    action_camera_perspective->setData(QVariant((int)CameraMode::PERSPECTIVE));
    action_camera_perspective->setShortcut(Qt::Key_5);
    action_camera_orthographic->setText(tr("Orthographic"));
    action_camera_orthographic->setData(QVariant((int)CameraMode::ORTHOGRAPHIC));
    action_camera_orthographic->setShortcut(Qt::CTRL | Qt::Key_5);
    action_projection_two_dimensional->setText(tr("Two-dimensional"));
    action_projection_anaglyph_red_cyan->setText(tr("Anaglyph (red/cyan)"));
    action_projection_interlaced_rows_lr->setText(tr("Interlaced rows (left first)"));
    action_projection_interlaced_rows_rl->setText(tr("Interlaced rows (right first)"));
    action_projection_interlaced_columns_lr->setText(tr("Interlaced columns (left first)"));
    action_projection_interlaced_columns_rl->setText(tr("Interlaced columns (right first)"));
    action_projection_interlaced_checkerboard_lr->setText(tr("Checkerboard (left first)"));
    action_projection_interlaced_checkerboard_rl->setText(tr("Checkerboard (right first)"));
    menu_projection_interlaced->setIcon(QIcon(":/assets/icon/interlaced_rows_lr_32.png"));
    action_projection_two_dimensional->setIcon(QIcon(":/assets/icon/two_dimensional_32.png"));
    action_projection_anaglyph_red_cyan->setIcon(QIcon(":/assets/icon/anaglyph_red_cyan_32.png"));
    action_projection_interlaced_rows_lr->setIcon(QIcon(":/assets/icon/interlaced_rows_lr_32.png"));
    action_projection_interlaced_rows_rl->setIcon(QIcon(":/assets/icon/interlaced_rows_rl_32.png"));
    action_projection_interlaced_columns_lr->setIcon(QIcon(":/assets/icon/interlaced_columns_lr_32.png"));
    action_projection_interlaced_columns_rl->setIcon(QIcon(":/assets/icon/interlaced_columns_rl_32.png"));
    action_projection_interlaced_checkerboard_lr->setIcon(QIcon(":/assets/icon/interlaced_checkerboard_lr_32.png"));
    action_projection_interlaced_checkerboard_rl->setIcon(QIcon(":/assets/icon/interlaced_checkerboard_rl_32.png"));

    // create actions for about menu
    action_about->setText(tr("About"));
    action_about->setIcon(QIcon(":/assets/icon/info.png"));

    // add actions to file menu
    menu_file->addAction(action_new);
    menu_file->addAction(action_quit);

    // add actions to projection menu
    menu_view->addMenu(menu_projection);
    menu_view->addMenu(menu_camera);
    menu_camera->addMenu(menu_camera_align);
    menu_camera_align->addAction(action_camera_default);
    menu_camera_align->addAction(action_camera_top);
    menu_camera_align->addAction(action_camera_bottom);
    menu_camera_align->addAction(action_camera_left);
    menu_camera_align->addAction(action_camera_right);
    menu_camera_align->addAction(action_camera_front);
    menu_camera_align->addAction(action_camera_back);
    menu_camera->addMenu(menu_camera_mode);
    menu_camera_mode->addAction(action_camera_perspective);
    menu_camera_mode->addAction(action_camera_orthographic);
    menu_projection->addAction(action_projection_two_dimensional);
    menu_projection->addAction(action_projection_anaglyph_red_cyan);
    menu_projection->addMenu(menu_projection_interlaced);
    menu_projection_interlaced->addAction(action_projection_interlaced_rows_lr);
    menu_projection_interlaced->addAction(action_projection_interlaced_rows_rl);
    menu_projection_interlaced->addAction(action_projection_interlaced_columns_lr);
    menu_projection_interlaced->addAction(action_projection_interlaced_columns_rl);
    menu_projection_interlaced->addAction(action_projection_interlaced_checkerboard_lr);
    menu_projection_interlaced->addAction(action_projection_interlaced_checkerboard_rl);

    // add actions to help menu
    menu_help->addAction(action_about);

    // connect actions file menu
    connect(action_new, &QAction::triggered, this, &MainWindow::new_simulation);
    connect(action_quit, &QAction::triggered, this, &MainWindow::exit);

    // connect actions projection menu (note; [this]{} is the idiomatic way by providing a functor - "this is the way")
    connect(action_projection_two_dimensional, &QAction::triggered, this, [this]{ MainWindow::set_stereo("no_stereo_flat"); });
    connect(action_projection_anaglyph_red_cyan, &QAction::triggered, this, [this]{ MainWindow::set_stereo("stereo_anaglyph_red_cyan"); });
    connect(action_projection_interlaced_rows_lr, &QAction::triggered, this, [this]{ MainWindow::set_stereo("stereo_interlaced_rows_lr"); });
    connect(action_projection_interlaced_rows_rl, &QAction::triggered, this, [this]{ MainWindow::set_stereo("stereo_interlaced_rows_rl"); });
    connect(action_projection_interlaced_columns_lr, &QAction::triggered, this, [this]{ MainWindow::set_stereo("stereo_interlaced_columns_lr"); });
    connect(action_projection_interlaced_columns_rl, &QAction::triggered, this, [this]{ MainWindow::set_stereo("stereo_interlaced_columns_rl"); });
    connect(action_projection_interlaced_checkerboard_lr, &QAction::triggered, this, [this]{ MainWindow::set_stereo("stereo_interlaced_checkerboard_lr"); });
    connect(action_projection_interlaced_checkerboard_rl, &QAction::triggered, this, [this]{ MainWindow::set_stereo("stereo_interlaced_checkerboard_rl"); });

    // connect actions camera menu
    // connect(menu_camera_align, SIGNAL(triggered(QAction*)), this->interface_window, SLOT(set_camera_align(QAction*)));
    // connect(menu_camera_mode, SIGNAL(triggered(QAction*)), this->interface_window, SLOT(set_camera_mode(QAction*)));

    // connect actions about menu
    connect(action_about, &QAction::triggered, this, &MainWindow::about);

    // build interface
    setMenuBar(menuBar);
}

/**
 * @brief Sets up the main interface layout for the application.
 *
 * This function initializes and arranges the widgets within the main window.
 * It creates a central container with a horizontal layout, which consists of:
 * - A left vertical container holding the AnaglyphWidget (maximized horizontally)
 *   and the control interface (minimized as much as possible).
 * - A GraphWidget on the right.
 */
void MainWindow::set_widgets() {
    // Create main container widget and layout
    QWidget *mainWidget = new QWidget;
    QHBoxLayout *mainLayout = new QHBoxLayout();
    mainWidget->setLayout(mainLayout);
    this->setCentralWidget(mainWidget);

    // Create vertical container for anaglyph and control interface
    QWidget *leftContainer = new QWidget;
    QVBoxLayout *leftLayout = new QVBoxLayout();
    leftContainer->setLayout(leftLayout);

    // Add the anaglyph widget
    this->anaglyph_widget = new AnaglyphWidget(this);
    leftLayout->addWidget(this->anaglyph_widget);
    this->anaglyph_widget->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

    // Build and add control interface
    this->build_control_interface();
    leftLayout->addWidget(this->control_interface);
    this->control_interface->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Minimum);

    // Add the left container to the main layout
    mainLayout->addWidget(leftContainer);

    // Create and add the graph widget
    this->graph_widget = new GraphWidget();
    mainLayout->addWidget(this->graph_widget);
}

/**
 * @brief Construct a new simulation
 */
 void MainWindow::new_simulation() {
    if(this->thread != nullptr) {
        this->thread->stop();
    }
    auto params = std::make_shared<LennardJonesParameters>();
    DialogNewSimulation dialog(params, this);
    int ret = dialog.exec();

    if(ret == QDialog::Accepted) {
        if(this->thread == nullptr) {
            qDebug() << "Building new thread object.";
            // construct LJ simulation
            this->build_simulation(params);

            // spawn thread and run simulation
            this->thread = new ThreadIntegrate(this->ljsim);
            this->thread->toggle_pause();
            connect(this->thread, SIGNAL(signal_integration_step()), this, SLOT(integration_step()));
            connect(this->thread, SIGNAL(signal_ekin(double, double)), this->graph_widget, SLOT(add_data_item_kinetic_energy(double, double)));
            connect(this->thread, SIGNAL(signal_epot(double, double)), this->graph_widget, SLOT(add_data_item_potential_energy(double, double)));
            connect(this->thread, SIGNAL(signal_etot(double, double)), this->graph_widget, SLOT(add_data_item_total_energy(double, double)));
            connect(this->thread, SIGNAL(signal_velocities()), this, SLOT(slot_transmit_velocities()));
            connect(this->graph_widget, SIGNAL(signal_simulation_continue()), this->thread, SLOT(slot_simulation_unlock()));
            connect(this->button_simulation_pause, SIGNAL(released()), this->thread, SLOT(toggle_pause()));
            this->thread->start();
        } else {
            qDebug() << "Constructing new simulation";

            // construct LJ simulation
            this->build_simulation(params);
            this->thread->set_simulation(this->ljsim);
        }
    }
}