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

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QMouseEvent>
#include <QCoreApplication>
#include <QString>
#include <QGuiApplication>
#include <QScreen>
#include <QSysInfo>
#include <QDebug>
#include <QTimer>
#include <QMenu>
#include <QtGlobal>

#include <QtCore/qmath.h>
#include <QtCore/qvariant.h>

#include <sstream>
#include <fstream>
#include <math.h>
#include <string>

#include "shader_program_manager.h"
#include "shader_program_types.h"
#include "structure_renderer.h"
#include "user_action.h"
#include "scene.h"

QT_FORWARD_DECLARE_CLASS(QOpenGLShaderProgram)

enum FrameBuffer {
    STRUCTURE_NORMAL,
    STRUCTURE_LEFT,
    STRUCTURE_RIGHT,
    ANAGLYPH_LEFT,
    ANAGLYPH_RIGHT,
    COORDINATE_AXES,

    NR_FRAMEBUFFERS
};

class AnaglyphWidget : public QOpenGLWidget, protected QOpenGLFunctions {
    Q_OBJECT

private:
    QPoint m_lastPos;

    QString root_path;

    QPoint top_left;

    // default background color
    const QVector4D background = {21.0/255.0, 
                                  21.0/255.0, 
                                  21.0/255.0, 
                                  1.0};

    unsigned int framebuffers[FrameBuffer::NR_FRAMEBUFFERS];
    unsigned int texture_color_buffers[FrameBuffer::NR_FRAMEBUFFERS];
    unsigned int rbo[FrameBuffer::NR_FRAMEBUFFERS];

    QOpenGLVertexArrayObject quad_vao;
    QOpenGLVertexArrayObject quad_vao_small;
    QOpenGLBuffer quad_vbo;

    // used for arcball rotation
    bool arcball_rotation_flag = false;             // whether arcball rotation is active

    QOpenGLVertexArrayObject vao;
    QOpenGLBuffer vbo[4];

    bool flag_axis_enabled = true;                  // whether to draw coordinate axes
    bool flag_draw_unitcell = true;                 // whether to draw the unitcell

    // stereographic projections
    bool flag_stereographic_projection = false;     // whether stereographic rendering is used
    QString stereographic_type_name = "NONE";

    std::shared_ptr<Scene> scene;
    std::unique_ptr<StructureRenderer> structure_renderer;
    std::shared_ptr<ShaderProgramManager> shader_manager;
    std::shared_ptr<Structure> structure;
    std::shared_ptr<UserAction> user_action;        // object that stores current action of the user on a structure

public:
    AnaglyphWidget(QWidget *parent = 0);

    /**
     * @brief Rotate the scene
     * @param angle
     */
     void rotate_scene(float angle);

    /**
     * @brief      Paint the models in the models vector to the screen
     */
    void draw_structure();

    /**
     * @brief      Sets a new structure without modifying camera settings.
     *
     * This function updates the internal structure reference with a new one.
     * It does not alter camera settings, ensuring the view remains consistent.
     *
     * @param[in]  _structure  A shared pointer to the new structure.
     */
    void set_structure(const std::shared_ptr<Structure>& _structure);

    /**
     * @brief      Retrieves the current structure.
     *
     * This function returns a constant reference to the currently set structure.
     *
     * @return     A constant reference to the structure.
     */
    const auto& get_structure() const {
        return this->structure;
    }

    /**
     * @brief      Handles window movement events.
     *
     * This function is triggered when the window is moved, allowing the 
     * application to respond accordingly (e.g., updating UI elements or
     * recalculating positions).
     */
    void window_move_event();

    /**
     * @brief      Sets the stereo mode.
     *
     * Configures the rendering mode for stereoscopic 3D display based on the 
     * provided stereo name.
     *
     * @param[in]  stereo_name  The name of the stereo mode to apply.
     */
    void set_stereo(QString stereo_name);

    /**
     * @brief      Provides a minimum size hint for the widget.
     *
     * This function returns the smallest recommended size for the widget,
     * ensuring proper layout and usability.
     *
     * @return     The minimum recommended widget size as a QSize object.
     */
    QSize minimumSizeHint() const Q_DECL_OVERRIDE;

    /**
     * @brief      Provides a preferred size hint for the widget.
     *
     * Returns the suggested size for the widget, helping layout managers 
     * determine an optimal default size.
     *
     * @return     The preferred widget size as a QSize object.
     */
    QSize sizeHint() const Q_DECL_OVERRIDE;

    /**
    * @brief      Destructor for AnaglyphWidget.
    *
    * Cleans up resources associated with the widget before destruction.
    */
    ~AnaglyphWidget();

    /**
     * @brief      Enables or disables the axis display.
     *
     * This function sets a flag to enable or disable the axis rendering
     * in the visualization.
     *
     * @param[in]  value  Boolean flag to enable (true) or disable (false) the axis.
     */
    inline void set_axis_enabled(bool value) {
        this->flag_axis_enabled = value;
    }

    /**
     * @brief      Checks if the axis display is enabled.
     *
     * This function returns whether the axis rendering is currently enabled.
     *
     * @return     True if axis rendering is enabled, otherwise false.
     */
    inline bool get_axis_enabled() const {
        return this->flag_axis_enabled;
    }


    /**
     * @brief      Gets the user action object.
     *
     * @return     The user action.
     */
    const auto& get_user_action() {
        return this->user_action;
    }

    /**
     * @brief      Disables the drawing of the unitcell
     */
    inline void disable_draw_unitcell() {
        this->flag_draw_unitcell = false;
        if(this->structure_renderer) {
            this->structure_renderer->disable_draw_unitcell();
        }
    }

public slots:
    /**
     * @brief      Clean the anaglyph class
     */
    void cleanup();

protected:
    /**
     * @brief      Initialize OpenGL environment
     */
    void initializeGL() Q_DECL_OVERRIDE;

    /**
     * @brief      Render scene
     */
    void paintGL() Q_DECL_OVERRIDE;

    /**
     * @brief      Resize window
     *
     * @param[in]  width   screen width
     * @param[in]  height  screen height
     */
    void resizeGL(int width, int height) Q_DECL_OVERRIDE;

    /**
     * @brief      Parse mouse press event
     *
     * @param      event  The event
     */
    void mousePressEvent(QMouseEvent *event) Q_DECL_OVERRIDE;

    /**
     * @brief      Parse mouse release event
     *
     * @param      event  The event
     */
    void mouseReleaseEvent(QMouseEvent *event) Q_DECL_OVERRIDE;

    /**
     * @brief      Parse mouse move event
     *
     * @param      event  The event
     */
    void mouseMoveEvent(QMouseEvent *event) Q_DECL_OVERRIDE;

    /**
     * @brief      Parse mouse wheel event
     *
     * @param      event  The event
     */
    void wheelEvent(QWheelEvent *event) Q_DECL_OVERRIDE;

    /**
     * @brief      Calculate the arcball vector for mouse rotation
     *
     * @param[in]  x, y  The mouse position
     * @param[out] P  The arcball vector
     */
    QVector3D get_arcball_vector(int x, int y);

    /**
     * @brief      Set the arcball vector rotation (angle and vector) to model and updates
     *
     * @param      arcball_angle, arcball_vector
     */
    void set_arcball_rotation(float arcball_angle, const QVector4D& arcball_vector);

private:
    /**
     * @brief      Build a canvas used for stereographic rendering
     */
    void build_framebuffers();

    /**
     * @brief      Load OpenGL shaders
     */
    void load_shaders();

    /**
     * @brief      Reset rotation matrices
     */
    void reset_matrices();

    /**
     * @brief      Regular draw call
     */
    void paint_regular();

    /**
     * @brief      Stereographic draw call
     */
    void paint_stereographic();

private slots:
    /**
     * @brief      Open menu for atom
     *
     * @param[in]  pos   The position
     */
    void custom_menu_requested(QPoint pos);

    /**
     * @brief      Update the scene
     */
    void call_update();

    /**
     * @brief      Transmit a message to the user
     */
    void transmit_message(const QString& text);

signals:
    /**
     * @brief      Send signal that opengl engine is ready
     */
    void opengl_ready();

    /**
     * @brief      Message received
     */
    void signal_interaction_message(const QString& text);

    /**
     * @brief      Message received
     */
    void signal_selection_message(const QString& text);
};
