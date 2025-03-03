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

#include "user_action.h"

/**
 * @brief      Constructs a new instance.
 */
UserAction::UserAction(const std::shared_ptr<Scene>& _scene):
    scene(_scene) {
}

/**
 * @brief      Updates the given cursor position.
 *
 * @param[in]  cursor_position  The cursor position
 */
void UserAction::update(QPoint _cursor_position) {
    this->cursor_position_now = _cursor_position;
}

/**
 * @brief      Handle a translation action
 */
void UserAction::handle_action_movement() {
}

/**
 * @brief      Handle a rotation action
 */
void UserAction::handle_action_rotation() {
}

/**
 * @brief      Handle left mouse click action
 */
void UserAction::handle_left_mouse_click() {
}

/**
 * @brief      Handle key stroke
 *
 * @param[in]  key        The key
 * @param[in]  modifiers  The modifiers
 */
void UserAction::handle_key(int key, Qt::KeyboardModifiers modifiers) {
}

/**
 * @brief      Sets the camera alignment.
 *
 * @param[in]  direction  The direction
 */
void UserAction::set_camera_alignment(int direction) {
    QVector3D dirvec;

    switch((CameraAlignment)direction) {
        case CameraAlignment::DEFAULT:
            this->scene->rotation_matrix.setToIdentity();
            this->scene->rotation_matrix.rotate(20.0, QVector3D(1,0,0));
            this->scene->rotation_matrix.rotate(30.0, QVector3D(0,0,1));
            emit(request_update());
        return;
        case CameraAlignment::TOP:
            dirvec = QVector3D(0.0f, 0.0f, 1.0f);
        break;
        case CameraAlignment::BOTTOM:
            dirvec = QVector3D(0.0f, 0.0f, -1.0f);
        break;
        case CameraAlignment::LEFT:
            dirvec = QVector3D(-1.0f, 0.0f, 0.0f);
        break;
        case CameraAlignment::RIGHT:
            dirvec = QVector3D(1.0f, 0.0f, 0.0f);
        break;
        case CameraAlignment::FRONT:
            dirvec = QVector3D(0.0f, 1.0f, 0.0f);
        break;
        case CameraAlignment::BACK:
            dirvec = QVector3D(0.0f, -1.0f, 0.0f);
        break;
    }

    QVector3D axis;
    float angle;

    // avoid gimball locking
    if (fabs(dirvec[1]) > .999) {
        if(dirvec[1] < 0.0) {
            axis = QVector3D(0.0, 0.0, 1.0);
            angle = -M_PI;
        } else {
            axis = QVector3D(0.0, 0.0, 1.0);
            angle = 0.0;
        }
    } else {
        axis = QVector3D::crossProduct(QVector3D(0.0, 1.0, 0.0), dirvec);
        angle = std::acos(dirvec[1]);
    }

    this->scene->rotation_matrix.setToIdentity();
    this->scene->rotation_matrix.rotate(qRadiansToDegrees(angle), axis);
    emit(signal_message_statusbar("Change camera alignment"));
    emit(request_update());
}

/**
 * @brief      Sets the camera mode.
 *
 * @param[in]  mode  The mode
 */
void UserAction::set_camera_mode(int mode) {
    float w = (float)this->scene->canvas_width;
    float h = (float)this->scene->canvas_height;

    switch((CameraMode)mode) {
        case CameraMode::PERSPECTIVE:
            this->scene->camera_mode = CameraMode::PERSPECTIVE;
            this->scene->projection.setToIdentity();
            this->scene->projection.perspective(45.0f, w / h, 0.01f, 1000.0f);
            emit(signal_message_statusbar("Set camera to perspective"));
        break;
        case CameraMode::ORTHOGRAPHIC:
            this->scene->camera_mode = CameraMode::ORTHOGRAPHIC;
            float w = (float)this->scene->canvas_width;
            float h = (float)this->scene->canvas_height;
            float ratio = w/h;
            float zoom = -this->scene->camera_position[1];
            this->scene->projection.setToIdentity();
            this->scene->projection.ortho(-zoom/2.0f, zoom/2.0f, -zoom / ratio /2.0f, zoom / ratio / 2.0f, 0.01f, 1000.0f);
            emit(signal_message_statusbar("Set camera to orthographic"));
        break;
    }

    emit(request_update());
}

/**
 * @brief      Sets the cursor position.
 */
void UserAction::set_cursor_position() {
    this->cursor_position_start = QCursor::pos();
    qDebug() << this->cursor_position_start;
}