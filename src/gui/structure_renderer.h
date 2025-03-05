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

#include <QOpenGLFunctions>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QDebug>
#include <QMatrix4x4>
#include <QtMath>
#include <QTemporaryDir>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/norm.hpp>

#include <vector>

#include "../data/model_loader.h"
#include "../data/structure.h"
#include "shader_program_manager.h"
#include "user_action.h"

class StructureRenderer {
private:
    // sphere facets
    std::vector<glm::vec3> sphere_vertices;
    std::vector<glm::vec3> sphere_normals;
    std::vector<unsigned int> sphere_indices;

    // cylinder facets
    std::vector<glm::vec3> cylinder_vertices;
    std::vector<glm::vec3> cylinder_normals;
    std::vector<unsigned int> cylinder_indices;

    // vao and vbo for rendering
    QOpenGLVertexArrayObject vao_sphere;
    QOpenGLBuffer vbo_sphere[3];

    QOpenGLVertexArrayObject vao_cylinder;
    QOpenGLBuffer vbo_cylinder[3];

    QOpenGLVertexArrayObject vao_unitcell;
    QOpenGLBuffer vbo_unitcell[2];

    QOpenGLVertexArrayObject vao_line;
    QOpenGLBuffer vbo_line[2];

    QOpenGLVertexArrayObject vao_plane;
    QOpenGLBuffer vbo_plane[2];

    // couple of pointers to important matrices
    std::shared_ptr<Scene> scene;

    std::shared_ptr<ShaderProgramManager> shader_manager;
    std::shared_ptr<UserAction> user_action;

    // models
    std::shared_ptr<Model> axis_model;

    bool flag_draw_unitcell = true;     // whether to draw the unitcell

    // used for velocity color map
    const QVector<QVector3D> color_scheme = {
        QVector3D(0x05 / 255.0f, 0x30 / 255.0f, 0x61 / 255.0f),
        QVector3D(0x21 / 255.0f, 0x66 / 255.0f, 0xac / 255.0f),
        QVector3D(0x43 / 255.0f, 0x93 / 255.0f, 0xc3 / 255.0f),
        QVector3D(0x92 / 255.0f, 0xc5 / 255.0f, 0xde / 255.0f),
        QVector3D(0xd1 / 255.0f, 0xe5 / 255.0f, 0xf0 / 255.0f),
        QVector3D(0xf7 / 255.0f, 0xf7 / 255.0f, 0xf7 / 255.0f),
        QVector3D(0xfd / 255.0f, 0xdb / 255.0f, 0xc7 / 255.0f),
        QVector3D(0xf4 / 255.0f, 0xa5 / 255.0f, 0x82 / 255.0f),
        QVector3D(0xd6 / 255.0f, 0x60 / 255.0f, 0x4d / 255.0f),
        QVector3D(0xb2 / 255.0f, 0x18 / 255.0f, 0x2b / 255.0f),
        QVector3D(0x67 / 255.0f, 0x00 / 255.0f, 0x1f / 255.0f)
    };

    bool flag_relative_coloring = false;

public:
    /**
     * @brief      Constructs a new instance.
     *
     * @param[in]  _scene           The scene
     * @param[in]  _shader_manager  The shader manager
     * @param[in]  _user_action     The user action
     */
    StructureRenderer(const std::shared_ptr<Scene>& _scene,
                      const std::shared_ptr<ShaderProgramManager>& _shader_manager,
                      const std::shared_ptr<UserAction>& _user_action);

    /**
     * @brief      Draw the structure
     *
     * @param[in]  structure     The structure
     * @param      model_shader  The model shader
     */
    void draw(const Structure *structure);

    /**
     * @brief      Draw the structure
     *
     * @param[in]  structure     The structure
     */
    void draw_silhouette(const Structure *structure);

    /**
     * @brief      Draws coordinate axes.
     */
    void draw_coordinate_axes();

    /**
     * @brief      Disables the drawing of the unitcell
     */
    inline void disable_draw_unitcell() {
        this->flag_draw_unitcell = false;
    }

private:
    /**
     * @brief      Draws atoms in the regular unit cell.
     *
     * @param[in]  structure  The structure
     */
    void draw_atoms(const Structure* structure);

    /**
     * @brief      Draws the unitcell.
     *
     * @param[in]  structure  The structure
     */
    void draw_unitcell(const Structure* structure);

    /**
     * @brief      Draw movement lines
     *
     * @param[in]  structure  The structure
     */
    void draw_movement_lines(const Structure* structure);

    /**
     * @brief      Draw movement plane
     *
     * @param[in]  structure  The structure
     */
    void draw_movement_plane(const Structure* structure);

    /**
     * @brief      Generate coordinates of a sphere
     *
     * @param[in]  tesselation_level  The tesselation level
     */
    void generate_sphere_coordinates(unsigned int tesselation_level);

    /**
     * @brief      Generate coordinates for a default cylinder (radius 1, height 1)
     *
     * @param[in]  stack_count  The stack count
     * @param[in]  slice_count  The slice count
     */
    void generate_cylinder_coordinates(unsigned int stack_count, unsigned int slice_count);

    /**
     * @brief      Generate the coordinates of the unitcell
     */
    void generate_coordinates_unitcell(const MatrixUnitcell& unitcell);

    /**
     * @brief      Update the unitcell vertices in the unit cell
     *
     * @param[in]  unitcell  The unitcell
     */
    void set_unitcell_vertices(const MatrixUnitcell& unitcell);

    /**
     * @brief      Load all data to a vertex array object
     */
    void load_sphere_to_vao();

    /**
     * @brief      Load all data to a vertex array object
     */
    void load_cylinder_to_vao();

    /**
     * @brief      Load simple line data to vertex array object
     */
    void load_line_to_vao();

    /**
     * @brief      Load simple plane data to vertex array object
     */
    void load_plane_to_vao();

    /**
     * @brief      Loads an arrow model.
     */
    void load_arrow_model();

    /**
     * @brief      Darken color
     *
     * @param[in]  color   The color
     * @param[in]  amount  The amount
     *
     * @return     The 3D vector.
     */
    QVector3D darken(const QVector3D& color, float amount) const;

    /**
     * @brief      Lighten color
     *
     * @param[in]  color   The color
     * @param[in]  amount  The amount
     *
     * @return     The 3D vector.
     */
    QVector3D lighten(const QVector3D& color, float amount) const;

    /**
     * @brief      Mix colors
     *
     * @param[in]  color1  First color
     * @param[in]  color2  Second color
     * @param[in]  amount  The amount
     *
     * @return     The 3D vector.
     */
    QVector3D mix(const QVector3D& color1, const QVector3D& color2, float amount) const;

    /**
     * @brief Get color from scale
     * @return QVector3D color
     */
    QVector3D get_color_from_scale(double low, double high, double val);

};
