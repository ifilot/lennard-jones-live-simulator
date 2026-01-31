cmake_minimum_required(VERSION 3.16)

# Set the project name
project(ljsim LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# get Git HASH
execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
add_compile_definitions(GIT_HASH="${GIT_HASH}")

# Find Qt5 packages
find_package(Qt5 REQUIRED COMPONENTS Widgets Charts)

# Use GLM
find_package(glm REQUIRED)

# Use Eigen
find_package(Eigen3 REQUIRED NO_MODULE)

add_executable(ljsim
    src/main.cpp
    src/gui/anaglyph_widget.cpp
    src/gui/mainwindow.cpp
    src/gui/shader_program.cpp
    src/gui/shader_program_manager.cpp
    src/gui/structure_renderer.cpp
    src/gui/user_action.cpp
    src/gui/scene.cpp
    src/gui/logwindow.cpp
    src/data/model.cpp
    src/data/model_loader.cpp
    src/data/structure.cpp
    src/simulator/dialognewsimulation.cpp
    src/simulator/graphwidget.cpp
    src/simulator/lennardjonesparameters.cpp
    src/simulator/lennardjonessimulation.cpp
    src/simulator/threadintegrate.cpp
    src/ljsimapp.cpp
    resources.qrc
)

target_link_libraries(ljsim PRIVATE 
    Qt5::Core 
    Qt5::Widgets 
    Qt5::Charts
    Eigen3::Eigen
    glm::glm
)

set_target_properties(ljsim PROPERTIES
    WIN32_EXECUTABLE ON
    MACOSX_BUNDLE ON
)

if (WIN32)
    target_sources(ljsim PRIVATE ljsim.rc)
endif()