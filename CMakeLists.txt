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

# Find Qt6 packages
find_package(Qt6 REQUIRED COMPONENTS Widgets Charts)

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
    src/data/atom_settings.cpp
    src/data/atom.cpp
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

include_directories("src/vendor/eigen-3.4.0"
                    "src/vendor/glm-1.0.1")

target_link_libraries(ljsim PRIVATE Qt6::Core Qt6::Widgets Qt6::Charts)

set_target_properties(ljsim PROPERTIES
    WIN32_EXECUTABLE ON
    MACOSX_BUNDLE ON
)