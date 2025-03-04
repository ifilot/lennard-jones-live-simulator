# Project name
TEMPLATE = app
TARGET = ljsim

# Qt Modules
QT += core widgets charts

# C++ Standard
CONFIG += c++17

# Include directories
INCLUDEPATH += \
    src/vendor/eigen-3.4.0 \
    src/vendor/glm-1.0.1

# Source files
SOURCES += \
    src/main.cpp \
    src/gui/anaglyph_widget.cpp \
    src/gui/mainwindow.cpp \
    src/gui/shader_program.cpp \
    src/gui/shader_program_manager.cpp \
    src/gui/structure_renderer.cpp \
    src/gui/user_action.cpp \
    src/gui/scene.cpp \
    src/gui/logwindow.cpp \
    src/data/model.cpp \
    src/data/model_loader.cpp \
    src/data/structure.cpp \
    src/simulator/dialognewsimulation.cpp \
    src/simulator/graphwidget.cpp \
    src/simulator/lennardjonesparameters.cpp \
    src/simulator/lennardjonessimulation.cpp \
    src/simulator/threadintegrate.cpp \
    src/ljsimapp.cpp

# Header files (optional, but good practice)
HEADERS += \
    src/gui/anaglyph_widget.h \
    src/gui/mainwindow.h \
    src/gui/shader_program.h \
    src/gui/shader_program_manager.h \
    src/gui/structure_renderer.h \
    src/gui/user_action.h \
    src/gui/scene.h \
    src/gui/logwindow.h \
    src/data/model.h \
    src/data/model_loader.h \
    src/data/structure.h \
    src/simulator/dialognewsimulation.h \
    src/simulator/graphwidget.h \
    src/simulator/lennardjonesparameters.h \
    src/simulator/lennardjonessimulation.h \
    src/simulator/threadintegrate.h

# Resource files
RESOURCES += resources.qrc

# Enable Qt's automatic MOC, UIC, and RCC handling
CONFIG += qt autogen openmp c++17

# Platform-specific settings
win32 {
    CONFIG += windows
    QMAKE_CXXFLAGS += -fopenmp
    QMAKE_LFLAGS += -fopenmp
}
macx {
    CONFIG += app_bundle
}

