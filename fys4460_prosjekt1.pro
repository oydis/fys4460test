TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    Lattice.cpp \
    Atom.cpp \
    Integrator.cpp \
    Cell.cpp \
    Cell_container.cpp

HEADERS += \
    Lattice.h \
    Atom.h \
    Integrator.h \
    Cell.h \
    Cell_container.h

release {
    DEFINES += ARMA_NO_DEBUG
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE += -O3
}
