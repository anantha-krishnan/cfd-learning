#-------------------------------------------------
#
# Project created by QtCreator 2016-03-09T15:34:31
#
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = couette_flow
CONFIG   += console
CONFIG   -= app_bundle


TEMPLATE = app


SOURCES += main.cpp \
    couette_flow.cpp \
    simple.cpp \
    numerical_methods.cpp \
    structured_domain.cpp \
    constants.cpp \
    simplepatankar.cpp

HEADERS += \
    couette_flow.h \
    simple.h \
    numerical_methods.h \
    structured_domain.h \
    constants.h \
    simplepatankar.h
