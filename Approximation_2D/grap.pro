QT += core gui

greaterThan(QT_MAJOR_VERSION, 4) : QT += widgets

QMAKE_CC = clang
QMAKE_CXX = clang++

QMAKE_LINK = clang++
QMAKE_LINK_SHLIB = clang++

QMAKE_CXXFLAGS += -Wall -Wextra -Wpedantic -fsanitize=thread

QMAKE_LFLAGS += -fsanitize=thread

SOURCES = main.cpp fill_msr.cpp residual.cpp solve.cpp thread.cpp window.cpp 
HEADERS = fill_msr.hpp residual.hpp solve.hpp thread.hpp window.hpp
TARGET = a.out
