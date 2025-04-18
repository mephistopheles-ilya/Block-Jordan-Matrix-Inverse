QT += core gui

greaterThan(QT_MAJOR_VERSION, 4) : QT += widgets

QMAKE_CC = g++
QMAKE_CXX = g++

QMAKE_LINK = g++
QMAKE_LINK_SHLIB = g++


QMAKE_CXXFLAGS += -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align\
                  -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security\
                  -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long\
                 -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format -pthread\
                 -fsanitize=address,undefined,integer-divide-by-zero,null,signed-integer-overflow,float-divide-by-zero,float-cast-overflow 

QMAKE_LFLAGS += -pthread -fsanitize=address,undefined,integer-divide-by-zero,null,signed-integer-overflow,float-divide-by-zero,float-cast-overflow 

SOURCES = main.cpp fill_msr.cpp residuals.cpp solve.cpp window.cpp 
HEADERS = fill_msr.hpp residuals.hpp solve.hpp window.hpp
TARGET = a.out
