QT += widgets

QMAKE_CXXFLAGS += -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align\
			-Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security\
			-Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long\
			-Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format\
            -fsanitize=leak,undefined,address

QMAKE_LFLAGS += -fsanitize=leak,undefined,address

HEADERS       = window.hpp approximations.hpp
SOURCES       = main.cpp window.cpp approximations.cpp
TARGET        = a.out
