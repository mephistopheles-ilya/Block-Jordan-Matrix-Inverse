ifeq ($(origin CXX),default)
  CXX = g++
endif

# mpirun -np 4 -x ASAN_OPTIONS=detect_leaks=0 ./a.out
# mpirun -np 2 -x ASAN_OPTIONS=detect_leaks=0 xterm -e gdb --args ./a.out
#CXXFLAGS ?=  -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align\
			-Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security\
			-Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long\
			-Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format\

CXXFLAGS ?=  -O0 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align\
			-pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security\
			-Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long\
			-Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format\
			-fsanitize=undefined,address

#CXXFLAGS ?= -fsanitize=leak,undefined,address -g
LDFLAGS ?=  -fsanitize=undefined,address
#CXXFLAGS ?= -O3 -mfpmath=sse

CSRC = main.cpp read_print_fill.cpp
COBJ = main.o read_print_fill.o

%.o : %.cpp
	mpicxx $(CXXFLAGS) -c $^ -o $@

.PHONY: all
all: a.out

a.out: $(COBJ)
	mpicxx $^ -o $@ $(LDFLAGS)

.PHONY: clean
clean:
	rm -rf *.o
