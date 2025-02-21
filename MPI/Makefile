ifeq ($(origin CXX),default)
  CXX = g++
endif

# mpirun -np 4 -x ASAN_OPTIONS=detect_leaks=0 ./a.out
# mpirun -np 2 -x ASAN_OPTIONS=detect_leaks=0 xterm -e gdb --args ./a.out
#/lib/x86_64-linux-gnu/mpich/include
#/opt/impi-5.1.3.223/intel64/include
#/usr/lib/x86_64-linux-gnu/openmpi/include 

CXXFLAGS ?= -isystem  /usr/lib/x86_64-linux-gnu/openmpi/include -O3 -fstack-protector-all\
			-g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal\
			-Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings\
			-Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual\
			-Wno-suggest-attribute=format -std=c++17\
			#-fsanitize=undefined,address

#CXXFLAGS ?= -fsanitize=leak,undefined,address -g
#LDFLAGS ?=  -fsanitize=undefined,address
#CXXFLAGS ?= -O3 -mfpmath=sse

CSRC = main.cpp read_print_fill.cpp matrix.cpp inverse.cpp discrepancy.cpp
COBJ = main.o read_print_fill.o matrix.o inverse.o discrepancy.o

%.o : %.cpp
	mpicxx $(CXXFLAGS) -c $^ -o $@

.PHONY: all
all: a.out

a.out: $(COBJ)
	mpicxx $^ -o $@ $(LDFLAGS)

.PHONY: clean
clean:
	rm -rf *.o
