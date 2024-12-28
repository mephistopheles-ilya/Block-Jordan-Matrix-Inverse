ifeq ($(origin CXX),default)
  CXX = g++
endif

CXXFLAGS ?=  -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align\
			-Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security\
			-Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long\
			-Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format\
			-pthread
			#-march=native -ffast-math
			#-ffast-math #-march=native
			#-fsanitize=leak,undefined,address

#CXXFLAGS ?= -fsanitize=leak,undefined,address -g
#LDFLAGS ?=  -fsanitize=leak,undefined,address
#CXXFLAGS ?= -O3 -mfpmath=sse


CSRC = main.cpp read_print_fill.cpp matrix.cpp inverse.cpp discrepancy.cpp thread.cpp
COBJ = main.o read_print_fill.o matrix.o inverse.o discrepancy.o thread.o

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@

.PHONY: all
all: a.out

a.out: $(COBJ)
	$(CXX) $^ -o $@ $(LDFLAGS)

.PHONY: clean
clean:
	rm -rf *.o
