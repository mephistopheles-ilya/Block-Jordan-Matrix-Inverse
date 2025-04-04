ifeq ($(origin CXX),default)
  CXX = g++
endif

CXXFLAGS ?=  -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align\
			-Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security\
			-Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long\
			-Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format\
			-pthread -std=c++14 #-fsanitize=thread#leak,undefined,address
			#-march=native -ffast-math
			#-ffast-math #-march=native
			#-fsanitize=leak,undefined,address

#CXXFLAGS ?= -fsanitize=leak,undefined,address -g
#LDFLAGS ?= -fsanitize=thread#leak,undefined,address
#CXXFLAGS ?= -O3 -mfpmath=sse
LDFLAGS ?= -pthread


CSRC = main.cpp thread.cpp fill_msr.cpp solve.cpp residuals.cpp
COBJ = main.o thread.o fill_msr.o solve.o residuals.o

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@

.PHONY: all
all: a.out

a.out: $(COBJ)
	$(CXX) $^ -o $@ $(LDFLAGS)

.PHONY: clean
clean:
	rm -rf *.o
