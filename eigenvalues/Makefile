ifeq ($(origin CXX),default)
  CXX = g++
endif

#CXXFLAGS ?=  -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align\
			-Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security\
			-Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long\
			-Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format\
			#-ffast-math #-march=native
			#-ffast-math #-march=native

#CXXFLAGS ?= -g -O0 -fsanitize=leak,undefined,address,null
#LDFLAGS ?=  -fsanitize=leak,undefined,address,null
CXXFLAGS ?= -O3

CSRC = main.cpp func.cpp eigenvalues.cpp
COBJ = main.o func.o eigenvalues.o

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@

.PHONY: all
all: a.out

a.out: $(COBJ)
	$(CXX) $^ -o $@ $(LDFLAGS)

.PHONY: clean
clean:
	rm -rf *.o
