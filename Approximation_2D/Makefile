ifeq ($(origin CXX),default)
  CXX = g++
endif

CXXFLAGS ?=  -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align\
			-Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security\
			-Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long\
			-Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format\
			-pthread -std=c++14 #-fsanitize=leak,undefined,address -fno-sanitize-recover=all

LDFLAGS ?= -pthread #-fsanitize=leak,undefined,address

#CXXFLAGS ?= -O3 -g -mfpmath=sse -Wall -Wextra -Wpedantic -fsanitize=memory
#LDFLAGS ?= -pthread -fsanitize=memory


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
