# See LICENSE.txt for license details.

ifneq ($(DEBUG),)
    CXX_FLAGS += -std=c++11 -O0 -g -Wall
else
    CXX_FLAGS += -std=c++11 -O3 -g -Wall
endif

PAR_FLAG = -fopenmp -march=native

ifneq (,$(findstring icpc,$(CXX)))
	PAR_FLAG = -openmp 
endif

ifneq (,$(findstring sunCC,$(CXX)))
	CXX_FLAGS = -std=c++11 -xO3 -m64 -xtarget=native
	PAR_FLAG = -xopenmp
endif

ifneq ($(SERIAL), 1)
	CXX_FLAGS += $(PAR_FLAG)
endif

KERNELS = pr cc_sv #pr_delta 
SUITE = $(KERNELS) #converter

.PHONY: all
all: $(SUITE)

% : %.cc *.h
	$(CXX) $(CXX_FLAGS) $< -o $@

# Testing
#include test/test.mk

# Benchmark Automation
#include benchmark/bench.mk
#include benchmark/bench_medium.mk


.PHONY: clean
clean:
	rm -f $(SUITE) test/out/*
