PIN_HOME=../../pin-2.14/

CFLAGS=-DBIGARRAY_MULTIPLIER=1 -Wall -Werror -Wno-unknown-pragmas -g -fno-stack-protector -DTARGET_IA32E -DHOST_IA32E -fPIC -DTARGET_LINUX  -I$(PIN_HOME)/source/include/pin -I$(PIN_HOME)/source/include/pin/gen -I$(PIN_HOME)/extras/components/include -I$(PIN_HOME)/extras/xed-intel64/include -I$(PIN_HOME)/source/tools/InstLib 

LFLAGS=-shared -Wl,--hash-style=sysv -Wl,-Bsymbolic -Wl,--version-script=$(PIN_HOME)/source/include/pin/pintool.ver  

PIN_LIB=-L$(PIN_HOME)/intel64/lib -L$(PIN_HOME)/intel64/lib-ext -L$(PIN_HOME)/intel64/runtime/glibc -L$(PIN_HOME)/extras/xed-intel64/lib -lpin -lxed -lpindwarf -ldl

ifndef DEBUG
OPT_FLAGS=-O2
else
OPT_FLAGS=-O0 -g
endif

all:
	g++-4.9 -std=c++11 $(CFLAGS) $(OPT_FLAGS) -c -o l1.o l1.cpp
	g++-4.9 -std=c++11 $(CFLAGS) $(OPT_FLAGS) -c -o l2.o l2.cpp
	g++-4.9 -std=c++11 $(CFLAGS) $(OPT_FLAGS) -c -o llc.o llc.cpp
	g++-4.9 -std=c++11 $(CFLAGS) $(OPT_FLAGS) -c -o cache_backend.o cache_backend.cpp 
	g++-4.9 -std=c++11 $(CFLAGS) $(OPT_FLAGS) -c -o cache_pinsim.o cache_pinsim.cpp
	g++-4.9 -std=c++11 $(LFLAGS) -o cache_pinsim.so cache_pinsim.o cache_backend.o l1.o l2.o llc.o $(PIN_LIB)

clean:
	rm *.o cache_pinsim.so
