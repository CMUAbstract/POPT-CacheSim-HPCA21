PIN_HOME=../../pin-3.21/

CFLAGS=-DBIGARRAY_MULTIPLIER=1 -Wall -Werror -Wno-unknown-pragmas -DPIN_CRT=1 -fno-stack-protector -fno-exceptions -funwind-tables -fasynchronous-unwind-tables -fno-rtti -DTARGET_IA32E -DHOST_IA32E -fPIC -DTARGET_LINUX -fabi-version=2 -faligned-new  -I$(PIN_HOME)/source/include/pin -I$(PIN_HOME)/source/include/pin/gen -isystem $(PIN_HOME)/extras/stlport/include -isystem $(PIN_HOME)/extras/libstdc++/include -isystem $(PIN_HOME)/extras/crt/include -isystem $(PIN_HOME)/extras/crt/include/arch-x86_64 -isystem $(PIN_HOME)/extras/crt/include/kernel/uapi -isystem $(PIN_HOME)/extras/crt/include/kernel/uapi/asm-x86 -I$(PIN_HOME)/extras/components/include -I$(PIN_HOME)/extras/xed-intel64/include/xed -I$(PIN_HOME)/source/tools/Utils -I$(PIN_HOME)/source/tools/InstLib

LFLAGS=-shared -Wl,--hash-style=sysv $(PIN_HOME)/intel64/runtime/pincrt/crtbeginS.o -Wl,-Bsymbolic -Wl,--version-script=$(PIN_HOME)/source/include/pin/pintool.ver -fabi-version=2  -g

PIN_LIB=-L$(PIN_HOME)/intel64/runtime/pincrt -L$(PIN_HOME)/intel64/lib -L$(PIN_HOME)/intel64/lib-ext -L$(PIN_HOME)/extras/xed-intel64/lib -lpin -lxed $(PIN_HOME)/intel64/runtime/pincrt/crtendS.o -lpin3dwarf  -ldl-dynamic -nostdlib -lstlport-dynamic -lm-dynamic -lc-dynamic -lunwind-dynamic

ifndef DEBUG
OPT_FLAGS=-O2
else
OPT_FLAGS=-O0 -g
endif

all:
	g++ -std=c++11 $(CFLAGS) $(OPT_FLAGS) -c -o l1.o l1.cpp
	g++ -std=c++11 $(CFLAGS) $(OPT_FLAGS) -c -o l2.o l2.cpp
	g++ -std=c++11 $(CFLAGS) $(OPT_FLAGS) -c -o llc.o llc.cpp
	g++ -std=c++11 $(CFLAGS) $(OPT_FLAGS) -c -o cache_backend.o cache_backend.cpp 
	g++ -std=c++11 $(CFLAGS) $(OPT_FLAGS) -c -o cache_pinsim.o cache_pinsim.cpp
	g++ -std=c++11 $(LFLAGS) -o cache_pinsim.so cache_pinsim.o cache_backend.o l1.o l2.o llc.o $(PIN_LIB)

clean:
	rm *.o cache_pinsim.so
