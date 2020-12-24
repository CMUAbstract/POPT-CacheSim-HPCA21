# Cache Simulators 

This directory contains cache simulators with different Last Level Cache (LLC) 
replacement policies:
* LRU       : Note that this is NOT an approximation of LRU
* DRRIP     : Uses Set-dueling to switch between SRRIP & BRRIP policies (Refer to [1] for more info)
* P-OPT     : A practical OPT implementation that reserves LLC capacity for Rereference Matrix Columns (Refer to [2] for more info)
* OPT-IDEAL : Transpose-driven Belady's OPT replacement policy (incurs zero overhead for accessing rereference information)

## Modifying the simulator

The following sections provide information on the most common modifications to the cache simulator pintool

### Changing parameters of the cache hierarchy

To specify different cache sizes, line sizes, and associativity modify parameters in `l{1,2,lc}.h`. 

The default hierarchy simulates a 3-level cache hierarchy with 8 cores. Each core
has a private L1 and L2 cache with PLRU (Pseudo-LRU) replacement policy and the 
replacement policy for the NUCA (with 8 banks) LLC defines the different simulator versions in this directory.

All levels of cache use 64B cache line size and are non-inclusive. 

### Adding cache replacement policies

We are mostly concerned with the replacement policy at the LLC (for 
most experiments the L1 & L2 replacement policies remain fixed at 
default PLRU policy).

Each new replacement policy must define the following functions:
* `getReplacementIndex(...)`    : Identify the set to be evicted on replacement
* `updateReplacementIndex(...)` : Identify what to do on a cache hit or insertion 

Checking the difference in the implementations of above functions across the 
different simulator versions will give a sense on how to define new policies

### Communicating information from an application to the simulator

To communicate information from the instrumented application to the simulator (pintool), 
we use the pin feature of "function replacement". As the name suggests,
function replacement replaces the function definition provided in the 
application with the definition within the simulator. This allows us
to place "hooks" in the instrumented application for which our simulator
behaves differently. Examples of function replacement are in `cache_pinsim.cpp` (search for `RTN_Replace(...)`).

### Adding new data structures to track

The cache simulator provides global hit/miss statistics for each cache 
level. Additionally, statistics can also be collected on a per data
structure level (NOTE: currently only tracking for arrays is supported)

To add tracking for a new data structure, modify the following:
* `numMainDataTypes` in llc.h
* Assign a new integer identifier for the data structure
* Make sure to use the same integer identifier in the instrumented application (when registering the datatype with `PIN_RegisterDataType(...)`
* Modify statistics-printing functions in the simulator (Mostly in `llc.cpp` and `cache_pinsim.cpp`)

### Modifying data flow between cache levels

Orchestrating the flow of information between different cache levels happens in `cache_backend.cpp`. 

Modifying how requests are sent between different levels of caches will allow implementing advanced
cache features such as selectively caching only some data-structures in an application.

### Multi-threaded runs

To support multiple threads, simply set `NUM_THREADS` in `llc.h` to the number of threads used in the 
multi-threaded applications.

Currently, the simulators will work for serial versions of graph applications.

## References

1. [RRIP-Paper](https://people.csail.mit.edu/emer/papers/2010.06.isca.rrip.pdf)
2. [P-OPT](https://users.ece.cmu.edu/~vigneshb/papers/POPT_HPCA21_CameraReady.pdf)
