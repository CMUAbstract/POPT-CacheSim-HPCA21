#ifndef LLC_BACKEND_H
#define LLC_BACKEND_H

#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <omp.h>
#include <cstring>
#include <cassert>
#include <vector>
#include "pin.H"

typedef uint64_t stat_t;
const int numThreads = 1; //Set this the same as the application threads in the ROI


/* Different data types being tracked */
const int numMainDataTypes = 5; //allow the application to specify these
const int IRREGDATA     {0};
const int REGDATA       {1};
const int CSR_OFFSETS   {2};
const int CSR_COORDS    {3};
const int FRONTIER      {4};
const int OTHER         {5};

class LLC
{
    public:
        LLC();
        ~LLC(); 
        
        void Init();
        void registerDataType(intptr_t addr, int dataTypeID, int numElements, size_t elemSz, int totalDataTypes);
        void insertData(intptr_t addr, int tid, bool isWrite, bool &hit, bool &writeback, intptr_t &evictedAddr, bool updateReplacement); 
        void reportTotalStats();
        void reportDataTypeMisses();
        void reportEvictionReasons();

        //finding accesses to special data types
        bool isSpecialDataType(intptr_t addr, int &dataTypeID);
        
        //stat update functions
		void updateHits(int tid);	
		void updateMisses(int tid);
        void updateWritebacks(int tid);

    private:
        /////////////// DATA ////////////
        /* LLC parameters */
        /* Modeling a 3MB/core, 16-way set associative LLC */
        const int numBanks  = 8; //To simulate a NUCA bank 
        const int numCores  = 8;
        const int m_numSets = 3072 * numCores; 
        const int m_numWays = 16;
        const int m_lineSz  = 64;   

        /* tag store */
        intptr_t** m_tagArray; //[numSets][numWays];
        uint16_t* m_dirty;      
        stat_t** m_referenceCtr; //keep track of no. of references between insertion & eviction
        PIN_LOCK* m_setLocks;  //[numSets]; //per-set lock
        
        /* stats */
        stat_t m_hits[numThreads];
        stat_t m_misses[numThreads];
        stat_t m_writebacks[numThreads];

        // per-datatype LLC hit/miss characterization
        int m_numDataTypes;
        stat_t* m_datatype_misses[numThreads];
        stat_t* m_datatype_hits[numThreads];
        stat_t* m_datatype_writebacks[numThreads];
        stat_t* m_datatype_deadlines[numThreads];
        stat_t* m_datatype_references[numThreads];
       
        // tracking for who evicted who (2d matrix; index 0 for evictee, index 1 for evicter) 
        stat_t** m_eviction_reason[numThreads];
        
        
        //address information for special dataTypes
        std::vector<std::vector<intptr_t> > m_dType_addrStart;
        std::vector<std::vector<intptr_t> > m_dType_addrEnd;
        std::vector<std::vector<int> > m_dType_elemSz;
        
        //LRU related
        uint8_t** m_lru_bits;

        /////////////// METHODS ////////////
        /* internal methods to implement inserting data */
        bool isCacheHit(intptr_t addr, int setID, bool isWrite, bool updateReplacementMetadata);
        bool installNewLine(intptr_t addr, int setID, intptr_t &evictedAddr, bool isWrite, int tid);
		int findSet(intptr_t addr);

        //replacement related (LRU)
        int getReplacementIndex(int setID, int tid);
        void updateReplacementState(int setID, int wayID);
        void moveToMRU(int setID, int wayID);

        // get total number of misses (to implement 1-out-of-2^n counter for BRRIP)
        uint64_t reduceNumMisses();
        
        int findVtxID(intptr_t addr);
};

#endif //LLC_BACKEND_H
