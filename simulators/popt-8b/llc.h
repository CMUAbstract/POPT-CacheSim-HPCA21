#ifndef LLC_BACKEND_H
#define LLC_BACKEND_H

#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <omp.h>
#include <cstring>
#include <cassert>
#include <vector>
#include <algorithm>
#include "pin.H"
#include "graph.h"

typedef uint64_t stat_t;
const int numThreads = 1; //Set this the same as the application threads in the ROI 

typedef CSRGraph<int32_t> Graph;

const uint8_t MAX_REREF {127}; //The maxReref value that the offsetMatrix can track


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
        void registerGraph(Graph& g, bool isPull);
        void registerOffsetMatrix(uint8_t* offsetMatrix, int32_t numEpochs, int32_t numCacheLines, int vtxPerLine, int dTypeID);
        void updateRegIndex(int32_t index, int tid);
        void updateIrregRanges(int32_t startID, int32_t endID);
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

        // DRRIP replacement related
        uint8_t** m_rrpv;
        uint8_t m_MAX_rrpv = 7; // distant re-reference = 2^M-1 (for M = 3)
        int* m_policySelector;
        const int m_pSEL_max {1024};
        const int m_pSEL_min {0};
        const int m_pSEL_init {512};
        const int FOLLOWER_SRRIP  {0};
        const int FOLLOWER_BRRIP  {1};
        const int DEDICATED_SRRIP {2};
        const int DEDICATED_BRRIP {3};
        const int m_SRRIP {1};
        const int m_BRRIP {2};

        //Graph for finding precise re-reference value
        Graph m_graph;
        bool m_isPull;
        int32_t m_regIndex[numThreads];
        PIN_MUTEX m_dstLock;

        //Data structure for POPT replacement
        pvector<uint8_t*> m_offsetMatrix;
        pvector<int32_t>  m_numEpochs;
        pvector<int32_t>  m_numCacheLines;
        pvector<int32_t>  m_vtxPerLine;
        int m_startWay;

        //Metadata to find when srcData values are live
        int32_t m_irregData_startID;
        int32_t m_irregData_endID;

        /////////////// METHODS ////////////
        /* internal methods to implement inserting data */
        bool isCacheHit(intptr_t addr, int setID, bool isWrite, bool updateReplacementMetadata);
        bool installNewLine(intptr_t addr, int setID, intptr_t &evictedAddr, bool isWrite, int setType, int tid);
		int findSet(intptr_t addr);

        //replacement related (LRU)
        int getReplacementIndex(intptr_t addr, int setID, int setType, int tid);
        void updateReplacementState(int setID, int wayID);
        void moveToMRU(int setID, int wayID);
        uint8_t findRereferenceVal(int vtxID, int currDst, int tid, int dTypeID);

        int findVtxID(intptr_t addr, int dTypeID);

        // get total number of misses (to implement 1-out-of-2^n counter for BRRIP)
        uint64_t reduceNumMisses();

        //DRRIP-specific
        void incPSEL(int bankID);
        void decPSEL(int bankID);
        bool winningPolicy(int bankID);
        int determineSetType(int setID, int bankID);
};

#endif //LLC_BACKEND_H
