#ifndef L1_BACKEND_H
#define L1_BACKEND_H

#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <omp.h>
#include <cstring>
#include <cassert>
#include "pin.H"

typedef uint64_t stat_t;

class L1
{
    public:
        L1();
        ~L1(); 
        
        void Init();
        void insertData(intptr_t addr, int tid, bool isWrite, bool &hit, bool &writeback, intptr_t &evictedAddr); 
        stat_t getHits();
        stat_t getMisses();
        stat_t getWritebacks();
        stat_t getBackInvalidates();

        void resetCacheStats();

        //inclusion related
        void backInvalidate(intptr_t addr, int tid);
		
    private:
        /////////////// DATA ////////////
        /* Modeling a 32KB, 8-way set associative L1 */
        int m_numSets = 64; 
        int m_numWays = 8;
        int m_lineSz  = 64;   

        /* tag store */
        intptr_t** m_tagArray; //[numSets][numWays];
        uint8_t* m_dirty;      //one hot encoding for dirty within ways (supports upto 256 ways)

        #if 0 //we needed locks to protect from backinvalidates from
              //LLC. Since LLC is non-inclusive now, we dont need no locks
        PIN_LOCK m_setLock;    //since we wont have a lot of parallel L1 accesses - using a giant lock for entire L1 (to save memory)
        #endif
        
        /* stats */
        stat_t m_hits;
        stat_t m_misses;
        stat_t m_writebacks;
        stat_t m_backInvalidates;
        
        // PLRU replacement related
        uint8_t** m_plru_bits;

        /////////////// METHODS ////////////
        /* internal methods to implement inserting data */
        bool isCacheHit(intptr_t addr, int setID, bool isWrite);
        bool installNewLine(intptr_t addr, int setID, intptr_t &evictedAddr, bool isWrite);
		int findSet(intptr_t addr);

        //replacement related (LRU)
        int getReplacementIndex(int setID);
        void updateReplacementState(int setID, int wayID);
        void moveToMRU(int setID, int wayID);
        void moveToLRU(int setID, int wayID);
};

#endif //L1_BACKEND_H
