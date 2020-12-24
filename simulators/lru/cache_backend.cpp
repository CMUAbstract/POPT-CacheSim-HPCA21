#include "cache_backend.h"

Cache::Cache() { }

Cache::~Cache() { }

// App programmers responsibility to call with the write dataTypeID args 
void Cache::registerDataType(intptr_t addr, int dataTypeID, int numElements, size_t elemSz, int totalDataTypes = 5)
{
    m_llc.registerDataType(addr, dataTypeID, numElements, elemSz, totalDataTypes);
}


void Cache::Init()
{
    for (int tid = 0; tid < numThreads; ++tid)
    {
        m_l1d[tid].Init();
        m_l1i[tid].Init();
        m_l2[tid].Init();
    }
    m_llc.Init();
}

void Cache::insertData(intptr_t addr, int tid, bool isWrite, bool isInstrn)
{
    bool l1Hit {false};
    bool l1WB {false}; //writeback dirty data from L1
    intptr_t l1EvictAddr {-1};
    
    if (isInstrn == true)
        m_l1i[tid].insertData(addr, tid, isWrite, l1Hit, l1WB, l1EvictAddr);
    else
        m_l1d[tid].insertData(addr, tid, isWrite, l1Hit, l1WB, l1EvictAddr);

    if (l1Hit == false)
    {
       //Missed in L1; access L2
       if (l1WB == true)
       {
           //service the writeback request first
           bool l2Hit {false};
           bool l2WB {false};
           intptr_t l2EvictAddr {-1};
           m_l2[tid].insertData(l1EvictAddr, tid, true, l2Hit, l2WB, l2EvictAddr, false);

           if (l2Hit == false)
           {
                //we tried to update value in L2 but data wasnt there. Forward writeback to LLC
                bool l3Hit {false};
                bool l3WB {false};
                intptr_t l3EvictAddr {-1};
                m_llc.insertData(l1EvictAddr, tid, true, l3Hit, l3WB, l3EvictAddr, false);
           }
       }

       // service actual (demand) miss request
       bool l2Hit {false};
       bool l2WB {false};
       intptr_t l2EvictAddr {-1};
       m_l2[tid].insertData(addr, tid, isWrite, l2Hit, l2WB, l2EvictAddr, true);

       if (l2Hit == false)
       {
           //Missed in L2; time to access LLC
           if (l2WB == true)
           {
               //service the writeback request first
               bool l3Hit {false};
               bool l3WB {false};
               intptr_t l3EvictAddr {-1};
               m_llc.insertData(l2EvictAddr, tid, true, l3Hit, l3WB, l3EvictAddr, false);
           }
            
           #if 0 //no level of cache is inclusive
           //enforce inclusion -> kick data out of L1 (if present)
           if (l2EvictAddr != -1)
           {
               /* NOTE: issuing backinvalidate to both D and I cache
               because we dont track instructions in the Unified L2
               cache */
               m_l1d[tid].backInvalidate(l2EvictAddr, tid);
               m_l1i[tid].backInvalidate(l2EvictAddr, tid);
               //std::cout << "[DEBUG] L2 eviction - recalled data from L1" << std::endl;
           }
           #endif

           //service actual miss request
           bool l3Hit {false};
           bool l3WB {false};
           intptr_t l3EvictAddr {-1};
           m_llc.insertData(addr, tid, isWrite, l3Hit, l3WB, l3EvictAddr, true);

           if (l3Hit == false && l3EvictAddr != -1)
           {
               // LLC is non-inclusive; no need to issue backinvalidates
               // to L1 and L2 on a LLC eviction
               #if 0
               // a line was evicted from LLC; enforce inclusion by removing from all private levels
               for (int thread = 0; thread < numThreads; ++thread)
               {
                   m_l1[thread].backInvalidate(l3EvictAddr, tid);
                   m_l2[thread].backInvalidate(l3EvictAddr, tid);
               }
               //std::cout << "[DEBUG] LLC eviction - recalled data from L1 & L2" << std::endl;
               #endif
           }
       }
    }
}

void Cache::reportTotalStats()
{
    //collect stats for L1D first
    {
        stat_t totalHits {0};
        stat_t totalMisses {0};
        stat_t totalWritebacks {0};
        stat_t totalBackInvalidates {0};

        for (int tid = 0; tid < numThreads; ++tid)
        {
            totalHits += m_l1d[tid].getHits();
            totalMisses += m_l1d[tid].getMisses();
            totalWritebacks += m_l1d[tid].getWritebacks();
            totalBackInvalidates += m_l1d[tid].getBackInvalidates();

            m_l1d[tid].resetCacheStats();
        }
        std::cout << "[L1D-STAT] Total Hits = " << totalHits << std::endl;
        std::cout << "[L1D-STAT] Total Misses = " << totalMisses << std::endl;
        std::cout << "[L1D-STAT] Total Writebacks = " << totalWritebacks << std::endl;
        std::cout << "[L1D-STAT] Total BackInvalidates = " << totalBackInvalidates << std::endl;
        std::cout << "~~~~~~~~~~~~~~~~~~~\n\n";
    }
    
    //collect stats for L1I 
    {
        stat_t totalHits {0};
        stat_t totalMisses {0};
        stat_t totalWritebacks {0};
        stat_t totalBackInvalidates {0};

        for (int tid = 0; tid < numThreads; ++tid)
        {
            totalHits += m_l1i[tid].getHits();
            totalMisses += m_l1i[tid].getMisses();
            totalWritebacks += m_l1i[tid].getWritebacks();
            totalBackInvalidates += m_l1i[tid].getBackInvalidates();
            
            m_l1i[tid].resetCacheStats();
        }
        std::cout << "[L1I-STAT] Total Hits = " << totalHits << std::endl;
        std::cout << "[L1I-STAT] Total Misses = " << totalMisses << std::endl;
        std::cout << "[L1I-STAT] Total Writebacks = " << totalWritebacks << std::endl;
        std::cout << "[L1I-STAT] Total BackInvalidates = " << totalBackInvalidates << std::endl;
        std::cout << "~~~~~~~~~~~~~~~~~~~\n\n";
    }
    
    //collect stats for L2 
    {
        stat_t totalHits {0};
        stat_t totalMisses {0};
        stat_t totalWritebacks {0};
        stat_t totalBackInvalidates {0};

        for (int tid = 0; tid < numThreads; ++tid)
        {
            totalHits += m_l2[tid].getHits();
            totalMisses += m_l2[tid].getMisses();
            totalWritebacks += m_l2[tid].getWritebacks();
            totalBackInvalidates += m_l2[tid].getBackInvalidates();
            
            m_l2[tid].resetCacheStats();
        }
        std::cout << "[L2-STAT] Total Hits = " << totalHits << std::endl;
        std::cout << "[L2-STAT] Total Misses = " << totalMisses << std::endl;
        std::cout << "[L2-STAT] Total Writebacks = " << totalWritebacks << std::endl;
        std::cout << "[L2-STAT] Total BackInvalidates = " << totalBackInvalidates << std::endl;
        std::cout << "~~~~~~~~~~~~~~~~~~~\n\n";
    }

    //finally report stats for L3
    {
        m_llc.reportTotalStats();
        m_llc.reportDataTypeMisses();
        m_llc.reportEvictionReasons();
    }
}

