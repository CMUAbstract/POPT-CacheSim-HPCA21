#include "llc.h"

LLC::LLC() { }

LLC::~LLC() 
{
    for (int set = 0; set < m_numSets; ++set)
    {
        delete[] m_tagArray[set];
        delete[] m_lru_bits[set];
        delete[] m_referenceCtr[set];
    }
    delete[] m_tagArray;
    delete[] m_dirty;
    delete[] m_referenceCtr;
    delete[] m_lru_bits;
    delete[] m_setLocks;

    for (int tid = 0; tid < numThreads; ++tid)
    {
        delete[] m_datatype_misses[tid];
        delete[] m_datatype_hits[tid];
        delete[] m_datatype_writebacks[tid];
        delete[] m_datatype_deadlines[tid];
        delete[] m_datatype_references[tid];

        for (int d = 0; d <= m_numDataTypes; ++d)
        {
            delete[] m_eviction_reason[tid][d];
        }
        delete[] m_eviction_reason[tid];
    }
}

void LLC::Init()
{
    for (int tid = 0; tid < numThreads; ++tid)
    {
        m_hits[tid]       = 0;
        m_misses[tid]     = 0;
        m_writebacks[tid] = 0;
        m_datatype_hits[tid] = nullptr;
        m_datatype_misses[tid] = nullptr;
        m_datatype_writebacks[tid] = nullptr;
        m_datatype_deadlines[tid] = nullptr;
        m_datatype_references[tid] = nullptr;
        m_eviction_reason[tid] = nullptr;
    }

    m_dType_addrStart.resize(numMainDataTypes);
    m_dType_addrEnd.resize(numMainDataTypes);
    m_dType_elemSz.resize(numMainDataTypes);

    m_tagArray = new intptr_t* [m_numSets];
    m_dirty    = new uint16_t [m_numSets];
    m_referenceCtr = new stat_t* [m_numSets];
    m_setLocks = new PIN_LOCK [m_numSets];
    m_lru_bits  = new uint8_t* [m_numSets];
    assert(m_tagArray != nullptr && m_dirty != nullptr && m_setLocks != nullptr && m_lru_bits != nullptr);
    for (int set = 0; set < m_numSets; ++set)
    {
        PIN_InitLock(&m_setLocks[set]);
        m_tagArray[set] = new intptr_t [m_numWays];
        m_lru_bits[set]  = new uint8_t [m_numWays];
        m_referenceCtr[set] = new stat_t [m_numWays];
        assert(m_tagArray[set] != nullptr && m_lru_bits[set] != nullptr);
        m_dirty[set] = 0;
        for (int way = 0; way < m_numWays; ++way)
        {
            m_tagArray[set][way] = -1; //the entire array is invalid initially
            m_lru_bits[set][way]     = way;
            m_referenceCtr[set][way] = 0;
        }
    }
}

// App programmers responsibility to call with the write dataTypeID args 
void LLC::registerDataType(intptr_t addr, int dataTypeID, int numElements, size_t elemSz, int totalDataTypes)
{
    if (m_datatype_misses[0] == nullptr)
    {
        /* data structures havent been initialized - initialize them */
        m_numDataTypes = totalDataTypes;
        for (int tid = 0; tid < numThreads; ++tid)
        {
            m_datatype_misses[tid]     = new stat_t [m_numDataTypes](); 
            m_datatype_hits[tid]       = new stat_t [m_numDataTypes](); 
            m_datatype_writebacks[tid] = new stat_t [m_numDataTypes](); 
            m_datatype_deadlines[tid]  = new stat_t [m_numDataTypes](); 
            m_datatype_references[tid] = new stat_t [m_numDataTypes](); 
            m_eviction_reason[tid]     = new stat_t* [m_numDataTypes+1];
            for (int d = 0; d <= m_numDataTypes; ++d)
            {
                m_eviction_reason[tid][d] = new stat_t[m_numDataTypes+1]();
            }
        }
    }

    assert(dataTypeID < m_numDataTypes);
    if (dataTypeID == FRONTIER)
    {
      // Slightly different tracking for frontiers because of double-buffering
      if (m_dType_addrStart[dataTypeID].size() == 0)
      {
        m_dType_addrStart[dataTypeID].push_back(addr);
        m_dType_addrEnd[dataTypeID].push_back(addr + (elemSz * numElements));
        m_dType_elemSz[dataTypeID].push_back(elemSz);
      }
      else 
      {
        m_dType_addrStart[dataTypeID][0] = addr;
        m_dType_addrEnd[dataTypeID][0]   = addr + (elemSz * numElements);
        m_dType_elemSz[dataTypeID][0]    = elemSz;
      }
    }
    else
    {
      m_dType_addrStart[dataTypeID].push_back(addr);
      m_dType_addrEnd[dataTypeID].push_back(addr + (elemSz * numElements));
      m_dType_elemSz[dataTypeID].push_back(elemSz);
    }


    if (dataTypeID == REGDATA || dataTypeID == IRREGDATA || dataTypeID == FRONTIER)
    {
        //std::cout << "[SIM] Alignment error for data structure no. " << dataTypeID << std::endl;
        assert(addr % 64 == 0);
    }
}
 

void LLC::insertData(intptr_t addr, int tid, bool isWrite, bool &hit, bool &writeback, intptr_t &evictedAddr, bool updateReplacement)
{
    int setID = findSet(addr);
    //omp_set_lock(m_setLocks[setID]);
    PIN_GetLock(&m_setLocks[setID], tid+1);
    {
        //entered critical section
        hit = isCacheHit(addr, setID, isWrite, updateReplacement);
        if (hit == false && updateReplacement == true)
        {
            // NOTE: Assumption that if a L2 WB misses in the LLC, then
            // we forward the request to DRAM and dont bring it to the LLC
            writeback = installNewLine(addr, setID, evictedAddr, isWrite, tid);
        }
        //leaving critical section
    }
    //omp_unset_lock(m_setLocks[setID]);
    PIN_ReleaseLock(&m_setLocks[setID]);

    if (updateReplacement == false)
      return; //dont update hit/miss for evictions

    if (hit == true)
    {
      updateHits(tid);
      int dTypeID {-1};
      // characterizing miss data
      bool specialData = isSpecialDataType(addr, dTypeID);
      if (specialData == true)
      {
          m_datatype_hits[tid][dTypeID] += 1;
      }
    }
    else
    {
        if (writeback == true)
        {
            updateWritebacks(tid);
            assert(evictedAddr != -1 * m_lineSz);
            int dTypeID {-1};
            // characterizing miss data
            bool specialData = isSpecialDataType(evictedAddr, dTypeID);
            if (specialData == true)
            {
                m_datatype_writebacks[tid][dTypeID] += 1;
            }
        }

        if (evictedAddr == -1 * m_lineSz)
            evictedAddr = -1;

        updateMisses(tid);
        int dTypeID {-1};
        // characterizing miss data
        bool specialData = isSpecialDataType(addr, dTypeID);
        if (specialData == true)
        {
            m_datatype_misses[tid][dTypeID] += 1;
        }
        
        if (dTypeID == -1)
            dTypeID = OTHER;
        
        // noting who caused the eviction
        if (evictedAddr != -1)
        {
            int dType_evicted {-1};
            isSpecialDataType(evictedAddr, dType_evicted); 
            if (dType_evicted == -1)
                dType_evicted = m_numDataTypes;
            if (dTypeID == -1)
                dTypeID = m_numDataTypes;
            m_eviction_reason[tid][dType_evicted][dTypeID] += 1;
        }
    }
    return;
}

/* NOTE: the reads here are racy */
uint64_t LLC::reduceNumMisses()
{
    uint64_t misses {0};
    for (int tid = 0; tid < numThreads; ++tid)
    {
        misses += m_misses[tid];
    }
    return misses;
} 

void LLC::reportTotalStats()
{
    stat_t totalHits {0};
    for (int tid = 0; tid < numThreads; ++tid)
    {
        totalHits += m_hits[tid];
        m_hits[tid] = 0;
    }
    std::cout << "[LLC-STAT] Total Hits = " << totalHits << std::endl;

    stat_t totalMisses {0};
    for (int tid = 0; tid < numThreads; ++tid)
    {
        totalMisses += m_misses[tid];
        m_misses[tid] = 0;
    }
    std::cout << "[LLC-STAT] Total Misses = " << totalMisses << std::endl;

    stat_t totalWritebacks {0};
    for (int tid = 0; tid < numThreads; ++tid)
    {
        totalWritebacks += m_writebacks[tid];
        m_writebacks[tid] = 0;
    }
    std::cout << "[LLC-STAT] Total Writebacks = " << totalWritebacks << std::endl;
}

void LLC::reportDataTypeMisses()
{
    std::string names[] = {"irregData", "regData", "CSR-offsets", "CSR-coords", "frontier"};
    for (int dType = 0; dType < m_numDataTypes; ++dType)
    {
        stat_t totalHits   = 0;
        stat_t totalMisses = 0;
        stat_t totalWritebacks = 0;
        stat_t totalDeadlines = 0;
        stat_t totalReferences = 0;
        for (int tid = 0; tid < numThreads; ++tid)
        {
            totalHits   += m_datatype_hits[tid][dType];
            totalMisses += m_datatype_misses[tid][dType];
            totalWritebacks += m_datatype_writebacks[tid][dType];
            totalDeadlines += m_datatype_deadlines[tid][dType];
            totalReferences += m_datatype_references[tid][dType];

            m_datatype_hits[tid][dType] = 0;
            m_datatype_misses[tid][dType] = 0;
            m_datatype_writebacks[tid][dType] = 0;
            m_datatype_deadlines[tid][dType] = 0;
            m_datatype_references[tid][dType] = 0;
        }
        if (dType < numMainDataTypes)
        {
            std::cout << "[LLC-STAT] Hits for DataType - " << names[dType] << " - " << totalHits << std::endl;
            std::cout << "[LLC-STAT] Misses for DataType - " << names[dType] << " - " << totalMisses << std::endl;
            std::cout << "[LLC-STAT] Writebacks for DataType - " << names[dType] << " - " << totalWritebacks << std::endl;
            std::cout << "[LLC-STAT] Deadlines for DataType - " << names[dType] << " - " << totalDeadlines << std::endl;
            std::cout << "[LLC-STAT] Avg. ReReference for DataType - " << names[dType] << " - " << static_cast<double>(totalReferences) / static_cast<double>(totalMisses) << std::endl;
        }
        else 
        {
            assert(false);
        }
    }
}

void LLC::reportEvictionReasons()
{
    std::string names[] = {"irregData", "regData", "CSR-offsets", "CSR-coords", "frontier", "others"};

    #if 0
    for (int d = 0; d <= numDataTypes; ++d)
    {
        stat_t evictCtr = 0;
        for (int tid = 0; tid < numThreads; ++tid)
        {
            evictCtr += m_eviction_reason[tid][d];
        }
        std::cout << "[LLC-STAT] No. of Evictions by DataType - " << names[d] << " - " << evictCtr << std::endl;
    }
    #endif

    for (int dEvictee = 0; dEvictee <= m_numDataTypes; ++dEvictee)
    {
        if (dEvictee < numMainDataTypes)
            std::cout << "[LLC-EVICT-STAT] Eviction stats for DataType - " << names[dEvictee] << std::endl;
        else
            std::cout << "[LLC-EVICT-STAT] Eviction stats for DataType - others" << std::endl;

        for (int dEvicter = 0; dEvicter <= m_numDataTypes; ++dEvicter)
        {
            stat_t evictCtr = 0;
            for (int tid = 0; tid < numThreads; ++tid)
            {
                evictCtr += m_eviction_reason[tid][dEvictee][dEvicter];
                m_eviction_reason[tid][dEvictee][dEvicter] = 0;
            }

            if (dEvicter < numMainDataTypes)
                std::cout << "[LLC-STAT] No. of Evictions by DataType - " << names[dEvicter] << " - " << evictCtr << std::endl;
            else
                std::cout << "[LLC-STAT] No. of Evictions by DataType - others - " << evictCtr << std::endl;
            
        }
        std::cout << "[LLC-EVICT-STAT] ~~~~~~~~~~~~~~~~~~~~~~~~\n";
    }
}


bool LLC::isCacheHit(intptr_t addr, int setID, bool isWrite, bool updateReplacementMetadata)
{
    intptr_t maskedAddr = addr / m_lineSz;
    for (int way = 0; way < m_numWays; ++way)
    {
        if ((m_tagArray[setID][way] / m_lineSz) == maskedAddr)
        {
            if (updateReplacementMetadata)
                updateReplacementState(setID, way);
            if (isWrite == true)
            {
                //m_dirty[setID][way] = 1;
                uint16_t mask = 1 << way;
                m_dirty[setID] = m_dirty[setID] | (mask);
            }
            m_referenceCtr[setID][way]++;
            return true;
        }
    }
    return false;
}

bool LLC::installNewLine(intptr_t addr, int setID, intptr_t &evictedAddr, bool isWrite, int tid)
{
    int index = getReplacementIndex(setID, tid);
    //evictedAddr = m_tagArray[setID][index] * m_lineSz; //line to be kicked out
    evictedAddr = m_tagArray[setID][index]; //line to be kicked out
    uint16_t mask = 1 << index;
    uint16_t retVal = m_dirty[setID] & mask; 
    assert(retVal == 0 || retVal == mask);
    if (evictedAddr == -1 * m_lineSz)
        assert(retVal == 0);

    //check how many times the line was rereferenced
    stat_t reReference = m_referenceCtr[setID][index]; 
    int dTypeID {-1}; 
    bool specialData = isSpecialDataType(evictedAddr, dTypeID);
    if (specialData == true)
    {
        if (reReference == 0)
        {
            // this is a dead-line (no reuse between insertion and eviction)
            m_datatype_deadlines[tid][dTypeID] += 1;
        }

        // update total number of references
        m_datatype_references[tid][dTypeID] += reReference;
    }
    //reset ctr for newly inserted line
    m_referenceCtr[setID][index] = 0;

    //m_tagArray[setID][index] = addr / m_lineSz; //new line inserted 
    m_tagArray[setID][index] = addr; //new line inserted 
    //m_dirty[setID][index]    = (isWrite == true) ? 1 : 0;
    if (isWrite == true)
        m_dirty[setID] = m_dirty[setID] | (mask);
    else
        m_dirty[setID] = m_dirty[setID] & (~mask);
    return (retVal != 0);
}


/* Set mapping policy from Sniper */ 
int LLC::findSet(intptr_t addr)
{
    // RNG parameters, defaults taken from drand48
    #define RNG_A __UINT64_C(0x5DEECE66D)
    #define RNG_C __UINT64_C(0xB)
    #define RNG_M ((__UINT64_C(1) << 48) - 1)

    assert(m_lineSz == 64);
    uint64_t block_num = addr >> 6;
    uint64_t state = (block_num << 16) + 0x330E;
    state          = (RNG_A * state + RNG_C) & RNG_M;
    int setIndex   = (state >> 16) % m_numSets;
    return setIndex;
}

int LLC::getReplacementIndex(int setID, int tid)
{
    // check if there is a way with an invalid line
    for (int way = 0; way < m_numWays; ++way)
    {
        if (m_tagArray[setID][way] == -1)
        {
            moveToMRU(setID, way);
            return way;
        }
    }
    
    // all ways contain valid data
    // unlike sniper code we are not implementing QBS
    uint32_t index   = 0;
    uint8_t max_bits = 0;
    for (int way = 0; way < m_numWays; ++way)
    {
        if (m_lru_bits[setID][way] > max_bits)
        {
            index    = way;
            max_bits = m_lru_bits[setID][way];
        }
    }
    moveToMRU(setID, index);
    return index;
}

void LLC::updateReplacementState(int setID, int wayID)
{
    moveToMRU(setID, wayID);
}

void LLC::moveToMRU(int setID, int wayID)
{
    for (int way = 0; way < m_numWays; ++way)
    {
       if (m_lru_bits[setID][way] < m_lru_bits[setID][wayID])
          m_lru_bits[setID][way]++;
    }
    m_lru_bits[setID][wayID] = 0;
}


void LLC::updateHits(int tid)
{
	m_hits[tid] += 1;
}


void LLC::updateMisses(int tid)
{
	m_misses[tid] += 1;
}

void LLC::updateWritebacks(int tid)
{
    m_writebacks[tid] += 1;
}

/* For srcData accesses (find if this was a hub or non-hub) */
int LLC::findVtxID(intptr_t addr)
{
    assert(m_dType_addrStart[IRREGDATA].size() == 1);
    assert(m_dType_addrEnd[IRREGDATA].size() == 1);
    assert(m_dType_elemSz[IRREGDATA].size() == 1);
    auto startAddr = m_dType_addrStart[IRREGDATA][0]; 
    auto endAddr   = m_dType_addrEnd[IRREGDATA][0]; 
    auto elemSz    = m_dType_elemSz[IRREGDATA][0];
    
    assert(addr >= startAddr && addr < endAddr);

    auto diff = addr - startAddr;
    int vtxID = (diff) / elemSz;
    return vtxID;
}

bool LLC::isSpecialDataType(intptr_t addr, int &dataTypeID)
{
    if (addr == -1)
        return false;
    assert(addr > 0);
    for (int d = 0; d < m_numDataTypes; ++d) 
    {
        for (size_t i = 0; i < m_dType_addrStart[d].size(); ++i)
        {
            if (addr >= m_dType_addrStart[d][i] && addr < m_dType_addrEnd[d][i])
            {
                dataTypeID = d;
                return true;
            }
        }
    }
    return false;
}

