#include "llc.h"

LLC::LLC() { }

LLC::~LLC()
{
    for (int set = 0; set < m_numSets; ++set)
    {
        delete[] m_tagArray[set];
        delete[] m_rrpv[set];
        delete[] m_referenceCtr[set];
    }
    delete[] m_tagArray;
    delete[] m_dirty;
    delete[] m_referenceCtr;
    delete[] m_rrpv;
    delete[] m_setLocks;
    delete[] m_policySelector;

    for (int d = 0; d < m_numDataTypes; ++d)
    { 
      if (m_offsetMatrix[d] != nullptr)
        delete[] m_offsetMatrix[d];
    }

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
         
        m_regIndex[tid]  = -1;
    }

    m_dType_addrStart.resize(numMainDataTypes);
    m_dType_addrEnd.resize(numMainDataTypes);
    m_dType_elemSz.resize(numMainDataTypes);
     
    m_tagArray = new intptr_t* [m_numSets];
    m_dirty    = new uint16_t [m_numSets];
    m_referenceCtr = new stat_t* [m_numSets];
    m_setLocks = new PIN_LOCK [m_numSets];
    PIN_MutexInit(&m_dstLock);
    m_rrpv  = new uint8_t* [m_numSets];
    assert(m_tagArray != nullptr && m_dirty != nullptr && m_setLocks != nullptr && m_rrpv != nullptr);
    for (int set = 0; set < m_numSets; ++set)
    {
        PIN_InitLock(&m_setLocks[set]);
        m_tagArray[set] = new intptr_t [m_numWays];
        m_rrpv[set] = new uint8_t [m_numWays];
        m_referenceCtr[set] = new stat_t [m_numWays];
        assert(m_tagArray[set] != nullptr && m_rrpv[set] != nullptr);
        m_dirty[set] = 0;
        for (int way = 0; way < m_numWays; ++way)
        {
            m_tagArray[set][way] = -1; //the entire array is invalid initially
            m_rrpv[set][way] = m_MAX_rrpv;
            m_referenceCtr[set][way] = 0;
        }
    }
    m_irregData_startID = -1;
    m_irregData_endID   = -1;
    m_policySelector = new int[numBanks];
    assert(m_policySelector != nullptr);
    for (int b = 0; b < numBanks; ++b)
        m_policySelector[b] = m_pSEL_init;
    
    for (int d = 0; d < numMainDataTypes; ++d)
    {
      m_offsetMatrix.push_back(nullptr);
      m_numEpochs.push_back(-1);
      m_numCacheLines.push_back(-1);
      m_vtxPerLine.push_back(-1);
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

void LLC::registerGraph(Graph &g, bool isPull)
{
    m_graph.setGraphProperties(g.num_nodes(), g.num_edges(), g.directed());
    m_graph.setGraphDatastructures(g.out_index(), g.out_neighbors(),
                                   g.in_index(), g.in_neighbors());
    m_isPull = isPull;

}

void LLC::registerOffsetMatrix(uint8_t* offsetMatrix, int32_t numEpochs, int32_t numCacheLines, int vtxPerLine, int dTypeID)
{
  m_numEpochs[dTypeID]     = numEpochs;
  m_numCacheLines[dTypeID] = numCacheLines;
  m_vtxPerLine[dTypeID]    = vtxPerLine;
  m_offsetMatrix[dTypeID]  = new uint8_t [numEpochs * numCacheLines];
  assert(m_offsetMatrix[dTypeID] != nullptr);
  
  for (int i = 0; i < numEpochs * numCacheLines; ++i)
    m_offsetMatrix[dTypeID][i] = offsetMatrix[i];
  
  assert(dTypeID == FRONTIER || dTypeID == IRREGDATA);

  int sizePerWay = (m_numSets * m_numWays * m_lineSz) / (m_numWays); //(3 * 1024 * 1024) / 2; 
  if (dTypeID == IRREGDATA)
  {
    if (m_offsetMatrix[FRONTIER] != nullptr)
    {
      //SpMSpV execution
      int frontExtra {0};
      int irregExtra {0};
      if (m_numCacheLines[FRONTIER] % numBanks != 0)
        frontExtra = numBanks - (m_numCacheLines[FRONTIER] % numBanks);
      if (m_numCacheLines[IRREGDATA] % numBanks != 0)
        irregExtra = numBanks - (m_numCacheLines[IRREGDATA] % numBanks);

      int frontierMetadata = m_numCacheLines[FRONTIER] + frontExtra;  //Tracking the matrices for each datatype & epoch will become easier
      int irregMetadata    = m_numCacheLines[IRREGDATA] + irregExtra;

      m_startWay = (((frontierMetadata + irregMetadata) * 2) + sizePerWay - 1) / sizePerWay; //need to store current & next Epoch for two datatypes
    }
    else
    {
      //SpMDV execution
      int irregExtra {0};
      if (m_numCacheLines[IRREGDATA] % numBanks != 0)
        irregExtra = numBanks - (m_numCacheLines[IRREGDATA] % numBanks);

      int irregMetadata = m_numCacheLines[IRREGDATA] + irregExtra;

      m_startWay = (((irregMetadata) * 2) + sizePerWay - 1) / sizePerWay; //need to store current & next Epoch for two datatypes
    }
  }
  std::cout << "[OPT-METADATA-SZ] Start Way = " << m_startWay << std::endl; 
}

void LLC::updateRegIndex(int32_t index, int tid)
{
    m_regIndex[tid] = index;
}

void LLC::updateIrregRanges(int32_t startID, int32_t endID)
{
    m_irregData_startID = startID;
    m_irregData_endID   = endID;
}


void LLC::insertData(intptr_t addr, int tid, bool isWrite, bool &hit, bool &writeback, intptr_t &evictedAddr, bool updateReplacement)
{
    int setID = findSet(addr);
    assert(m_lineSz == 64);
    int bankID = (addr >> 6) % numBanks; //S-NUCA policy
    int setType   = determineSetType(setID, bankID);
    //omp_set_lock(m_setLocks[setID]);
    PIN_GetLock(&m_setLocks[setID], tid+1);
    {
        //entered critical section
        hit = isCacheHit(addr, setID, isWrite, updateReplacement);
        if (hit == false && updateReplacement == true)
        {
            // NOTE: Assumption that if a L2 WB misses in the LLC, then
            // we forward the request to DRAM and dont bring it to the LLC
            writeback = installNewLine(addr, setID, evictedAddr, isWrite, setType, tid);
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

bool LLC::installNewLine(intptr_t addr, int setID, intptr_t &evictedAddr, bool isWrite, int setType, int tid)
{
    int index = getReplacementIndex(addr, setID, setType, tid);
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
    assert(addr != 0);
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

int LLC::getReplacementIndex(intptr_t addr, int setID, int setType, int tid)
{
    if (m_irregData_startID != -1) //should execute OPT
    {
        int evictionWay {-1};
        // check if there is a way with an invalid line
        for (int way = m_startWay; way < m_numWays; ++way)
        {
            if (m_tagArray[setID][way] == -1)
            {
                evictionWay = way;
                break;
            }
        }
        
        if (evictionWay == -1)
        {
          // all ways have valid data
          // check if a way contains data besides irregData & frontier
          for (int way = m_startWay; way < m_numWays; ++way)
          {
              int dTypeID {-1};
              isSpecialDataType(m_tagArray[setID][way], dTypeID);
              if (dTypeID != IRREGDATA && dTypeID != FRONTIER)
              {
                  evictionWay = way;
                  break;
              }
          }
        }

        if (evictionWay == -1)
        {
          // all ways contain irregData & frontier 
          // find the line that is going to be accessed farthest into the future
          uint8_t maxRerefDist {0};
          int currRegIndex {-1};
          currRegIndex = m_regIndex[0];  //TODO: Fix this racy update
          uint8_t wayRerefDists [m_numWays] {0};
          for (int way = m_startWay; way < m_numWays; ++way)
          {
              uint8_t rerefDist {MAX_REREF};
              intptr_t wayAddr = m_tagArray[setID][way] / m_lineSz;
              int dTypeID {-1};
              isSpecialDataType(wayAddr * m_lineSz, dTypeID);
              assert(dTypeID == FRONTIER || dTypeID == IRREGDATA);
              int baseVtx = findVtxID(wayAddr * m_lineSz, dTypeID);
              rerefDist   = findRereferenceVal(baseVtx, currRegIndex, tid, dTypeID);
              assert(rerefDist <= MAX_REREF);

              if (maxRerefDist < rerefDist)
              {
                  maxRerefDist = rerefDist;
                  evictionWay  = way;
              }
              wayRerefDists[way] = rerefDist;
          }
        
          //To handle ties in reref dist we also factor in DRRIP information
          while (true)
          {
             //search for a way with distant re-reference
             for (int way = m_startWay; way < m_numWays; ++way)
             {
                 if (wayRerefDists[way] == maxRerefDist && m_rrpv[setID][way] == m_MAX_rrpv)
                 {
                     evictionWay = way;
                     break;
                 }
             }

             if (evictionWay != -1)
               break;

             //if we reached here then no line had distance reference. 
             // age all rrpv among competing reref values 
             for (int way = m_startWay; way < m_numWays; ++way)
             {
               if (wayRerefDists[way] == maxRerefDist)
               {
                 ++m_rrpv[setID][way];
                 if (m_rrpv[setID][way] > m_MAX_rrpv)
                     m_rrpv[setID][way] = m_MAX_rrpv;
               }
             }
          }
        }

        assert(evictionWay != -1);
        if (setType == FOLLOWER_BRRIP || setType == DEDICATED_BRRIP)
        {
            if (reduceNumMisses() % 32 == 0)
                m_rrpv[setID][evictionWay] = m_MAX_rrpv-1; //long re-reference
            else
                m_rrpv[setID][evictionWay] = m_MAX_rrpv;  //long re-reference
        }
        else if (setType == FOLLOWER_SRRIP || setType == DEDICATED_SRRIP)
        {
            m_rrpv[setID][evictionWay] = m_MAX_rrpv-1;  //long re-reference
        }
        else 
        {
            assert(false && "Messed up DRRIP set classification\n");
        }

        return evictionWay;
    }
    else //execute DRRIP as normal
    {
        // check if there is a way with an invalid line
        for (int way = 0; way < m_numWays; ++way)
        {
            if (m_tagArray[setID][way] == -1)
            {
                if (setType == FOLLOWER_BRRIP || setType == DEDICATED_BRRIP)
                {
                    if (reduceNumMisses() % 32 == 0)
                        m_rrpv[setID][way] = m_MAX_rrpv-1; //long re-reference
                    else
                        m_rrpv[setID][way] = m_MAX_rrpv;  //distant re-reference
                }
                else if (setType == FOLLOWER_SRRIP || setType == DEDICATED_SRRIP)
                {
                    m_rrpv[setID][way] = m_MAX_rrpv-1;  //long re-reference
                }
                else
                {
                    assert(false && "Messed up DRRIP set classification\n");
                }
                return way;
            }
        }

        // no empty way have to kick out a way with distant rereference
        while (true)
        {
            //search for a way with distant re-reference
            for (int way = 0; way < m_numWays; ++way)
            {
                if (m_rrpv[setID][way] == m_MAX_rrpv)
                {
                    if (setType == FOLLOWER_BRRIP || setType == DEDICATED_BRRIP)
                    {
                        if (reduceNumMisses() % 32 == 0)
                            m_rrpv[setID][way] = m_MAX_rrpv-1; //long re-reference
                        else
                            m_rrpv[setID][way] = m_MAX_rrpv;  //long re-reference
                    }
                    else if (setType == FOLLOWER_SRRIP || setType == DEDICATED_SRRIP)
                    {
                        m_rrpv[setID][way] = m_MAX_rrpv-1;  //long re-reference
                    }
                    else
                    {
                        assert(false && "Messed up DRRIP set classification\n");
                    }
                    return way;
                }
            }

            //if we reached here then no line had distance reference.
            // age all rrpv
            for (int way = 0; way < m_numWays; ++way)
                ++m_rrpv[setID][way];
        }
    }
}

void LLC::updateReplacementState(int setID, int wayID)
{
    assert(m_rrpv[setID][wayID] <= m_MAX_rrpv); //checking for overflows
    m_rrpv[setID][wayID] = 0; //Hit Promotion policy
}

void LLC::moveToMRU(int setID, int wayID)
{
    ; // may be required for aggressive pinning of hubs
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
int LLC::findVtxID(intptr_t addr, int dTypeID)
{
    assert(m_dType_addrStart[dTypeID].size() == 1);
    assert(m_dType_addrEnd[dTypeID].size() == 1);
    assert(m_dType_elemSz[dTypeID].size() == 1);
    auto startAddr = m_dType_addrStart[dTypeID][0]; 
    auto endAddr   = m_dType_addrEnd[dTypeID][0]; 
    auto elemSz    = m_dType_elemSz[dTypeID][0];
    
    assert(addr >= startAddr && addr < endAddr);

    auto diff = addr - startAddr;
    int vtxID;
    if (dTypeID == FRONTIER)  //Frontiers are bitvectors
      vtxID = (diff) * 8; 
    else
      vtxID = (diff) / elemSz;
    return vtxID;
}

uint8_t LLC::findRereferenceVal(int irregInd, int regInd, int tid, int dTypeID)
{
  int epochSz = (m_graph.num_nodes() + m_numEpochs[dTypeID] - 1) / m_numEpochs[dTypeID];
  int numSubEpochs = m_numEpochs[dTypeID] / 2; 
  int subEpochSz = (epochSz + numSubEpochs - 1) / numSubEpochs;
  int VTX_PER_LINE;
  if (dTypeID == FRONTIER)
    VTX_PER_LINE = m_lineSz * 8; 
  else
    VTX_PER_LINE = m_lineSz / m_dType_elemSz[IRREGDATA][0];
  

  // get vertex/cache line number
  int cacheLineID = irregInd / VTX_PER_LINE;
  int currEpoch   = regInd / epochSz;
  uint8_t mask    = 1;
  uint8_t orMask  = (mask << 7);
  uint8_t andMask = (~orMask);
  assert(orMask == 128 && andMask == 127);
    
  if ((m_offsetMatrix[dTypeID][(currEpoch * m_numCacheLines[dTypeID]) + cacheLineID] & orMask) != 0)
  {
    //We dont have a reference this epoch - Find next Epoch access
    uint8_t reref = m_offsetMatrix[dTypeID][(currEpoch * m_numCacheLines[dTypeID]) + cacheLineID] & andMask;
    assert(reref <= MAX_REREF);
    return reref;
  }
  else
  {
    //We have a reference this epoch
    uint8_t lastRefQ = m_offsetMatrix[dTypeID][(currEpoch * m_numCacheLines[dTypeID]) + cacheLineID] & andMask;
    int delta = regInd - (currEpoch * epochSz);
    assert(delta >= 0);
    uint8_t subEpochDistQ = delta / subEpochSz;
    assert(subEpochDistQ <= MAX_REREF);
    if (subEpochDistQ <= lastRefQ)
    {
      //We are yet to be accessed within this epoch
      return 0;
    }
    else
    {
      //No further accesses this epoch; Look beyond to the next Epoch
      uint8_t nextRef = m_offsetMatrix[dTypeID][((currEpoch + 1) * m_numCacheLines[dTypeID]) + cacheLineID]; //this is why we store current & next Epoch information
      if ((nextRef & orMask) != 0)
      {
        uint8_t reref = nextRef & andMask; 
        return std::min(++reref, MAX_REREF); 
      }
      else
      {
        //Going to be accessed next epoch
        return 1;
      }
    }
  }
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

//DRRIP related functions
void LLC::incPSEL(int bankID)
{
  int localVal;
  do
  {
    localVal = m_policySelector[bankID];
    if (localVal >= m_pSEL_max)
      return;
  } while (!__sync_val_compare_and_swap(&m_policySelector[bankID], localVal, localVal + 1));
  return;
}

void LLC::decPSEL(int bankID)
{
  int localVal;
  do
  {
    localVal = m_policySelector[bankID];
    if (localVal <= m_pSEL_min)
      return;
  } while (!__sync_val_compare_and_swap(&m_policySelector[bankID], localVal, localVal - 1));
  return;
}

bool LLC::winningPolicy(int bankID)
{
  auto pSEL_val = m_policySelector[bankID];
  if (pSEL_val > m_pSEL_init)
    return m_BRRIP;
  else
    return m_SRRIP;
}

int LLC::determineSetType (int setID, int bankID)
{
  uint32_t constituencySz = m_numSets / 32;
  uint32_t constituency   = static_cast<uint32_t>(setID) / constituencySz;
  uint32_t offset         = static_cast<uint32_t>(setID) % constituencySz;
  assert(constituency < 32 && offset < constituencySz);
  if (constituency == offset)
  {
    incPSEL(bankID);
    return DEDICATED_SRRIP;
  }
  else if (constituency == ((~offset) % constituencySz))
  {
    decPSEL(bankID);
    return DEDICATED_BRRIP;
  }

  //We are dealing with a follower set
  if (winningPolicy(bankID) == m_SRRIP)
  {
    return FOLLOWER_SRRIP;
  }
  else
  {
    return FOLLOWER_BRRIP;
  }
}
