#include "l1.h"

L1::L1() { }

L1::~L1() 
{
    for (int set = 0; set < m_numSets; ++set)
    {
        delete[] m_tagArray[set];
        delete[] m_plru_bits[set];
    }
    delete[] m_tagArray;
    delete[] m_dirty;
    delete[] m_plru_bits;
    //delete[] m_setLocks;
}

void L1::Init()
{
    m_hits = 0;
    m_misses = 0;
    m_writebacks = 0;

    assert(m_numWays == 8); //PLRU implementation assumption

    m_tagArray = new intptr_t* [m_numSets];
    m_dirty    = new uint8_t [m_numSets];
    //m_setLocks = new PIN_LOCK [m_numSets];
    m_plru_bits  = new uint8_t* [m_numSets];
    //assert(m_tagArray != nullptr && m_dirty  != nullptr && m_setLocks != nullptr && m_plru_bits != nullptr);
    assert(m_tagArray != nullptr && m_dirty  != nullptr && m_plru_bits != nullptr);
    //PIN_InitLock(&m_setLock);
    for (int set = 0; set < m_numSets; ++set)
    {
        //PIN_InitLock(&m_setLocks[set]);
        m_tagArray[set] = new intptr_t [m_numWays];
        m_plru_bits[set]  = new uint8_t [m_numWays];
        assert(m_tagArray[set] != nullptr && m_plru_bits[set] != nullptr);
        m_dirty[set] = 0;
        for (int way = 0; way < m_numWays; ++way)
        {
            m_tagArray[set][way]  = -1; //the entire array is invalid initially
            m_plru_bits[set][way] = 0;
        }
    }
}

void L1::insertData(intptr_t addr, int tid, bool isWrite, bool &hit, bool &writeback, intptr_t &evictedAddr)
{
    int setID = findSet(addr);
    //PIN_GetLock(&m_setLocks[setID], tid+1);
    //PIN_GetLock(&m_setLock, tid+1); //Since LLC is non-inclusive, no need for locking
    hit = isCacheHit(addr, setID, isWrite);
    if (hit == false)
    {
        writeback = installNewLine(addr, setID, evictedAddr, isWrite);
    }
    //PIN_ReleaseLock(&m_setLock);
    //PIN_ReleaseLock(&m_setLocks[setID]);

    if (hit == true)
        ++m_hits;
    else
        ++m_misses;

    if (writeback == true)
    {
        ++m_writebacks;
        assert(evictedAddr != -1 * m_lineSz);
    }

    if (evictedAddr == -1 * m_lineSz)
        evictedAddr = -1;

    return;
}

void L1::backInvalidate(intptr_t addr, int tid)
{
    assert(false); //assuming zero inclusion in hierarchy
    #if 0
    int setID = findSet(addr);
    bool invalidated {false};
    assert(setID >= 0 && setID < m_numSets);
    //PIN_GetLock(&m_setLocks[setID], tid+1);
    //PIN_GetLock(&m_setLock, tid+1);
    /* NOTE: we dont need the above lock because only a private L2
    can backinvalidate a private L1; these operations cannot happen
    concurrently */
    intptr_t maskedAddr = addr / m_lineSz;
    for (int way = 0; way < m_numWays; ++way)
    {
        if (m_tagArray[setID][way] == maskedAddr)
        {
            moveToLRU(setID, way);
            m_tagArray[setID][way] = -1;
            //m_dirty[setID][way]    = 0;
            uint8_t mask   = 1 << way;
            m_dirty[setID] = m_dirty[setID] & (~mask);
            invalidated    = true;
            break;
        }
    }
    //PIN_ReleaseLock(&m_setLock);
    //PIN_ReleaseLock(&m_setLocks[setID]);
    if (invalidated == true)
        ++m_backInvalidates;
    #endif
}

stat_t L1::getHits()
{
    return m_hits;
}

stat_t L1::getMisses()
{
    return m_misses;
}

stat_t L1::getWritebacks()
{
    return m_writebacks;
}

stat_t L1::getBackInvalidates()
{
    return m_backInvalidates;
}

void L1::resetCacheStats()
{
    m_hits = 0;
    m_misses = 0;
    m_writebacks = 0;
    m_backInvalidates = 0;
}

bool L1::isCacheHit(intptr_t addr, int setID, bool isWrite)
{
    intptr_t maskedAddr = addr / m_lineSz;
    for (int way = 0; way < m_numWays; ++way)
    {
        if ((m_tagArray[setID][way] / m_lineSz) == maskedAddr)
        {
            updateReplacementState(setID, way);
            if (isWrite == true)
            {
                //m_dirty[setID][way] = 1;
                uint8_t mask = 1 << way;
                m_dirty[setID] = m_dirty[setID] | (mask);
            }
            return true;
        }
    }
    return false;
}

bool L1::installNewLine(intptr_t addr, int setID, intptr_t &evictedAddr, bool isWrite)
{
    int index = getReplacementIndex(setID);
    //evictedAddr = m_tagArray[setID][index] * m_lineSz; //line to be kicked out
    evictedAddr = m_tagArray[setID][index]; //line to be kicked out
    uint8_t mask = 1 << index;
    uint8_t retVal = m_dirty[setID] & mask; 
    assert(retVal == 0 || retVal == mask);
    if (evictedAddr == -1 * m_lineSz)
        assert(retVal == 0);
    //m_tagArray[setID][index] = addr / m_lineSz; //new line inserted 
    m_tagArray[setID][index] = addr; //new line inserted 
    //m_dirty[setID][index]    = (isWrite == true) ? 1 : 0;
    if (isWrite == true)
        m_dirty[setID] = m_dirty[setID] | (mask);
    else
        m_dirty[setID] = m_dirty[setID] & (~mask);
    return (retVal != 0);
}

int L1::findSet(intptr_t addr)
{
	intptr_t maskedAddr = addr / m_lineSz;
    int setIndex = maskedAddr % m_numSets;
    return setIndex;
}


//Tree-based PLRU implementation
int L1::getReplacementIndex(int setID)
{
    // check if there is a way with an invalid line
    for (int way = 0; way < m_numWays; ++way)
    {
        if (m_tagArray[setID][way] == -1)
        {
            updateReplacementState(setID, way);
            return way;
        }
    }

    // no empty way have to kick out the oldest way
    int retValue {-1};
    if (m_plru_bits[setID][0] == 0)
    {
       if (m_plru_bits[setID][1] == 0)
       {
          if (m_plru_bits[setID][2] == 0) retValue = 0;
          else                            retValue = 1;   // b2==1
       }
       else
       {                                                  // b1==1
          if (m_plru_bits[setID][3] == 0) retValue = 2;
          else                            retValue = 3;   // b3==1
       }
    }
    else
    {                                                     // b0==1
       if (m_plru_bits[setID][4] == 0)
       {
          if (m_plru_bits[setID][5] == 0) retValue = 4;
          else                            retValue = 5;   // b5==1
       }
       else
       {                                                  // b4==1
          if (m_plru_bits[setID][6] == 0) retValue = 6;
          else                            retValue = 7;   // b6==1
       }
    }
    assert(retValue != -1);
    updateReplacementState(setID, retValue);
    return retValue;
}

void L1::updateReplacementState(int setID, int wayID)
{
      if      (wayID==0) { m_plru_bits[setID][0]=1; m_plru_bits[setID][1]=1; m_plru_bits[setID][2]=1; }
      else if (wayID==1) { m_plru_bits[setID][0]=1; m_plru_bits[setID][1]=1; m_plru_bits[setID][2]=0; }
      else if (wayID==2) { m_plru_bits[setID][0]=1; m_plru_bits[setID][1]=0; m_plru_bits[setID][3]=1; }
      else if (wayID==3) { m_plru_bits[setID][0]=1; m_plru_bits[setID][1]=0; m_plru_bits[setID][3]=0; }
      else if (wayID==4) { m_plru_bits[setID][0]=0; m_plru_bits[setID][4]=1; m_plru_bits[setID][5]=1; }
      else if (wayID==5) { m_plru_bits[setID][0]=0; m_plru_bits[setID][4]=1; m_plru_bits[setID][5]=0; }
      else if (wayID==6) { m_plru_bits[setID][0]=0; m_plru_bits[setID][4]=0; m_plru_bits[setID][6]=1; }
      else if (wayID==7) { m_plru_bits[setID][0]=0; m_plru_bits[setID][4]=0; m_plru_bits[setID][6]=0; }
}

/* used during backinvalidates, could also be used for a
   aggressive thrash-resistence. 
   NOTE: update replacement seems to set the current way as far away 
   from the LRU position. This function sets the invalidated way at
   the LRU position
*/
void L1::moveToLRU (int setID, int wayID)
{
      if      (wayID==0) { m_plru_bits[setID][0]=0; m_plru_bits[setID][1]=0; m_plru_bits[setID][2]=0; }
      else if (wayID==1) { m_plru_bits[setID][0]=0; m_plru_bits[setID][1]=0; m_plru_bits[setID][2]=1; }
      else if (wayID==2) { m_plru_bits[setID][0]=0; m_plru_bits[setID][1]=1; m_plru_bits[setID][3]=0; }
      else if (wayID==3) { m_plru_bits[setID][0]=0; m_plru_bits[setID][1]=1; m_plru_bits[setID][3]=1; }
      else if (wayID==4) { m_plru_bits[setID][0]=1; m_plru_bits[setID][4]=0; m_plru_bits[setID][5]=0; }
      else if (wayID==5) { m_plru_bits[setID][0]=1; m_plru_bits[setID][4]=0; m_plru_bits[setID][5]=1; }
      else if (wayID==6) { m_plru_bits[setID][0]=1; m_plru_bits[setID][4]=1; m_plru_bits[setID][6]=0; }
      else if (wayID==7) { m_plru_bits[setID][0]=1; m_plru_bits[setID][4]=1; m_plru_bits[setID][6]=1; }
    
}

void L1::moveToMRU (int setID, int wayID)
{
    // could be useful for pinning-hubs
    ;
}

