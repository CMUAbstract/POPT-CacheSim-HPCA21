#ifndef CACHE_BACKEND_H
#define CACHE_BACKEND_H

#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <omp.h>
#include <cstring>
#include <cassert>
#include <cstring>
#include "l1.h"
#include "l2.h"
#include "llc.h"

class Cache
{
    public:
        Cache();
        ~Cache(); 
        
        void Init();
        void registerDataType(intptr_t addr, int dataTypeID, int numElements, size_t elemSz, int totalDataTypes);
        void insertData(intptr_t addr, int tid, bool isWrite, bool isInstrn = false); 
        void reportTotalStats();
		
    private:
        /////////////// DATA ////////////
        /* cache hierarchy */
        L1 m_l1d[numThreads]; //d-cache
        L1 m_l1i[numThreads]; //i-cache
        L2 m_l2[numThreads];
        LLC m_llc;
};

#endif //CACHE_BACKEND_H
