/*! @file
 *  This is an example of the PIN tool that demonstrates some basic PIN APIs 
 *  and could serve as the starting point for developing your first PIN tool
 */

#include "pin.H"
#include "cache_backend.h"
#include <iostream>
#include <fstream>

/* ================================================================== */
// Global variables 
/* ================================================================== */
Cache cache; 
bool simulateAccesses {false};
uint64_t numInsns {0};
int* m_tidMapAddrStart {nullptr};
int m_numThreads {-1};
int* m_threadMap {nullptr};
PIN_LOCK m_tidLock;

/* ===================================================================== */
// Command line switches
/* ===================================================================== */
KNOB<string> KnobOutputFile(KNOB_MODE_WRITEONCE,  "pintool",
    "o", "", "specify file name for MyPinTool output");

KNOB<BOOL>   KnobCount(KNOB_MODE_WRITEONCE,  "pintool",
    "count", "1", "count instructions, basic blocks and threads in the application");


/* ===================================================================== */
// Utilities
/* ===================================================================== */

/*!
 *  Print out help message.
 */
INT32 Usage()
{
    cerr << "This tool prints out the number of dynamically executed " << endl <<
            "instructions, basic blocks and threads in the application." << endl << endl;

    cerr << KNOB_BASE::StringKnobSummary() << endl;

    return -1;
}

/* ===================================================================== */
// Analysis routines
/* ===================================================================== */

VOID RecordIp(VOID* ip, THREADID tid)
{
    ++numInsns;
    //assert(false); //not doing I-Cache modeling to save sim time
                     //(Graph application code footprint is often small)
    #if 0
    if (simulateAccesses == true)
        cache.insertData(reinterpret_cast<intptr_t>(ip), static_cast<int>(tid), false, true); 
    #endif
}

VOID RecordMemAddrRead(VOID* addr, THREADID tid)
{
    if (simulateAccesses == true)
        cache.insertData(reinterpret_cast<intptr_t>(addr), static_cast<int>(tid), false);
    else
    {
        intptr_t addrStart = reinterpret_cast<intptr_t>(m_tidMapAddrStart);
        intptr_t addrEnd   = addrStart + (m_numThreads * sizeof(int)); 
        intptr_t tidAddr   = reinterpret_cast<intptr_t>(addr);
        assert(tidAddr < addrStart || tidAddr > addrEnd);
    }
}

VOID RecordMemAddrWrite(VOID* addr, THREADID tid)
{
    if (simulateAccesses == true)
        cache.insertData(reinterpret_cast<intptr_t>(addr), static_cast<int>(tid), true);
    else
    {
        intptr_t addrStart = reinterpret_cast<intptr_t>(m_tidMapAddrStart);
        intptr_t addrEnd   = addrStart + (m_numThreads * sizeof(int)); 
        intptr_t tidAddr   = reinterpret_cast<intptr_t>(addr);
        if (tidAddr >= addrStart && tidAddr < addrEnd)
        {
          int diff  = tidAddr - addrStart;
          int index = diff / sizeof(int);
          assert(index >= 0 && index < m_numThreads);
          m_threadMap[index] = static_cast<int>(tid);
          PIN_GetLock(&m_tidLock, tid + 1);
          std::cout << "[PIN] Remapped OS Tid-" << index << " to PIN Tid-" << tid << std::endl;
          PIN_ReleaseLock(&m_tidLock);
        }
    }
}

void replacementFunction(intptr_t addr, int dType, int numElements, size_t elemSz, int totalDataTypes)
{
    std::string names[] = {"irregData", "regData", "CSR-offsets", "CSR-coords", "frontier"};

    cache.registerDataType(addr, dType, numElements, elemSz, totalDataTypes);
    if (dType < numMainDataTypes)
    {
        std::cout << "[PIN-FUNC-REPLACEMENT] registered data type - " << names[dType] << std::endl;
    }
    else
    {
        assert(false);
    }
}

void mapThreads(int* threadMap, int numThreads)
{
    m_tidMapAddrStart = threadMap;
    m_numThreads      = numThreads;
    m_threadMap = new int [numThreads];
    for (int tid = 0; tid < m_numThreads; ++tid)
        m_threadMap[tid] = -1;
}


void registerGraph(Graph &g, bool isPull)
{
    cache.registerGraph(g, isPull);
}

void updateRegIndex(int32_t index, int tid)
{
    cache.updateRegIndex(index, m_threadMap[tid]);
}

void updateIrregRanges(int32_t startID, int32_t endID)
{
    cache.updateIrregRanges(startID, endID);
}


void startLogging()
{
    simulateAccesses = true;
}

void stopLogging()
{
    simulateAccesses = false;
}

void printStats()
{
    cache.reportTotalStats();
}

void initFunc()
{
    cache.Init();
}

/* ===================================================================== */
// Instrumentation callbacks
/* ===================================================================== */

void routineCallback(RTN rtn, void* v)
{
   //std::string rtn_name = RTN_Name(rtn).c_str();
   std::string rtn_name = RTN_Name(rtn);

   if (rtn_name.find("PIN_RegisterDataType") != std::string::npos)
   {
      RTN_Replace(rtn, AFUNPTR(replacementFunction));
   }
   
   if (rtn_name.find("PIN_Init") != std::string::npos)
   {
      RTN_Replace(rtn, AFUNPTR(initFunc));
   }

   if (rtn_name.find("PIN_Start") != std::string::npos)
   {
      RTN_Replace(rtn, AFUNPTR(startLogging));
   }

   if (rtn_name.find("PIN_Stop") != std::string::npos)
   {
      RTN_Replace(rtn, AFUNPTR(stopLogging));
   }
   if (rtn_name.find("PIN_RegisterGraph") != std::string::npos)
   {
      RTN_Replace(rtn, AFUNPTR(registerGraph));
   }
   if (rtn_name.find("PIN_MapThreads") != std::string::npos)
   {
      RTN_Replace(rtn, AFUNPTR(mapThreads));
   }
   if (rtn_name.find("PIN_UpdateRegIndex") != std::string::npos)
   {
      RTN_Replace(rtn, AFUNPTR(updateRegIndex));
   }
   if (rtn_name.find("PIN_UpdateIrregRanges") != std::string::npos)
   {
      RTN_Replace(rtn, AFUNPTR(updateIrregRanges));
   }

   if (rtn_name.find("PIN_DumpStats") != std::string::npos)
   {
      RTN_Replace(rtn, AFUNPTR(printStats));
   }
}

VOID Instruction(INS ins, VOID *v)
{
    if (!INS_IsSyscall(ins))
    {
        //pass instruction pointer to llc
        //#if 0
        INS_InsertPredicatedCall (
            ins, IPOINT_BEFORE, (AFUNPTR)RecordIp, 
            IARG_INST_PTR,
            IARG_THREAD_ID,
            IARG_END);
        //#endif

        //send addresses for memory instructions
        if (INS_IsMemoryRead(ins) || INS_IsMemoryWrite(ins))
        {
            for (unsigned int i = 0; i < INS_MemoryOperandCount(ins); ++i)
            {
                //pass load address 
                if (INS_MemoryOperandIsRead(ins, i))
                {
                    INS_InsertPredicatedCall (
                        ins, IPOINT_BEFORE, (AFUNPTR)RecordMemAddrRead, 
                        IARG_MEMORYOP_EA, i,
                        IARG_THREAD_ID,
                        IARG_END);
                }
                //pass store address 
                if (INS_MemoryOperandIsWritten(ins, i))
                {
                    INS_InsertPredicatedCall (
                        ins, IPOINT_BEFORE, (AFUNPTR)RecordMemAddrWrite, 
                        IARG_MEMORYOP_EA, i,
                        IARG_THREAD_ID,
                        IARG_END);
                }
            }
        }
    }
}

/*!
 * Increase counter of threads in the application.
 * This function is called for every thread created by the application when it is
 * about to start running (including the root thread).
 * @param[in]   threadIndex     ID assigned by PIN to the new thread
 * @param[in]   ctxt            initial register state for the new thread
 * @param[in]   flags           thread creation flags (OS specific)
 * @param[in]   v               value specified by the tool in the 
 *                              PIN_AddThreadStartFunction function call
 */
VOID ThreadStart(THREADID threadIndex, CONTEXT *ctxt, INT32 flags, VOID *v)
{
    assert(static_cast<int>(threadIndex) < numThreads);
}

/*!
 * Print out analysis results.
 * This function is called when the application exits.
 * @param[in]   code            exit code of the application
 * @param[in]   v               value specified by the tool in the 
 *                              PIN_AddFiniFunction function call
 */
VOID Fini(INT32 code, VOID *v)
{
    std::cout << "[PINTOOL] No. of Instructions = " << numInsns << std::endl;
    cache.reportTotalStats();
}

/*!
 * The main procedure of the tool.
 * This function is called when the application image is loaded but not yet started.
 * @param[in]   argc            total number of elements in the argv array
 * @param[in]   argv            array of command line arguments, 
 *                              including pin -t <toolname> -- ...
 */
int main(int argc, char *argv[])
{
    PIN_InitLock(&m_tidLock);

    PIN_InitSymbols();

    // Initialize PIN library. Print help message if -h(elp) is specified
    // in the command line or the command line is invalid 
    if( PIN_Init(argc,argv) )
    {
        return Usage();
    }
    
    #if 0
    string fileName = KnobOutputFile.Value();

    if (!fileName.empty()) { out = new std::ofstream(fileName.c_str());}
    #endif

    if (KnobCount)
    {
        // Register function to be called for function replacement
        RTN_AddInstrumentFunction(routineCallback, 0);

        // Register function to be called to instrument traces
        INS_AddInstrumentFunction(Instruction, 0);

        // Register function to be called for every thread before it starts running
        PIN_AddThreadStartFunction(ThreadStart, 0);

        // Register function to be called when the application exits
        PIN_AddFiniFunction(Fini, 0);
    }
    
    cerr <<  "===============================================" << endl;
    cerr <<  "This application is instrumented by MyPinTool" << endl;
    if (!KnobOutputFile.Value().empty()) 
    {
        cerr << "See file " << KnobOutputFile.Value() << " for analysis results" << endl;
    }
    cerr <<  "===============================================" << endl;

    // Start the program, never returns
    PIN_StartProgram();
    
    return 0;
}

/* ===================================================================== */
/* eof */
/* ===================================================================== */
