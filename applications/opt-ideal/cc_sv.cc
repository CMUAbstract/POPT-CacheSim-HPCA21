// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <algorithm>
#include <cinttypes>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <omp.h>

#include "benchmark.h"
#include "bitmap.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "pvector.h"
#include "timer.h"


/*
GAP Benchmark Suite
Kernel: Connected Components (CC)
Author: Scott Beamer

Will return comp array labelling each vertex with a connected component ID

This CC implementation makes use of the Shiloach-Vishkin [2] algorithm with
implementation optimizations from Bader et al. [1]. Michael Sutton contributed
a fix for directed graphs using the min-max swap from [3], and it also produces
more consistent performance for undirected graphs.

[1] David A Bader, Guojing Cong, and John Feo. "On the architectural
    requirements for efficient execution of graph algorithms." International
    Conference on Parallel Processing, Jul 2005.

[2] Yossi Shiloach and Uzi Vishkin. "An o(logn) parallel connectivity algorithm"
    Journal of Algorithms, 3(1):57–67, 1982.

[3] Kishore Kothapalli, Jyothish Soman, and P. J. Narayanan. "Fast GPU
    algorithms for graph connectivity." Workshop on Large Scale Parallel
    Processing, 2010.
*/


using namespace std;

const int numDataTypes = 5;    
const int IRREGDATA   {0};
const int REGDATA     {1};
const int CSR_OFFSETS {2};
const int CSR_COORDS  {3};
const int FRONTIER    {4};
const int OTHERS      {5};



/* dummy functions used by pin */
void __attribute__((noinline)) __attribute__((optimize("O0"))) PIN_Start()
{
    //assert(false); //on successful replacement - this shouldnt be executed
    asm("");
}

void __attribute__((noinline)) __attribute__((optimize("O0"))) PIN_Stop()
{
    //assert(false);
    asm("");
}

void __attribute__((noinline)) __attribute__((optimize("O0"))) PIN_MapThreads(int* threadMapping, int numThreads)
{
    //assert(false);
    asm("");
}


void __attribute__((noinline)) __attribute__((optimize("O0"))) PIN_RegisterDataType(intptr_t addr, int dType, int numElements, size_t elemSz, int totalDataTypes)
{
    //assert(false);
    asm("");
}

void __attribute__((noinline)) __attribute__((optimize("O0"))) PIN_Init()
{
    //assert(false); //on successful replacement - this shouldnt be executed
    asm("");
}

void __attribute__((noinline)) __attribute__((optimize("O0"))) PIN_RegisterGraph(Graph& g, bool isPull)
{
    asm("");
}

void __attribute__((noinline)) __attribute__((optimize("O0"))) PIN_RegisterOffsetMatrix(uint8_t* offsetMatrix, int32_t numEpochs, int32_t numCachelines, int vtxPerLine, int dTypeID)
{
    asm("");
}

void __attribute__((noinline)) __attribute__((optimize("O0"))) PIN_UpdateRegIndex(int32_t index, int tid)
{
    asm("");
}

void __attribute__((noinline)) __attribute__((optimize("O0"))) PIN_UpdateIrregRanges(int32_t startID, int32_t endID)
{
    asm("");
}

void __attribute__((noinline)) __attribute__((optimize("O0"))) PIN_DumpStats()
{
    asm("");
}




// The hooking condition (comp_u < comp_v) may not coincide with the edge's
// direction, so we use a min-max swap such that lower component IDs propagate
// independent of the edge's direction.
pvector<NodeID> ShiloachVishkin(const Graph &g, int roiIter) 
{
  pvector<int> threadMapping (omp_get_max_threads());
  PIN_MapThreads(threadMapping.data(), omp_get_max_threads());
  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    threadMapping[tid] = tid;
  }

  pvector<NodeID> comp(g.num_nodes());
  PIN_RegisterDataType(reinterpret_cast<intptr_t>(comp.data()), IRREGDATA, g.num_nodes(), sizeof(NodeID), numDataTypes);

  /* Avoiding any parallelism before the ROI */
  //#pragma omp parallel for
  for (NodeID n=0; n < g.num_nodes(); n++)
    comp[n] = n;

  bool change = true;
  int num_iter = 0;

  while (change) 
  {
    change = false;

    if (roiIter == num_iter)
    {
      PIN_Start();
      PIN_UpdateIrregRanges(0, g.num_nodes()); 
    }

    #pragma omp parallel for schedule(dynamic, 64)
    for (NodeID u = 0; u < g.num_nodes(); u++) 
    {
      int tid = omp_get_thread_num();
      PIN_UpdateRegIndex(u, tid);

      for (NodeID v : g.out_neigh(u)) 
      {
        NodeID comp_u = comp[u];
        NodeID comp_v = comp[v];
        if (comp_u == comp_v) continue;
        // Hooking condition so lower component ID wins independent of direction
        NodeID high_comp = comp_u > comp_v ? comp_u : comp_v;
        NodeID low_comp = comp_u + (comp_v - high_comp);
        if (high_comp == comp[high_comp]) {
          change = true;
          comp[high_comp] = low_comp;
        }
      }
    }
    

    if (roiIter == num_iter)
    {
      PIN_UpdateIrregRanges(-1, -1);
      PIN_Stop();
  
      std::cout << "~~~ PINTOOL STATS BEGIN ~~~" << std::endl;
      PIN_DumpStats();
      std::cout << "~~~ PINTOOL STATS END ~~~" << std::endl;

      std::exit(0);
    }
    
    /* Avoiding parallelism outside the ROI */
    #pragma omp parallel for 
    for (NodeID n=0; n < g.num_nodes(); n++) {
      while (comp[n] != comp[comp[n]]) {
        comp[n] = comp[comp[n]];
      }
    }

    num_iter++;
  }
  cout << "Shiloach-Vishkin took " << num_iter << " iterations" << endl;
  return comp;
}


void PrintCompStats(const Graph &g, const pvector<NodeID> &comp) {
  cout << endl;
  unordered_map<NodeID, NodeID> count;
  for (NodeID comp_i : comp)
    count[comp_i] += 1;
  int k = 5;
  vector<pair<NodeID, NodeID>> count_vector;
  count_vector.reserve(count.size());
  for (auto kvp : count)
    count_vector.push_back(kvp);
  vector<pair<NodeID, NodeID>> top_k = TopK(count_vector, k);
  k = min(k, static_cast<int>(top_k.size()));
  cout << k << " biggest clusters" << endl;
  for (auto kvp : top_k)
    cout << kvp.second << ":" << kvp.first << endl;
  cout << "There are " << count.size() << " components" << endl;
}


// Verifies CC result by performing a BFS from a vertex in each component
// - Asserts search does not reach a vertex with a different component label
// - If the graph is directed, it performs the search as if it was undirected
// - Asserts every vertex is visited (degree-0 vertex should have own label)
bool CCVerifier(const Graph &g, const pvector<NodeID> &comp) {
  unordered_map<NodeID, NodeID> label_to_source;
  for (NodeID n : g.vertices())
    label_to_source[comp[n]] = n;
  Bitmap visited(g.num_nodes());
  visited.reset();
  vector<NodeID> frontier;
  frontier.reserve(g.num_nodes());
  for (auto label_source_pair : label_to_source) {
    NodeID curr_label = label_source_pair.first;
    NodeID source = label_source_pair.second;
    frontier.clear();
    frontier.push_back(source);
    visited.set_bit(source);
    for (auto it = frontier.begin(); it != frontier.end(); it++) {
      NodeID u = *it;
      for (NodeID v : g.out_neigh(u)) {
        if (comp[v] != curr_label)
          return false;
        if (!visited.get_bit(v)) {
          visited.set_bit(v);
          frontier.push_back(v);
        }
      }
      if (g.directed()) {
        for (NodeID v : g.in_neigh(u)) {
          if (comp[v] != curr_label)
            return false;
          if (!visited.get_bit(v)) {
            visited.set_bit(v);
            frontier.push_back(v);
          }
        }
      }
    }
  }
  for (NodeID n=0; n < g.num_nodes(); n++)
    if (!visited.get_bit(n))
      return false;
  return true;
}


int main(int argc, char* argv[]) {
  CLApp cli(argc, argv, "connected-components");
  if (!cli.ParseArgs())
    return -1;
  omp_set_num_threads(1); 
  PIN_Init();
  Builder b(cli);
  Graph g = b.MakeGraph();
  std::cout << "[GRAPH-STATS] Nodes = " << g.num_nodes() << std::endl;
  std::cout << "[GRAPH-STATS] Edges = " << g.num_edges_directed() << std::endl;
    
  auto offsetsPtr = g.returnOffsetsArrayForCSR();
  auto coordsPtr  = g.returnCoordsArrayForCSR();
  PIN_RegisterDataType(reinterpret_cast<intptr_t>(offsetsPtr), CSR_OFFSETS, g.num_nodes(), sizeof(NodeID*), numDataTypes);
  PIN_RegisterDataType(reinterpret_cast<intptr_t>(coordsPtr), CSR_COORDS, g.num_edges_directed(), sizeof(NodeID), numDataTypes);

  PIN_RegisterGraph(g, false);

  //BenchmarkKernel(cli, g, ShiloachVishkin, PrintCompStats, CCVerifier);
  ShiloachVishkin(g, cli.roiIter());
  return 0;
}
