// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef BUILDER_H_
#define BUILDER_H_

#include <algorithm>
#include <parallel/algorithm>
#include <cinttypes>
#include <fstream>
#include <functional>
#include <type_traits>
#include <utility>
#include <set>
#include <omp.h>
#include <cassert>

#include "command_line.h"
#include "generator.h"
#include "graph.h"
#include "platform_atomics.h"
#include "pvector.h"
#include "reader.h"
#include "timer.h"
#include "util.h"
#include "sliding_queue.h"


/*
GAP Benchmark Suite
Class:  BuilderBase
Author: Scott Beamer

Given arguements from the command line (cli), returns a built graph
 - MakeGraph() will parse cli and obtain edgelist and call
   MakeGraphFromEL(edgelist) to perform actual graph construction
 - edgelist can be from file (reader) or synthetically generated (generator)
 - Common case: BuilderBase typedef'd (w/ params) to be Builder (benchmark.h)
*/


template <typename NodeID_, typename DestID_ = NodeID_,
          typename WeightT_ = NodeID_, bool invert = true>
class BuilderBase {
  typedef EdgePair<NodeID_, DestID_> Edge;
  typedef pvector<Edge> EdgeList;

  const CLBase &cli_;
  bool symmetrize_;
  bool needs_weights_;
  int64_t num_nodes_ = -1;

 public:
  explicit BuilderBase(const CLBase &cli) : cli_(cli) {
    symmetrize_ = cli_.symmetrize();
    needs_weights_ = !std::is_same<NodeID_, DestID_>::value;
  }

  DestID_ GetSource(EdgePair<NodeID_, NodeID_> e) {
    return e.u;
  }

  DestID_ GetSource(EdgePair<NodeID_, NodeWeight<NodeID_, WeightT_>> e) {
    return NodeWeight<NodeID_, WeightT_>(e.u, e.v.w);
  }

  NodeID_ FindMaxNodeID(const EdgeList &el) {
    NodeID_ max_seen = 0;
    #pragma omp parallel for reduction(max : max_seen)
    for (auto it = el.begin(); it < el.end(); it++) {
      Edge e = *it;
      max_seen = std::max(max_seen, e.u);
      max_seen = std::max(max_seen, (NodeID_) e.v);
    }
    return max_seen;
  }

  pvector<NodeID_> CountDegrees(const EdgeList &el, bool transpose) {
    pvector<NodeID_> degrees(num_nodes_, 0);
    #pragma omp parallel for
    for (auto it = el.begin(); it < el.end(); it++) {
      Edge e = *it;
      if (symmetrize_ || (!symmetrize_ && !transpose))
        fetch_and_add(degrees[e.u], 1);
      if (symmetrize_ || (!symmetrize_ && transpose))
        fetch_and_add(degrees[(NodeID_) e.v], 1);
    }
    return degrees;
  }

  static
  pvector<SGOffset> PrefixSum(const pvector<NodeID_> &degrees) {
    pvector<SGOffset> sums(degrees.size() + 1);
    SGOffset total = 0;
    for (size_t n=0; n < degrees.size(); n++) {
      sums[n] = total;
      total += degrees[n];
    }
    sums[degrees.size()] = total;
    return sums;
  }

  static
  pvector<SGOffset> ParallelPrefixSum(const pvector<NodeID_> &degrees) {
    const size_t block_size = 1<<20;
    const size_t num_blocks = (degrees.size() + block_size - 1) / block_size;
    pvector<SGOffset> local_sums(num_blocks);
    #pragma omp parallel for
    for (size_t block=0; block < num_blocks; block++) {
      SGOffset lsum = 0;
      size_t block_end = std::min((block + 1) * block_size, degrees.size());
      for (size_t i=block * block_size; i < block_end; i++)
        lsum += degrees[i];
      local_sums[block] = lsum;
    }
    pvector<SGOffset> bulk_prefix(num_blocks+1);
    SGOffset total = 0;
    for (size_t block=0; block < num_blocks; block++) {
      bulk_prefix[block] = total;
      total += local_sums[block];
    }
    bulk_prefix[num_blocks] = total;
    pvector<SGOffset> prefix(degrees.size() + 1);
    #pragma omp parallel for
    for (size_t block=0; block < num_blocks; block++) {
      SGOffset local_total = bulk_prefix[block];
      size_t block_end = std::min((block + 1) * block_size, degrees.size());
      for (size_t i=block * block_size; i < block_end; i++) {
        prefix[i] = local_total;
        local_total += degrees[i];
      }
    }
    prefix[degrees.size()] = bulk_prefix[num_blocks];
    return prefix;
  }

  // Removes self-loops and redundant edges
  // Side effect: neighbor IDs will be sorted
  void SquishCSR(const CSRGraph<NodeID_, DestID_, invert> &g, bool transpose,
                 DestID_*** sq_index, DestID_** sq_neighs) {
    pvector<NodeID_> diffs(g.num_nodes());
    DestID_ *n_start, *n_end;
    #pragma omp parallel for private(n_start, n_end)
    for (NodeID_ n=0; n < g.num_nodes(); n++) {
      if (transpose) {
        n_start = g.in_neigh(n).begin();
        n_end = g.in_neigh(n).end();
      } else {
        n_start = g.out_neigh(n).begin();
        n_end = g.out_neigh(n).end();
      }
      std::sort(n_start, n_end);
      DestID_ *new_end = std::unique(n_start, n_end);
      new_end = std::remove(n_start, new_end, n);
      diffs[n] = new_end - n_start;
    }
    pvector<SGOffset> sq_offsets = ParallelPrefixSum(diffs);
    *sq_neighs = new DestID_[sq_offsets[g.num_nodes()]];
    *sq_index = CSRGraph<NodeID_, DestID_>::GenIndex(sq_offsets, *sq_neighs);
    #pragma omp parallel for private(n_start)
    for (NodeID_ n=0; n < g.num_nodes(); n++) {
      if (transpose)
        n_start = g.in_neigh(n).begin();
      else
        n_start = g.out_neigh(n).begin();
      std::copy(n_start, n_start+diffs[n], (*sq_index)[n]);
    }
  }

  CSRGraph<NodeID_, DestID_, invert> SquishGraph(
      const CSRGraph<NodeID_, DestID_, invert> &g) {
    DestID_ **out_index, *out_neighs, **in_index, *in_neighs;
    SquishCSR(g, false, &out_index, &out_neighs);
    if (g.directed()) {
      if (invert)
        SquishCSR(g, true, &in_index, &in_neighs);
      return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), out_index,
                                                out_neighs, in_index,
                                                in_neighs);
    } else {
      return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), out_index,
                                                out_neighs);
    }
  }

  /*
  Graph Bulding Steps (for CSR):
    - Read edgelist once to determine vertex degrees (CountDegrees)
    - Determine vertex offsets by a prefix sum (ParallelPrefixSum)
    - Allocate storage and set points according to offsets (GenIndex)
    - Copy edges into storage
  */
  void MakeCSR(const EdgeList &el, bool transpose, DestID_*** index,
               DestID_** neighs) {
    pvector<NodeID_> degrees = CountDegrees(el, transpose);
    pvector<SGOffset> offsets = ParallelPrefixSum(degrees);
    *neighs = new DestID_[offsets[num_nodes_]];
    *index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, *neighs);
    #pragma omp parallel for
    for (auto it = el.begin(); it < el.end(); it++) {
      Edge e = *it;
      if (symmetrize_ || (!symmetrize_ && !transpose))
        (*neighs)[fetch_and_add(offsets[e.u], 1)] = e.v;
      if (symmetrize_ || (!symmetrize_ && transpose))
        (*neighs)[fetch_and_add(offsets[static_cast<NodeID_>(e.v)], 1)] =
            GetSource(e);
    }
  }

  CSRGraph<NodeID_, DestID_, invert> MakeGraphFromEL(EdgeList &el) {
    DestID_ **index = nullptr, **inv_index = nullptr;
    DestID_ *neighs = nullptr, *inv_neighs = nullptr;
    Timer t;
    t.Start();
    if (num_nodes_ == -1)
      num_nodes_ = FindMaxNodeID(el)+1;
    //#if 0 //TEMP
    if (needs_weights_)
      Generator<NodeID_, DestID_, WeightT_>::InsertWeights(el);
    //#endif
    MakeCSR(el, false, &index, &neighs);
    if (!symmetrize_ && invert)
      MakeCSR(el, true, &inv_index, &inv_neighs);
    t.Stop();
    PrintTime("Build Time", t.Seconds());
    if (symmetrize_)
      return CSRGraph<NodeID_, DestID_, invert>(num_nodes_, index, neighs);
    else
      return CSRGraph<NodeID_, DestID_, invert>(num_nodes_, index, neighs,
                                                inv_index, inv_neighs);
  }
  
  #if 0 //will complete the code later
  CSRGraph<NodeID_, DestID_, invert> relabelForSpatialLocality(
      const CSRGraph<NodeID_, DestID_, invert> &g) {
    if (g.directed()) {
        // Will add support soon
    }
    else {
        Timer t;
        t.start();
        
        /* STEP I: make a map between new and old vertex labels */

        long long counter = 0; //keep track of local counts
        std::map<NodeID_, int64_t> reMap[128]; //Conservatively assuming we will never use more than 128 threads

        /* relabel vertices in parallel (using local counter) */
        #pragma omp parallel for firstprivate(count)
        for (NodeID_ v = 0; v < g.num_nodes(); v++) {
            if (reMap[omp_get_thread_num()].find(v) == reMap.end()) {
                // vertex hasn't been labelled
                reMap.insert(std::pair<NodeID_, int64_t>(v, counter));
                counter++;
            }
            for (NodeID_ u : g.in_neigh(v)) {
                if (reMap[omp_get_thread_num()].find(u) == reMap.end()) {
                    // vertex hasn't been labelled
                    reMap.insert(std::pair<NodeID_, int64_t>(u, counter));
                    counter++;
                }
            }
        }

        /* Update counts based on maximum count for each thread */
        int64_t offset = 0;
        for (int i = 0; i < 128; i++) {
            if (reMap[i].size() != 0) {
                // adding offset to all counts of current map
                std::map<NodeID_, int64_t>::iterator it, it_end;
                #pragma omp parallel for 
                for (it = reMap[i].begin(), it_end = reMap[i].end(); it != it_end; it++) {
                    it->second += offset; 
                }
                
                // finding maximum value of current set 
                int64_t maxVal = 0;
                #pragma omp parallel for reduction(max: maxVal)
                for (it = reMap[i].begin(), it_end = reMap[i].end(); it != it_end; it++) {
                    if (it->second > maxVal) {
                        maxVal = it->second;
                    }
                }
                offset = maxVal;
            }
        }
        
        /* Merge local containers */
        std::map <NodeID_, int64_t> merged_reMap; 
        for (int i = 0; i < 128; i++) {
            if (reMap[i].size() != 0) {
                merged_reMap.insert(reMap[i].begin(), reMap[i].end());
            }
        }

        /* STEP II: rewrite CSR based on this reMap */
        DestID_* neighs = new DestID_[2 * g.num_edges()];
        DestID_** index = CSRGraph<NodeID_, DestID_>::relabelIndex(offsets, neighs, reMap);


    }
  }
  #endif

  CSRGraph<NodeID_, DestID_, invert> MakeGraph() {
    CSRGraph<NodeID_, DestID_, invert> g;
    {  // extra scope to trigger earlier deletion of el (save memory)
      EdgeList el;
      if (cli_.filename() != "") {
        Reader<NodeID_, DestID_, WeightT_, invert> r(cli_.filename());
        if ((r.GetSuffix() == ".sg") || (r.GetSuffix() == ".wsg")) {
          return r.ReadSerializedGraph();
        } else {
          el = r.ReadFile(needs_weights_);
        }
      } else if (cli_.scale() != -1) {
        Generator<NodeID_, DestID_> gen(cli_.scale(), cli_.degree());
        el = gen.GenerateEL(cli_.uniform());
      }
      g = MakeGraphFromEL(el);
    }
    #if 0
    if (cli_.relabel() == 1) {
        g_new = relabelForSpatialLocality(g); 
    }
    #endif
    return SquishGraph(g);
  }

  // Relabels (and rebuilds) graph by order of decreasing degree
  static
  CSRGraph<NodeID_, DestID_, invert> RelabelByDegree(
      const CSRGraph<NodeID_, DestID_, invert> &g) {
    if (g.directed()) {
      std::cout << "Cannot relabel directed graph" << std::endl;
      std::exit(-11);
    }
    Timer t;
    t.Start();
    typedef std::pair<int64_t, NodeID_> degree_node_p;
    pvector<degree_node_p> degree_id_pairs(g.num_nodes());
    #pragma omp parallel for
    for (NodeID_ n=0; n < g.num_nodes(); n++)
      degree_id_pairs[n] = std::make_pair(g.out_degree(n), n);
    std::sort(degree_id_pairs.begin(), degree_id_pairs.end(),
              std::greater<degree_node_p>());
    pvector<NodeID_> degrees(g.num_nodes());
    pvector<NodeID_> new_ids(g.num_nodes());
    #pragma omp parallel for
    for (NodeID_ n=0; n < g.num_nodes(); n++) {
      degrees[n] = degree_id_pairs[n].first;
      new_ids[degree_id_pairs[n].second] = n;
    }
    pvector<SGOffset> offsets = ParallelPrefixSum(degrees);
    DestID_* neighs = new DestID_[offsets[g.num_nodes()]];
    DestID_** index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, neighs);
    #pragma omp parallel for schedule (dynamic, 1024)
    for (NodeID_ u=0; u < g.num_nodes(); u++) {
      for (NodeID_ v : g.out_neigh(u))
        neighs[offsets[new_ids[u]]++] = new_ids[v];
      std::sort(index[new_ids[u]], index[new_ids[u]+1]);
    }
    t.Stop();
    PrintTime("Relabel", t.Seconds());
    return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), index, neighs);
  }

  /* 
    CSR-segmenting as proposed in the Cagra paper

    Partitions the graphs and produces a sub-graph within a specified range of vertices. 

    The following implementation assumes use in pull implementation and that only the 
    partitioned CSC is required
  */ 
  static
  CSRGraph<NodeID_, DestID_, invert> graphSlicer(
      const CSRGraph<NodeID_, DestID_, invert> &g, NodeID_ startID, NodeID_ stopID, 
      bool outDegree = false,  bool modifyBothDestlists = false) { 
      /* create a partition of a graph in the range [startID, stopID) */
      Timer t;
      t.Start();

      //NOTE: that pull implementation should specify outDegree == false
      //      and push implementations should use outDegree == true

      if (g.directed() == true) {
        /* Step I : For the requested range [startID, stopID), construct the reduced degree per vertex */
        pvector<NodeID_> degrees(g.num_nodes()); //note that stopID is not included in the range
        #pragma omp parallel for schedule(dynamic,1024)
        for (NodeID_ n = 0; n < g.num_nodes(); ++n) {
          if (outDegree == true) {
            NodeID_ newDegree(0);
            for (NodeID_ m : g.out_neigh(n)) {
              if(m >= startID && m < stopID) {
                ++newDegree;
              }
            }
            degrees[n] = newDegree;
          }
          else { 
            NodeID_ newDegree(0);
            for (NodeID_ m : g.in_neigh(n)) {
              if(m >= startID && m < stopID) {
                ++newDegree;
              }
            }
            degrees[n] = newDegree;
          }
        }

        /* Step II : Construct a trimmed offset list */
        pvector<SGOffset> offsets = ParallelPrefixSum(degrees);
        DestID_* neighs = new DestID_[offsets[g.num_nodes()]];
        DestID_** index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, neighs);
        if (outDegree == true) {
          #pragma omp parallel for schedule(dynamic, 1024) 
          for (NodeID_ u = 0; u < g.num_nodes(); ++u) {
            for (NodeID_ v : g.out_neigh(u)) {
              if (v >= startID && v < stopID) {
                neighs[offsets[u]++] = v;
              }
            }
          }
        }
        else {
          #pragma omp parallel for schedule(dynamic, 1024) 
          for (NodeID_ u = 0; u < g.num_nodes(); ++u) {
            for (NodeID_ v : g.in_neigh(u)) {
              if (v >= startID && v < stopID) {
                neighs[offsets[u]++] = v;
              }
            }
          }
        }
        
        /* Step III : Populate the inv dest lists (for push-pull implementations) */
        DestID_* inv_neighs(nullptr);
        DestID_** inv_index(nullptr);
        if (modifyBothDestlists == true) {
          //allocate space
          pvector<NodeID_> inv_degrees(g.num_nodes());
          #pragma omp parallel for 
          for (NodeID_ u = 0; u < g.num_nodes(); ++u) {
            if (outDegree == true) {
              inv_degrees[u] = g.in_degree(u);
            }
            else {
              inv_degrees[u] = g.out_degree(u);
            }
          }
          pvector<SGOffset> inv_offsets = ParallelPrefixSum(inv_degrees);
          inv_neighs = new DestID_[inv_offsets[g.num_nodes()]];
          inv_index = CSRGraph<NodeID_, DestID_>::GenIndex(inv_offsets, inv_neighs);

          //populate the inv dest list
          #pragma omp parallel for schedule(dynamic, 1024) 
          for (NodeID_ u = 0; u < g.num_nodes(); ++u) {
            if (outDegree == true) {
              for (NodeID_ v : g.in_neigh(u)) {
                inv_neighs[inv_offsets[u]++] = v;
              }
            }
            else {
              for (NodeID_ v : g.out_neigh(u)) {
                inv_neighs[inv_offsets[u]++] = v;
              }
            }
          }
        }
        #if 0
        else {
          inv_neighs = neighs;
          inv_index  = index;
        }
        #endif

        /* Step IV : return the appropriate graph */
        if (outDegree == true) {
          t.Stop();
          PrintTime("Slice-time", t.Seconds());
          return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), index, neighs, inv_index, inv_neighs);
        }
        else {
          t.Stop();
          PrintTime("Slice-time", t.Seconds());
          return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), inv_index, inv_neighs, index, neighs);
        }
      }
      else {
        /* Step I : For the requested range [startID, stopID), construct the reduced degree per vertex */
        pvector<NodeID_> degrees(g.num_nodes()); //note that stopID is not included in the range
        #pragma omp parallel for schedule(dynamic,1024)
        for (NodeID_ n = 0; n < g.num_nodes(); ++n) {
          NodeID_ newDegree(0);
          for (NodeID_ m : g.out_neigh(n)) {
            if(m >= startID && m < stopID) {
              ++newDegree; //if neighbor is in current partition
            }
          }
          degrees[n] = newDegree;
        }

        /* Step II : Construct a trimmed offset list */
        pvector<SGOffset> offsets = ParallelPrefixSum(degrees);
        DestID_* neighs = new DestID_[offsets[g.num_nodes()]];
        DestID_** index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, neighs);
        #pragma omp parallel for schedule(dynamic, 1024) 
        for (NodeID_ u = 0; u < g.num_nodes(); ++u) {
          for (NodeID_ v : g.out_neigh(u)) {
            if (v >= startID && v < stopID) {
              neighs[offsets[u]++] = v;
            }
          }
        }
        
        /* Step III : return the appropriate graph */
        t.Stop();
        PrintTime("Slice-time", t.Seconds());
        return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), index, neighs);
      }
    }
    
    static
    CSRGraph<NodeID_, DestID_, invert> quantizeGraph(
        const CSRGraph<NodeID_, DestID_, invert> &g, NodeID_ numTiles) {

        NodeID_ tileSz = g.num_nodes() / numTiles;
        if (numTiles > g.num_nodes())
            tileSz = 1;
        else if (g.num_nodes() % numTiles != 0)
            tileSz += 1;
        
        pvector<NodeID_> degrees(g.num_nodes(), 0);
        #pragma omp parallel for
        for (NodeID_ n = 0; n < g.num_nodes(); ++n)
        {
          std::set<NodeID_> uniqNghs;
          for (NodeID_ ngh : g.out_neigh(n))
          {
            uniqNghs.insert(ngh / tileSz);
          }
          degrees[n] = uniqNghs.size();
          assert(degrees[n] <= numTiles);
        }
        
        pvector<SGOffset> offsets = ParallelPrefixSum(degrees);
        DestID_* neighs = new DestID_[offsets[g.num_nodes()]];
        DestID_** index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, neighs);
        #pragma omp parallel for schedule (dynamic, 1024)
        for (NodeID_ u=0; u < g.num_nodes(); u++) {
          std::set<NodeID_> uniqNghs;
          for (NodeID_ ngh : g.out_neigh(u))
          {
            uniqNghs.insert(ngh / tileSz);
          }
          
          auto it = uniqNghs.begin();
          for (NodeID_ i = 0; i < static_cast<NodeID_>(uniqNghs.size()); ++i)
          {
            neighs[offsets[u]++] = *it;
            it++;
          }
        }
        return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), index, neighs);
    }
  
  // Proper degree sorting
  static
  CSRGraph<NodeID_, DestID_, invert> DegSort (
      const CSRGraph<NodeID_, DestID_, invert> &g, bool outDegree, pvector<NodeID_> &new_ids, bool createOnlyDegList, bool createBothCSRs) {
      Timer t;
      t.Start();
       
      typedef std::pair<int64_t, NodeID_> degree_node_p;
      pvector<degree_node_p> degree_id_pairs(g.num_nodes());
      if (g.directed() == true) {
        /* Step I: Create a list of degrees */
        #pragma omp parallel for
        for (NodeID_ n=0; n < g.num_nodes(); n++) { 
          if (outDegree == true) {
            degree_id_pairs[n] = std::make_pair(g.out_degree(n), n);
          }
          else {
            degree_id_pairs[n] = std::make_pair(g.in_degree(n), n);
          }
        }

        /* Step II: Sort based on degree order */
        __gnu_parallel::sort(degree_id_pairs.begin(), degree_id_pairs.end(),
                             std::greater<degree_node_p>()); //TODO:Use parallel sort

        /* Step III: assigned remap for the hub vertices */
        #pragma omp parallel for 
        for (NodeID_ n=0; n < g.num_nodes(); n++) {
          new_ids[degree_id_pairs[n].second] = n;
        }

        /* Step VI: generate degree to build a new graph */
        pvector<NodeID_> degrees(g.num_nodes());
        pvector<NodeID_> inv_degrees(g.num_nodes());
        if (outDegree == true) {
          #pragma omp parallel for 
          for (NodeID_ n=0; n < g.num_nodes(); n++) {
            degrees[new_ids[n]]     = g.out_degree(n); 
            inv_degrees[new_ids[n]] = g.in_degree(n);
          }
        }
        else {
          #pragma omp parallel for 
          for (NodeID_ n=0; n < g.num_nodes(); n++) {
            degrees[new_ids[n]]     = g.in_degree(n); 
            inv_degrees[new_ids[n]] = g.out_degree(n);
          }
        }
        
        /* Graph building phase */
		pvector<SGOffset> offsets = ParallelPrefixSum(inv_degrees);
        DestID_* neighs = new DestID_[offsets[g.num_nodes()]];
        DestID_** index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, neighs);
        #pragma omp parallel for schedule (dynamic, 1024)
        for (NodeID_ u=0; u < g.num_nodes(); u++) {
          if (outDegree == true) {
            for (NodeID_ v : g.in_neigh(u))
              neighs[offsets[new_ids[u]]++] = new_ids[v];
          }
          else {
            for (NodeID_ v : g.out_neigh(u))
              neighs[offsets[new_ids[u]]++] = new_ids[v];
          }
          std::sort(index[new_ids[u]], index[new_ids[u]+1]); //sort neighbors of each vertex
        }
        DestID_* inv_neighs(nullptr);
        DestID_** inv_index(nullptr);
        if (createOnlyDegList == true || createBothCSRs == true) {
          // making the inverse list (in-degrees in this case)
          pvector<SGOffset> inv_offsets = ParallelPrefixSum(degrees);
          inv_neighs = new DestID_[inv_offsets[g.num_nodes()]]; inv_index = CSRGraph<NodeID_, DestID_>::GenIndex(inv_offsets, inv_neighs);
          if (createBothCSRs == true) {
            #pragma omp parallel for schedule(dynamic, 1024)
            for (NodeID_ u=0; u < g.num_nodes(); u++) {
              if (outDegree == true) {
                for (NodeID_ v : g.out_neigh(u))
                  inv_neighs[inv_offsets[new_ids[u]]++] = new_ids[v];
              }
              else {
                for (NodeID_ v : g.in_neigh(u))
                  inv_neighs[inv_offsets[new_ids[u]]++] = new_ids[v];
              }
              std::sort(inv_index[new_ids[u]], inv_index[new_ids[u]+1]); //sort neighbors of each vertex
            }
          }
        }
        t.Stop();
        PrintTime("DegSort time", t.Seconds());
        if (outDegree == true) {
          return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), inv_index, inv_neighs, index, neighs);
        }
        else {
          return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), index, neighs, inv_index, inv_neighs);
        }
      }
      else {
        /* Undirected graphs - no need to make separate lists for in and out degree */
        /* Step I: Create a list of degrees */
        #pragma omp parallel for 
        for (NodeID_ n=0; n < g.num_nodes(); n++) { 
          degree_id_pairs[n] = std::make_pair(g.out_degree(n), n);
        }

        /* Step II: Sort based on degree order */
        __gnu_parallel::sort(degree_id_pairs.begin(), degree_id_pairs.end(),
                             std::greater<degree_node_p>()); //TODO:Use parallel sort

        /* Step III: assigned remap for the hub vertices */
        #pragma omp parallel for 
        for (NodeID_ n=0; n < g.num_nodes(); n++) {
          new_ids[degree_id_pairs[n].second] = n;
        }

        /* Step VI: generate degree to build a new graph */
        pvector<NodeID_> degrees(g.num_nodes());
        #pragma omp parallel for 
        for (NodeID_ n=0; n < g.num_nodes(); n++) {
          degrees[new_ids[n]] = g.out_degree(n); 
        }
        
        /* Graph building phase */
        pvector<SGOffset> offsets = ParallelPrefixSum(degrees);
        DestID_* neighs = new DestID_[offsets[g.num_nodes()]];
        DestID_** index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, neighs);
        #pragma omp parallel for schedule (dynamic, 1024)
        for (NodeID_ u=0; u < g.num_nodes(); u++) {
          for (NodeID_ v : g.out_neigh(u))
            neighs[offsets[new_ids[u]]++] = new_ids[v];
          std::sort(index[new_ids[u]], index[new_ids[u]+1]);
        }
        t.Stop();
        PrintTime("DegSort time", t.Seconds());
        return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), index, neighs);
      }
    }
  
  static
  CSRGraph<NodeID_, DestID_, invert> RandOrder (
      const CSRGraph<NodeID_, DestID_, invert> &g, pvector<NodeID_> &new_ids, bool createOnlyDegList, bool createBothCSRs) {
      Timer t;
      t.Start();
      std::srand(0); //so that the random graph generated is the same everytime
      bool outDegree = true;
       
      if (g.directed() == true) {
        //Step I: create a random permutation - SLOW implementation
        pvector<NodeID_> claimedVtxs(g.num_nodes(), 0);
        
        //#pragma omp parallel for
        for (NodeID_ v = 0; v < g.num_nodes(); ++v)
        {
          while (true) 
          {
            NodeID_ randID = std::rand() % g.num_nodes();
            if (claimedVtxs[randID] != 1)
            {
              if (compare_and_swap(claimedVtxs[randID], 0, 1) == true)
              {
                new_ids[v] = randID;
                break;
              }
              else
                continue;
            }
          }
        }

        #pragma omp parallel for
        for (NodeID_ v = 0; v < g.num_nodes(); ++v)
            assert(new_ids[v] != -1);

        /* Step VI: generate degree to build a new graph */
        pvector<NodeID_> degrees(g.num_nodes());
        pvector<NodeID_> inv_degrees(g.num_nodes());
        if (outDegree == true) {
          #pragma omp parallel for 
          for (NodeID_ n=0; n < g.num_nodes(); n++) {
            degrees[new_ids[n]]     = g.out_degree(n); 
            inv_degrees[new_ids[n]] = g.in_degree(n);
          }
        }
        else {
          #pragma omp parallel for 
          for (NodeID_ n=0; n < g.num_nodes(); n++) {
            degrees[new_ids[n]]     = g.in_degree(n); 
            inv_degrees[new_ids[n]] = g.out_degree(n);
          }
        }
        
        /* Graph building phase */
		pvector<SGOffset> offsets = ParallelPrefixSum(inv_degrees);
        DestID_* neighs = new DestID_[offsets[g.num_nodes()]];
        DestID_** index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, neighs);
        #pragma omp parallel for schedule (dynamic, 1024)
        for (NodeID_ u=0; u < g.num_nodes(); u++) {
          if (outDegree == true) {
            for (NodeID_ v : g.in_neigh(u))
              neighs[offsets[new_ids[u]]++] = new_ids[v];
          }
          else {
            for (NodeID_ v : g.out_neigh(u))
              neighs[offsets[new_ids[u]]++] = new_ids[v];
          }
          std::sort(index[new_ids[u]], index[new_ids[u]+1]); //sort neighbors of each vertex
        }
        DestID_* inv_neighs(nullptr);
        DestID_** inv_index(nullptr);
        if (createOnlyDegList == true || createBothCSRs == true) {
          // making the inverse list (in-degrees in this case)
          pvector<SGOffset> inv_offsets = ParallelPrefixSum(degrees);
          inv_neighs = new DestID_[inv_offsets[g.num_nodes()]];
          inv_index = CSRGraph<NodeID_, DestID_>::GenIndex(inv_offsets, inv_neighs);
          if (createBothCSRs == true) {
            #pragma omp parallel for schedule(dynamic, 1024)
            for (NodeID_ u=0; u < g.num_nodes(); u++) {
              if (outDegree == true) {
                for (NodeID_ v : g.out_neigh(u))
                  inv_neighs[inv_offsets[new_ids[u]]++] = new_ids[v];
              }
              else {
                for (NodeID_ v : g.in_neigh(u))
                  inv_neighs[inv_offsets[new_ids[u]]++] = new_ids[v];
              }
              std::sort(inv_index[new_ids[u]], inv_index[new_ids[u]+1]); //sort neighbors of each vertex
            }
          }
        }
        t.Stop();
        PrintTime("RandOrder time", t.Seconds());
        if (outDegree == true) {
          return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), inv_index, inv_neighs, index, neighs);
        }
        else {
          return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), index, neighs, inv_index, inv_neighs);
        }
      }
      else {
        /* Undirected graphs - no need to make separate lists for in and out degree */
        //Step I: create a random permutation - SLOW implementation
        pvector<NodeID_> claimedVtxs(g.num_nodes(), 0);
        
        //#pragma omp parallel for 
        for (NodeID_ v = 0; v < g.num_nodes(); ++v)
        {
          while (true) 
          {
            NodeID_ randID = std::rand() % g.num_nodes();
            if (claimedVtxs[randID] != 1)
            {
              if (compare_and_swap(claimedVtxs[randID], 0, 1) == true)
              {
                new_ids[v] = randID;
                break;
              }
              else
                continue;
            }
          }
        }
        
        #pragma omp parallel for
        for (NodeID_ v = 0; v < g.num_nodes(); ++v)
            assert(new_ids[v] != -1);

        /* Step VI: generate degree to build a new graph */
        pvector<NodeID_> degrees(g.num_nodes());
        #pragma omp parallel for 
        for (NodeID_ n=0; n < g.num_nodes(); n++) {
          degrees[new_ids[n]] = g.out_degree(n); 
        }
        
        /* Graph building phase */
        pvector<SGOffset> offsets = ParallelPrefixSum(degrees);
        DestID_* neighs = new DestID_[offsets[g.num_nodes()]];
        DestID_** index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, neighs);
        #pragma omp parallel for schedule (dynamic, 1024)
        for (NodeID_ u=0; u < g.num_nodes(); u++) {
          for (NodeID_ v : g.out_neigh(u))
            neighs[offsets[new_ids[u]]++] = new_ids[v];
          std::sort(index[new_ids[u]], index[new_ids[u]+1]);
        }
        t.Stop();
        PrintTime("RandOrder time", t.Seconds());
        return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), index, neighs);
      }
    }
  
  /*
    Return a compressed transpose matrix (Rereference Matrix)
  */
  static void makeOffsetMatrix(const CSRGraph<NodeID_, DestID_, invert> &g, 
                                      pvector<uint8_t> &offsetMatrix, int numVtxPerLine, int numEpochs, bool traverseCSR = true)
    {
      if (g.directed() == false)
        traverseCSR = true;

      Timer tm;
      #if 0
      /* Step 0: Sanity check that neighborhoods are sorted */
      #pragma omp parallel for schedule(dynamic, 64)
      for (NodeID_ v = 0; v < g.num_nodes(); ++v)
      {
        NodeID_ maxNgh {-1};
        for (NodeID ngh: g.out_neigh(v))
        {
          assert(ngh > maxNgh);
          maxNgh = ngh;
        }
      }
      #endif

      /* Step I: Collect quantized edges & Compact vertices into "super vertices" */ 
      tm.Start();
      NodeID_ numCacheLines = (g.num_nodes() + numVtxPerLine - 1) / numVtxPerLine;
      NodeID_ epochSz       = (g.num_nodes() + numEpochs - 1) / numEpochs;
      pvector<NodeID_> lastRef(numCacheLines * numEpochs, -1);
      NodeID_ chunkSz = 64 / numVtxPerLine;
      if (chunkSz == 0)
        chunkSz = 1;
        
      #pragma omp parallel for schedule(dynamic, chunkSz) 
      for (NodeID_ c = 0; c < numCacheLines; ++c)
      {
        NodeID_ startVtx = c * numVtxPerLine;
        NodeID_ endVtx   = (c+1) * numVtxPerLine;
        if (c == numCacheLines - 1)
          endVtx = g.num_nodes();

        for (NodeID_ v = startVtx; v < endVtx; ++v)
        {
          if (traverseCSR == true)
          {
            for (NodeID_ ngh : g.out_neigh(v))
            {
              NodeID_ nghEpoch = ngh / epochSz; 
              lastRef[(c * numEpochs) + nghEpoch] = std::max(ngh, lastRef[(c * numEpochs) + nghEpoch]);
            }
          }
          else
          {
            for (NodeID_ ngh : g.in_neigh(v))
            {
              NodeID_ nghEpoch = ngh / epochSz; 
              lastRef[(c * numEpochs) + nghEpoch] = std::max(ngh, lastRef[(c * numEpochs) + nghEpoch]);
            }
          }
        }
      }
      tm.Stop();
      std::cout << "[CSR-HYBRID-PREPROCESSING] Time to quantize nghs and compact vertices = " << tm.Seconds() << std::endl;
      assert(numEpochs == 256);

      /* Step II: Converting adjacency matrix into offsets */
      tm.Start();
      uint8_t maxReref = 127; //because MSB is reserved for identifying between reref val (1) & switch point (0)
      NodeID_ subEpochSz = (epochSz + 127) / 128; //Using remaining 7 bits to identify intra-epoch information
      pvector<uint8_t> compressedOffsets(numCacheLines * numEpochs);
      uint8_t mask    = 1;
      uint8_t orMask  = mask << 7;
      uint8_t andMask = ~(orMask);
      assert(orMask == 128 && andMask == 127);
      #pragma omp parallel for schedule (static)
      for (NodeID_ c = 0; c < numCacheLines; ++c)
      {
        { // first set values for the last epoch 
          NodeID_ e = numEpochs - 1;
          if (lastRef[(c * numEpochs) + e] != -1)
          {
            compressedOffsets[(c * numEpochs) + e] = maxReref; 
            compressedOffsets[(c * numEpochs) + e] &= andMask; 
          }
          else
          {
            compressedOffsets[(c * numEpochs) + e] = maxReref; 
            compressedOffsets[(c * numEpochs) + e] |= orMask; 
          }
        }
        
        //Now back track and set values for all epochs
        for (NodeID_ e = numEpochs - 2; e >= 0; --e)
        {
          if (lastRef[(c * numEpochs) + e] != -1)
          {
            // There was a ref this epoch - store the quantized val of the lastRef 
            NodeID_ subEpochDist = lastRef[(c * numEpochs) + e] - (e * epochSz);
            assert(subEpochDist >= 0);
            NodeID_ lastRefQ = (subEpochDist / subEpochSz);
            assert(lastRefQ <= maxReref);
            compressedOffsets[(c * numEpochs) + e] = static_cast<uint8_t>(lastRefQ);
            compressedOffsets[(c * numEpochs) + e] &= andMask; 
          }
          else
          {
            if ((compressedOffsets[(c * numEpochs) + e + 1] & orMask) != 0)
            {
              //No access next epoch as well - add inter-epoch distance
              uint8_t nextRef = compressedOffsets[(c * numEpochs) + e + 1] & andMask;
              if (nextRef == maxReref)
                compressedOffsets[(c * numEpochs) + e] = maxReref;
              else
                compressedOffsets[(c * numEpochs) + e] = nextRef + 1;
            }
            else
            {
              //There is an access next epoch - so inter-epoch distance is set to next epoch
              compressedOffsets[(c * numEpochs) + e] = 1; 
            }
            compressedOffsets[(c * numEpochs) + e] |= orMask;
          }
        }
      }
      tm.Stop();
      std::cout << "[CSR-HYBRID-PREPROCESSING] Time to convert to offsets matrix = " << tm.Seconds() << std::endl;

      /* Step III: Transpose edgePresent*/  
      tm.Start();
      #pragma omp parallel for schedule (static)
      for (NodeID_ c = 0; c < numCacheLines; ++c)
      {
        for (NodeID_ e = 0; e < numEpochs; ++e)
        {
          offsetMatrix[(e * numCacheLines) + c] = compressedOffsets[(c * numEpochs) + e];
        }
      }
      tm.Stop();
      std::cout << "[CSR-HYBRID-PREPROCESSING] Time to transpose offsets matrix =  " << tm.Seconds() << std::endl;
      
      /*
      // Sanity check - The following conditions will be checked
      // 1. If the MSB bit in an epoch is zero, then there should be a reference to the cacheline in that epoch
      // 2. If the MSB bit in an epoch is zero, then the remaining 7 bits should capture the last reference within the epoch
      // 3. If the MSB bit in an epoch is 1, then there should be no reference to the cacheline in that epoch
      // 4. If the MSB bit in an epoch is 1, then use the information in the remaining 7 bits to ascertain which epoch a line will be accessed next
      if (traverseCSR == true)
      {
        //Step I: compute per-cacheline neighborhoods
        std::vector<std::set<NodeID_> > lineNghs;
        lineNghs.resize(numCacheLines);
        #pragma omp parallel for schedule(dynamic, 8)
        for (NodeID_ c = 0; c < numCacheLines; ++c)
        {
          NodeID_ startVtx = c * numVtxPerLine;
          NodeID_ endVtx   = (c+1) * numVtxPerLine;
          if (c == numCacheLines - 1)
            endVtx = g.num_nodes();
          for (NodeID_ v = startVtx; v < endVtx; ++v)
          {
            for (NodeID_ ngh : g.out_neigh(v))
            {
              lineNghs[c].insert(ngh);
            }
          }
        }

        //Step II: Check conditions 1 & 3 for epoch data structure
        #pragma omp parallel for schedule(static)
        for (NodeID_ c = 0; c < numCacheLines; ++c)
        {
          for (NodeID_ e = 0; e < numEpochs; ++e)
          {
            if ((compressedOffsets[(c * numEpochs) + e] & orMask) != 0)
            {
              // Cond 3: There should be no reference this epoch
              NodeID_ epochStart = e * epochSz;
              NodeID_ epochEnd   = (e + 1) * epochSz;
              auto it = std::lower_bound(lineNghs[c].begin(), lineNghs[c].end(), epochStart);
              if (it != lineNghs[c].end())
              {
                NodeID_ lastRef = *it;
                assert(lastRef >= epochEnd);
              }
            }
            else
            {
              // Cond 1: There must be a reference this epoch
              NodeID_ epochStart = e * epochSz;
              NodeID_ epochEnd   = (e + 1) * epochSz;
              auto it = std::lower_bound(lineNghs[c].begin(), lineNghs[c].end(), epochStart);
              assert(it != lineNghs[c].end());
              NodeID_ lastRef = *it;
              assert(lastRef < epochEnd);
            }
          }
        }
        
        //Step III: Check conditions 2 & 4 for epoch data structure
        #pragma omp parallel for schedule(static)
        for (NodeID_ c = 0; c < numCacheLines; ++c)
        {
          for (NodeID_ e = 0; e < numEpochs; ++e)
          {
            if ((compressedOffsets[(c * numEpochs) + e] & orMask) != 0)
            {
              // Cond 4: The inter-epoch info points to the next epoch 
              uint8_t interEpochDist = compressedOffsets[(c * numEpochs) + e] & andMask;
              NodeID_ epochStart = (e + interEpochDist) * epochSz;
              auto it = std::lower_bound(lineNghs[c].begin(), lineNghs[c].end(), epochStart);
              if (it != lineNghs[c].end())
              {
                NodeID_ lastRef = *it;
                assert(lastRef >= epochStart); //minimum condition
                                           //because of quantization we cant be sure exactly which epoch (or if ever) the line will be accessed
              }
            }
            else
            {
              // Cond 2: We track the correct intra-Epoch lastRef 
              NodeID_ epochStart = e * epochSz;
              NodeID_ epochEnd   = (e + 1) * epochSz;
              auto it = std::lower_bound(lineNghs[c].begin(), lineNghs[c].end(), epochEnd);
              it--;
              NodeID_ lastRef = *it;
              assert(lastRef < epochEnd && lastRef >= epochStart);
              NodeID_ subEpochDist   = (*it) - epochStart;
              uint8_t subEpochDistQ  = subEpochDist / subEpochSz;
              uint8_t intraEpochDist = compressedOffsets[(c * numEpochs) + e] & andMask;
              if (e != numEpochs - 1)
                assert(intraEpochDist == subEpochDistQ);
            }
          }
        }
      }
      else 
      {
        //Step I: compute per-cacheline neighborhoods
        std::vector<std::set<NodeID_> > lineNghs;
        lineNghs.resize(numCacheLines);
        #pragma omp parallel for schedule(dynamic, 8)
        for (NodeID_ c = 0; c < numCacheLines; ++c)
        {
          NodeID_ startVtx = c * numVtxPerLine;
          NodeID_ endVtx   = (c+1) * numVtxPerLine;
          if (c == numCacheLines - 1)
            endVtx = g.num_nodes();
          for (NodeID_ v = startVtx; v < endVtx; ++v)
          {
            for (NodeID_ ngh : g.in_neigh(v))
            {
              lineNghs[c].insert(ngh);
            }
          }
        }

        //Step II: Check conditions 1 & 3 for epoch data structure
        #pragma omp parallel for schedule(static)
        for (NodeID_ c = 0; c < numCacheLines; ++c)
        {
          for (NodeID_ e = 0; e < numEpochs; ++e)
          {
            if ((compressedOffsets[(c * numEpochs) + e] & orMask) != 0)
            {
              // Cond 3: There should be no reference this epoch
              NodeID_ epochStart = e * epochSz;
              NodeID_ epochEnd   = (e + 1) * epochSz;
              auto it = std::lower_bound(lineNghs[c].begin(), lineNghs[c].end(), epochStart);
              if (it != lineNghs[c].end())
              {
                NodeID_ lastRef = *it;
                assert(lastRef >= epochEnd);
              }
            }
            else
            {
              // Cond 1: There must be a reference this epoch
              NodeID_ epochStart = e * epochSz;
              NodeID_ epochEnd   = (e + 1) * epochSz;
              auto it = std::lower_bound(lineNghs[c].begin(), lineNghs[c].end(), epochStart);
              assert(it != lineNghs[c].end());
              NodeID_ lastRef = *it;
              assert(lastRef < epochEnd);
            }
          }
        }
        
        //Step III: Check conditions 2 & 4 for epoch data structure
        #pragma omp parallel for schedule(static)
        for (NodeID_ c = 0; c < numCacheLines; ++c)
        {
          for (NodeID_ e = 0; e < numEpochs; ++e)
          {
            if ((compressedOffsets[(c * numEpochs) + e] & orMask) != 0)
            {
              // Cond 4: The inter-epoch info points to the next epoch 
              uint8_t interEpochDist = compressedOffsets[(c * numEpochs) + e] & andMask;
              NodeID_ epochStart = (e + interEpochDist) * epochSz;
              auto it = std::lower_bound(lineNghs[c].begin(), lineNghs[c].end(), epochStart);
              if (it != lineNghs[c].end())
              {
                NodeID_ lastRef = *it;
                assert(lastRef >= epochStart); //minimum condition
                                           //because of quantization we cant be sure exactly which epoch (or if ever) the line will be accessed
              }
            }
            else
            {
              // Cond 2: We track the correct intra-Epoch lastRef 
              NodeID_ epochStart = e * epochSz;
              NodeID_ epochEnd   = (e + 1) * epochSz;
              auto it = std::lower_bound(lineNghs[c].begin(), lineNghs[c].end(), epochEnd);
              it--;
              NodeID_ lastRef = *it;
              assert(lastRef < epochEnd && lastRef >= epochStart);
              NodeID_ subEpochDist   = (*it) - epochStart;
              uint8_t subEpochDistQ  = subEpochDist / subEpochSz;
              uint8_t intraEpochDist = compressedOffsets[(c * numEpochs) + e] & andMask;
              if (e != numEpochs - 1)
                assert(intraEpochDist == subEpochDistQ);
            }
          }
        }
      }
      */
    }
};

#endif  // BUILDER_H_
