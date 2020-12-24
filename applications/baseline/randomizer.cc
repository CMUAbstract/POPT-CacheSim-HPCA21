// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <iostream>
#include <omp.h>

#include "benchmark.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "reader.h"
#include "writer.h"
#include "pvector.h"

using namespace std;

int main(int argc, char* argv[]) {
  CLConvert cli(argc, argv, "converter");
  cli.ParseArgs();
  if (cli.out_weighted()) {
    assert(false); //TODO: add randomization support for weighted graphs
    #if 0
    WeightedBuilder bw(cli);
    WGraph wg = bw.MakeGraph();
    wg.PrintStats();
    WeightedWriter ww(wg);
    ww.WriteGraph(cli.out_filename(), cli.out_sg());
    #endif
  } else {
    Builder b(cli);
    Graph g = b.MakeGraph();
    pvector<NodeID> newIds(g.num_nodes(), -1);
    Graph rand_g = Builder::RandOrder(g, newIds, false, true);
    rand_g.PrintStats();
    Writer w(rand_g);
    w.WriteGraph(cli.out_filename(), cli.out_sg());
  }
  return 0;
}
