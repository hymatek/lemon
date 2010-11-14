/* -*- mode: C++; indent-tabs-mode: nil; -*-
 *
 * This file is a part of LEMON, a generic C++ optimization library.
 *
 * Copyright (C) 2003-2010
 * Egervary Jeno Kombinatorikus Optimalizalasi Kutatocsoport
 * (Egervary Research Group on Combinatorial Optimization, EGRES).
 *
 * Permission to use, modify and distribute this software is granted
 * provided that this copyright notice appears in all copies. For
 * precise terms see the accompanying LICENSE file.
 *
 * This software is provided "AS IS" with no warranty of any kind,
 * express or implied, and with no claim as to its suitability for any
 * purpose.
 *
 */

#include <lemon/concepts/bpgraph.h>
//#include <lemon/list_graph.h>
//#include <lemon/smart_graph.h>
//#include <lemon/full_graph.h>
//#include <lemon/grid_graph.h>
//#include <lemon/hypercube_graph.h>

#include "test_tools.h"
#include "graph_test.h"

using namespace lemon;
using namespace lemon::concepts;

void checkConcepts() {
  { // Checking graph components
    checkConcept<BaseBpGraphComponent, BaseBpGraphComponent >();

    checkConcept<IDableBpGraphComponent<>,
      IDableBpGraphComponent<> >();

    checkConcept<IterableBpGraphComponent<>,
      IterableBpGraphComponent<> >();

    checkConcept<AlterableBpGraphComponent<>,
      AlterableBpGraphComponent<> >();

    checkConcept<MappableBpGraphComponent<>,
      MappableBpGraphComponent<> >();

    checkConcept<ExtendableBpGraphComponent<>,
      ExtendableBpGraphComponent<> >();

    checkConcept<ErasableBpGraphComponent<>,
      ErasableBpGraphComponent<> >();

    checkConcept<ClearableGraphComponent<>,
      ClearableGraphComponent<> >();

  }
  { // Checking skeleton graph
    checkConcept<BpGraph, BpGraph>();
  }
}

void checkGraphs() {
}

int main() {
  checkConcepts();
  checkGraphs();
  return 0;
}
