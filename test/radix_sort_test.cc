/* -*- C++ -*-
 *
 * This file is a part of LEMON, a generic C++ optimization library
 *
 * Copyright (C) 2003-2008
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

#include <lemon/time_measure.h>
#include <lemon/smart_graph.h>
#include <lemon/maps.h>
#include <lemon/radix_sort.h>
#include <lemon/math.h>

#include "test_tools.h"

#include <vector>
#include <algorithm>

using namespace lemon;

static const int n = 10000;

struct Negate {
  typedef int argument_type;
  typedef int result_type;
  int operator()(int a) { return - a; }
};

int negate(int a) { return - a; }


void generateIntSequence(int n, std::vector<int>& data) {
  int prime = 9973;
  int root = 136, value = 1;
  for (int i = 0; i < n; ++i) {
    data.push_back(value - prime / 2);
    value = (value * root) % prime;
  }
}

void generateCharSequence(int n, std::vector<unsigned char>& data) {
  int prime = 251;
  int root = 3, value = root;
  for (int i = 0; i < n; ++i) {
    data.push_back(static_cast<unsigned char>(value));
    value = (value * root) % prime;
  }
}

void checkRadixSort() {
  {
    std::vector<int> data1;
    generateIntSequence(n, data1);

    std::vector<int> data2(data1);
    std::sort(data1.begin(), data1.end());

    radixSort(data2.begin(), data2.end());
    for (int i = 0; i < n; ++i) {
      check(data1[i] == data2[i], "Test failed");
    }

    radixSort(data2.begin(), data2.end(), Negate());
    for (int i = 0; i < n; ++i) {
      check(data1[i] == data2[n - 1 - i], "Test failed");
    }

    radixSort(data2.begin(), data2.end(), negate);
    for (int i = 0; i < n; ++i) {
      check(data1[i] == data2[n - 1 - i], "Test failed");
    }

  } 
  
  {
    std::vector<unsigned char> data1(n);
    generateCharSequence(n, data1);

    std::vector<unsigned char> data2(data1);
    std::sort(data1.begin(), data1.end());

    radixSort(data2.begin(), data2.end());
    for (int i = 0; i < n; ++i) {
      check(data1[i] == data2[i], "Test failed");
    }

  }
}


void checkCounterSort() {
  {
    std::vector<int> data1;
    generateIntSequence(n, data1);

    std::vector<int> data2(data1);
    std::sort(data1.begin(), data1.end());

    counterSort(data2.begin(), data2.end());
    for (int i = 0; i < n; ++i) {
      check(data1[i] == data2[i], "Test failed");
    }

    counterSort(data2.begin(), data2.end(), Negate());
    for (int i = 0; i < n; ++i) {
      check(data1[i] == data2[n - 1 - i], "Test failed");
    }

    counterSort(data2.begin(), data2.end(), negate);
    for (int i = 0; i < n; ++i) {
      check(data1[i] == data2[n - 1 - i], "Test failed");
    }
  } 

  {
    std::vector<unsigned char> data1(n);
    generateCharSequence(n, data1);

    std::vector<unsigned char> data2(data1);
    std::sort(data1.begin(), data1.end());

    radixSort(data2.begin(), data2.end());
    for (int i = 0; i < n; ++i) {
      check(data1[i] == data2[i], "Test failed");
    }

  }
}

int main() {

  checkRadixSort();  
  checkCounterSort();

  return 0;
}
