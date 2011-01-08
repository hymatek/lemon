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

#ifndef LEMON_NEAREST_NEIGHBOUR_TSP_H
#define LEMON_NEAREST_NEIGHBOUR_TSP_H

/// \ingroup tsp
/// \file
/// \brief Nearest neighbor algorithm for symmetric TSP

#include <deque>
#include <limits>
#include <lemon/full_graph.h>
#include <lemon/maps.h>

namespace lemon {

  /// \brief Nearest neighbor algorithm for symmetric TSP.
  ///
  /// NearestNeighborTsp implements the nearest neighbor heuristic for solving
  /// symmetric \ref tsp "TSP".
  ///
  /// This is probably the simplest TSP heuristic.
  /// It starts with a minimum cost edge and at each step, it connects the
  /// nearest unvisited node to the current path.
  /// Finally, it connects the two end points of the path to form a tour.
  ///
  /// This method runs in O(n<sup>2</sup>) time.
  /// It quickly finds an effectively short tour for most TSP
  /// instances, but in special cases, it could yield a really bad
  /// (or even the worst) solution.
  ///
  /// \tparam CM Type of the cost map.
  template <typename CM>
  class NearestNeighborTsp
  {
    public:

      /// Type of the cost map
      typedef CM CostMap;
      /// Type of the edge costs
      typedef typename CM::Value Cost;

    private:

      GRAPH_TYPEDEFS(FullGraph);

      const FullGraph &_gr;
      const CostMap &_cost;
      Cost _sum;
      std::deque<Node> _path;

    public:

      /// \brief Constructor
      ///
      /// Constructor.
      /// \param gr The \ref FullGraph "full graph" the algorithm runs on.
      /// \param cost The cost map.
      NearestNeighborTsp(const FullGraph &gr, const CostMap &cost)
        : _gr(gr), _cost(cost) {}

      /// \name Execution Control
      /// @{

      /// \brief Runs the algorithm.
      ///
      /// This function runs the algorithm.
      ///
      /// \return The total cost of the found tour.
      Cost run() {
        _path.clear();

        if (_gr.nodeNum() == 0) return _sum = 0;
        else if (_gr.nodeNum() == 1) {
          _path.push_back(_gr(0));
          return _sum = 0;
        }

        Edge min_edge1 = INVALID,
             min_edge2 = INVALID;

        min_edge1 = mapMin(_gr, _cost);
        Node n1 = _gr.u(min_edge1),
             n2 = _gr.v(min_edge1);
        _path.push_back(n1);
        _path.push_back(n2);

        FullGraph::NodeMap<bool> used(_gr, false);
        used[n1] = true;
        used[n2] = true;

        min_edge1 = INVALID;
        while (int(_path.size()) != _gr.nodeNum()) {
          if (min_edge1 == INVALID) {
            for (IncEdgeIt e(_gr, n1); e != INVALID; ++e) {
              if (!used[_gr.runningNode(e)] &&
                  (_cost[e] < _cost[min_edge1] || min_edge1 == INVALID)) {
                min_edge1 = e;
              }
            }
          }

          if (min_edge2 == INVALID) {
            for (IncEdgeIt e(_gr, n2); e != INVALID; ++e) {
              if (!used[_gr.runningNode(e)] &&
                  (_cost[e] < _cost[min_edge2] || min_edge2 == INVALID)) {
                min_edge2 = e;
              }
            }
          }

          if (_cost[min_edge1] < _cost[min_edge2]) {
            n1 = _gr.oppositeNode(n1, min_edge1);
            _path.push_front(n1);

            used[n1] = true;
            min_edge1 = INVALID;

            if (_gr.u(min_edge2) == n1 || _gr.v(min_edge2) == n1)
              min_edge2 = INVALID;
          } else {
            n2 = _gr.oppositeNode(n2, min_edge2);
            _path.push_back(n2);

            used[n2] = true;
            min_edge2 = INVALID;

            if (_gr.u(min_edge1) == n2 || _gr.v(min_edge1) == n2)
              min_edge1 = INVALID;
          }
        }

        _sum = _cost[_gr.edge(_path.back(), _path.front())];
        for (int i = 0; i < int(_path.size())-1; ++i) {
          _sum += _cost[_gr.edge(_path[i], _path[i+1])];
        }

        return _sum;
      }

      /// @}

      /// \name Query Functions
      /// @{

      /// \brief The total cost of the found tour.
      ///
      /// This function returns the total cost of the found tour.
      ///
      /// \pre run() must be called before using this function.
      Cost tourCost() const {
        return _sum;
      }

      /// \brief Returns a const reference to the node sequence of the
      /// found tour.
      ///
      /// This function returns a const reference to the internal structure
      /// that stores the node sequence of the found tour.
      ///
      /// \pre run() must be called before using this function.
      const std::deque<Node>& tourNodes() const {
        return _path;
      }

      /// \brief Gives back the node sequence of the found tour.
      ///
      /// This function copies the node sequence of the found tour into
      /// the given standard container.
      ///
      /// \pre run() must be called before using this function.
      template <typename Container>
      void tourNodes(Container &container) const {
        container.assign(_path.begin(), _path.end());
      }

      /// \brief Gives back the found tour as a path.
      ///
      /// This function copies the found tour as a list of arcs/edges into
      /// the given \ref concept::Path "path structure".
      ///
      /// \pre run() must be called before using this function.
      template <typename Path>
      void tour(Path &path) const {
        path.clear();
        for (int i = 0; i < int(_path.size()) - 1; ++i) {
          path.addBack(_gr.arc(_path[i], _path[i+1]));
        }
        if (int(_path.size()) >= 2) {
          path.addBack(_gr.arc(_path.back(), _path.front()));
        }
      }

      /// @}

  };

}; // namespace lemon

#endif
