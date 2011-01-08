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

#ifndef LEMON_INSERTION_TSP_H
#define LEMON_INSERTION_TSP_H

/// \ingroup tsp
/// \file
/// \brief Insertion algorithm for symmetric TSP

#include <vector>
#include <lemon/full_graph.h>
#include <lemon/maps.h>
#include <lemon/random.h>

namespace lemon {

  /// \ingroup tsp
  ///
  /// \brief Insertion algorithm for symmetric TSP.
  ///
  /// InsertionTsp implements the insertion heuristic for solving
  /// symmetric \ref tsp "TSP".
  ///
  /// This is a basic TSP heuristic that has many variants.
  /// It starts with a subtour containing a few nodes of the graph and it
  /// iteratively inserts the other nodes into this subtour according to a
  /// certain node selection rule.
  ///
  /// This implementation provides four different node selection rules,
  /// from which the most powerful one is used by default.
  /// For more information, see \ref SelectionRule.
  ///
  /// \tparam CM Type of the cost map.
  template <typename CM>
  class InsertionTsp
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
      std::vector<Node> _notused;
      std::vector<Node> _path;
      Cost _sum;

    public:

      /// \brief Constants for specifying the node selection rule.
      ///
      /// Enum type containing constants for specifying the node selection
      /// rule for the \ref run() function.
      ///
      /// During the algorithm, nodes are selected for addition to the current
      /// subtour according to the applied rule.
      /// In general, the FARTHEST method yields the best tours, thus it is the
      /// default option. The RANDOM rule usually gives somewhat worse results,
      /// but it is much faster than the others and it is the most robust.
      ///
      /// The desired selection rule can be specified as a parameter of the
      /// \ref run() function.
      enum SelectionRule {

        /// An unvisited node having minimum distance from the current
        /// subtour is selected at each step.
        /// The algorithm runs in O(n<sup>3</sup>) time using this
        /// selection rule.
        NEAREST,

        /// An unvisited node having maximum distance from the current
        /// subtour is selected at each step.
        /// The algorithm runs in O(n<sup>3</sup>) time using this
        /// selection rule.
        FARTHEST,

        /// An unvisited node whose insertion results in the least
        /// increase of the subtour's total cost is selected at each step.
        /// The algorithm runs in O(n<sup>3</sup>) time using this
        /// selection rule.
        CHEAPEST,

        /// An unvisited node is selected randomly without any evaluation
        /// at each step.
        /// The global \ref rnd "random number generator instance" is used.
        /// You can seed it before executing the algorithm, if you
        /// would like to.
        /// The algorithm runs in O(n<sup>2</sup>) time using this
        /// selection rule.
        RANDOM
      };

    public:

      /// \brief Constructor
      ///
      /// Constructor.
      /// \param gr The \ref FullGraph "full graph" the algorithm runs on.
      /// \param cost The cost map.
      InsertionTsp(const FullGraph &gr, const CostMap &cost)
        : _gr(gr), _cost(cost) {}

      /// \name Execution Control
      /// @{

      /// \brief Runs the algorithm.
      ///
      /// This function runs the algorithm.
      ///
      /// \param rule The node selection rule. For more information, see
      /// \ref SelectionRule.
      ///
      /// \return The total cost of the found tour.
      Cost run(SelectionRule rule = FARTHEST) {
        _path.clear();

        if (_gr.nodeNum() == 0) return _sum = 0;
        else if (_gr.nodeNum() == 1) {
          _path.push_back(_gr(0));
          return _sum = 0;
        }

        switch (rule) {
          case NEAREST:
            init(true);
            start<NearestSelection, DefaultInsertion>();
            break;
          case FARTHEST:
            init(false);
            start<FarthestSelection, DefaultInsertion>();
            break;
          case CHEAPEST:
            init(true);
            start<CheapestSelection, CheapestInsertion>();
            break;
          case RANDOM:
            init(true);
            start<RandomSelection, DefaultInsertion>();
            break;
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
      /// This function returns a const reference to a vector
      /// that stores the node sequence of the found tour.
      ///
      /// \pre run() must be called before using this function.
      const std::vector<Node>& tourNodes() const {
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

    private:

      // Initializes the algorithm
      void init(bool min) {
        Edge min_edge = min ? mapMin(_gr, _cost) : mapMax(_gr, _cost);

        _path.clear();
        _path.push_back(_gr.u(min_edge));
        _path.push_back(_gr.v(min_edge));

        _notused.clear();
        for (NodeIt n(_gr); n!=INVALID; ++n) {
          if (n != _gr.u(min_edge) && n != _gr.v(min_edge)) {
            _notused.push_back(n);
          }
        }

        _sum = _cost[min_edge] * 2;
      }

      // Executes the algorithm
      template <class SelectionFunctor, class InsertionFunctor>
      void start() {
        SelectionFunctor selectNode(_gr, _cost, _path, _notused);
        InsertionFunctor insertNode(_gr, _cost, _path, _sum);

        for (int i=0; i<_gr.nodeNum()-2; ++i) {
          insertNode.insert(selectNode.select());
        }

        _sum = _cost[_gr.edge(_path.back(), _path.front())];
        for (int i = 0; i < int(_path.size())-1; ++i) {
          _sum += _cost[_gr.edge(_path[i], _path[i+1])];
        }
      }


      // Implementation of the nearest selection rule
      class NearestSelection {
        public:
          NearestSelection(const FullGraph &gr, const CostMap &cost,
                           std::vector<Node> &path, std::vector<Node> &notused)
            : _gr(gr), _cost(cost), _path(path), _notused(notused) {}

          Node select() const {
            Cost insert_val = 0;
            int insert_node = -1;

            for (unsigned int i=0; i<_notused.size(); ++i) {
              Cost min_val = _cost[_gr.edge(_notused[i], _path[0])];
              int min_node = 0;

              for (unsigned int j=1; j<_path.size(); ++j) {
                Cost curr = _cost[_gr.edge(_notused[i], _path[j])];
                if (min_val > curr) {
                    min_val = curr;
                    min_node = j;
                }
              }

              if (insert_val > min_val || insert_node == -1) {
                insert_val = min_val;
                insert_node = i;
              }
            }

            Node n = _notused[insert_node];
            _notused.erase(_notused.begin()+insert_node);

            return n;
          }

        private:
          const FullGraph &_gr;
          const CostMap &_cost;
          std::vector<Node> &_path;
          std::vector<Node> &_notused;
      };


      // Implementation of the farthest selection rule
      class FarthestSelection {
        public:
          FarthestSelection(const FullGraph &gr, const CostMap &cost,
                            std::vector<Node> &path, std::vector<Node> &notused)
            : _gr(gr), _cost(cost), _path(path), _notused(notused) {}

          Node select() const {
            Cost insert_val = 0;
            int insert_node = -1;

            for (unsigned int i=0; i<_notused.size(); ++i) {
              Cost min_val = _cost[_gr.edge(_notused[i], _path[0])];
              int min_node = 0;

              for (unsigned int j=1; j<_path.size(); ++j) {
                Cost curr = _cost[_gr.edge(_notused[i], _path[j])];
                if (min_val > curr) {
                  min_val = curr;
                  min_node = j;
                }
              }

              if (insert_val < min_val || insert_node == -1) {
                insert_val = min_val;
                insert_node = i;
              }
            }

            Node n = _notused[insert_node];
            _notused.erase(_notused.begin()+insert_node);

            return n;
          }

        private:
          const FullGraph &_gr;
          const CostMap &_cost;
          std::vector<Node> &_path;
          std::vector<Node> &_notused;
      };


      // Implementation of the cheapest selection rule
      class CheapestSelection {
        private:
          Cost costDiff(Node u, Node v, Node w) const {
            return
              _cost[_gr.edge(u, w)] +
              _cost[_gr.edge(v, w)] -
              _cost[_gr.edge(u, v)];
          }

        public:
          CheapestSelection(const FullGraph &gr, const CostMap &cost,
                            std::vector<Node> &path, std::vector<Node> &notused)
            : _gr(gr), _cost(cost), _path(path), _notused(notused) {}

          Cost select() const {
            int insert_notused = -1;
            int best_insert_index = -1;
            Cost insert_val = 0;

            for (unsigned int i=0; i<_notused.size(); ++i) {
              int min = i;
              int best_insert_tmp = 0;
              Cost min_val =
                costDiff(_path.front(), _path.back(), _notused[i]);

              for (unsigned int j=1; j<_path.size(); ++j) {
                Cost tmp =
                  costDiff(_path[j-1], _path[j], _notused[i]);

                if (min_val > tmp) {
                  min = i;
                  min_val = tmp;
                  best_insert_tmp = j;
                }
              }

              if (insert_val > min_val || insert_notused == -1) {
                insert_notused = min;
                insert_val = min_val;
                best_insert_index = best_insert_tmp;
              }
            }

            _path.insert(_path.begin()+best_insert_index,
                         _notused[insert_notused]);
            _notused.erase(_notused.begin()+insert_notused);

            return insert_val;
          }

        private:
          const FullGraph &_gr;
          const CostMap &_cost;
          std::vector<Node> &_path;
          std::vector<Node> &_notused;
      };

      // Implementation of the random selection rule
      class RandomSelection {
        public:
          RandomSelection(const FullGraph &, const CostMap &,
                          std::vector<Node> &, std::vector<Node> &notused)
            : _notused(notused) {}

          Node select() const {
            const int index = rnd[_notused.size()];
            Node n = _notused[index];
            _notused.erase(_notused.begin()+index);
            return n;
          }
        private:
          std::vector<Node> &_notused;
      };


      // Implementation of the default insertion method
      class DefaultInsertion {
        private:
          Cost costDiff(Node u, Node v, Node w) const {
            return
              _cost[_gr.edge(u, w)] +
              _cost[_gr.edge(v, w)] -
              _cost[_gr.edge(u, v)];
          }

        public:
          DefaultInsertion(const FullGraph &gr, const CostMap &cost,
                           std::vector<Node> &path, Cost &total_cost) :
            _gr(gr), _cost(cost), _path(path), _total(total_cost) {}

          void insert(Node n) const {
            int min = 0;
            Cost min_val =
              costDiff(_path.front(), _path.back(), n);

            for (unsigned int i=1; i<_path.size(); ++i) {
              Cost tmp = costDiff(_path[i-1], _path[i], n);
              if (tmp < min_val) {
                min = i;
                min_val = tmp;
              }
            }

            _path.insert(_path.begin()+min, n);
            _total += min_val;
          }

        private:
          const FullGraph &_gr;
          const CostMap &_cost;
          std::vector<Node> &_path;
          Cost &_total;
      };

      // Implementation of a special insertion method for the cheapest
      // selection rule
      class CheapestInsertion {
        TEMPLATE_GRAPH_TYPEDEFS(FullGraph);
        public:
          CheapestInsertion(const FullGraph &, const CostMap &,
                            std::vector<Node> &, Cost &total_cost) :
            _total(total_cost) {}

          void insert(Cost diff) const {
            _total += diff;
          }

        private:
          Cost &_total;
      };

  };

}; // namespace lemon

#endif
