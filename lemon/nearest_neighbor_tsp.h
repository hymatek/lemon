#ifndef LEMON_NEAREST_NEIGHBOUR_TSP_H
#define LEMON_NEAREST_NEIGHBOUR_TSP_H

#include <deque>
#include <lemon/full_graph.h>
#include <lemon/path.h>
#include <lemon/maps.h>

namespace lemon {
  
  namespace nn_helper {
    template <class L>
    L dequeConvert(const std::deque<FullGraph::Node> &_path) {
      return L(_path.begin(), _path.end());
    }

    template <>
    std::deque<FullGraph::Node> dequeConvert(const std::deque<FullGraph::Node> &_path) {
      return _path;
    }
  }
  
  template <typename CM>
  class NearestNeighborTsp {
    private:
      GRAPH_TYPEDEFS(FullGraph);

    public:
      typedef CM CostMap;
      typedef typename CM::Value Cost;
    
      NearestNeighborTsp(const FullGraph &gr, const CostMap &cost) : _gr(gr), _cost(cost) {}

      Cost run() {
        _path.clear();

        Edge min_edge1 = INVALID,
             min_edge2 = INVALID;
        
        min_edge1 = mapMin(_gr, _cost);

        FullGraph::NodeMap<bool> used(_gr, false);

        Node n1 = _gr.u(min_edge1), 
             n2 = _gr.v(min_edge1);
        
        _path.push_back(n1);
        _path.push_back(n2);

        used[n1] = true;
        used[n2] = true;

        min_edge1 = INVALID;

        while (int(_path.size()) != _gr.nodeNum()) {
          if (min_edge1 == INVALID) {
            for (IncEdgeIt e(_gr, n1); e!=INVALID; ++e) {
              if (!used[_gr.runningNode(e)]) {
                if (min_edge1==INVALID || _cost[min_edge1] > _cost[e])
                  min_edge1 = e;
              }
            }
          }

          if (min_edge2 == INVALID) {
            for (IncEdgeIt e(_gr, n2); e!=INVALID; ++e) {
              if (!used[_gr.runningNode(e)]) {
                if (min_edge2==INVALID || _cost[min_edge2] > _cost[e])
                  min_edge2 = e;
              }
            }
          }

          if ( _cost[min_edge1] < _cost[min_edge2] ) {
            n1 = (_gr.u(min_edge1) == n1) ? _gr.v(min_edge1) : _gr.u(min_edge1);
            _path.push_front(n1);

            used[n1] = true;
            min_edge1 = INVALID;

            if (_gr.u(min_edge2)==n1 || _gr.v(min_edge2)==n1)
              min_edge2 = INVALID;
          } else {
            n2 = (_gr.u(min_edge2) == n2) ? _gr.v(min_edge2) : _gr.u(min_edge2);
            _path.push_back(n2);

            used[n2] = true;
            min_edge2 = INVALID;

            if (_gr.u(min_edge1)==n2 || _gr.v(min_edge1)==n2)
              min_edge1 = INVALID;
          }
        }

        _sum = _cost[ _gr.edge(_path.back(), _path.front()) ];
        for (unsigned int i=0; i<_path.size()-1; ++i)
          _sum += _cost[ _gr.edge(_path[i], _path[i+1]) ];

        return _sum;
      }

      
      template <typename L>
      void tourNodes(L &container) {
        container(nn_helper::dequeConvert<L>(_path));
      }
      
      template <template <typename> class L>
      L<Node> tourNodes() {
        return nn_helper::dequeConvert<L<Node> >(_path);
      }      

      const std::deque<Node>& tourNodes() {
        return _path;
      }
      
      Path<FullGraph> tour() {
        Path<FullGraph> p;
        if (_path.size()<2)
          return p;
        
        for (unsigned int i=0; i<_path.size()-1; ++i) {
          p.addBack(_gr.arc(_path[i], _path[i+1]));
        }
        p.addBack(_gr.arc(_path.back(), _path.front()));
        
        return p;
      }
      
      Cost tourCost() {
        return _sum;
      }
      

  private:
    const FullGraph &_gr;
    const CostMap &_cost;
    Cost _sum;
    std::deque<Node> _path;
  };


}; // namespace lemon

#endif
