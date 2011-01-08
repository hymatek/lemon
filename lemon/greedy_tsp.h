#ifndef LEMON_GREEDY_TSP_H
#define LEMON_GREEDY_TSP_H

#include <lemon/path.h>
#include <lemon/full_graph.h>
#include <lemon/unionfind.h>
#include <algorithm>

namespace lemon {

  namespace greedy_tsp_helper {

    template <typename CostMap>
    class KeyComp {
      typedef typename CostMap::Key Key;
      const CostMap &cost;
      
      public:
        KeyComp(const CostMap &_cost) : cost(_cost) {}
      
        bool operator() (const Key &a, const Key &b) const {
          return cost[a] < cost[b];
        }
    };

  
    template <class L>
    L vectorConvert(const std::vector<FullGraph::Node> &_path) {
      return L(_path.begin(), _path.end());
    }
  
    template <>
    std::vector<FullGraph::Node> vectorConvert(const std::vector<FullGraph::Node> &_path) {
      return _path;
    }
    
  }


  template <typename CM>
  class GreedyTsp {
    private:
      GRAPH_TYPEDEFS(FullGraph);
    
    public:
      typedef CM CostMap;
      typedef typename CM::Value Cost;
      
      GreedyTsp(const FullGraph &gr, const CostMap &cost) : _gr(gr), _cost(cost) {}  

      Cost run() {
        typedef UnionFind<FullGraph::NodeMap<int> > Union;
        _nodes.clear();
        
        std::vector<int> path;
        path.resize(_gr.nodeNum()*2, -1);
        
        std::vector<typename CostMap::Key> sorted_edges;
        sorted_edges.reserve(_gr.edgeNum());
        for (EdgeIt n(_gr); n != INVALID; ++n)
          sorted_edges.push_back(n);

        std::sort(sorted_edges.begin(), sorted_edges.end(), greedy_tsp_helper::KeyComp<CostMap>(_cost));

        FullGraph::NodeMap<int> nodemap(_gr);
        Union unionfind(nodemap);

        for (NodeIt n(_gr); n != INVALID; ++n)
          unionfind.insert(n);
        
        FullGraph::NodeMap<int> degree(_gr, 0);

        int nodesNum = 0, i = 0;

        while ( nodesNum != _gr.nodeNum()-1 ) {
          const Edge &e = sorted_edges[i];
          
          const Node u = _gr.u(e),
                     v = _gr.v(e);
          
          if (degree[u]<=1 && degree[v]<=1) {
            if (unionfind.join(u, v)) {
              ++degree[u];
              ++degree[v];
              ++nodesNum;
              
              const int uid = _gr.id(u),
                        vid = _gr.id(v);
              
              
              path[uid*2 + (path[uid*2]==-1 ? 0 : 1)] = vid;
              path[vid*2 + (path[vid*2]==-1 ? 0 : 1)] = uid;
            }
          }

          ++i;
        }


        for (int i=0, n=-1; i<_gr.nodeNum()*2; ++i) {
          if (path[i] == -1) {
            if (n==-1) {
              n = i;
            } else {
              path[n] = i/2;
              path[i] = n/2;
              break;
            }
          }
        }


        for (int i=0, j=0, last=-1; i!=_gr.nodeNum(); ++i) {
          _nodes.push_back(_gr.nodeFromId(j));
          
          if (path[2*j] != last) {
            last = j;
            j = path[2*j];
          } else {
            last = j;
            j = path[2*j+1];
          }
        }

        _sum = _cost[_gr.edge(_nodes.back(), _nodes.front())];
        for (unsigned int i=0; i<_nodes.size()-1; ++i)
          _sum += _cost[_gr.edge(_nodes[i], _nodes[i+1])];

        return _sum;
      }



      template <typename L>
      void tourNodes(L &container) {
        container(greedy_tsp_helper::vectorConvert<L>(_nodes));
      }
      
      template <template <typename> class L>
      L<Node> tourNodes() {
        return greedy_tsp_helper::vectorConvert<L<Node> >(_nodes);
      }
      
      const std::vector<Node>& tourNodes() {
        return _nodes;
      }
      
      Path<FullGraph> tour() {
        Path<FullGraph> p;
        if (_nodes.size()<2)
          return p;
        
        for (unsigned int i=0; i<_nodes.size()-1; ++i) {
          p.addBack(_gr.arc(_nodes[i], _nodes[i+1]));
        }
        
        p.addBack(_gr.arc(_nodes.back(), _nodes.front()));
        
        return p;
      }
      
      Cost tourCost() {
        return _sum;
      }
      

    private:
      const FullGraph &_gr;
      const CostMap &_cost;
      Cost _sum;
      std::vector<Node> _nodes;
  };

}; // namespace lemon

#endif
