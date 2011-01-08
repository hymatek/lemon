#ifndef LEMON_CHRISTOFIDES_TSP_H
#define LEMON_CHRISTOFIDES_TSP_H

#include <lemon/full_graph.h>
#include <lemon/smart_graph.h>
#include <lemon/path.h>
#include <lemon/kruskal.h>
#include <lemon/matching.h>
#include <lemon/adaptors.h>
#include <lemon/maps.h>
#include <lemon/euler.h>

namespace lemon {
  
  namespace christofides_helper {
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
  class ChristofidesTsp {
    private:
      GRAPH_TYPEDEFS(SmartGraph);

    public:
      typedef typename CM::Value Cost;
      typedef SmartGraph::EdgeMap<Cost> CostMap;
    
      ChristofidesTsp(const FullGraph &gr, const CM &cost) : _cost(_gr), fullg(gr), fullcost(cost), nr(_gr) {
        graphCopy(gr, _gr).nodeCrossRef(nr).edgeMap(cost, _cost).run();
      }

      Cost run() {
        _path.clear();
        
        SmartGraph::EdgeMap<bool> tree(_gr);
        kruskal(_gr, _cost, tree);
        
        FilterEdges<SmartGraph> treefiltered(_gr, tree);
        InDegMap<FilterEdges<SmartGraph> > deg(treefiltered);
        
        SmartGraph::NodeMap<bool> oddfilter(_gr, false);
        FilterNodes<SmartGraph> oddNodes(_gr, oddfilter);
        
        for (NodeIt n(_gr); n!=INVALID; ++n) {
          if (deg[n]%2 == 1) {
            oddNodes.enable(n);
          }
        }
        
        NegMap<CostMap> negmap(_cost);
        MaxWeightedPerfectMatching<FilterNodes<SmartGraph>,
                  NegMap<CostMap> > perfmatch(oddNodes, negmap);
        perfmatch.run();
        
        for (FilterNodes<SmartGraph>::EdgeIt e(oddNodes); e!=INVALID; ++e) {
          if (perfmatch.matching(e)) {
            treefiltered.enable(_gr.addEdge(_gr.u(e), _gr.v(e)));
          }
        }
        
        FilterEdges<SmartGraph>::NodeMap<bool> seen(treefiltered, false);
        for (EulerIt<FilterEdges<SmartGraph> > e(treefiltered); e!=INVALID; ++e) {
          if (seen[treefiltered.target(e)]==false) {
            _path.push_back(nr[treefiltered.target(e)]);
            seen[treefiltered.target(e)] = true;
          }
        }

        _sum = fullcost[ fullg.edge(_path.back(), _path.front()) ];
        for (unsigned int i=0; i<_path.size()-1; ++i)
          _sum += fullcost[ fullg.edge(_path[i], _path[i+1]) ];

        return _sum;
      }

      template <typename L>
      void tourNodes(L &container) {
        container(christofides_helper::vectorConvert<L>(_path));
      }
      
      template <template <typename> class L>
      L<FullGraph::Node> tourNodes() {
        return christofides_helper::vectorConvert<L<FullGraph::Node> >(_path);
      }

      const std::vector<Node>& tourNodes() {
        return _path;
      }
      
      Path<FullGraph> tour() {
        Path<FullGraph> p;
        if (_path.size()<2)
          return p;

        for (unsigned int i=0; i<_path.size()-1; ++i) {
          p.addBack(fullg.arc(_path[i], _path[i+1]));
        }
        p.addBack(fullg.arc(_path.back(), _path.front()));
        
        return p;
      }
      
      Cost tourCost() {
        return _sum;
      }
      

  private:
    SmartGraph _gr;
    CostMap _cost;
    Cost _sum;
    const FullGraph &fullg;
    const CM &fullcost;
    std::vector<FullGraph::Node> _path;
    SmartGraph::NodeMap<FullGraph::Node> nr;
  };


}; // namespace lemon

#endif
