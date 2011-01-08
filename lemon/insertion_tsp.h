#ifndef LEMON_INSERTION_TSP_H
#define LEMON_INSERTION_TSP_H

#include <lemon/full_graph.h>
#include <lemon/path.h>
#include <lemon/maps.h>
#include <lemon/random.h>
#include <vector>

namespace lemon {

  namespace insertion_tsp_helper {
  
    template <class L>
    L vectorConvert(const std::vector<FullGraph::Node> &_path) {
      return L(_path.begin(), _path.end());
    };
  
    template <>
    std::vector<FullGraph::Node> vectorConvert(
        const std::vector<FullGraph::Node> &_path) {
      return _path;
    };
    
  };

  template <typename CM>
  class InsertionTsp {
    private:
      GRAPH_TYPEDEFS(FullGraph);

    public:      
      typedef CM CostMap;
      typedef typename CM::Value Cost;
      
      InsertionTsp(const FullGraph &gr, const CostMap &cost) : 
            _gr(gr), _cost(cost) {}
      
      enum InsertionMethod {
        INSERT_NEAREST,
        INSERT_FARTHEST,
        INSERT_CHEAPEST,
        INSERT_RANDOM
      };     
      
      Cost run(InsertionMethod method = INSERT_FARTHEST) {
        switch (method) {
          case INSERT_NEAREST:
            start<InitTour<true>, NearestSelection, DefaultInsert>();
            break;
          case INSERT_FARTHEST:
            start<InitTour<false>, FarthestSelection, DefaultInsert>();
            break;
          case INSERT_CHEAPEST:
            start<InitTour<true>, CheapestSelection, CheapestInsert>();
            break;
          case INSERT_RANDOM:
            start<InitTour<true>, RandomSelection, DefaultInsert>();
            break;
        }
        return sum;
      }

      template <typename L>
      void tourNodes(L &container) {
        container(insertion_tsp_helper::vectorConvert<L>(nodesPath));
      }
      
      template <template <typename> class L>
      L<Node> tourNodes() {
        return insertion_tsp_helper::vectorConvert<L<Node> >(nodesPath);
      }
      
      const std::vector<Node>& tourNodes() {
        return nodesPath;
      }
      
      Path<FullGraph> tour() {
        Path<FullGraph> p;
        if (nodesPath.size()<2)
          return p;
        
        for (unsigned int i=0; i<nodesPath.size()-1; ++i)
          p.addBack(_gr.arc(nodesPath[i], nodesPath[i+1]));
        
        p.addBack(_gr.arc(nodesPath.back(), nodesPath.front()));
        return p;
      }
      
      Cost tourCost() {
        return sum;
      }


    private:

      template <bool MIN>
      class InitTour {
        public:
          InitTour(const FullGraph &gr, const CostMap &cost,
                   std::vector<Node> &used, std::vector<Node> &notused,
                   Cost &fullcost) :
              _gr(gr), _cost(cost), _used(used), _notused(notused),
              _fullcost(fullcost) {}

          std::vector<Node> init() const {
            Edge min = (MIN) ? mapMin(_gr, _cost) : mapMax(_gr, _cost);
            std::vector<Node> path;
            path.push_back(_gr.u(min));
            path.push_back(_gr.v(min));
            
            _used.clear();
            _notused.clear();
            for (NodeIt n(_gr); n!=INVALID; ++n) {
              if (n != _gr.u(min) && n != _gr.v(min)) {
                _notused.push_back(n);
              }
            }
            _used.push_back(_gr.u(min));
            _used.push_back(_gr.v(min));
            
            _fullcost = _cost[min]*2;
            return path;
          }

        private:
          const FullGraph &_gr;
          const CostMap &_cost;
          std::vector<Node> &_used;
          std::vector<Node> &_notused;
          Cost &_fullcost;
      };

      class NearestSelection {
        public:
          NearestSelection(const FullGraph &gr, const CostMap &cost,
                           std::vector<Node> &used, std::vector<Node> &notused) : 
              _gr(gr), _cost(cost), _used(used), _notused(notused) {}
      
          Node select() const {

            Cost insert_val =
              std::numeric_limits<Cost>::max();
            int insert_node = -1;
            
            for (unsigned int i=0; i<_notused.size(); ++i) {
              Cost min_val = _cost[_gr.edge(_notused[i], _used[0])];
              int min_node = 0;
            
              for (unsigned int j=1; j<_used.size(); ++j) {
                if (min_val > _cost[_gr.edge(_notused[i], _used[j])]) {
                    min_val = _cost[_gr.edge(_notused[i], _used[j])];
                    min_node = j;
                }
              }
              
              if (insert_val > min_val) {
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
          std::vector<Node> &_used;
          std::vector<Node> &_notused;
      };

      class FarthestSelection {
        public:
          FarthestSelection(const FullGraph &gr, const CostMap &cost,
                            std::vector<Node> &used, std::vector<Node> &notused) : 
              _gr(gr), _cost(cost), _used(used), _notused(notused) {}
      
          Node select() const {

            Cost insert_val =
              std::numeric_limits<Cost>::min();
            int insert_node = -1;
            
            for (unsigned int i=0; i<_notused.size(); ++i) {
              Cost min_val = _cost[_gr.edge(_notused[i], _used[0])];
              int min_node = 0;
              
              for (unsigned int j=1; j<_used.size(); ++j) {
                if (min_val > _cost[_gr.edge(_notused[i], _used[j])]) {
                  min_val = _cost[_gr.edge(_notused[i], _used[j])];
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
          std::vector<Node> &_used;
          std::vector<Node> &_notused;
      };


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
                            std::vector<Node> &used, std::vector<Node> &notused) : 
              _gr(gr), _cost(cost), _used(used), _notused(notused) {}
        
          Cost select() const {
            int insert_notused = -1;
            int best_insert_index = -1;
            Cost insert_val =
              std::numeric_limits<Cost>::max();
            
            for (unsigned int i=0; i<_notused.size(); ++i) {
              int min = i;
              int best_insert_tmp = 0;
              Cost min_val =
                costDiff(_used.front(), _used.back(), _notused[i]);
              
              for (unsigned int j=1; j<_used.size(); ++j) {
                Cost tmp =
                  costDiff(_used[j-1], _used[j], _notused[i]);

                if (min_val > tmp) {
                  min = i;
                  min_val = tmp;
                  best_insert_tmp = j;
                }
              }

              if (insert_val > min_val) {
                insert_notused = min;
                insert_val = min_val;
                best_insert_index = best_insert_tmp;
              }
            }

            _used.insert(_used.begin()+best_insert_index, _notused[insert_notused]);
            _notused.erase(_notused.begin()+insert_notused);

            return insert_val;
          }
          
        private:
          const FullGraph &_gr;
          const CostMap &_cost;
          std::vector<Node> &_used;
          std::vector<Node> &_notused;
      };

      class RandomSelection {
        public:
          RandomSelection(const FullGraph &, const CostMap &,
                            std::vector<Node> &, std::vector<Node> &notused) : 
                           _notused(notused) {
            rnd.seedFromTime();
          }
          
          Node select() const {
            const int index = rnd[_notused.size()];
            Node n = _notused[index];
            _notused.erase(_notused.begin()+index);
            return n;
          }
        private:
          std::vector<Node> &_notused;
      };


      class DefaultInsert {
        private:
          Cost costDiff(Node u, Node v, Node w) const {
            return 
              _cost[_gr.edge(u, w)] +
              _cost[_gr.edge(v, w)] -
              _cost[_gr.edge(u, v)];
          }
  
        public:
          DefaultInsert(const FullGraph &gr, const CostMap &cost,
                        std::vector<Node> &path, Cost &fullcost) : 
            _gr(gr), _cost(cost), _path(path), _fullcost(fullcost) {}
          
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
            _fullcost += min_val;
          }
        
        private:
          const FullGraph &_gr;
          const CostMap &_cost;
          std::vector<Node> &_path;
          Cost &_fullcost;
      };
  
      class CheapestInsert {
        TEMPLATE_GRAPH_TYPEDEFS(FullGraph);
        public:
          CheapestInsert(const FullGraph &, const CostMap &,
                         std::vector<Node> &, Cost &fullcost) : 
            _fullcost(fullcost) {}
          
          void insert(Cost diff) const {
            _fullcost += diff;
          }

        private:
          Cost &_fullcost;
      };  
      
      
      template <class InitFunctor, class SelectionFunctor, class InsertFunctor>
      void start() {
        InitFunctor init(_gr, _cost, nodesPath, notused, sum);
        SelectionFunctor selectNode(_gr, _cost, nodesPath, notused);
        InsertFunctor insertNode(_gr, _cost, nodesPath, sum);
        
        nodesPath = init.init();
        
        for (int i=0; i<_gr.nodeNum()-2; ++i) {
          insertNode.insert(selectNode.select());
        }
        
        sum = _cost[ _gr.edge(nodesPath.front(), nodesPath.back()) ];
        for (unsigned int i=0; i<nodesPath.size()-1; ++i)
          sum += _cost[ _gr.edge(nodesPath[i], nodesPath[i+1]) ];
      }
    
      const FullGraph &_gr;
      const CostMap &_cost;
      std::vector<Node> notused;
      std::vector<Node> nodesPath;
      Cost sum;
  };
  
}; // namespace lemon

#endif
