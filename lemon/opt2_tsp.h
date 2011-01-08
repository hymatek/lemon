#ifndef LEMON_OPT2_TSP_H
#define LEMON_OPT2_TSP_H

#include <vector>
#include <lemon/full_graph.h>
#include <lemon/path.h>

namespace lemon {
  
  namespace opt2_helper {
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
  class Opt2Tsp {
    private:
      GRAPH_TYPEDEFS(FullGraph);

    public:
      typedef CM CostMap;
      typedef typename CM::Value Cost;
      
    
      Opt2Tsp(const FullGraph &gr, const CostMap &cost) : _gr(gr), _cost(cost), 
            tmppath(_gr.nodeNum()*2) {
            
        for (int i=1; i<_gr.nodeNum()-1; ++i) {
          tmppath[2*i] = i-1;
          tmppath[2*i+1] = i+1;
        }
        tmppath[0] = _gr.nodeNum()-1;
        tmppath[1] = 1;
        tmppath[2*(_gr.nodeNum()-1)] = _gr.nodeNum()-2;
        tmppath[2*(_gr.nodeNum()-1)+1] = 0;
      }
      
      Opt2Tsp(const FullGraph &gr, const CostMap &cost, const std::vector<Node> &path) : 
              _gr(gr), _cost(cost), tmppath(_gr.nodeNum()*2) {

        for (unsigned int i=1; i<path.size()-1; ++i) {
          tmppath[2*_gr.id(path[i])] = _gr.id(path[i-1]);
          tmppath[2*_gr.id(path[i])+1] = _gr.id(path[i+1]);
        }
        
        tmppath[2*_gr.id(path[0])] = _gr.id(path.back());
        tmppath[2*_gr.id(path[0])+1] = _gr.id(path[1]);
        tmppath[2*_gr.id(path.back())] = _gr.id(path[path.size()-2]);
        tmppath[2*_gr.id(path.back())+1] = _gr.id(path.front());
      }

    private:
      Cost c(int u, int v) {
        return _cost[_gr.edge(_gr.nodeFromId(u), _gr.nodeFromId(v))];
      }
      
      class It {
        public:
          It(const std::vector<int> &path, int i=0) : tmppath(path), act(i), last(tmppath[2*act]) {}
          It(const std::vector<int> &path, int i, int l) : tmppath(path), act(i), last(l) {}

          int next_index() const {
            return (tmppath[2*act]==last)? 2*act+1 : 2*act;
          }
          
          int prev_index() const {
            return (tmppath[2*act]==last)? 2*act : 2*act+1;
          }
          
          int next() const {
            return tmppath[next_index()];
          }

          int prev() const {
            return tmppath[prev_index()];
          }
          
          It& operator++() {
            int tmp = act;
            act = next();
            last = tmp;
            return *this;
          }
          
          operator int() const {
            return act;
          }
          
        private:
          const std::vector<int> &tmppath;
          int act;
          int last;
      };

      bool check(std::vector<int> &path, It i, It j) {
        if (c(i, i.next()) + c(j, j.next()) > 
            c(i, j) + c(j.next(), i.next())) {

            path[ It(path, i.next(), i).prev_index() ] = j.next();
            path[ It(path, j.next(), j).prev_index() ] = i.next();

            path[i.next_index()] = j;
            path[j.next_index()] = i;

            return true;
        }
        return false;
      }
      
    public:
      
      Cost run() {
        _path.clear();

        if (_gr.nodeNum()>3) {

opt2_tsp_label:
          It i(tmppath);
          It j(tmppath, i, i.prev());
          ++j; ++j;
          for (; j.next()!=0; ++j) {
            if (check(tmppath, i, j))
              goto opt2_tsp_label;
          }
          
          for (++i; i.next()!=0; ++i) {
            It j(tmppath, i, i.prev());
            if (++j==0)
              break;
            if (++j==0)
              break;
            
            for (; j!=0; ++j) {
              if (check(tmppath, i, j))
                goto opt2_tsp_label;
            }
          }
        }

        It k(tmppath);
        _path.push_back(_gr.nodeFromId(k));
        for (++k; k!=0; ++k)
          _path.push_back(_gr.nodeFromId(k));

        

        _sum = _cost[ _gr.edge(_path.back(), _path.front()) ];
        for (unsigned int i=0; i<_path.size()-1; ++i)
          _sum += _cost[ _gr.edge(_path[i], _path[i+1]) ];
        return _sum;
      }

      
      template <typename L>
      void tourNodes(L &container) {
        container(opt2_helper::vectorConvert<L>(_path));
      }
      
      template <template <typename> class L>
      L<Node> tourNodes() {
        return opt2_helper::vectorConvert<L<Node> >(_path);
      }

      const std::vector<Node>& tourNodes() {
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
    std::vector<int> tmppath;
    std::vector<Node> _path;
  };


}; // namespace lemon

#endif
