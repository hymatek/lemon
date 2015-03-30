/* -*- mode: C++; indent-tabs-mode: nil; -*-
 *
 * This file is a part of LEMON, a generic C++ optimization library.
 *
 * Copyright (C) 2015
 * EMAXA Kutato-fejleszto Kft. (EMAXA Research Ltd.)
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

#ifndef LEMON_VF2_H
#define LEMON_VF2_H

#include <lemon/core.h>
#include <lemon/concepts/graph.h>
#include <lemon/dfs.h>
#include <lemon/bfs.h>
#include <test/test_tools.h>

#include <vector>
#include <set>

namespace lemon
{
  namespace bits
  {
    namespace vf2
    {
      class AlwaysEq
      {
      public:
        template<class T1, class T2>
        bool operator()(T1, T2) const
        {
          return true;
        }
      };

      template<class M1, class M2>
      class MapEq
      {
        const M1 &_m1;
        const M2 &_m2;
      public:
        MapEq(const M1 &m1, const M2 &m2) : _m1(m1), _m2(m2) {}
        bool operator()(typename M1::Key k1, typename M2::Key k2) const
        {
          return _m1[k1] == _m2[k2];
        }
      };

      template <class G>
      class DfsLeaveOrder : public DfsVisitor<G>
      {
        const G &_g;
        std::vector<typename G::Node> &_order;
        int i;
      public:
        DfsLeaveOrder(const G &g, std::vector<typename G::Node> &order)
          : i(countNodes(g)), _g(g), _order(order)
        {}
        void leave(const typename G::Node &node)
        {
          _order[--i]=node;
        }
      };

      template <class G>
      class BfsLeaveOrder : public BfsVisitor<G>
      {
        int i;
        const G &_g;
        std::vector<typename G::Node> &_order;
      public:
        BfsLeaveOrder(const G &g, std::vector<typename G::Node> &order)
          : i(0), _g(g), _order(order)
        {}
        void process(const typename G::Node &node)
        {
          _order[i++]=node;
        }
      };
    }
  }

  enum MappingType {
    SUBGRAPH = 0,
    INDUCED = 1,
    ISOMORPH = 2
  };

  template<class G1, class G2, class I, class NEq = bits::vf2::AlwaysEq >
  class Vf2
  {
    //Current depth in the DFS tree.
    int _depth;
    //Functor with bool operator()(G1::Node,G2::Node), which returns 1
    //if and only if the 2 nodes are equivalent.
    NEq _nEq;

    typename G2::template NodeMap<int> _conn;
    //Current matching. We index it by the nodes of g1, and match[v] is
    //a node of g2.
    I &_match;
    //order[i] is the node of g1, for which we find a pair if depth=i
    std::vector<typename G1::Node> order;
    //currEdgeIts[i] is an edge iterator, witch is last used in the ith
    //depth to find a pair for order[i].
    std::vector<typename G2::IncEdgeIt> currEdgeIts;
    //The small graph.
    const G1 &_g1;
    //The big graph.
    const G2 &_g2;
    //lookup tables for cut the searchtree
    typename G1::template NodeMap<int> rNew1t,rInOut1t;

    MappingType _mapping_type;

    //cut the search tree
    template<MappingType MT>
    bool cut(const typename G1::Node n1,const typename G2::Node n2) const
    {
      int rNew2=0,rInOut2=0;
      for(typename G2::IncEdgeIt e2(_g2,n2); e2!=INVALID; ++e2)
        {
          const typename G2::Node currNode=_g2.oppositeNode(n2,e2);
          if(_conn[currNode]>0)
            ++rInOut2;
          else if(MT!=SUBGRAPH&&_conn[currNode]==0)
            ++rNew2;
        }
      switch(MT)
        {
        case INDUCED:
          return rInOut1t[n1]<=rInOut2&&rNew1t[n1]<=rNew2;
        case SUBGRAPH:
          return rInOut1t[n1]<=rInOut2;
        case ISOMORPH:
          return rInOut1t[n1]==rInOut2&&rNew1t[n1]==rNew2;
        default:
          return false;
        }
    }

    template<MappingType MT>
    bool feas(const typename G1::Node n1,const typename G2::Node n2)
    {
      if(!_nEq(n1,n2))
        return 0;

      for(typename G1::IncEdgeIt e1(_g1,n1); e1!=INVALID; ++e1)
        {
          const typename G1::Node currNode=_g1.oppositeNode(n1,e1);
          if(_match[currNode]!=INVALID)
            --_conn[_match[currNode]];
        }
      bool isIso=1;
      for(typename G2::IncEdgeIt e2(_g2,n2); e2!=INVALID; ++e2)
        {
          const typename G2::Node currNode=_g2.oppositeNode(n2,e2);
          if(_conn[currNode]<-1)
            ++_conn[currNode];
          else if(MT!=SUBGRAPH&&_conn[currNode]==-1)
            {
              isIso=0;
              break;
            }
        }

      for(typename G1::IncEdgeIt e1(_g1,n1); e1!=INVALID; ++e1)
        {
          const typename G1::Node currNode=_g1.oppositeNode(n1,e1);
          if(_match[currNode]!=INVALID&&_conn[_match[currNode]]!=-1)
            {
              switch(MT)
                {
                case INDUCED:
                case ISOMORPH:
                  isIso=0;
                  break;
                case SUBGRAPH:
                  if(_conn[_match[currNode]]<-1)
                    isIso=0;
                  break;
                }
              _conn[_match[currNode]]=-1;
            }
        }
      return isIso&&cut<MT>(n1,n2);
    }

    void addPair(const typename G1::Node n1,const typename G2::Node n2)
    {
      _conn[n2]=-1;
      _match.set(n1,n2);
      for(typename G2::IncEdgeIt e2(_g2,n2); e2!=INVALID; ++e2)
        if(_conn[_g2.oppositeNode(n2,e2)]!=-1)
          ++_conn[_g2.oppositeNode(n2,e2)];
    }

    void subPair(const typename G1::Node n1,const typename G2::Node n2)
    {
      _conn[n2]=0;
      _match.set(n1,INVALID);
      for(typename G2::IncEdgeIt e2(_g2,n2); e2!=INVALID; ++e2)
        {
          const typename G2::Node currNode=_g2.oppositeNode(n2,e2);
          if(_conn[currNode]>0)
            --_conn[currNode];
          else if(_conn[currNode]==-1)
            ++_conn[n2];
        }
    }

    void setOrder()//we will find pairs for the nodes of g1 in this order
    {
      // bits::vf2::DfsLeaveOrder<G1> v(_g1,order);
      //   DfsVisit<G1,bits::vf2::DfsLeaveOrder<G1> >dfs(_g1, v);
      //   dfs.run();

      //it is more efficient in practice than DFS
      bits::vf2::BfsLeaveOrder<G1> v(_g1,order);
      BfsVisit<G1,bits::vf2::BfsLeaveOrder<G1> >bfs(_g1, v);
      bfs.run();
    }

  public:

    template<MappingType MT>
    bool extMatch()
    {
      while(_depth>=0)
        {
          //there are not nodes in g1, which has not pair in g2.
          if(_depth==static_cast<int>(order.size()))
            {
              --_depth;
              return true;
            }
          //the node of g2, which neighbours are the candidates for
          //the pair of order[_depth]
          typename G2::Node currPNode;
          if(currEdgeIts[_depth]==INVALID)
            {
              typename G1::IncEdgeIt fstMatchedE(_g1,order[_depth]);
              //if _match[order[_depth]]!=INVALID, we dont use
              //fstMatchedE
              if(_match[order[_depth]]==INVALID)
                for(; fstMatchedE!=INVALID &&
                      _match[_g1.oppositeNode(order[_depth],
                                              fstMatchedE)]==INVALID;
                    ++fstMatchedE) ; //find fstMatchedE
              if(fstMatchedE==INVALID||_match[order[_depth]]!=INVALID)
                {
                  //We did not find an covered neighbour, this means
                  //the graph is not connected(or _depth==0).  Every
                  //uncovered(and there are some other properties due
                  //to the spec. problem types) node of g2 is
                  //candidate.  We can read the iterator of the last
                  //tryed node from the match if it is not the first
                  //try(match[order[_depth]]!=INVALID)
                  typename G2::NodeIt n2(_g2);
                  //if its not the first try
                  if(_match[order[_depth]]!=INVALID)
                    {
                      n2=++typename G2::NodeIt(_g2,_match[order[_depth]]);
                      subPair(order[_depth],_match[order[_depth]]);
                    }
                  for(; n2!=INVALID; ++n2)
                    if(MT!=SUBGRAPH&&_conn[n2]==0)
                      {
                        if(feas<MT>(order[_depth],n2))
                          break;
                      }
                    else if(MT==SUBGRAPH&&_conn[n2]>=0)
                      if(feas<MT>(order[_depth],n2))
                        break;
                  // n2 is the next candidate
                  if(n2!=INVALID)
                    {
                      addPair(order[_depth],n2);
                      ++_depth;
                    }
                  else // there is no more candidate
                    --_depth;
                  continue;
                }
              else
                {
                  currPNode=_match[_g1.oppositeNode(order[_depth],fstMatchedE)];
                  currEdgeIts[_depth]=typename G2::IncEdgeIt(_g2,currPNode);
                }
            }
          else
            {
              currPNode=_g2.oppositeNode(_match[order[_depth]],
                                         currEdgeIts[_depth]);
              subPair(order[_depth],_match[order[_depth]]);
              ++currEdgeIts[_depth];
            }
          for(; currEdgeIts[_depth]!=INVALID; ++currEdgeIts[_depth])
            {
              const typename G2::Node currNode =
                _g2.oppositeNode(currPNode, currEdgeIts[_depth]);
              //if currNode is uncovered
              if(_conn[currNode]>0&&feas<MT>(order[_depth],currNode))
                {
                  addPair(order[_depth],currNode);
                  break;
                }
            }
          currEdgeIts[_depth]==INVALID?--_depth:++_depth;
        }
      return false;
    }

    //calc. the lookup table for cut the searchtree
    void setRNew1tRInOut1t()
    {
      typename G1::template NodeMap<int> tmp(_g1,0);
      for(unsigned int i=0; i<order.size(); ++i)
        {
          tmp[order[i]]=-1;
          for(typename G1::IncEdgeIt e1(_g1,order[i]); e1!=INVALID; ++e1)
            {
              const typename G1::Node currNode=_g1.oppositeNode(order[i],e1);
              if(tmp[currNode]>0)
                ++rInOut1t[order[i]];
              else if(tmp[currNode]==0)
                ++rNew1t[order[i]];
            }
          for(typename G1::IncEdgeIt e1(_g1,order[i]); e1!=INVALID; ++e1)
            {
              const typename G1::Node currNode=_g1.oppositeNode(order[i],e1);
              if(tmp[currNode]!=-1)
                ++tmp[currNode];
            }
        }
    }
  public:
    Vf2(const G1 &g1, const G2 &g2,I &i, const NEq &nEq = NEq() ) :
      _nEq(nEq), _conn(g2,0), _match(i), order(countNodes(g1)),
      currEdgeIts(countNodes(g1),INVALID), _g1(g1), _g2(g2), rNew1t(g1,0),
      rInOut1t(g1,0), _mapping_type(SUBGRAPH)
    {
      _depth=0;
      setOrder();
      setRNew1tRInOut1t();
    }

    MappingType mappingType() const { return _mapping_type; }
    void mappingType(MappingType m_type) { _mapping_type = m_type; }

    bool find()
    {
      switch(_mapping_type)
        {
        case SUBGRAPH:
          return extMatch<SUBGRAPH>();
        case INDUCED:
          return extMatch<INDUCED>();
        case ISOMORPH:
          return extMatch<ISOMORPH>();
        default:
          return false;
        }
    }
  };

  template<class G1, class G2>
  class Vf2WizardBase
  {
  protected:
    typedef G1 Graph1;
    typedef G2 Graph2;

    const G1 &_g1;
    const G2 &_g2;

    MappingType _mapping_type;

    typedef typename G1::template NodeMap<typename G2::Node> Mapping;
    bool _local_mapping;
    void *_mapping;
    void createMapping()
    {
      _mapping = new Mapping(_g1);
    }

    typedef bits::vf2::AlwaysEq NodeEq;
    NodeEq _node_eq;

    Vf2WizardBase(const G1 &g1,const G2 &g2)
      : _g1(g1), _g2(g2), _mapping_type(SUBGRAPH), _local_mapping(true) {}
  };

  template<class TR>
  class Vf2Wizard : public TR
  {
    typedef TR Base;
    typedef typename TR::Graph1 Graph1;
    typedef typename TR::Graph2 Graph2;

    typedef typename TR::Mapping Mapping;
    typedef typename TR::NodeEq NodeEq;

    using TR::_g1;
    using TR::_g2;
    using TR::_mapping_type;
    using TR::_mapping;
    using TR::_node_eq;

  public:
    ///Copy constructor
    Vf2Wizard(const Graph1 &g1,const Graph2 &g2) : Base(g1,g2) {}

    ///Copy constructor
    Vf2Wizard(const Base &b) : Base(b) {}


    template<class T>
    struct SetMappingBase : public Base {
      typedef T Mapping;
      SetMappingBase(const Base &b) : Base(b) {}
    };

    ///\brief \ref named-templ-param "Named parameter" for setting
    ///the mapping.
    ///
    ///\ref named-templ-param "Named parameter" function for setting
    ///the map that stores the found embedding.
    template<class T>
    Vf2Wizard< SetMappingBase<T> > mapping(const T &t)
    {
      Base::_mapping=reinterpret_cast<void*>(const_cast<T*>(&t));
      Base::_local_mapping = false;
      return Vf2Wizard<SetMappingBase<T> >(*this);
    }

    template<class NE>
    struct SetNodeEqBase : public Base {
      typedef NE NodeEq;
      NodeEq _node_eq;
      SetNodeEqBase(const Base &b, const NE &node_eq)
        : Base(b), _node_eq(node_eq) {}
    };

    ///\brief \ref named-templ-param "Named parameter" for setting
    ///the node equivalence relation.
    ///
    ///\ref named-templ-param "Named parameter" function for setting
    ///the equivalence relation between the nodes.
    template<class T>
    Vf2Wizard< SetNodeEqBase<T> > nodeEq(const T &node_eq)
    {
      return Vf2Wizard<SetNodeEqBase<T> >(SetNodeEqBase<T>(*this,node_eq));
    }

    ///\brief \ref named-templ-param "Named parameter" for setting
    ///the node labels.
    ///
    ///\ref named-templ-param "Named parameter" function for setting
    ///the node labels defining equivalence relation between them.
    template<class M1, class M2>
    Vf2Wizard< SetNodeEqBase<bits::vf2::MapEq<M1,M2> > >
    nodeLabels(const M1 &m1,const M2 &m2)
    {
      return nodeEq(bits::vf2::MapEq<M1,M2>(m1,m2));
    }

    Vf2Wizard<Base> &mappingType(MappingType m_type)
    {
      _mapping_type = m_type;
      return *this;
    }

    Vf2Wizard<Base> &induced()
    {
      _mapping_type = INDUCED;
      return *this;
    }

    Vf2Wizard<Base> &iso()
    {
      _mapping_type = ISOMORPH;
      return *this;
    }

    bool run()
    {
      if(Base::_local_mapping)
        Base::createMapping();

      Vf2<Graph1, Graph2, Mapping, NodeEq >
        alg(_g1, _g2, *reinterpret_cast<Mapping*>(_mapping), _node_eq);

      alg.mappingType(_mapping_type);

      bool ret = alg.find();

      if(Base::_local_mapping)
        delete reinterpret_cast<Mapping*>(_mapping);

      return ret;
    }
  };

  template<class G1, class G2>
  Vf2Wizard<Vf2WizardBase<G1,G2> > vf2(const G1 &g1, const G2 &g2)
  {
    return Vf2Wizard<Vf2WizardBase<G1,G2> >(g1,g2);
  }

}

#endif
