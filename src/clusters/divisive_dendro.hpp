#ifndef divisive_dendro_hpp_defined
#define divisive_dendro_hpp_defined

#include <queue>
#include <set>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/property_map/property_map.hpp>
#include <iostream> 
namespace clusters {

  template <typename V, typename T>
  struct dendro_vertex {
    T level;                    // the lowest density level, at which the cluster enters into the hierarchy. Should we specify also at which it leaves?
    std::set<V> elts;           // set of data element indices in a given cluster.
  };

  /**
   * @brief hierarchy graph type traits
   * 
   */
  template<typename V, typename T>
  struct dendro_graph_traits {
    typedef dendro_vertex<V,T> vertex_type;
    typedef boost::adjacency_list
    <boost::listS               // vertices are listS
     , boost::listS
     , boost::directedS
     , dendro_vertex<V,T>       // vertex type
     , boost::no_property
     , boost::listS
     > graph_type;
    typedef typename boost::graph_traits<graph_type>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<graph_type>::edge_descriptor edge_descriptor;
  };

  
  /** 
   * 
   * @note What do we need from data_graph: 1. copy-constructible, has vertex_descriptor, edge_descriptor, vertices(), edges(), edge_index_t (!that's important). So, basically, a boost graph with edge_index_t; the rest is minimum graph_traits. 
   * @note Dendro_graph is a bit more specific that data_graph. For example, it relies on a special vertex type, called dendro_vertex. That is because we are doing something quite specific here: collecting the elements of data into sets, and specifying a numeric value (pair of values?) that are specific for the collection. So it makes sense to define a type.
   *
   * @param[in] e_beg iterator to start of the edge mass data. It must be an iterator passable to the boost::iterator_property_map. @note Constraint: the range of e_beg should be at least as large as boost::num_edges(G).
   * @param G 
   * @param Kmin 
   * 
   * @return 
   */
  template
  <class Iterator,
   class data_graph,
   class dendro_graph = typename dendro_graph_traits<typename boost::graph_traits<data_graph>::vertex_descriptor,typename Iterator::value_type>::graph_type
   >
  dendro_graph
  divisive_dendro(const Iterator e_beg, const data_graph& G, const int Kmin) // : input
  {

    data_graph Gc(G);
    dendro_graph D;

    typedef typename Iterator::value_type real_type;
    typedef typename boost::graph_traits<data_graph>::vertex_descriptor density_vertex_index;
    typedef typename boost::graph_traits<data_graph>::edge_descriptor   density_edge_index;
    typedef typename boost::graph_traits<dendro_graph>::vertex_descriptor vertex_descriptor;
    
    /*
     * brief Extract connected component indices. 
     *
     * This is a depth first visitor that we use to collect the indices
     * of data nodes that are within a connected component of a graph,
     * starting from a given node.
     *
     * param G data graph type
     *
     */
    struct component_extractor : boost::default_dfs_visitor {
      component_extractor(std::set<density_vertex_index> &elts) : m_elts(elts) {}
      // add a discovered vertex;
      void discover_vertex(const density_vertex_index u, const data_graph& g) {
        m_elts.insert(u);
      }
    private:
      std::set<density_vertex_index>& m_elts;
    };
    
    // --------------------------------------------------------------------------------------------
    auto vv=boost::vertices(Gc);
    auto ee=boost::edges(Gc);
    // --------------------------------------------------------------------------------------------
    std::vector<density_edge_index> eind(ee.first,ee.second);
    //    RealEdgeMap density_edge_map(e_beg, boost::get(boost::edge_index,Gc));
    auto density_edge_map = make_iterator_property_map(e_beg, boost::get(boost::edge_index,Gc));

    // for(auto ee=boost::edges(Gc);ee.first!=ee.second;ee.first++) {
    //   std::cout << density_edge_map[*ee.first] << "\n";
    // }

    auto edge_order_lambda = [&density_edge_map](density_edge_index const& lhs, density_edge_index const& rhs) -> bool { return (density_edge_map[lhs]>=density_edge_map[rhs]); };
    
    std::priority_queue
      <density_edge_index
       , std::vector<density_edge_index>
       , decltype(edge_order_lambda)
       > Q(edge_order_lambda, eind);

    // ---------------------------------------------------------------------------------------------
    std::vector<vertex_descriptor> owner_data(boost::num_vertices(Gc));
    std::vector<boost::default_color_type> color_data(boost::num_vertices(Gc));
    auto color_map = boost::make_iterator_property_map(color_data.begin(), get(boost::vertex_index, Gc));
    auto owner_map = boost::make_iterator_property_map(owner_data.begin(), get(boost::vertex_index, Gc));
    //       OwnerMap owner_map(owner_data.begin(), boost::get(boost::vertex_index,Gc));
    //============================================================================================
    // let's add everyting to the root vertex:
    vertex_descriptor v
      , u=boost::add_vertex(D);
    //    D[u]=std::set<density_vertex_index>(vv.first,vv.second);
    D[u].elts=std::set<density_vertex_index>(vv.first,vv.second);
    D[u].level=0;
    for(auto &i : owner_data) i=u;
    // MASS[u]=riemann_sum(vertex_mass_map,H[u].elts,0.0);
    // ===========================================================================================
    //    std::cerr << "priority ordering:\n";
    std::set<density_vertex_index> elts;
    while(!Q.empty()) {
      // std::cerr << "Before edge processing: V=" << boost::num_vertices(H)
      //           << " E=" << boost::num_edges(H) << "\n";
      const auto& q = Q.top();
      //      std::cerr << "first edge to consider: E" << q << " w=" << density_edge_map[q] << "\n";
      density_vertex_index s = boost::source(q,Gc);
      density_vertex_index t = boost::target(q,Gc);
      real_type z = density_edge_map[q];
      boost::remove_edge(q,Gc);  // this modifies G! That's why we're using Gc
      // The "s" (source) branch
      u=owner_map[s];
      std::fill(color_data.begin(),color_data.end(),boost::white_color);
      elts.clear();
      boost::depth_first_visit
        (Gc
         , s
         , component_extractor(elts)
         , color_map
         );
      // assign to v if large, enough, assign to NULL otherwise
      if(u!=NULL && elts.size()>=size_t(Kmin)) {
        v=boost::add_vertex(D);
        //        D[v]=std::set<density_vertex_index>(elts);
        D[v].elts=std::set<density_vertex_index>(elts);
        D[v].level=z;
        boost::add_edge(u,v,D);
        for(auto &i: elts) owner_map[i]=v;
        //        MASS[v]=riemann_sum(vertex_mass_map,elts,z);
      } else if(u!=NULL) for(auto &i: elts) owner_map[i]=NULL;
      // --------------------------------------------------------------
      // The "t" (target) branch
      u=owner_map[t];
      std::fill(color_data.begin(),color_data.end(),boost::white_color);
      elts.clear();
      depth_first_visit
        (Gc
         , t
         , component_extractor(elts)
         , color_map
         );
      // assign to v if large, enough, assign to NULL otherwise
      if(u!=NULL && elts.size()>=size_t(Kmin)) {
        v=boost::add_vertex(D);
        // D[v]=std::set<density_vertex_index>(elts);
        D[v].elts=std::set<density_vertex_index>(elts);
        D[v].level=z;
        boost::add_edge(u,v,D);
        for(auto &i: elts) owner_map[i]=v;
        //        MASS[v]=riemann_sum(vertex_mass_map,elts,z);
      } else if(u!=NULL) for(auto& i: elts) owner_map[i]=NULL;
      Q.pop();
    }
    // ---------------------------------------------------------------
    // concatenation of thin links' step
    bool concat=true;
    while(concat) {
      concat=false;
      //      vv=boost::vertices(H); // <= assuming this is valid
      for(auto vv=boost::vertices(D);(!concat) && (vv.first!=vv.second);vv.first++) {
        auto const& o=*vv.first;
        if(boost::out_degree(o,D)==1) { // concatenate
          concat=true;
          auto const t = boost::target(*out_edges(o,D).first,D);
          // this is done zero times is out_degree(v)=0
          for(auto ee=boost::out_edges(t,D);ee.first!=ee.second;++ee.first)
            if(boost::add_edge(o,boost::target(*ee.first,D),D).second==false) throw std::runtime_error("concatenation: failed to add edge\n");
          //          MASS.erase(MASS.find(t)); // assuming find does not return NULL.
          clear_vertex(t,D);        // may not be necessary (because we took care to eliminate out edges of t), but we do it anyway.
          remove_vertex(t,D);
        }
      }
    }
    //    return outcome; // {H,MASS}
    return D;
  } // divisive_dendro


  /** 
   * A snippet that is used a lot.
   * 
   * @param u 
   * @param G 
   * @param m 
   * 
   * @return 
   */
  template <class Graph, class vertex, class Map, typename K = typename Map::mapped_type>
  K self_minus_children(vertex u, Graph& G, Map m) {
    K res = m.at(u);
    for(auto e=boost::out_edges(u,G);e.first!=e.second;e.first++) res -= m.at(boost::target(*e.first,G));
    return res;
  }


  template <class graph_type, // = typename dendro_graph_traits<typename boost::graph_traits<data_graph>::vertex_descriptor,typename Iterator::value_type>::graph_type,
            class MassMap,
            typename real_type = typename MassMap::mapped_type
            >
  graph_type
  max_contrast_transform(const graph_type& H, const MassMap& mass)
  {
    typedef typename boost::graph_traits<graph_type>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<graph_type>::edge_descriptor edge_descriptor;
    // design decision: we'll use a std::map for mapping vertices (see below)
    // The following typedefs are the result of this decision. If we replace map with something else, these also should be modified.
    typedef std::map<vertex_descriptor,size_t> index_map_container;
    typedef boost::associative_property_map<index_map_container> IndexMap;
    typedef boost::iterator_property_map<typename std::vector<vertex_descriptor>::iterator,IndexMap> VertexMap;
    typedef boost::iterator_property_map<typename std::vector<bool>::iterator, IndexMap> BoolMap;
    graph_type Hc(H);           // a copy
    // {
    //   auto a=boost::vertices(H);
    //   auto b=boost::vertices(Hc);
    //   for(;a.first!=a.second;a.first++,b.first++) {
    //     auto h = H[*a.first];
    //     auto hc= Hc[*b.first];
    //     std::cout << h.level << "," << hc.level
    //               << "\n";
    //   }
    // }
    // this is the aforementioned decision:
    index_map_container index; // need to fill it
    {
      size_t idx=0;
      for(auto vv=boost::vertices(Hc);vv.first!=vv.second;vv.first++,idx++) index[*vv.first]=idx;
    }
    bool is_reduced = false;
    // containers....
    std::vector<real_type> contrast_data(boost::num_vertices(Hc)); // here we'll store relative masses; the total mass should remain invariate!
    std::vector<bool> leaf_data(boost::num_vertices(Hc));
    std::vector<vertex_descriptor> predecessor_data(boost::num_vertices(Hc));
    std::vector<boost::default_color_type> color_data(boost::num_vertices(Hc));
    // and their property maps
    IndexMap index_map(index);
    auto predecessor_map = boost::make_iterator_property_map(predecessor_data.begin(), index_map); // actually a VertexMap
    auto leaf_map = boost::make_iterator_property_map(leaf_data.begin(),index_map); // actually a BoolMap
    auto color_map = boost::make_iterator_property_map(color_data.begin(), index_map);
    auto contrast_map = boost::make_iterator_property_map(contrast_data.begin(),index_map);
    // marker selects: leaves and predecessors
    struct leaf_marker : public boost::default_dfs_visitor
    {
      typedef typename std::list<vertex_descriptor> VertexListType;
      leaf_marker(const graph_type& g, VertexListType& l, VertexMap& p, BoolMap& s) : m_graph(g), m_list(l), m_predecessor(p), m_select(s) {}
      void start_vertex(const vertex_descriptor u, const graph_type& g) const {
        //        m_predecessor[u]=u;
        boost::put(m_predecessor,u,u);
      }
      void tree_edge(const edge_descriptor e, graph_type const& g) const {
        m_predecessor[boost::target(e,g)]=boost::source(e,g);
      }
      void discover_vertex(const vertex_descriptor u, graph_type const& g) const {
        if(boost::out_degree(u,g)==0) {
          m_list.push_back(u);
        }
      }
      void finish_vertex(const vertex_descriptor u, const graph_type& g) const {
        m_select[u]=true;
      }
    private:
      graph_type const& m_graph;
      VertexListType& m_list;
      VertexMap& m_predecessor;
      BoolMap& m_select;
    };
    std::list<vertex_descriptor> leaf;
    std::set<vertex_descriptor> S; // we use set to avoid duplicates
    //    size_t idx=0;
    for(auto uu=vertices(H),c=vertices(Hc);uu.first!=uu.second;++uu.first,++c.first) {
      auto& u=*uu.first;
      // real_type tmp = mass.at(u); // <== const mass, no const mass ==> mass[u];
      // for(auto ee=out_edges(u,H);ee.first!=ee.second;ee.first++) tmp -= mass.at(boost::target(*ee.first,H));
      real_type tmp = self_minus_children(u,H,mass);
      contrast_map[*c.first]=tmp;
    }
    // do we assume that the first vertex is the "root"?
    // the main loop
    while(true) {
      std::fill(color_data.begin(),color_data.end(),boost::white_color);
      std::fill(leaf_data.begin(),leaf_data.end(),false);
      leaf.clear();
      S.clear();
      auto vv=boost::vertices(Hc);
      is_reduced=true;
      boost::depth_first_visit
        (Hc
         , *vv.first
         , leaf_marker(H,leaf,predecessor_map,leaf_map)
         , color_map
         );
      // end condition: if all predecessors are equal to predecessors of the root, then get out
      for(auto const& ref=predecessor_map[*vv.first];is_reduced && (vv.first!=vv.second);++vv.first)
        is_reduced = is_reduced && (ref == predecessor_map[*vv.first]);
      if(is_reduced) break;
      // verify leaves: if all children of a parent are leaves, push it to S.
      for(auto li=leaf.begin();li!=leaf.end();++li) {
        auto const&l = *li;
        if(boost::get(leaf_map,l)) {
          bool p_terminal = true;
          auto& p = predecessor_map[l];
          for(auto ee=boost::out_edges(p,Hc); p_terminal && (ee.first!=ee.second); ++ee.first) {
            auto const t=boost::target(*ee.first,Hc);
            p_terminal = (boost::out_degree(t,Hc)==0);
            leaf_map[t]=false; // to avoid returning to already verified leaves;
          }
          if(p_terminal) S.insert(p);
        }
      }

      // the pruning step:
      // here we are removing light leaves according to the maximization algo.
      for(auto& s: S) {
        real_type rel = 0;
        for(auto ee=boost::out_edges(s,Hc);ee.first!=ee.second;++ee.first) rel += contrast_map[target(*ee.first,Hc)];
        auto ee=boost::out_edges(s,Hc);
        if (contrast_map[s] > rel) {   // parent is stronger
          for(auto next=ee.first;ee.first!=ee.second;ee.first=next) {
            ++next;
            auto t=boost::target(*ee.first,Hc);
            boost::remove_edge(*ee.first,Hc);
            boost::clear_vertex(t,Hc); // doubt if that's needed...
            boost::remove_vertex(t,Hc); // this is safe for VertexList=listS
          }
        } else {                // children are stronger
          auto& p=predecessor_map[s];
          for(auto next=ee.first;ee.first!=ee.second;ee.first=next) {
            ++next;
            auto t=boost::target(*ee.first,Hc);
            boost::add_edge(p,t,Hc); // should be safe because it is outside the range of [ee.first,ee.second)
            boost::remove_edge(*ee.first,Hc);
          }
          remove_edge(p,s,Hc);
          clear_vertex(s,Hc);    // not sure it's needed
          remove_vertex(s,Hc); // safe for VertexList=listS
        }
      }
    }
    return Hc;
  } // max_contrast_reduction

  template <typename V, typename T>
  std::set<V> parent_minus_children(typename dendro_graph_traits<V,T>::vertex_descriptor const u, typename dendro_graph_traits<V,T>::graph_type const& H)
  {
    std::set<V> source(H[u].elts);
    std::set<V> target;
    std::set<V> *s=&source;
    std::set<V> *t=&target;
    for(auto ee=boost::out_edges(u,H); ee.first!=ee.second;ee.first++) {
      auto &m = H[boost::target(*ee.first,H)].elts;
      std::set_difference((*s).begin(),(*s).end(),m.begin(),m.end(),std::inserter((*t),(*t).begin()));
      std::swap(s,t);
      (*t).clear();
    }
    return *s;
  }
  
}

#endif
