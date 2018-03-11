/**
 * @file   knn.hpp
 * @author  <rytis@casamia>
 * @date   Sun Mar 11 17:16:07 2018
 * 
 * @brief  
 * 
 * 
 */

#ifndef knn_hpp_defined
#define knn_hpp_defined

#include <algorithm>            // transform,nth_element
#include <iterator>             // iterator_traits
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/lexical_cast.hpp>

//! Hierarchical Density Clustering interfaces and implementations
namespace clusters {
  /*!
   *  \addtogroup clusters
   *  @{
   */

  /** 
   * std algo-like \f$K\f$-nearest neighbor density estimator
   * @param[in] d_beg iterator
   * @param[in] d_end iterator
   * @param[in,out] m_beg iterator
   * @param[in] volume function
   * @param[in] Knn 
   * @todo How about we allow the data to be circular?
   */
  template<class Iterator1, class Iterator2, class VolumeFunction>
  void knn_density(Iterator1 d_beg, Iterator1 d_end, Iterator2 m_beg, VolumeFunction volume, int const Knn = 1) {
    typedef typename std::iterator_traits<Iterator1>::value_type data_type;
    typedef typename std::iterator_traits<Iterator2>::value_type real_type;
    if(d_beg==d_end) return; // trivial
    std::vector<real_type> dij(d_end-d_beg);
    const real_type scale = real_type(Knn)/real_type(d_end-d_beg);
    data_type ref;
    const auto lambda = [&volume,&ref] (data_type& x) -> real_type { return 1.0/volume(x,ref); };
    for (auto cur=d_beg;cur!=d_end;cur++,m_beg++) {
      ref=*cur;
      std::transform(d_beg, d_end, dij.begin(), lambda);
      std::nth_element(dij.begin(),dij.begin() + Knn, dij.end(), std::greater<double>());
      *m_beg = scale*dij[Knn]; // mass = scale/V
    }
  }


  /** 
   * Spanning tree based on the k-th nearest neighbor density estimation
   * 
   * Let's build a MST starting from given data. 
   * Assumptions: the range of data can be traversed using random
   * iterator access. Therefore, it could be also indexed by an
   * integer. The resulting graph should also be indexable by size_t.
   * 
   * We are assuming that the graph has both vertex map convertible to size_t and an edge map convertible to size_t.
   * The default return type can be just copied from the third template argument. 
   * @param[in] d_beg range iterator to data. It could be a pointer, a std::vector<T>::iterator, etc... we are quite flexible!
   * @param[in] d_end range iterator to data. Same as above.
   * @param[in] vm_beg Iterator to vertex weights. The Iterator type must be convertible to floating number (real_type). The underlying data container should be large enough to hold at least as many items as are in the range of data.
   * @param[in] em_beg Iterator to edge weights. Same as above, with the exception of the capacity: we only require it to hold one less items than in the data range.  
   * @param[out] G The resulting graph. Note that G must be a tree, with all the consequences (such as, n. of edges = n. of nodes -1, etc.) G can be initialized and non-pure, we'll clean it before using it. On the outcome, the edges will be randomly ordered. 
   * @param[in] volume Tells what is the volume of a ball with a center at the reference point and a boundary touching the second point. 
   * @param[in] Knn the "k" in the k-th nearest neighbor.
   */
  template
  <class Iterator1, class Iterator2,
   class graph_type = typename boost::adjacency_list<
     boost::listS,
     boost::vecS,
     boost::undirectedS,
     boost::no_property,
     boost::property<boost::edge_index_t, std::size_t>, // edge index can be used to relate Edges(G) <-> edge mass vector
     boost::no_property,
     boost::listS
     >,
   class VolumeFunction,
   class data_type = typename std::iterator_traits<Iterator1>::value_type,
   class real_type = typename std::iterator_traits<Iterator2>::value_type
   >
  void
  knn_density_tree (Iterator1 d_beg,
                    Iterator1 d_end,
                    Iterator2 vm_beg,
                    Iterator2 em_beg,
                    graph_type& G,
                    VolumeFunction volume,
                    const int Knn)
  {
    typedef typename boost::graph_traits<graph_type>::vertex_descriptor vertex_descriptor;
    typedef typename boost::property_map<graph_type, boost::edge_index_t>::type EdgeIndexMap;
    typedef typename boost::property_map<graph_type, boost::vertex_index_t>::type VertexIndexMap;
    typedef typename boost::iterator_property_map<typename std::vector<vertex_descriptor>::iterator,VertexIndexMap> IndexVertexMap;
    typedef typename boost::iterator_property_map<Iterator2,VertexIndexMap> RealVertexMap;
    typedef typename boost::iterator_property_map<Iterator2, EdgeIndexMap> RealEdgeMap;
    typedef typename boost::iterator_property_map<std::vector<bool>::iterator, EdgeIndexMap> BoolEdgeMap;

    struct mst_marker : boost::default_dijkstra_visitor {
      mst_marker(IndexVertexMap const& p, BoolEdgeMap& m, vertex_descriptor const u) : m_pred(p), m_mark(m), start(u) {}
      /*
       * Mark maximum spanning tree edges.
       * note According to boost description of dijkstra visitor, at examine_vertex event, the
       * dijkstra algorithm knows that (p[u],u) a MST edge.
       */
      void examine_vertex(const vertex_descriptor& v, const graph_type&g) const {
        if(v==start) return;    // do nothing;
        vertex_descriptor u = m_pred[v];
        auto e = boost::edge(u,v,g); // assuming there is exactly one edge from u to v.
        // most typical problem here are orphaned vertices,
        // i.e. vertices that have no edges, with an error message
        // like this: "error mst_marker: 648 - 648 : edge does not
        // exist."
        if(!e.second) {
          throw std::runtime_error(std::string("mst_marker: ")+boost::lexical_cast<std::string>(u) + " - " + boost::lexical_cast<std::string>(v)+std::string(" : edge does not exist. Attempting to mark a self-edge may be a symptom of an orphaned vertex, which in turn may indicate a defective metric of event similarity comparison."));
        }
        boost::put(m_mark,e.first,true);
        return;
      }
    private:
      const IndexVertexMap& m_pred;
      BoolEdgeMap& m_mark;
      vertex_descriptor start;
    };

    // 1st step: get vertex mass:
    knn_density(d_beg,d_end,vm_beg,volume,Knn);

    // 2nd step: get edges and the tree:
    G.clear();
    G = graph_type(d_end-d_beg);

    // let's inflate G by adding all-to-all edges:
    auto vv = boost::vertices(G);
    for(auto idx=vv.first;idx!=vv.second;idx++)
      for(auto idy=vv.first;idy<idx;++idy) {
        auto ee=boost::add_edge(*idy,*idx,G);
        if(!ee.second) throw std::runtime_error(std::string("... error. Adding edge failed :") + boost::lexical_cast<std::string>(*idy) + std::string(" - ") + boost::lexical_cast<std::string>(*idx));
      }
    
    std::vector<bool> mst_marker_data(boost::num_edges(G),false);
    std::vector<vertex_descriptor> predecessor_data(boost::num_vertices(G));
    std::vector<real_type> distance_data(boost::num_vertices(G));
    std::vector<real_type> em_temporary(boost::num_edges(G));
    Iterator2 em = em_beg;

    EdgeIndexMap edge_index_map = boost::get(boost::edge_index,G);
    VertexIndexMap vertex_index_map = boost::get(boost::vertex_index,G);

    RealEdgeMap edge_weight_map(em_temporary.begin(), edge_index_map);
    BoolEdgeMap mst_marker_map(mst_marker_data.begin(), edge_index_map);
    IndexVertexMap predecessor_map(predecessor_data.begin(), vertex_index_map);
    RealVertexMap distance_map(distance_data.begin(), vertex_index_map);
    RealVertexMap mass_map(vm_beg, vertex_index_map);
    auto data_map = make_iterator_property_map(d_beg, vertex_index_map);

    // fill up edges:
    {
      size_t eid=0;
      for(auto ee=boost::edges(G);ee.first!=ee.second;ee.first++,eid++) boost::put(edge_index_map,*ee.first,eid);
    }

    real_type scale = double(Knn)/double(d_end-d_beg);

    for(auto ee=boost::edges(G);ee.first!=ee.second;ee.first++) {
      auto s=boost::source(*ee.first,G);
      auto t=boost::target(*ee.first,G);
      // In this way, only the factor would be hard-defined.
      real_type v = std::min({mass_map[s],mass_map[t],scale/volume(data_map[s],data_map[t])});
      boost::put(edge_weight_map,*ee.first,v);
    }
    
    auto start_vertex = *(boost::vertices)(G).first;
    boost::dijkstra_shortest_paths
      (G
       , start_vertex
       , predecessor_map
       , distance_map
       , edge_weight_map
       , vertex_index_map
       , // compare:
       std::less<double>()
       , // combine:
       // assuming that similarity is monotone decreasing, whereas
       // dijkstra priority is monotone increasing, we can use 1/w as
       // a map between the two. Note that 1/w is not related to the
       // originally used 1/distance similarity measure.
       // @todo Should I worry about w=0?
       [](const real_type d, const real_type sim) -> real_type { return 1.0/sim; }
       , // infinity element:
       (std::numeric_limits<double>::max)()
       , // zero element:
       0
       , mst_marker(predecessor_map,mst_marker_map,start_vertex)
       );
    // erase edges that aren't marked by mst_edge_map
    {
      auto ee = boost::edges(G);
      for (auto next = ee.first; ee.first != ee.second; ee.first = next) {
        next++;
        if(mst_marker_map[*ee.first]==false) boost::remove_edge(*ee.first, G);
        else *em++ = edge_weight_map[*ee.first];
      }
    }
    // relabel edges:
    {
      size_t eid=0;
      for(auto ee=boost::edges(G);ee.first!=ee.second;ee.first++,eid++) boost::put(edge_index_map,*ee.first,eid);
    }
    assert((em==em_beg+boost::num_edges(G)) && (boost::num_edges(G)+1==boost::num_vertices(G)));
  } // done!
  
  /** 
   * knn_density variant which can be called with boost::property_map and similar types.
   * 
   * @param[in] d Any property map that maps to data
   * @param[in] m Any property map that can be converted to a real type. 
   * @param[in] range The range of iterators. The iterators can be
   * arbitrary keys, but they must be convertible to at least forward
   * iterators. A typical simplest example would be to use
   * boost::counting_iterator.
   * @param[in] volume Volume of a ball with the reference point in center and the referred point on the boundary.
   * @param[in] Knn the "k" in k-nearest neighbors.
   */
  template<class DataMap, class MassMap, class KeyIterator, class VolumeFunction>
  void knn_density_map(DataMap d, MassMap m, typename std::pair<KeyIterator,KeyIterator> range, VolumeFunction volume, int Knn)
  {
    typedef typename boost::property_traits<MassMap>::value_type real_type;
    if(range.second==range.first) return; // trivial
    std::vector<real_type> dij(range.second-range.first);
    real_type ref;
    const auto lambda = [&d,&volume,&ref] (size_t i) -> real_type { return 1.0/volume(d[i],ref); };
    const real_type scale = real_type(Knn)/real_type(range.second-range.first);
    for (auto cur=range.first;cur!=range.second;cur++) {
      ref=d[*cur];
      std::transform(range.first, range.second, dij.begin(), lambda);
      std::nth_element(dij.begin(), dij.begin() + Knn, dij.end(), std::greater<real_type>());
      m[*cur] = scale*dij[Knn];
    }
  }

}      // namespace clusters
#endif // knn_hpp_defined
