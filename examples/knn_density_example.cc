/**
 * @file   knn_histogram.cc
 * @author  <rytis@casamia>
 * @date   Fri Mar  2 14:51:48 2018
 * 
 * @brief  Demonstrate several modes of use of knn_density
 * 
 * 
 */
#include "clusters/knn.hpp"

#include <iostream>
#include <utility>
#include <boost/iterator/counting_iterator.hpp>

//#include <boost/property_map/property_map.hpp>

using std::abs;
double volume (const double& a, const double& b) {return 2.0*abs(a-b);}

int main (int ac, char*av[])
{
  std::vector<double> data((std::istream_iterator<double>(std::cin)),std::istream_iterator<double>());
  size_t n=data.size();
  int K = int(std::sqrt(double(n)));

  std::vector<double> mass(n);
  std::sort(data.begin(),data.end());
  
  // density using iterators:
  clusters::knn_density(data.begin(),data.end(),mass.begin(),volume,K);
  clusters::knn_density(&data[0],&data[0]+n, &mass[0], volume, K);
  // using boost's property maps:
  std::pair<boost::counting_iterator<int>, boost::counting_iterator<int> > range{0,n};
  // separately constructed
  boost::identity_property_map ident;
  boost::iterator_property_map<std::vector<double>::iterator,boost::identity_property_map> data_map(data.begin(),ident);
  boost::iterator_property_map<std::vector<double>::iterator,boost::identity_property_map> mass_map(mass.begin(),ident);
  clusters::knn_density_map(data_map,mass_map,range,volume,K);
  // or on the fly
  clusters::knn_density_map(boost::make_iterator_property_map(data.begin(),ident),
                            boost::make_iterator_property_map(mass.begin(),ident),
                            range,volume,K);
  
  std::copy(mass.begin(),mass.end(),std::ostream_iterator<double>(std::cout, "\n"));
  return 0;
}
