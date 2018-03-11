/**
 * @file   KL_statistic_example.cc
 * @author  <rytis@casamia>
 * @date   Tue Jan 23 15:26:46 2018
 * 
 * @brief  Test of MSE::KL_statistic
 * 
 * Show various overloads of MSE::KL_statistic in action, in particular showing how it works with C++ and GSL formatted distribution functions. Does not take arguments.
 * 
 */
#include "spacings/KL_statistic.hpp"
#include "spacings/distros.hpp"
#include <iostream>
#include <algorithm>
#include <boost/random.hpp>

// #include <stdlib.h>
// #include <math.h>
// #include <vector>
// #include <gsl/gsl_vector.h>
// #include <gsl/gsl_multimin.h>


typedef std::pair<double*,double*> sample_with_pt;
typedef std::vector<double>::iterator iterator;
typedef std::pair<iterator,iterator> sample_with_it;

int main (void) {
  size_t N=10;                  // sample size

  // this is the container for our sample:
  std::vector<double> d(N);

  // you can use almost any type of parameter
  double par = 1.0;
  std::vector<double> v_par{par}; // parameter; using vector only for the purpose of example2.
  double *p_par = v_par.data();   // parameter; using for example1

  sample_with_pt ptr{d.data(),d.data()+N};
  sample_with_it itr{d.begin(),d.end()};

  //  boost::random::exponential_distribution<double> variate(par[0]);
  printf("Compare various methods of invoking \"spacings::KL_statistic\" using C++ and GSL conventions.\n");
  printf("Below we print hopefully the same number, computed using different calls to KL_statistic:\n");

  // let's init the data:
  boost::random::mt19937 rng(time(NULL));
  boost::random::exponential_distribution<double> variate(par);
  for (auto &s:d) s=variate(rng);
  // don't forget to sort data! That's because KL_statistic wants an order statistic of a sample
  std::sort(d.begin(),d.end());

  // An example using generic C++ (flexible)
  std::cout << spacings::KL_statistic<spacings::Exponential,double,sample_with_pt>(spacings::Exponential(),par,ptr) << "\n";

  // The same using iterators:
  std::cout << spacings::KL_statistic<spacings::Exponential,double,sample_with_it>(spacings::Exponential(),par,itr) << "\n";

  // with automatic type deduction:
  std::cout << spacings::KL_statistic(spacings::Exponential(),par,ptr) << "\n";
  std::cout << spacings::KL_statistic(spacings::Exponential(),par,itr) << "\n";

  // How about we create a distro on the fly and use a non-conventional parameter :
  auto expon_with_vec = [] (double x, const std::vector<double>& p) -> double { double arg = p[0]*x; return (arg>0?1.0-exp(-arg):0.0); };
  std::cout << spacings::KL_statistic(expon_with_vec,v_par,ptr) << "\n";
  auto expon_with_ptr = [](double x, double *p) -> double {double arg=(*p)*x; return (arg>0?1.0-exp(-arg):0.0);};
  std::cout << spacings::KL_statistic(expon_with_ptr,p_par,itr) << "\n";


  // Here is a specialization for GSL (note a different signature: two arguments!)
  // The point with GSL is that KL_statistic must look like a GSL function (with signature=(gsl_vector*,void*)->double
  // Therefore, to comply with GSL, we need to do some masquerading.
  // For the first argument, the following two will work (gsl_vector_view will also work):
  gsl_vector_const_view gp = gsl_vector_const_view_array (&par,1);

  // Then we call KL_statistic with a pre-fabricated GSL-complying distribution function (they also have to comply!), the GSL_Exp (you can make your own easily!)
  // Note: you cannot ask the compiler to do type-deduction. Why? Because the first par is not a type (a concrete function), and the second argument cannot be deduced!
  std::cout << spacings::KL_statistic<spacings::GSL_Exp,sample_with_it>(&gp.vector,(void*)&ptr) << "\n";

  // // Can we fabricate a GSL-distribution function on the fly?
  // auto gsl_expon = [] (double x, void*p) -> double { double arg = x* *(double*)p; return (arg>0?1.0-exp(-arg):0.0); };
  // std::cout << spacings::KL_statistic<gsl_expon,sample_with_pt>(&gp.vector,(void*)&ptr) << "\n";
  // // ......
  // // This seems to be out-of-reach for now: that's because the GSL variant of KL_statistic uses a different kind of template mechanism. The first argument must be a name of a concrete function of the specific type (univariate_distribution_function). So, I'll have to find out how does one construct a lambda and mask it as some other type.

  // do not free the vector view!
  //  gsl_vector_free(gp); // NO!
  return 0;
}
