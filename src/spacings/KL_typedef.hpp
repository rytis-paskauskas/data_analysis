#ifndef KL_typedef_hpp_defined
#define KL_typedef_hpp_defined

namespace spacings {
  /*!
   *  \addtogroup spacings
   *  @{
   */
  /// Parametric test distribution function \f$F(x;\theta)\f$ in GSL format
  typedef double(*univariate_distribution_function)(double, void*);
  // typedef double(*univariate_distribution_function)(double, double*);
  /// Derivative of parametric test df w.r.t. parameter: \f$\frac{d F(x;\theta)}{d\theta}\f$
  typedef void(*univariate_distribution_function_d)(double, void*, double*);
  /// Model for discrete alternative probabilities \f$P_i(\theta)\f$
  typedef double(*distribution_multi)(size_t, double, void*);
  typedef double(*probability_multi)(size_t, void*);
  /*! @} End of Doxygen Groups*/
  
}

#endif
