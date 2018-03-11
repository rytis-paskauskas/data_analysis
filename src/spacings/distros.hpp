/**
 * @file   distros.hpp
 * @author  <rytis@casamia>
 * @date   Sun Mar 11 17:40:19 2018
 * 
 * @brief  
 * 
 */
#include <math.h>
#include <iterator>             // for iterator_traits
#include <type_traits>          // static_assert
#include <limits>               // numeric_limits
//#include <stdexcept>
#include <assert.h>
#include <gsl/gsl_vector.h>

#ifndef spacings_distros_hpp_defined
#define spacings_distros_hpp_defined

#include "KL_typedef.hpp"

//! Maximum Spacing Estimation interfaces and implementations
namespace spacings {
  /*!
   *  \addtogroup spacings
   *  @{
   */

  /**
   * The exponential distribution \f$\mathrm{Exp}(\lambda)\f$
   *
   * Some distributions will be used more frequently: we provide some elementary helpers. 
   * This one is the exponential distribution:
   * \f{equation*}{X\sim\mathrm{Exp}(\lambda)~:\quad X\ge 0,\quad \mathbb{P}(X>x)=\mathrm{e}^{-\lambda x}\f}
   *
   * @note that when this is used with KL_statistic, there is an
   * implied Parameter=double specification.
   * 
   */
  struct Exponential {
    typedef double value_type;
    /** 
     * The operator used by the spacings functions
     * @return Value of the distribution \f$F(x,p)\f$.
     */
    value_type operator()(value_type x, const value_type& p) const {
      value_type arg = p*x;
      return (arg>0?1.0-exp(-arg):0.0);
    }
  };

  /**
   * Helper distribution for the purely exponential multiple alternative model.
   * This structure can be used as a prototype for general multi-alternative distribution construction for each alternative.
   */
  struct Markov_continuous {
    typedef double value_type;
    Markov_continuous(size_t n) : num(n) {} /**< Constructor with the number of alternatives */
    /** 
     * The operator used by the spacings' functions. 
     * @note The operator signature implies Parameter = double*, when it is used with mKL_statistic
     * @return The distribution function of the continuous variable, conditional to the occurrence of the \f$a\f$'th alternative. In the min-model, it is given by \f$\displaystyle \mathbb{P}(X\le x\vert A=a)=\exp{\left(-\left\{\sum_{i=1}^n p_i\right\}x\right)}\f$
     */
    value_type operator()(size_t a, value_type x, const value_type* par) const {
      value_type arg=value_type(0);
      const value_type* tmp=par;
      for(size_t j=0;j!=num;j++) arg += *(tmp++);
      arg=x*arg;
      return (arg>0?1.0-exp(-arg):0.0);
    }
  private:
    size_t num;
  };
  
  /**
   * Helper function for a purely exponential multiple alternative model.
   * 
   * This is the discrete probability model for the minimum of the multiple exponential random variables.
   */
  struct Markov_discrete {
    typedef double value_type;
    Markov_discrete(size_t n) : num(n) {} /**< Constructor with the number of alternatives */
    /** 
     * The operator used by the KL_statistic
     * @note The operator signature implies Parameter = double*, when it is used with mKL_statistic
     * @return The probability of occurrence of the \f$a\f$'th alternative which, in the min-model, is  \f$\displaystyle \mathbb{P}(a)=\frac{n_a}{\sum_{i=1}^n n_i}\f$
     */
    value_type operator()(size_t a, const value_type* par) const {
      value_type sum=value_type(0);
      const double* tmp=par;
      assert(a<num);
      for(size_t j=0;j!=num;j++) sum += *(tmp++);
      return *(par+a)/sum;
    }
  private:
    size_t num;
  };

  /// GSL SPECIFICATIONS FOR FREQUENTLY USED DISTROS

  /** 
   * This one satisfies the signature requirements of GSL. 
   * See KL_statistic_example.cc and KL_MSE_example.cc. 
   * 
   * @param x 
   * @param p 
   * 
   * @return 
   */  double GSL_Exp(double x, void *p) {
    double arg = x * *(double*)p;
    return (arg>0?1.0-exp(-arg):0.0);
  }

  /** 
   * Multi-alternative distribution for the multi-exponential model.
   * @note Make sure that p={N,par1,...parN}.
   * @param a 
   * @param x 
   * @param p 
   * 
   * @return 
   */
  double GSL_Markov_continuous(size_t a, double x, void *p) {
    size_t n=*(size_t*)p;
    //    p+=sizeof(size_t);
    p = (size_t*)p + 1;
    double*last=(double*)p+n;
    assert((double*)p + a < last);
    double arg=0;
    while((last--)!=(double*)p) arg+=*last;
    arg=x*arg;
    // assert(a<*(size_t*)p);
    // p+=sizeof(size_t);
    // double arg=*((double*)p+a)*x;
    return (arg>0?1.0-exp(-arg):0.0);
  }

  /** 
   * Multi-alternative discrete probability model following the signature requirements of the GSL.
   * @note Make sure that p={N,par1,...parN}.
   * @param a 
   * @param p 
   * 
   * @return 
   */
  double GSL_Markov_discrete(size_t a, void *p) {
    size_t n=*(size_t*)p;
    //    p+=sizeof(size_t);
    p = (size_t*)p+1;
    double*last=(double*)p+n;
    assert((double*)p + a < last);
    double sum=0;
    while((last--)!=(double*)p) sum+=*last;
    return *((double*)p+a)/sum;
  }
  /*! @} End of Doxygen Groups*/
}

#endif

