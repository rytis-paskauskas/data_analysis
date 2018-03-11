/**
 * @file   MSE.hpp
 * @author  <rytis@casamia>
 * @date   Tue Jan 23 15:45:44 2018
 * 
 * @brief  
 * 
 * 
 */
#include <math.h>
#include <iterator>             // for iterator_traits
#include <type_traits>          // static_assert
#include <limits>               // numeric_limits
//#include <stdexcept>
#include <assert.h>
#include <gsl/gsl_vector.h>

#ifndef KL_statistic_hpp_defined
#define KL_statistic_hpp_defined

#include "KL_typedef.hpp"


#ifndef M_EULER
#define M_EULER      0.57721566490153286060
#endif

//! Maximum Spacing Estimation interfaces and implementations
namespace spacings {
  /*!
   *  \addtogroup MSE
   *  @{
   */

  /// \private
  template<class Iterator> typename std::iterator_traits<Iterator>::difference_type look_ahead(std::pair<Iterator,Iterator>&);

  /** 
   * The Kullback-Leibler entropy statistic
   * \f{equation*}{ S_n(\theta)=-\frac{1}{(n+1)}\sum_{i=0}^n\ln{\big[(n+1)( F_\theta(x_{i+1})-F_\theta(x_i))\big]}-\gamma \f}
   * - \f$x_1 \le x_2 \dots \le x_n\f$ is an ordered sample of size \f$n\f$, possibly with ties (i.e. equalities are allowed). In the above formula, \f$x_0=-\infty\f$, \f$x_{n+1}=\infty\f$.
   * - \f$F_\theta(x)=F(x,\theta)\f$ is a distribution function depending on parameter \f$\theta\f$.
   * - \f$\gamma\f$ Euler's constant
   *
   * This statistic is related to the Kullback-Leibler entropy \f$\mathcal{D}(f\vert\vert g)\f$, to the differential entropy \f$H[f]\f$, to the Moran's statistic of uniformity \f$M_n\f$, and to the spacings parametric modelling and goodness-of-fit tests. We use it in the Maximum Spacings Estimation and goodness-of-fit.

   * Relation to the Kullback-Leibler entropy: 
   \f{equation*}{ \lim_{n\to\infty} S_n(\theta) - \gamma = \mathcal{D}( F_X \vert\vert F_{\theta} ) \f}
   * 
   * Relation to Moran's statistic for uniformity:
   * If \f$X_i\in [0,1]\f$, and \f$F_\theta(x) = x\f$ on \f$[0,1\f$ then \f$M_n=(n+1)S_n(\theta)\f$ is the Moran's statistic of uniformity. 
   *
   * Maximum Spacings Estimation:
   * Our main interest is to use \f$S_n(\theta)\f$ as the Maximum Spacings estimation (MSE) and the goodness-of-fit testing of the parametric hypothesis test:
   \f{equation*}{ \text{H: each } X_i \text{ is i.i.d. rv} \sim X_\theta : \mathbb{P}(X_\theta\le x)=F_\theta(x)\f}
   I.e. the data set comes from a family of distribution functions \f$F_\theta(x)\f$.

   *
   * Ties in data are taken care of as follows.
   *


   * The \f$n\f$ scaling has been tested and shows a finite \f$\lim_{n\to\infty}S_n(\theta)\f$ limit.
   *
   * We define the Kullback-Leibler entropy statistic as 

   * The normalization convention is the same as in Bo Ranneby. But we
   * use the negative sign: Then it becomes KL entropy (non-negative)
   * and we can use minimization procedures of GSL.

   
   * Relation to Kullback-Leibler entropy

   * This normalization convention is, I think, the original Ronneby's normalization. And also that of Moran's statistic. 

   * @param[in] Iterator type of iterator. Typically would be \c double* or \c vector<double>::iterator
   * @param[in] theta range for the the parameter \f$\theta\f$ 
   * @param[in] F template for the distribution function. Distribution must have the following member function:
   * \code{.hpp}value_type Distribution::operator()(value_type,iterator_type) const\endcode to return values from [0,1] (the distribution function).
   * A simplest implementation of the exponential distribution \f$\mathsf{Exp}(\lambda)\f$: 
   * \snippet MSE.hpp Exp
   * @param[in] samp data sample of class Sample. Sample must have a member function \code{.hpp} std::pair<Iterator,Iterator> Sample::operator()() const \endcode 
   * An example of a simple structure using pointers for the range:
   * \snippet KL_statistic_test.cc sample with pointer
   * 
   * @return Value of \f$S_n(\theta)\f$
   */
  template<class Distribution,class Parameter,class SampleRange,typename value_type=typename std::iterator_traits<typename SampleRange::first_type>::value_type >
  value_type
  KL_statistic(Distribution F, Parameter theta, SampleRange range)
  {
    static_assert(std::is_same<typename SampleRange::first_type, typename SampleRange::second_type>::value, "KL_statistic: SampleRange must be a pair of iterator to the range of a sample");
    using std::log;
    const value_type m=value_type(range.second-range.first+1);
    value_type result(0),next(0),prev(0),tie(0);
    for(;range.first!=range.second;++range.first) {
      next = F(*range.first,theta); // THe distribution F should know how to operate with theta!
      if((next<0)||(next>1.0)) throw std::range_error (std::string(__FUNCTION__) + std::string(": distribution function out of range"));
      if(next==NAN) throw std::range_error (std::string(__FUNCTION__) + std::string(": distribution function NaN"));      
      tie=value_type(look_ahead(range));
      result+=tie*log((next-prev)/tie);
      if(result==+std::numeric_limits<value_type>::infinity() || result==-std::numeric_limits<value_type>::infinity())
        throw std::range_error( std::string(__FUNCTION__) + std::string(" : infinite function value (most likely, bad data or bad parameters)"));
      prev = next;
    }
    next = value_type(1);
    result += log(next-prev);
    return -log(m)-(result/m) - M_EULER; // we follow Beirlant to add M_EULER.
  }
  // this scaling has been tested and shows a finite N->infty limit (close to zero).
  // it correspondes to the target function
  // \f$ \frac{1}{N}\sum_i^N \ln{\big( n( F(x_i)-F(x_{i-1}))\big)} \f$
  // of which we take the negative for minimization.
  // This normalization convention is, I think, the original Ronneby's normalization. 


  /** 
   * A wrapper to KL_statistic made specifically fo the Gnu Scientific
   * Library functions like gsl_fminimizer, gsl_fdfminimizer, etc. The
   * GSL requires a specific signature for its functions (in order to
   * be usable with GSL minimization etc. procedures). This
   * restriction forces us to define a standard signature for the
   * allowed distribution function: univariate_distribution_function,
   * and to use this kind of template, so that the distribution is
   * specified at compile time. In particular this allows to take a
   * pointer in the form \code{.hpp}&KL_statistic<F>\endcode

   * @param[in] p the vector containing distribution parameter(s)
   * @param[in] s the pointer to the data sample on which the spacings
   * calculations are made. It must be a pointer to
   * \code{.hpp}struct GSL_sample\endcode
   *
   * @note after a lot of trials, this is a really short and nifty
   * solution. See at the end for older implementations.
   * 
   * @return 
   */
  template <univariate_distribution_function F,class SampleRange>
  double KL_statistic(const gsl_vector*p, void*s)
  {
    assert(p->stride==1);       // we don't want to handle strides.
    SampleRange *rp = static_cast<SampleRange*>(s);
    return KL_statistic<univariate_distribution_function,void*,SampleRange,double>(F, (void*)(p->data),*rp);
  }
  
  /**
   * The exponential distribution \f$\mathrm{Exp}(\lambda)\f$
   *
   * Some distributions will be used more frequently: we provide the
   * elementary exponential distribution:
   * \f{equation*}{X\sim\mathrm{Exp}(\lambda)~:\quad X\ge 0,\quad \mathbb{P}(X>x)=\mathrm{e}^{-\lambda x}\f}
   *
   * @note that when this is used with KL_statistic, there is an
   * implied Parameter=double specification.
   * 
   */
  /// [Exp]
  struct Exponential {
    typedef double value_type;
    /** 
     * The operator used by the spacings functions
     * 
     * 
     * @return Value of the distribution \f$F(x,p)\f$.
     */
    value_type operator()(value_type x, const value_type& p) const {
      value_type arg = p*x;
      return (arg>0?1.0-exp(-arg):0.0);
    }
  };
  /// [Exp]

  /**
   * Helper distribution for the purely exponential multiple alternative model.
   * 
   * This structure shows how a continuous distribution for each alternative is constructed.
   */
  /// [multi Exp]
  struct Markov_continuous {
    typedef double value_type;
    Markov_continuous(size_t n) : num(n) {} /**< Constructor with the number of alternatives */
    /** 
     * 
     * @note The operator signature implies Parameter = double*, when it is used with mKL_statistic
     * 
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
   * probability model for the minimum of the multiple exponential random variables
   * the signature of operator() implies that Parameter=double*
   */
  /// [minExp model]
  struct Markov_discrete {
    typedef double value_type;
    Markov_discrete(size_t n) : num(n) {} /**< Constructor with the number of alternatives */
    /** 
     * 
     * @note The operator signature implies Parameter = double*, when it is used with mKL_statistic
     * 
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
  
  double GSL_Exp(double x, void *p) {
    double arg = x * *(double*)p;
    return (arg>0?1.0-exp(-arg):0.0);
  }

  // Make sure that p={N,par1,...parN}.
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

  // Make sure that p={N,par1,...parN}.
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




  
  /** 
   * Counts the number of ties and moves the left boundary of the range to the first element past the tie
   * 
   * @param[in] Iterator 
   * @param[in] range 
   * 
   * @return number of ties
   */
  template<class Iterator=double*>
  typename std::iterator_traits<Iterator>::difference_type
  look_ahead(std::pair<Iterator,Iterator>& range) {
    typename std::iterator_traits<Iterator>::difference_type count(1);
    typename std::iterator_traits<Iterator>::value_type t = *range.first;
    while((range.first+1)!=range.second && *(range.first+1)==t) {
      ++range.first;
      ++count;
    }
    return count;
  }

  int look_ahead(double const *& x, double const * y) {
    int m=1;
    double t = *x;
    while((x+1)!=y && *(x+1) == t) {
      x++;
      m++;
    }
    return m;
  }




  /*! @} End of Doxygen Groups*/
}

#endif





  
  //  // old implementation: 
  // template <class Distribution,class Sample, class Iterator=double*>
  // typename std::iterator_traits<Iterator>::value_type
  // KL_statistic(const std::pair<Iterator,Iterator>& theta, Distribution F, const Sample* samp)
  // {
  //   using std::log;
  //   typedef typename std::iterator_traits<Iterator>::value_type value_type;
  //   auto input = (*samp)();
  //   const value_type m=value_type(input.second-input.first+1);
  //   value_type result(0),next(0),prev(0),tie(0);
  //   for(;input.first!=input.second;++input.first) {
  //     next = F(*input.first,theta.first); // This one should know how to operate with theta!
  //     if((next<0)||(next>1.0)) throw std::range_error (std::string(__FUNCTION__) + std::string(": distribution function out of range"));
  //     if(next==NAN) throw std::range_error (std::string(__FUNCTION__) + std::string(": distribution function NaN"));      
  //     tie=value_type(look_ahead(input));
  //     result+=tie*log((next-prev)/tie);
  //     if(result==+std::numeric_limits<value_type>::infinity() || result==-std::numeric_limits<value_type>::infinity())
  //       throw std::range_error( std::string(__FUNCTION__) + std::string(" : infinite function value (most likely, bad data or bad parameters)"));
  //     prev = next;
  //   }
  //   next = value_type(1);
  //   result += log(next-prev);
  //   // this scaling has been tested and shows a finite N->infty limit (close to zero).
  //   // it correspondes to the target function
  //   // \f$ \frac{1}{N}\sum_i^N \ln{\big( n( F(x_i)-F(x_{i-1}))\big)} \f$
  //   // of which we take the negative for minimization.
  //   // This normalization convention is, I think, the original Ronneby's normalization. 
  //   return -log(m)-(result/m) - M_EULER; // we follow Beirlant to add M_EULER.
  // }



  // // Older implementation:
  // template <class Continuous, class Discrete, class Sample, class Iterator=double*>
  // typename std::iterator_traits<Iterator>::value_type
  // mKL_statistic(const std::pair<Iterator,Iterator>& theta, Continuous F, Discrete P, const Sample* samp)
  // {
  //   using std::log;
  //   typedef typename std::iterator_traits<Iterator>::value_type value_type;
  //   typedef typename std::iterator_traits<Iterator>::difference_type integral_type;
  //   value_type ret = 0, tmp;
  //   value_type sum = value_type(0);
  //   size_t n=(*samp)();
  //   std::vector<integral_type> freq(n);
  //   typename std::vector<integral_type>::iterator f=freq.begin();
  //   for(size_t a=0;a!=n;a++,f++) {
  //     auto r = (*samp)(a);
  //     *f=r.second - r.first;
  //     sum += value_type(*f);
  //   }
  //   ret = 0.0;
  //   for(size_t a=0;a!=n;a++) {
  //     auto Fa = [&F,a] (value_type x,Iterator ptr) -> value_type { return F(a,x,ptr); };
  //     auto Sa = [samp,a](void) -> std::pair<Iterator,Iterator> { return (*samp)(a); };
  //     if(freq[a]>0) {
  //       tmp = value_type(freq[a])/sum;
  //       ret += tmp *(log(tmp/P(a,theta.first)) + 0); // KL_statistic<decltype(Fa),decltype(Sa),Iterator>(theta,Fa,&Sa));
  //     }
  //   }
  //   return ret;
  // }

  // // This is an older implementation:
  // template <univariate_distribution_function F>
  // double KL_statistic(const gsl_vector*p, void*s)
  // {
  //   typedef double* ptr_type;
  //   //    typedef double  value_type;
  //   assert(p->stride==1);
  //   std::pair<ptr_type,ptr_type> theta = {p->data,p->data + p->size};
  //   //   // Lambda wrapper to F (unnecessary):
  //   //   auto GSL_distribution = [] (value_type x, ptr_type ptr) -> value_type { return F(x, (void*)ptr); };
  //   return KL_statistic<univariate_distribution_function,GSL_sample,ptr_type>(theta, F,(const GSL_sample*)s);
  //   //   return KL_statistic<decltype(GSL_distribution),GSL_sample,ptr_type>(theta, GSL_distribution,(const GSL_sample*)s);
  // }


  // // An older implementation, high on lambda sugar.
  // template <distribution_multi C, probability_multi D>
  // double mKL_statistic(const gsl_vector*p, void*s)
  // {
  //   typedef double* ptr_type;
  //   typedef double  value_type;
  //   assert(p->stride==1);
  //   // Lambda wrapper to F:
  //   auto GSL_cont = [p] (size_t a, value_type x, ptr_type ptr) -> value_type { return C(a, x, (void*)ptr, size_t(p->size)); };
  //   auto GSL_disc = [p] (size_t a, ptr_type ptr) -> value_type { return D(a, (void*)ptr, size_t(p->size)); };
  //   std::pair<ptr_type,ptr_type> theta = { p->data, p->data+p->size };
  //   return mKL_statistic<decltype(GSL_cont),decltype(GSL_disc),GSL_msample,ptr_type>(theta, GSL_cont,GSL_disc,(const GSL_msample*)s);
  // }

  // template <class Iterator>
  // struct Exponential {
  //   typedef typename std::iterator_traits<Iterator>::value_type value_type;
  //   /// [distribution signature]
  //   value_type operator()(value_type x, Iterator lambda_ptr) const
  //   /// [distribution signature]
  //   {
  //     value_type arg = *lambda_ptr*x;
  //     return (arg>0?1.0-exp(-arg):0.0);
  //   }
  // };
  /// [multi Exp]
  // template <class Iterator>
  // struct mExponential {
  //   typedef typename std::iterator_traits<Iterator>::value_type value_type;
  //   value_type operator()(size_t a,value_type x, Iterator lambda_ptr) const {
  //     value_type arg = *(lambda_ptr+a)*x;
  //     return (arg>0?1.0-exp(-arg):0.0);
  //   }
  // };

  // template<class Iterator>
  // struct minexp_model {
  //   typedef typename std::iterator_traits<Iterator>::value_type value_type;
  //   minexp_model(size_t n) : num(n) {}
  //   value_type operator()(size_t a, Iterator lambda_ptr) const {
  //     assert(a<num);
  //     value_type sum=value_type(0);
  //     Iterator tmp=lambda_ptr;
  //     for(size_t j=0;j!=num;j++) sum += *(tmp++);
  //     return *(lambda_ptr+a)/sum;
  //   }
  // private:
  //   size_t num;
  // };
  /// [minExp model]

  // /** 
  //  * The Maximum Spacings Estimate statistic
  //  * 
  //  * This function implements the Maximum Spacings Estimator (MSE) function, related to the hypothesis test H: data \f${\bf X}_N=(X_1,\dots, X_N)\f$ has a distribution function \f$F(x;\theta)\f$ where \f$\theta\f$ is given (i.e. \f$\mathbb{P}(X\le x)=F(x;\theta)\f$). The value of the MSE function is defined as 
  //  * \f{equation}{ S_n(\theta) = -\frac{1}{(n+1)}\sum_{i=0}^n \ln{\big( (n+1)( F(x_{i+1})-F(x_i))\big)} \f}
  //   this scaling has been tested and shows a finite N->infty limit (close to zero).
  //   it correspondes to the target function


  //   version with ties.


  //   of which we take the negative for minimization.
  //   This normalization convention is, I think, the original Ronneby's normalization. And also that of Moran's statistic. 

  //  * @param[in] v is the vector of distribution parameters
  //  * @param[in] p is the void pointer to the array of random samples (which are parameters to the MSE function)
  //  * 
  //  * @return value of the MSE estimator
  //  */
  // template <univariate_distribution_function F>
  // double
  // KLentropy_tie(const gsl_vector* v, void*p)
  // {
  //   using std::log;
  //   info* i = (info*)p;
  //   double const * x_end=i->data + i->size;
  //   double const * x = i->data; // pointers to the first and past-last datums.
  //   double result=0,next,prev=0.0; // F(*x,v->data);
  //   for(;x!=x_end;x++) {
  //     next = F(*x,v->data);
  //     if(next<0 || next > 1.0) throw std::range_error (std::string(__FUNCTION__) + std::string(": distribution function out of range"));
  //     if(next==NAN) throw std::range_error (std::string(__FUNCTION__) + std::string(": distribution function NaN"));      
  //     // look-ahead step
  //     //      double const * y = x;
  //     int m=look_ahead(x,x_end);
  //     result += double(m)*log((next-prev)/double(m));
  //     if(result==+std::numeric_limits<double>::infinity() || result==-std::numeric_limits<double>::infinity())
  //       throw std::range_error( std::string(__FUNCTION__) + std::string(" : infinite function value (most likely, bad data or bad parameters)"));
  //     prev = next;
  //   }
  //   next = 1.0;
  //   result += log(next-prev); 
  //   // this scaling has been tested and shows a finite N->infty limit (close to zero).
  //   // it correspondes to the target function
  //   // \f$ \frac{1}{N}\sum_i^N \ln{\big( n( F(x_i)-F(x_{i-1}))\big)} \f$
  //   // of which we take the negative for minimization.
  //   // This normalization convention is, I think, the original Ronneby's normalization. 
  //   return -log(double(1 + i->size)) - (result/double(1 + i->size));
  // }


  // /** 
  //  * Older version of @ref KLentropy_tie, without testing for ties.
  //  * 
  //  * @param v 
  //  * @param p 
  //  * 
  //  * @return 
  //  */
  // template <univariate_distribution_function F>
  // double
  // KLentropy(const gsl_vector* v, void*p)
  // {
  //   using std::log;
  //   info* i = (info*)p;
  //   double const * x = i->data, *x_end=i->data + i->size; // pointers to the first and past-last datums.
  //   double result=0,next,prev=0.0; // F(*x,v->data);
  //   for(;x!=x_end;x++) {
  //     next = F(*x,v->data);
  //     if(next==std::numeric_limits<double>::infinity()) {
  //       std::cerr << "infinity encountered for x=" << *x << '\n';
  //     }
  //     result += log(next-prev);
  //     if(result==std::numeric_limits<double>::infinity()) {
  //       std::cerr << "infinity encountered F[i]=" << next << " F[i-1]=" << prev << '\n';
  //     }
  //     prev = next;
  //   }
  //   next = 1.0;
  //   result += log(next-prev); 
  //   // this scaling has been tested and shows a finite N->infty limit (close to zero).
  //   // it correspondes to the target function
  //   // \f$ -\frac{1}{(n+1)}\sum_{i=0}^n \ln{\big( (n+1)( F(x_{i+1})-F(x_i))\big)} \f$
  //   // of which we take the negative for minimization.
  //   // This normalization convention is, I think, the original Ronneby's normalization. And also that of Moran's statistic. 
  //   return -log(double(1 + i->size)) - (result/double(1 + i->size));
  // }

