/**
 * @file   MSE.hpp
 * @author  <rytis@casamia>
 * @date   Tue Jan 23 15:45:44 2018
 * 
 * @brief  
 * 
 * 
 */
//#include <math.h>
//#include <utility>              // for std::pair
#include <iterator>             // for iterator_traits, distance
#include <cmath>
#include <vector>
#include <limits>               // numeric_limits
#include <stdexcept>
#include <assert.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>

#ifndef KL_MSE_hpp_defined
#define KL_MSE_hpp_defined

#include "KL_typedef.hpp"
#include "KL_statistic.hpp"

//! Maximum Spacing Estimation interfaces and implementations
namespace spacings {
  /*!
   *  \addtogroup spacings
   *  @{
   */


  /** 
   * Minimizer of the KL statistic.
   * 
   * Find the minimum of Kullback-Leibler statistic using GSL's minimum search with function value only.

   * @param[out] statistic The minimum KL statistic \f$\displaystyle \hat{S}_n=S_n(\hat{\theta})=\min_{\theta}S_n(\theta)\f$
   * @param[out] theta The estimated minimizer \f$\displaystyle \hat{\theta} = \arg\min_{\theta}S_n(\theta)\f$
   * @param[in] step The initial step size used by GSL search routine. dim(step)=dim(theta).
   * @param[in] dim Dimension of the parameter \f$\theta\f$
   * @param[in] samp Pointer to the first datum of the sample.
   * @param[in] npt Sample size
   * @param[in] T The GSL search algorithm to be used. The choice is from: 
   * @param[out] iter Number of iterations performed
   * @param[in] max_iter Maximum number of iterations allowed
   * @param[out] error Granularity of the \f$\hat{\theta}\f$ solution.
   * @param[in] tol The largest granularity allowed.
   *
   * @note It is assumed that theta, samp, step are contiguous double arrays.

   * @note Don't forget to switch off the standard gsl error handler
   * before calling this routine if you don't wish the program to
   * abort upon failure. This can be done by remembering the old
   * handler as follows: \code{.cc} gsl_error_handler_t * oh = gsl_set_error_handler_off();\endcode
   * If needed, return control to the old handler after the routine, do
   * \code{.cc} gsl_set_error_handler (oh);\endcode
   *
   * @return the GSL status variable. If
   */
  template<
    univariate_distribution_function Distribution
    >
  int 
  KL_MSE_fmin(double& statistic,
              double *theta, const double *step, const size_t dim,
              const double *samp, const size_t npt,
              gsl_multimin_fminimizer_type const *T,
              size_t& iter, const size_t max_iter,
              double& error, const double tol
              )
  {
    typedef std::pair<const double*,const double*> sample;
    int gsl_status = GSL_SUCCESS;
    sample p={samp,samp+npt};
    gsl_vector_const_view gtheta = gsl_vector_const_view_array (theta, dim);
    gsl_vector_const_view gstep =  gsl_vector_const_view_array (step, dim);
    gsl_multimin_function mf;
    gsl_multimin_fminimizer *state = NULL;
    try {
      mf.n=dim;
      mf.f=&KL_statistic<Distribution,sample>;
      mf.params=(void*)&p;
      state = gsl_multimin_fminimizer_alloc (T, dim);
      gsl_multimin_fminimizer_set (state, &mf, &gtheta.vector, &gstep.vector);
      // out.name = gsl_multimin_fminimizer_name(state);
      // out.iter = 0;
      do {
        gsl_status = gsl_multimin_fminimizer_iterate(state);
        if (gsl_status) break;
        error = gsl_multimin_fminimizer_size (state);
        gsl_status = gsl_multimin_test_size (error, tol);
        //        iter++;
      } while (gsl_status == GSL_CONTINUE && (++iter) < max_iter);
      // if(gsl_status) {                 // failure
      //   if(gsl_status==GSL_CONTINUE) { // MILD FAILURE
      //     //          out.status = false;
      //   } else {                  // HARD FAILURE
      //     statistic = error = std::numeric_limits<double>::infinity();
      //     //          out.status = false;
      //   }
      // } else {
      //   //        out.status = true;
      // }
      // if(gsl_status==GSL_SUCCESS || gsl_status==GSL_CONTINUE) {
      //   statistic=gsl_multimin_fminimizer_minimum(state);
      //   //        for(size_t i=0;i!=par.size();i++) par[i] = gsl_vector_get(state->x,i);
      //   //        error = extent;
      // }
    }
    // this is going to catch infinities
    catch(std::range_error const& e) {
      // we will exit as normally.
      // values are set to infinities to signal failed procedure
      fprintf(stderr,"%s %s : setting failure values;\n",e.what(),__FUNCTION__);
      //      out.status=false;
      //      out.gsl_code=;
      //      out.err=;
      //      out.statistic=std::numeric_limits<double>::infinity();
      //      for(auto& p: out.p) p= NAN ; // std::numeric_limits<double>::infinity();
    }
    statistic=gsl_multimin_fminimizer_minimum(state);
    for(size_t i=0;i!=dim;i++) theta[i]=gsl_vector_get(state->x,i);
    if(state!=NULL) gsl_multimin_fminimizer_free (state);
    //    gsl_vector_free(step);
    //    gsl_vector_free(var);
    // gsl_set_error_handler (oh);
    //    return out;
    return gsl_status;
  }

  // template<
  //   univariate_distribution_function Distribution
  //   >
  // int 
  // KL_mse_fmin(double& statistic, double *theta, const double step_size, const size_t dim,
  //             const double *samp, const size_t npt,
  //             gsl_multimin_fminimizer_type const *T,
  //             size_t& iter, const size_t max_iter,
  //             double& error, const double tol
  //             )
  // {
  //   assert(step_size>0);
  //   std::vector<double> step(dim);
  //   for(size_t i=0;i!=dim;i++) step[i]=step_size*theta[i];
  //   return KL_mse_fmin<Distribution>(statistic,theta,step,dim,samp,npt,T,iter,max_iter,error,tol);
  // }

    /** 
   * Find the minimizer for KL statistic. 
   * 
   * Version for single-parameter distributions and fixed step size multiplier.
   * 
   * @param[out] statistic 
   * @param[out] error 
   * @param[in/out] theta 
   * @param[in] step_size 
   * @param[in] samp 
   * @param[in] npt 
   * @param[in] T 
   * @param[in] tol 
   * @param[in] max_iter 
   * 
   * @return 
   */
  // template<
  //   univariate_distribution_function Distribution
  //   >
  // int 
  // KL_mse_fmin(double& statistic,double& theta, const double step_size,
  //             const double *samp, const size_t npt,
  //             gsl_multimin_fminimizer_type const *T,
  //             size_t& iter, const size_t max_iter,
  //             double& error, const double tol
  //             )
  // {
  //   assert(step_size>0);
  //   double step=theta*step_size;
  //   return KL_mse_fmin<Distribution>(statistic,&theta,&step,1,samp,npt,T,iter,max_iter,error,tol);
  // }


  // template<
  //   univariate_distribution_function Distribution
  //   >
  // int 
  // mKL_mse_fmin(double& statistic,
  //              double *theta, const double *step, const size_t dim,
  //              const double *samp, const size_t npt,
  //              gsl_multimin_fminimizer_type const *T,
  //              size_t& iter, const size_t max_iter,
  //              double& error, const double tol
  //              )
  // {
  //   int gsl_status;
  //   const GSL_msample p(a);
  //   gsl_vector_const_view gtheta = gsl_vector_const_view_array (theta, dim);
  //   gsl_vector_const_view gstep =  gsl_vector_const_view_array (step, dim);
  //   gsl_multimin_function mf;
  //   gsl_multimin_fminimizer *state;
  //   //    gsl_error_handler_t * oh = gsl_set_error_handler_off();
  //   try {
  //     mf.n=dim;
  //     mf.f=&KL_statistic<Distribution>;
  //     mf.params=(void*)&p;
  //     state = gsl_multimin_fminimizer_alloc (T, dim);
  //     gsl_multimin_fminimizer_set (state, &mf, &gtheta.vector, &gstep.vector);
  //     // out.name = gsl_multimin_fminimizer_name(state);
  //     // out.iter = 0;
  //     do {
  //       gsl_status = gsl_multimin_fminimizer_iterate(state);
  //       if (gsl_status) break;
  //       error = gsl_multimin_fminimizer_size (state);
  //       gsl_status = gsl_multimin_test_size (error, tol);
  //       //        iter++;
  //     } while (gsl_status == GSL_CONTINUE && (++iter) < max_iter);
  //     // if(gsl_status) {                 // failure
  //     //   if(gsl_status==GSL_CONTINUE) { // MILD FAILURE
  //     //     //          out.status = false;
  //     //   } else {                  // HARD FAILURE
  //     //     statistic = error = std::numeric_limits<double>::infinity();
  //     //     //          out.status = false;
  //     //   }
  //     // } else {
  //     //   //        out.status = true;
  //     // }
  //     // if(gsl_status==GSL_SUCCESS || gsl_status==GSL_CONTINUE) {
  //     //   statistic=gsl_multimin_fminimizer_minimum(state);
  //     //   //        for(size_t i=0;i!=par.size();i++) par[i] = gsl_vector_get(state->x,i);
  //     //   //        error = extent;
  //     // }
  //   }
  //   // this is going to catch infinities
  //   catch(std::range_error const& e) {
  //     // we will exit as normally.
  //     // values are set to infinities to signal failed procedure
  //     // std::cerr << e.what() << "\n"
  //     //           << __FUNCTION__
  //     //           <<  ": setting failure values;\n";
  //     //      out.status=false;
  //     //      out.gsl_code=;
  //     //      out.err=;
  //     //      out.statistic=std::numeric_limits<double>::infinity();
  //     //      for(auto& p: out.p) p= NAN ; // std::numeric_limits<double>::infinity();
  //   }
  //   statistic=gsl_multimin_fminimizer_minimum(state);
  //   for(size_t i=0;i!=dim;i++) theta[i]=gsl_vector_get(state->x,i);
  //   gsl_multimin_fminimizer_free (state);
  //   //    gsl_vector_free(step);
  //   //    gsl_vector_free(var);
  //   // gsl_set_error_handler (oh);
  //   //    return out;
  //   return gsl_status;
  // }









  
  
  

  // template <
  //   univariate_distribution_function F
  //   , univariate_distribution_function_d D
  //   >
  // void
  // KLentropy_tie_d(const gsl_vector *v, void *p, gsl_vector *g)
  // {
  //   info* i = (info*)p;
  //   double const * x_end=i->data + i->size;
  //   double const * x = i->data;
  //   double next,prev=0.0,tmp;
  //   double *prev_d = new double[g->size]{0};
  //   double *next_d = new double[g->size];
  //   //    prev=0; // F(*x,v->data);
  //   std::fill(prev_d,prev_d+g->size,0.0);
  //   //    D(*x,v->data,prev_d);
  //   gsl_vector_set_all (g,0.0);
  //   for(;x!=x_end;x++) {    // loop over data x+=1
  //     next = F(*x,v->data);
  //     if(next<0 || next > 1.0) throw std::range_error (std::string(__FUNCTION__) + std::string(": distribution function out of range"));
  //     if(next==NAN) throw std::range_error (std::string(__FUNCTION__) + std::string(": distribution function NaN"));
  //     D(*x,v->data,next_d);
  //     // look-ahead step
  //     //      int m=look_ahead(x,next,x_end);
  //     // look-ahead step
  //     int m=look_ahead(x,x_end);
  //     tmp  = double(m)/(next-prev);
  //     if(tmp==+std::numeric_limits<double>::infinity() || tmp==-std::numeric_limits<double>::infinity())
  //       throw std::range_error( std::string(__FUNCTION__) + std::string(" : infinite function value (most likely, bad data or bad parameters)"));
  //     for(size_t i=0;i!=g->size;i++) gsl_vector_set(g,i,(gsl_vector_get(g,i) + tmp*(next_d[i]-prev_d[i])));
  //     prev = next;
  //     std::swap(next_d,prev_d);
  //   }
  //   next = 1.0;
  //   tmp  = 1.0/(next-prev);
  //   if(tmp==+std::numeric_limits<double>::infinity() || tmp==-std::numeric_limits<double>::infinity())
  //     throw std::range_error( std::string(__FUNCTION__) + std::string(" : infinite function value (most likely, bad data or bad parameters)"));
  //   for(size_t i=0;i!=g->size;i++) gsl_vector_set(g,i,(gsl_vector_get(g,i) + tmp*(0 - prev_d[i])));
  //   // note the minus sign (confront with KLentropy and KLentropy_df)
  //   gsl_vector_scale(g,-1.0/double(1+i->size));
  //   delete[] next_d;
  //   delete[] prev_d;
  //   return;
  // }

  
  // template <
  //   univariate_distribution_function F
  //   , univariate_distribution_function_d D
  //   >
  // void
  // KLentropy_d(const gsl_vector *v, void *p, gsl_vector *g)
  // {
  //   info* i = (info*)p;
  //   auto * x = i->data, *x_end=i->data + i->size;
  //   double next,prev=0.0,tmp;
  //   double *prev_d = new double[g->size]{0};
  //   double *next_d = new double[g->size];
  //   //    prev=0; // F(*x,v->data);
  //   std::fill(prev_d,prev_d+g->size,0.0);
  //   //    D(*x,v->data,prev_d);
  //   gsl_vector_set_all (g,0.0);
  //   for(;x!=x_end;x++) {    // loop over data x+=1
  //     next = F(*x,v->data);
  //     D(*x,v->data,next_d);
  //     tmp  = 1.0/(next-prev);
  //     for(size_t i=0;i!=g->size;i++) gsl_vector_set(g,i,(gsl_vector_get(g,i) + tmp*(next_d[i]-prev_d[i])));
  //     prev = next;
  //     std::swap(next_d,prev_d);
  //   }
  //   next = 1.0;
  //   tmp  = 1.0/(next-prev);
  //   for(size_t i=0;i!=g->size;i++) gsl_vector_set(g,i,(gsl_vector_get(g,i) + tmp*(0 - prev_d[i])));    
  //   // note the minus sign (confront with KLentropy and KLentropy_fd)
  //   gsl_vector_scale(g,-1.0/double(1+i->size));
  //   delete[] next_d;
  //   delete[] prev_d;
  //   return;
  // }


  
  // template <
  //   univariate_distribution_function F
  //   , univariate_distribution_function_d D
  //   >
  // void
  // KLentropy_fd(const gsl_vector * v, void *p, double *f, gsl_vector *g) {
  //   assert(v->stride==1);
  //   info* i = (info*)p;
  //   auto * x = i->data, *x_end=i->data + i->size;
  //   double next,prev=0.0,tmp;
  //   double *prev_d = new double[g->size]{0};
  //   double *next_d = new double[g->size];
  //   std::fill(prev_d,prev_d+g->size,0.0);
  //   //    prev=F(*x,v->data);
  //   //    D(*x,v->data,prev_d);
  //   *f = 0;    
  //   gsl_vector_set_all (g, 0.0);
  //   for(;x!=x_end;x++) {    // loop over data x+=1
  //     next = F(*x,v->data);
  //     D(*x,v->data,next_d);
  //     tmp = 1.0/(next-prev);
  //     *f += log(next-prev);
  //     for(size_t i=0;i!=g->size;i++) gsl_vector_set(g,i,gsl_vector_get(g,i) + tmp*(next_d[i]-prev_d[i]));
  //     prev = next;
  //     std::swap(next_d,prev_d);
  //   }
  //   next = 1.0;
  //   tmp = 1.0/(next-prev);
  //   *f += log(next-prev);
  //   for(size_t i=0;i!=g->size;i++) gsl_vector_set(g,i,(gsl_vector_get(g,i) + tmp*(0 - prev_d[i])));
  //   // note the minus sign!  (confront with KLentropy and KLentropy_d)
  //   gsl_vector_scale(g,-1.0/double(1 + i->size));
  //   *f = - log(double(1 + i->size)) - ((*f)/double(1 + i->size));
  //   delete[] next_d;
  //   delete[] prev_d;
  //   return;
  // }

  
  // template <
  //   univariate_distribution_function F
  //   , univariate_distribution_function_d D
  //   >
  // void
  // KLentropy_tie_fd(const gsl_vector * v, void *p, double *f, gsl_vector *g) {
  //   assert(v->stride==1);
  //   info* i = (info*)p;
  //   double const * x_end=i->data + i->size;
  //   double const * x = i->data;
  //   double next,prev=0.0,tmp;
  //   double *prev_d = new double[g->size]{0};
  //   double *next_d = new double[g->size];
  //   std::fill(prev_d,prev_d+g->size,0.0);
  //   *f = 0;
  //   gsl_vector_set_all (g, 0.0);
  //   for(;x!=x_end;x++) {    // loop over data x+=1
  //     next = F(*x,v->data);
  //     if(next<0 || next > 1.0) throw std::range_error (std::string(__FUNCTION__) + std::string(": distribution function out of range"));
  //     if(next==NAN) throw std::range_error (std::string(__FUNCTION__) + std::string(": distribution function NaN"));
  //     D(*x,v->data,next_d);
  //     //      int m = look_ahead(x,next,x_end);
  //           // look-ahead step
  //     int m=look_ahead(x,x_end);
  //     tmp = double(m)/(next-prev);
  //     *f += double(m)*log((next-prev)/double(m));
  //     if(tmp==+std::numeric_limits<double>::infinity() || tmp==-std::numeric_limits<double>::infinity() ||
  //        *f ==+std::numeric_limits<double>::infinity() ||  *f==-std::numeric_limits<double>::infinity())
  //       throw std::range_error( std::string(__FUNCTION__) + std::string(" : infinite function value (resulting from, most likely, bad data or bad parameters"));
  //     for(size_t i=0;i!=g->size;i++) gsl_vector_set(g,i,gsl_vector_get(g,i) + tmp*(next_d[i]-prev_d[i]));
  //     prev = next;
  //     std::swap(next_d,prev_d);
  //   }
  //   next = 1.0;
  //   tmp = 1.0/(next-prev);
  //   *f += log(next-prev);
  //   if(tmp==+std::numeric_limits<double>::infinity() || tmp==-std::numeric_limits<double>::infinity() ||
  //      *f ==+std::numeric_limits<double>::infinity() ||  *f==-std::numeric_limits<double>::infinity())
  //     throw std::range_error( std::string(__FUNCTION__) + std::string(" : infinite function value (resulting from, most likely, bad data or bad parameters"));
  //   for(size_t i=0;i!=g->size;i++) gsl_vector_set(g,i,(gsl_vector_get(g,i) + tmp*(0 - prev_d[i])));
  //   // note the minus sign!  (confront with KLentropy and KLentropy_d)
  //   gsl_vector_scale(g,-1.0/double(1 + i->size));
  //   *f = - log(double(1 + i->size)) - ((*f)/double(1 + i->size));
  //   delete[] next_d;
  //   delete[] prev_d;
  //   return;
  // }







  
  // // template <univariate_distribution_function F
  // //           , discrete_probability G
  // //           , class D
  // //           >
  // // double
  // // multimode_KLentropy(const gsl_vector* v, void*p) {
  // template <class D>
  // double
  // multimode_KLentropy(void*p) {
  //   D* sample =*((D**)p);
  //   // std::vector<double> freq(sample->alternatives());
  //   // for(size_t i=0;i!=freq.size();i++) freq[i] = sample->frequency(i);
  //   double ans = 0.0, freq;
  //   for(size_t i=0;i!=sample->alternatives();i++) {
  //     freq = sample->frequency(i);
  //     //      std::pair<double*,double*> r = sample.range(i);
  //     // double tmp = freq * (log( freq/G(i,p)) +
  //     //                      KLentropy<F,G>(const gsl_vector* v, void*p)
  //     //                      tmp +=   template <univariate_distribution_function F>
  //     //                      double
                           

  //   }
  //   return ans;
  // }

  



  

  // /** 
  //  * Maximum Spacing Estimator. 
  //  * 
  //  * @param init Parameter initial values
  //  * @param data Observation data
  //  * 
  //  * @return 
  //  */
  // template<
  //   univariate_distribution_function Distribution
  //   >
  // outcome 
  // find_fminimizer(std::vector<double> const& init, const double *d, const size_t n, gsl_multimin_fminimizer_type const *T, double const tol=1.e-3)
  // {
  //   info p(d,n);
  //   outcome out(init.size());
  //   std::vector<double> &par = out.p;
  //   gsl_vector *var = gsl_vector_alloc (init.size());
  //   gsl_vector *step = gsl_vector_alloc (init.size());
  //   gsl_multimin_function mf;
  //   gsl_multimin_fminimizer *state;
  //   //    gsl_error_handler_t * oh = gsl_set_error_handler_off();
  //   try {
  //     for(size_t i=0;i!=init.size();++i) gsl_vector_set (var, i, init[i]);
  //     // seems to be a reasonable initial step
  //     for(size_t i=0;i!=init.size();i++) gsl_vector_set(step,i, 0.25*init[i]);
  //     mf.n=init.size();
  //     mf.f=&KLentropy_tie<Distribution>;
  //     mf.params=(void*)&p;
  //     state = gsl_multimin_fminimizer_alloc (T, init.size());
  //     gsl_multimin_fminimizer_set (state, &mf, var, step);
  //     out.name = gsl_multimin_fminimizer_name(state);
  //     out.iter = 0;
  //     double extent;
  //     do {
  //       out.gsl_status = gsl_multimin_fminimizer_iterate(state);
  //       if (out.gsl_status) break;
  //       extent = gsl_multimin_fminimizer_size (state);
  //       out.gsl_status = gsl_multimin_test_size (out.err, tol);
  //       out.iter++;
  //     } while (out.gsl_status == GSL_CONTINUE && out.iter < 100);
  //     if(out.gsl_status) {                 // failure
  //       if(out.gsl_status==GSL_CONTINUE) { // MILD FAILURE
  //         out.status = false;
  //       } else {                  // HARD FAILURE
  //         out.statistic = out.err = std::numeric_limits<double>::infinity();
  //         out.status = false;
  //       }
  //     } else {
  //       out.status = true;
  //     }
  //     if(out.gsl_status==GSL_SUCCESS || out.gsl_status==GSL_CONTINUE) {
  //       out.statistic=gsl_multimin_fminimizer_minimum(state);
  //       for(size_t i=0;i!=par.size();i++) par[i] = gsl_vector_get(state->x,i);
  //       out.err = extent;
  //     }
  //   }
  //   // this is going to catch infinities
  //   catch(std::range_error const& e) {
  //     // we will exit as normally.
  //     // values are set to infinities to signal failed procedure
  //     std::cerr << e.what() << "\n"
  //               << __FUNCTION__
  //               <<  ": setting failure values;\n";
  //     out.status=false;
  //     //      out.gsl_code=;
  //     //      out.err=;
  //     out.statistic=std::numeric_limits<double>::infinity();
  //     for(auto& p: out.p) p= NAN ; // std::numeric_limits<double>::infinity();
  //   }
  //   gsl_multimin_fminimizer_free (state);
  //   gsl_vector_free(step);
  //   gsl_vector_free(var);
  //   // gsl_set_error_handler (oh);
  //   return out;
  // }

  // template<
  //   univariate_distribution_function Distribution
  //   >
  // outcome 
  // find_fminimizer(std::vector<double> const& init, std::vector<double> const& data, gsl_multimin_fminimizer_type const *T, double const tol=1.e-3)
  // {
  //   return find_fminimizer<Distribution>(init,data.data(),data.size(),T,tol);
  // }
  
  // template<
  //   univariate_distribution_function Distribution
  //   , univariate_distribution_function_d Derivative
  //   >
  // outcome
  // find_fdfminimizer(std::vector<double> const& init, std::vector<double> const& data, const gsl_multimin_fdfminimizer_type * T, double const tol=1.e-4)
  // {
  //   info p(data);
  //   outcome out(init.size());
  //   std::vector<double>& par = out.p;
  //   gsl_vector *var = gsl_vector_alloc (init.size());
  //   gsl_multimin_function_fdf mf;
  //   gsl_multimin_fdfminimizer *state;
  //   try{ 
  //     //    gsl_error_handler_t * oh = gsl_set_error_handler_off();    
  //     for(size_t i=0;i!=init.size();++i) gsl_vector_set (var, i, init[i]);
  //     mf.n=init.size();
  //     mf.f=KLentropy_tie<Distribution>;
  //     mf.df=KLentropy_tie_d<Distribution,Derivative>;
  //     mf.fdf=KLentropy_tie_fd<Distribution,Derivative>; 
  //     mf.params=(void*)&p;
  //     state = gsl_multimin_fdfminimizer_alloc (T, init.size());
  //     gsl_multimin_fdfminimizer_set (state, &mf, var, 0.1, tol);
  //     out.name = gsl_multimin_fdfminimizer_name(state);
  //     //      int iter = 0, status;
  //     double extent;
  //     do {
  //       out.gsl_status = gsl_multimin_fdfminimizer_iterate(state);
  //       if(out.gsl_status) break;
  //       out.gsl_status = gsl_multimin_test_gradient (state->gradient,tol);
  //       out.iter++;
  //     } while (out.gsl_status == GSL_CONTINUE && out.iter < 100);
  //     if(out.gsl_status) { // failure
  //       if(out.gsl_status==GSL_CONTINUE) { // Mild failure
  //         out.status=false;
  //       } else {
  //         out.statistic=out.err=std::numeric_limits<double>::infinity();
  //         out.status=false;
  //       }
  //     } else {
  //       out.status=false;
  //     }
  //     if(out.gsl_status==GSL_SUCCESS || out.gsl_status==GSL_CONTINUE) {
  //       out.statistic=gsl_multimin_fdfminimizer_minimum(state);
  //       for(size_t i=0;i!=par.size();i++) par[i] = gsl_vector_get(state->x,i);
  //       out.err=0.0;
  //       for(size_t i=0;i!=par.size();++i) out.err += fabs(gsl_vector_get(state->gradient,i));
  //     }
  //   }
  //   // this is going to catch infinities
  //   catch(std::range_error const& e) {
  //     // we will exit as normally.
  //     // values are set to infinities to signal failed procedure
  //     std::cerr << e.what() << "\n"
  //               << __FUNCTION__
  //               <<  ": setting failure values;\n";
  //     out.status=false;
  //     //      out.gsl_code=;
  //     //      out.err=;
  //     out.statistic=std::numeric_limits<double>::infinity();
  //     for(auto& p: out.p) p=NAN; // std::numeric_limits<double>::infinity();
  //   }
  //   gsl_multimin_fdfminimizer_free (state);
  //   gsl_vector_free(var);
  //   //    gsl_set_error_handler (oh);
  //   return out;
  // }





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

