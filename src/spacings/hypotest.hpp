
// pi^2/6 - 1
#define MSE_SIG2_INF 0.64493406684822643647
// SQRT of the above
#define MSE_ASYMPTOTIC_ROOTVAR  0.80307787097405842818

namespace spacings {
  
  /** 
   * Large sample central limit.
   * The large sample central limit transformation of \f$X_n\f$ is \f$Y_n = \frac{\sqrt{n}X_n}{\sigma}\f$. In theory, \f$\displaystyle\lim_{n\to\infty}Y_n=\mathsf{N}(0,1)\f$.
   *
   * @param[in] f KL statistic \f$S_n(\theta)\f$
   * @param[in] n sample size (number of data; the \f$n\f$ in \f$S_n(\theta)\f$).
   * 
   * @return A transformed input random variable, distributed approximately as \f$\displaystyle\mathrm{N}(0,1)\f$.
   */
  double normal_cl_statistic(double f, size_t n) {
    double sm = MSE_ASYMPTOTIC_ROOTVAR/sqrt(double(n+1));
    //    return (f - M_EULER)/sm;
    return f/sm;
  }

  double normal_pvalue(double f, size_t n) {
    return gsl_cdf_gaussian_Q(f,1);
  }

  /** 
   * Small sample central limit.
   *
   * The small sample central limit transformation of \f$X_n\f$ is 
   * \f$\displaystyle Y_n = \frac{\sqrt{n+1}(X_n-\gamma_n)}{\sigma_n}\f$.
   * In theory, \f$\displaystyle Y_n\approx\frac{\chi^2_n-n}{\sqrt{2n}}\f$.
   * 
   * @param[in] f KL statistic \f$S_n(\theta)\f$.
   * @param[in] n size of the sample (number of data).
   * @param[in] k \f$=0\f$ if \f$\theta\f$ is known, otherwise \f$k=\dim{\theta}\f$.
   * 
   * @return A transformed input random variable, distributed approximately as \f$\displaystyle \frac{\chi^2_n-n}{\sqrt{2n}}\f$.
   */
  double chisq_cl_statistic(double f, size_t n, int const k=0) {
    if(f==+std::numeric_limits<double>::infinity() ||
       f==-std::numeric_limits<double>::infinity()) 
      return std::numeric_limits<double>::infinity();
    double m = 1.0/double(n+1);
    //    double gm = M_EULER - m * (0.5 + (0.08333333333333333333 * m)); // asymptotic value + O(m^3) of KL_statistic
    double gm = - m * (0.5*double(k) + 0.5 + (0.08333333333333333333 * m)); // asymptotic value + O(m^3) of KL_statistic
    double sm = sqrt( m *( MSE_SIG2_INF   - m * ( 0.5 +  0.16666666666666666667 * m )) ); // asymptotic variance 
    //    double c1 = gm - sqrt(0.5*double(n)) * sm;
    //     double c2 = sqrt (0.5/double(n)) * sm;
    // recycle m;
    //    m = (f - c1 + 0.5*double(k) * m)/c2;
    //    return (f - gm + (m * 0.5*double(k)))/sm; // should yield (chi^2(n) - n)/sqrt(2n)
    return (f - gm)/sm; // should yield (chi^2(n) - n)/sqrt(2n)
  }
  
  /** 
   * Returns the p-value assuming that input is a statistic given by chisq_cl
   * 
   * @param f 
   * @param n 
   * @param g 
   * 
   * @return 
   */
  double chisq_pvalue(double f, size_t n) {
    if(f==+std::numeric_limits<double>::infinity() ||
       f==-std::numeric_limits<double>::infinity()) return 0;
    // double m = 1.0/double(n+1);
    // double gm = M_EULER - m * (0.5 + (0.08333333333333333333 * m));
    // double sm = sqrt( m *( MSE_SIG2_INF   - m * ( 0.5 +  0.16666666666666666667 * m )) );
    // double c1 = gm - sqrt(0.5*double(n)) * sm;
    // double c2 = sqrt (0.5/double(n)) * sm;
    // // recycle m;
    // m = (f - c1 + 0.5*double(k) * m)/c2;
    // return gsl_cdf_chisq_Q(m,n);
    return gsl_cdf_chisq_Q(sqrt(2*n)*f + n, n);
  }

}
