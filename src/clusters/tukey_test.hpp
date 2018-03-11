/**
 * @file   tukey_test.hpp
 * @author Rytis Paskauskas <rytis@W530>
 * @date   Thu Jul  6 19:26:23 2017
 *
 * (C) Rytis Paskauskas rytis.paskauskas@gmail.com		2017-2018
 *
 * @brief  Tukey's outlier test.
 * https://en.wikipedia.org/wiki/Outlier#Tukey.27s_test
 * https://en.wikipedia.org/wiki/Tukey%27s_range_test

 * Tukey's range test, also known as the Tukey's test, Tukey method,
 * Tukey's honest significance test, Tukey's HSD (honest significant
 * difference) test, or the Tukey–Kramer method, is a single-step
 * multiple comparison procedure and statistical test. Tukey's test
 * compares the means of every treatment to the means of every other
 * treatment; that is, it applies simultaneously to the set of all
 * pairwise comparisons μi − μj.

 * In statistics, an outlier is an observation point that is distant
 * from other observations. An outlier may be due to variability in
 * the measurement or it may indicate experimental error; the latter
 * are sometimes excluded from the data set.

 * Some methods flag observations based on measures such as the
 * interquartile range. For example, if Q1 and Q3 are the lower and
 * upper quartiles respectively, then one could define an outlier to
 * be any observation outside the range: [ Q1 − k ( Q3 − Q1 ) , Q3 + k
 * ( Q3 − Q1 ) ] for some nonnegative constant k. John Tukey proposed
 * this test, where k = 1.5 indicates an "outlier", and k = 3
 * indicates data that is "far out".  I use a modification of Tukey's
 * test with larger constant k because the distribution of deviations
 * is heavy tailed. As default value, k=10.
 *
 */

#ifndef tukey_test_hpp_defined
#define tukey_test_hpp_defined
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/framework/accumulator_set.hpp>
#include <boost/accumulators/statistics/p_square_quantile.hpp>
#include <type_traits>
#include <iterator>
#include <map>

/** 
 * Let's define an "floating type equal" function
 * This function is used to capture uninitialized quantiles in the tukey test.
 * 
 * @param[in] x a float to be compared
 * @param[in] y a float to be compared against
 * @param[in] ulp the machine epsilon has to be scaled to the magnitude of the values used and multiplied by the desired precision in ULPs (units in the last place)
 * 
 * @return true if relative error is less than ULP last siginificant digits.
 */
template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
almost_equal(T x, T y, int ulp)
{
  return std::abs(x-y) <= std::numeric_limits<T>::epsilon() * std::abs(x+y) * ulp
    // unless the result is subnormal
    || std::abs(x-y) < std::numeric_limits<T>::min();
}

/** 
 * Tukey test that returns extreme values' map to data iterators.
 * See tukey_test_example.cc for an example.
 * @note Skip must be at least 5 but possibly even more because of how
 * p_square_quantile is constructed. I believe it does not initialize
 * unless at least 5 data are accumulated. It may require more data
 * actually. For this reason, we monitor the quantiles and may force
 * additional data to be skipped. This is a hardened version.
 * @param[in] first data range
 * @param[in] last  data range
 * @param[in] value this tells us how to extract a number from the
 * iterator.
 * @param[in] farout This is the parameter of the Tukey test. For
 * example, if farout=3, data that are three times the interquantile
 * range are considered as outliers.
 * @param[in] skip how many initial data to skip. Note that skipping
 * more data than specified can be imposed if the quantiles aren't
 * correctly initialized.
 * @return A map to that is of the form [extreme value] -> [data
 * iterator]. This function returns an empty map if first==last
 */
template <
  typename K
  , class Iterator
  , class Extractor
  , typename V=typename std::iterator_traits<Iterator>::value_type
  , class Compare = std::greater<K>
  >
std::map< K, Iterator, Compare >
tukey_test(Iterator first, Iterator last, Extractor value, K farout, size_t skip = 5) // use farout = 3.0 as a reasonable default option...
{
  namespace ba = boost::accumulators;
  static_assert(std::is_same<std::random_access_iterator_tag, typename std::iterator_traits<Iterator>::iterator_category>::value
                , "outlier_turkey_test error: Iterator must be a random access iterator or a raw pointers to an array.\n");
  assert(first + skip < last);
  std::map<K,Iterator,Compare> res;
  ba::accumulator_set
    <
      K
    , ba::stats<ba::tag::p_square_quantile> // p_square_quantile may need at least 5 readings before it can compute quantiles.
    > Q1(ba::quantile_probability = 0.25)
    , Q3(ba::quantile_probability = 0.75);
  K v,vo,r,q1,q3,iqr;
  vo=value(*first++);
  while(skip!=0)
    {
      skip--;
      v=value(*first++);
      r = v-vo;
      Q1(r);
      Q3(r);
      vo = v;
    }
  // this is an additional safeguard to avoid uninitialized quantiles
  while(almost_equal(ba::p_square_quantile(Q1),ba::p_square_quantile(Q3),2) && first!=last)
    {
      v=value(*first++);
      r = v-vo;
      Q1(r);
      Q3(r);
      vo = v;
    }
  // main part of the loop
  // returns an empty map if first==last
  for(;first!=last;first++)
    {
      v=value(*first);
      r = v-vo;
      q1 = ba::p_square_quantile(Q1);
      q3 = ba::p_square_quantile(Q3);
      iqr = farout*(q3-q1);
      if( r<=q1-iqr || r>=q3+iqr )
        //        res[(r-((q1+q3)/2))/farout/(q3-q1)] = *(first-1);
        res[(r-((q1+q3)/2))/iqr] = first;
      else {
        vo = v;
        Q1(r);
        Q3(r);
      }
    }
  return res;
}

#endif
