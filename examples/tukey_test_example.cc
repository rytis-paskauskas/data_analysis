#include "clusters/tukey_test.hpp"

#include <utility>
#include <vector>

typedef std::pair<int,double> event;

double smooth_fcn(double x) {
  return x*x;
}

int main(void) {

  size_t n=30;
  std::vector<event> data(n);
  for(size_t i=0;i!=n;i++) data[i]={i,smooth_fcn( double(i)/double(n) - 0.5 )};
  data[n/2].second += 10.0;

  auto extractor = [] ( const event& e) -> double {
    return e.second;
  };

  auto outlier_test = tukey_test(data.begin(),data.end(),extractor,3.0,10);
  for(auto&m:outlier_test) {
    event& d=*(m.second-1);
    printf("One before the outlier %d %f\n",d.first,d.second);
  }
  return 0;
}
