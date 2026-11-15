#include <array>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
using std::int64_t;
constexpr int B=2855,Q=5760,T=4034,S=3914,h=45,M=16;
struct Counts{int64_t total=0;std::array<int64_t,M> c{};};
Counts enumerate(int64_t X,int64_t limit_x,int64_t limit_y,const std::vector<int>& p){
  Counts r;
  for(int64_t x=0;x<=std::min(limit_x,(X-B)/T);++x){
    int64_t b=B+T*x+h*p[x], a=Q*x+S;
    if(b>X) continue;
    for(int64_t y=0;y<=std::min(limit_y,(X-b)/a);++y){
      int64_t n=b+a*y+h*p[y];
      if(n>X) continue;
      if(n%2!=1)throw std::runtime_error("parity");
      ++r.total;++r.c[(n/2)%M];
    }
  }
  return r;
}
int main(){
  const int64_t bound=1000000000;
  std::vector<int> p(bound/std::min(T,S)+2);
  for(std::size_t x=1;x<p.size();++x){int u=x/2;p[x]=p[u]+(x%2?0:u+u%2);}
  auto c=[](double a) {return (std::lgamma(a/Q)-std::lgamma((a+h)/Q))/h;};
  double K=c(S)+c(T)-(1+std::log(double(Q)))/Q;
  std::cout << std::fixed << std::setprecision(10) << "K " << K << '\n';
  for(int64_t X : {1000000, 10000000, 100000000, 1000000000}){
    int64_t L=std::sqrt(X/double(Q));
    while(Q*(L+1)*(L+1)<=X)++L;
    while(Q*L*L>X)--L;
    auto full=enumerate(X,X,X,p), v=enumerate(X,L,X,p), hor=enumerate(X,X,L,p), sq=enumerate(X,L,L,p);
    if(full.total!=v.total+hor.total-sq.total || (L+1)*(L+1)-sq.total>16)
      throw std::runtime_error("hyperbola decomposition");
    for(int r=0;r<M;++r)if(full.c[r]!=v.c[r]+hor.c[r]-sq.c[r])throw std::runtime_error("modular decomposition");
    double normalized_error=(full.total-X/double(Q)*std::log(double(X))-K*X)/std::sqrt(double(X));
    auto lohi=std::minmax_element(full.c.begin(),full.c.end());
    std::cout << X << " " << full.total << " " << *lohi.first << " " << *lohi.second
       << " " << std::setprecision(6) << normalized_error << '\n';
  }
}
