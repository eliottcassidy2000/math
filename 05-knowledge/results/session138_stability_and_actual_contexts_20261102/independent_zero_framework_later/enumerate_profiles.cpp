// Independent complete small-tournament profile census.
// Forward endpoint-extension DP counts Hamilton paths on each subset;
// an unrestricted first-block convolution counts ordered path partitions.
// No permutation/descent census, OCF formula, or producer code is used.
#include <array>
#include <cstdint>
#include <iostream>
#include <set>
#include <stdexcept>
#include <vector>
using I=std::int64_t;
using V=std::vector<I>;
void need(bool b,const char* s){if(!b)throw std::runtime_error(s);}
V profile(const std::vector<unsigned>& a){
 const int n=static_cast<int>(a.size()),N=1<<n;
 std::vector<std::array<I,6>> ends(N);std::vector<I> h(N);
 for(int j=0;j<n;++j)ends[1<<j][j]=1;
 for(int S=1;S<N;++S){
  for(int j=0;j<n;++j){h[S]+=ends[S][j];unsigned next=a[j]&unsigned((N-1)^S);
   while(next){unsigned bit=next&(~next+1u);next-=bit;int k=__builtin_ctz(bit);ends[S|int(bit)][k]+=ends[S][j];}
  }
 }
 std::vector<std::array<I,7>> g(N);g[0][0]=1;
 for(int S=1;S<N;++S)for(int T=S;T;T=(T-1)&S)
  for(int r=1;r<=n;++r)g[S][r]+=h[T]*g[S^T][r-1];
 V out;for(int r=1;r<=n;++r)out.push_back(g[N-1][r]);return out;
}
int main(){I total=0,pairsum=0,profilecount=0,impliedwords=0;std::cout<<"{\"status\":\"PASS\",\"method\":\"forward endpoint-extension and ordered first-block convolution\",\"orders\":[";
 for(int n=1;n<=6;++n){unsigned labels=1u<<(n*(n-1)/2);std::set<V> fs;
  for(unsigned code=0;code<labels;++code){std::vector<unsigned>a(n);unsigned bit=1;
   for(int j=1;j<n;++j)for(int i=0;i<j;++i){if(code&bit)a[j]|=1u<<i;else a[i]|=1u<<j;bit<<=1;}
   fs.insert(profile(a));
  }
  I pairs=0;
  for(auto it=fs.begin();it!=fs.end();++it){for(int i=1;i+1<n;++i)need((*it)[i]*(*it)[i]>=(*it)[i-1]*(*it)[i+1],"LC failure");
   for(auto jt=std::next(it);jt!=fs.end();++jt){V d(n+3);for(int i=0;i<n;++i)d[i+1]=(*jt)[i]-(*it)[i];bool positive=false,negative=false;
    for(int i=1;i<=n+1;++i){I diagonal=d[i-1]+3*d[i]+d[i+1];positive|=diagonal>0;negative|=diagonal<0;if(i<=n){I off=d[i]+d[i+1];positive|=off>0;negative|=off<0;}}
    need(!(positive&&negative),"mixed kernel found");++pairs;
   }
  }
  I factorial=1;for(int j=2;j<=n;++j)factorial*=j;impliedwords+=labels*factorial;
  total+=labels;pairsum+=pairs;profilecount+=static_cast<I>(fs.size());
  if(n>1)std::cout<<',';
  std::cout<<"{\"n\":"<<n<<",\"labelled_tournaments\":"<<labels<<",\"profile_count\":"<<fs.size()<<",\"unordered_profile_pairs\":"<<pairs<<",\"all_kernels_one_signed\":true,\"all_profiles_log_concave\":true,\"profiles\":[";
  bool first=true;for(const auto& f:fs){if(!first)std::cout<<',';first=false;std::cout<<'[';for(int i=0;i<n;++i){if(i)std::cout<<',';std::cout<<f[i];}std::cout<<']';}std::cout<<"]}";
 }
 need(total==33867&&pairsum==308&&profilecount==39&&impliedwords==23717429,"total mismatch");
 std::cout<<"],\"total_tournaments\":"<<total<<",\"total_profiles\":"<<profilecount<<",\"total_profile_pairs\":"<<pairsum<<",\"implied_permutation_universe\":"<<impliedwords<<"}\n";
}
