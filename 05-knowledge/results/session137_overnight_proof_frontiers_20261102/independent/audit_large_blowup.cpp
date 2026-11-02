// Independent selected k6 row: forward subset endpoint DP, reverse endpoints,
// and direct endpoint insertion polynomial, NOT orientation-capacity evaluation.
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>
using I=std::int64_t;using Graph=std::vector<unsigned>;
void require(bool b,const char*m){if(!b)throw std::runtime_error(m);}
std::vector<I> endpoints(const Graph&a){
  int n=a.size();unsigned N=1u<<n;std::vector<I> dp(std::size_t(N)*n);
  for(int i=0;i<n;++i)dp[std::size_t(1u<<i)*n+i]=1;
  for(unsigned s=1;s<N;++s)for(int last=0;last<n;++last){I value=dp[std::size_t(s)*n+last];if(!value)continue;unsigned next=a[last]&~s;
    while(next){int j=__builtin_ctz(next);next&=next-1;dp[std::size_t(s|1u<<j)*n+j]+=value;}}
  return dp;
}
I hp(const Graph&a){auto dp=endpoints(a);int n=a.size();I total=0;for(int i=0;i<n;++i)total+=dp[std::size_t((1u<<n)-1)*n+i];return total;}
Graph add(const Graph&a,unsigned cut){int n=a.size();Graph b=a;for(int i=0;i<n;++i)if(!(cut>>i&1))b[i]|=1u<<n;b.push_back(cut);return b;}
int main(){
  Graph q(6);unsigned code=16;for(int i=0;i<6;++i)for(int j=i+1;j<6;++j){if(code&1)q[i]|=1u<<j;else q[j]|=1u<<i;code>>=1;}
  Graph T(18);for(int b=0;b<6;++b)for(int i=0;i<3;++i){T[3*b+i]|=1u<<(3*b+(i+1)%3);for(int c=0;c<6;++c)if(q[b]>>c&1)T[3*b+i]|=7u<<(3*c);}
  I target=hp(T);require(target==1049305833,"full order18 target");
  int n=17;unsigned N=1u<<n,full=N-1;Graph a(n),reverse(n);unsigned original=0;
  for(int i=0;i<n;++i){if(T[0]>>(i+1)&1)original|=1u<<i;for(int j=0;j<n;++j)if(T[i+1]>>(j+1)&1)a[i]|=1u<<j;}
  for(int i=0;i<n;++i)for(int j=0;j<n;++j)if(a[j]>>i&1)reverse[i]|=1u<<j;
  auto E=endpoints(a),S=endpoints(reverse);std::vector<I> Q(n*n),linear(n),values(N);I h=0;
  for(int i=0;i<n;++i)h+=E[std::size_t(full)*n+i];
  for(unsigned s=1;s<full;++s){unsigned r=full^s;
    for(int i=0;i<n;++i)if(s>>i&1){I e=E[std::size_t(s)*n+i];if(!e)continue;
      for(int j=0;j<n;++j)if(r>>j&1)Q[i*n+j]+=e*S[std::size_t(r)*n+j];}}
  for(int i=0;i<n;++i){linear[i]=S[std::size_t(full)*n+i]-E[std::size_t(full)*n+i];for(int j=0;j<n;++j)linear[i]+=Q[j*n+i];}
  values[0]=h;I divisor=0;std::vector<unsigned> zero;
  for(unsigned cut=1;cut<N;++cut){int i=__builtin_ctz(cut);unsigned rest=cut^(1u<<i);I value=values[rest]+linear[i];
    for(int j=0;j<n;++j) {
      if(rest>>j&1)value-=Q[i*n+j]+Q[j*n+i];
    }
    values[cut]=value;
  }
  for(unsigned cut=0;cut<N;++cut){divisor=std::gcd(divisor,values[cut]-h);if(values[cut]==target)zero.push_back(cut);require(std::abs(values[cut]-target)!=2,"no unit replacement");}
  require(h==28273779&&divisor==54&&original==114689&&zero==std::vector<unsigned>{original},"selected complete row");
  for(unsigned cut:std::vector<unsigned>{0,full,original,original^1u})require(hp(add(a,cut))==values[cut],"direct order18 child control");
  std::cout<<"{\n  \"scope\": \"one independently recomputed k6 quotient code16, deleted vertex0; all131072 literal cuts, no orbit compression\",\n  \"source_order\": 18,\n  \"parent_order\": 17,\n  \"target_H\": "<<target<<",\n  \"parent_H\": "<<h<<",\n  \"response_gcd\": "<<divisor<<",\n  \"zero_cuts\": ["<<original<<"],\n  \"direct_child_controls\": 4,\n  \"verdict\": \"PASS unique original zero and no unit response; not the entire k6 census\"\n}\n";
}
