// Exact two-vertex induced extensions of the specified determinant modules.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <vector>
using I=std::int64_t; using U=std::uint64_t;
struct P {I h,p; unsigned a,b,e;};
using Bank=std::unordered_map<U,P>;
const std::array<unsigned,8> L={4,9,82,21,35,143,187,31};
const std::array<unsigned,8> R={6,12,48,21,131,219,159,15};
U key(I h,I p){if(h<0||p<0||h>0xffffffffLL||p>0xffffffffLL)throw std::runtime_error("key range");return(U(h)<<32)|U(p);}
std::pair<I,I> profile(const unsigned* rows,unsigned n){
 std::array<std::array<I,10>,1024> dp{};std::array<I,1024> hh{};
 unsigned full=(1u<<n)-1;hh[0]=1;
 for(unsigned i=0;i<n;i++)dp[1u<<i][i]=1;
 for(unsigned s=1;s<=full;s++)for(unsigned v=0;v<n;v++){
  I z=dp[s][v];hh[s]+=z;if(!z)continue;
  for(unsigned rem=rows[v]&~s;rem;rem&=rem-1){unsigned b=rem&-rem;dp[s|b][__builtin_ctz(b)]+=z;}
 }
 I p=0;for(unsigned s=0;s<=full;s++)p+=hh[s]*hh[full^s];
 return {hh[full],p};
}
Bank build(const std::array<unsigned,8>& base){
 Bank out;out.reserve(40000);unsigned rows[10];
 for(unsigned a=0;a<256;a++)for(unsigned b=0;b<256;b++)for(unsigned e=0;e<2;e++){
  for(unsigned i=0;i<8;i++)rows[i]=base[i]|((a>>i&1)<<8)|((b>>i&1)<<9);
  rows[8]=(255^a)|(e<<9);rows[9]=(255^b)|((1-e)<<8);
  auto [h,p]=profile(rows,10);out.emplace(key(h,p),P{h,p,a,b,e});
 }
 return out;
}
I inverse(I a,I m){I r=a,s=m,x=1,y=0;while(s){I q=r/s,t=r-q*s;r=s;s=t;t=x-q*y;x=y;y=t;}if(r!=1)throw std::runtime_error("inverse");x%=m;return x<0?x+m:x;}
int main(int argc,char**argv){
 try{
  if(profile(L.data(),8)!=std::make_pair<I,I>(181,3424)||profile(R.data(),8)!=std::make_pair<I,I>(179,3416))throw std::runtime_error("base controls");
  auto left=build(L);std::cerr<<"left profiles "<<left.size()<<"\n";
  auto right=build(R);std::cerr<<"right profiles "<<right.size()<<"\n";
  I maxj=0;for(auto&[k,t]:right)maxj=std::max(maxj,t.h);
  std::vector<std::pair<P,P>> hits;
  for(auto&[k,m]:left){
   if(m.h%3==0||std::gcd(m.h,m.p)!=1)continue;
   I x=6*m.p+m.h,mod=6*m.h,inv=inverse(x%mod,mod);
   for(I e:{-1,1}){
    I j=(e==1?inv:mod-inv);if(!j)j=mod;
    for(;j<=maxj;j+=mod){I q=(x*j-e)/mod;if(q%2)continue;auto it=right.find(key(j,q));if(it!=right.end())hits.emplace_back(m,it->second);}
   }
  }
  std::sort(hits.begin(),hits.end(),[](auto&a,auto&b){return std::tie(a.first.h,a.first.p,a.second.h,a.second.p)<std::tie(b.first.h,b.first.p,b.second.h,b.second.p);});
  std::ofstream file(argc>1?argv[1]:"pair_growth2.json");
  file<<"{\n\"status\":\"FINITE-EXACT\",\n\"universe\":\"All 131072 labelled induced two-vertex extensions per supplied order-eight module; all ordered response pairs via exact modular inverse\",\n\"left_profiles\":"<<left.size()<<",\n\"right_profiles\":"<<right.size()<<",\n\"hit_count\":"<<hits.size()<<",\n\"hits\":[";
  bool first=true;for(auto&pair:hits){auto m=pair.first,n=pair.second;I e=(6*m.p+m.h)*n.h-6*m.h*n.p;if(std::abs(e)!=1)throw std::runtime_error("hit");
   if(!first)file<<",";first=false;
   file<<"\n{\"M\":["<<m.h<<","<<m.p<<"],\"N\":["<<n.h<<","<<n.p<<"],\"left_attachment\":["<<m.a<<","<<m.b<<","<<m.e<<"],\"right_attachment\":["<<n.a<<","<<n.b<<","<<n.e<<"],\"e\":"<<e<<"}";
  }file<<"\n]}\n";std::cout<<"PASS exact two-vertex universes; hits "<<hits.size()<<"\n";
 }catch(std::exception&e){std::cerr<<e.what()<<"\n";return 1;}
}
