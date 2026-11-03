#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>
using I=std::int64_t;
I hp(const std::vector<unsigned>& a){
 int n=int(a.size());std::vector<I> dp(std::size_t(1u<<n)*n,0);
 for(int v=0;v<n;v++)dp[(1u<<v)*n+v]=1;
 for(unsigned m=1;m<(1u<<n);m++)for(int v=0;v<n;v++){
  I c=dp[m*n+v];if(!c)continue;unsigned next=a[v]&~m;
  while(next){int w=__builtin_ctz(next);next&=next-1;dp[(m|(1u<<w))*n+w]+=c;}
 }
 I h=0;for(int v=0;v<n;v++)h+=dp[((1u<<n)-1)*n+v];return h;
}
int main(){
 std::vector<std::vector<unsigned>> contexts={{}};
 for(unsigned a=0;a<64;a++)contexts.push_back({a});
 for(unsigned a=0;a<16;a++)contexts.push_back({(13*a)%64,(37*a+9)%128});
 for(unsigned a=0;a<8;a++)contexts.push_back({(7*a+3)%64,(19*a+5)%128,(53*a+11)%256});
 bool first=true;std::cout<<"{\"contexts\":[";
 for(const auto &ctx:contexts){
  I sum=0,c1=0,c2=0,q1=0,q2=0;
  for(unsigned code=0;code<32768;code++){
   std::vector<unsigned>a(6,0);int e=0;
   for(int i=0;i<6;i++)for(int j=i+1;j<6;j++,e++)
    if(code&(1u<<e))a[i]|=1u<<j;else a[j]|=1u<<i;
   for(unsigned cut:ctx){int n=int(a.size());a.push_back(cut);
    for(int i=0;i<n;i++)if(!(cut&(1u<<i)))a[i]|=1u<<n;}
   I h=hp(a);if(h<1 || !(h&1))return 2;
   I z1=(__builtin_popcount(code&10376)&1)?-1:1;
   I z2=(__builtin_popcount(code&5392)&1)?-1:1;
   sum+=h;c1+=z1*h;c2+=z2*h;q1+=z1*h*h;q2+=z2*h*h;
  }
  if(c1!=c2)return 3;
  if(!first)std::cout<<",";first=false;
  std::cout<<"{\"cuts\":[";for(std::size_t i=0;i<ctx.size();i++){if(i)std::cout<<",";std::cout<<ctx[i];}
  std::cout<<"],\"sum_H\":"<<sum<<",\"sum_chi1_H\":"<<c1<<",\"sum_chi2_H\":"<<c2
           <<",\"sum_chi1_H2\":"<<q1<<",\"sum_chi2_H2\":"<<q2<<"}";
 }
 std::cout<<"],\"all_old_orientations\":32768,\"graphs_counted\":"<<32768*contexts.size()
          <<",\"algorithm\":\"literal Boolean Hamiltonian DP; independent of forest contraction\"}\n";
}
