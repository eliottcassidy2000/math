// Exact second counted-complement moment of the labelled 81/27 hostile.
// No canonical-labelling library: normalized path presentations are literal bits.
#include <array>
#include <cstdint>
#include <iostream>
#include <set>
#include <stdexcept>
#include <vector>
using U=uint64_t; using T=std::array<unsigned,12>;
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
U paths(const T&t){
 std::vector<std::array<U,12>> dp(4096);
 for(int v=0;v<12;v++)dp[1<<v][v]=1;
 for(unsigned s=1;s<4096;s++)for(int v=0;v<12;v++)if(dp[s][v]){
   unsigned avail=t[v]&(4095^s);
   while(avail){unsigned b=avail&-avail;avail-=b;int w=__builtin_ctz(b);dp[s|b][w]+=dp[s][v];}}
 U h=0;for(auto x:dp[4095])h+=x;return h;
}
T build(bool minus){T t{};
 for(int b=0;b<3;b++){
  for(int i=0;i<3;i++){t[4*b+i]|=1u<<(4*b+(i+1)%3);if(minus&&b==2)t[4*b+3]|=1u<<(4*b+i);else t[4*b+i]|=1u<<(4*b+3);}
  for(int i=0;i<4;i++)t[4*b+i]|=15u<<(4*((b+1)%3));}
 return t;
}
int main(){std::cout<<"[\n";
 for(int minus=0;minus<2;minus++){
  auto t=build(minus);std::array<int,12> p{};std::set<U> codes;U raw=0;
  auto dfs=[&](auto&&self,unsigned used,int len)->void{
   if(len==12){raw++;U code=0;int k=0;for(int i=0;i<12;i++)for(int j=i+2;j<12;j++)code|=U((t[p[i]]>>p[j])&1)<<k++;codes.insert(code);return;}
   unsigned avail=(len?t[p[len-1]]:4095)&(4095^used);
   while(avail){unsigned b=avail&-avail;avail-=b;p[len]=__builtin_ctz(b);self(self,used|b,len+1);}
  };dfs(dfs,0,0);
  need(raw==paths(t)&&raw==462915,"source paths");
  unsigned aut=minus?27:81;need(codes.size()*aut==raw,"presentation mass");
  U sum=0,minh=~U(0),maxh=0;std::vector<std::pair<U,U>> sample;
  for(U c:codes){T partner{};int k=0;
   for(int i=0;i<12;i++)for(int j=i+1;j<12;j++){
    bool forward=j==i+1 || !((c>>k)&1);if(j>i+1)k++;
    if(forward)partner[i]|=1u<<j;else partner[j]|=1u<<i;}
   U h=paths(partner);sum+=h;if(h<minh)minh=h;if(h>maxh)maxh=h;
   if(sample.size()<6)sample.push_back({c,h});
  }
  std::cout<<"{\"word\":\""<<(minus?"001":"000")<<"\",\"H\":"<<raw<<",\"Aut\":"<<aut<<",\"presentations\":"<<codes.size()<<",\"sum_partner_H_over_presentations\":"<<sum<<",\"L2_one\":"<<sum*aut<<",\"min_partner_H\":"<<minh<<",\"max_partner_H\":"<<maxh<<",\"sample\":[";
  bool first=true;for(auto[c,h]:sample){if(!first)std::cout<<",";first=false;std::cout<<"["<<c<<","<<h<<"]";}
  std::cout<<"]}"<<(minus?"\n":" ,\n");std::cout.flush();
 }std::cout<<"]\n";
}
