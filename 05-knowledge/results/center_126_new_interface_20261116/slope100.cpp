#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using A=std::array<unsigned,10>;
using S=std::array<int,1024>;
S subsets(const A& a,int n,int root=-1) {
 int d[1024][10]={}; S sum{};
 for(int i=0;i<n;++i) if(root<0 || root==i)d[1<<i][i]=1;
 for(int s=1;s<(1<<n);++s)for(int i=0;i<n;++i)if(int v=d[s][i]){
  sum[s]+=v;
  for(unsigned E=a[i]&~s;E;E&=E-1){unsigned b=E&-E;d[s|b][__builtin_ctz(b)]+=v;}
 }
 sum[0]=root<0; return sum;
}
A append(A a,int n,unsigned v){
 a[n]=v;
 for(int i=0;i<n;++i) if(!(v>>i&1))a[i]|=1u<<n;
 return a;
}


int main(int argc,char**argv){
 if(argc!=2)return 1;
 std::ifstream f(argv[1]); std::string s; std::getline(f,s);
 unsigned bases=0,gated=0,wins=0;
 while(std::getline(f,s)){
  std::istringstream ss(s);unsigned long long code;unsigned h;char c;
  if(!(ss>>code>>c>>h)||c!=',')return 2;
  if(h>127)continue;
  A a{};for(int j=0;j<8;++j)a[j]=(code>>(8*j))&255;
  for(unsigned x=0;x<256;++x){
    auto b=append(a,8,x);
    auto h=subsets(b,9);
    if(h[511]<27||h[511]>127)continue;
    ++bases;
    // Ending paths of the converse are starts in the original base.
    int dp[512][9]={};unsigned R[9]={},D[9]={};
    for(int u=0;u<9;++u)dp[1<<u][u]=1;
    for(int v=1;v<512;++v)for(int u=0;u<9;++u)if(int w=dp[v][u]){
     ((__builtin_popcount(v)&1)?R:D)[u]+=w*h[511^v];
     for(unsigned z=(511^b[u]^(1u<<u))&~v;z;z&=z-1){
      unsigned bit=z&-z;dp[v|bit][__builtin_ctz(bit)]+=w;
     }
    }
    unsigned r[512]={},d[512]={};
    for(unsigned y=1;y<512;++y){
      unsigned bit=y&-y,u=__builtin_ctz(bit);
      r[y]=r[y^bit]+R[u];d[y]=d[y^bit]+D[u];
      if(r[y]!=50||d[y]!=50)continue;
      ++gated;
      auto A10=append(b,9,y);
      if(subsets(A10,10)[1023]!=127)continue;
      ++wins;std::cout<<"hits g="<<h[511]<<" root9 rows ";
      for(auto w:A10)std::cout<<w<<" ";
      std::cout<<std::endl;
    }
  }
 }
 std::cout<<bases<<" base presentations, "<<gated<<" masks with slope 100 and stable, "<<wins<<" 127 hits\n";
 return (bases==17090 && gated==84 && wins==0)?0:3;
}
