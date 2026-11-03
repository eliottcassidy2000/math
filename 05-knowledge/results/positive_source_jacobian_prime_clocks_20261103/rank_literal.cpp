// Independent lower bound: actual tournaments, no forest/contraction code.
#include <cstdint>
#include <iostream>
#include <set>
#include <vector>
using I=std::int64_t;
I count(const std::vector<unsigned>&a){
 int n=int(a.size());std::vector<I>dp(std::size_t(1u<<n)*n,0);
 for(int v=0;v<n;v++)dp[(1u<<v)*n+v]=1;
 for(unsigned m=1;m<(1u<<n);m++)for(int v=0;v<n;v++){
  I c=dp[m*n+v];if(!c)continue;unsigned next=a[v]&~m;
  while(next){int w=__builtin_ctz(next);next&=next-1;dp[(m|(1u<<w))*n+w]+=c;}
 }
 I h=0;for(int v=0;v<n;v++)h+=dp[((1u<<n)-1)*n+v];return h;
}
int rank(std::vector<std::vector<int>>a,int p){
 int r=0;for(auto&row:a)for(int&v:row)v%=p;
 for(int j=0;j<int(a[0].size())&&r<int(a.size());j++){
  int pivot=r;while(pivot<int(a.size())&&!a[pivot][j])pivot++;
  if(pivot==int(a.size()))continue;std::swap(a[r],a[pivot]);
  int inv=1;while(inv*a[r][j]%p!=1)inv++;
  for(int k=j;k<int(a[r].size());k++)a[r][k]=a[r][k]*inv%p;
  for(int i=r+1;i<int(a.size());i++)if(a[i][j]){
   int t=a[i][j];for(int k=j;k<int(a[i].size());k++)a[i][k]=(a[i][k]-t*a[r][k]%p+p)%p;
  }r++;
 }return r;
}
int main(){
 int q;if(!(std::cin>>q)||q!=571)return 2;std::vector<std::vector<int>> words(q);
 for(auto&w:words){int d;std::cin>>d;if(d<0||d>5)return 3;w.resize(d);
  for(int&z:w)if(!(std::cin>>z)||z<0||z>6)return 4;}
 std::vector<unsigned>sources;std::set<unsigned>seen;std::uint32_t x=137;
 while(sources.size()<768){x^=x<<13;x^=x>>17;x^=x<<5;unsigned z=x&32767;
  if(seen.insert(z).second)sources.push_back(z);}
 std::vector<std::vector<int>>values(q,std::vector<int>(sources.size()));
 for(int j=0;j<q;j++){
  for(std::size_t s=0;s<sources.size();s++){
   unsigned code=sources[s];std::vector<unsigned>a(6,0);int e=0;
   for(int u=0;u<6;u++)for(int v=u+1;v<6;v++,e++)
    if(code&(1u<<e))a[u]|=1u<<v;else a[v]|=1u<<u;
   for(int op:words[j]){int n=int(a.size());a.push_back(op<6?(1u<<op):0);
    for(int v=0;v<n;v++)if(v!=op||op==6)a[v]|=1u<<n;}
   I h=count(a);if(h<1||!(h&1)||h>39916800)return 5;values[j][s]=int(h);
  }
  if(j%100==99)std::cerr<<"literal rows "<<j+1<<"/"<<q<<"\n";
 }
 int r101=rank(values,101),r103=rank(values,103);
 std::cout<<"{\"actual_contexts\":"<<q<<",\"actual_sources\":"<<sources.size()
          <<",\"literal_graph_counts\":"<<q*sources.size()<<",\"rank_mod101\":"<<r101
          <<",\"rank_mod103\":"<<r103<<",\"max_births\":5,\"source_codes\":[";
 for(std::size_t i=0;i<sources.size();i++){if(i)std::cout<<",";std::cout<<sources[i];}
 std::cout<<"],\"algorithm\":\"Boolean subset DP then full evaluation-matrix elimination; no forest formulas\"}\n";
 return r101==571&&r103==571?0:6;
}
