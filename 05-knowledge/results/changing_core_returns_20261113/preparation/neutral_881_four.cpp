#include <array>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <functional>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>
constexpr int N=11,FULL=(1<<N)-1;
using Rows=std::array<unsigned,N>;
static uint32_t dp[1<<N][N],hs[1<<N];
uint32_t fill(const Rows&t,int start=-1){
 std::memset(dp,0,sizeof(dp));
 for(int u=0;u<N;u++)if(start<0||start==u)dp[1<<u][u]=1;
 for(unsigned mask=1;mask<=FULL;mask++){
  unsigned occupied=mask;
  while(occupied){unsigned bit=occupied&-occupied;occupied-=bit;int u=__builtin_ctz(bit);unsigned value=dp[mask][u];if(!value)continue;
   unsigned rem=t[u]&~mask;
   while(rem){unsigned b=rem&-rem;rem-=b;dp[mask|b][__builtin_ctz(b)]+=value;}
  }
 }
 uint32_t h=0;for(int u=0;u<N;u++)h+=dp[FULL][u];return h;
}
bool strong(const Rows&t){
 for(int u=0;u<N;u++){unsigned seen=1u<<u;while(true){unsigned more=seen;for(int v=0;v<N;v++)if(seen>>v&1)more|=t[v];if(more==seen)break;seen=more;}if(seen!=FULL)return false;}return true;
}
int main(int argc,char**argv){
 int width=argc>1?std::atoi(argv[1]):4;if(width<3||width>4)return 2;
 Rows source;for(auto &row:source)if(!(std::cin>>row))return 3;
 std::vector<std::pair<int,int>> arcs;for(int a=0;a<N;a++)for(int b=a+1;b<N;b++)arcs.emplace_back(a,b);
 uint64_t considered=0,filtered=0,neutral=0;const int root=10;std::vector<int> patch(width);
 auto process=[&](){
  considered++;Rows t=source;for(int i:patch){auto[u,v]=arcs[i];t[u]^=1u<<v;t[v]^=1u<<u;}
  int indegree=0;for(int v=0;v<N;v++)indegree+=(t[v]>>root)&1u;if(indegree<2)return;filtered++;
  auto h=fill(t);if(h!=881)return;neutral++;
  hs[0]=1;for(int mask=1;mask<=FULL;mask++){hs[mask]=0;for(int v=0;v<N;v++)hs[mask]+=dp[mask][v];}
  uint64_t F=0;for(int mask=1;mask<=FULL;mask++)F+=uint64_t(dp[mask][root])*hs[FULL^mask];uint32_t q=hs[FULL^(1<<root)],s=fill(t,root);
  std::cout<<"{\"arcs\":[";for(int i=0;i<width;i++){if(i)std::cout<<',';auto[u,v]=arcs[patch[i]];std::cout<<'['<<u<<','<<v<<']';}
  std::cout<<"],\"state\":[881,"<<s<<','<<q<<','<<F<<"],\"indegree\":"<<indegree<<",\"strong\":"<<(strong(t)?"true":"false")<<"}\n";
 };
 std::function<void(int,int)> choose=[&](int depth,int lower){
  if(depth==width){process();return;}
  for(int i=lower;i<=55-(width-depth);i++){patch[depth]=i;choose(depth+1,i+1);}
 };
 choose(0,0);
 std::cout<<"{\"status\":\"PASS\",\"width\":"<<width<<",\"all_arc_patches\":"<<considered<<",\"indegree_at_least_two\":"<<filtered<<",\"H881\":"<<neutral<<"}\n";
}
