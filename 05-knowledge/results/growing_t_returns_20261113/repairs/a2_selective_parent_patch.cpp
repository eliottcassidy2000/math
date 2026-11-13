#include <algorithm>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iostream>
#include <vector>
#include <utility>
constexpr int N=10,FULL=(1<<N)-1,ROOT=9;
using Rows=std::vector<unsigned>;
static uint32_t dp[1<<N][N],hs[1<<N];
uint32_t fill(const Rows&t,int start=-1){
 std::memset(dp,0,sizeof(dp));for(int u=0;u<N;u++)if(start<0||start==u)dp[1<<u][u]=1;
 for(unsigned mask=1;mask<=FULL;mask++){unsigned occ=mask;while(occ){unsigned bit=occ&-occ;occ-=bit;int u=__builtin_ctz(bit);auto value=dp[mask][u];if(!value)continue;unsigned rem=t[u]&~mask;while(rem){unsigned b=rem&-rem;rem-=b;dp[mask|b][__builtin_ctz(b)]+=value;}}}
 uint32_t h=0;for(int u=0;u<N;u++)h+=dp[FULL][u];return h;
}
bool strong(const Rows&t){for(int u=0;u<N;u++){unsigned seen=1u<<u;while(true){unsigned more=seen;for(int v=0;v<N;v++)if(seen>>v&1)more|=t[v];if(more==seen)break;seen=more;}if(seen!=FULL)return false;}return true;}
int main(){Rows source(N);for(auto&x:source)if(!(std::cin>>x))return 2;std::vector<std::pair<int,int>> arcs;for(int a=0;a<N;a++)for(int b=a+1;b<N;b++)arcs.emplace_back(a,b);
 for(int width=0;width<=4;width++){uint64_t checked=0,one_start=0,hits=0;std::vector<int> patch(width);
  auto process=[&](){checked++;Rows t=source;for(int i:patch){auto[a,b]=arcs[i];t[a]^=1u<<b;t[b]^=1u<<a;}
   Rows reverse=t;for(int u=0;u<N;u++)reverse[u]=FULL^t[u]^(1u<<u);auto h=fill(reverse);std::vector<int> roots;for(int u=0;u<N;u++)if(dp[FULL][u]==1)roots.push_back(u);if(roots.empty())return;one_start++;if(h>846)return;fill(t);hs[0]=1;for(int mask=1;mask<=FULL;mask++){hs[mask]=0;for(int v=0;v<N;v++)hs[mask]+=dp[mask][v];}for(int root:roots){uint64_t F=0;for(int mask=1;mask<=FULL;mask++)F+=uint64_t(dp[mask][root])*hs[FULL^mask];
   if(h+F-1!=1691)continue;hits++;std::cout<<"{\"width\":"<<width<<",\"patch\":[";for(int i=0;i<width;i++){if(i)std::cout<<',';auto[a,b]=arcs[patch[i]];std::cout<<'['<<a<<','<<b<<']';}std::cout<<"],\"state\":["<<h<<",1,"<<hs[FULL^(1<<root)]<<','<<F<<"],\"root\":"<<root<<",\"strong\":"<<(strong(t)?"true":"false")<<"}\n";}
  };
  std::function<void(int,int)> choose=[&](int depth,int lower){if(depth==width){process();return;}for(int i=lower;i<=45-(width-depth);i++){patch[depth]=i;choose(depth+1,i+1);}};choose(0,0);std::cout<<"{\"status\":\"PASS\",\"width\":"<<width<<",\"checked\":"<<checked<<",\"one_start\":"<<one_start<<",\"hits\":"<<hits<<"}\n";
 }
}
