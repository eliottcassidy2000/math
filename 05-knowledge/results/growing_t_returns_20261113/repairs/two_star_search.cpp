#include <array>
#include <vector>
#include <iostream>
#include <cstring>
#include <cstdint>
#include <algorithm>
#include <string>
#ifndef ORDER
#define ORDER 11
#endif
#ifndef TARGET
#define TARGET 3311
#endif
#ifndef REQUIRED_H
#define REQUIRED_H 0
#endif
constexpr int SMALL=ORDER-1,SMALLFULL=(1<<SMALL)-1;
template<int N> struct Counter{
 uint32_t dp[1<<N][N],h[1<<N];
 uint32_t fill(const std::array<unsigned,N>&t){std::memset(dp,0,sizeof(dp));for(int v=0;v<N;v++)dp[1<<v][v]=1;
 for(unsigned mask=1;mask<(1<<N);mask++){unsigned bits=mask;while(bits){unsigned b=bits&-bits;bits-=b;int u=__builtin_ctz(b);auto value=dp[mask][u];if(!value)continue;unsigned rem=t[u]&~mask;while(rem){unsigned c=rem&-rem;rem-=c;dp[mask|c][__builtin_ctz(c)]+=value;}}h[mask]=0;for(int v=0;v<N;v++)h[mask]+=dp[mask][v];}h[0]=1;return h[(1<<N)-1];}
};
static Counter<SMALL> small;static Counter<ORDER> big;
int main(){std::array<unsigned,ORDER>T;for(auto&x:T)if(!(std::cin>>x))return 2;uint64_t choices=0,inner_stars=0;bool found=false;
 for(int root=0;root<ORDER&&!found;root++){
  std::vector<int> old;for(int v=0;v<ORDER;v++)if(v!=root)old.push_back(v);std::array<unsigned,SMALL> base;
  for(int i=0;i<SMALL;i++){base[i]=0;for(int j=0;j<SMALL;j++)if(T[old[i]]>>old[j]&1)base[i]|=1u<<j;}
  for(int w=0;w<SMALL&&!found;w++){
   std::vector<int> other;for(int i=0;i<SMALL;i++)if(i!=w)other.push_back(i);
   for(unsigned orient=0;orient<(1u<<(SMALL-1))&&!found;orient++){
    inner_stars++;auto U=base;U[w]=0;for(int i=0;i<(SMALL-1);i++){int v=other[i];if(orient>>i&1){U[w]|=1u<<v;U[v]&=~(1u<<w);}else U[v]|=1u<<w;}
    auto rev=U;for(int i=0;i<SMALL;i++)rev[i]=unsigned(SMALLFULL)^U[i]^(1u<<i);small.fill(rev);std::vector<int> ones,zeros;
    for(int i=0;i<SMALL;i++){auto s=small.dp[SMALLFULL][i];if(s==1)ones.push_back(i);if(s==0)zeros.push_back(i);}if(ones.empty())continue;
    auto hU=small.fill(U);if(hU>((TARGET+1)/2))continue;std::array<uint64_t,SMALL> F{};for(int mask=1;mask<(1u<<SMALL);mask++)for(int i=0;i<SMALL;i++)F[i]+=uint64_t(small.dp[mask][i])*small.h[SMALLFULL^mask];
    for(int entry:ones)for(unsigned z=0;z<(1u<<zeros.size())&&!found;z++){
     unsigned O=1u<<entry;for(unsigned i=0;i<zeros.size();i++)if(z>>i&1)O|=1u<<zeros[i];uint64_t FF=hU;for(int i=0;i<SMALL;i++)if(!(O>>i&1))FF+=F[i];if(FF>TARGET)continue;
     std::array<unsigned,ORDER>V{};for(int i=0;i<SMALL;i++){for(int j=0;j<SMALL;j++)if(U[i]>>j&1)V[old[i]]|=1u<<old[j];if(O>>i&1)V[root]|=1u<<old[i];else V[old[i]]|=1u<<root;}
     choices++;auto hV=big.fill(V);if(REQUIRED_H && hV!=REQUIRED_H)continue;if(hV+FF-1!=TARGET)continue;found=true;
     std::cout<<"{\"status\":\"HIT\",\"root\":"<<root<<",\"other_star\":"<<old[w]<<",\"parent_state\":["<<hV<<",1,"<<hU<<','<<FF<<"],\"rows\":[";for(int i=0;i<ORDER;i++){if(i)std::cout<<',';std::cout<<V[i];}std::cout<<"],\"changed_arcs\":[";bool first=true;for(int i=0;i<ORDER;i++)for(int j=i+1;j<ORDER;j++)if((T[i]^V[i])>>j&1){if(!first)std::cout<<',';first=false;std::cout<<'['<<i<<','<<j<<']';}std::cout<<"]}\n";
    }
   }
  }
  std::cerr<<"completed root "<<root<<", choices "<<choices<<"\n";
 }
 std::cout<<"{\"status\":\"PASS\",\"stopped_at_first_hit\":"<<(found?"true":"false")<<",\"inner_stars_checked\":"<<inner_stars<<",\"target_face_choices_counted\":"<<choices<<"}\n";
}
