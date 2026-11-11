#include <array>
#include <cstdint>
#include <iostream>
int main(){
  std::uint64_t code,dual; int aut,oldh; std::array<int,9> r;
  std::size_t seen=0,hits=0;
  while(std::cin>>code>>dual>>aut>>oldh){
    for(auto &x:r)std::cin>>x;
    int dp[512][9]{};
    for(int v=0;v<9;++v)dp[1<<v][v]=1;
    for(int mask=1;mask<512;++mask)for(int v=0;v<9;++v){
      if(!dp[mask][v])continue;
      int avail=r[v]&~mask;
      while(avail){int w=__builtin_ctz(avail);avail&=avail-1;dp[mask|(1<<w)][w]+=dp[mask][v];}
    }
    int h=0;for(int v=0;v<9;++v)h+=dp[511][v];
    if(h==221||h==223){std::cout<<code<<' '<<dual<<' '<<h;for(int x:r)std::cout<<' '<<x;std::cout<<'\n';++hits;}
    ++seen;
  }
  std::cerr<<"converse representatives "<<seen<<", target hits "<<hits<<'\n';
  return seen==97144?0:1;
}
