#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>
#include <cstdint>
int main(int argc,char**){
 int n,sr,si,tr,ti;std::cin>>n>>sr>>si>>tr>>ti;
 std::vector<unsigned>A(n),B(n);for(auto &x:A)std::cin>>x;for(auto &x:B)std::cin>>x;
 std::vector<int> p(n,-1),bestp;std::vector<int> assigned{sr,si};p[sr]=tr;p[si]=ti;
 int best=n*n;unsigned used=(1<<tr)|(1<<ti);uint64_t nodes=0;
 auto edge=[&](int a,int b,int x,int y){return int(((A[a]>>b)^(B[x]>>y))&1);};
 std::function<void(int,unsigned)> dfs=[&](int cost,unsigned mask){
  nodes++;if(cost>=best)return;
  if((int)assigned.size()==n){best=cost;bestp=p;return;}
  std::vector<int> S,J;unsigned sm=0;
  for(int i=0;i<n;i++){if(p[i]<0){S.push_back(i);sm|=1<<i;}if(!(mask>>i&1))J.push_back(i);}
  int lb=0,pick=-1,contrast=-1;unsigned jm=((1<<n)-1)^mask;
  for(int i:S){
   int mn=100000;int da=__builtin_popcount(A[i]&sm);
   int score=std::abs(2*da-int(S.size())+1);
   if(score>contrast){contrast=score;pick=i;}
   for(int j:J){int z=0;for(int a:assigned)z+=edge(i,a,j,p[a]);int db=__builtin_popcount(B[j]&jm);mn=std::min(mn,2*z+std::abs(da-db));}
   lb+=mn;
  }
  if(argc==1 && 2*cost+lb>=2*best)return;
  std::vector<std::pair<int,int>> choices;
  for(int j:J){int z=0;for(int a:assigned)z+=edge(pick,a,j,p[a]);choices.emplace_back(z,j);}
  std::sort(choices.begin(),choices.end());
  for(auto [z,j]:choices){p[pick]=j;assigned.push_back(pick);dfs(cost+z,mask|(1<<j));assigned.pop_back();p[pick]=-1;}
 };
 dfs(edge(sr,si,tr,ti),used);
 std::cout<<best<<" "<<nodes<<"\n";for(int x:bestp)std::cout<<x<<" ";std::cout<<"\n";
}
