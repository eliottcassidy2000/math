#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>
using I=std::int64_t;
int n;std::vector<unsigned> adjacency;std::vector<I> memo;
I paths_from(unsigned set,int v){
 I &cached=memo[std::size_t(set)*n+v];if(cached>=0)return cached;
 unsigned rem=set^(1u<<v);if(!rem)return cached=1;
 cached=0;unsigned next=adjacency[v]&rem;
 while(next){int w=__builtin_ctz(next);next&=next-1;cached+=paths_from(rem,w);}
 return cached;
}
I paths(){memo.assign(std::size_t(1u<<n)*n,-1);I h=0;for(int v=0;v<n;v++)h+=paths_from((1u<<n)-1,v);return h;}
int rank(std::vector<std::vector<int>> a,int p){
 int r=0;if(a.empty())return 0;
 for(auto &row:a)for(int &x:row)x%=p;
 for(int j=0;j<int(a[0].size())&&r<int(a.size());j++){
  int pivot=r;while(pivot<int(a.size())&&a[pivot][j]==0)pivot++;
  if(pivot==int(a.size()))continue;std::swap(a[pivot],a[r]);
  int inv=1;while((inv*a[r][j])%p!=1)inv++;
  for(int k=j;k<int(a[r].size());k++)a[r][k]=(a[r][k]*inv)%p;
  for(int i=r+1;i<int(a.size());i++){int t=a[i][j];if(!t)continue;
   for(int k=j;k<int(a[i].size());k++)a[i][k]=(a[i][k]-t*a[r][k]%p+p)%p;
  }r++;
 }return r;
}
int main(){
 int old,q;bool first=true;std::uint64_t total=0;std::cout<<"{\"audits\":[";
 while(std::cin>>old>>q){
  if(old<3||old>5||q<1||q>2000)return 2;
  std::vector<std::vector<int>> words(q);
  for(auto &w:words){int d;if(!(std::cin>>d)||d<0||d>4)return 3;w.resize(d);
   for(int &op:w)if(!(std::cin>>op)||op<0||op>old)return 4;
  }
  int m=old*(old-1)/2,states=1<<m;
  std::vector<std::vector<int>> values(q,std::vector<int>(states));
  for(int t=0;t<q;t++)for(int mask=0;mask<states;mask++){
   adjacency.assign(old+words[t].size(),0);int e=0;
   for(int i=0;i<old;i++)for(int j=i+1;j<old;j++,e++)
    if(mask&(1<<e))adjacency[i]|=1u<<j;else adjacency[j]|=1u<<i;
   n=old;
   for(int op:words[t]){for(int i=0;i<n;i++)
     if(op<old&&i==op)adjacency[n]|=1u<<i;else adjacency[i]|=1u<<n;
    n++;
   }
   I h=paths();if(h<1||h%2!=1||h>362880)return 5;
   values[t][mask]=int(h);total++;
  }
  if(!first)std::cout<<",";first=false;
  std::cout<<"{\"n\":"<<old<<",\"contexts\":"<<q<<",\"orientations\":"<<states
   <<",\"rank_mod101\":"<<rank(values,101)<<",\"rank_mod103\":"<<rank(values,103)<<"}";
 }
 std::cout<<"],\"literal_graph_counts\":"<<total<<",\"algorithm\":\"independent starting-endpoint DP; elimination of full Boolean evaluation rows\"}\n";
}
