// Complete exact Hamming-distance-two arc-flip shells. No external libraries.
#include <cstdint>
#include <iostream>
#include <vector>
#include <utility>
#include <cstdlib>
#include <array>
#include <set>
#include <algorithm>
using namespace std;
uint64_t paths(const vector<unsigned>& adj,vector<uint64_t>& dp){
 int n=adj.size();unsigned N=1u<<n;fill(dp.begin(),dp.end(),0);
 for(int i=0;i<n;i++)dp[(1u<<i)*n+i]=1;
 for(unsigned S=1;S<N;S++){
  unsigned vs=S;while(vs){int v=__builtin_ctz(vs);vs&=vs-1;uint64_t x=dp[S*n+v];if(!x)continue;
   unsigned choices=adj[v]&~S;while(choices){int w=__builtin_ctz(choices);choices&=choices-1;dp[(S|(1u<<w))*n+w]+=x;}
  }
 }
 uint64_t ans=0;for(int v=0;v<n;v++)ans+=dp[(N-1)*n+v];return ans;
}
void flip(vector<unsigned>& adj,pair<int,int> e){adj[e.first]^=1u<<e.second;adj[e.second]^=1u<<e.first;}
int main(int argc,char**argv){int lo=argc>1?atoi(argv[1]):2,hi=argc>2?atoi(argv[2]):5;int mode=argc>3?atoi(argv[3]):2;cout<<"[";bool first=true;
 for(int a=lo;a<=hi;a++){
  int n=2*a+2;vector<unsigned> adj(n);vector<pair<int,int>> edges;
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){edges.push_back({i,j});if(j==n-1&&i<a)adj[j]|=1u<<i;else adj[i]|=1u<<j;}
  vector<uint64_t> dp((1u<<n)*n);uint64_t base=paths(adj,dp),target=base-2;long long tested=0,hits=0;uint64_t minloss=base;vector<pair<int,int>> hitidx;
  if(mode==3){
   set<array<int,3>> triples;vector<int> apex,cross;
   for(int e=0;e<(int)edges.size();e++){auto z=edges[e];if(z.second==n-1)apex.push_back(e);else if(z.first<a&&z.second>=a)cross.push_back(e);}
   auto put=[&](int e,int f,int g){array<int,3>x={e,f,g};sort(x.begin(),x.end());triples.insert(x);};
   for(int i=0;i<(int)apex.size();i++)for(int j=i+1;j<(int)apex.size();j++){
    for(int e:cross)put(apex[i],apex[j],e);
    for(int k=j+1;k<(int)apex.size();k++)put(apex[i],apex[j],apex[k]);
   }
   for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)for(int k=j+1;k<n;k++){
    int es[3],p=0;for(int e=0;e<(int)edges.size();e++){auto z=edges[e];if((z.first==i&&z.second==j)||(z.first==i&&z.second==k)||(z.first==j&&z.second==k))es[p++]=e;}put(es[0],es[1],es[2]);
   }
   vector<array<int,3>> good,mins;
   for(auto z:triples){for(int e:z)flip(adj,edges[e]);auto h=paths(adj,dp);if(h==target)good.push_back(z);if(h<base&&base-h<=minloss){if(base-h<minloss)mins.clear();minloss=base-h;mins.push_back(z);}for(int e:z)flip(adj,edges[e]);}
   if(!first)cout<<",";first=false;
   cout<<"{\"a\":"<<a<<",\"order\":"<<n<<",\"base\":"<<base<<",\"target\":"<<target<<",\"tested\":"<<triples.size()<<",\"hits\":"<<good.size()<<",\"minimum_positive_decrement\":"<<minloss<<",\"edge_triples\":[";
   bool one=true;for(auto z:good){if(!one)cout<<",";one=false;cout<<"[";for(int i=0;i<3;i++){if(i)cout<<",";auto e=edges[z[i]];cout<<"["<<e.first<<","<<e.second<<"]";}cout<<"]";}cout<<"],\"minimum_triples\":[";
   one=true;for(auto z:mins){if(!one)cout<<",";one=false;cout<<"[";for(int i=0;i<3;i++){if(i)cout<<",";auto e=edges[z[i]];cout<<"["<<e.first<<","<<e.second<<"]";}cout<<"]";}
   set<array<int,4>> repairs;long long fourthtests=0;
   if(minloss>2){for(auto z:mins){for(int e:z)flip(adj,edges[e]);for(int f=0;f<(int)edges.size();f++){
    if(find(z.begin(),z.end(),f)!=z.end())continue;flip(adj,edges[f]);fourthtests++;if(paths(adj,dp)==target){array<int,4>rr={z[0],z[1],z[2],f};sort(rr.begin(),rr.end());repairs.insert(rr);}flip(adj,edges[f]);
   }for(int e:z)flip(adj,edges[e]);}}
   cout<<"],\"fourth_edge_tests\":"<<fourthtests<<",\"four_edge_repairs\":[";
   one=true;for(auto z:repairs){if(!one)cout<<",";one=false;cout<<"[";for(int i=0;i<4;i++){if(i)cout<<",";auto e=edges[z[i]];cout<<"["<<e.first<<","<<e.second<<"]";}cout<<"]";}cout<<"]}"<<flush;continue;
  }
  for(int e=0;e<(int)edges.size();e++){flip(adj,edges[e]);for(int f=e+1;f<(int)edges.size();f++){
   flip(adj,edges[f]);uint64_t h=paths(adj,dp);tested++;if(h==target){hits++;hitidx.push_back({e,f});}if(h<base&&base-h<minloss)minloss=base-h;flip(adj,edges[f]);
  }flip(adj,edges[e]);}
  if(!first)cout<<",";first=false;
  cout<<"{\"a\":"<<a<<",\"order\":"<<n<<",\"base\":"<<base<<",\"target\":"<<target<<",\"tested\":"<<tested<<",\"hits\":"<<hits<<",\"minimum_positive_decrement\":"<<minloss<<",\"edge_pairs\":[";
  bool one=true;for(auto z:hitidx){if(!one)cout<<",";one=false;auto e=edges[z.first],f=edges[z.second];cout<<"[["<<e.first<<","<<e.second<<"],["<<f.first<<","<<f.second<<"]]";}
  cout<<"]}"<<flush;
 }
 cout<<"]\n";
}
