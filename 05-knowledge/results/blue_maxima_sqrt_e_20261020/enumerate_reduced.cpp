#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <vector>
using namespace std;
struct Rep { int a; int orbit; vector<array<int,7>> aut; };
int main(int argc,char**argv){
 int k=argc>1?stoi(argv[1]):6; if(k<1||k>6) return 2;
 int m=k*(k-1)/2,sz=1<<m,n=2*k,N=1<<k;
 int id[7][7]; vector<pair<int,int>> pairs;
 for(int i=0;i<k;i++)for(int j=i+1;j<k;j++){id[i][j]=id[j][i]=pairs.size();pairs.emplace_back(i,j);}
 vector<array<int,7>> perms; array<int,7> p{}; iota(p.begin(),p.begin()+k,0);
 do{perms.push_back(p);}while(next_permutation(p.begin(),p.begin()+k));
 vector<Rep> reps;vector<bool> seen(sz,false);
 for(int a=0;a<sz;a++)if(!seen[a]){
  Rep r{a,0,{}};set<int> orb;
  for(auto p:perms){int t=0;for(int z=0;z<m;z++){auto [i,j]=pairs[z];int u=p[i],v=p[j],bit=(a>>z)&1;if(u>v)bit^=1;t|=bit<<id[u][v];} orb.insert(t);if(t==a)r.aut.push_back(p);}
  for(int x:orb)seen[x]=true;r.orbit=orb.size();if(r.orbit*r.aut.size()!=perms.size())return 3;reps.push_back(r);
 }
 int avail[128]={};for(int mask=1;mask<N;mask++)for(int j=0;j<k;j++)if(!(mask>>j&1))avail[mask]|=3<<(2*j);
 vector<uint32_t> all;all.reserve(reps.size()*sz);vector<uint64_t> hist(100000,0);int best=0;vector<pair<int,int>> maxima;uint64_t sum=0,sq=0,count=0;
 auto start=chrono::steady_clock::now();
 for(auto &r:reps){
  int rb=0; vector<int> rmax;
  for(int b=0;b<sz;b++){
   int adj[14]={};auto edge=[&](int u,int v,int yes){if(yes)adj[u]|=1<<v;else adj[v]|=1<<u;};
   for(int z=0;z<m;z++){auto [i,j]=pairs[z];int a=r.a>>z&1,c=b>>z&1;edge(2*i,2*j,a);edge(2*j+1,2*i+1,a);edge(2*i,2*j+1,c);edge(2*j,2*i+1,c);}
   uint32_t dp[128][14]={};for(int i=0;i<n;i++)dp[1<<(i/2)][i]=1;
   for(int mask=1;mask<N;mask++)for(int i=0;i<n;i++){int v=dp[mask][i];if(v){int next=adj[i]&avail[mask];while(next){int j=__builtin_ctz(next);next&=next-1;dp[mask|(1<<(j/2))][j]+=v;}}}
   int h=0;for(int i=0;i<n;i+=2)h+=dp[N-1][i];if(!(h&1))return 4;all.push_back(h);
   hist[h]+=r.orbit;sum+=(uint64_t)h*r.orbit;sq+=(uint64_t)h*h*r.orbit;count+=r.orbit;
   if(h>best){best=h;maxima.clear();}if(h==best)maxima.emplace_back(r.a,b);
   if(h>rb){rb=h;rmax.clear();}if(h==rb)rmax.push_back(b);
  }
  cerr<<"A="<<r.a<<" orbit="<<r.orbit<<" maximum="<<rb<<" b_count="<<rmax.size()<<" elapsed="<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<endl;
 }
 ofstream bin("reduced_k"+to_string(k)+".bin",ios::binary);bin.write((char*)all.data(),all.size()*sizeof(uint32_t));
 cout<<"{\"k\":"<<k<<",\"tournament_representatives\":[";bool first=true;
 for(auto&r:reps){if(!first)cout<<",";first=false;cout<<"{\"A\":"<<r.a<<",\"orbit_size\":"<<r.orbit<<",\"aut_order\":"<<r.aut.size()<<"}";}cout<<"],\"reduced_frames_checked\":"<<all.size()<<",\"full_labelled_frames\":"<<count<<",\"sum_h\":"<<sum<<",\"sum_h_squared\":"<<sq<<",\"maximum\":"<<best<<",\"maximizing_reduced_frames\":[";
 first=true;for(auto[a,b]:maxima){if(!first)cout<<",";first=false;cout<<"["<<a<<","<<b<<"]";}cout<<"],\"histogram\":{";first=true;for(int h=0;h<(int)hist.size();h++)if(hist[h]){if(!first)cout<<",";first=false;cout<<"\""<<h<<"\":"<<hist[h];}cout<<"},\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
}
