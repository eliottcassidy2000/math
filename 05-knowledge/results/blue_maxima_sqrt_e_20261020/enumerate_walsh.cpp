#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <vector>
using namespace std;
struct Rep {int a,orbit,aut;};
string decimal(unsigned __int128 n){if(!n)return "0";string s;while(n){s+=char('0'+n%10);n/=10;}reverse(s.begin(),s.end());return s;}
int countdp(int k,int a,int b){int n=2*k,N=1<<k,adj[14]={},z=0;auto edge=[&](int u,int v,int x){if(x)adj[u]|=1<<v;else adj[v]|=1<<u;};for(int i=0;i<k;i++)for(int j=i+1;j<k;j++,z++){int x=a>>z&1,y=b>>z&1;edge(2*i,2*j,x);edge(2*j+1,2*i+1,x);edge(2*i,2*j+1,y);edge(2*j,2*i+1,y);}int dp[128][14]={};for(int i=0;i<n;i++)dp[1<<(i/2)][i]=1;for(int mask=1;mask<N;mask++){int avail=0;for(int j=0;j<k;j++)if(!(mask>>j&1))avail|=3<<(2*j);for(int i=0;i<n;i++)if(dp[mask][i]){int next=adj[i]&avail;while(next){int j=__builtin_ctz(next);next&=next-1;dp[mask|1<<(j/2)][j]+=dp[mask][i];}}}int h=0;for(int i=0;i<n;i+=2)h+=dp[N-1][i];return h;}
int main(int argc,char**argv){int k=argc>1?stoi(argv[1]):7;if(k<1||k>7)return 2;int m=k*(k-1)/2,sz=1<<m,scale=1<<(k-1);int id[7][7];vector<pair<int,int>> pairs;for(int i=0;i<k;i++)for(int j=i+1;j<k;j++){id[i][j]=id[j][i]=pairs.size();pairs.emplace_back(i,j);}vector<array<int,7>> perms;array<int,7> p{};iota(p.begin(),p.begin()+k,0);do{perms.push_back(p);}while(next_permutation(p.begin(),p.begin()+k));vector<Rep>reps;vector<bool>seen(sz,false);for(int a=0;a<sz;a++)if(!seen[a]){set<int>orb;int aut=0;for(auto p:perms){int t=0;for(int z=0;z<m;z++){auto[i,j]=pairs[z];int u=p[i],v=p[j],bit=(a>>z)&1;if(u>v)bit^=1;t|=bit<<id[u][v];}orb.insert(t);if(t==a)aut++;}for(int t:orb)seen[t]=true;if(orb.size()*aut!=perms.size())return 3;reps.push_back({a,(int)orb.size(),aut});}
 vector<uint32_t>expected;if(argc>2 && string(argv[2])!="-"){ifstream f(argv[2],ios::binary);expected.resize(sz*reps.size());f.read((char*)expected.data(),expected.size()*4);if(f.gcount()!=expected.size()*4)return 4;}
 ofstream stream;if(argc>3)stream.open(argv[3],ios::binary);uint64_t total=0,sum=0,updates=0,sampled=0;unsigned __int128 squares=0;vector<uint64_t>hist(10000000,0);int best=0;vector<pair<int,int>> maxima;mt19937 rng(20261020);auto start=chrono::steady_clock::now();
 cout<<"{\"k\":"<<k<<",\"method\":\"signed-prefix Walsh transform\",\"tournament_representatives\":[";bool first=true;
 for(size_t r=0;r<reps.size();r++){auto R=reps[r];vector<int32_t>coeff(sz,0);
  for(auto p:perms)for(int signs=0;signs<scale;signs++){
   bool good=true;int sel=0,req=0;
   for(int t=0;t<k-1;t++){int i=p[t],j=p[t+1],s=signs>>t&1,u=signs>>(t+1)&1,q=id[i][j];if(s==u){int arc=(R.a>>q&1)^(i>j);if(arc==s){good=false;break;}}else{sel|=1<<q;if(!s)req|=1<<q;}}
   if(!good)continue;int weight=scale>>__builtin_popcount((unsigned)sel),sub=sel;
   while(true){coeff[sub]+=(__builtin_parity((unsigned)(req&sub))?-weight:weight);updates++;if(!sub)break;sub=(sub-1)&sel;}
  }
  for(int len=1;len<sz;len*=2)for(int base=0;base<sz;base+=2*len)for(int j=0;j<len;j++){int x=coeff[base+j],y=coeff[base+len+j];coeff[base+j]=x+y;coeff[base+len+j]=x-y;}
  int rb=0,rc=0;uint64_t rsum=0,rsq=0;vector<int>rmax;
  for(int b=0;b<sz;b++){int v=coeff[b];if(v%scale){cerr<<"nonintegral"<<endl;return 5;}int h=v/scale;if(h<1||!(h&1)||(h>=(int)hist.size())){cerr<<"bad h"<<endl;return 6;}if(!expected.empty()&&expected[r*sz+b]!=(uint32_t)h){cerr<<"mismatch"<<endl;return 7;}coeff[b]=h;hist[h]+=R.orbit;sum+=(uint64_t)h*R.orbit;squares+=(unsigned __int128)h*h*R.orbit;total+=R.orbit;rsum+=h;rsq+=(uint64_t)h*h;if(h>best){best=h;maxima.clear();}if(h==best)maxima.emplace_back(R.a,b);if(h>rb){rb=h;rc=0;rmax.clear();}if(h==rb){rc++;rmax.push_back(b);}}
  if(stream.is_open())stream.write((char*)coeff.data(),coeff.size()*sizeof(int32_t));for(int t=0;t<256;t++){int b=rng()&(sz-1);if(countdp(k,R.a,b)!=coeff[b]){cerr<<"sample mismatch"<<endl;return 8;}sampled++;}
  for(int b:rmax){if(countdp(k,R.a,b)!=coeff[b]){cerr<<"maximum mismatch"<<endl;return 9;}sampled++;}
  if(!first)cout<<",";first=false;cout<<"{\"A\":"<<R.a<<",\"orbit_size\":"<<R.orbit<<",\"aut_order\":"<<R.aut<<",\"max_h\":"<<rb<<",\"max_count\":"<<rc<<",\"sum_h\":"<<rsum<<",\"sum_h_squared\":"<<rsq<<"}";
  cerr<<r+1<<"/"<<reps.size()<<" A="<<R.a<<" max="<<rb<<" global="<<best<<" elapsed="<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<endl;
 }
 cout<<"],\"reduced_frames_checked\":"<<(uint64_t)reps.size()*sz<<",\"full_labelled_frames\":"<<total<<",\"sum_h\":"<<sum<<",\"sum_h_squared\":"<<decimal(squares)<<",\"maximum\":"<<best<<",\"maximizing_reduced_frames\":[";first=true;for(auto[a,b]:maxima){if(!first)cout<<",";first=false;cout<<"["<<a<<","<<b<<"]";}cout<<"],\"histogram\":{";first=true;for(int h=0;h<(int)hist.size();h++)if(hist[h]){if(!first)cout<<",";first=false;cout<<"\""<<h<<"\":"<<hist[h];}cout<<"},\"coefficient_updates\":"<<updates<<",\"sampled_DP_checks\":"<<sampled<<",\"pointwise_DP_checks\":"<<expected.size()<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
}
