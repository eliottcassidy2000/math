#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
using namespace std;
using U=uint64_t;
int n,goal; U adj[64],full; vector<int> path; map<U,pair<vector<int>,uint64_t>> witnesses; uint64_t calls=0,countpaths=0;
int bit(U x){return __builtin_ctzll(x);}
void dfs(int x,U unused,U odd){
 ++calls;
 if(!unused){if(x==goal){++countpaths;auto &entry=witnesses[odd]; if(entry.second==0)entry.first=path; ++entry.second;}return;}
 if(x==goal)return;
 // Every unused internal vertex needs two available neighbors; the terminal needs one.
 U available=unused|(1ULL<<x),scan=unused;
 while(scan){int v=bit(scan);scan&=scan-1;if(__builtin_popcountll(adj[v]&available)<(v==goal?1:2))return;}
 // The remaining graph, including the current endpoint, must be connected.
 U seen=1ULL<<x,front=seen;
 while(front){int v=bit(front);front&=front-1;U fresh=adj[v]&available&~seen;seen|=fresh;front|=fresh;}
 if(seen!=available)return;
 U opts=adj[x]&unused;
 while(opts){int v=bit(opts);opts&=opts-1;
  int pos=path.size()+1;
  // Any compatible pair has 1 and 3 on odd positions and 8 on even positions.
  if((v==2&&pos%2==0)||(v==7&&pos%2==1))continue;
  path.push_back(v+1);dfs(v,unused^(1ULL<<v),odd|(pos%2?(1ULL<<v):0));path.pop_back();
 }
}
int main(int argc,char**argv){
 int lo=argc>1?stoi(argv[1]):25,hi=argc>2?stoi(argv[2]):42;
 if(lo<8 || hi>42 || lo>hi){cerr<<"require 8 <= low <= high <= 42\n";return 2;}
 for(n=lo;n<=hi;n++){
  fill(adj,adj+64,0);for(int a=1;a<=n;a++)for(int b=a+1;b<=n;b++)for(int s=2;s*s<=2*n;s++)if(a+b==s*s){adj[a-1]|=1ULL<<(b-1);adj[b-1]|=1ULL<<(a-1);}
  goal=(n%2?3:8)-1;full=(1ULL<<n)-1;path={1};witnesses.clear();calls=countpaths=0;dfs(0,full^1,1);
  cerr<<n<<" paths="<<countpaths<<" masks="<<witnesses.size()<<" calls="<<calls<<endl;
  cout<<"{\"n\":"<<n<<",\"path_count\":"<<countpaths<<",\"calls\":"<<calls<<",\"witnesses\":[";
  bool comma=false;for(auto &kv:witnesses){if(comma)cout<<",";comma=true;cout<<"{\"odd_mask\":"<<kv.first<<",\"count\":"<<kv.second.second<<",\"path\":[";for(int i=0;i<n;i++){if(i)cout<<",";cout<<kv.second.first[i];}cout<<"]}";}cout<<"]}"<<endl;
 }
}
