#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>
using namespace std;
array<int,28> source,XX,X,Y,W;
array<vector<uint32_t>,2048> groups;
uint32_t three=0,four=0;
uint64_t checks=0;
void gen(int start,int left,uint32_t m,int x,int y) {
 if(!left){groups[((x&127)<<4)|y].push_back(m);return;}
 for(int i=start;i<=28-left;i++)gen(i+1,left-1,m|(1u<<i),x+X[i],y^Y[i]);
}
int wt(uint32_t m){return 2*__builtin_popcount(m)+__builtin_popcount(m&three)+2*__builtin_popcount(m&four);}
int main(int argc,char**argv){
 if(argc!=5){cerr<<"input margin_first margin_last budget\n";return 1;}
 ifstream f(argv[1]);for(int i=0;i<28;i++){if(!(f>>source[i]>>XX[i]>>X[i]>>Y[i]>>W[i]))abort();if(W[i]==3)three|=1u<<i;else if(W[i]==4)four|=1u<<i;else if(W[i]!=2)abort();}
 int budget=atoi(argv[4]);
 for(int k=atoi(argv[2]);k<=atoi(argv[3]);k++){
  auto begin=chrono::steady_clock::now();for(auto &g:groups){vector<uint32_t>().swap(g);}gen(0,k,0,0,0);
  uint64_t total=0;int fibers=0,bad=0,maxc=0;map<int,int>hist;checks=0;
  for(int q=0;q<2048;q++){
   auto &ss=groups[q];if(ss.empty())continue;fibers++;total+=ss.size();vector<int> p(ss.size());iota(p.begin(),p.end(),0);int active=0;
   auto root=[&](int i){while(p[i]!=i){p[i]=p[p[i]];i=p[i];}return i;};
   for(int i=0;i<(int)ss.size();i++){
    active++;
    for(int j=0;j<i;j++){
     checks++;if(wt(ss[i]^ss[j])>budget)continue;
     int a=root(i),b=root(j);if(a!=b){p[a]=b;active--;}
     if(active==1)break;
    }
   }
   hist[active]++;maxc=max(maxc,active);if(active>1){bad++;cerr<<"BAD k="<<k<<" q="<<(q>>4)<<","<<(q&15)<<" states="<<ss.size()<<" components="<<active<<"\n";}
  }
  double sec=chrono::duration<double>(chrono::steady_clock::now()-begin).count();
  cout<<"{\"margin\":"<<k<<",\"budget\":"<<budget<<",\"states\":"<<total<<",\"fibres\":"<<fibers<<",\"disconnected_fibres\":"<<bad<<",\"edge_tests\":"<<checks<<",\"seconds\":"<<sec<<",\"component_counts\":{";bool first=true;for(auto [n,c]:hist){if(!first)cout<<",";first=false;cout<<"\""<<n<<"\":"<<c;}cout<<"}}"<<endl;
 }
}
