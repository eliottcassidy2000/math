#include <algorithm>
#include <array>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <unordered_set>
#include <vector>
using namespace std;
using V=array<int,4>;
struct Hash{
 size_t operator()(const V &v) const{
  size_t h=0x9e3779b97f4a7c15ULL;
  for(auto x:v)h^=(x+0x9e3779b97f4a7c15ULL)+(h<<6)+(h>>2);
  return h;
 }
};
int main(int argc,char**argv){
 assert(argc==2); int X=atoi(argv[1]);
 unordered_set<V,Hash> seen;
 vector<V> queue;
 vector<int> parent,edge, best(X+1,-1);
 auto add=[&](V w,int from,int op){
  auto [h,l,u,g]=w;
  if(h>X)return;
  assert(0<=l&&l<=h&&u>=l&&g<=h);
  if(!seen.insert(w).second)return;
  if(best[h]<0)best[h]=queue.size();
  queue.push_back(w);parent.push_back(from);edge.push_back(op);
 };
 // g=H(T-a); empty deletion has one HP.
 add({1,1,1,1},-1,-1);
 for(int i=0;i<(int)queue.size();i++){
  auto [h,l,u,g]=queue[i];
  if(h-l+u<=X){
   add({h-l+u,g,2*u-g,g},i,0);
   add({h-l+u,h-l,h+u,h},i,1);
  }
  if(h+2*u<=X) add({h,0,2*u,g},i,2);
 }
 cerr<<"States "<<queue.size()<<"\n";
 int count=0;
 for(int h=1;h<=X;h+=2)if(best[h]>=0){
  ++count;
  int p=best[h];string ops;
  while(parent[p]>=0){ops.push_back("KMP"[edge[p]]);p=parent[p];}
  reverse(ops.begin(),ops.end());
  cout<<h<<" "<<ops<<"\n";
 }
 cerr<<"Distinct "<<count<<"\n";
}
