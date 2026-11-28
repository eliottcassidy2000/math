#include <algorithm>
#include <array>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <unordered_set>
#include <vector>
#include <set>
#include <string>
using namespace std;
using V=array<int,4>;
struct Hash{
 size_t operator()(const V &v) const{
  size_t h=0x9e3779b97f4a7c15ULL;
  for(auto x:v)h^=(x+0x9e3779b97f4a7c15ULL)+(h<<6)+(h>>2);
  return h;
 }
};
int main(){
 const int X=2647;
 unordered_set<V,Hash> seen;
 vector<V> queue;
 vector<int> parent,edge, best(X+1,-1);
 auto add=[&](V w,int from,int op){
  auto [h,l,u,g]=w;
  if(h>X || u>X)return;
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
  if(2*u<=X) add({h,0,2*u,g},i,2);
 }
 set<int> missing={13,41,43,55,67,77,83,163,173,231,233,247,251,277,299,421,449,503,533,569,583,637,821,1067,1339,1861,2647};
 auto word=[&](int p){
  string w;
  while(parent[p]>=0){w.push_back("KMP"[edge[p]]);p=parent[p];}
  reverse(w.begin(),w.end());
  return w.empty()?string("-"):w;
 };
 for(int i=0;i<(int)queue.size();i++){
  auto [a,e,u,g]=queue[i];
  if(a*a>X)continue;
  for(int j=0;j<(int)queue.size();j++){
   auto [b,f,v,h]=queue[j];
   if(a>b || a*b>X)continue;
   int c=a*b-e*f+u*v;
   if(missing.erase(c))cout<<c<<" "<<word(i)<<" "<<word(j)<<"\n";
  }
 }
 cerr<<"Bridge states "<<queue.size()<<" residual:";
 for(int c:missing)cerr<<" "<<c;
 cerr<<"\n";
}
