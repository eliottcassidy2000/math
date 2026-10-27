#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>
#include <stdexcept>
using namespace std;
using U=uint64_t;
static vector<uint32_t> A;
static int n;
static U hamilton(){
 size_t sz=1u<<n; vector<U>d(sz*n,0);
 for(uint32_t s=1;s<sz;s++)for(uint32_t vs=s;vs;vs&=vs-1){
  int v=__builtin_ctz(vs); uint32_t prev=s^(1u<<v); U total=0;
  if(!prev) total=1;
  else for(uint32_t ps=prev&~A[v];ps;ps&=ps-1)total+=d[(size_t)prev*n+__builtin_ctz(ps)];
  d[(size_t)s*n+v]=total;
 }
 U total=0;for(int v=0;v<n;v++)total+=d[(sz-1)*n+v];return total;
}
static U automorphisms(){
 vector<int>p(n,-1),degree(n);for(int i=0;i<n;i++)degree[i]=__builtin_popcount(A[i]);
 auto rec=[&](auto&&self,int i,uint32_t used)->U{
  if(i==n)return 1;
  U total=0;
  for(int v=0;v<n;v++)if(!(used>>v&1)&&degree[i]==degree[v]){
   bool good=true;for(int j=0;j<i;j++)if(((A[i]>>j)&1)!=((A[v]>>p[j])&1)){good=false;break;}
   if(good){p[i]=v;total+=self(self,i+1,used|(1u<<v));}
  }return total;
 };return rec(rec,0,0);
}
static U bruteH(){vector<int>p(n);iota(p.begin(),p.end(),0);U total=0;do{bool good=true;for(int j=1;j<n;j++)if(!(A[p[j-1]]>>p[j]&1)){good=false;break;}if(good)total++;}while(next_permutation(p.begin(),p.end()));return total;}
int main(int argc, char**){
 if(argc>1){
  while(cin>>n){ A.resize(n);for(auto &x:A)cin>>x;U H=hamilton(),G=automorphisms();if(H%G)throw runtime_error("nonintegral tc");cout<<n<<" H="<<H<<" Aut="<<G<<" tc="<<H/G<<endl;}return 0;
 }
 for(auto [a,b]:vector<pair<int,int>>{{3,1},{5,1},{7,1},{9,1},{3,3},{11,1},{13,1},{15,1},{17,1}}){
  n=a*b;A.resize(n);auto id=[&](int x,int y){return ((x%a+a)%a)*b+(y%b+b)%b;};
  vector<pair<int,int>>pairs;vector<bool>seen(n,false);seen[0]=true;
  for(int x=0;x<a;x++)for(int y=0;y<b;y++)if(!seen[id(x,y)]){seen[id(x,y)]=seen[id(-x,-y)]=true;pairs.push_back({x,y});}
  U mn=~U(0),sets=1u<<pairs.size();map<U,U>small;
  for(U s=0;s<sets;s++){
   fill(A.begin(),A.end(),0);for(int x=0;x<a;x++)for(int y=0;y<b;y++)for(size_t j=0;j<pairs.size();j++){auto [dx,dy]=pairs[j];if(s>>j&1){dx=-dx;dy=-dy;}A[id(x,y)]|=1u<<id(x+dx,y+dy);}
   U H=hamilton(),G=automorphisms();if(H%G)throw runtime_error("nonintegral tc");if(n<=7&&bruteH()!=H)throw runtime_error("brute H mismatch");U tc=H/G;mn=min(mn,tc);if(tc<=21)small[tc]++;
  }
  cout<<"Z"<<a<<"xZ"<<b<<" sets="<<sets<<" min_tc="<<mn<<" small";for(auto[k,v]:small)cout<<" "<<k<<":"<<v;cout<<endl;
 }
}
