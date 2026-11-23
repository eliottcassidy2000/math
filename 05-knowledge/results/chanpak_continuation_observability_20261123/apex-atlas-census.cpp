#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <numeric>
using U=uint64_t;
struct T {int n; std::array<int,8> r{};};
std::vector<int> deg(T const&t){std::vector<int>D;for(int i=0;i<t.n;i++)D.push_back(__builtin_popcount(t.r[i]));return D;}
U canon(T const&t){
 auto dg=deg(t);std::vector<int> p(t.n);
 std::iota(p.begin(),p.end(),0);std::stable_sort(p.begin(),p.end(),[&](int a,int b){return dg[a]<dg[b];});
 U best=~U(0);
 auto f=[&](auto&& self,int i)->void{
   if(i==t.n){U m=0;int l=0;for(int a=0;a<t.n;a++)for(int b=a+1;b<t.n;b++,l++)if((t.r[p[a]]>>p[b])&1)m|=1ULL<<l;best=std::min(best,m);return;}
   for(int j=i;j<t.n&&dg[p[j]]==dg[p[i]];j++){std::swap(p[i],p[j]);self(self,i+1);std::swap(p[i],p[j]);}
 };
 f(f,0);return best;
}
T decode(int n,U m){T t{n};for(int a=0,l=0;a<n;a++)for(int b=a+1;b<n;b++,l++){if((m>>l)&1)t.r[a]|=1<<b;else t.r[b]|=1<<a;}return t;}
std::string canon_prof(T const&t){
 const int n=t.n,M=(1<<n)-1;
 int in[8]={},h[256]={1},V[256][8]={};int B[8][8]{};
 for(int i=0;i<n;i++)for(int v=0;v<n;v++)if((t.r[i]>>v)&1)in[v]|=1<<i;
 for(int s=1;s<=M;s++)for(unsigned m=s;m;m&=m-1) {
    auto i=__builtin_ctz(m);int prev=s^(1<<i);
    if(!prev)V[s][i]=1;
    else for(unsigned u=prev&in[i];u;u&=u-1)V[s][i]+=V[prev][__builtin_ctz(u)];
    h[s]+=V[s][i];
 }
 for(int i=0;i<n;i++){
   int A[256][8]={};A[1<<i][i]=1;
   for(int s=(1<<i)+1;s<=M;s++){
     if(!((s>>i)&1))continue;
     for(unsigned m=s^(1<<i);m;m&=m-1){
       int j=__builtin_ctz(m),p=s^(1<<j);
       for(unsigned u=p&in[j];u;u&=u-1)A[s][j]+=A[p][__builtin_ctz(u)];
       if(!(__builtin_popcount(s)&1))B[i][j]+=A[s][j]*h[M^s];
     }
   }
 }
 std::vector<int> a(n),b(n);int W[8][8]{};
 for(int i=0;i<n;i++)for(int j=0;j<n;j++){a[i]+=B[i][j];b[i]+=B[j][i];W[i][j]=B[i][j]+B[j][i];}
 std::vector<int> p(n);std::iota(p.begin(),p.end(),0);
 auto ty=[&](int i){return std::pair(a[i],b[i]);};
 std::stable_sort(p.begin(),p.end(),[&](int i,int j){return ty(i)<ty(j);});
 std::string best;
 auto add=[](std::string&s,int x){s.push_back(char((x>>8)&255));s.push_back(char(x&255));};
 auto f=[&](auto&&self,int i)->void{
  if(i==n){std::string s;add(s,h[M]);for(int j:p){add(s,a[j]);add(s,b[j]);}for(int j=0;j<n;j++)for(int k=j+1;k<n;k++)add(s,W[p[j]][p[k]]);
     if(best.empty()||s<best)best=s;
     return;
  }
  for(int j=i;j<n&&ty(p[i])==ty(p[j]);j++){std::swap(p[i],p[j]);self(self,i+1);std::swap(p[i],p[j]);}
 };
 f(f,0);return best;
}
int main(int argc,char**argv){
 int max=argc>1?std::stoi(argv[1]):7;
 if(max<1||max>8)return 1;
 std::vector<U> L{0};
 for(int n=1;n<=max;n++){
  if(n>1){
   std::unordered_map<U,int> N;
   for(U m:L)for(int j=0;j<(1<<(n-1));j++){
      auto t=decode(n-1,m);t.n=n;
      for(int a=0;a<n-1;a++)if((j>>a)&1)t.r[a]|=1<<(n-1);else t.r[n-1]|=1<<a;
      N.emplace(canon(t),0);
   }
   L.clear();for(auto [m,k]:N)L.push_back(m);
  }
  std::sort(L.begin(),L.end());
  std::unordered_map<std::string,U> k;int coll=0;
  for(U m:L){
    auto s=canon_prof(decode(n,m));
    auto [it,newval]=k.emplace(s,m);
    if(!newval){if(coll<12)std::cout<<"collision n="<<n<<": "<<it->second<<" "<<m<<"\n";++coll;}
  }
  std::cout<<"n="<<n<<" types "<<L.size()<<" profiles "<<k.size()<<" collisions "<<coll<<std::endl;
 }
}
