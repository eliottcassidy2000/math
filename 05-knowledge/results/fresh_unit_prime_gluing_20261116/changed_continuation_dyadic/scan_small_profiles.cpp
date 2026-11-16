#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>
using namespace std;
using V=vector<int>; using ll=long long;
int dp[8192][13],h[8192],route[8192][13];
void require(bool b){if(!b)throw runtime_error("certificate mismatch");}
int total(const V& a){
 int n=a.size(),N=1<<n;h[0]=1;
 for(int s=1;s<N;s++){
  int H=0;
  for(int x=s;x;x&=x-1){int i=__builtin_ctz(x),r=s^(1<<i),z=!r; for(int y=r&a[i];y;y&=y-1)z+=dp[r][__builtin_ctz(y)];dp[s][i]=z;H+=z;}
  h[s]=H;
 }
 return h[N-1];
}
bool strong(const V&a){
 int n=a.size(),F=(1<<n)-1;
 for(int rev=0;rev<2;rev++){
  int m=1,q=1;while(q){
   int i=__builtin_ctz(q);q&=q-1;int b=0;if(!rev)b=a[i];else for(int j=0;j<n;j++)if(a[j]>>i&1)b|=1<<j;
   q|=b&~m;m|=b;}
  if(m!=F)return false;
 }
 return true;
}
uint32_t canon(const V& a){
 V p(a.size());for(int i=0;i<int(p.size());i++)p[i]=i;
 uint32_t best=UINT32_MAX;
 do{uint32_t c=0; for(int i=0;i<int(p.size());i++)for(int j=i+1;j<int(p.size());j++)c=2*c+(a[p[i]]>>p[j]&1); best=min(best,c);}while(next_permutation(p.begin(),p.end()));
 return best;
}
map<int,map<pair<int,uint32_t>,V>> bank;
map<pair<int,int>,int> labeled;
void generate(V a){
 int q=total(a),n=a.size();if(q>45)return;
 if((q==3||q==5||q==9||q==15||q==45)&&strong(a)){labeled[{q,n}]++;auto& b=bank[q]; b.emplace(make_pair(n,canon(a)),a);}
 if(n==8)return;
 for(int mask=0;mask<(1<<(n-1));mask++){
  V c=a;for(int i=0;i<n;i++)if(!(mask>>i&1))c[i]|=1<<n;c.push_back(mask);generate(move(c));
 }
}
int hosts=0,G,a;
map<array<int,4>,int> table; map<array<int,5>,V> witnesses;
void inspect(const V& adj,int fst,int last){
 int n=adj.size(),N=1<<n,F=N-1,first=(1<<fst)-1,tail=F^((1<<(n-last))-1);hosts++;
 require(total(adj)==G);
 ll E[13]={},O[13]={},P[13][13]={};
 for(int i=0;i<n;i++){
   for(int s=1<<i;s<N;s++)if(s>>i&1){
    int r=s^(1<<i),odd=__builtin_popcount(s)&1;
    for(int bits=s;bits;bits&=bits-1){
     int j=__builtin_ctz(bits),v=0;
     if(j==i)v=!r;
     else{int t=s^(1<<j);for(int z=t&~adj[j];z;z&=z-1)v+=route[t][__builtin_ctz(z)];}
     route[s][j]=v;
     if(v){ll w=ll(v)*h[F^s];if(odd)E[i]+=w;else{O[i]+=w;P[i][j]+=w;}}
   }}
 }
 for(int c=1;c<F;c++)if((c&first)&&((F^c)&tail)){
   ll z=G;for(int x=c;x;x&=x-1)for(int y=F^c;y;y&=y-1)z+=2*P[__builtin_ctz(x)][__builtin_ctz(y)];
   ll y=a*z;
   int q=(711<=y&&y<=727&&n+1<=13)?79:(1467<=y&&y<=1483&&n+1<=14)?163:0;
   if(!q)continue;
   ll e=0,o=G;for(int x=c;x;x&=x-1){int i=__builtin_ctz(x);e+=E[i];o+=O[i];}
   array<int,4> key={q,a*G,int(a*e),int(a*o)};
   if((q==163&&a*G==15&&a*e==1107&&a*o==1122)||(q==79&&a*G==5&&a*e==545&&a*o==550)){
    V v=adj;for(int i=0;i<n;i++)if(!(c>>i&1))v[i]|=1<<n;v.push_back(c);
    cerr<<"candidate "<<q<<" "<<y<<" "<<a;for(auto r:v)cerr<<" "<<r;cerr<<"\n";
   }
   table[key]|=1<<((y-9*q)/2);
   array<int,5> w={q,a*G,int(a*e),int(a*o),int(y)};
   if(!witnesses.count(w)){V a2=adj;for(int i=0;i<n;i++)if(!(c>>i&1))a2[i]|=1<<n;a2.push_back(c);witnesses[w]=a2;}
 }

}
V join(const V&a,const V&b){
 int m=a.size(),n=b.size();V c=a;
 for(int &r:c)r|=((1<<n)-1)<<m;for(int r:b)c.push_back(r<<m);return c;
}
void layouts(const vector<V>& chosen,int j,V cur,int f,int l){
 if(j==(int)chosen.size()){
  for(int s=0;s+(int)cur.size()<=13;s++){ V a=cur;for(int i=0;i<s;i++)a=join(a,{0});inspect(a,f,s?1:l); }
  return;
 }
 int n=chosen[j].size();
 for(int s=0;s+(int)cur.size()+n<=13;s++){
  V next=cur;for(int i=0;i<s;i++)next=join(next,{0});
  layouts(chosen,j+1,join(next,chosen[j]),f?f:(s?1:n),n);
 }
}
void factors(const V& f,int j,vector<V> cur){
 if(j==(int)f.size()){layouts(cur,0,{},0,0);return;}
 for(const auto& [key,a]:bank[f[j]]){cur.push_back(a);factors(f,j+1,cur);cur.pop_back();}
}
int main(){
 generate({0});
 for(auto [g0,fs]:map<int,vector<V>>{{1,{{}}},{3,{{3}}},{5,{{5}}},{9,{{9},{3,3}}},{15,{{15},{3,5}}},{45,{{45},{3,15},{5,9},{3,3,5}}}}){
  G=g0;
  for(int m:{1,3,5,9,15}){
   if(m>1 && G*m>45)continue;
   if(m>1 && G*m!=3&&G*m!=5&&G*m!=9&&G*m!=15&&G*m!=45)continue;
   if(m==45)continue; 
   a=m; hosts=0;
   if(G==1){V t={0};for(int n=1;n<=13;n++){inspect(t,1,1);t=join(t,{0});}}
   else for(V f:fs){do{factors(f,0,{});}while(next_permutation(f.begin(),f.end()));}
   cerr<<G<<" * "<<m<<" hosts "<<hosts<<"\n";
  }
 }
 for(auto [k,b]:table){
  auto [q,g,e,o]=k;
  cout<<q<<" "<<g<<" "<<e<<" "<<o<<" "<<b<<"\n";
 }
 for(auto [k,adj]:witnesses){auto [q,g,e,o,h]=k; if(__builtin_popcount(table[{q,g,e,o}])<3)continue;
  cerr<<"witness "<<q<<" "<<g<<" "<<e<<" "<<o<<" "<<h<<" rows";for(int r:adj)cerr<<" "<<r;cerr<<"\n";
 }
}
