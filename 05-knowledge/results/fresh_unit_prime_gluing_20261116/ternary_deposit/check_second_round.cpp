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
int hosts=0, mult=1, g=45; map<pair<int,int>,map<int,int>> report; map<array<int,3>,V> examples;
void inspect(const V& a,int fst,int last){
 int n=a.size(),N=1<<n,F=N-1,first=(1<<fst)-1,tail=F^((1<<(n-last))-1);hosts++;
 require(total(a)==g);
 ll E[13]={},O[13]={},P[13][13]={};
 for(int i=0;i<n;i++){
   for(int s=1<<i;s<N;s++)if(s>>i&1){
    int r=s^(1<<i),odd=__builtin_popcount(s)&1;
    for(int bits=s;bits;bits&=bits-1){
     int j=__builtin_ctz(bits),v=0;
     if(j==i)v=!r;
     else{int t=s^(1<<j);for(int z=t&~a[j];z;z&=z-1)v+=route[t][__builtin_ctz(z)];}
     route[s][j]=v;
     if(v){ll w=ll(v)*h[F^s];if(odd)E[i]+=w;else{O[i]+=w;P[i][j]+=w;}}
   }}
 }
 for(int c=1;c<F;c++)if((c&first)&&((F^c)&tail)){
   ll e=0,o=g;for(int x=c;x;x&=x-1){int i=__builtin_ctz(x);e+=E[i];o+=O[i];}
   int q=(mult*e==378&&mult*o==423)?79:(mult*e==756&&mult*o==801)?163:0;
   if(!q)continue;
   int f=3, idx=3, up=0;
   const int base[]={3,5,9};
   for(;idx<=14;idx++){
    f=base[(idx-3)%3]; for(int j=0;j<(idx-3)/3;j++)f*=5;
    if(mult*f>2*mult*o-45)break;
    up=idx;
   }
   if(n+1>up)continue;
   ll z=g;for(int x=c;x;x&=x-1)for(int y=F^c;y;y&=y-1)z+=2*P[__builtin_ctz(x)][__builtin_ctz(y)];
   z*=mult;
   int val=static_cast<int>(z);
   report[{q,mult}][val]++;
   if(!examples.count({q,mult,val})){
    V a2=a;for(int i=0;i<n;i++)if(!(c>>i&1))a2[i]|=1<<n;a2.push_back(c);require(strong(a2));examples[{q,mult,val}]=a2;
   }
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
 for(int q:{3,5,9,15,45}){cout<<"type "<<q<<" "<<bank[q].size()<<" labeled";for(int n=3;n<=8;n++)cout<<" "<<labeled[{q,n}];cout<<"\n";}
 require(bank[3].size()==1&&bank[5].size()==1&&bank[9].size()==2&&bank[15].size()==4&&bank[45].size()==12);
 for(auto [a,fs]:map<int,vector<V>>{{1,{{45},{3,15},{5,9},{3,3,5}}},{3,{{15},{3,5}}},{9,{{5}}}}){
   mult=a;g=45/a;hosts=0;
   for(V f:fs){do{factors(f,0,{});}while(next_permutation(f.begin(),f.end()));}
   cout<<"hosts "<<a<<" "<<hosts<<"\n";
   require(hosts==(a==1?926:(a==3?330:55)));
 }
 for(auto [qm,d]:report)for(auto [z,c]:d){auto [q,m]=qm;cout<<q<<" "<<m<<" "<<z<<" "<<c<<" example";for(int r:examples.at({q,m,z}))cout<<" "<<r;cout<<"\n";}
 map<pair<int,int>,set<int>> expected={
  {{79,1},{573,717,723,735,801}}, {{79,3},{609,711,717,723,801}}, {{79,9},{711}},
  {{163,1},{1027,1053,1151,1205,1263,1305,1329,1365,1401,1467,1473,1479,1491,1557}},
  {{163,3},{999,1053,1305,1365,1467,1473,1479,1557}}, {{163,9},{1053,1467}}
 };
 map<pair<int,int>,set<int>> support;
 for(auto [qm,sub]:report)for(auto [z,c]:sub)support[qm].insert(z);
 require(support==expected);

 map<int,set<int>> actual;
 for(const auto& [key,m]:report)for(const auto& [H,_]:m)actual[key.first].insert(H);
 const map<int,set<int>> expected_union{
   {79,{573,609,711,717,723,735,801}},
   {163,{999,1027,1053,1151,1205,1263,1305,1329,1365,1401,1467,1473,1479,1491,1557}}
 };
 require(actual==expected_union);
 const map<int,map<int,int>> target{{79,{{717,5},{723,8}}},
                                   {163,{{1467,2},{1473,10},{1479,13}}}};
 for(auto [q,m]:target){
  map<int,int> rows;
  for(auto [H,num]:report[{q,1}])if(H>=9*q&&H<=9*q+16)rows.emplace(H,num);
  require(rows==m);
 }
}
