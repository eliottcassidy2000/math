#include<bits/stdc++.h>
using namespace std;
using V=vector<int>;using LL=long long;
V d={54,32,2,119,6,20,55};
V make(int u,int w){
 int a=32-__builtin_clz(u),b=32-__builtin_clz(w),m=a+7+b,n=m+1;
 V p(n);int bits=26<<a;
 for(int i=0;i<a;i++)if(u>>(a-i-1)&1)bits|=1<<i;
 for(int i=0;i<b;i++)if(!(w>>i&1))bits|=1<<(a+7+i);
 for(int i=0;i<m;i++){
  for(int j=i+1;j<m;j++)if(i<a||j>=a+7||(d[i-a]>>(j-a)&1))p[i]|=1<<j;else p[j]|=1<<i;
  if(!(bits>>i&1))p[i]|=1<<m;
 }
 p[m]=bits;return p;
}
LL literal(V a){
 int n=a.size(),M=1<<n;vector<array<LL,22>> p(M);
 for(int i=0;i<n;i++)p[1<<i][i]=1;
 for(int m=1;m<M;m++)for(int i=0;i<n;i++)if(p[m][i])
   for(int t=a[i]&~m;t;t&=t-1)p[m|(t&-t)][__builtin_ctz(t)]+=p[m][i];
 return accumulate(p.back().begin(),p.back().end(),0LL);
}
LL F(int n){if(!n)return 1;int i=n/2;return F(i)+(n%2?0:i+i%2);}
vector<array<LL,22>> form(V a){
 int n=a.size(),M=1<<n;
 vector<array<LL,22>> dp(M),Q(n);
 for(int i=0;i<n;i++)dp[1<<i][i]=1;
 auto go=[&](){for(int s=1;s<M;s++)for(int j=0;j<n;j++)if(dp[s][j])
   for(int t=a[j]&~s;t;t&=t-1)dp[s|(t&-t)][__builtin_ctz(t)]+=dp[s][j];};
 go();vector<LL> h(M);h[0]=1;for(int s=1;s<M;s++)h[s]=accumulate(dp[s].begin(),dp[s].end(),0ll);
 for(int i=0;i<n;i++){
  fill(dp.begin(),dp.end(),array<LL,22>{});dp[1<<i][i]=1;go();
  for(int s=1;s<M;s++)if(!(__builtin_popcount(s)&1))
   for(int j=0;j<n;j++)Q[i][j]+=dp[s][j]*h[M-1-s];
 }
 return Q;
}
vector<LL> table(int u,int w){
 auto G=make(u,w);auto Q=form(G);
 int a=32-__builtin_clz(u),b=32-__builtin_clz(w),n=G.size();
 vector<LL> ret(4096);
 LL h=91+5*(F(u)+F(w)-2)+142*(u+w)+256*u*w+2*(w%2);
 if(literal(G)!=h)exit(2);
 for(int m=0;m<4096;m++){
  int mask=(m>>2&127)<<a;
  auto yes=[](int code,int bit){return code==0?0:code==1?1:code==2?bit:1-bit;};
  for(int i=0;i<a;i++)if(yes(m&3,G[n-1]>>i&1))mask|=1<<i;
  for(int i=a+7;i<n-1;i++)if(yes(m>>9&3,G[n-1]>>i&1))mask|=1<<i;
  if(m>>11&1)mask|=1<<(n-1);
  for(int i=0;i<n;i++)if(mask>>i&1)
   for(int j=0;j<n;j++)if(!(mask>>j&1))ret[m]+=2*Q[i][j];
  if((u<=4 && w<=4)&&(m%257==0||m==4095)){
   auto T=G;for(int i=0;i<n;i++)if(!(mask>>i&1))T[i]|=1<<n;T.push_back(mask);
   if(literal(T)!=h+ret[m])exit(3);
  }
 }
 return ret;
}
int main(){
 vector<vector<LL>> S(4096),T(4096);
 for(auto [u,w]:vector<pair<int,int>>{{1,2},{1,4},{3,2},{3,4},{5,4},{7,6},{9,2},{1,10},{11,4},{15,8},{17,2},{3,20},{21,14},{31,16}}){
   auto a=table(u,w),b=table(w,u);
   for(int j=0;j<4096;j++){S[j].push_back(a[j]);T[j].push_back(b[j]);}
 }
 map<vector<LL>,vector<int>> ta;for(int j=0;j<4096;j++)ta[T[j]].push_back(j);
 int hits=0;
 for(int i=0;i<4096;i++)if(ta.count(S[i])){
  for(int j:ta[S[i]]){
   if((i&3)>=2||(i>>9&3)>=2||(j&3)>=2||(j>>9&3)>=2)exit(4);
   ++hits;cout<<i<<" "<<j<<" masks "<<(i>>2&127)<<" "<<(j>>2&127)<<" shift ";
   for(auto a:S[i])cout<<a<<",";
   cout<<"\n";
  }
 }
 if(hits!=9)return 5;
 cout<<"PASS: 9 surviving pairs, all constant on both flanks; no digital callback\n";
}
