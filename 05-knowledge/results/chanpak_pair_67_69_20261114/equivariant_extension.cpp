#include <bits/stdc++.h>
using namespace std;using LL=long long;using A=vector<int>;
const A P{2,3,5,6,4,0,1,8,9,7};
void put(A& a,int i,int j){a[i]|=1<<j;}
vector<A> dp(const A& a){
 int n=a.size(),F=1<<n;vector<A>d(F,A(n));
 for(int i=0;i<n;i++)d[1<<i][i]=1;
 for(int s=1;s<F;s++)for(int i=0;i<n;i++)if(d[s][i])
  for(int b=a[i]&~s;b;b&=b-1){
   int j=__builtin_ctz(b);d[s|(1<<j)][j]+=d[s][i];
  }
 return d;
}
bool strong(A a){
 int n=a.size(),F=(1<<n)-1;
 for(int i=0;i<n;i++){int x=1<<i,y=0;while(x!=y){y=x;for(int j=0;j<n;j++)if(x>>j&1)x|=a[j];}if(x!=F)return false;}
 return true;
}
int main(){
 vector<pair<int,int>> rep{{7,0},{7,2},{7,5},{7,1},{7,3},{7,6},{7,4}};
 int candidates=0,strong_count=0;
 for(int z=0;z<256;z++){
  A a{14,28,104,112,37,67,19,0,0,0};
  for(int i=0;i<7;i++){
   auto [x,y]=rep[i];for(int j=0;j<3;j++){if(z>>i&1)put(a,x,y);else put(a,y,x);x=P[x];y=P[y];}
  }
  for(int i=7;i<10;i++){if(z>>7&1)put(a,i,P[i]);else put(a,P[i],i);}
  if(!strong(a))continue;strong_count++;
  auto b=a;for(int i=0;i<10;i++)b[i]=1023^(1<<i)^a[i];
  auto e=dp(a),s=dp(b);LL M[10][10]={};
  for(int t=1;t<1023;t++)for(int i=0;i<10;i++)if(e[t][i])
    for(int j=0;j<10;j++)M[i][j]+=LL(e[t][i])*s[1023^t][j];
  for(int X=1;X<1023;X++){
    int U=X,V=0;
    for(int i=0;i<10;i++)if(U>>i&1)V|=1<<P[i];
    LL d=0;
    for(int i=0;i<10;i++)if(U>>i&1)for(int j=0;j<10;j++)if(V>>j&1)d+=M[j][i]-M[i][j];
    if(abs(d)==2){
     candidates++;cout<<z<<" h "<<accumulate(e[1023].begin(),e[1023].end(),0LL)<<" U,V "<<U<<","<<V<<" skew "<<d<<" rows";
     for(auto v:a)cout<<" "<<v;
     cout<<"\n";
    }
  }
 }
 cout<<"strong "<<strong_count<<" survivors "<<candidates<<"\n";
}
