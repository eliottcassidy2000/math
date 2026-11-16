#include <bits/stdc++.h>
using namespace std;
int main(){
 map<int,set<pair<int,int>>> result;
 for(int n=1;n<=6;n++){
  vector<pair<int,int>> edges;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)edges.emplace_back(i,j);
  for(unsigned u=0;u<(1u<<edges.size());u++){
   int a[6]={};for(unsigned e=0;e<edges.size();e++){auto [i,j]=edges[e];if(u>>e&1)a[i]|=1<<j;else a[j]|=1<<i;}
   int F=(1<<n)-1;int p[64][6]={},h[64]={};h[0]=1;
   for(int S=1;S<=F;S++)for(int b=S;b;b&=b-1){int i=__builtin_ctz((unsigned)b),R=S^(1<<i),v=!R;for(int c=R&a[i];c;c&=c-1)v+=p[R][__builtin_ctz((unsigned)c)];p[S][i]=v;h[S]+=v;}
   if(h[F]!=3 && h[F]!=5 && h[F]!=15 && h[F]!=9)continue;
   bool strong=true;for(int i=0;i<n;i++){int S=1<<i,old=0;while(old!=S){old=S;for(int b=S;b;b&=b-1)S|=a[__builtin_ctz((unsigned)b)];}strong &= S==F;}
   if(!strong)continue;
   int X=0;for(int S=0;S<=F;S++)X+=h[S]*h[F^S];result[h[F]].insert({n,X});
   if(h[F]==9){
    int mask=0;for(int i=0;i<n;i++){
      int W=0;for(int S=1;S<=F;S++)if(S>>i&1)W+=p[S][i]*h[S^F];
      mask|=1<<(W%3);
    }
    static set<pair<int,int>> seen;
    if(seen.insert({X,mask}).second)cout<<"H9 mod3 starts "<<X<<" "<<mask<<"\n";
   }
  }
 }
 for(auto &[h,ss]:result){cout<<h<<":";for(auto[n,x]:ss)cout<<" "<<n<<","<<x;cout<<"\n";}
}
