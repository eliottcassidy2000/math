#include <bits/stdc++.h>
using namespace std;
int h(const vector<int>& a) {
 int n=a.size(),S=(1<<n)-1,dp[1024][10]={};
 for(int i=0;i<n;i++)dp[1<<i][i]=1;
 for(int x=1;x<=S;x++)for(int i=0;i<n;i++)if(x>>i&1)
   for(int b=a[i]&(S^x);b;b&=b-1){int j=__builtin_ctz(b);dp[x|1<<j][j]+=dp[x][i];}
 return accumulate(dp[S],dp[S]+n,0);
}
vector<int> grow(vector<int> a,int mask) {
 int n=a.size();for(int i=0;i<n;i++)if(!((mask>>i)&1))a[i]|=1<<n;
 a.push_back(mask);return a;
}
void print(const vector<int>&a){for(int z:a)cout<<z<<" ";cout<<"\n";}

int ready(const vector<int>&a) {
 int n=a.size(),F=(1<<n)-1,dp[1024][10]={},H[1024]={};
 H[0]=1;
 for(int i=0;i<n;i++)dp[1<<i][i]=1;
 for(int M=1;M<=F;M++)for(int i=0;i<n;i++)if(M>>i&1) {
  if(M!=(1<<i))for(int b=a[i]&M;b;b&=b-1)dp[M][i]+=dp[M^(1<<i)][__builtin_ctz(b)];
  H[M]+=dp[M][i];
 }
 if(H[F]!=209)throw 1;
 for(int x=0;x<n;x++){
  int r=0,s=0;
  for(int M=1;M<=F;M++)if(M>>x&1)((__builtin_popcount(M)%2)?s:r)+=dp[M][x]*H[F^M];
  if(2*r==208 && s-r==H[F^(1<<x)]){cout<<"prof "<<r<<" "<<s<<" "<<H[F^(1<<x)]<<" ";return x;}
 }
 return -1;
}

int main(){
 vector<vector<int>> A={{6,12,24,49,99,71,15,383,127},{188,509,248,496,32,192,401,16,181}};
 int counts=0;
 for (auto a:A) {
   for(int i=-1;i<9;i++)for(int j=i+1;j<9;j++){
    if(i<0&&j>0)break;
    auto k=a;
    if(i>=0){k[i]^=1<<j;k[j]^=1<<i;}
    for(int b=0;b<512;b++){
      auto m=grow(k,b);
      if(h(m)!=209)continue;
      counts++;
      int x=ready(m);
      if(x>=0){cout<<"root "<<x<<" edit "<<i<<","<<j<<"\n"; print(m);return 0;}
    }
   }
 }
 cout<<"no seed; hits "<<counts<<"\n";
}
