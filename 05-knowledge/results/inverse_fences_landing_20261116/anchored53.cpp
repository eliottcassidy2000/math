#include <bits/stdc++.h>
using namespace std;
int h(const vector<int>& a) {
 int n=a.size(),S=(1<<n)-1,dp[512][9]={};
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

int main(int argc,char**argv){
 if(argc!=2)return 2;
 ifstream f(argv[1]);vector<int> a(7);
 int cand=0,h105=0,matched=0;
 while(f>>a[0]){
  for(int i=1;i<7;i++)f>>a[i];
  if(h(a)!=53)continue;
  cand++;
  for(int j=0;j<128;j++){
   auto b=grow(a,j); if(h(b)>105)continue;
   for(int k=0;k<256;k++){
    auto c=grow(b,k);++matched;
    if(h(c)==105){h105++;cout<<"hit "<<j<<" "<<k<<" ";print(a);print(c);return 0;}
   }
  }
 }
 cout<<"seeds="<<cand<<" tested="<<matched<<" results="<<h105<<"\n";
}
