#include<bits/stdc++.h>
using namespace std;
int main(int argc,char**argv){
 ifstream f(argv[1]);string line;getline(f,line);
 const int n=8,F=255;
 map<pair<int,int>,map<int,vector<int>>> bank;
 while(getline(f,line)){
  auto z=line.find(',');auto zz=line.find(',',z+1);
  if(stoi(line.substr(z+1,zz-z))>121)continue;
  uint64_t c=stoull(line);vector<int>a(n);
  for(int i=0;i<n;i++)a[i]=(c>>(8*i))&255;
  int dp[256][8]={},H[256]={};H[0]=1;
  for(int i=0;i<n;i++)dp[1<<i][i]=1;
  for(int M=1;M<=F;M++)for(int i=0;i<n;i++)if(M>>i&1){
    if(M!=(1<<i))for(int b=a[i]&M;b;b&=b-1)dp[M][i]+=dp[M^(1<<i)][__builtin_ctz(b)];
    H[M]+=dp[M][i];
  }
  int S[8]={},R[8]={},W[8][8]={},g=H[F];
  for(int i=0;i<n;i++)for(int M=1;M<=F;M++)if(M>>i&1)
    ((__builtin_popcount(M)%2)?S[i]:R[i])+=dp[M][i]*H[F^M];
  for(int i=0;i<n;i++){
   int pf[256][8]={};pf[1<<i][i]=1;
   for(int M=1;M<=F;M++)if(M>>i&1)for(int j=0;j<n;j++)if(pf[M][j]){
     if(!(__builtin_popcount(M)%2))W[i][j]+=pf[M][j]*H[F^M];
     for(int b=a[j]&~M;b;b&=b-1){int q=__builtin_ctz(b);pf[M|1<<q][q]+=pf[M][j];}
   }
  }
  for(int q=1;q<=F;q++){
   int A=0,B=g;
   for(int i=0;i<n;i++)if(q>>i&1){
     A+=R[i]+S[i];
     for(int j=0;j<n;j++)if(q>>j&1)B-=2*W[i][j];
   }
   auto &t=bank[{g,A}];
   if(t.count(B))continue;
   auto b=a;
   for(int i=0;i<n;i++)if(!(q>>i&1))b[i]|=1<<8;
   b.push_back(q);t[B]=b;
  }
 }
 for(int len=14;len>=3;len--){
  int hits=0;
  for(auto& [p,tab]:bank){
    auto [g,z]=p;
    for(auto& [B,b]:tab){
     if(tab.count(B-2))continue;
     int k=B;while(tab.count(k))k+=2;
     if(k-B<2*len)continue;
     if(hits++>=5)return 0;
     cout<<"len "<<len<<" g "<<g<<" Z "<<z<<" minB "<<B<<"\n";
     for(int i=B;i<k;i+=2){
      for(int a:tab[i])cout<<a<<" ";
      cout<<"\n";
     }
    }
  }
  if(hits)return 0;
 }
}
