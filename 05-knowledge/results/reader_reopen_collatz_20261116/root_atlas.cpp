#include<bits/stdc++.h>
using namespace std;
int main(int argc,char**argv){
 if(argc!=2)return 1;
 ifstream f(argv[1]);string line;getline(f,line);
 const int n=8,F=(1<<8)-1;
 map<int,tuple<int,int,int,int,int,vector<int>>> best;
 int read=0;
 while(getline(f,line)){
  read++;uint64_t c=stoull(line);
  vector<int>a(n);for(int i=0;i<n;i++)a[i]=(c>>(8*i))&255;
  int dp[256][8]={},H[256]={};H[0]=1;
  for(int i=0;i<n;i++)dp[1<<i][i]=1;
  for(int M=1;M<=F;M++)for(int i=0;i<n;i++)if(M>>i&1){
    if(M!=(1<<i))for(int b=a[i]&M;b;b&=b-1)dp[M][i]+=dp[M^(1<<i)][__builtin_ctz(b)];
    H[M]+=dp[M][i];
  }
  int S[8]={},R[8]={},W[8][8]={},g=H[F];
  for(int i=0;i<n;i++)for(int M=1;M<=F;M++)if(M>>i&1)
    ((__builtin_popcount(M)%2)?S[i]:R[i])+=dp[M][i]*H[F^M];
  vector<int>cand;
  for(int q=1;q<=F;q++){
    int A=0;for(int i=0;i<n;i++)if(q>>i&1)A+=R[i]+S[i];
    if(A>0&&(A>>__builtin_ctz(A))==121)cand.push_back(q);
  }
  if(cand.empty())continue;
  for(int i=0;i<n;i++){
   int pf[256][8]={};pf[1<<i][i]=1;
   for(int M=1;M<=F;M++)if(M>>i&1)for(int j=0;j<n;j++)if(pf[M][j]){
     if(!(__builtin_popcount(M)%2))W[i][j]+=pf[M][j]*H[F^M];
     for(int b=a[j]&~M;b;b&=b-1){int q=__builtin_ctz(b);pf[M|1<<q][q]+=pf[M][j];}
   }
  }
  for(int q:cand){
   int h=g,r=0,s=g;
   for(int i=0;i<n;i++)if(q>>i&1){
     s+=R[i];r+=S[i];
     for(int j=0;j<n;j++)if(!(q>>j&1))h+=2*W[i][j];
   }
   int A=r+s-g,B=h-2*s+2*g;
   if(abs(B)>400||(best.count(B)&&get<0>(best[B])<=A)||A>121*32)continue;
   auto b=a;
   for(int i=0;i<n;i++)if(!(q>>i&1))b[i]|=1<<8;
   b.push_back(q);
   best[B]=make_tuple(A,h,g,r,s,b);
  }
 }
 cout<<"Atlas parent rows "<<read<<"\n";
 for(auto [B,z]:best){
   auto [A,h,g,r,s,a]=z;
   cout<<B<<" A"<<A<<" x8 h"<<h<<" g"<<g<<" r"<<r<<" s"<<s<<" :";
   for(auto j:a)cout<<" "<<j;cout<<"\n";
 }
}
