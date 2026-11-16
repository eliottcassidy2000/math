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

int ready(const vector<int>&a) {
 int n=a.size(),F=(1<<n)-1,dp[512][9]={},H[512]={};
 H[0]=1;
 for(int i=0;i<n;i++)dp[1<<i][i]=1;
 for(int M=1;M<=F;M++)for(int i=0;i<n;i++)if(M>>i&1) {
  if(M!=(1<<i))for(int b=a[i]&M;b;b&=b-1)dp[M][i]+=dp[M^(1<<i)][__builtin_ctz(b)];
  H[M]+=dp[M][i];
 }
 if(H[F]!=105)throw 1;
 for(int x=0;x<n;x++){
  int r=0,s=0;
  for(int M=1;M<=F;M++)if(M>>x&1)((__builtin_popcount(M)%2)?s:r)+=dp[M][x]*H[F^M];
  if(2*r==104 && s-r==H[F^(1<<x)]){cout<<"prof "<<r<<" "<<s<<" "<<H[F^(1<<x)]<<" ";return x;}
 }
 return -1;
}

int main(int argc,char **argv){
 ifstream f(argv[1]);string line;getline(f,line);int sources=0,found=0,readyfound=0;
 while(getline(f,line)){
  auto z=line.find(',');auto zz=line.find(',',z+1);
  if(stoi(line.substr(z+1,zz-z))>105)continue;
  ++sources;
  uint64_t c=stoull(line);
  vector<int>a(8);for(int i=0;i<8;i++)a[i]=(c>>(8*i))&255;
  for(int j=0;j<256;j++){
   auto b=grow(a,j);if(h(b)!=105)continue;found++;
   int x=ready(b);
   if(x>=0){cout<<"ready port "<<x<<" grown by "<<j<<" from h="<<h(a)<<"\n";print(b);return 0;}
  }
 }
 cout<<"none sources="<<sources<<" hits="<<found<<"\n";
}
