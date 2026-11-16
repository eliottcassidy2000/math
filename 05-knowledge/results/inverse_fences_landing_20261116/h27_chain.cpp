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
 vector<int> a={62,124,120,112,32,64,17};
 cout<<"h:" <<h(a)<<"\n";
 auto b0=grow(a,16);
 cout<<"strict seed chain start "<<h(b0)<<"\n";print(a);print(b0);
 for(int i=0;i<256;i++){auto b=grow(b0,i);if(h(b)==105){cout<<"105 extends literal first hit by mask "<<i<<" ";print(b);break;}}
 for(int i=0;i<128;i++){auto b=grow(a,i);if(h(b)!=53)continue;
  for(int j=0;j<256;j++){auto c=grow(b,j);if(h(c)==105){cout<<"two-step: "<<i<<" "<<j<<" ";print(b);print(c);goto done;}}
 }done:
 if(argc==2){
  ifstream f(argv[1]);string line;getline(f,line);
  int cnt=0,cls=0;
  while(getline(f,line)){
   auto z=line.find(',');auto zz=line.find(',',z+1);
   if(stoi(line.substr(z+1,zz-z))!=53)continue;
   ++cls;
   uint64_t c=stoull(line);
   a.resize(8);for(int i=0;i<8;i++)a[i]=(c>>(8*i))&255;
   for(int j=0;j<256;j++)if(h(grow(a,j))==105){
    if(!cnt){cout<<"first from any atlas: "<<j<<" ";print(a);print(grow(a,j));}cnt++;
   }
  }
  cout<<"all "<<cls<<" classes H53 n8 had "<<cnt<<" labelled extension masks\n";
 }
}
