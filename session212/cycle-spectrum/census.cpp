#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <cstdint>
#include <cstring>
using namespace std;
int cycles(const vector<int>&a){
 const int n=a.size(),m=n-1,S=1<<m; int dp[256][8];memset(dp,0,sizeof(dp));
 for(int v=1;v<n;v++)if(a[0]>>v&1)dp[1<<(v-1)][v-1]=1;
 for(int mask=1;mask<S;mask++){
  int tmp=mask;
  while(tmp){int v=__builtin_ctz(tmp);tmp&=tmp-1;int prior=mask^(1<<v);if(!prior)continue;
   int possible=prior & ((~(a[v+1]>>1))& (S-1));int s=0;
   while(possible){int u=__builtin_ctz(possible);possible&=possible-1;s+=dp[prior][u];}
   dp[mask][v]=s;
  }
 }
 int s=0;for(int v=1;v<n;v++)if(a[v]&1)s+=dp[S-1][v-1];return s;
}
void output(int n,const map<int,uint64_t>&cnt,const map<int,vector<int>>&w){
 cout<<"{\"order\":"<<n<<",\"labelled_or_extension_presentations\":";uint64_t total=0;for(auto&p:cnt)total+=p.second;cout<<total<<",\"spectrum\":[";bool first=true;for(auto&p:cnt){if(!first)cout<<",";first=false;cout<<p.first;}cout<<"],\"gaps\":[";first=true;for(int i=0;i<=cnt.rbegin()->first;i++)if(!cnt.count(i)){if(!first)cout<<",";first=false;cout<<i;}cout<<"],\"witnesses\":{";first=true;for(auto&p:w){if(!first)cout<<",";first=false;cout<<"\""<<p.first<<"\":[";for(int i=0;i<n;i++){if(i)cout<<",";cout<<p.second[i];}cout<<"]";}cout<<"}}\n";
}
int main(int argc,char**argv){
 map<int,uint64_t>cnt;map<int,vector<int>>w;
 if(argc==3&&string(argv[1])=="--extend8"){
  ifstream in(argv[2]);if(!in)return 2;vector<int>b(8),a(9);int bases=0;
  while(in>>b[0]){for(int i=1;i<8;i++)if(!(in>>b[i]))return 3;bases++;
   for(int m=0;m<256;m++){for(int i=0;i<8;i++)a[i]=b[i]|(((m>>i)&1)<<8);a[8]=m^255;int c=cycles(a);cnt[c]++;if(!w.count(c))w[c]=a;}
  }if(bases!=6880)return 4;output(9,cnt,w);
 }else if(argc==3&&string(argv[1])=="--fixedcycle"){
  int n=stoi(argv[2]);if(n<3||n>8)return 5;vector<int>base(n),a(n);vector<pair<int,int>>free;
  for(int i=0;i<n;i++)base[i]|=1<<((i+1)%n);
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)if(!((base[i]>>j&1)||(base[j]>>i&1)))free.push_back({i,j});
  for(uint64_t bits=0;bits<(1ULL<<free.size());bits++){a=base;for(size_t k=0;k<free.size();k++){auto [i,j]=free[k];if(bits>>k&1)a[i]|=1<<j;else a[j]|=1<<i;}int c=cycles(a);if(!c)return 6;cnt[c]++;if(!w.count(c))w[c]=a;}output(n,cnt,w);
 }else return 1;
}
