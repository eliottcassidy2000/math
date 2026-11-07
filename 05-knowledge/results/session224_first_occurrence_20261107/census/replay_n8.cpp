#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <vector>
using namespace std;using U=uint64_t;
pair<U,int> canon(array<unsigned,8> a){
 vector<vector<int>> cells;for(int d=0;d<8;d++){vector<int> c;for(int i=0;i<8;i++)if(__builtin_popcount(a[i])==d)c.push_back(i);if(!c.empty())cells.push_back(c);}
 array<int,8> p{};U best=~U(0);int multiplicity=0;
 auto rec=[&](auto&& self,int k,int start)->void{
   if(k==int(cells.size())){U code=0;for(int i=0;i<8;i++)for(int j=0;j<8;j++)if((a[p[i]]>>p[j])&1)code|=U(1)<<(8*i+j);if(code<best){best=code;multiplicity=1;}else if(code==best)multiplicity++;return;}
   auto c=cells[k];do{for(int j=0;j<int(c.size());j++)p[start+j]=c[j];self(self,k+1,start+c.size());}while(next_permutation(c.begin(),c.end()));
 };rec(rec,0,0);return {best,multiplicity};
}
array<unsigned,8> decode(U x){array<unsigned,8>a{};for(int i=0;i<8;i++)a[i]=(x>>(8*i))&255;return a;}
int hamilton(array<unsigned,8>a){int dp[256][8]={};for(int i=0;i<8;i++)dp[1<<i][i]=1;for(int s=1;s<256;s++)for(int i=0;i<8;i++)if((s>>i)&1){unsigned t=a[i]&~s;while(t){int j=__builtin_ctz(t);t&=t-1;dp[s|(1<<j)][j]+=dp[s][i];}}int sum=0;for(int i=0;i<8;i++)sum+=dp[255][i];return sum;}
bool strong(array<unsigned,8>a){for(int i=0;i<8;i++)a[i]|=1<<i;for(int k=0;k<8;k++)for(int i=0;i<8;i++)if((a[i]>>k)&1)a[i]|=a[k];for(int i=0;i<8;i++)if(a[i]!=255)return false;return true;}
int main(int argc,char**argv){ifstream f(argv[1]);map<U,int> classes;array<unsigned,8>a{};int n7=0;while(f>>a[0]){for(int i=1;i<7;i++)f>>a[i];n7++;for(int x=0;x<128;x++){auto b=a;b[7]=x;for(int i=0;i<7;i++)if(!((x>>i)&1))b[i]|=128;auto [code,aut]=canon(b);classes.emplace(code,aut);}}cerr<<"n7="<<n7<<" classes="<<classes.size()<<"\n";cout<<"code,H,aut,self_converse,strong\n";for(auto [code,aut]:classes){a=decode(code);auto b=a;for(int i=0;i<8;i++)b[i]^=255^(1<<i);auto dual=canon(b).first;cout<<code<<","<<hamilton(a)<<","<<aut<<","<<(dual==code)<<","<<strong(a)<<"\n";} }
