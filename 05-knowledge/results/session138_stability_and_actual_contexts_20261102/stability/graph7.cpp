#include <iostream>
#include <vector>
#include <cstdint>
#include <array>
using namespace std;
int main(){
 const int n=7,m=21;vector<pair<int,int>>e;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)e.push_back({i,j});
 array<unsigned,21>disj{};for(int i=0;i<m;i++)for(int j=0;j<m;j++)if(e[i].first!=e[j].first&&e[i].first!=e[j].second&&e[i].second!=e[j].first&&e[i].second!=e[j].second)disj[i]|=1u<<j;
 int mn=100000,zeros=0;unsigned witness=0;uint64_t count=0;
 for(unsigned mask=0;mask<(1u<<m);mask++){
  int E=__builtin_popcount(mask),M2=0,M3=0;
  for(unsigned a=mask;a;a&=a-1){int i=__builtin_ctz(a);unsigned rest=mask&disj[i]&(~0u<<(i+1));M2+=__builtin_popcount(rest);
   for(unsigned b=rest;b;b&=b-1){int j=__builtin_ctz(b);M3+=__builtin_popcount(rest&disj[j]&(~0u<<(j+1)));}
  }
  int W=E*(E-1)/2-M2;int D=22*E-24*M2+6*M3+15*W;
  if(D<mn){mn=D;witness=mask;}if(D==0)zeros++;count++;
 }
 cout<<"{\"n\":7,\"graphs\":"<<count<<",\"minimum\":"<<mn<<",\"zero_graphs\":"<<zeros<<",\"witness\":"<<witness<<"}\n";
}
