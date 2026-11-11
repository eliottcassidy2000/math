// Exact order-nine screen. Hashes are only lookup accelerators: keys are full vectors.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <vector>
using namespace std;
using I=int64_t;
struct Hash {size_t operator()(const vector<I>& a)const{size_t h=0;for(auto x:a)h=(h^(x+0x9e3779b97f4a7c15ULL+(h<<6)+(h>>2)));return h;}};
void need(bool b,const char*s){if(!b)throw runtime_error(s);}
int main(int argc,char**argv){
 need(argc==2,"screen9 input");ifstream f(argv[1]);uint64_t code,dual;int aut,scH;array<int,9> rows;
 unordered_map<vector<I>,uint64_t,Hash> supports,distributions;uint64_t ordinary=0,merged=0,mass=0,unitSources=0,unitSets=0,unitOrbits=0,maxWidth=0,sumWidth=0;
 set<I> hvalues,tcvalues,coefficients;vector<array<uint64_t,3>> collisions;vector<array<uint64_t,3>> distributionCollisions;
 while(f>>code>>dual>>aut>>scH){for(int &r:rows)f>>r;bool sc=code==dual;int weight=sc?1:2;ordinary+=weight;merged++;
  constexpr int n=9,N=512,full=511;array<array<I,9>,N> end{},start{};
  for(int v=0;v<n;v++)end[1<<v][v]=start[1<<v][v]=1;
  for(int S=1;S<N;S++)if(S&(S-1))for(int v=0;v<n;v++)if(S>>v&1){int R=S^(1<<v);for(int u=0;u<n;u++)if(R>>u&1){if(rows[u]>>v&1)end[S][v]+=end[R][u];if(rows[v]>>u&1)start[S][v]+=start[R][u];}}
  I H=accumulate(end[full].begin(),end[full].end(),I(0));need(H%aut==0,"automorphism divisibility");if(sc)need(H==scH,"inherited SC Hamilton count");mass+=weight*(H/aut);hvalues.insert(H);tcvalues.insert(H/aut);
  array<array<I,9>,9> C{};for(int S=1;S<full;S++)for(int i=0;i<n;i++)if(S>>i&1)for(int j=0;j<n;j++)if(!(S>>j&1))C[i][j]+=end[S][i]*start[full^S][j];
  array<I,9> uv{},dualuv{};for(int j=0;j<n;j++){uv[j]=start[full][j];dualuv[j]=end[full][j];for(int i=0;i<n;i++){uv[j]+=C[i][j];dualuv[j]+=C[j][i];}}
  array<I,N> U{},dualU{},W{},E{};vector<I> spectrum(N);spectrum[0]=H;int seeds=0;
  for(int S=1;S<N;S++){int j=__builtin_ctz((unsigned)S),R=S^(1<<j);U[S]=U[R]+uv[j];dualU[S]=dualU[R]+dualuv[j];W[S]=W[R];E[S]=E[R]+end[full][j];for(int i=0;i<n;i++)if(R>>i&1)W[S]+=C[i][j]+C[j][i];spectrum[S]=H+U[S]-E[S]-W[S];
   if(W[S]==H-1&&U[S]>0){seeds++;I c=U[S];while(c%2==0)c/=2;coefficients.insert(c);need(dualU[S]>0,"positive dual unit slope");c=dualU[S];while(c%2==0)c/=2;coefficients.insert(c);}
  }
  need(seeds%aut==0,"unit orbits");unitSources+=weight*(seeds>0);unitSets+=weight*seeds;unitOrbits+=weight*(seeds/aut);
  sort(spectrum.begin(),spectrum.end());
  auto [itd,newd]=distributions.emplace(spectrum,code);if(!newd)distributionCollisions.push_back({itd->second,code,(uint64_t)H});
  spectrum.erase(unique(spectrum.begin(),spectrum.end()),spectrum.end());maxWidth=max(maxWidth,(uint64_t)spectrum.size());sumWidth+=weight*spectrum.size();
  auto [its,news]=supports.emplace(spectrum,code);if(!news)collisions.push_back({its->second,code,(uint64_t)H});
  if(merged%10000==0)cerr<<merged<<" merged classes\n";
 }
 need(ordinary==191536&&merged==97144&&mass==(1ULL<<28),"complete input/mass controls");
 cout<<"{\"n\":9,\"ordinary\":"<<ordinary<<",\"merged\":"<<merged<<",\"path_presentations\":"<<mass<<",\"H_distinct\":"<<hvalues.size()<<",\"tc_distinct\":"<<tcvalues.size()<<",\"extension_support_distinct\":"<<supports.size()<<",\"extension_distribution_distinct\":"<<distributions.size()<<",\"unit_sources\":"<<unitSources<<",\"unit_selections_sum\":"<<unitSets<<",\"unit_orbits_sum\":"<<unitOrbits<<",\"unit_odd_coefficients\":"<<coefficients.size()<<",\"extension_width_max\":"<<maxWidth<<",\"extension_width_sum\":"<<sumWidth;
 auto output=[&](const char*key,const vector<array<uint64_t,3>>&v){cout<<",\""<<key<<"\":[";for(size_t i=0;i<v.size();i++){if(i)cout<<',';cout<<'['<<v[i][0]<<','<<v[i][1]<<','<<v[i][2]<<']';}cout<<']';};
 output("support_collisions",collisions);output("distribution_collisions",distributionCollisions);cout<<",\"unit_coefficients\":[";bool first=true;for(I c:coefficients){if(!first)cout<<',';first=false;cout<<c;}cout<<"]}\n";
}
