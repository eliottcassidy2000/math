// Deterministic bounded witness search; absence is never a spectrum exclusion.
// g++ -O3 -std=c++17 search.cpp -o search.exe
// search.exe bases.txt 20000 > search.jsonl
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <vector>
using I=long long;
static I en[1024][10],st[1024][10],hs[1024],ms[1024],ws[1024];
static const std::set<I> targets={245,853,909,1049,1407,1467,1515,1577,1629,1829,1887,2023,2037,2085,2343,2393,2703,3001,3033,3123,3313,3411,3427,3545,3713,3733,3753,3773,3907,3909,3985,3991,4023};
int main(int argc,char**argv){
 if(argc!=3 && argc!=4)return 2;
 bool exhaustive=argc==4 && std::string(argv[3])=="exhaustive";
 std::ifstream f(argv[1]); std::vector<std::array<int,9>> bases;
 std::array<int,9>b;while(f>>b[0]){for(int j=1;j<9;j++)f>>b[j];bases.push_back(b);}
 if(bases.empty())return 3;
 std::mt19937_64 rng(20261110);int total=exhaustive?512*bases.size():std::stoi(argv[2]);std::set<I>found;
 I units=0;
 for(int it=0;it<total;it++){
  int out[10]={},inn[10]={};
  // Alternate uniform labelled bases and a random extension of an inherited base.
  if(exhaustive || it%2){auto base=bases[exhaustive?it/512:rng()%bases.size()];for(int i=0;i<9;i++)out[i]=base[i];
    auto cut=exhaustive?it%512:rng()&511;for(int i=0;i<9;i++)if(cut>>i&1)out[9]|=1<<i;else out[i]|=512;
  }else{for(int i=0;i<10;i++)for(int j=i+1;j<10;j++)if(rng()&1)out[i]|=1<<j;else out[j]|=1<<i;}
  for(int i=0;i<10;i++)inn[i]=1023^out[i]^(1<<i);
  for(int s=0;s<1024;s++){std::fill(en[s],en[s]+10,0);std::fill(st[s],st[s]+10,0);}
  for(int v=0;v<10;v++)en[1<<v][v]=st[1<<v][v]=1;
  for(int s=1;s<1024;s++)for(int v=0;v<10;v++){
    if(en[s][v])for(int nx=out[v]&~s;nx;nx&=nx-1){int u=__builtin_ctz(nx);en[s|1<<u][u]+=en[s][v];}
    if(st[s][v])for(int nx=inn[v]&~s;nx;nx&=nx-1){int u=__builtin_ctz(nx);st[s|1<<u][u]+=st[s][v];}
  }
  hs[0]=1;for(int s=1;s<1024;s++){hs[s]=0;for(int v=0;v<10;v++)hs[s]+=en[s][v];}
  I U[10]={},W[10][10]={};
  for(int s=1;s<1024;s++)for(int v=0;v<10;v++)U[v]+=st[s][v]*hs[1023^s];
  for(int s=1;s<1023;s++)for(int v=0;v<10;v++)if(en[s][v])for(int w=0;w<10;w++)if(st[1023^s][w]){
   I z=en[s][v]*st[1023^s][w];W[v][w]+=z;W[w][v]+=z;
  }
  ms[0]=ws[0]=0;
  // Every selected set; no cardinality filter is inherited from the old search.
  for(int O=1;O<1024;O++){
   int v=__builtin_ctz(O),R=O^(1<<v);ms[O]=ms[R]+U[v];ws[O]=ws[R];
   for(int q=R;q;q&=q-1)ws[O]+=W[v][__builtin_ctz(q)];
   if(hs[1023]-ws[O]!=1)continue;
   units++;I c=ms[O];while(c&&!(c&1))c/=2;
   if(!targets.count(c)||found.count(c))continue;
   found.insert(c);std::cout<<"{\"mode\":\""<<(exhaustive?"exhaustive":"sample")<<"\",\"iteration\":"<<it<<",\"c\":"<<c<<",\"h\":"<<hs[1023]<<",\"M\":"<<ms[O]<<",\"O\":"<<O<<",\"rows\":[";
   for(int i=0;i<10;i++){if(i)std::cout<<',';std::cout<<out[i];}std::cout<<"]}\n"<<std::flush;
  }
 }
 std::cerr<<"mode="<<(exhaustive?"exhaustive":"sample")<<" graphs="<<total<<" unit_sets="<<units<<" targets_found="<<found.size()<<" seed=20261110 inherited_bases="<<bases.size()<<"\n";
}
