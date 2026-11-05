#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <unordered_set>
#include <vector>
using namespace std;
int main(int argc,char **argv){
 if(argc!=3)throw runtime_error("usage: census bank.txt bound");
 int N=stoi(argv[2]);if(N<1||N>65535)throw runtime_error("bound must be in 1..65535");
 vector<array<int64_t,8>>bank;ifstream f(argv[1]);if(!f)throw runtime_error("cannot open bank");array<int64_t,8>y;
 while(f>>y[0]>>y[1]>>y[2]>>y[3]>>y[4]>>y[5]>>y[6]>>y[7])bank.push_back(y);
 auto pack=[](int64_t h,int64_t e,int64_t r,int64_t d){return uint64_t(h)|(uint64_t(e)<<16)|(uint64_t(r)<<32)|(uint64_t(d)<<48);};
 unordered_set<uint64_t>seen;vector<uint64_t>q;vector<bool>vals(N+1);seen.insert(pack(1,1,1,1));q.push_back(pack(1,1,1,1));vals[1]=true;
 for(size_t i=0;i<q.size();i++){
  uint64_t s=q[i];int64_t h=s&65535,e=(s>>16)&65535,r=(s>>32)&65535,d=(s>>48)&65535;
  auto add=[&](int64_t H,int64_t E,int64_t R,int64_t D){
   if(H<h||E<0||E>H||R<H||D<1||D>H)throw runtime_error("invalid physical invariant");
   if(H>N)return;
   vals[H]=true;
   if(R>N)return;
   uint64_t t=pack(H,E,R,D);if(seen.insert(t).second)q.push_back(t);
  };
  add(h+r-e,d,2*r-d,d);add(h,0,2*r,d);
  for(auto y:bank)add(y[0]*h-y[1]*e+y[2]*r,y[3]*h-y[6]*e,y[4]*h+y[5]*r,y[7]*h);
 }
 cout<<"{\"target_bound\":"<<N<<",\"bank_entries\":"<<bank.size()<<",\"retained_states\":"<<q.size()<<",\"missing_odds\":[";
 bool first=true;for(int h=1;h<=N;h+=2)if(!vals[h]){if(!first)cout<<",";cout<<h;first=false;}
 cout<<"]}\n";
}
