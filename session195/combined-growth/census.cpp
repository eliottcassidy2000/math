// Exact target-bounded arbitrary-depth census. Standard C++17 only.
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>
using namespace std;
int main(int argc,char**argv){try{
 if(argc!=4)throw runtime_error("usage: census inventory.txt bound combined|near");
 int N=stoi(argv[2]);string mode=argv[3];
 if(N<1||N>65535)throw runtime_error("bound must be 1..65535 (16-bit state packing)");
 if(mode!="combined"&&mode!="near")throw runtime_error("unknown mode");
 vector<array<int64_t,8>> bank;ifstream f(argv[1]);if(!f)throw runtime_error("cannot open inventory");
 array<int64_t,8> y;while(f>>y[0]>>y[1]>>y[2]>>y[3]>>y[4]>>y[5]>>y[6]>>y[7])bank.push_back(y);
 if(!f.eof()||bank.empty())throw runtime_error("invalid/empty inventory");
 auto pack=[](int64_t h,int64_t e,int64_t r,int64_t d){return (uint64_t)h|((uint64_t)e<<16)|((uint64_t)r<<32)|((uint64_t)d<<48);};
 unordered_set<uint64_t> seen;vector<uint64_t> q;vector<bool> vals(N+1);
 auto seed=pack(1,1,1,1);seen.insert(seed);q.push_back(seed);vals[1]=true;
 for(size_t i=0;i<q.size();i++){
  uint64_t z=q[i];int64_t h=z&65535,e=(z>>16)&65535,r=(z>>32)&65535,d=(z>>48)&65535;
  auto add=[&](int64_t H,int64_t E,int64_t R,int64_t D){
   if(!(0<=E&&E<=D&&D<=H&&H<=R&&R+E>=H+D&&H>=h&&R>=r&&H%2==1))throw runtime_error("physical cone/parity/monotonicity failed");
   if(H>N)return;vals[H]=true;if(R>N)return;
   uint64_t v=pack(H,E,R,D);if(seen.insert(v).second)q.push_back(v);
  };
  add(h+r-e,h-e,h+r,h);add(h+r-e,d,2*r-d,d);add(h,0,2*r,d);
  if(mode=="combined")for(auto y:bank)add(y[0]*h-y[1]*e+y[2]*r,y[3]*h-y[6]*e,y[4]*h+y[5]*r,y[7]*h);
 }
 size_t nv=0;for(bool v:vals)nv+=v;
 cout<<"{\n  \"status\": \"PASS\",\n  \"mode\": \""<<mode<<"\",\n  \"target_bound\": "<<N<<",\n  \"gadget_matrices\": "<<(mode=="combined"?bank.size():1)<<",\n  \"retained_states\": "<<q.size()<<",\n  \"attained_odd_counts\": "<<nv<<",\n  \"missing_odd_counts\": [";
 bool first=true;for(int h=1;h<=N;h+=2)if(!vals[h]){if(!first)cout<<", ";first=false;cout<<h;}
 cout<<"],\n  \"scope\": \"Complete arbitrary-depth grammar support at the stated target bound; not all tournaments.\"\n}\n";
}catch(const exception&e){cerr<<e.what()<<"\n";return 1;}}
