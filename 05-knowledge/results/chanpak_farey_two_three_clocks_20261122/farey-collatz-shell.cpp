// Search all positive flank pairs linked by <= d steps of a signed shortened Collatz map.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
using U=__uint128_t; using u=uint64_t;
struct Chart{u h,q,c,t,e,r,k;};
constexpr Chart charts[]={{2465,2293760,1344270,712944,1740,986,423385},
                          {37265,144834560,84839910,57602544,24360,14906,33838015}};
U f(U x){U z=0;for(;x;x>>=1)if(!(x&1))z+=(x>>1)+((x>>1)&1);return z;}
U P(Chart a,U x,U y){return a.k+a.t*x+a.c*y+a.q*x*y;}
U H(Chart a,U x,U y){return P(a,x,y)+a.r*(x&1)+a.e*(y&1)+a.h*(f(x)+f(y));}
struct Stats{u pieces=0,phase=0,band=0;std::vector<std::pair<u,u>> hits;};
class Search{
  Chart a; U N; int sgn,limit; bool swap; public:std::vector<Stats> stats;
  Search(Chart a,U N,int sgn,int limit,bool sw):a(a),N(N),sgn(sgn),limit(limit),swap(sw),stats(limit+1){}
  void test(int d,u A,u B,u C,u D) {
    // x=A*Z+B, y=C*Z+D, possibly after swapping physically.
    if(swap){std::swap(A,C);std::swap(B,D);}
    U Q=U(a.q)*A*C, L=U(a.q)*(U(A)*D+U(B)*C)+U(a.t)*A+U(a.c)*C;
    U K=P(a,B,D);
    if(K>N)return;
    U M=N-K;
    u Z=2*(double)M/(sqrt((double)L*(double)L+4*(double)Q*(double)M)+(double)L);
    auto le=[&](u z){return z==0||Q*z+L<=M/z;}; // exact, no product of three terms
    // Floating point is only a seed. Below the anchor cutoff, all addresses
    // are small enough that even after adding one to Z no quantity wraps.
    if(!le(Z)) {
      if(Z&&le(Z-1))--Z;
      else {u lo=0,hi=Z;while(hi-lo>1){u m=lo+(hi-lo)/2;if(le(m))lo=m;else hi=m;}Z=lo;}
    } else if(le(Z+1)) {
      if(!le(Z+2))++Z;
      else {u lo=Z+1,hi=2*Z+4;while(le(hi))hi*=2;
            while(hi-lo>1){u m=lo+(hi-lo)/2;if(le(m))lo=m;else hi=m;}Z=lo;}
    }
    U x=U(A)*Z+B,y=U(C)*Z+D;
    if(!x||!y)return;
    ++stats[d].pieces;
    U M0=M-Z*(Q*Z+L), parity=a.r*(x&1)+a.e*(y&1);
    if(M0<parity||M0-parity>a.h*(x+y))return;
    ++stats[d].band;
    M0-=parity;
    if(M0%(2*a.h))return;
    ++stats[d].phase;
    if(M0==a.h*(f(x)+f(y)))stats[d].hits.emplace_back(u(x),u(y));
  }
  void step(int d,u B,u C,u D) {
    // x=2^d Z+B and T_sgn^d(x)=C Z+D, T(n)=n/2 (even) or (3n+sgn)/2 (odd).
    if(d)test(d,1ULL<<d,B,C,D);
    if(d==limit)return;
    for(u b=0;b<2;++b){
      u E=D+C*b; bool odd=E&1;
      step(d+1,B+(b<<d),C*(odd?3:1),odd?(sgn>0?3*E+1:3*E-1)/2:E/2);
    }
  }
  void run(){step(0,0,1,0);}
};
int main(int argc,char**argv) {
 const U N=((U(1)<<60)*((U(1)<<61)-1))-1;
 int maxdepth=argc>1?std::stoi(argv[1]):23;
 if(maxdepth<1||maxdepth>27)throw std::runtime_error("depth range 1..27");
 if(argc>2&&std::string(argv[2])=="check"){
   for(auto a:charts)for(int sign:{-1,1})for(bool sw:{false,true}){
     for(u x:{1,2,27,153,999,1024,76543,999999,43626901}){
       u y=x;
       for(int d=1;d<=8;++d){
         y=(y&1)?(sign==1?3*y+1:3*y-1)/2:y/2;
         if(!y)continue;
         u p=sw?y:x,q=sw?x:y;
         Search se(a,H(a,p,q),sign,d,sw);se.run();
         auto h=se.stats[d].hits;
         if(std::find(h.begin(),h.end(),std::make_pair(p,q))==h.end()){std::cerr<<a.h<<" "<<sign<<" "<<sw<<" "<<d<<" "<<p<<" "<<q<<"\n";throw std::runtime_error("lost planted pair");}
       }
     }
   }
   for(u x=0;x<=2048;++x)if(f(x)>x||(f(x)&1))throw std::runtime_error("f control");
   std::cout<<"planted signed/oriented controls PASS\n";return 0;
 }
 for(int j:{0,1})for(int sign:{-1,1})for(bool sw:{false,true}){
   Search se(charts[j],N,sign,maxdepth,sw);se.run();
   u pieces=0,band=0,phase=0,hits=0;
   for(int d=1;d<=maxdepth;++d){
     auto &v=se.stats[d];pieces+=v.pieces;band+=v.band;phase+=v.phase;hits+=v.hits.size();
     for(auto [x,y]:v.hits)std::cout<<"HIT depth="<<d<<" x="<<x<<" y="<<y<<"\n";
   }
   std::cout<<"spine="<<4+4*j<<" sign="<<sign<<" swapped="<<sw<<" depth="<<maxdepth
     <<" pieces="<<pieces<<" band="<<band<<" phase="<<phase<<" hits="<<hits<<std::endl;
 }
}
