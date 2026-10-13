// Packed exact parity census. C++17; unsigned-128 overflow guard.
// Usage: binary MAX_EXP EDGES_TXT OUTPUT_JSON [EVEN_JUMP_DEPTH]
// Every source stops at a smaller source; no unproved convergence input.
// Even fixed jumps preserve total-stopping parity, including through the
// terminal 1<->2 cycle. This is a parity computation, not a strict ROOT word.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
using u64=uint64_t; using u128=__uint128_t;
void require(bool b,const char* s){if(!b)throw std::runtime_error(s);}
std::string decimal(u128 n){std::string s;do{s.push_back('0'+n%10);n/=10;}while(n);std::reverse(s.begin(),s.end());return s;}
struct Jump{u64 p,t;};
int main(int argc,char**argv){try{
 require(argc>=4,"arguments");int E=std::stoi(argv[1]),K=argc>4?std::stoi(argv[4]):12;
 require(E>=10&&E<=34&&K>=4&&K<=16&&K%2==0,"declared exponent/depth scope");
 u64 N=1ULL<<E, mask=(1ULL<<K)-1;std::vector<Jump> jumps(1ULL<<K);
 for(u64 r=0;r<=mask;r++){u64 m=r,p=1;for(int j=0;j<K;j++){if(m&1){m=(3*m+1)/2;p*=3;}else m/=2;}jumps[r]={p,m};}
 std::vector<u64> bits((N+64)/64,0);auto parity=[&](u64 n){return (bits[n>>6]>>(n&63))&1ULL;};
 auto start=std::chrono::steady_clock::now();u64 max_jumps=0,wide_jumps=0;u128 max_state=0;
 for(u64 n=2;n<=N;n++){
  u64 m=n,loops=0,sgn=0;
  if(!(n&1)){m=n/2;sgn=1;}
  else{
   do{
    Jump j=jumps[m&mask];u128 z=u128(m>>K)*j.p+j.t;
    if(z>max_state)max_state=z;
    ++loops;require(loops<1000000,"descent guard exceeded; no completeness claimed");
    if(z>std::numeric_limits<u64>::max()){
     ++wide_jumps;
     while(z>std::numeric_limits<u64>::max()){
      require(z<(u128(1)<<120),"128-bit safety bound exceeded");
      for(int t=0;t<K;t++)z=(z&1)?(3*z+1)/2:z/2;
      ++loops;require(loops<1000000,"wide descent guard exceeded");
     }
    }
    m=u64(z);
   }while(m>=n);
  }
  require(m>=1&&m<n,"parent must be an earlier positive source");
  if((parity(m)^sgn)!=0)bits[n>>6]|=1ULL<<(n&63);
  max_jumps=std::max(max_jumps,loops);
  if(n%(1ULL<<28)==0||n==N){double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();std::cerr<<"sources "<<n<<" / "<<N<<" seconds "<<sec<<"\n";}
 }
 // Independent strict first-hit recurrence, using single Terras steps.
 u64 controls=0,max_tau=0;
 auto control=[&](u64 n){u128 m=n;u64 tau=0;while(m!=1){require(m<(u128(1)<<120),"control overflow");m=(m&1)?(3*m+1)/2:m/2;++tau;require(tau<1000000,"control cap");}require(parity(n)==tau%2,"strict ROOT parity control");max_tau=std::max(max_tau,tau);++controls;};
 for(u64 n=1;n<=std::min(N,1ULL<<18);n++)control(n);
 u64 rng=0x9e3779b97f4a7c15ULL;
 for(int i=0;i<16384;i++){rng^=rng<<13;rng^=rng>>7;rng^=rng<<17;control(1+rng%N);}
 for(int e=0;e<=E;e++){control(1ULL<<e);require(parity(1ULL<<e)==u64(e%2),"power-of-two phase");}
 struct Sum{u64 count,ones,even_count,even_ones;};
 auto interval=[&](u64 lo,u64 hi){Sum s{hi-lo,0,0,0};
  for(u64 w=lo>>6;w<=((hi-1)>>6);w++){
   unsigned a=w==(lo>>6)?lo%64:0,b=w==((hi-1)>>6)?unsigned((hi-1)%64+1):64;
   u64 take=(~0ULL<<a)&(b==64?~0ULL:((1ULL<<b)-1));u64 v=bits[w]&take;
   s.ones+=__builtin_popcountll(v);s.even_count+=__builtin_popcountll(take&0x5555555555555555ULL);s.even_ones+=__builtin_popcountll(v&0x5555555555555555ULL);
  }return s;
 };
 std::ifstream in(argv[2]);require(bool(in),"edges input");std::vector<u64> edges;u64 x;
 while(in>>x)if(x<=N)edges.push_back(x);
 require(edges.size()>1&&edges.back()==N,"edges terminate at N");
 for(size_t j=1;j<edges.size();j++)require(edges[j]>edges[j-1],"strict edges");
 std::ofstream out(argv[3]);require(bool(out),"output file");
 out<<"{\n\"max_exp\":"<<E<<",\"N\":"<<N<<",\"jump_depth\":"<<K<<",\"packed_bytes\":"<<bits.size()*8<<",\"strict_controls\":"<<controls<<",\"max_control_tau\":"<<max_tau<<",\"max_jump_count\":"<<max_jumps<<",\"wide_jump_count\":"<<wide_jumps<<",\"max_jump_endpoint\":\""<<decimal(max_state)<<"\",\n\"bins\":[\n";
 for(size_t j=1;j<edges.size();j++){
  u64 lo=edges[j-1],hi=edges[j];auto s=interval(lo,hi);auto ev=interval((lo+1)/2,(hi+1)/2);
  require(s.even_count==ev.count&&s.even_ones==ev.count-ev.ones,"exact even-bin doubling identity");
  out<<"["<<lo<<","<<hi<<","<<s.count<<","<<int64_t(s.count)-2*int64_t(s.ones)<<","<<s.even_count<<","<<int64_t(s.even_count)-2*int64_t(s.even_ones)<<"]"<<(j+1==edges.size()?"\n":",\n");
 }
 auto cuts=[&](u64 lo,u64 hi){u64 count=0;
  for(u64 w=lo>>6;w<=((hi-1)>>6);w++){
   unsigned a=w==(lo>>6)?lo%64:0,b=w==((hi-1)>>6)?unsigned((hi-1)%64+1):64;
   u64 take=(~0ULL<<a)&(b==64?~0ULL:((1ULL<<b)-1));
   u64 adjacent=(bits[w]>>1)|((w+1<bits.size()?bits[w+1]:0ULL)<<63);
   count+=__builtin_popcountll((bits[w]^adjacent)&take);
  }return count;
 };
 for(u64 hi=2;hi<200;hi++){
  u64 direct=0;for(u64 j=1;j<hi;j++)direct+=(parity(j)!=parity(j+1));
  require(cuts(1,hi)==direct,"packed adjacent-cut control");
 }
 out<<"],\"dyadic_transfers\":[\n";
 bool comma=false;
 for(u64 M=1;3*M<=N;M*=2){
  auto signed_sum=[&](u64 end){auto ss=interval(1,end+1);return int64_t(ss.count)-2*int64_t(ss.ones);};
  int64_t s1=signed_sum(M),s2=signed_sum(2*M),s3=signed_sum(3*M);
  int64_t numerator=3*s2+3*s1+s3;u64 changes=cuts(1,3*M);
  require(numerator<=int64_t(2*changes)&&numerator>=-int64_t(2*changes),"exact signed-transfer variation bound");
  if(comma)out<<",\n";comma=true;
  out<<"["<<M<<","<<s1<<","<<s2<<","<<s3<<","<<changes<<","<<numerator<<"]";
 }
 auto total=interval(1,N+1);out<<"\n],\"prefix_sum\":"<<int64_t(total.count)-2*int64_t(total.ones)<<"}\n";
 std::cerr<<"complete; strict controls "<<controls<<" max tau "<<max_tau<<"\n";
 }catch(std::exception const&e){std::cerr<<"ERROR "<<e.what()<<"\n";return 1;}}
