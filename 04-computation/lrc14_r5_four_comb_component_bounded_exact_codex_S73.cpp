// Exact bounded certificate for THM-1123.
//
// For each eight-speed core P in {1,...,12}, let lo=13*max(P)+1 and scan
// every ordered quadruple in the forty-speed window [lo,lo+39].  Construct
//
//   S(P) \ (D_k1 union D_k2 union D_k3 union D_k4)
//
// by exact endpoint subtraction and decide 7*k4*L>1 by integer cross-
// multiplication.  This is a bounded bank, not an all-scale r=5 theorem.
//
// Tournament/alternate-carrier audit: runner, comb, core-component, tooth,
// endpoint, residue, and proof-obligation vertices were considered.  The
// exact endpoint word is the useful carrier.  Pairwise rational order gives
// a transitive tournament (ties are exact coincidences; coalesce them before
// taking the unique sorted Hamiltonian path).  Order alone destroys metric
// gap lengths and endpoint owners, so the interval word retains coordinate
// and owner sidecars.  Scalar half-coverage and comb-overlap MST carriers are
// audited independently in the Fraction replay and are not faithful here.

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

using i64 = std::int64_t;
using i128 = __int128_t;

struct Rat { i64 n; i64 d; };

static int compare(Rat a, Rat b) {
  const i128 x = static_cast<i128>(a.n) * b.d;
  const i128 y = static_cast<i128>(b.n) * a.d;
  return (x > y) - (x < y);
}
static Rat rmin(Rat a, Rat b) { return compare(a,b)<=0?a:b; }
static Rat rmax(Rat a, Rat b) { return compare(a,b)>=0?a:b; }
static Rat difference(Rat b, Rat a) {
  return {b.n*a.d-a.n*b.d,b.d*a.d};
}
static Rat reduced(Rat x) {
  if (x.d<0) { x.n=-x.n; x.d=-x.d; }
  const i64 g=std::gcd(std::llabs(x.n),x.d);
  return {x.n/g,x.d/g};
}
static std::string show(Rat x) {
  x=reduced(x); return std::to_string(x.n)+"/"+std::to_string(x.d);
}

using Interval=std::pair<Rat,Rat>;
using Region=std::vector<Interval>;

static Region remove_bad(const Region& input,int k) {
  Region output; output.reserve(input.size()+32);
  const i64 den=14LL*k;
  for (const auto& [a,b]:input) {
    int lo=static_cast<int>((static_cast<i128>(a.n)*k)/a.d)-1;
    int hi=static_cast<int>((static_cast<i128>(b.n)*k)/b.d)+1;
    lo=std::max(lo,0); hi=std::min(hi,k);
    Rat cursor=a;
    for (int j=lo;j<=hi;++j) {
      Rat left{14LL*j-1,den},right{14LL*j+1,den};
      if (compare(right,a)<=0 || compare(left,b)>=0) continue;
      left=rmax(left,a); right=rmin(right,b);
      if (compare(left,cursor)>0) output.push_back({cursor,left});
      cursor=rmax(cursor,right);
    }
    if (compare(cursor,b)<0) output.push_back({cursor,b});
  }
  return output;
}

static Region safe_set(const std::vector<int>& core) {
  Region r{{Rat{0,1},Rat{1,1}}};
  for (int p:core) r=remove_bad(r,p);
  return r;
}
static Rat longest(const Region& r) {
  Rat best{0,1};
  for (const auto& [a,b]:r) {
    Rat x=difference(b,a); if (compare(x,best)>0) best=x;
  }
  return best;
}
static std::vector<Interval> maximizers(const Region& r) {
  const Rat L=longest(r); std::vector<Interval> out;
  for (const auto& x:r)
    if (compare(difference(x.second,x.first),L)==0) out.push_back(x);
  return out;
}
static std::string show_core(const std::vector<int>& p) {
  std::ostringstream s;s<<"[";
  for (std::size_t i=0;i<p.size();++i){if(i)s<<",";s<<p[i];}
  return s.str()+"]";
}

struct Candidate {
  bool set=false; Rat L{0,1}; int a=0,b=0,c=0,d=0;
};
static int compare_metric(const Candidate& x,const Candidate& y) {
  const i128 lhs=static_cast<i128>(7)*x.d*x.L.n*y.L.d;
  const i128 rhs=static_cast<i128>(7)*y.d*y.L.n*x.L.d;
  return (lhs>rhs)-(lhs<rhs);
}
static bool harder(const Candidate& x,const Candidate& y) {
  if(!x.set)return false;if(!y.set)return true;
  int q=compare_metric(x,y);if(q)return q<0;
  return std::tie(x.a,x.b,x.c,x.d)<std::tie(y.a,y.b,y.c,y.d);
}

struct CoreResult {
  std::vector<int> core; Rat ell{0,1}; std::uint64_t rows=0,failures=0;
  Candidate hardest;
};

static CoreResult scan(const std::vector<int>& core) {
  CoreResult z;z.core=core;Region r0=safe_set(core);z.ell=longest(r0);
  const int lo=13*core.back()+1,hi=lo+39;
  for(int a=lo;a<=hi-3;++a){Region r1=remove_bad(r0,a);
   for(int b=a+1;b<=hi-2;++b){Region r2=remove_bad(r1,b);
    for(int c=b+1;c<=hi-1;++c){Region r3=remove_bad(r2,c);
     for(int d=c+1;d<=hi;++d){
      ++z.rows;Region r4=remove_bad(r3,d);Rat L=longest(r4);
      Candidate x{true,L,a,b,c,d};
      if(static_cast<i128>(7)*d*L.n<=L.d)++z.failures;
      if(harder(x,z.hardest))z.hardest=x;
  }}}}
  return z;
}

int main(){
  std::vector<std::vector<int>> cores;
  for(int mask=0;mask<(1<<12);++mask){
    if(__builtin_popcount(static_cast<unsigned>(mask))!=8)continue;
    std::vector<int> p;for(int i=0;i<12;++i)if((mask>>i)&1)p.push_back(i+1);
    cores.push_back(p);
  }
  std::vector<CoreResult> rows(cores.size());std::atomic<std::size_t> next{0};
  const unsigned detected=std::thread::hardware_concurrency();
  const unsigned workers=std::max(1u,std::min(8u,detected?detected:8u));
  std::vector<std::thread> pool;
  for(unsigned w=0;w<workers;++w)pool.emplace_back([&]{for(;;){
    std::size_t i=next.fetch_add(1);if(i>=cores.size())return;rows[i]=scan(cores[i]);
  }});
  for(auto& t:pool)t.join();

  std::uint64_t total=0,failures=0;Candidate hardest;std::size_t hard_core=0;
  Rat minell{1,1};std::vector<std::vector<int>> mincores;
  std::cout<<"THM-1123 bounded exact sharp four-comb certificate\n";
  std::cout<<"arithmetic=integer-rational; decisions=__int128 cross-products\n";
  std::cout<<"scope=all 495 eight-speed cores; first 40 legal speeds per core\n";
  std::cout<<"window=[13*max(P)+1,13*max(P)+40]\n";
  std::cout<<"target=7*k4*L>1\ncores="<<rows.size()<<"\n";
  for(std::size_t i=0;i<rows.size();++i){const auto& r=rows[i];
    total+=r.rows;failures+=r.failures;
    if(compare(r.ell,minell)<0){minell=r.ell;mincores={r.core};}
    else if(compare(r.ell,minell)==0)mincores.push_back(r.core);
    if(harder(r.hardest,hardest)){hardest=r.hardest;hard_core=i;}
  }
  const auto& hc=rows[hard_core];Region rr=safe_set(hc.core);
  rr=remove_bad(rr,hardest.a);rr=remove_bad(rr,hardest.b);
  rr=remove_bad(rr,hardest.c);rr=remove_bad(rr,hardest.d);
  Rat metric{7LL*hardest.d*hardest.L.n,hardest.L.d};
  Rat ratio{metric.d,metric.n};
  std::cout<<"rows_per_core=91390\n";
  std::cout<<"rows_total="<<total<<"\nfailures_total="<<failures<<"\n";
  std::cout<<"min_core_ell="<<show(minell)<<"\n";
  for(const auto& p:mincores)std::cout<<"min_ell_core="<<show_core(p)<<"\n";
  std::cout<<"hardest_core="<<show_core(hc.core)<<"\n";
  std::cout<<"hardest_quadruple=("<<hardest.a<<","<<hardest.b<<","<<hardest.c<<","<<hardest.d<<")\n";
  std::cout<<"hardest_longest_component="<<show(hardest.L)<<"\n";
  std::cout<<"min_7k4L="<<show(metric)<<"\nmax_sharp_R="<<show(ratio)<<"\n";
  for(const auto& [a,b]:maximizers(rr))
    std::cout<<"hardest_longest_interval=["<<show(a)<<","<<show(b)<<"]\n";
  std::cout<<"tournament_vertices=exact surviving endpoints\n";
  std::cout<<"pairwise_observable=exact rational order; tournament=transitive\n";
  std::cout<<"order_only_destroyed=metric lengths|endpoint owners|removal stages\n";
  std::cout<<"challenged_vertices=runners|combs|core components|teeth|endpoints|residues|proof obligations\n";
  std::cout<<"scope_warning=bounded bank only; uniform r=5 remains open\n";
  std::cout<<"certificate="<<(failures==0?"PASS":"FAIL")<<"\n";
  return failures==0?0:1;
}
