#include <algorithm>
#include <gmpxx.h>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using cpp_int = mpz_class;

namespace {

cpp_int ab(const cpp_int& x) { return x < 0 ? -x : x; }
cpp_int gc(cpp_int a, cpp_int b) {
  a = ab(a); b = ab(b);
  while (b != 0) { cpp_int r = a % b; a = b; b = r; }
  return a;
}

struct Rat {
  cpp_int n = 0, d = 1;
  Rat() = default;
  Rat(cpp_int num, cpp_int den = 1) : n(std::move(num)), d(std::move(den)) {
    if (d == 0) throw std::runtime_error("zero denominator");
    if (d < 0) { n = -n; d = -d; }
    cpp_int g = gc(n,d); n /= g; d /= g;
  }
};
int cmp(const Rat& a, const Rat& b) {
  cpp_int x=a.n*b.d, y=b.n*a.d; return (x>y)-(x<y);
}
std::string str(const cpp_int& x) { return x.get_str(); }
std::string str(const Rat& x) { return str(x.n)+"/"+str(x.d); }
cpp_int z(uint64_t x) { return cpp_int(std::to_string(x)); }
cpp_int ceil_div(cpp_int a, const cpp_int& b) {
  if (b <= 0) throw std::runtime_error("ceil denominator");
  cpp_int q=a/b, r=a%b;
  if (r>0) ++q;
  return q;
}

struct Tournament {
  int n;
  std::vector<uint32_t> out, in;
};

void arc(Tournament& t, int u, int v) { t.out[u] |= uint32_t(1u<<v); }

Tournament construction(int r) {
  if (r<1 || !(r&1)) throw std::runtime_error("odd positive block size");
  const int n=r+4, a=0, c=r+1, d=r+2, e=r+3;
  Tournament t{n,std::vector<uint32_t>(n),std::vector<uint32_t>(n)};
  // B is the cyclic regular tournament: i -> i+1,...,i+(r-1)/2 mod r.
  for(int i=0;i<r;++i) for(int step=1;step<=(r-1)/2;++step)
    arc(t,1+i,1+(i+step)%r);
  for(int b=1;b<=r;++b) {
    arc(t,a,b); arc(t,b,c); arc(t,b,d); arc(t,b,e);
  }
  arc(t,a,c); arc(t,d,a); arc(t,a,e);
  arc(t,c,d); arc(t,e,c); arc(t,d,e);
  for(int v=0;v<n;++v) for(int u=0;u<n;++u)
    if(t.out[u]&(1u<<v))t.in[v]|=uint32_t(1u<<u);
  for(int u=0;u<n;++u) for(int v=u+1;v<n;++v)
    if(bool(t.out[u]&(1u<<v))==bool(t.out[v]&(1u<<u)))
      throw std::runtime_error("construction missing/double arc");
  return t;
}

Tournament deleted_cyclic_construction(int R) {
  if (R<3 || !(R&1)) throw std::runtime_error("odd cyclic parent size");
  const int bsize=R-1,n=bsize+4,a=0,c=bsize+1,d=bsize+2,e=bsize+3;
  Tournament t{n,std::vector<uint32_t>(n),std::vector<uint32_t>(n)};
  // Vertex b in 1,...,R-1 represents the same vertex of R_R; vertex 0 is deleted.
  for(int i=1;i<R;++i)for(int step=1;step<=(R-1)/2;++step){int j=(i+step)%R;if(j!=0)arc(t,i,j);}
  for(int b=1;b<=bsize;++b){arc(t,a,b);arc(t,b,c);arc(t,b,d);arc(t,b,e);}
  arc(t,a,c);arc(t,d,a);arc(t,a,e);arc(t,c,d);arc(t,e,c);arc(t,d,e);
  for(int v=0;v<n;++v)for(int u=0;u<n;++u)if(t.out[u]&(1u<<v))t.in[v]|=uint32_t(1u<<u);
  for(int u=0;u<n;++u)for(int v=u+1;v<n;++v)if(bool(t.out[u]&(1u<<v))==bool(t.out[v]&(1u<<u)))throw std::runtime_error("deleted construction missing/double arc");
  return t;
}

bool strong(const Tournament& t) {
  auto run=[&](bool rev){uint32_t seen=1,front=1;while(front){int v=__builtin_ctz(front);front&=front-1;uint32_t fresh=(rev?t.in[v]:t.out[v])&~seen;seen|=fresh;front|=fresh;}return __builtin_popcount(seen)==t.n;};
  return run(false)&&run(true);
}

struct Layer { int m,t; Rat J; uint64_t L,A,d; };
struct Result {
  int r,n; bool is_strong; uint64_t H,W2; cpp_int D4x4,Chdx4; Rat theta,rho;
  std::vector<Layer> layers; std::vector<int> rt,ct,at; uint64_t response_hash=0;
};

Result evaluate_tournament(Tournament t, int r, bool direct_cut_check) {
  const int n=t.n, size=1<<n, full=size-1;
  auto idx=[&](int mask,int v){return size_t(mask)*n+v;};
  std::vector<uint64_t> ending(size_t(size)*n), starting(size_t(size)*n);
  for(int v=0;v<n;++v) ending[idx(1<<v,v)]=starting[idx(1<<v,v)]=1;
  for(int mask=1;mask<size;++mask){
    uint32_t vs=mask;
    while(vs){int v=__builtin_ctz(vs);vs&=vs-1;int rest=mask^(1<<v);if(!rest)continue;
      uint32_t bits=uint32_t(rest)&t.in[v];while(bits){int u=__builtin_ctz(bits);bits&=bits-1;ending[idx(mask,v)]+=ending[idx(rest,u)];}
      bits=uint32_t(rest)&t.out[v];while(bits){int u=__builtin_ctz(bits);bits&=bits-1;starting[idx(mask,v)]+=starting[idx(rest,u)];}
    }
  }
  uint64_t H=0;for(int v=0;v<n;++v)H+=ending[idx(full,v)];
  std::vector<uint64_t> before(size_t(size)*n),after(size_t(size)*n);
  for(int v=0;v<n;++v)before[idx(0,v)]=after[idx(0,v)]=1;
  for(int mask=1;mask<size;++mask)for(int v=0;v<n;++v){
    uint32_t bits=uint32_t(mask)&t.in[v];while(bits){int u=__builtin_ctz(bits);bits&=bits-1;before[idx(mask,v)]+=ending[idx(mask,u)];}
    bits=uint32_t(mask)&t.out[v];while(bits){int u=__builtin_ctz(bits);bits&=bits-1;after[idx(mask,v)]+=starting[idx(mask,u)];}
  }
  std::vector<uint64_t> cap(size_t(n)*n), edge(size_t(n)*n), outcap(n);
  uint64_t W2=0;
  for(int x=0;x<n;++x)for(int y=x+1;y<n;++y){int u=x,v=y;if(!(t.out[u]&(1u<<v)))std::swap(u,v);int rem=full^(1<<u)^(1<<v);unsigned __int128 total=0;
    for(int left=rem;;left=(left-1)&rem){int right=rem^left;total+=(unsigned __int128)before[idx(left,u)]*after[idx(right,v)];total+=(unsigned __int128)before[idx(left,v)]*after[idx(right,u)];if(!left)break;}
    if(total>std::numeric_limits<uint64_t>::max())throw std::runtime_error("capacity overflow");uint64_t q=(uint64_t)total;cap[size_t(u)*n+v]=q;edge[size_t(u)*n+v]=edge[size_t(v)*n+u]=q;outcap[u]+=q;W2+=q;
  }
  std::vector<cpp_int>d2(n),h2(n);
  for(int v=0;v<n;++v)for(int u=0;u<n;++u)if(u!=v){cpp_int q=z(edge[size_t(v)*n+u]);d2[v]+=q;if(t.out[v]&(1u<<u))h2[v]+=q;else h2[v]-=q;}
  cpp_int Chdx4=0;for(int v=0;v<n;++v)Chdx4+=h2[v]*d2[v];
  cpp_int D4x4=0;
  for(int a=0;a<n;++a)for(int b=a+1;b<n;++b)for(int c=a+1;c<n;++c)if(c!=b)for(int d=c+1;d<n;++d)if(d!=b)D4x4+=z(edge[size_t(a)*n+b])*z(edge[size_t(c)*n+d]);
  Rat theta(cpp_int(n-3)*Chdx4,2*D4x4);
  Rat rho(cpp_int(n-3)*ab(Chdx4),(n%2?4:2)*D4x4);

  std::vector<uint64_t> es(size_t(n)*size);
  for(int v=0;v<n;++v)for(int mask=1;mask<size;++mask){int u=__builtin_ctz(mask),rest=mask^(1<<u);es[size_t(v)*size+mask]=es[size_t(v)*size+rest]+edge[size_t(v)*n+u];}
  std::vector<int64_t> cut(size);std::vector<uint64_t> values(size);values[0]=H;
  for(int mask=1;mask<size;++mask){int v=__builtin_ctz(mask),rest=mask^(1<<v);cut[mask]=cut[rest]+int64_t(outcap[v])-int64_t(es[size_t(v)*size+rest]);if(cut[mask]<0)throw std::runtime_error("negative cut");values[mask]=H+uint64_t(cut[mask]);if(values[mask]<H||!(values[mask]&1))throw std::runtime_error("support/parity");}
  if(values[full]!=H)throw std::runtime_error("full constant cut");
  if(direct_cut_check){for(int mask=0;mask<size;++mask){uint64_t q=H;for(int u=0;u<n;++u)if(mask&(1<<u))for(int v=0;v<n;++v)if(!(mask&(1<<v)))q+=cap[size_t(u)*n+v];if(q!=values[mask])throw std::runtime_error("direct cut replay");}}

  std::vector<Layer> layers; Rat bestJ(cpp_int(-1)<<100); uint64_t bestL=0,bestA=0;std::vector<int>rt,ct,at;
  for(int m=1;m<n;++m){uint64_t count=0,anchor=0,A=0,d=0;bool first=true;cpp_int sum=0,sumsq=0;for(int mask=0;mask<size;++mask)if(__builtin_popcount(mask)==m){uint64_t x=values[mask];if(first){anchor=x;first=false;}++count;cpp_int zx=z(x);sum+=zx;sumsq+=zx*zx;A=std::max(A,x);d=std::gcd(d,x>anchor?x-anchor:anchor-x);}
    Rat J(sumsq-sum*z(H),sum-z(count)*z(H));cpp_int k=d?ceil_div(J.n-z(anchor)*J.d,z(d)*J.d):0;cpp_int bigL=z(anchor)+z(d)*k;cpp_int umax(std::to_string(std::numeric_limits<uint64_t>::max()));if(bigL<0||bigL>umax)throw std::runtime_error("coset floor range");uint64_t L=bigL.get_ui();if(L>A)throw std::runtime_error("coset floor support");int tt=n-2*m;layers.push_back({m,tt,J,L,A,d});int cj=cmp(J,bestJ);if(cj>0){bestJ=J;rt={tt};}else if(cj==0)rt.push_back(tt);if(L>bestL){bestL=L;ct={tt};}else if(L==bestL)ct.push_back(tt);if(A>bestA){bestA=A;at={tt};}else if(A==bestA)at.push_back(tt);
  }
  // FNV-1a over mask-ordered uint64 responses, a compact replay checksum.
  uint64_t vh=1469598103934665603ULL;for(uint64_t x:values)for(int i=0;i<8;++i){vh^=(x>>(8*i))&255;vh*=1099511628211ULL;}
  return {r,n,strong(t),H,W2,D4x4,Chdx4,theta,rho,std::move(layers),std::move(rt),std::move(ct),std::move(at),vh};
}

Result evaluate(int r, bool direct_cut_check) {
  return evaluate_tournament(construction(r),r,direct_cut_check);
}

std::string tuple(const std::vector<int>& xs){std::ostringstream o;o<<'(';for(size_t i=0;i<xs.size();++i){if(i)o<<',';o<<xs[i];}o<<')';return o.str();}

} // namespace

int main(){try{
  int first_rational=0,first_coset=0;
  for(int r:{1,3,5,7,9,11,13}){
    Result q=evaluate(r,r==13);bool rc=std::all_of(q.rt.begin(),q.rt.end(),[](int t){return t==-1||t==1;});bool cc=std::all_of(q.ct.begin(),q.ct.end(),[](int t){return t==-1||t==1;});if(!rc&&!first_rational)first_rational=r;if(!cc&&!first_coset)first_coset=r;
    std::cout<<"r="<<r<<",n="<<q.n<<",strong="<<(q.is_strong?"yes":"no")<<",H="<<q.H<<",W="<<str(Rat(z(q.W2),2))<<",D4="<<str(Rat(q.D4x4,4))<<",Chd="<<str(Rat(q.Chdx4,4))<<",theta="<<str(q.theta)<<",rho="<<str(q.rho)<<",rational_t="<<tuple(q.rt)<<",coset_t="<<tuple(q.ct)<<",actual_t="<<tuple(q.at)<<",response_fnv64="<<std::hex<<q.response_hash<<std::dec<<"\n";
    if(r==13){
      if(!q.is_strong||q.H!=264761757ULL||cmp(Rat(z(q.W2),2),Rat(cpp_int("7641744513")))||cmp(Rat(q.D4x4,4),Rat(cpp_int("20351251042445697517")))||cmp(Rat(q.Chdx4,4),Rat(cpp_int("-6742645661951046276")))||cmp(q.rho,Rat(cpp_int("165029788928871762"),cpp_int("142316440856263619")))||cmp(q.theta,Rat(cpp_int("-330059577857743524"),cpp_int("142316440856263619")))||q.rt!=std::vector<int>{-3}||q.ct!=std::vector<int>{-3})throw std::runtime_error("n17 frozen packet mismatch");
      std::cout<<"n17_layers=(m,t,J,L,A,lattice)\n";for(const Layer&x:q.layers)std::cout<<x.m<<","<<x.t<<","<<str(x.J)<<","<<x.L<<","<<x.A<<","<<x.d<<"\n";
      std::cout<<"n17_direct_cut_recurrence_replay=PASS\n";
    }
  }
  std::cout<<"first_rational_noncentral_block_r="<<first_rational<<"\n"<<"first_coset_noncentral_block_r="<<first_coset<<"\n";
  if(first_rational!=13||first_coset!=13)throw std::runtime_error("unexpected first family failure");

  int first_deleted_rational=0,first_deleted_coset=0;
  for(int R:{3,5,7,9}){
    Tournament t=deleted_cyclic_construction(R);
    if(R==9){
      const std::vector<uint32_t> frozen={3070,3644,3704,3824,4064,4032,3970,3846,3598,1024,2049,512};
      if(t.out!=frozen)throw std::runtime_error("n12 frozen adjacency mismatch");
    }
    Result q=evaluate_tournament(std::move(t),R,R==9);
    bool rc=std::all_of(q.rt.begin(),q.rt.end(),[](int x){return x==0;});
    bool cc=std::all_of(q.ct.begin(),q.ct.end(),[](int x){return x==0;});
    if(!rc&&!first_deleted_rational)first_deleted_rational=R;
    if(!cc&&!first_deleted_coset)first_deleted_coset=R;
    std::cout<<"deleted_R="<<R<<",n="<<q.n<<",strong="<<(q.is_strong?"yes":"no")<<",H="<<q.H<<",W="<<str(Rat(z(q.W2),2))<<",D4="<<str(Rat(q.D4x4,4))<<",Chd="<<str(Rat(q.Chdx4,4))<<",theta="<<str(q.theta)<<",rho="<<str(q.rho)<<",rational_t="<<tuple(q.rt)<<",coset_t="<<tuple(q.ct)<<",actual_t="<<tuple(q.at)<<",response_fnv64="<<std::hex<<q.response_hash<<std::dec<<"\n";
    if(R==9){
      if(!q.is_strong||q.H!=27759||cmp(Rat(z(q.W2),2),Rat(cpp_int("506085")))||cmp(Rat(q.D4x4,4),Rat(cpp_int("80871049732")))||cmp(Rat(q.Chdx4,4),Rat(cpp_int("-23596773036")))||cmp(q.rho,Rat(cpp_int("53092739331"),cpp_int("40435524866")))||q.rt!=std::vector<int>{-2}||q.ct!=std::vector<int>{-2}||q.at!=std::vector<int>{-6})throw std::runtime_error("n12 frozen packet mismatch");
      int64_t central=0,outer=0;for(const Layer&x:q.layers){if(x.t==0)central=x.L;if(x.t==-2)outer=x.L;}
      if(central!=350727||outer!=352951||central-outer!=-2224)throw std::runtime_error("n12 central/coset margin mismatch");
      std::cout<<"n12_frozen_adjacency=PASS\n"<<"n12_layers=(m,t,J,L,A,lattice)\n";for(const Layer&x:q.layers)std::cout<<x.m<<","<<x.t<<","<<str(x.J)<<","<<x.L<<","<<x.A<<","<<x.d<<"\n";
      std::cout<<"n12_direct_cut_recurrence_replay=PASS\n";
    }
  }
  std::cout<<"first_deleted_rational_noncentral_parent_R="<<first_deleted_rational<<"\n"<<"first_deleted_coset_noncentral_parent_R="<<first_deleted_coset<<"\n";
  if(first_deleted_rational!=9||first_deleted_coset!=9)throw std::runtime_error("unexpected first deleted-family failure");
  std::cout<<"status=ACCEPT\n";
}catch(const std::exception&e){std::cerr<<"status=REJECT error="<<e.what()<<"\n";return 1;}}
