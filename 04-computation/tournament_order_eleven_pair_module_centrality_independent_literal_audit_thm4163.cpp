#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {
constexpr int N=11, SZ=1<<N, FULL=SZ-1, E=N*(N-1)/2;
struct P { bool strong=false; int64_t H=0,W2=0,D4x4=0,Chdx4=0; double rho=0; int64_t coset_margin=std::numeric_limits<int64_t>::max(); uint16_t inner_bad_mask=0; int64_t max_inner_lattice=0; std::array<int64_t,N+1> layer_lattice{}; std::array<int64_t,N+1> layer_floor{}; };
static int64_t gcd64(int64_t a,int64_t b){return std::gcd(std::llabs(a),std::llabs(b));}
static int64_t ceil_div(__int128 a,__int128 b){__int128 q=a/b,r=a%b;if(r>0)++q;return int64_t(q);}

static std::array<uint16_t,N> adj_from_label(uint64_t z) {
  std::array<uint16_t,N> out{}; int k=0;
  for(int i=0;i<N;i++) for(int j=i+1;j<N;j++,k++) {
    if(z&(1ull<<k)) out[i]|=uint16_t(1u<<j); else out[j]|=uint16_t(1u<<i);
  }
  return out;
}
static bool strong(const std::array<uint16_t,N>& out) {
  std::array<uint16_t,N> in{};
  for(int u=0;u<N;u++) for(int v=0;v<N;v++) if(out[u]&(1u<<v)) in[v]|=uint16_t(1u<<u);
  for(int reverse=0;reverse<2;reverse++) {
    uint16_t seen=1, todo=1;
    while(todo) { int v=__builtin_ctz(todo); todo&=uint16_t(todo-1); uint16_t q=(reverse?in[v]:out[v])&uint16_t(~seen); seen|=q; todo|=q; }
    if((seen&FULL)!=FULL) return false;
  }
  return true;
}

static P eval(uint64_t z) {
  auto out=adj_from_label(z); P p; p.strong=strong(out);
  std::array<uint16_t,N> in{};
  for(int u=0;u<N;u++) for(int v=0;v<N;v++) if(out[u]&(1u<<v)) in[v]|=uint16_t(1u<<u);
  static uint32_t ending[SZ][N], starting[SZ][N], before[SZ][N], after[SZ][N];
  std::fill(&ending[0][0],&ending[0][0]+SZ*N,0); std::fill(&starting[0][0],&starting[0][0]+SZ*N,0);
  std::fill(&before[0][0],&before[0][0]+SZ*N,0); std::fill(&after[0][0],&after[0][0]+SZ*N,0);
  for(int v=0;v<N;v++) ending[1<<v][v]=starting[1<<v][v]=1;
  for(int mask=1;mask<SZ;mask++) {
    uint16_t vs=mask;
    while(vs) { int v=__builtin_ctz(vs); vs&=uint16_t(vs-1); int rest=mask^(1<<v); if(!rest) continue;
      uint16_t q=uint16_t(rest)&in[v]; while(q){int u=__builtin_ctz(q);q&=uint16_t(q-1);ending[mask][v]+=ending[rest][u];}
      q=uint16_t(rest)&out[v]; while(q){int u=__builtin_ctz(q);q&=uint16_t(q-1);starting[mask][v]+=starting[rest][u];}
    }
  }
  for(int v=0;v<N;v++) p.H+=ending[FULL][v];
  for(int v=0;v<N;v++) before[0][v]=after[0][v]=1;
  for(int mask=1;mask<SZ;mask++) for(int v=0;v<N;v++) {
    uint16_t q=uint16_t(mask)&in[v]; while(q){int u=__builtin_ctz(q);q&=uint16_t(q-1);before[mask][v]+=ending[mask][u];}
    q=uint16_t(mask)&out[v]; while(q){int u=__builtin_ctz(q);q&=uint16_t(q-1);after[mask][v]+=starting[mask][u];}
  }
  uint32_t cap[N][N]{}; int64_t d2[N]{},h2[N]{};
  for(int a=0;a<N;a++) for(int b=a+1;b<N;b++) {
    int u=a,v=b; if(!(out[u]&(1u<<v))) std::swap(u,v); int rem=FULL^(1<<u)^(1<<v); uint64_t c=0;
    for(int L=rem;;L=(L-1)&rem){int R=rem^L;c+=uint64_t(before[L][u])*after[R][v]+uint64_t(before[L][v])*after[R][u];if(!L)break;}
    cap[u][v]=cap[v][u]=uint32_t(c); p.W2+=c;
  }
  for(int v=0;v<N;v++) for(int u=0;u<N;u++) if(u!=v){d2[v]+=cap[v][u];h2[v]+=(out[v]&(1u<<u))?cap[v][u]:-int64_t(cap[v][u]);}
  for(int v=0;v<N;v++) p.Chdx4+=h2[v]*d2[v];
  for(int a=0;a<N;a++) for(int b=a+1;b<N;b++) for(int c=a+1;c<N;c++) if(c!=b) for(int d=c+1;d<N;d++) if(d!=b) p.D4x4+=int64_t(cap[a][b])*cap[c][d];
  p.rho=2.0*std::abs(double(p.Chdx4))/double(p.D4x4);
  std::array<int64_t,SZ> values{};
  for(int mask=0;mask<SZ;mask++) {
    int64_t value=p.H;
    for(int u=0;u<N;u++) if(mask&(1<<u)) for(int v=0;v<N;v++) if(!(mask&(1<<v)) && (out[u]&(1<<v))) value+=cap[u][v];
    values[mask]=value;
  }
  int64_t bestc=std::numeric_limits<int64_t>::min(), besto=bestc;
  for(int m=1;m<N;m++) {
    int64_t cnt=0,sum=0,sq=0,anchor=-1,g=0;
    for(int mask=0;mask<SZ;mask++) if(__builtin_popcount(mask)==m){int64_t v=values[mask];if(anchor<0)anchor=v;++cnt;sum+=v;sq+=v*v;g=gcd64(g,v-anchor);}
    int64_t delta=sum-cnt*p.H, L;
    __int128 jn=__int128(sq)-__int128(sum)*p.H, jd=delta;
    if(g==0)L=anchor;else L=anchor+g*ceil_div(jn-__int128(anchor)*jd,jd*g);
    p.layer_lattice[m]=g;p.layer_floor[m]=L;
    if(m>=2&&m<=N-2&&g!=2){p.inner_bad_mask|=uint16_t(1u<<m);p.max_inner_lattice=std::max(p.max_inner_lattice,g);}
    if(m==5||m==6)bestc=std::max(bestc,L);else besto=std::max(besto,L);
  }
  p.coset_margin=bestc-besto;
  return p;
}

static void detail(uint64_t z) {
  auto out=adj_from_label(z);
  std::array<uint16_t,N> in{};
  for(int u=0;u<N;u++) for(int v=0;v<N;v++) if(out[u]&(1u<<v)) in[v]|=uint16_t(1u<<u);
  static uint32_t ending[SZ][N], starting[SZ][N], before[SZ][N], after[SZ][N];
  std::fill(&ending[0][0],&ending[0][0]+SZ*N,0); std::fill(&starting[0][0],&starting[0][0]+SZ*N,0);
  std::fill(&before[0][0],&before[0][0]+SZ*N,0); std::fill(&after[0][0],&after[0][0]+SZ*N,0);
  for(int v=0;v<N;v++) ending[1<<v][v]=starting[1<<v][v]=1;
  for(int mask=1;mask<SZ;mask++){uint16_t vs=mask;while(vs){int v=__builtin_ctz(vs);vs&=uint16_t(vs-1);int rest=mask^(1<<v);uint16_t q=uint16_t(rest)&in[v];while(q){int u=__builtin_ctz(q);q&=uint16_t(q-1);ending[mask][v]+=ending[rest][u];}q=uint16_t(rest)&out[v];while(q){int u=__builtin_ctz(q);q&=uint16_t(q-1);starting[mask][v]+=starting[rest][u];}}}
  for(int v=0;v<N;v++) before[0][v]=after[0][v]=1;
  for(int mask=1;mask<SZ;mask++)for(int v=0;v<N;v++){uint16_t q=uint16_t(mask)&in[v];while(q){int u=__builtin_ctz(q);q&=uint16_t(q-1);before[mask][v]+=ending[mask][u];}q=uint16_t(mask)&out[v];while(q){int u=__builtin_ctz(q);q&=uint16_t(q-1);after[mask][v]+=starting[mask][u];}}
  uint32_t cap[N][N]{};
  for(int a=0;a<N;a++)for(int b=a+1;b<N;b++){int u=a,v=b;if(!(out[u]&(1u<<v)))std::swap(u,v);int rem=FULL^(1<<u)^(1<<v);uint64_t c=0;for(int L=rem;;L=(L-1)&rem){int R=rem^L;c+=uint64_t(before[L][u])*after[R][v]+uint64_t(before[L][v])*after[R][u];if(!L)break;}cap[u][v]=cap[v][u]=uint32_t(c);}
  for(int i=0;i<N;i++){int64_t d=0,h=0;for(int j=0;j<N;j++)if(i!=j){d+=cap[i][j];h+=(out[i]&(1u<<j))?cap[i][j]:-int64_t(cap[i][j]);}std::cout<<"v="<<i<<" d2="<<d<<" h2="<<h<<" caps=";for(int j=0;j<N;j++)if(i!=j)std::cout<<j<<":"<<cap[i][j]<<",";std::cout<<"\n";}
  int64_t hamilton=0;for(int v=0;v<N;v++)hamilton+=ending[FULL][v];
  int64_t values[SZ]{};
  for(int mask=0;mask<SZ;mask++){
    int64_t value=hamilton;
    for(int u=0;u<N;u++)if(mask&(1<<u))for(int v=0;v<N;v++)if(!(mask&(1<<v))&&(out[u]&(1<<v)))value+=cap[u][v];
    values[mask]=value;
  }
  for(int m=1;m<N;m++){
    int64_t anchor=-1,g=0;
    for(int mask=0;mask<SZ;mask++)if(__builtin_popcount(mask)==m){if(anchor<0)anchor=values[mask];g=gcd64(g,values[mask]-anchor);}
    std::cout<<"layer="<<m<<" gcd="<<g<<" anchor="<<anchor<<"\n";
  }
}

static std::string bits(uint64_t z) { std::string s; for(int k=0;k<E;k++) s.push_back(z&(1ull<<k)?'1':'0'); return s; }
static void show(const char* tag,uint64_t z,const P&p){
  std::cout<<tag<<" rho="<<std::setprecision(15)<<p.rho<<" label="<<bits(z)<<" H="<<p.H<<" W2="<<p.W2<<" D4x4="<<p.D4x4<<" Chdx4="<<p.Chdx4<<" coset_margin="<<p.coset_margin<<" strong="<<p.strong<<"\n";
}
static uint64_t parse_bits(const std::string&s){uint64_t z=0;for(int k=0;k<(int)s.size();k++)if(s[k]=='1')z|=1ull<<k;return z;}
static void fnv_u64(uint64_t& h,uint64_t x){for(int i=0;i<8;i++){h^=(x>>(8*i))&255u;h*=1099511628211ull;}}
static uint64_t mix64(uint64_t x){x+=0x9e3779b97f4a7c15ull;x=(x^(x>>30))*0xbf58476d1ce4e5b9ull;x=(x^(x>>27))*0x94d049bb133111ebull;return x^(x>>31);}
static std::pair<uint64_t,uint64_t> row_semantic(uint64_t qz,int root,uint64_t z,const P&p){
  uint64_t first=14695981039346656037ull,second=0x243f6a8885a308d3ull;
  const auto feed=[&](uint64_t field){fnv_u64(first,field);second=mix64(second^mix64(field));};
  feed(0x4f31315041495232ull);
  feed(qz);feed(uint64_t(root));feed(z);feed(uint64_t(p.H));feed(uint64_t(p.W2));
  feed(uint64_t(p.D4x4));feed(uint64_t(p.Chdx4));feed(uint64_t(p.coset_margin));
  feed(uint64_t(p.inner_bad_mask));feed(uint64_t(p.max_inner_lattice));
  for(int m=1;m<N;m++){feed(uint64_t(p.layer_lattice[m]));feed(uint64_t(p.layer_floor[m]));}
  return {first,second};
}

static uint64_t expand_pair(uint64_t qz, int root) {
  bool q[10][10] = {};
  int k = 0;
  for (int i = 0; i < 10; ++i) for (int j = i + 1; j < 10; ++j, ++k) {
    q[i][j] = (qz >> k) & 1ull;
    q[j][i] = !q[i][j];
  }
  bool t[11][11] = {};
  for (int i = 0; i < 10; ++i) for (int j = 0; j < 10; ++j) t[i][j] = q[i][j];
  for (int u = 0; u < 10; ++u) if (u != root) {
    t[10][u] = q[root][u];
    t[u][10] = q[u][root];
  }
  t[root][10] = true;
  uint64_t z = 0; k = 0;
  for (int i = 0; i < 11; ++i) for (int j = i + 1; j < 11; ++j, ++k)
    if (t[i][j]) z |= 1ull << k;
  return z;
}

} // namespace

int main(int argc,char**argv){
  if(argc>1 && std::string(argv[1])=="--quotient-stdin") {
    std::string s; uint64_t quotients=0, rooted=0, strong_count=0, failures=0, ties=0, coset_failures=0;
    uint64_t semantic=14695981039346656037ull,semantic_sum64=0,semantic_xor64=0;
    uint64_t gz=0,gq=0,mz=0,mq=0; int groot=-1,mroot=-1; P global,minp; minp.coset_margin=std::numeric_limits<int64_t>::max();
    while(std::cin>>s) {
      if((int)s.size()!=45){std::cerr<<"bad quotient label length\n";return 2;}
      ++quotients;uint64_t qz=parse_bits(s);
      for(int root=0;root<10;root++) {
        ++rooted;uint64_t z=expand_pair(qz,root);P p=eval(z);strong_count+=p.strong;
        __int128 lhs=2*__int128(std::llabs(p.Chdx4)),rhs=p.D4x4;
        failures += lhs >= rhs;
        ties += lhs == rhs;
        coset_failures += p.coset_margin <= 0;
        fnv_u64(semantic,qz);fnv_u64(semantic,root);fnv_u64(semantic,z);
        fnv_u64(semantic,p.strong);fnv_u64(semantic,p.H);fnv_u64(semantic,p.W2);
        fnv_u64(semantic,p.D4x4);fnv_u64(semantic,p.Chdx4);fnv_u64(semantic,p.coset_margin);
        const auto row=row_semantic(qz,root,z,p);semantic_sum64+=row.first;semantic_xor64^=row.second;
        if(!gz || (__int128)std::llabs(p.Chdx4)*global.D4x4 > (__int128)std::llabs(global.Chdx4)*p.D4x4){global=p;gz=z;gq=qz;groot=root;}
        if(p.coset_margin<minp.coset_margin){minp=p;mz=z;mq=qz;mroot=root;}
      }
    }
    std::cout<<"quotients="<<quotients<<" rooted_presentations="<<rooted<<" strong_children="<<strong_count<<" rational_failures="<<failures<<" rho_equal_one="<<ties<<" coset_failures="<<coset_failures<<"\n";
    show("ROOTED_MAX_RHO",gz,global);std::cout<<"q="<<bits(gq).substr(0,45)<<" root="<<groot<<" exact_rho="<<2*std::llabs(global.Chdx4)<<"/"<<global.D4x4<<"\n";
    show("ROOTED_MIN_COSET",mz,minp);std::cout<<"q="<<bits(mq).substr(0,45)<<" root="<<mroot<<"\n";
    std::cout<<"semantic_fnv64="<<std::hex<<semantic<<std::dec<<"\n";
    std::cout<<std::hex<<std::setfill('0')<<"semantic_sum64="<<std::setw(16)<<semantic_sum64<<" semantic_xor64="<<std::setw(16)<<semantic_xor64<<std::dec<<"\n";
    return 0;
  }
  if(argc>1 && std::string(argv[1])=="--pair-any-search") {
    uint64_t seed=argc>2?std::strtoull(argv[2],nullptr,10):1;
    int restarts=argc>3?std::atoi(argv[3]):20;
    std::mt19937_64 rng(seed); P global; uint64_t gz=0,gq=0; int groot=-1;
    const uint64_t qmask=(1ull<<45)-1;
    for(int rr=0;rr<restarts;rr++) {
      uint64_t qz=rng()&qmask; int root=rng()%10;
      uint64_t z=expand_pair(qz,root); P cur=eval(z);
      for(int it=0;it<80;it++) {
        P best=cur; uint64_t bq=qz;
        for(int e=0;e<45;e++) {uint64_t yq=qz^(1ull<<e);P p=eval(expand_pair(yq,root));if(p.rho>best.rho+1e-15){best=p;bq=yq;}}
        if(bq==qz) break;
        qz=bq;cur=best;
      }
      if(cur.rho>global.rho){global=cur;gq=qz;groot=root;gz=expand_pair(gq,groot);show("PAIR_ANY_BEST",gz,global);std::cout<<"q="<<std::bitset<45>(gq)<<" root="<<groot<<"\n";}
    }
    show("PAIR_ANY_FINAL",gz,global);std::cout<<"q="<<std::bitset<45>(gq)<<" root="<<groot<<"\n";return 0;
  }
  if(argc>2 && std::string(argv[1])=="--detail") {uint64_t z=parse_bits(argv[2]);show("DETAIL",z,eval(z));detail(z);return 0;}
  if(argc>1 && std::string(argv[1])=="--pair-search") {
    uint64_t seed=argc>2?std::strtoull(argv[2],nullptr,10):1;
    int restarts=argc>3?std::atoi(argv[3]):20;
    std::mt19937_64 rng(seed); P global; uint64_t gz=0,gq=0; int groot=-1;
    const uint64_t qmask=(1ull<<45)-1;
    for(int rr=0;rr<restarts;rr++) {
      uint64_t qz=rng()&qmask; int root=rng()%10;
      uint64_t z=expand_pair(qz,root); P cur=eval(z);
      if(!cur.strong){--rr;continue;}
      for(int it=0;it<80;it++) {
        P best=cur; uint64_t bq=qz;
        for(int e=0;e<45;e++) {uint64_t yq=qz^(1ull<<e);P p=eval(expand_pair(yq,root));if(p.strong&&p.rho>best.rho+1e-15){best=p;bq=yq;}}
        if(bq==qz) break;
        qz=bq;cur=best;
      }
      if(cur.rho>global.rho){global=cur;gq=qz;groot=root;gz=expand_pair(gq,groot);show("PAIR_BEST",gz,global);std::cout<<"q="<<std::bitset<45>(gq)<<" root="<<groot<<"\n";}
    }
    show("PAIR_FINAL",gz,global);std::cout<<"q="<<std::bitset<45>(gq)<<" root="<<groot<<"\n";return 0;
  }
  if(argc>1 && std::string(argv[1])=="--stdin") {
    std::string s; uint64_t count=0, strong_count=0, failures=0, coset_failures=0, gz=0,mz=0; P global,minp; minp.coset_margin=std::numeric_limits<int64_t>::max();
    while(std::cin>>s){ if((int)s.size()!=E) { std::cerr<<"bad label length\n"; return 2; } ++count; uint64_t z=parse_bits(s); P p=eval(z); strong_count+=p.strong; if(!p.strong) continue; failures+=p.rho>1.0; coset_failures+=p.coset_margin<=0; if(p.rho>global.rho){global=p;gz=z;} if(p.coset_margin<minp.coset_margin){minp=p;mz=z;} }
    std::cout<<"count="<<count<<" strong="<<strong_count<<" rational_failures="<<failures<<" coset_failures="<<coset_failures<<"\n"; show("MAX_RHO",gz,global); show("MIN_COSET",mz,minp); return 0;
  }
  uint64_t seed=argc>1?std::strtoull(argv[1],nullptr,10):1; int restarts=argc>2?std::atoi(argv[2]):2000; std::mt19937_64 rng(seed);
  P global; uint64_t gz=0; uint64_t mask=(1ull<<E)-1;
  // Named seeds: the order-11 regular-block substitution and its converse.
  std::vector<uint64_t> starts;
  starts.push_back(parse_bits("1111111011100111111101111101111111111111111011111111111"));
  for(int rr=0;rr<restarts;rr++) {
    uint64_t z;
    if(rr<(int)starts.size()) z=starts[rr]&mask; else do z=rng()&mask; while(!strong(adj_from_label(z)));
    P cur=eval(z);
    // Greedy best-improvement, with occasional accepted two-flip escape.
    for(int it=0;it<80;it++) {
      P best=cur; uint64_t bz=z;
      for(int e=0;e<E;e++){uint64_t y=z^(1ull<<e);P q=eval(y);if(q.strong&&q.rho>best.rho+1e-15){best=q;bz=y;}}
      if(bz!=z){z=bz;cur=best;continue;}
      bool escaped=false;
      for(int k=0;k<100;k++){int e=rng()%E,f=rng()%E;if(e==f)continue;uint64_t y=z^(1ull<<e)^(1ull<<f);P q=eval(y);if(q.strong&&q.rho>cur.rho+1e-15){z=y;cur=q;escaped=true;break;}}
      if(!escaped)break;
    }
    if(cur.rho>global.rho){global=cur;gz=z;show("BEST",gz,global);}
    if(global.rho>1.0) break;
  }
  show("FINAL",gz,global);
}
