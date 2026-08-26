#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace {
constexpr int QN=10, N=11, QSZ=1<<QN, SZ=1<<N;
constexpr int QFULL=QSZ-1, FULL=SZ-1, E=N*(N-1)/2;

static uint32_t qe[QSZ][QN], qs[QSZ][QN];
static uint32_t ce[SZ][N], cs[SZ][N];
static uint32_t before_[SZ][N], after_[SZ][N];

static constexpr std::size_t ix(int value){return static_cast<std::size_t>(value);}

struct Profile {
  int64_t H=0,W2=0,D4x4=0,Chdx4=0;
  int64_t margin=std::numeric_limits<int64_t>::max();
  bool inner_lattice_non2=false;
  uint16_t inner_bad_mask=0;
  int64_t max_inner_lattice=0;
  std::array<int64_t,N+1> layer_lattice{};
  std::array<int64_t,N+1> layer_floor{};
};

static int64_t gcd64(int64_t a,int64_t b){return std::gcd(std::llabs(a),std::llabs(b));}
static int64_t ceil_div(__int128 a,__int128 b){
  __int128 q=a/b,r=a%b;
  if(r>0)++q;
  return int64_t(q);
}
static int64_t choose(int n,int k){
  if(k<0||k>n)return 0;
  int64_t z=1;
  for(int i=1;i<=k;i++)z=z*(n-k+i)/i;
  return z;
}
static uint64_t parse_bits(const std::string&s){
  uint64_t z=0;for(int k=0;k<(int)s.size();k++)if(s[ix(k)]=='1')z|=1ull<<k;return z;
}
static std::string bits(uint64_t z){
  std::string s;for(int k=0;k<E;k++)s.push_back((z>>k)&1?'1':'0');return s;
}
static std::array<uint16_t,QN> qadj(uint64_t z){
  std::array<uint16_t,QN> out{};int k=0;
  for(int i=0;i<QN;i++)for(int j=i+1;j<QN;++j,++k){
    if((z>>k)&1)out[ix(i)]|=uint16_t(1u<<j);else out[ix(j)]|=uint16_t(1u<<i);
  }
  return out;
}
static std::array<uint16_t,N> child_adj(const std::array<uint16_t,QN>&q,int root){
  std::array<uint16_t,N> out{};
  for(int u=0;u<QN;u++)out[ix(u)]=q[ix(u)];
  for(int u=0;u<QN;u++)if(u!=root){
    if(q[ix(root)]&(1u<<u))out[ix(QN)]|=uint16_t(1u<<u);else out[ix(u)]|=uint16_t(1u<<QN);
  }
  out[ix(root)]|=uint16_t(1u<<QN);
  return out;
}
static uint64_t child_label(const std::array<uint16_t,N>&out){
  uint64_t z=0;int k=0;
  for(int i=0;i<N;i++)for(int j=i+1;j<N;j++,k++)if(out[ix(i)]&(1u<<j))z|=1ull<<k;
  return z;
}

static void base_endpoint_dp(const std::array<uint16_t,QN>&out){
  std::fill(&qe[0][0],&qe[0][0]+QSZ*QN,0);
  std::fill(&qs[0][0],&qs[0][0]+QSZ*QN,0);
  std::array<uint16_t,QN> in{};
  for(int u=0;u<QN;u++)for(int v=0;v<QN;v++)if(out[ix(u)]&(1u<<v))in[ix(v)]|=uint16_t(1u<<u);
  for(int v=0;v<QN;v++)qe[1<<v][v]=qs[1<<v][v]=1;
  for(int mask=1;mask<QSZ;mask++){
    uint16_t vs=uint16_t(mask);
    while(vs){
      int v=__builtin_ctz(vs);vs&=uint16_t(vs-1);int rest=mask^(1<<v);
      uint16_t p=uint16_t(rest)&in[ix(v)];uint32_t x=0;
      while(p){int u=__builtin_ctz(p);p&=uint16_t(p-1);x+=qe[rest][u];}qe[mask][v]=x+(rest?0:1);
      p=uint16_t(rest)&out[ix(v)];x=0;
      while(p){int u=__builtin_ctz(p);p&=uint16_t(p-1);x+=qs[rest][u];}qs[mask][v]=x+(rest?0:1);
    }
  }
}

// Materialize the virtual child endpoint tables.  The no-clone and one-clone
// slices are inherited verbatim from Q.  Only entries with both clones and an
// endpoint in U are new stored state: 2*(QN-1)*2^(QN-2)=4608 integers.
static int64_t rooted_endpoints(const std::array<uint16_t,QN>&q,int root){
  const int rb=1<<root, bb=1<<QN;
  const int umask=QFULL^rb;
  for(int A=0;A<QSZ;A++){
    if(A&rb)continue;
    const int qm=A|rb, ma=A|rb, mb=A|bb;
    uint16_t p=uint16_t(A);
    while(p){
      int u=__builtin_ctz(p);p&=uint16_t(p-1);
      ce[A][u]=qe[A][u];cs[A][u]=qs[A][u];
      ce[ma][u]=ce[mb][u]=qe[qm][u];
      cs[ma][u]=cs[mb][u]=qs[qm][u];
    }
    ce[ma][root]=ce[mb][QN]=qe[qm][root];
    cs[ma][root]=cs[mb][QN]=qs[qm][root];

    uint32_t ea=0,sb=0;
    p=uint16_t(A);
    while(p){
      int u=__builtin_ctz(p);p&=uint16_t(p-1);
      if(q[ix(u)]&rb)ea+=qe[qm][u];
      if(q[ix(root)]&(1u<<u))sb+=qs[qm][u];
    }
    const uint32_t eb=ea+qe[qm][root];
    const uint32_t sa=sb+qs[qm][root];
    const int both=A|rb|bb;
    ce[both][root]=ea;ce[both][QN]=eb;
    cs[both][root]=sa;cs[both][QN]=sb;

    p=uint16_t(A);
    while(p){
      int u=__builtin_ctz(p);p&=uint16_t(p-1);const int rest=A^(1<<u), restboth=rest|rb|bb;
      uint32_t ev=0,sv=0;uint16_t z=uint16_t(rest);
      while(z){
        int v=__builtin_ctz(z);z&=uint16_t(z-1);
        if(q[ix(v)]&(1u<<u))ev+=ce[restboth][v];
        if(q[ix(u)]&(1u<<v))sv+=cs[restboth][v];
      }
      if(q[ix(root)]&(1u<<u))ev+=ce[restboth][root]+ce[restboth][QN];
      else sv+=cs[restboth][root]+cs[restboth][QN];
      ce[both][u]=ev;cs[both][u]=sv;
    }
  }
  int64_t H=0;for(int v=0;v<N;v++)H+=ce[FULL][v];
  (void)umask;
  return H;
}

static Profile evaluate_root(const std::array<uint16_t,QN>&q,int root){
  const auto out=child_adj(q,root);
  Profile p;p.H=rooted_endpoints(q,root);

  // Only mask/target pairs with target outside the mask are ever consumed.
  for(int mask=0;mask<SZ;mask++)for(int x=0;x<N;x++)if(!(mask&(1<<x))){
    if(!mask){before_[mask][x]=after_[mask][x]=1;continue;}
    uint32_t bv=0,av=0;uint16_t z=uint16_t(mask);
    while(z){
      int u=__builtin_ctz(z);z&=uint16_t(z-1);
      if(out[ix(u)]&(1u<<x))bv+=ce[mask][u];
      if(out[ix(x)]&(1u<<u))av+=cs[mask][u];
    }
    before_[mask][x]=bv;after_[mask][x]=av;
  }

  int64_t cap[N][N]{};
  for(int x=0;x<N;x++)for(int y=x+1;y<N;y++){
    int u=x,v=y;if(!(out[ix(u)]&(1u<<v)))std::swap(u,v);
    const int rem=FULL^(1<<u)^(1<<v);int left=rem;uint64_t c=0;
    while(true){
      const int right=rem^left;
      c+=uint64_t(before_[left][u])*after_[right][v]
        +uint64_t(before_[left][v])*after_[right][u];
      if(!left)break;left=(left-1)&rem;
    }
    cap[x][y]=cap[y][x]=int64_t(c);p.W2+=int64_t(c);
  }
  int64_t parentH=0;for(int v=0;v<QN;v++)parentH+=qe[QFULL][v];
  if(cap[root][QN]!=2*parentH){std::cerr<<"internal exposure failure\n";std::exit(3);}

  int64_t degree[N]{},signed_degree[N]{},outcap[N]{},incap[N]{};
  int64_t qsum=0,outstar=0,instar=0;
  for(int i=0;i<N;i++)for(int j=i+1;j<N;j++)qsum+=cap[i][j]*cap[i][j];
  for(int v=0;v<N;v++)for(int u=0;u<N;u++)if(u!=v){
    degree[v]+=cap[v][u];
    if(out[ix(v)]&(1u<<u)){signed_degree[v]+=cap[v][u];outcap[v]+=cap[v][u];}
    else{signed_degree[v]-=cap[v][u];incap[v]+=cap[v][u];}
  }
  int64_t degree_sq=0;
  for(int v=0;v<N;v++){
    p.Chdx4+=degree[v]*signed_degree[v];degree_sq+=degree[v]*degree[v];
    int64_t oq=0,iq=0;
    for(int u=0;u<N;u++)if(u!=v){if(out[ix(v)]&(1u<<u))oq+=cap[v][u]*cap[v][u];else iq+=cap[v][u]*cap[v][u];}
    outstar+=(outcap[v]*outcap[v]-oq)/2;
    instar+=(incap[v]*incap[v]-iq)/2;
  }
  p.D4x4=(p.W2*p.W2+qsum-degree_sq)/2;
  if(p.D4x4<=0){std::cerr<<"D4 failure\n";std::exit(4);}

  int64_t anchors[N+1]{};anchors[0]=p.H;
  for(int m=1;m<=N;m++){
    int mask=(1<<m)-1;int64_t f=p.H;
    for(int u=0;u<N;u++)if(mask&(1<<u))for(int v=0;v<N;v++)if(!(mask&(1<<v))&&(out[ix(u)]&(1u<<v)))f+=cap[u][v];
    anchors[m]=f;
  }

  int64_t latt[N+1]{};
  for(int i=0;i<N;i++)for(int j=i+1;j<N;j++){
    int rest[N-2],bn[N-2],len=0;
    for(int u=0;u<N;u++)if(u!=i&&u!=j){rest[len]=u;bn[len]=int(cap[j][u]-cap[i][u]);len++;}
    int64_t gd=0;for(int z=1;z<len;z++)gd=gcd64(gd,int64_t(bn[z])-bn[0]);
    int64_t prefix=0,base=outcap[j]-outcap[i];
    for(int k=0;k<=len;k++){
      int64_t pg=std::llabs(base-prefix);
      if(k>0&&k<len)pg=gcd64(pg,gd);
      latt[k+1]=gcd64(latt[k+1],pg);
      if(k<len)prefix+=bn[k];
    }
  }
  for(int m=2;m<=N-2;m++)if(latt[m]!=2){
    p.inner_lattice_non2=true;p.inner_bad_mask|=uint16_t(1u<<m);
    p.max_inner_lattice=std::max(p.max_inner_lattice,latt[m]);
  }

  int64_t bestc=std::numeric_limits<int64_t>::min(),besto=bestc;
  for(int m=1;m<N;m++){
    const int64_t count=choose(N,m),p1=choose(N-2,m-1);
    const int64_t sum=count*p.H+p1*p.W2;
    const int64_t sq=count*p.H*p.H+p1*(2*p.H*p.W2+qsum)
      +2*(choose(N-3,m-1)*outstar+choose(N-3,m-2)*instar+choose(N-4,m-2)*p.D4x4);
    const int64_t jd=sum-count*p.H;
    const __int128 jn=__int128(sq)-__int128(sum)*p.H;
    int64_t floor=anchors[m];
    if(latt[m])floor+=latt[m]*ceil_div(jn-__int128(anchors[m])*jd,__int128(jd)*latt[m]);
    else if(jn!=__int128(anchors[m])*jd){std::cerr<<"constant layer failure\n";std::exit(5);}
    p.layer_lattice[ix(m)]=latt[m];p.layer_floor[ix(m)]=floor;
    if(m==5||m==6)bestc=std::max(bestc,floor);else besto=std::max(besto,floor);
  }
  p.margin=bestc-besto;
  return p;
}

static void show(const char*tag,uint64_t z,const Profile&p){
  std::cout<<tag<<" label="<<bits(z)<<" H="<<p.H<<" W2="<<p.W2
           <<" D4x4="<<p.D4x4<<" Chdx4="<<p.Chdx4
           <<" exact_rho="<<2*std::llabs(p.Chdx4)<<"/"<<p.D4x4
           <<" coset_margin="<<p.margin<<"\n";
}

static uint64_t mix64(uint64_t x){
  x+=0x9e3779b97f4a7c15ull;
  x=(x^(x>>30))*0xbf58476d1ce4e5b9ull;
  x=(x^(x>>27))*0x94d049bb133111ebull;
  return x^(x>>31);
}
static void fnv_u64(uint64_t&h,uint64_t x){
  for(int byte=0;byte<8;byte++){h^=(x>>(8*byte))&0xffull;h*=1099511628211ull;}
}
static std::pair<uint64_t,uint64_t> row_semantic(
    uint64_t qz,int root,uint64_t childz,const Profile&p){
  uint64_t first=14695981039346656037ull;
  uint64_t second=0x243f6a8885a308d3ull;
  const auto feed=[&](uint64_t field){
    fnv_u64(first,field);second=mix64(second^mix64(field));
  };
  feed(0x4f31315041495232ull); // "O11PAIR2": semantic schema/domain word.
  feed(qz);feed(static_cast<uint64_t>(root));feed(childz);
  feed(static_cast<uint64_t>(p.H));feed(static_cast<uint64_t>(p.W2));
  feed(static_cast<uint64_t>(p.D4x4));feed(static_cast<uint64_t>(p.Chdx4));
  feed(static_cast<uint64_t>(p.margin));feed(static_cast<uint64_t>(p.inner_bad_mask));
  feed(static_cast<uint64_t>(p.max_inner_lattice));
  for(int m=1;m<N;m++){
    feed(static_cast<uint64_t>(p.layer_lattice[ix(m)]));
    feed(static_cast<uint64_t>(p.layer_floor[ix(m)]));
  }
  return {first,second};
}
}

int main(){
  std::string s;uint64_t quotients=0,rooted=0,rational_failures=0,coset_failures=0,inner_lattice_non2=0;
  uint64_t semantic_sum64=0,semantic_xor64=0;
  uint64_t maxz=0,minz=0,maxq=0,minq=0;int maxr=-1,minr=-1;Profile maxp,minp;minp.margin=std::numeric_limits<int64_t>::max();
  uint64_t badz=0,badq=0;int badr=-1;Profile badp;
  while(std::cin>>s){
    if(s.size()!=45){std::cerr<<"bad quotient label length\n";return 2;}
    const uint64_t qz=parse_bits(s);const auto q=qadj(qz);base_endpoint_dp(q);quotients++;
    for(int r=0;r<QN;r++){
      Profile p=evaluate_root(q,r);const uint64_t childz=child_label(child_adj(q,r));rooted++;
      rational_failures+=2*std::llabs(p.Chdx4)>=p.D4x4;
      coset_failures+=p.margin<=0;
      inner_lattice_non2+=p.inner_lattice_non2;
      const auto semantic=row_semantic(qz,r,childz,p);
      semantic_sum64+=semantic.first;semantic_xor64^=semantic.second;
      if(p.inner_lattice_non2&&!badz){badp=p;badz=childz;badq=qz;badr=r;}
      if(!maxz||__int128(std::llabs(p.Chdx4))*maxp.D4x4>__int128(std::llabs(maxp.Chdx4))*p.D4x4){maxp=p;maxz=childz;maxq=qz;maxr=r;}
      if(p.margin<minp.margin){minp=p;minz=childz;minq=qz;minr=r;}
    }
  }
  std::cout<<"quotients="<<quotients<<" rooted_presentations="<<rooted
           <<" rational_failures="<<rational_failures<<" coset_failures="<<coset_failures
           <<" inner_lattice_non2="<<inner_lattice_non2<<"\n";
  show("ROOTED_MAX_RHO",maxz,maxp);std::cout<<"q="<<bits(maxq).substr(0,45)<<" root="<<maxr<<"\n";
  show("ROOTED_MIN_COSET",minz,minp);std::cout<<"q="<<bits(minq).substr(0,45)<<" root="<<minr<<"\n";
  if(badz){show("FIRST_INNER_LATTICE_NON2",badz,badp);std::cout<<"q="<<bits(badq).substr(0,45)<<" root="<<badr<<" bad_mask="<<badp.inner_bad_mask<<" max_gcd="<<badp.max_inner_lattice<<"\n";}
  std::cout<<std::hex<<std::setfill('0')
           <<"semantic_sum64="<<std::setw(16)<<semantic_sum64
           <<" semantic_xor64="<<std::setw(16)<<semantic_xor64<<std::dec<<"\n";
}
