#include <algorithm>
#include <array>
#include <atomic>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

namespace {

constexpr int N = 11;
constexpr int SZ = 1 << N;
constexpr int FULL = SZ - 1;
constexpr int SHARDS = 16;

struct Sha256 {
  std::array<uint32_t,8> h{{0x6a09e667u,0xbb67ae85u,0x3c6ef372u,0xa54ff53au,
                             0x510e527fu,0x9b05688cu,0x1f83d9abu,0x5be0cd19u}};
  std::array<unsigned char,64> buf{};
  uint64_t bytes=0;
  size_t used=0;
  static uint32_t rotr(uint32_t x,int n){return (x>>n)|(x<<(32-n));}
  void block(const unsigned char* p){
    static constexpr uint32_t k[64]={
      0x428a2f98u,0x71374491u,0xb5c0fbcfu,0xe9b5dba5u,0x3956c25bu,0x59f111f1u,0x923f82a4u,0xab1c5ed5u,
      0xd807aa98u,0x12835b01u,0x243185beu,0x550c7dc3u,0x72be5d74u,0x80deb1feu,0x9bdc06a7u,0xc19bf174u,
      0xe49b69c1u,0xefbe4786u,0x0fc19dc6u,0x240ca1ccu,0x2de92c6fu,0x4a7484aau,0x5cb0a9dcu,0x76f988dau,
      0x983e5152u,0xa831c66du,0xb00327c8u,0xbf597fc7u,0xc6e00bf3u,0xd5a79147u,0x06ca6351u,0x14292967u,
      0x27b70a85u,0x2e1b2138u,0x4d2c6dfcu,0x53380d13u,0x650a7354u,0x766a0abbu,0x81c2c92eu,0x92722c85u,
      0xa2bfe8a1u,0xa81a664bu,0xc24b8b70u,0xc76c51a3u,0xd192e819u,0xd6990624u,0xf40e3585u,0x106aa070u,
      0x19a4c116u,0x1e376c08u,0x2748774cu,0x34b0bcb5u,0x391c0cb3u,0x4ed8aa4au,0x5b9cca4fu,0x682e6ff3u,
      0x748f82eeu,0x78a5636fu,0x84c87814u,0x8cc70208u,0x90befffau,0xa4506cebu,0xbef9a3f7u,0xc67178f2u};
    uint32_t w[64];
    for(int i=0;i<16;i++) w[i]=(uint32_t(p[4*i])<<24)|(uint32_t(p[4*i+1])<<16)|(uint32_t(p[4*i+2])<<8)|p[4*i+3];
    for(int i=16;i<64;i++){uint32_t s0=rotr(w[i-15],7)^rotr(w[i-15],18)^(w[i-15]>>3);uint32_t s1=rotr(w[i-2],17)^rotr(w[i-2],19)^(w[i-2]>>10);w[i]=w[i-16]+s0+w[i-7]+s1;}
    uint32_t a=h[0],b=h[1],c=h[2],d=h[3],e=h[4],f=h[5],g=h[6],hh=h[7];
    for(int i=0;i<64;i++){uint32_t S1=rotr(e,6)^rotr(e,11)^rotr(e,25);uint32_t ch=(e&f)^((~e)&g);uint32_t t1=hh+S1+ch+k[i]+w[i];uint32_t S0=rotr(a,2)^rotr(a,13)^rotr(a,22);uint32_t maj=(a&b)^(a&c)^(b&c);uint32_t t2=S0+maj;hh=g;g=f;f=e;e=d+t1;d=c;c=b;b=a;a=t1+t2;}
    h[0]+=a;h[1]+=b;h[2]+=c;h[3]+=d;h[4]+=e;h[5]+=f;h[6]+=g;h[7]+=hh;
  }
  void update(const void* vp,size_t n){const auto*p=static_cast<const unsigned char*>(vp);bytes+=n;while(n){size_t take=std::min(n,64-used);std::memcpy(buf.data()+used,p,take);used+=take;p+=take;n-=take;if(used==64){block(buf.data());used=0;}}}
  std::string finish(){uint64_t bits=bytes*8;unsigned char one=0x80;update(&one,1);unsigned char zero=0;while(used!=56)update(&zero,1);unsigned char len[8];for(int i=0;i<8;i++)len[7-i]=static_cast<unsigned char>(bits>>(8*i));update(len,8);std::ostringstream o;for(uint32_t x:h)o<<std::hex<<std::setw(8)<<std::setfill('0')<<x;return o.str();}
};

template<class U> void sha_le(Sha256& s,U x){static_assert(std::is_unsigned<U>::value,"unsigned");unsigned char b[sizeof(U)];for(size_t i=0;i<sizeof(U);i++){b[i]=x&255u;x>>=8;}s.update(b,sizeof(b));}
void sha_i64(Sha256&s,int64_t x){sha_le<uint64_t>(s,static_cast<uint64_t>(x));}
std::string sha_text(const std::string& x){Sha256 s;s.update(x.data(),x.size());return s.finish();}
std::string sha_file(const std::string& path){std::ifstream f(path,std::ios::binary);if(!f)throw std::runtime_error("open hash file: "+path);Sha256 s;std::array<char,1<<16>b{};while(f){f.read(b.data(),b.size());std::streamsize n=f.gcount();if(n>0)s.update(b.data(),size_t(n));}return s.finish();}

[[noreturn]] void fail(const std::string& x){throw std::runtime_error(x);}
void require(bool x,const std::string& why){if(!x)fail(why);}
int64_t gcd64(int64_t a,int64_t b){return std::gcd(std::llabs(a),std::llabs(b));}
int64_t ceil_div(__int128 a,__int128 b){require(b>0,"ceil denominator");__int128 q=a/b,r=a%b;if(r>0)++q;require(q>=std::numeric_limits<int64_t>::min()&&q<=std::numeric_limits<int64_t>::max(),"ceil overflow");return static_cast<int64_t>(q);}

struct Rat { int64_t n=0,d=1; };
Rat rat(int64_t n,int64_t d){require(d!=0,"rat denominator");if(d<0){n=-n;d=-d;}int64_t g=gcd64(n,d);return {n/g,d/g};}
int cmp(Rat a,Rat b){__int128 x=__int128(a.n)*b.d,y=__int128(b.n)*a.d;return (x>y)-(x<y);}
std::string rs(Rat x){return std::to_string(x.n)+"/"+std::to_string(x.d);}

struct Stream { int n=0; bool strong=false; std::vector<std::string> labels; std::vector<std::vector<uint16_t>> decoded; std::string digest; };

Stream read_stream(const std::string& exe,int n,bool strong){
  std::string cmd=exe+" -q "+(strong?"-c ":"")+std::to_string(n);
  FILE* f=popen(cmd.c_str(),"r");if(!f)fail("popen: "+cmd);
  Stream s; s.n=n;s.strong=strong;Sha256 hash;char buf[256];
  while(std::fgets(buf,sizeof(buf),f)){std::string line(buf);while(!line.empty()&&(line.back()=='\n'||line.back()=='\r'))line.pop_back();if(line.empty())continue;require(int(line.size())==n*(n-1)/2,"gentourng label length");s.labels.push_back(line);hash.update(line.data(),line.size());char nl='\n';hash.update(&nl,1);}
  int rc=pclose(f);require(rc==0,"gentourng exit: "+cmd);s.digest=hash.finish();return s;
}

std::vector<uint16_t> decode(const std::string& s,int n){
  std::vector<uint16_t> out(n);int z=0;
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,z++){char c=s[z];require(c=='0'||c=='1',"binary label");if(c=='1')out[i]|=uint16_t(1u<<j);else out[j]|=uint16_t(1u<<i);}return out;
}

uint64_t encode11(const std::array<uint16_t,N>& out){uint64_t z=0;int e=0;for(int i=0;i<N;i++)for(int j=i+1;j<N;j++,e++)if(out[i]&(1u<<j))z|=uint64_t(1)<<e;return z;}
std::string bits11(uint64_t z){std::string s;for(int e=0;e<55;e++)s.push_back(z&(uint64_t(1)<<e)?'1':'0');return s;}

std::array<uint16_t,N> substitute(const std::vector<uint16_t>& Q,int k,const std::vector<uint16_t>& B){
  int q=Q.size(),r=B.size();require(q+r-1==N,"substitution order");std::vector<int> size(q,1),start(q);size[k]=r;int cur=0;for(int i=0;i<q;i++){start[i]=cur;cur+=size[i];}
  std::array<uint16_t,N> out{};
  for(int i=0;i<r;i++)for(int j=0;j<r;j++)if(B[i]&(1u<<j))out[start[k]+i]|=uint16_t(1u<<(start[k]+j));
  for(int i=0;i<q;i++)for(int j=0;j<q;j++)if(Q[i]&(1u<<j))for(int u=start[i];u<start[i]+size[i];u++)for(int v=start[j];v<start[j]+size[j];v++)out[u]|=uint16_t(1u<<v);
  return out;
}

bool strong11(const std::array<uint16_t,N>& out){
  std::array<uint16_t,N> in{};for(int u=0;u<N;u++)for(int v=0;v<N;v++)if(out[u]&(1u<<v))in[v]|=uint16_t(1u<<u);
  for(int rev=0;rev<2;rev++){uint16_t seen=1,todo=1;while(todo){int v=__builtin_ctz(todo);todo&=uint16_t(todo-1);uint16_t fresh=(rev?in[v]:out[v])&uint16_t(~seen);seen|=fresh;todo|=fresh;}if((seen&FULL)!=FULL)return false;}return true;
}

struct Layer { Rat J; int64_t L=0,A=0,g=0; };
struct Profile {
  uint64_t label=0,response_hash=0;
  int64_t H=0,W2=0,D4x4=0,C4=0,margin=0;
  bool strong=false,rational_failure=false,coset_failure=false;
  std::array<Layer,N-1> layer{};
};

uint64_t fnv_values(const std::array<int64_t,SZ>& values){uint64_t h=1469598103934665603ull;for(int64_t v:values)for(int i=0;i<8;i++){h^=(uint64_t(v)>>(8*i))&255u;h*=1099511628211ull;}return h;}

Profile evaluate(const std::array<uint16_t,N>& out){
  Profile p;p.label=encode11(out);p.strong=strong11(out);if(!p.strong)return p;
  std::array<uint16_t,N> in{};for(int u=0;u<N;u++)for(int v=0;v<N;v++)if(out[u]&(1u<<v))in[v]|=uint16_t(1u<<u);
  uint32_t ending[SZ][N]{};uint32_t starting[SZ][N]{};
  for(int v=0;v<N;v++)ending[1<<v][v]=starting[1<<v][v]=1;
  for(int mask=1;mask<SZ;mask++){
    uint16_t vs=mask;while(vs){int v=__builtin_ctz(vs);vs&=uint16_t(vs-1);int rest=mask^(1<<v);if(!rest)continue;
      uint16_t q=uint16_t(rest)&in[v];while(q){int u=__builtin_ctz(q);q&=uint16_t(q-1);ending[mask][v]+=ending[rest][u];}
      q=uint16_t(rest)&out[v];while(q){int u=__builtin_ctz(q);q&=uint16_t(q-1);starting[mask][v]+=starting[rest][u];}
    }
  }
  std::array<uint32_t,N> starts{},ends{};for(int v=0;v<N;v++){starts[v]=starting[FULL][v];ends[v]=ending[FULL][v];p.H+=ends[v];}
  uint32_t exposed[N][N]{};
  for(int a=0;a<N;a++)for(int b=0;b<N;b++)if(a!=b){int rem=FULL^(1<<a)^(1<<b);uint64_t q=0;for(int mid=rem;;mid=(mid-1)&rem){int left=mid|(1<<a);int right=(rem^mid)|(1<<b);q+=uint64_t(ending[left][a])*starting[right][b];if(!mid)break;}require(q<=std::numeric_limits<uint32_t>::max(),"exposed overflow");exposed[a][b]=q;}
  std::array<int64_t,SZ> values{};std::array<int64_t,N> linear{};
  int64_t base=std::accumulate(ends.begin(),ends.end(),int64_t(0));require(base==p.H,"empty response");
  for(int v=0;v<N;v++){linear[v]=int64_t(starts[v])-ends[v];for(int a=0;a<N;a++)if(a!=v)linear[v]+=exposed[a][v];}
  values[0]=base;
  for(int mask=1;mask<SZ;mask++){int bit=mask&-mask,v=__builtin_ctz(bit),old=mask^bit;int64_t penalty=0;for(int m=old;m;m&=m-1){int u=__builtin_ctz(m);penalty+=exposed[u][v]+exposed[v][u];}values[mask]=values[old]+linear[v]-penalty;}
  require(values[FULL]==p.H,"full response");for(auto v:values)require(v>=p.H&&(v&1),"response support/parity");p.response_hash=fnv_values(values);

  uint32_t cap[N][N]{};int64_t d2[N]{},h2[N]{};
  for(int i=0;i<N;i++)for(int j=i+1;j<N;j++){
    int64_t c=values[1<<i]+values[1<<j]-p.H-values[(1<<i)|(1<<j)];require(c>=0&&(c%2==0),"curvature");require(c<=std::numeric_limits<uint32_t>::max(),"capacity overflow");cap[i][j]=cap[j][i]=uint32_t(c);p.W2+=c;
  }
  for(int v=0;v<N;v++)for(int u=0;u<N;u++)if(u!=v){d2[v]+=cap[v][u];h2[v]+=(out[v]&(1u<<u))?cap[v][u]:-int64_t(cap[v][u]);}
  for(int v=0;v<N;v++){int64_t from_single=2*(values[1<<v]-p.H)-d2[v];require(from_single==h2[v],"field/divergence");p.C4+=h2[v]*d2[v];}
  for(int a=0;a<N;a++)for(int b=a+1;b<N;b++)for(int c=a+1;c<N;c++)if(c!=b)for(int d=c+1;d<N;d++)if(d!=b)p.D4x4+=int64_t(cap[a][b])*cap[c][d];
  require(p.D4x4>0,"D4");p.rational_failure=2*std::llabs(p.C4)>=p.D4x4;

  Rat bestc{std::numeric_limits<int64_t>::min(),1},besto=bestc;int64_t bestcL=std::numeric_limits<int64_t>::min(),bestoL=bestcL;
  for(int m=1;m<N;m++){
    int64_t cnt=0,sum=0,sq=0,anchor=-1,g=0,A=0;
    for(int mask=0;mask<SZ;mask++)if(__builtin_popcount(mask)==m){int64_t v=values[mask];if(anchor<0)anchor=v;++cnt;sum+=v;sq+=v*v;A=std::max(A,v);g=gcd64(g,v-anchor);}
    int64_t delta=sum-cnt*p.H;require(delta>0,"layer mean");Rat J=rat(sq-sum*p.H,delta);int64_t L;if(g==0){require(cmp(J,{anchor,1})==0,"constant layer");L=anchor;}else L=anchor+g*ceil_div(__int128(J.n)-__int128(anchor)*J.d,__int128(J.d)*g);require(L<=A,"coset support");p.layer[m-1]={J,L,A,g};
    bool central=(m==5||m==6);if(central){if(cmp(J,bestc)>0)bestc=J;bestcL=std::max(bestcL,L);}else{if(cmp(J,besto)>0)besto=J;bestoL=std::max(bestoL,L);}
  }
  require((cmp(besto,bestc)>=0)==p.rational_failure,"criterion/layer comparison");p.margin=bestcL-bestoL;p.coset_failure=p.margin<=0;return p;
}

struct Extreme { bool set=false;Profile p;int q=0,r=0;size_t qi=0,bi=0;int k=0;uint64_t multiplicity=0; };
bool better_rho(const Profile&a,const Profile&b){return __int128(std::llabs(a.C4))*b.D4x4>__int128(std::llabs(b.C4))*a.D4x4;}
bool same_rho(const Profile&a,const Profile&b){return __int128(std::llabs(a.C4))*b.D4x4==__int128(std::llabs(b.C4))*a.D4x4;}
bool row_less(const Extreme&a,const Extreme&b){return std::tie(a.p.label,a.q,a.r,a.qi,a.bi,a.k)<std::tie(b.p.label,b.q,b.r,b.qi,b.bi,b.k);}
void add_max(Extreme& e,const Extreme& x){if(!e.set||better_rho(x.p,e.p)){e=x;e.set=true;}else if(same_rho(x.p,e.p)){uint64_t m=e.multiplicity+x.multiplicity;if(row_less(x,e))e=x;e.multiplicity=m;}}
void add_min(Extreme& e,const Extreme& x){if(!e.set||x.p.margin<e.p.margin){e=x;e.set=true;}else if(x.p.margin==e.p.margin){uint64_t m=e.multiplicity+x.multiplicity;if(row_less(x,e))e=x;e.multiplicity=m;}}

struct Task {int q=0,r=0,shard=0;};
struct Result {Task t;uint64_t rows=0,strong=0,rfail=0,cfail=0;Extreme mx,mn;std::string profile_digest;};

void hash_profile(Sha256& h,int q,int r,size_t qi,size_t bi,int k,const Profile&p){
  sha_le<uint32_t>(h,q);sha_le<uint32_t>(h,r);sha_le<uint64_t>(h,qi);sha_le<uint64_t>(h,bi);sha_le<uint32_t>(h,k);sha_le<uint64_t>(h,p.label);sha_i64(h,p.H);sha_i64(h,p.W2);sha_i64(h,p.D4x4);sha_i64(h,p.C4);sha_le<uint64_t>(h,p.response_hash);sha_i64(h,p.margin);
  for(const auto& L:p.layer){sha_i64(h,L.J.n);sha_i64(h,L.J.d);sha_i64(h,L.L);sha_i64(h,L.A);sha_i64(h,L.g);}
}

Result run_task(const Task&t,const Stream& qs,const Stream& bs){
  Result z;z.t=t;Sha256 profile;
  uint64_t total=uint64_t(qs.labels.size())*bs.labels.size()*t.q;
  for(uint64_t idx=t.shard;idx<total;idx+=SHARDS){uint64_t x=idx;int k=x%t.q;x/=t.q;size_t bi=x%bs.labels.size();size_t qi=x/bs.labels.size();auto out=substitute(qs.decoded[qi],k,bs.decoded[bi]);Profile p=evaluate(out);++z.rows;z.strong+=p.strong;z.rfail+=p.rational_failure;z.cfail+=p.coset_failure;hash_profile(profile,t.q,t.r,qi,bi,k,p);Extreme e{true,p,t.q,t.r,qi,bi,k,1};add_max(z.mx,e);add_min(z.mn,e);}
  z.profile_digest=profile.finish();return z;
}

std::string extreme_text(const Extreme&e){Rat rho=rat(2*std::llabs(e.p.C4),e.p.D4x4);std::ostringstream o;o<<"q="<<e.q<<",r="<<e.r<<",qi="<<e.qi<<",bi="<<e.bi<<",k="<<e.k<<",label="<<bits11(e.p.label)<<",H="<<e.p.H<<",W="<<rs(rat(e.p.W2,2))<<",D4="<<rs(rat(e.p.D4x4,4))<<",Chd="<<rs(rat(e.p.C4,4))<<",rho="<<rs(rho)<<",margin="<<e.p.margin<<",response_fnv="<<std::hex<<std::setw(16)<<std::setfill('0')<<e.p.response_hash<<std::dec<<",multiplicity="<<e.multiplicity;return o.str();}

} // namespace

int main(int argc,char**argv){try{
  std::string exe="/opt/homebrew/bin/gentourng";int threads=8,q_only=0;
  for(int i=1;i<argc;i++){std::string a=argv[i];if(a=="--gentourng"&&i+1<argc)exe=argv[++i];else if(a=="--threads"&&i+1<argc)threads=std::stoi(argv[++i]);else if(a=="--q-only"&&i+1<argc)q_only=std::stoi(argv[++i]);else fail("unknown argument: "+a);}require(q_only==0||(q_only>=3&&q_only<=9),"q-only range");
  const std::array<uint64_t,10> all_expected{{0,0,0,2,4,12,56,456,6880,191536}};
  const std::array<uint64_t,10> strong_expected{{0,0,0,1,1,6,35,353,6008,178133}};
  const std::array<std::string,10> all_hash{{"","","",
    "b8bd9f18dbaa6e48bf08aea7368585a4e3e001c4914a4012b138c7aa8b1bb6a0",
    "91dde84a6b6286a5e0b7b6295ccb1d16955dc4d74b8b93163d50d8a9cd7c7921",
    "5d16cab0f65c58f402ef49c3442e3e8fedbb5498a8cb09f58a5d8d733784f6e4",
    "814e0ed10e5d3aaae92d809b9eb915f1eb5542a65a33fe4d47463d22159a525b",
    "164260b94960af0cc63faf3f178ceb95f4dd23bbca376ec23872c33b30d94261",
    "fc96c6997724e54ccea3bd166f4117d9e27925d85f568e31b0623527e5139dad",
    "4f7d6c43cfed87e1e5293dc751736efe2d7bc1554946cdc83f4026a575fbbbf8"}};
  const std::array<std::string,10> strong_hash{{"","","",
    "39b8dc3fc8b44765c8e6f1adee04c5b465e555ab791cc42d0d9e810d5b64297c",
    "76650554f5ab120115e47364bbe1822257753e89b72c0731b69edd77c0cd9404",
    "61a69b4844a4f1ec611b88250fb68654626be205c3c9245e30ad341eca745972",
    "3d9ab51665d390c367ced4940d6b80a1233ac8d1a8e158d07134f9abd7bb9ab2",
    "9b96ef048acfddb3b5b1ea29a3964f0987052cdf140f9dc1e150cbcd255c21bf",
    "6900758e8f64444dd2a75450f35e05a1f1f2e00bdf6a77f1395fecd26d689e5a",
    "3073d5ecf5f34345aa5f35c349c51b35f4c244e687db25c16627fcb602b019a1"}};
  const std::string gentourng_expected="89df605922cc574b28688248b7c256d24342cc615f887e89b2d096038970c110";
  std::string gentourng_hash=sha_file(exe);require(gentourng_hash==gentourng_expected,"gentourng binary hash");std::cout<<"gentourng_sha256="<<gentourng_hash<<"\n";
  std::array<Stream,10> all,strong;
  for(int n=3;n<=9;n++){all[n]=read_stream(exe,n,false);strong[n]=read_stream(exe,n,true);require(all[n].labels.size()==all_expected[n],"all count n="+std::to_string(n));require(strong[n].labels.size()==strong_expected[n],"strong count n="+std::to_string(n));require(all[n].digest==all_hash[n],"all stream hash n="+std::to_string(n));require(strong[n].digest==strong_hash[n],"strong stream hash n="+std::to_string(n));for(const auto&s:all[n].labels)all[n].decoded.push_back(decode(s,n));for(const auto&s:strong[n].labels)strong[n].decoded.push_back(decode(s,n));}
  std::vector<Task> tasks;for(int q=3;q<=9;q++)if(!q_only||q==q_only){int r=12-q;for(int s=0;s<SHARDS;s++)tasks.push_back({q,r,s});}
  std::vector<Result> result(tasks.size());std::atomic<size_t> next{0};std::vector<std::thread> pool;
  for(int th=0;th<threads;th++)pool.emplace_back([&]{for(;;){size_t i=next.fetch_add(1);if(i>=tasks.size())break;result[i]=run_task(tasks[i],strong[tasks[i].q],all[tasks[i].r]);}});for(auto&th:pool)th.join();

  uint64_t total=0,strong_rows=0,rfail=0,cfail=0;Extreme global_max,global_min;Sha256 semantic;
  for(int n=3;n<=9;n++){std::cout<<"stream n="<<n<<" all="<<all[n].labels.size()<<" all_sha256="<<all[n].digest<<" strong="<<strong[n].labels.size()<<" strong_sha256="<<strong[n].digest<<"\n";semantic.update(all[n].digest.data(),all[n].digest.size());semantic.update(strong[n].digest.data(),strong[n].digest.size());}
  size_t cursor=0;
  for(int q=3;q<=9;q++)if(!q_only||q==q_only){int r=12-q;uint64_t rows=0,srows=0,rf=0,cf=0;Extreme mx,mn;std::ostringstream shards;
    for(int s=0;s<SHARDS;s++,cursor++){const Result&z=result[cursor];require(z.t.q==q&&z.t.r==r&&z.t.shard==s,"task order");rows+=z.rows;srows+=z.strong;rf+=z.rfail;cf+=z.cfail;add_max(mx,z.mx);add_min(mn,z.mn);if(s)shards<<',';shards<<z.profile_digest;semantic.update(z.profile_digest.data(),z.profile_digest.size());}
    uint64_t expected=uint64_t(strong[q].labels.size())*all[r].labels.size()*q;require(rows==expected,"regime row count");require(srows==rows,"regime strong count");require(rf==0&&cf==0,"regime centrality");static const int64_t min_margin[10]={0,0,0,598,598,468,430,418,388,380};require(mn.p.margin==min_margin[q],"regime minimum margin");total+=rows;strong_rows+=srows;rfail+=rf;cfail+=cf;add_max(global_max,mx);add_min(global_min,mn);
    std::cout<<"regime q="<<q<<" r="<<r<<" rows="<<rows<<" strong="<<srows<<" rational_failures="<<rf<<" coset_failures="<<cf<<"\n";
    std::cout<<"regime_max "<<extreme_text(mx)<<"\n";std::cout<<"regime_min "<<extreme_text(mn)<<"\n";std::cout<<"regime_shard_sha256="<<shards.str()<<"\n";
  }
  uint64_t expected_total=0;for(int q=3;q<=9;q++)if(!q_only||q==q_only)expected_total+=uint64_t(strong[q].labels.size())*all[12-q].labels.size()*q;require(total==expected_total,"global row count");
  std::cout<<"global rows="<<total<<" strong="<<strong_rows<<" rational_failures="<<rfail<<" coset_failures="<<cfail<<"\n";
  std::cout<<"global_max "<<extreme_text(global_max)<<"\n";std::cout<<"global_min "<<extreme_text(global_min)<<"\n";
  std::string sem=semantic.finish();if(!q_only){require(sem=="6516491b7097782ed50e200bdde26b032109f7275e95c5d5cd8c506c3a292829","semantic digest");require(global_max.q==5&&global_max.r==7&&bits11(global_max.p.label)=="1111111101110100111110101111101111110111111111111111101"&&global_max.p.H==6615&&global_max.p.W2==204246&&global_max.p.D4x4==12741213492&&global_max.p.C4==-3629965752&&global_max.p.margin==1494&&global_max.multiplicity==2,"global maximum row");require(global_min.q==9&&global_min.r==3&&bits11(global_min.p.label)=="1011101111111111111111111111111111111111111111101111110"&&global_min.p.H==243&&global_min.p.W2==23742&&global_min.p.D4x4==169212908&&global_min.p.C4==3515272&&global_min.p.margin==380&&global_min.multiplicity==4,"global minimum row");}std::cout<<"semantic_sha256="<<sem<<"\n";
  std::cout<<"status="<<((strong_rows==total&&rfail==0&&cfail==0)?"PASS":"FAIL")<<"\n";
  return (strong_rows==total&&rfail==0&&cfail==0)?0:1;
}catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<"\n";return 2;}}
