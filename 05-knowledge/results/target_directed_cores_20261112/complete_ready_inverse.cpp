// Complete target-specific core search using P_C(2) <= 2*target.
// The inherited unique-start grammar is used, with safe positive pruning.
#define main finite_bank_main
#include "target_inverse.cpp"
#undef main
#include <array>
#include <limits>
struct Node {std::array<W,6> z; int parent_order,parent_index; I word;};
struct Coeff {W l[3],d[3],b[3];};
I ipow(I x,int n) {I r=1;while(n--)r*=x;return r;}
Coeff coefficients(I word,int t) {
  Coeff c;
  for(int m=1;m<=3;++m) {
    W L=ipow(m-1,t),D=L,vec[6]={1,0,0,0,0,0};
    for(int i=0;i<t;++i) {
      int bit=(word>>i)&1;
      if(bit)L+=W(ipow(m,i))*ipow(m-1,t-i-1);
      else D+=W(ipow(m-1,i))*ipow(m,t-i-1);
      W next[6]={0,0,0,0,0,0};
      for(int last=0;last<3;++last)for(int started=0;started<2;++started) {
        W a=vec[2*last+started];
        next[2*last+started]+=(m-1)*a;
        next[2*(bit?1:2)+started]+=a;
        if(started||!bit)next[2*last+1]+=a;
      }
      std::copy(next,next+6,vec);
    }
    c.l[m-1]=L;c.d[m-1]=D;c.b[m-1]=vec[0]+vec[1]+vec[2]+vec[3];
  }
  return c;
}
std::array<W,6> advance(const std::array<W,6>&z,const Coeff &c,int t) {
  W h=z[0],e=z[1],s=z[2],p=z[3],f=z[4],v=z[5];
  return {c.l[1]*s+c.b[0]*h-c.d[0]*e,c.l[0]*h-(t?0:e),
          W(ipow(2,t))*s+c.d[1]*h,
          2*(c.l[2]*v+c.b[1]*p-c.d[1]*f),c.l[1]*p-f,
          W(ipow(3,t))*v+c.d[2]*p};
}
void write_address(std::ostream&out,const std::vector<std::vector<Node>>&bank,int n,int index) {
  if(n==1)return;
  const auto &v=bank[n][index];
  write_address(out,bank,v.parent_order,v.parent_index);
  int t=n-v.parent_order-1;out<<"\"";
  for(int i=0;i<t;++i)out<<((v.word>>i)&1);
  out<<"\", ";
}
int main(int argc,char**argv) {
  if(argc!=3){std::cerr<<"usage: complete_ready_inverse TARGET receipt.json\n";return 2;}
  I N=std::stoll(argv[1]);if(N<3||N>1000000000000LL)return 3;
  I bound=2*N;int cap=0;while((I(1)<<(cap+1))<=bound)++cap;
  std::vector<std::vector<Node>> bank(cap+1);
  bank[1].push_back({{1,1,1,2,1,1},0,0,0});
  std::vector<std::vector<Coeff>> cache(cap+1);
  std::ofstream out(argv[2]);out<<"{\n\"target\": "<<N<<", \"core_p2_bound\": "<<bound<<", \"core_order_bound\": "<<cap<<",\n\"layers\": [\n";
  bool first=true;uint64_t total=0;visits=0;
  for(int n=1;n<=cap;++n) {
    if(n==2)continue;
    uint64_t raw=0,pruned=0;
    if(n>1)for(int m=1;m<n;++m) {
      if(bank[m].empty())continue;
      int t=n-m-1;
      W lower=t?(W(ipow(2,t))+ipow(3,t-1)+1):2;
      bool any=false;for(const auto&v:bank[m])if(lower*v.z[3]<=bound){any=true;break;}
      if(!any)continue;
      I words=t?(I(1)<<(t-1)):1;
      if(cache[t].empty()) {
        cache[t].reserve(words);
        for(I j=0;j<words;++j)cache[t].push_back(coefficients(t?j+(I(1)<<(t-1)):0,t));
      }
      for(I j=0;j<words;++j)for(size_t index=0;index<bank[m].size();++index) {
        const auto &v=bank[m][index];
        if(lower*v.z[3]>bound)continue;
        ++raw;auto z=advance(v.z,cache[t][j],t);
        if(z[3]>bound){++pruned;continue;}
        if(z[3]<2*v.z[3]){std::cerr<<"monotonicity failure\n";return 5;}
        bank[n].push_back({z,m,int(index),t?j+(I(1)<<(t-1)):0});
      }
    }
    total+=bank[n].size();
    if(!first)out<<",\n";first=false;
    out<<"{\"order\": "<<n<<", \"retained_cores\": "<<bank[n].size()<<", \"evaluated\": "<<raw<<", \"p2_pruned\": "<<pruned<<"}";out.flush();
    std::cout<<"order="<<n<<" retained="<<bank[n].size()<<" evaluated="<<raw<<" pruned="<<pruned<<std::endl;
    for(size_t i=0;i<bank[n].size();++i) {
      auto z=bank[n][i].z;I h=I(z[0]),G=I(z[2]+z[3]),e=I(z[1]),d=I(z[2]+z[3]-z[1]-z[4]-1),hit=-1;
      if(inverse(N,h,G,e,d,1,0,hit)&&(n>1||hit>0)) {
        out<<"\n],\n\"found\": true, \"core_order\": "<<n<<", \"filtered_index\": "<<i<<", \"core_state\": [";
        for(int j=0;j<6;++j){if(j)out<<", ";if(z[j]>std::numeric_limits<I>::max())return 6;out<<I(z[j]);}
        out<<"], \"word_integer\": "<<hit<<", \"address\": [";
        write_address(out,bank,n,int(i));out<<"\"";
        for(I x=hit;x;x>>=1)out<<(x&1);
        out<<"\"], \"inverse_visits\": "<<visits<<"\n}\n";
        std::cout<<"FOUND target="<<N<<" core_order="<<n<<" index="<<i<<" word="<<hit<<" cores="<<total<<std::endl;return 0;
      }
    }
  }
  out<<"\n], \"found\": false, \"complete_core_universe\": true, \"inverse_visits\": "<<visits<<"\n}\n";
  std::cout<<"EXCLUDED from complete unit-ready image; cores="<<total<<std::endl;
}
