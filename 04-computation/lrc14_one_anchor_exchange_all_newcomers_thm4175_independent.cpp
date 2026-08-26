#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {
using i128 = __int128_t;
using u128 = __uint128_t;
constexpr std::int64_t L = 18241159416480LL;
const std::array<int,3> A{120,126,143};
const std::vector<int> P{8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
  120,126,132,143,145,168,170,176,190,193,240,252,264,286,290};
const std::vector<int> O{8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
  132,145,168,170,176,190,193,240,252,264,286,290};
long long g_comparisons=0;
long long g_equalities=0;

void require(bool x,const std::string& s){if(!x)throw std::runtime_error(s);}
std::string text128(i128 x){
  if(x==0)return "0"; bool neg=x<0;u128 y=neg?u128(-x):u128(x);std::string s;
  while(y){s.push_back(char('0'+y%10));y/=10;}if(neg)s.push_back('-');
  std::reverse(s.begin(),s.end());return s;
}
std::string vec_text(const std::vector<int>& v){
  std::string s="(";for(std::size_t i=0;i<v.size();++i){if(i)s+=',';s+=std::to_string(v[i]);}return s+")";
}
std::vector<int> mask_labels(std::uint32_t m){
  std::vector<int> z;for(int i=0;i<27;++i)if((m>>i)&1U)z.push_back(O[i]);return z;
}
bool safe_mid(std::int64_t left,std::int64_t right,int v){
  const std::int64_t denom=2*L;
  const std::int64_t residue=std::int64_t((i128(v)*(left+right))%denom);
  return i128(14)*residue>=denom && i128(14)*residue<=i128(13)*denom;
}
i128 prefix(std::int64_t tick,int q){
  const i128 prod=i128(q)*tick;const i128 whole=prod/L;const i128 rem=prod-whole*L;
  const i128 scaled=14*rem;i128 partial=0;
  if(scaled<=L)partial=0;else if(scaled>=i128(13)*L)partial=i128(12)*L;else partial=scaled-L;
  return i128(12)*whole*L+partial;
}

struct Deletion{std::uint32_t full;std::uint32_t projected;};
struct Hyper{std::vector<std::uint32_t> edges;int raw=0;int equalities=0;i128 minpos=0,maxneg=0;};

std::vector<Deletion> deletions(int omitted,int t,const std::array<int,291>& pidx,
                                const std::array<int,291>& oidx){
  std::vector<int> allowed{omitted};allowed.insert(allowed.end(),O.begin(),O.end());
  std::vector<Deletion> out;std::vector<int> choice;
  auto rec=[&](auto&& self,int start,int left)->void{
    if(left==0){
      std::uint32_t fm=0,pm=0;
      for(int v:choice){fm|=1U<<pidx[v];if(v!=omitted)pm|=1U<<oidx[v];}
      out.push_back({fm,pm});return;
    }
    for(int i=start;i<=int(allowed.size())-left;++i){choice.push_back(allowed[i]);self(self,i+1,left-1);choice.pop_back();}
  };
  rec(rec,0,t);return out;
}

std::vector<std::uint32_t> minimal_edges(const std::vector<std::uint32_t>& input){
  std::unordered_set<std::uint32_t> all(input.begin(),input.end());
  std::vector<std::uint32_t> rows(all.begin(),all.end());
  std::sort(rows.begin(),rows.end(),[](auto x,auto y){
    if(std::popcount(x)!=std::popcount(y))return std::popcount(x)<std::popcount(y);return x<y;});
  std::vector<std::uint32_t> out;
  for(auto e:rows){
    bool redundant=false;
    for(auto old:out){if((old&e)==old){redundant=true;break;}}
    if(!redundant)out.push_back(e);
  }
  return out;
}

int greedy_matching(const std::vector<std::uint32_t>& edges){
  std::uint32_t used=0;int n=0;for(auto e:edges)if(!(e&used)){used|=e;++n;}return n;
}
struct SearchKey{std::uint32_t chosen;std::uint8_t left;bool operator==(const SearchKey&o)const{return chosen==o.chosen&&left==o.left;}};
struct SearchHash{std::size_t operator()(const SearchKey&k)const{return (std::size_t(k.left)<<32)^k.chosen;}};
bool cover_rec(const std::vector<std::uint32_t>& edges,std::uint32_t chosen,int left,
               std::uint32_t& witness,std::unordered_set<SearchKey,SearchHash>& dead){
  SearchKey key{chosen,std::uint8_t(left)};if(dead.contains(key))return false;
  std::uint32_t uncovered=0,used=0;int matching=0;
  for(auto e:edges){if(e&chosen)continue;if(!uncovered)uncovered=e;if(!(e&used)){used|=e;if(++matching>left){dead.insert(key);return false;}}}
  if(!uncovered){witness=chosen;return true;}if(!left){dead.insert(key);return false;}
  while(uncovered){auto bit=uncovered&(0U-uncovered);if(cover_rec(edges,chosen|bit,left-1,witness,dead))return true;uncovered&=uncovered-1;}
  dead.insert(key);return false;
}
bool cover(const std::vector<std::uint32_t>& edges,int budget,std::uint32_t& witness){
  std::unordered_set<SearchKey,SearchHash> dead;return cover_rec(edges,0,budget,witness,dead);
}
int tau_through8(const std::vector<std::uint32_t>& edges,std::uint32_t& witness){
  if(greedy_matching(edges)>=9)return -1;
  std::uint32_t w=0;if(!cover(edges,8,w))return -1;
  for(int b=0;b<=8;++b)if(cover(edges,b,w)){witness=w;return b;}
  throw std::runtime_error("tau logic");
}

struct Engine{
  std::vector<std::int64_t> walls;
  std::vector<std::uint32_t> failures;
  std::vector<int> group;
  std::vector<std::uint32_t> group_masks;
  std::unordered_map<std::uint32_t,int> group_index;
  std::array<int,291> pidx{},oidx{};
  std::map<std::pair<int,int>,std::vector<Deletion>> bank;
  Engine(){
    pidx.fill(-1);oidx.fill(-1);for(int i=0;i<30;++i)pidx[P[i]]=i;for(int i=0;i<27;++i)oidx[O[i]]=i;
    std::set<std::int64_t> ws{0,L};
    for(int v:P){require(L%(14*v)==0,"lattice divisor");auto unit=L/(14*v);for(int k=0;k<v;++k){ws.insert((14LL*k+1)*unit);ws.insert((14LL*k+13)*unit);}}
    walls.assign(ws.begin(),ws.end());require(walls.size()==7134,"wall count");
    for(std::size_t i=0;i+1<walls.size();++i){
      std::uint32_t f=0;for(int j=0;j<30;++j)if(!safe_mid(walls[i],walls[i+1],P[j]))f|=1U<<j;
      failures.push_back(f);auto it=group_index.find(f);int g;
      if(it==group_index.end()){g=int(group_masks.size());group_index[f]=g;group_masks.push_back(f);}else g=it->second;
      group.push_back(g);
    }
    require(failures.size()==7133,"cell count");
    for(int a:A)for(int t=3;t<=6;++t)bank[{a,t}]=deletions(a,t,pidx,oidx);
  }
  std::vector<i128> masses(int q)const{
    std::vector<i128> m(group_masks.size());i128 old=prefix(walls[0],q);
    for(std::size_t i=0;i<failures.size();++i){i128 cur=prefix(walls[i+1],q);require(cur>=old,"prefix monotone");m[group[i]]+=cur-old;old=cur;}
    return m;
  }
  i128 zeta(const std::vector<i128>& m,std::uint32_t d)const{
    i128 n=0,s=d;while(true){auto it=group_index.find(std::uint32_t(s));if(it!=group_index.end())n+=m[it->second];if(!s)break;s=(s-1)&d;}return n;
  }
  Hyper hyper(int a,int q,int t,const std::vector<i128>& m)const{
    Hyper h;bool haspos=false,hasneg=false;
    for(const auto& d:bank.at({a,t})){
      i128 diff=9*zeta(m,d.full)-i128(8)*q*L;
      if(diff>=0){++h.raw;h.edges.push_back(d.projected);if(diff==0)++h.equalities;else if(!haspos||diff<h.minpos){h.minpos=diff;haspos=true;}}
      else if(!hasneg||diff>h.maxneg){h.maxneg=diff;hasneg=true;}
    }
    h.edges=minimal_edges(h.edges);return h;
  }
  std::pair<i128,int> limit_mass_components(std::uint32_t deleted)const{
    i128 mass=0;int components=0;bool prev=false;
    for(std::size_t i=0;i<failures.size();++i){bool ok=(failures[i]&~deleted)==0;if(ok){mass+=walls[i+1]-walls[i];if(!prev)++components;}prev=ok;}
    return {mass,components};
  }
};

std::vector<int> B3(int a){
 if(a==120)return {3,6,22,24,25,31,46,48,50,55,59,64,70,72,75,83,93,96,100,103,104,105,110,116,122,127,128,130,140,147,153,158,166,172,173,183,186,192,206,208,210,220,244,256,260,270,271,282,294,306,313,320,325,332,346,361,366,372,378,381,383,384,416,420,437,440,462,512,516,519,520,540,550,567,626,650,704,722,744,768,924,1134};
 if(a==126)return {3,6,22,24,25,31,46,48,50,55,59,64,70,72,75,83,93,96,100,103,104,105,110,116,122,127,128,130,140,147,148,153,158,166,172,173,183,186,192,206,208,210,220,244,256,258,260,270,271,282,294,306,313,320,325,332,346,361,366,372,378,381,383,384,416,420,437,440,462,512,516,519,520,540,550,567,626,650,704,722,744,768,924,1134};
 return {3,6,22,24,25,46,48,50,55,64,70,72,75,83,93,96,100,103,105,110,122,128,140,147,153,158,166,172,173,183,186,192,206,208,210,220,256,260,270,282,294,306,313,320,325,332,346,366,372,384,416,420,440,462,512,516,519,520,550,567,704,744,768,924};
}
std::vector<int> B4(int a){
 if(a==120)return {6,24,25,50,96,100,105,128,140,183,192,210,256};
 if(a==126)return {6,24,25,50,96,100,105,128,140,183,192,210,256,366};
 return {6,25,50,96,100,105,128,192,210,256};
}

struct MatchRow{bool delete_anchor;std::vector<int> r;};
std::vector<MatchRow> tail_rows(int a){
 if(a==120)return {{true,{42,88}},{true,{80,252}},{true,{85,170}},{true,{95,190}},{true,{145,290}},{true,{168,286}},{true,{176,193}},{true,{240,264}},{false,{8,16,132}}};
 if(a==126)return {{true,{42,88}},{true,{60,252}},{true,{145,290}},{true,{176,193}},{false,{8,16,286}},{false,{30,85,170}},{false,{63,132,264}},{false,{80,95,190}},{false,{84,168,240}}};
 return {{true,{15,252}},{true,{16,286}},{true,{85,170}},{true,{95,190}},{true,{132,264}},{true,{145,290}},{false,{8,168,176}},{false,{42,63,88}},{false,{80,193,240}}};
}
int cutoff(int a){return a==120?9376:a==126?14511:7883;}

std::uint32_t full_mask(const Engine&e,int a,const MatchRow&r){
 std::uint32_t m=r.delete_anchor?(1U<<e.pidx[a]):0;for(int v:r.r)m|=1U<<e.pidx[v];return m;
}
std::uint32_t projected_mask(const Engine&e,const MatchRow&r){
 std::uint32_t m=0;for(int v:r.r)m|=1U<<e.oidx[v];return m;
}

void check_tail(const Engine&e,int a){
 std::uint32_t used=0;int maximum=0;std::cout<<"TAIL "<<a<<" rows=";
 for(const auto&r:tail_rows(a)){
  auto pm=projected_mask(e,r);require(!(pm&used),"tail projected disjoint");used|=pm;
  auto [mass,c]=e.limit_mass_components(full_mask(e,a,r));i128 slack=54*mass-i128(4)*L;require(slack>0,"tail slack");
  i128 num=i128(54)*c*L,den=i128(7)*slack;int ce=int((num+den-1)/den);maximum=std::max(maximum,ce);
  std::cout<<vec_text(r.r)<<':'<<(r.delete_anchor?'A':'P')<<':'<<c<<':'<<text128(mass)<<'/'<<L<<':'<<ce<<';';
 }
 require(maximum==cutoff(a),"tail cutoff");std::cout<<"cutoff="<<maximum<<"\n";
}

void finite_t3(const Engine&e,int a){
 std::set<int> pool(P.begin(),P.end());std::vector<int> bad;std::map<int,int> tauhist;long long equalities=0;int qualifiers=0;
 for(int q=1;q<cutoff(a);++q){if(pool.contains(q))continue;auto m=e.masses(q);auto h=e.hyper(a,q,3,m);equalities+=h.equalities;std::uint32_t w=0;int tau=tau_through8(h.edges,w);if(tau<0)++qualifiers;else{bad.push_back(q);++tauhist[tau];}}
 g_comparisons+=static_cast<long long>(cutoff(a)-1-static_cast<int>(P.size()))*static_cast<long long>(e.bank.at({a,3}).size());g_equalities+=equalities;
 require(bad==B3(a),"finite t3 failure ledger");require(equalities==0,"finite t3 equality");
 std::cout<<"T3 "<<a<<" qualifiers="<<qualifiers<<" failures="<<bad.size()<<" q="<<vec_text(bad)<<" tau=";
 for(auto [t,n]:tauhist)std::cout<<t<<':'<<n<<',';std::cout<<" equalities="<<equalities<<"\n";
}

void ladder(const Engine&e,int a){
 auto current=B3(a);for(int t=4;t<=6;++t){std::vector<int> next;long long eq=0;
  for(int q:current){auto m=e.masses(q);auto h=e.hyper(a,q,t,m);eq+=h.equalities;g_comparisons+=static_cast<long long>(e.bank.at({a,t}).size());std::uint32_t w=0;int tau=tau_through8(h.edges,w);if(tau>=0)next.push_back(q);
    if(q==50)std::cout<<"HOSTILE "<<a<<" t="<<t<<" raw="<<h.raw<<" minimal="<<h.edges.size()<<" tau="<<tau<<" cover="<<vec_text(mask_labels(w))<<"\n";
  }
  std::vector<int> expected=t==4?B4(a):t==5?std::vector<int>{50}:std::vector<int>{};require(next==expected,"ladder residual");require(eq==0,"ladder equality");
  g_equalities+=eq;std::cout<<"LADDER "<<a<<" t="<<t<<" input="<<current.size()<<" failure="<<vec_text(next)<<" equalities="<<eq<<"\n";current=next;if(current.empty())break;
 }}

void count_uniform_bodies(){
 long long c120=0,c126=0,c143dc=0,c143p=0,c143g2=0;std::vector<int> choose;
 auto gcd_body=[](const std::vector<int>&v){int g=0;for(int x:v)g=std::gcd(g,x);return g;};
 auto dc=[](const std::vector<int>&v){for(int d=2;d<=14;++d){bool ok=false;for(int x:v)if(x%d==0){ok=true;break;}if(!ok)return false;}return true;};
 auto rec=[&](auto&&self,int start,int left)->void{if(!left){
    std::vector<int>b120{126,143};b120.insert(b120.end(),choose.begin(),choose.end());if(dc(b120)){require(gcd_body(b120)==1,"120 primitive");++c120;}
    std::vector<int>b126{120,143};b126.insert(b126.end(),choose.begin(),choose.end());if(dc(b126)){require(gcd_body(b126)==1,"126 primitive");++c126;}
    std::vector<int>b143{120,126};b143.insert(b143.end(),choose.begin(),choose.end());if(dc(b143)){++c143dc;int g=gcd_body(b143);if(g==1)++c143p;else{require(g==2,"143 gcd hostile");++c143g2;}}
    return;}for(int i=start;i<=27-left;++i){choose.push_back(O[i]);self(self,i+1,left-1);choose.pop_back();}};
 rec(rec,0,8);require(c120==2029699&&c126==967956&&c143dc==657800&&c143p==580280&&c143g2==77520,"body counts");
 std::cout<<"BODY_COUNTS 120="<<c120<<" 126="<<c126<<" 143dc="<<c143dc<<" 143primitive="<<c143p<<" 143gcd2="<<c143g2<<" total_primitive="<<(c120+c126+c143p)<<"\n";
}
} // namespace

int main(){
 Engine e;std::cout<<"LRC14_ONE_ANCHOR_EXCHANGE_INDEPENDENT_CPP walls="<<e.walls.size()<<" cells="<<e.failures.size()<<" atoms="<<e.group_masks.size()<<"\n";
 for(int a:A)check_tail(e,a);
 for(int a:A)finite_t3(e,a);
 for(int a:A)ladder(e,a);
 count_uniform_bodies();require(g_comparisons==113249682,"comparison total");require(g_equalities==0,"global equality total");
 std::cout<<"AUDIT comparisons="<<g_comparisons<<" threshold_equalities="<<g_equalities<<"\nPASS\n";
}
