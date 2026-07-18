#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {
constexpr int P=13, C=27;
constexpr uint32_t FULL=(uint32_t{1}<<C)-1;
constexpr std::array<int,4> O{1,3,9,27};
using Support=std::array<uint8_t,6>;
using Word=std::array<uint8_t,6>;

[[noreturn]] void fail(const std::string&m){std::cerr<<"FAIL "<<m<<'\n';std::exit(1);}
void req(bool b,const std::string&m){if(!b)fail(m);}
int inv13(int x){for(int y=1;y<13;y++)if(x*y%13==1)return y;fail("inv");}
int ctr(int x,int m){x%=m;if(x<0)x+=m;return 2*x>m?x-m:x;}
std::vector<uint8_t> units(int d){if(d==1)return {0};std::vector<uint8_t>v;for(int u=1;u<d;u++)if(std::gcd(u,d)==1)v.push_back(u);return v;}
int crt(int r,int d,int u){for(int x=0;x<P*d;x++)if(x%P==d*r%P&&x%d==u%d)return x;fail("crt");}
uint32_t mask(int r,int d,int u){int b=crt(r,d,u);uint32_t m=0;for(int s=0;s<C;s++){int x=ctr(b*(1+P*s),P*d);if(-d<x&&x<=d)m|=uint32_t{1}<<s;}return m;}
int lc(int a,int b){return std::lcm(a,b);}
bool hered(const Word&w){for(int z=0;z<6;z++){int x=1;for(int i=0;i<6;i++)if(i!=z)x=lc(x,O[w[i]]);if(x!=C)return false;}return true;}

std::string key_for(const Support&s,const Word&w,int owner){
 std::vector<std::pair<int,int>> v; int iv=inv13(owner);
 for(int i=0;i<6;i++)v.push_back({s[i]*iv%P,w[i]});
 std::sort(v.begin(),v.end());std::string k;
 for(auto [r,d]:v){k.push_back(char(r));k.push_back(char(d));}return k;
}
std::string desc_key(const std::string&k){std::ostringstream o;for(int i=0;i<12;i+=2)o<<int((uint8_t)k[i])<<':'<<O[(uint8_t)k[i+1]]<<(i==10?"":" ");return o.str();}
struct Sum{int max=0;bool full=false;size_t bank=0;std::array<int,3> bestrow{};};
struct FlagState { uint16_t full=0; uint32_t hits=0; auto operator<=>(const FlagState&) const = default; };
uint16_t fibre_signature(uint32_t m){uint16_t z=0;for(int a=0;a<9;a++){int n=0;for(int s=a;s<C;s+=9)n+=(m>>s)&1;if(n)z|=uint16_t{1}<<a;}return z;}
uint16_t full_fibres(uint32_t m){uint16_t z=0;for(int a=0;a<9;a++){int n=0;for(int s=a;s<C;s+=9)n+=(m>>s)&1;if(n==3)z|=uint16_t{1}<<a;}return z;}
uint32_t add_hits(uint32_t h,uint16_t sig){for(int a=0;a<9;a++)if(sig>>a&1){int n=(h>>(2*a))&3;if(n<3)h+=uint32_t{1}<<(2*a);}return h;}
}

int main(){
 std::array<std::array<std::vector<uint32_t>,4>,13> fib;
 std::array<std::array<int,4>,13> card{};
 for(int r=1;r<13;r++)for(int d=0;d<4;d++){
  for(int u:units(O[d]))fib[r][d].push_back(mask(r,O[d],u));
  std::sort(fib[r][d].begin(),fib[r][d].end());fib[r][d].erase(std::unique(fib[r][d].begin(),fib[r][d].end()),fib[r][d].end());
  card[r][d]=std::popcount(fib[r][d][0]);for(auto m:fib[r][d])req(std::popcount(m)==card[r][d],"card");
 }
 std::vector<Word> words;uint64_t weighted=0;
 for(int n=0;n<4096;n++){int x=n;Word w{};uint64_t f=1;for(int i=0;i<6;i++){w[i]=x%4;x/=4;f*=units(O[w[i]]).size();}if(hered(w)){words.push_back(w);weighted+=f;}}
 std::cout<<"words "<<words.size()<<" weighted "<<weighted<<'\n';
 struct Row{Support s;Word w;std::array<uint8_t,6> cap;};std::vector<Row> rows;
 std::array<int,4> mult{};std::map<std::array<int,4>,int> prof;
 for(int a=1;a<=7;a++)for(int b=a+1;b<=8;b++)for(int c=b+1;c<=9;c++)for(int d=c+1;d<=10;d++)for(int e=d+1;e<=11;e++)for(int f=e+1;f<=12;f++){
  Support s{(uint8_t)a,(uint8_t)b,(uint8_t)c,(uint8_t)d,(uint8_t)e,(uint8_t)f};
  for(auto&w:words){Row row{s,w,{}};bool ok=true;for(int oi=0;oi<6;oi++){int iv=inv13(s[oi]),z=0;for(int i=0;i<6;i++)z+=card[s[i]*iv%P][w[i]];row.cap[oi]=z;if(z<C)ok=false;}if(ok){rows.push_back(row);std::array<int,4>p{};for(auto q:w)p[q]++;prof[p]++;}}
 }
 std::set<Support> supports;for(auto&r:rows)supports.insert(r.s);
 std::cout<<"rows "<<rows.size()<<" supports "<<supports.size()<<" profiles";for(auto&[p,n]:prof)std::cout<<" ["<<p[0]<<','<<p[1]<<','<<p[2]<<','<<p[3]<<"]="<<n;std::cout<<'\n';
 std::unordered_map<std::string,Sum> memo;
 std::map<int,int> feas_hist,max_hist,bank_hist,keys_by_full;
 std::map<std::pair<int,int>,int> order_max;
 std::map<std::array<int,4>,std::map<int,int>> profile_feas;
 std::map<std::pair<int,int>,int> self_order_feas;
 std::map<std::pair<int,int>,int> self_order_flag;
 std::map<std::pair<int,int>,int> self_order_flagmax;
 std::map<int,std::set<std::string>> keys_by_order;
 std::map<std::pair<std::array<int,4>,int>,int> d27_profile_flagmax;
 std::map<std::string,int> key_mult;
 int ri=0;
 for(auto&row:rows){int nf=0;std::array<int,6> yes{},mx{},bn{};for(int oi=0;oi<6;oi++){
   auto k=key_for(row.s,row.w,row.s[oi]);key_mult[k]++;
   bool new_owner_key=keys_by_order[O[row.w[oi]]].insert(k).second;
   if(!memo.contains(k)){
    std::vector<std::pair<int,int>> v;for(int i=0;i<12;i+=2)v.push_back({(uint8_t)k[i],(uint8_t)k[i+1]});
    std::unordered_set<uint32_t> reach{0};
    for(auto [ratio,di]:v){std::unordered_set<uint32_t> next;next.reserve(reach.size()*fib[ratio][di].size());for(auto x:reach)for(auto m:fib[ratio][di])next.insert(x|m);reach=std::move(next);}
    Sum z;z.bank=reach.size();for(auto x:reach){z.max=std::max(z.max,std::popcount(x));z.full|=x==FULL;std::array<int,3>cnt{};for(int s=0;s<C;s++)if(x>>s&1)cnt[s%3]++;if(std::min(cnt[0],std::min(cnt[1],cnt[2]))>std::min(z.bestrow[0],std::min(z.bestrow[1],z.bestrow[2])))z.bestrow=cnt;}
    memo.emplace(k,z);
   }
   auto z=memo.at(k);yes[oi]=z.full;mx[oi]=z.max;bn[oi]=z.bank;nf+=z.full;max_hist[z.max]++;order_max[{O[row.w[oi]],z.max}]++;bank_hist[z.bank]++;keys_by_full[z.full]++;self_order_feas[{O[row.w[oi]],z.full}]++;
   std::set<FlagState> flags{{}};for(int j=0;j<12;j+=2){int ratio=(uint8_t)k[j],di=(uint8_t)k[j+1];std::set<FlagState> nx;for(auto st:flags)for(auto m:fib[ratio][di]){auto q=st;if(di<3)q.full|=full_fibres(m);else q.hits=add_hits(q.hits,fibre_signature(m));nx.insert(q);}flags=std::move(nx);}bool flagok=false;int flagmax=0;for(auto st:flags){bool q=true;int score=0;for(int a=0;a<9;a++){int h=(st.hits>>(2*a))&3;q&=(st.full>>a&1)||(h==3);score+=(st.full>>a&1)?3:h;}flagok|=q;flagmax=std::max(flagmax,score);}req(flagok==z.full,"flag quotient not exact");self_order_flag[{O[row.w[oi]],flagok}]++;self_order_flagmax[{O[row.w[oi]],flagmax}]++;if(O[row.w[oi]]==27){std::array<int,4> pp{};for(auto q:row.w)pp[q]++;d27_profile_flagmax[{pp,flagmax}]++;}if(new_owner_key&&O[row.w[oi]]==27&&flagmax>=26)std::cout<<"D27-FLAG26 "<<desc_key(k)<<" actual="<<z.max<<" bank="<<z.bank<<'\n';
  }
  std::array<int,4>p{};for(auto q:row.w)p[q]++;profile_feas[p][nf]++;feas_hist[nf]++;
  if(nf>=1){std::cout<<"ROW "<<ri<<" F"<<nf<<" ";for(int i=0;i<6;i++)std::cout<<int(row.s[i])<<':'<<O[row.w[i]]<<(yes[i]?"+":"-")<<"("<<mx[i]<<","<<int(row.cap[i])<<","<<bn[i]<<") ";std::cout<<'\n';}
  ri++;
 }
 std::cout<<"memo "<<memo.size()<<" feasible";for(auto[x,n]:feas_hist)std::cout<<' '<<x<<':'<<n;std::cout<<" max";for(auto[x,n]:max_hist)std::cout<<' '<<x<<':'<<n;std::cout<<'\n';
 std::cout<<"profile-feasible\n";for(auto&[p,h]:profile_feas){std::cout<<p[0]<<p[1]<<p[2]<<p[3];for(auto[x,n]:h)std::cout<<' '<<x<<':'<<n;std::cout<<'\n';}
 std::cout<<"self-order-feas";for(auto[x,n]:self_order_feas)std::cout<<' '<<x.first<<','<<x.second<<':'<<n;std::cout<<'\n';
 std::cout<<"self-order-max";for(auto[x,n]:order_max)std::cout<<' '<<x.first<<','<<x.second<<':'<<n;std::cout<<'\n';
 std::cout<<"self-order-flag";for(auto[x,n]:self_order_flag)std::cout<<' '<<x.first<<','<<x.second<<':'<<n;std::cout<<'\n';
 std::cout<<"self-order-flagmax";for(auto[x,n]:self_order_flagmax)std::cout<<' '<<x.first<<','<<x.second<<':'<<n;std::cout<<'\n';
 std::cout<<"owner-key-count";for(auto&[d,s]:keys_by_order)std::cout<<' '<<d<<':'<<s.size();std::cout<<'\n';
 std::cout<<"d27-profile-flagmax";for(auto&[x,n]:d27_profile_flagmax){auto p=x.first;std::cout<<" ["<<p[0]<<p[1]<<p[2]<<p[3]<<","<<x.second<<"]="<<n;}std::cout<<'\n';
}
