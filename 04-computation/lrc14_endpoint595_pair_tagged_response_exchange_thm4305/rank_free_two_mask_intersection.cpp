// Independent rank-free two-mask audit for THM-4303's 145 tagged failures.
// Formal concepts are enumerated as the closure system of row intersections,
// not by NextClosure.  Unlike rank-8/9 response code, this geometry retains
// every literal failure class (including ranks 10 through 25).

#define main lrc595_rankfree_hidden_main
#include "../lrc14_mixed_rank_depth_recursive_signatures_thm4296/r632_detached_hostile_survivor.cpp"
#undef main

namespace {
constexpr u32 kFull=(u32{1}<<30)-1;
constexpr std::array<int,3> kQ={96,100,210};

struct FullGeometry {i64 grid=0;u64 cells=0;i64 pair_ticks=0;std::vector<std::pair<u32,i64>> classes;};

FullGeometry full_geometry(int q,int r){
 i64 grid=fixed_grid();grid=checked_lcm(grid,14LL*q);grid=checked_lcm(grid,14LL*r);std::vector<i64>walls={0,grid};
 auto add=[&](int speed){require(grid%(14LL*speed)==0,"wall grid");i64 unit=grid/(14LL*speed);for(int t=0;t<speed;++t){walls.push_back((14LL*t+1)*unit);walls.push_back((14LL*t+13)*unit);}};
 for(int s:kPool)add(s);add(q);add(r);std::sort(walls.begin(),walls.end());walls.erase(std::unique(walls.begin(),walls.end()),walls.end());
 std::map<u32,i64> by;FullGeometry g;g.grid=grid;g.cells=walls.size()-1;
 for(std::size_t i=1;i<walls.size();++i){i64 l=walls[i-1],rr=walls[i];if(!safe_midpoint(q,grid,l,rr)||!safe_midpoint(r,grid,l,rr))continue;u32 f=0;for(unsigned v=0;v<kPool.size();++v)if(!safe_midpoint(kPool[v],grid,l,rr))f|=u32{1}<<v;i64 w=rr-l;g.pair_ticks+=w;by[f]+=w;}
 for(auto e:by)g.classes.push_back(e);return g;
}

i128 full_ticks(const FullGeometry&g,u32 mask){i64 mass=0;for(auto[f,w]:g.classes)if((f&~mask)==0)mass+=w;return static_cast<i128>(63)*mass-static_cast<i128>(4)*g.grid;}

std::array<std::vector<u32>,3> read_tagged_full(const std::filesystem::path&p){
 std::ifstream in(p);require(bool(in),"open failures");std::string line;require(std::getline(in,line)&&line=="q,r,body_hex","header");std::array<std::vector<u32>,3> out;std::array<Fnv,3> f;
 while(std::getline(in,line)){if(line.empty())continue;std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);int q,r;std::string h;row>>q>>r>>h;require(bool(row)&&r==595,"row");int j=q==96?0:q==100?1:q==210?2:-1;require(j>=0,"pair");u32 b=std::stoul(h,nullptr,16);require(std::popcount(b)==9,"rank");out[j].push_back(b);f[j].add(b);}
 require(out[0].size()==116&&out[1].size()==13&&out[2].size()==16,"counts");require(f[0].state==UINT64_C(0xfedacdbff3f31981)&&f[1].state==UINT64_C(0x3ac9ac8b4b9ad93f)&&f[2].state==UINT64_C(0xa6a226f12c168d3a),"FNV");return out;
}

u32 family_union(const std::vector<u32>&v){u32 x=0;for(u32 b:v)x|=b;return x;}
u32 derive(u32 x,const std::array<u32,30>&rel){u32 y=kFull;while(x){unsigned b=std::countr_zero(x);x&=x-1;y&=rel[b];}return y;}
std::array<bool,3> sig(u32 m,const std::array<FullGeometry,3>&g){return {full_ticks(g[0],m)>=0,full_ticks(g[1],m)>=0,full_ticks(g[2],m)>=0};}
std::string sigstr(const std::array<bool,3>&s){return std::string{s[0]?'1':'0',s[1]?'1':'0',s[2]?'1':'0'};}
std::array<unsigned,3> cover(u32 a,u32 b,const std::array<std::vector<u32>,3>&f,const std::array<FullGeometry,3>&g){auto sa=sig(a,g),sb=sig(b,g);std::array<unsigned,3> c{};for(int p=0;p<3;++p)for(u32 body:f[p])c[p]+=((sa[p]&&!(a&body))||(sb[p]&&!(b&body)));return c;}
bool exact(const std::array<unsigned,3>&c){return c==std::array<unsigned,3>{116,13,16};}
}

int main(int argc,char**argv){
 try{require(argc==2,"usage: audit FAILURE_CSV");auto failures=read_tagged_full(argv[1]);std::array<FullGeometry,3> g={full_geometry(96,595),full_geometry(100,595),full_geometry(210,595)};
  std::array<u32,30> co{};for(const auto&fam:failures)for(u32 body:fam){u32 z=body;while(z){unsigned b=std::countr_zero(z);z&=z-1;co[b]|=body;}}
  std::array<u32,30> rel;for(unsigned b=0;b<30;++b)rel[b]=kFull^co[b];
  std::set<u32> intents={kFull};for(unsigned b=0;b<30;++b){std::vector<u32> old(intents.begin(),intents.end());for(u32 y:old)intents.insert(y&rel[b]);}
  std::cout<<"LRC595_TWO_MASK_INTERSECTION_AUDIT_V1\n";for(int p=0;p<3;++p){u32 allowed=kFull^family_union(failures[p]);auto s=sig(allowed,g);std::cout<<"WHOLE_MAX PAIR "<<kQ[p]<<",595 MASK "<<hex8(allowed)<<" RANK "<<std::popcount(allowed)<<" SIG "<<sigstr(s)<<" OWN_TICKS "<<decimal(full_ticks(g[p],allowed))<<" GRID "<<g[p].grid<<" FULL_CLASSES "<<g[p].classes.size()<<'\n';}
  u64 checked=0,active96210=0,activeall=0,exactn=0;Fnv ledger,activeledger;std::array<unsigned,3> best{};u32 besta=0,bestb=0;
  for(u32 y:intents){u32 x=derive(y,rel);require(derive(x,rel)==y,"intent not closed");++checked;ledger.add(x);ledger.add(y);auto sx=sig(x,g),sy=sig(y,g);if(sx[0]&&sx[2]&&sy[0]&&sy[2]){++active96210;activeledger.add(x);activeledger.add(y);activeall+=sx[1]&&sy[1];auto c=cover(x,y,failures,g);if(std::accumulate(c.begin(),c.end(),0u)>std::accumulate(best.begin(),best.end(),0u)){best=c;besta=x;bestb=y;}if(exact(c)){++exactn;std::cout<<"EXACT_CONCEPT "<<hex8(x)<<' '<<hex8(y)<<" SIG "<<sigstr(sx)<<' '<<sigstr(sy)<<'\n';}}}
  u32 allowed100=kFull^family_union(failures[1]);u64 subsets=0,candidates=0,onesided=0;Fnv oneledger;for(u32 x=allowed100;;x=(x-1)&allowed100){++subsets;auto sx=sig(x,g);if(sx[0]&&sx[1]&&sx[2]){u32 y=derive(x,rel);auto sy=sig(y,g);if(sy[0]&&sy[2]){++candidates;oneledger.add(x);oneledger.add(y);if(exact(cover(x,y,failures,g))){++onesided;std::cout<<"EXACT_ONE_SIDED "<<hex8(x)<<' '<<hex8(y)<<" SIG "<<sigstr(sx)<<' '<<sigstr(sy)<<'\n';}}}if(x==0)break;}
  std::cout<<"INTERSECTION_CONCEPTS "<<checked<<" SORTED_FNV "<<std::hex<<ledger.state<<std::dec<<" ACTIVE_96_210 "<<active96210<<" ACTIVE_FNV "<<std::hex<<activeledger.state<<std::dec<<" ACTIVE_ALL "<<activeall<<" EXACT "<<exactn<<'\n';
  std::cout<<"BEST_ACTIVE_CONCEPT A "<<hex8(besta)<<" B "<<hex8(bestb)<<" COVER "<<best[0]<<','<<best[1]<<','<<best[2]<<" TOTAL "<<std::accumulate(best.begin(),best.end(),0u)<<'\n';
  std::cout<<"ONE_SIDED_Q100 SUBSETS "<<subsets<<" CANDIDATES "<<candidates<<" FNV "<<std::hex<<oneledger.state<<std::dec<<" EXACT "<<onesided<<'\n';
  std::cout<<"REDUCTION if WHOLE_MAX q96 and q210 are inactive, both masks in any two-cover must be active on those rows; every no-cooccurrence pair extends monotonically to a formal concept; q100 both-active concepts and one-active subsets are audited separately\n";
  std::cout<<"SCOPE RANK_FREE_FULL_LITERAL_CLASSES_EXACT_TWO_MASK_OBSTRUCTION_ON_FIXED_145_FAILURES_ONLY_NO_PHYSICAL_ENTRY_NO_LRC14\nVERDICT "<<(exactn==0&&onesided==0?"PASS_NO_TWO_MASK_COVER":"FAIL_TWO_MASK_COVER_FOUND")<<'\n';return exactn||onesided;
 }catch(const std::exception&e){std::cerr<<"INDEPENDENT_TWO_MASK_ERROR "<<e.what()<<'\n';return 1;}
}
