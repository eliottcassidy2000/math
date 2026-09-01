// Complete pair-tagged rank-8/rank-9 response atlas for the three
// THM-4303 endpoint-595 carrier-failure rows.  This imports only the audited
// literal-wall primitives and reconstructs C596 from frozen component ledgers.

#define ENDPOINT626_EXCHANGE_MAIN lrc595_hidden_main
#include "../lrc14_size_preserving_response_staircase_thm4300/endpoint_exchange_primitives.cpp"
#undef ENDPOINT626_EXCHANGE_MAIN

#include <chrono>
#include <unordered_map>
#include <unordered_set>

namespace {

constexpr std::array<std::pair<int,int>,3> kPairs = {{{96,595},{100,595},{210,595}}};
constexpr std::array<std::size_t,3> kFailureCounts = {116,13,16};

struct Reply {
  u64 q96lo=0, q96hi=0, q100=0, q210=0;
  auto operator<=>(const Reply&) const = default;
  bool empty() const { return !(q96lo|q96hi|q100|q210); }
};

struct ReplyHash {
  std::size_t operator()(const Reply& x) const noexcept {
    u64 h=UINT64_C(0x9e3779b97f4a7c15);
    for(u64 v:{x.q96lo,x.q96hi,x.q100,x.q210}) {
      v ^= v>>30; v*=UINT64_C(0xbf58476d1ce4e5b9);
      v ^= v>>27; v*=UINT64_C(0x94d049bb133111eb);
      v ^= v>>31; h ^= v + UINT64_C(0x9e3779b97f4a7c15) + (h<<6) + (h>>2);
    }
    return static_cast<std::size_t>(h);
  }
};

struct Stats { u64 count8=0,count9=0; u32 least8=UINT32_MAX,least9=UINT32_MAX; };

std::vector<u32> read_mixed_local(const std::filesystem::path& path,
                                  std::size_t expected) {
  std::ifstream input(path); require(bool(input),"cannot open mixed ledger");
  std::vector<u32> out; std::set<u32> seen; std::string token;
  while(input>>token){u32 m=parse_mask_agent(token);require((std::popcount(m)==8||std::popcount(m)==9)&&seen.insert(m).second,"bad mixed mask");out.push_back(m);}
  require(input.eof()&&out.size()==expected,"mixed ledger count changed"); return out;
}

std::array<std::vector<u32>,3> read_failures_local(const std::filesystem::path& path) {
  std::ifstream input(path); require(bool(input),"cannot open failures");
  std::string line; require(std::getline(input,line)&&line=="q,r,body_hex","failure header changed");
  std::array<std::vector<u32>,3> out; std::array<Fnv,3> ledgers; Fnv global;
  while(std::getline(input,line)){
    if(line.empty())continue; std::replace(line.begin(),line.end(),',',' ');
    std::istringstream row(line); int q=0,r=0; std::string token; row>>q>>r>>token; require(bool(row),"bad failure row");
    int p=q==96&&r==595?0:q==100&&r==595?1:q==210&&r==595?2:-1; require(p>=0,"unexpected failure pair");
    u32 body=parse_mask_agent(token); require(std::popcount(body)==9,"failure body rank changed");
    out[p].push_back(body); ledgers[p].add(body); global.add(q);global.add(r);global.add(body);
  }
  for(int p=0;p<3;++p) require(out[p].size()==kFailureCounts[p]&&std::set<u32>(out[p].begin(),out[p].end()).size()==out[p].size(),"failure count/distinctness changed");
  require(ledgers[0].state==UINT64_C(0xfedacdbff3f31981)&&ledgers[1].state==UINT64_C(0x3ac9ac8b4b9ad93f)&&ledgers[2].state==UINT64_C(0xa6a226f12c168d3a),"failure FNV changed");
  return out;
}

Reply raw_reply(u32 mask,const std::array<std::vector<u32>,3>& failures){
  Reply x;
  for(std::size_t i=0;i<failures[0].size();++i) if(!(mask&failures[0][i])) {if(i<64)x.q96lo|=u64{1}<<i;else x.q96hi|=u64{1}<<(i-64);}
  for(std::size_t i=0;i<failures[1].size();++i) if(!(mask&failures[1][i])) x.q100|=u64{1}<<i;
  for(std::size_t i=0;i<failures[2].size();++i) if(!(mask&failures[2][i])) x.q210|=u64{1}<<i;
  return x;
}

struct TrieNode { std::array<int,30> child; std::array<i64,3> weight{}; TrieNode(){child.fill(-1);} };

struct MassTrie {
  std::vector<TrieNode> nodes{1};
  void insert(u32 mask,const std::array<i64,3>& weight){
    int n=0; for(unsigned bit=0;bit<30;++bit)if(mask&(u32{1}<<bit)){int c=nodes[n].child[bit];if(c<0){c=static_cast<int>(nodes.size());nodes[n].child[bit]=c;nodes.emplace_back();}n=c;}
    for(int p=0;p<3;++p)nodes[n].weight[p]+=weight[p];
  }
  std::array<i64,3> query(u32 mask,u64& visits,u64& max_stack)const{
    struct Frame{int node;u32 remaining;}; std::array<Frame,4096> stack;std::size_t top=0;stack[top++]={0,mask};std::array<i64,3> mass{};u64 local=0;
    while(top){Frame f=stack[--top];++local;for(int p=0;p<3;++p)mass[p]+=nodes[f.node].weight[p];u32 rem=f.remaining;while(rem){unsigned bit=std::countr_zero(rem);rem&=rem-1;int c=nodes[f.node].child[bit];if(c>=0){require(top<stack.size(),"trie stack overflow");stack[top++]={c,rem};max_stack=std::max<u64>(max_stack,top);}}}
    visits+=local;return mass;
  }
};

std::string h16(u64 x){std::ostringstream o;o<<std::hex<<std::setw(16)<<std::setfill('0')<<x;return o.str();}

void add_reply(Fnv& f,const Reply& x){f.add(x.q96lo);f.add(x.q96hi);f.add(x.q100);f.add(x.q210);}

} // namespace

int main(int argc,char**argv){
 try{
  require(argc==9,"usage: atlas BASE8951 ADD45 SUFFIX9 REPAIRS76 DELETE73 ADDITIONS4 FAILURES145 ATLAS_CSV");
  auto failures=read_failures_local(argv[7]);
  std::vector<u32> carrier=build_mixed_carrier(argv[1],argv[2],argv[3]);
  auto repairs=read_mixed_local(argv[4],76);std::set<u32> distinct(carrier.begin(),carrier.end());for(u32 m:repairs){require(distinct.insert(m).second,"repair overlap");carrier.push_back(m);}
  auto deletes=read_mixed_local(argv[5],73);std::set<u32> del(deletes.begin(),deletes.end());std::vector<u32> c596;for(u32 m:carrier)if(!del.contains(m))c596.push_back(m);
  auto additions=read_mixed_local(argv[6],4);for(u32 m:additions)c596.push_back(m);
  require(c596.size()==9019&&masks_fnv_agent(c596)==UINT64_C(0x892fef44a9e6b37e),"C596 changed");std::unordered_set<u32> c596set(c596.begin(),c596.end());

  std::array<Geometry,3> geometry;std::map<u32,std::array<i64,3>> union_weights;
  for(int p=0;p<3;++p){geometry[p]=build_geometry(kPairs[p].first,kPairs[p].second);for(auto [m,w]:geometry[p].classes)union_weights[m][p]+=w;}
  MassTrie trie;for(const auto&[m,w]:union_weights)trie.insert(m,w);
  std::cerr<<"TRIE classes="<<union_weights.size()<<" nodes="<<trie.nodes.size()<<'\n';

  std::unordered_map<Reply,Stats,ReplyHash> atlas;atlas.reserve(1000000);
  std::array<std::array<u64,8>,2> potential_hist{};
  std::array<std::array<u64,3>,2> active{},equal{},pair_responders{};
  std::array<std::array<Fnv,3>,2> pair_responder_fnv;
  std::array<u64,2> universes{},responders{},existing{};std::array<Fnv,2> responder_fnv,response_fnv;
  std::array<std::array<u64,3>,2> whole_family{};std::array<std::array<u32,3>,2> whole_least;for(auto&r:whole_least)r.fill(UINT32_MAX);
  u64 trie_visits=0,max_stack=0;
  const Reply full{UINT64_MAX,((u64{1}<<52)-1),((u64{1}<<13)-1),((u64{1}<<16)-1)};
  for(unsigned rank:{8u,9u}){
    int ri=rank-8;auto start=std::chrono::steady_clock::now();
    for(u32 mask=(u32{1}<<rank)-1;mask<(u32{1}<<30);mask=next_combination(mask)){
      ++universes[ri];Reply raw=raw_reply(mask,failures);unsigned pattern=(!(!(raw.q96lo|raw.q96hi))?1:0)|(!(!raw.q100)?2:0)|(!(!raw.q210)?4:0);++potential_hist[ri][pattern];
      auto mass=trie.query(mask,trie_visits,max_stack);std::array<i128,3> ticks;
      for(int p=0;p<3;++p){ticks[p]=static_cast<i128>(63)*mass[p]-static_cast<i128>(4)*geometry[p].grid;active[ri][p]+=ticks[p]>=0;equal[ri][p]+=ticks[p]==0;}
      Reply reply=raw;if(ticks[0]<0){reply.q96lo=reply.q96hi=0;}if(ticks[1]<0)reply.q100=0;if(ticks[2]<0)reply.q210=0;
      if(reply.empty())continue;++responders[ri];responder_fnv[ri].add(mask);response_fnv[ri].add(mask);add_reply(response_fnv[ri],reply);existing[ri]+=c596set.contains(mask);
      if(reply.q96lo||reply.q96hi){++pair_responders[ri][0];pair_responder_fnv[ri][0].add(mask);}if(reply.q100){++pair_responders[ri][1];pair_responder_fnv[ri][1].add(mask);}if(reply.q210){++pair_responders[ri][2];pair_responder_fnv[ri][2].add(mask);}
      std::array<bool,3> whole={reply.q96lo==full.q96lo&&reply.q96hi==full.q96hi,reply.q100==full.q100,reply.q210==full.q210};for(int p=0;p<3;++p)if(whole[p]){++whole_family[ri][p];whole_least[ri][p]=std::min(whole_least[ri][p],mask);}
      Stats&s=atlas[reply];if(rank==8){++s.count8;s.least8=std::min(s.least8,mask);}else{++s.count9;s.least9=std::min(s.least9,mask);}
    }
    auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();std::cerr<<"RANK "<<rank<<" done seconds="<<seconds<<" responders="<<responders[ri]<<" types="<<atlas.size()<<'\n';
  }
  require(universes==std::array<u64,2>{UINT64_C(5852925),UINT64_C(14307150)},"rank universes changed");require(existing==std::array<u64,2>{0,0},"recorded C596 failure has an existing responder");

  std::vector<std::pair<Reply,Stats>> types(atlas.begin(),atlas.end());std::sort(types.begin(),types.end(),[](const auto&a,const auto&b){return a.first<b.first;});
  Fnv atlas_ledger;std::ofstream csv(argv[8]);require(bool(csv),"cannot create atlas csv");
  csv<<"q96_lo,q96_hi,q100,q210,count8,count9,least8,least9,least8_tick96,least8_tick100,least8_tick210,least9_tick96,least9_tick100,least9_tick210\n";
  for(const auto&[reply,s]:types){add_reply(atlas_ledger,reply);atlas_ledger.add(s.count8);atlas_ledger.add(s.count9);atlas_ledger.add(s.least8);atlas_ledger.add(s.least9);
    csv<<h16(reply.q96lo)<<','<<h16(reply.q96hi)<<','<<h16(reply.q100)<<','<<h16(reply.q210)<<','<<s.count8<<','<<s.count9<<','<<(s.least8==UINT32_MAX?"-":hex8(s.least8))<<','<<(s.least9==UINT32_MAX?"-":hex8(s.least9));
    for(u32 rep:{s.least8,s.least9}){if(rep==UINT32_MAX){csv<<",-,-,-";continue;}u64 v=0,ms=0;auto m=trie.query(rep,v,ms);for(int p=0;p<3;++p)csv<<','<<decimal(static_cast<i128>(63)*m[p]-static_cast<i128>(4)*geometry[p].grid);}csv<<'\n';
  }
  csv.close();
  std::cout<<"LRC595_PAIR_TAGGED_RESPONSE_ATLAS_V1\nC596 9019 FNV 892fef44a9e6b37e\nFAILURES 145 SEGMENTS 96,595:116 100,595:13 210,595:16\n";
  for(int p=0;p<3;++p)std::cout<<"GEOMETRY "<<kPairs[p].first<<','<<kPairs[p].second<<" GRID "<<geometry[p].grid<<" CELLS "<<geometry[p].cells<<" PAIR_TICKS "<<geometry[p].pair_ticks<<" CLASSES "<<geometry[p].classes.size()<<" THRESHOLD_NUMERATOR 4 GRID_DENOMINATOR 63\n";
  std::cout<<"UNION_CLASSES "<<union_weights.size()<<" TRIE_NODES "<<trie.nodes.size()<<" TRIE_VISITS "<<trie_visits<<" MAX_STACK "<<max_stack<<'\n';
  for(int ri=0;ri<2;++ri){int rank=ri+8;std::cout<<"RANK "<<rank<<" UNIVERSE "<<universes[ri]<<" POTENTIAL_HIST";for(int x=0;x<8;++x)std::cout<<' '<<x<<':'<<potential_hist[ri][x];std::cout<<" ACTIVE";for(int p=0;p<3;++p)std::cout<<' '<<kPairs[p].first<<':'<<active[ri][p];std::cout<<" EQUAL";for(int p=0;p<3;++p)std::cout<<' '<<kPairs[p].first<<':'<<equal[ri][p];std::cout<<" RESPONDERS "<<responders[ri]<<" MASK_FNV "<<std::hex<<responder_fnv[ri].state<<" RESPONSE_FNV "<<response_fnv[ri].state<<std::dec<<" EXISTING_C596 "<<existing[ri]<<'\n';for(int p=0;p<3;++p)std::cout<<"PAIR_RESPONDERS RANK "<<rank<<" PAIR "<<kPairs[p].first<<",595 COUNT "<<pair_responders[ri][p]<<" FNV "<<std::hex<<pair_responder_fnv[ri][p].state<<std::dec<<" WHOLE_FAMILY "<<whole_family[ri][p]<<" LEAST "<<(whole_least[ri][p]==UINT32_MAX?"-":hex8(whole_least[ri][p]))<<'\n';}
  std::cout<<"RESPONSE_TYPES "<<types.size()<<" ATLAS_FNV "<<std::hex<<atlas_ledger.state<<std::dec<<" CSV "<<argv[8]<<"\nSCOPE COMPLETE_RANK8_RANK9_FIXED_POOL_PAIR_TAGGED_ACTIVE_RESPONSES_NO_CROSS_GRID_NUMERATOR_COMPARISON_NO_PHYSICAL_ENTRY_NO_LRC14\nVERDICT PASS COMPLETE_PAIR_TAGGED_ATLAS\n";
  return 0;
 }catch(const std::exception&e){std::cerr<<"ATLAS_ERROR "<<e.what()<<'\n';return 1;}
}
