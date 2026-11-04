// Independent tournament census / prime / endpoint audit, orders at most 9.
// Does not call or include the producer's code. Canonicalization uses only
// degree buckets, all bucket permutations, and row-major (not source) encoding.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
using U = uint64_t;
using R = std::array<unsigned,9>;
void require(bool b,const std::string&s){if(!b)throw std::runtime_error(s);}
unsigned pc(unsigned x){return __builtin_popcount(x);}
R decode(int n,U code){R r{};int k=n*(n-1)/2-1;for(int j=1;j<n;++j)for(int i=0;i<j;++i,--k){if((code>>k)&1)r[i]|=1u<<j;else r[j]|=1u<<i;}return r;}
R opposite(int n,const R&r){R s{};for(int i=0;i<n;++i)s[i]=((1u<<n)-1)^r[i]^(1u<<i);return s;}
struct Canon {U key=~U(0), aut=0, visited=0;};
Canon canonical(int n,const R&r){
 std::array<int,9> p{};std::iota(p.begin(),p.begin()+n,0);
 std::sort(p.begin(),p.begin()+n,[&](int a,int b){return std::pair(pc(r[a]),a)<std::pair(pc(r[b]),b);});
 std::vector<std::pair<int,int>> blocks;
 for(int a=0;a<n;){int b=a+1;while(b<n&&pc(r[p[b]])==pc(r[p[a]]))++b;blocks.emplace_back(a,b);a=b;}
 U expected_permutations=1;
 for(auto [a,b]:blocks) for(int k=2;k<=b-a;++k) expected_permutations*=k;
 Canon z;
 auto visit=[&](auto&&self,size_t c)->void {
  if(c==blocks.size()){
   U key=0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j)key=(key<<1)|((r[p[i]]>>p[j])&1u);
   ++z.visited;if(key<z.key){z.key=key;z.aut=1;}else if(key==z.key)++z.aut;return;
  }
  auto [a,b]=blocks[c];
  do{self(self,c+1);}while(std::next_permutation(p.begin()+a,p.begin()+b));
 };
 visit(visit,0);
 require(z.visited==expected_permutations,"degree bucket permutation coverage mismatch");
 return z;
}
// A module containing a seed pair must also contain every outside vertex
// distinguishing the current set. Closure reaches a proper set iff a proper
// nontrivial module containing that pair exists.
unsigned proper_module(int n,const R&r){
 unsigned F=(1u<<n)-1;
 for(int a=0;a<n;++a)for(int b=a+1;b<n;++b){
  unsigned m=(1u<<a)|(1u<<b);
  for(;;){unsigned old=m;for(int x=0;x<n;++x)if(!(m>>x&1u)){unsigned t=r[x]&m;if(t&&t!=m)m|=1u<<x;}if(m==F)break;if(m==old)return m;}
 }
 return 0;
}
struct Paths{std::array<U,9> starts{},ends{};U visited=0;};
Paths literal_paths(int n,const R&r){
 Paths p;unsigned F=(1u<<n)-1;
 auto visit=[&](auto&&self,int first,int last,unsigned used)->void{
  ++p.visited;if(used==F){++p.starts[first];++p.ends[last];return;}
  unsigned choices=r[last]&~used;while(choices){int next=__builtin_ctz(choices);choices&=choices-1;self(self,first,next,used|(1u<<next));}
 };
 for(int i=0;i<n;++i) visit(visit,i,i,1u<<i);
 return p;
}
std::string list(const std::array<U,9>&a,int n){std::ostringstream s;for(int i=0;i<n;++i){if(i)s<<',';s<<a[i];}return s.str();}
std::string rows(const R&a,int n){std::ostringstream s;for(int i=0;i<n;++i){if(i)s<<',';s<<a[i];}return s.str();}
int main(int argc,char**argv){try{
 bool compact=(argc==4&&std::string(argv[1])=="--compact");
 require(argc==3||compact,"usage: audit DATA_DIR OUT_DIR OR audit --compact CENSUS_TSV OUT_DIR");
 std::string dir=argv[compact?2:1],out=argv[compact?3:2];
 std::ofstream summary(out+"/summary.tsv"),audit(out+"/all_rows.tsv"),candidates(out+"/candidates.tsv");
 require(bool(summary)&&bool(audit)&&bool(candidates),"cannot open output files");
 summary<<"n\tmerged\tordinary\tself_converse\tlabeled_mass\texpected_mass\tpermutations\tprime_merged\tp_nonzero\tp_zero\tprime_ordinary\tp_zero_ordinary\tpath_prefixes\n";
 audit<<"n\tnode\tsource_code\tkey\tconverse_key\taut\tself_converse\tproper_module\tprime\tp_mask\tstarts\tends\n";
 candidates<<"node\tcode\trows\tstarts\tends\n";
 for(int n=compact?9:2;n<=9;++n){
  std::ifstream in(compact?dir:dir+"/tower_"+std::to_string(n)+".txt");require(bool(in),"input missing");std::string line;std::getline(in,line);
  std::map<U,int> owners;U factorial=1;for(int i=2;i<=n;++i)factorial*=i;
  U mass=0,perms=0,path_nodes=0;unsigned merged=0,sc_count=0,primes=0,p_nonzero=0,p_zero=0,prime_ord=0,zero_ord=0;
  while(std::getline(in,line)){
   require(!line.empty(),"blank input row");std::istringstream head(line);unsigned id,sc=0,aut=0;U H=0,tilings=0,source=0;std::string scores;
   if(compact)head>>id>>source;else{head>>id>>sc>>aut>>H>>tilings>>scores;source=std::stoull(line.substr(line.rfind('|')+1));}
   require(bool(head),"malformed row");require(id==merged,"nonsequential ID");require(sc<=1,"bad sc");
   require(source<(U(1)<<(n*(n-1)/2)),"code out of range");R r=decode(n,source),op=opposite(n,r);
   for(int i=0;i<n;++i){require(!(r[i]>>i&1u),"loop");for(int j=i+1;j<n;++j)require(((r[i]>>j)&1u)+((r[j]>>i)&1u)==1,"not tournament");}
   Canon c=canonical(n,r),d=canonical(n,op);perms+=c.visited+d.visited;
   require(c.aut==d.aut,"converse aut mismatch");
   if(compact)sc=(c.key==d.key);else{require(c.aut==aut,"source aut mismatch");require((c.key==d.key)==bool(sc),"source sc mismatch");}
   auto insert=[&](U key){require(owners.emplace(key,id).second,"ordinary duplicate or converse collision");require(factorial%c.aut==0,"invalid orbit divisor");mass+=factorial/c.aut;};
   insert(c.key);if(c.key!=d.key)insert(d.key);else ++sc_count;
   std::vector<unsigned> deg;for(int i=0;i<n;++i)deg.push_back(pc(r[i]));std::sort(deg.begin(),deg.end());std::ostringstream expected;for(int i=0;i<n;++i){if(i)expected<<',';expected<<deg[i];}if(!compact)require(scores==expected.str(),"source score mismatch");
   unsigned module=proper_module(n,r),p_mask=0;bool prime=!module;
   audit<<n<<'\t'<<id<<'\t'<<source<<'\t'<<c.key<<'\t'<<d.key<<'\t'<<c.aut<<'\t'<<sc<<'\t'<<module<<'\t'<<prime<<'\t';
   if(prime){
    ++primes;prime_ord+=sc?1:2;Paths p=literal_paths(n,r);path_nodes+=p.visited;U hs=0,he=0;
    for(int i=0;i<n;++i){hs+=p.starts[i];he+=p.ends[i];if((p.starts[i]+p.ends[i])&1u)p_mask|=1u<<i;}
    require(hs==he,"start/end sum mismatch");
    if(!compact){require(hs==H,"source H mismatch on prime");require((sc?hs:2*hs)/c.aut==tilings&&hs%c.aut==0,"source tilings mismatch on prime");}
    if(p_mask)++p_nonzero;else{++p_zero;zero_ord+=sc?1:2;if(n==9)candidates<<id<<'\t'<<source<<'\t'<<rows(r,n)<<'\t'<<list(p.starts,n)<<'\t'<<list(p.ends,n)<<'\n';}
    audit<<p_mask<<'\t'<<list(p.starts,n)<<'\t'<<list(p.ends,n)<<'\n';
   }else audit<<"-\t-\t-\n";
   ++merged;
  }
  U expected=U(1)<<(n*(n-1)/2);require(mass==expected,"labeled orbit mass not complete");
  summary<<n<<'\t'<<merged<<'\t'<<owners.size()<<'\t'<<sc_count<<'\t'<<mass<<'\t'<<expected<<'\t'<<perms<<'\t'<<primes<<'\t'<<p_nonzero<<'\t'<<p_zero<<'\t'<<prime_ord<<'\t'<<zero_ord<<'\t'<<path_nodes<<'\n';summary.flush();
  std::cerr<<"n="<<n<<" merged="<<merged<<" ordinary="<<owners.size()<<" orbit mass="<<mass<<" prime="<<primes<<" p-zero="<<p_zero<<" PASS\n";
 }
 std::cout<<"PASS: independent census coverage, automorphisms, merged partition, prime filter, literal endpoints "<<(compact?"at n=9":"through n=9")<<"\n";
 return 0;
 }catch(const std::exception&e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}}
