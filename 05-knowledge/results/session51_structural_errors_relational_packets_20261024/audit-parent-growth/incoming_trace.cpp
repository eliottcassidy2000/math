// Independent incoming-arc refinement. Numeric packed signatures, insertion
// grouping and reversed color order differ from the production array sorting.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;using Row=array<unsigned,10>;using Perm=array<unsigned,9>;
uint64_t checks=0,passes=0,early=0;void need(bool b,const char*s){++checks;if(!b)throw runtime_error(s);}unsigned bits(unsigned v){return __builtin_popcount(v);}
unsigned image(unsigned m,const Perm&p){unsigned result=0;while(m){unsigned v=__builtin_ctz(m);m&=m-1;result|=1u<<p[v];}return result;}
unsigned root(const Row&in){array<unsigned,10>c{},masks{},next{};unsigned count=1;
 for(;;){masks.fill(0);for(unsigned i=0;i<10;++i)masks[c[i]]|=1u<<i;
#ifndef FULL_REFINEMENT
 if(bits(masks[c[9]])==1){++early;return 1;}
#endif
 ++passes;array<uint64_t,10> sig{};for(unsigned v=0;v<10;++v){sig[v]=uint64_t(c[v])<<40;for(unsigned k=0;k<count;++k)sig[v]|=uint64_t(bits(in[v]&masks[k]))<<(4*k);}
 unsigned distinct=0;for(unsigned v=0;v<10;++v){unsigned w=0;for(;w<v;++w)if(sig[w]==sig[v])break;if(w<v)next[v]=next[w];else next[v]=distinct++;}
 if(distinct==count){unsigned size=0;for(unsigned v=0;v<10;++v)size+=next[v]==next[9];return size;}count=distinct;c=next;
 }
}
int main(int argc,char**argv){try{need(argc==5,"usage: incoming_trace atlas expected.tsv actual.tsv receipt.json");ifstream in(argv[1]),expected(argv[2]);ofstream detail(argv[3]);need(bool(in)&&bool(expected)&&bool(detail),"opened files");string header;getline(expected,header);need(header=="parent\tstates\tupper\tsize1\tsize3\tsize5\tsize7\tsize9","expected header");detail<<header<<"\n";unsigned n=0;in>>n;need(bool(in)&&n>0&&n<1000000,"atlas size");array<uint64_t,7>global{};uint64_t mass=0;
 for(unsigned i=0;i<n;++i){Row outgoing{},incoming{};for(unsigned v=0;v<9;++v)in>>outgoing[v];unsigned ng=0;in>>ng;need(bool(in)&&ng>0&&ng<100000,"action size");vector<Perm>g(ng);for(auto&p:g){for(auto&v:p){in>>v;need(bool(in)&&v<9,"valid action entry");}auto q=p;sort(q.begin(),q.end());for(unsigned v=0;v<9;++v)need(q[v]==v,"action permutation");}need(bool(in),"complete record");for(unsigned v=0;v<9;++v)for(unsigned w=0;w<9;++w)if(outgoing[w]>>v&1)incoming[v]|=1u<<w;need(362880%ng==0,"mass integral");mass+=362880/ng;
 array<uint64_t,7>subtotal{};for(unsigned mask=0;mask<512;++mask){unsigned least=mask;for(auto&p:g)least=min(least,image(mask,p));if(least!=mask)continue;Row child=incoming;child[9]=mask;for(unsigned v=0;v<9;++v)if(!(mask>>v&1))child[v]|=512;unsigned b=root(child);need(b>=1&&b<=9&&(b&1),"odd stable cell");++subtotal[0];subtotal[1]+=b;++subtotal[2+(b-1)/2];}
 unsigned ei=0;expected>>ei;need(bool(expected)&&ei==i,"expected parent index");detail<<i;for(unsigned k=0;k<7;++k){uint64_t x=0;expected>>x;need(bool(expected)&&x==subtotal[k],"per-parent incoming histogram matches");global[k]+=subtotal[k];detail<<"\t"<<subtotal[k];}detail<<"\n";need(bool(expected),"expected complete row");if((i+1)%20000==0)cerr<<"audited_parents="<<i+1<<"\n";}
 string extra;in>>extra;need(extra.empty(),"no atlas trailing data");expected>>extra;need(extra.empty(),"no expected trailing data");need(mass==(1ull<<36),"parent labelled mass");need(global[1]>5508&&(global[1]-5508)%10==0,"structurally corrected divisibility");ofstream out(argv[4]);out<<"{\n\"status\":\"PASS\",\n\"parent_classes\":"<<n<<",\n\"parent_mass\":"<<mass<<",\n\"pointed_states\":"<<global[0]<<",\n\"upper_trace\":"<<global[1]<<",\n\"structural_error_correction\":5508,\n\"prediction\":"<<(global[1]-5508)/10<<",\n\"stable_root_histogram\":[";for(unsigned k=0;k<5;++k){if(k)out<<",";out<<"["<<2*k+1<<","<<global[k+2]<<"]";}out<<"],\n\"per_parent_records_compared\":"<<n<<",\n\"passes\":"<<passes<<",\n\"singleton_early_exits\":"<<early<<",\n\"checks\":"<<checks<<"\n}\n";need(bool(out),"receipt written");cout<<"PASS incoming upper="<<global[1]<<" predicted="<<(global[1]-5508)/10<<"\n";
 }catch(exception&e){cerr<<e.what()<<"\n";return 1;}}
