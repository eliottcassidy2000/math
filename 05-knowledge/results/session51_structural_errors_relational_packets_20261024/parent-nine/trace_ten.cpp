// Refinement-only extension trace from authenticated nine-vertex parents.
// No deletion, canonicalization, or target automorphism computation.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
using Rows=array<unsigned,10>;using Perm=array<unsigned,9>;
uint64_t checks=0,rounds=0,actions_used=0,early_singletons=0;
void need(bool c,const char*m){++checks;if(!c)throw runtime_error(m);}
int bits(unsigned x){return __builtin_popcount(x);}
unsigned image(unsigned mask,const Perm&p){unsigned result=0;for(unsigned v=0;v<9;++v)if(mask>>v&1)result|=1u<<p[v];return result;}
int rootcell(const Rows&r){
 array<int,10> colors{},next{},order{};array<unsigned,10> masks{};
 for(int v=0;v<10;++v){colors[v]=bits(r[v]);order[v]=v;}
 for(;;){
  masks.fill(0);for(int v=0;v<10;++v)masks[colors[v]]|=1u<<v;
  #ifndef FULL_REFINEMENT
  if(bits(masks[colors[9]])==1){++early_singletons;return 1;}
  #endif
  ++rounds;int oldcount=0;for(auto m:masks)oldcount+=(m!=0);
  array<array<int,11>,10> sig{};
  for(int v=0;v<10;++v){sig[v][0]=colors[v];for(int c=0;c<10;++c)sig[v][c+1]=bits(r[v]&masks[c]);}
  sort(order.begin(),order.end(),[&](int a,int b){return sig[a]<sig[b];});
  int count=0;next[order[0]]=0;
  for(int j=1;j<10;++j){if(sig[order[j]]!=sig[order[j-1]])++count;next[order[j]]=count;}
  if(count+1==oldcount){int size=0;for(int v=0;v<10;++v)size+=(next[v]==next[9]);return size;}
  colors=next;
 }
}
int main(int argc,char**argv){try{
 need(argc==4,"usage: trace_ten parents-nine.txt receipt.json per-parent.tsv");ifstream input(argv[1]);need(bool(input),"input exists");int nparents=0;input>>nparents;need(bool(input)&&nparents>0,"valid nonempty atlas header");
 uint64_t states=0,upper=0,parent_mass=0;array<uint64_t,11>histogram{};
 ofstream detail(argv[3]);need(bool(detail),"per-parent output opened");detail<<"parent\tstates\tupper\tsize1\tsize3\tsize5\tsize7\tsize9\n";
 for(int i=0;i<nparents;++i){
  Rows rows{};for(int v=0;v<9;++v)input>>rows[v];int na=0;input>>na;need(bool(input)&&na>0&&na<=362880,"valid nonempty parent action list");vector<Perm>group(na);for(auto&p:group)for(auto&v:p)input>>v;need(bool(input),"complete parent record");set<Perm>unique(group.begin(),group.end());need(unique.size()==group.size(),"distinct parent actions");
  Perm identity{};for(int v=0;v<9;++v)identity[v]=v;need(unique.count(identity),"identity supplied");need(362880u%group.size()==0,"parent automorphism order divides nine factorial");parent_mass+=362880u/group.size();
  for(int v=0;v<9;++v){need(rows[v]<512,"parent width");need(!(rows[v]>>v&1),"no loops");for(int w=0;w<v;++w)need(((rows[v]>>w&1)+(rows[w]>>v&1))==1,"parent tournament");}
  for(auto p:group){auto q=p;sort(q.begin(),q.end());need(q==identity,"action is permutation");for(int v=0;v<9;++v)need(image(rows[v],p)==rows[p[v]],"action preserves arcs");for(auto h:group){Perm composed{};for(int v=0;v<9;++v)composed[v]=p[h[v]];need(unique.count(composed),"actions closed");}}
  array<bool,512>seen{};array<uint64_t,7>summary{};
  for(unsigned mask=0;mask<512;++mask){if(seen[mask])continue;
   if(na==1){seen[mask]=true;++actions_used;}else{set<unsigned>orbit;for(auto p:group){orbit.insert(image(mask,p));++actions_used;}need(*orbit.begin()==mask,"minimum subset orbit representative");for(auto m:orbit){need(!seen[m],"disjoint subset orbit");seen[m]=true;}}
   Rows child=rows;child[9]=511u^mask;for(int v=0;v<9;++v)if(mask>>v&1)child[v]|=1u<<9;
   int cell=rootcell(child);need(cell>=1&&cell<=9&&(cell&1),"odd stable root cell");++states;upper+=cell;++histogram[cell];++summary[0];summary[1]+=cell;++summary[2+(cell-1)/2];
  }
  need(all_of(seen.begin(),seen.end(),[](bool x){return x;}),"subset orbit coverage");detail<<i;for(auto x:summary)detail<<"\t"<<x;detail<<"\n";
  if((i+1)%20000==0)cerr<<"completed_parents="<<i+1<<" states="<<states<<"\n";
 }
 string trailing;input>>trailing;need(trailing.empty(),"no trailing input");need(parent_mass==(1ull<<36),"authenticated parent labelled orbit mass");
 const uint64_t correction=5508; // Proven structural order-ten error, not a target count.
 need(upper>correction&&(upper-correction)%10==0,"corrected trace positive multiple of ten");uint64_t predicted=(upper-correction)/10;
 ofstream out(argv[2]);need(bool(out),"receipt opened");out<<"{\n\"status\":\"PASS\",\n\"target_order\":10,\n\"parent_classes\":"<<nparents<<",\n\"parent_labelled_mass\":"<<parent_mass<<",\n\"pointed_extension_states\":"<<states<<",\n\"upper_refinement_trace\":"<<upper<<",\n\"structural_error_correction\":"<<correction<<",\n\"predicted_tournament_count\":"<<predicted<<",\n\"stable_root_cell_histogram\":[";
 for(int i=1;i<=9;i+=2){if(i>1)out<<",";out<<"["<<i<<","<<histogram[i]<<"]";}out<<"],\n\"refinement_rounds\":"<<rounds<<",\n\"singleton_early_exits\":"<<early_singletons<<",\n\"subset_action_evaluations\":"<<actions_used<<",\n\"checks\":"<<checks<<",\n\"target_canonicalizations\":0,\n\"target_deletions\":0,\n\"target_automorphism_computations\":0\n}\n";
 cout<<"states="<<states<<" upper="<<upper<<" predicted="<<predicted<<"\n";
 }catch(exception&e){cerr<<e.what()<<"\n";return 1;}}
