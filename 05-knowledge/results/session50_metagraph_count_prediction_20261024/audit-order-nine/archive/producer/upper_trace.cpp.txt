// Refinement-only target evaluator. No deletion, isomorphism, or child
// automorphism code. Parent actions are supplied in an authenticated atlas.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <vector>
using namespace std;
using Rows=array<unsigned,9>;
using Perm=array<unsigned,8>;
uint64_t checks=0,rounds=0,actions_used=0;
void need(bool c,const char*m){++checks;if(!c)throw runtime_error(m);}
int bits(unsigned x){return __builtin_popcount(x);}
int rootcell(const Rows &rows){
 array<int,9> colors{},next{},order{}; array<unsigned,9> masks{};
 for(int v=0;v<9;++v){colors[v]=bits(rows[v]);order[v]=v;}
 for(;;){
  ++rounds;masks.fill(0);int oldcount=0;
  for(int v=0;v<9;++v)masks[colors[v]]|=1u<<v;
  for(auto m:masks)oldcount+=(m!=0);
  array<array<int,10>,9> sig{};
  for(int v=0;v<9;++v){sig[v][0]=colors[v];for(int c=0;c<9;++c)sig[v][c+1]=bits(rows[v]&masks[c]);}
  sort(order.begin(),order.end(),[&](int a,int b){return sig[a]<sig[b];});
  int count=0;next[order[0]]=0;
  for(int j=1;j<9;++j){if(sig[order[j]]!=sig[order[j-1]])++count;next[order[j]]=count;}
  if(count+1==oldcount){int size=0;for(int v=0;v<9;++v)size+=(next[v]==next[8]);return size;}
  colors=next;
 }
}
int main(int argc,char**argv){
 try{
  need(argc==3,"usage: upper_trace parents-eight.txt output.json");ifstream input(argv[1]);need(bool(input),"input exists");
  int nparents;input>>nparents;need(nparents>0,"nonempty atlas");
  uint64_t states=0,upper=0;array<uint64_t,10> histogram{};vector<array<uint64_t,7>> perparent;
  for(int i=0;i<nparents;++i){
   Rows rows{};for(int v=0;v<8;++v)input>>rows[v];int na;input>>na;need(na>0,"parent nonempty action list");
   vector<Perm> group(na);for(auto &p:group)for(auto &v:p)input>>v;
   need(bool(input),"complete parent record");set<Perm> unique(group.begin(),group.end());need(unique.size()==group.size(),"distinct parent actions");
   Perm identity{};for(int v=0;v<8;++v)identity[v]=v;need(unique.count(identity),"identity supplied");
   for(int v=0;v<8;++v){need(rows[v]<256,"parent width");need(!(rows[v]>>v&1),"no loops");for(int w=0;w<v;++w)need(((rows[v]>>w&1)+(rows[w]>>v&1))==1,"parent tournament");}
   for(auto p:group){
    auto q=p;sort(q.begin(),q.end());need(q==identity,"action is permutation");
    for(int v=0;v<8;++v)for(int w=0;w<8;++w)need((rows[v]>>w&1)==(rows[p[v]]>>p[w]&1),"action preserves arcs");
    for(auto h:group){Perm composed{};for(int v=0;v<8;++v)composed[v]=p[h[v]];need(unique.count(composed),"actions closed");}
   }
   array<bool,256> seen{};array<uint64_t,7> summary{};
   for(unsigned mask=0;mask<256;++mask){
    if(seen[mask])continue;
    set<unsigned> orbit;
    for(auto p:group){unsigned image=0;for(int v=0;v<8;++v)if(mask>>v&1)image|=1u<<p[v];orbit.insert(image);++actions_used;}
    need(*orbit.begin()==mask,"minimum subset orbit representative");for(auto image:orbit){need(!seen[image],"disjoint subset orbit");seen[image]=true;}
    Rows child=rows;child[8]=255u^mask;for(int v=0;v<8;++v)if(mask>>v&1)child[v]|=1u<<8;
    int cell=rootcell(child);need(cell>=1&&cell<=9&&(cell&1),"odd stable root cell");
    ++states;upper+=cell;++histogram[cell];++summary[0];summary[1]+=cell;++summary[2+(cell-1)/2];
   }
   need(all_of(seen.begin(),seen.end(),[](bool x){return x;}),"subset orbits cover all attachments");perparent.push_back(summary);
  }
  string leftover;input>>leftover;need(leftover.empty(),"no trailing input");
  ofstream output(argv[2]);need(bool(output),"output exists");
  output<<"{\n  \"status\": \"PASS\",\n  \"target_order\": 9,\n  \"parent_classes\": "<<nparents<<",\n  \"pointed_extension_states\": "<<states<<",\n  \"upper_refinement_trace\": "<<upper<<",\n  \"stable_root_cell_histogram\": [";
  for(int i=1;i<=9;i+=2){if(i>1)output<<",";output<<"["<<i<<","<<histogram[i]<<"]";}
  output<<"],\n  \"refinement_rounds\": "<<rounds<<",\n  \"subset_action_evaluations\": "<<actions_used<<",\n  \"checks\": "<<checks<<",\n  \"child_canonicalizations\": 0,\n  \"child_deletions\": 0,\n  \"per_parent\": [\n";
  for(size_t i=0;i<perparent.size();++i){if(i)output<<",\n";output<<"[";for(int j=0;j<7;++j){if(j)output<<",";output<<perparent[i][j];}output<<"]";}
  output<<"\n  ]\n}\n";cout<<"states="<<states<<" upper="<<upper<<" rounds="<<rounds<<"\n";
 }catch(exception&e){cerr<<e.what()<<"\n";return 1;}
}
