// Independent audit only. Incoming-edge refinement and full color-preserving
// canonicalization at orders eight/nine; production does NOT use this path.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
using Vec=vector<unsigned>; using Graph=Vec;
uint64_t checks=0,candidates=0;
void require(bool b,const string&s){++checks;if(!b)throw runtime_error(s);}
unsigned pc(unsigned x){return __builtin_popcount(x);}
unsigned act(unsigned mask,const Vec&p){unsigned y=0;while(mask){unsigned v=__builtin_ctz(mask);mask&=mask-1;y|=1u<<p[v];}return y;}
vector<Vec> refine(const Graph&r){
 unsigned n=r.size();Graph incoming(n);for(unsigned i=0;i<n;++i)for(unsigned j=0;j<n;++j)if(r[j]>>i&1)incoming[i]|=1u<<j;
 vector<Vec> cells(1);for(unsigned v=0;v<n;++v)cells[0].push_back(v);
 for(;;){vector<unsigned> masks;for(auto&c:cells){unsigned m=0;for(auto v:c)m|=1u<<v;masks.push_back(m);}vector<Vec>next;
  for(auto&c:cells){map<vector<unsigned>,Vec>split;for(auto v:c){vector<unsigned>sig;for(auto m:masks)sig.push_back(pc(incoming[v]&m));split[sig].push_back(v);}for(auto&kv:split)next.push_back(kv.second);}
  if(next.size()==cells.size())return cells;cells=move(next);
 }
}
struct Canon{uint64_t key=UINT64_MAX;Vec first;vector<Vec>autos;vector<Vec>cells;};
void visit(const Graph&r,const vector<Vec>&cells,unsigned level,Vec&order,Canon&result){
 if(level==cells.size()){
  ++candidates;uint64_t key=0;unsigned bit=0;for(unsigned i=0;i<r.size();++i)for(unsigned j=i+1;j<r.size();++j)key|=uint64_t((r[order[i]]>>order[j])&1)<<bit++;
  if(key<result.key){result.key=key;result.first=order;result.autos.clear();}
  if(key==result.key){Vec g(r.size());for(unsigned i=0;i<r.size();++i)g[order[i]]=result.first[i];result.autos.push_back(g);}return;
 }
 Vec block=cells[level];sort(block.begin(),block.end());do{for(auto v:block)order.push_back(v);visit(r,cells,level+1,order,result);order.resize(order.size()-block.size());}while(next_permutation(block.begin(),block.end()));
}
Canon canon(const Graph&r){Canon c;c.cells=refine(r);Vec order;visit(r,c.cells,0,order,c);sort(c.autos.begin(),c.autos.end());return c;}
struct Item{Graph rows;vector<Vec>autos;};
vector<Item> read(const string&path,unsigned n){ifstream f(path);require(bool(f),"input opened "+path);unsigned cnt=0;f>>cnt;require(bool(f)&&cnt>0&&cnt<1000000,"input count");vector<Item>atlas(cnt);for(auto&it:atlas){it.rows.resize(n);for(auto&v:it.rows)f>>v;unsigned a=0;f>>a;require(bool(f)&&a>0&&a<100000,"automorphism count");it.autos.resize(a,Vec(n));for(auto&p:it.autos)for(auto&v:p)f>>v;require(bool(f),"record complete");for(unsigned v=0;v<n;++v){require(it.rows[v]<(1u<<n)&&!(it.rows[v]>>v&1),"width loop");for(unsigned w=0;w<v;++w)require(((it.rows[v]>>w&1)+(it.rows[w]>>v&1))==1,"arcs");}sort(it.autos.begin(),it.autos.end());}string tail;f>>tail;require(tail.empty(),"trailing content");return atlas;}
unsigned defect(const Canon&c){unsigned e=0;set<unsigned>done;for(auto&cell:c.cells)for(auto v:cell)if(!done.count(v)){set<unsigned>orb;for(auto&p:c.autos)orb.insert(p[v]);for(auto w:orb)done.insert(w);e+=cell.size()-orb.size();}return e;}
bool strong(const Graph&r){for(unsigned v=0;v<r.size();++v){unsigned seen=1u<<v;for(;;){unsigned next=seen;for(unsigned w=0;w<r.size();++w)if(seen>>w&1)next|=r[w];if(next==seen)break;seen=next;}if(seen!=(1u<<r.size())-1)return false;}return true;}
void cell_defects(const Canon&c,map<unsigned,uint64_t>&out){for(auto&cell:c.cells){set<unsigned>seen;for(auto v:cell)if(!seen.count(v)){set<unsigned>orbit;for(auto&p:c.autos)orbit.insert(p[v]);seen.insert(orbit.begin(),orbit.end());unsigned e=cell.size()-orbit.size();if(e)out[cell.size()]+=e;}}}
void hist(ostream&o,const map<unsigned,uint64_t>&m){o<<"[";bool first=true;for(auto [x,y]:m){if(!first)o<<",";first=false;o<<"["<<x<<","<<y<<"]";}o<<"]";}
int main(int argc,char**argv){try{
 require(argc==4,"usage: independent parents-eight parents-nine receipt");auto p8=read(argv[1],8),p9=read(argv[2],9);set<uint64_t>keys8,keys9;map<uint64_t,pair<unsigned,unsigned>>expected;uint64_t mass8=0,mass9=0,e9=0;map<unsigned,uint64_t>autHist,errorHist,totalCellError,strongCellError;
 for(auto&it:p8){auto c=canon(it.rows);require(c.autos==it.autos,"full independent order-eight action equality");require(keys8.insert(c.key).second,"nonisomorphic eight parents");require(40320%c.autos.size()==0,"eight orbit integrality");mass8+=40320/c.autos.size();}require(mass8==(1ull<<28),"complete eight labelled mass");
 for(auto&it:p9){auto c=canon(it.rows);require(c.autos==it.autos,"full independent order-nine action equality");require(keys9.insert(c.key).second,"nonisomorphic nine parents");require(362880%c.autos.size()==0,"nine orbit integrality");mass9+=362880/c.autos.size();auto e=defect(c);e9+=e;cell_defects(c,totalCellError);if(strong(it.rows))cell_defects(c,strongCellError);++autHist[c.autos.size()];++errorHist[e];expected[c.key]={c.autos.size(),e};}require(mass9==(1ull<<36),"complete nine labelled mass");
 // Independent constructive path: first incoming-refinement singleton, and
 // canonical set deduplication for the small singleton-free residue.
 set<uint64_t>rebuilt;uint64_t states=0,single=0,residualStates=0,residualClasses=0,upper=0;
 map<unsigned,uint64_t>rootHist;
 for(auto&it:p8){set<unsigned>attachments;for(unsigned mask=0;mask<256;++mask){unsigned minimum=mask;for(auto&p:it.autos)minimum=min(minimum,act(mask,p));if(mask!=minimum)continue;require(attachments.insert(mask).second,"attachment disjoint");++states;Graph r=it.rows;r.push_back(255^mask);for(unsigned v=0;v<8;++v)if(mask>>v&1)r[v]|=256;auto cells=refine(r);unsigned chosen=9,rootSize=0;for(auto&cell:cells){if(cell.size()==1&&chosen==9)chosen=cell[0];if(find(cell.begin(),cell.end(),8)!=cell.end())rootSize=cell.size();}require(rootSize>0,"root present");upper+=rootSize;++rootHist[rootSize];if(chosen!=9&&chosen!=8)continue;auto c=canon(r);auto pos=expected.find(c.key);require(pos!=expected.end(),"rebuilt key present in production atlas");require(pos->second==make_pair(unsigned(c.autos.size()),defect(c)),"rebuilt actions and defect");if(chosen==9){++residualStates;if(rebuilt.insert(c.key).second)++residualClasses;}else{++single;require(rebuilt.insert(c.key).second,"singleton construction distinct");}}
 }
 require(rebuilt==keys9,"independent reconstructed class set equals production");require(9*keys9.size()+e9==upper,"independent trace identity");
 ofstream out(argv[3]);out<<"{\n\"status\":\"PASS\",\n\"parent_classes\":"<<keys8.size()<<",\n\"constructed_classes\":"<<keys9.size()<<",\n\"parent_mass\":"<<mass8<<",\n\"constructed_mass\":"<<mass9<<",\n\"pointed_states\":"<<states<<",\n\"independent_singleton_classes\":"<<single<<",\n\"residual_states\":"<<residualStates<<",\n\"residual_classes\":"<<residualClasses<<",\n\"refinement_error\":"<<e9<<",\n\"upper_trace_order_nine\":"<<upper<<",\n\"automorphism_order_histogram\":";hist(out,autHist);out<<",\n\"error_histogram\":";hist(out,errorHist);out<<",\n\"root_cell_histogram\":";hist(out,rootHist);out<<",\n\"order_nine_error_by_stable_cell_size\":";hist(out,totalCellError);out<<",\n\"strong_order_nine_error_by_stable_cell_size\":";hist(out,strongCellError);out<<",\n\"audit_only_order_nine_canonicalization\":true,\n\"canonical_labelling_candidates\":"<<candidates<<",\n\"checks\":"<<checks<<"\n}\n";require(bool(out),"receipt complete");cout<<"PASS classes="<<keys9.size()<<" mass="<<mass9<<" error="<<e9<<" residual="<<residualClasses<<" candidates="<<candidates<<"\n";
 }catch(exception&e){cerr<<e.what()<<"\n";return 1;}}
