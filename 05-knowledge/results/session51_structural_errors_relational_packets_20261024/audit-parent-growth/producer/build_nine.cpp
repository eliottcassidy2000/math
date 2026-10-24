// Canonical marked-root construction from the complete eight-vertex atlas.
// Only eight-vertex parents are canonicalized. No nine-vertex canonicalizer.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
template<int N> using Rows=array<unsigned,N>;
template<int N> using Perm=array<unsigned,N>;
using R8=Rows<8>;using R9=Rows<9>;using P8=Perm<8>;using P9=Perm<9>;
uint64_t checks=0,refinements=0,rounds=0,canon_calls=0,canon_candidates=0,reroots=0;
void need(bool b,const char*msg){++checks;if(!b)throw runtime_error(msg);}
int bits(unsigned x){return __builtin_popcount(x);}
template<size_t N> array<int,N> refine(const array<unsigned,N>&r){
 ++refinements;array<int,N> c{},nxt{},order{};array<unsigned,N> masks{};
 for(unsigned v=0;v<N;++v){c[v]=bits(r[v]);order[v]=v;}
 for(;;){
  ++rounds;masks.fill(0);for(unsigned v=0;v<N;++v)masks[c[v]]|=1u<<v;
  int oldcount=0;for(auto m:masks)oldcount+=(m!=0);
  array<array<int,N+1>,N> sig{};
  for(unsigned v=0;v<N;++v){sig[v][0]=c[v];for(unsigned k=0;k<N;++k)sig[v][k+1]=bits(r[v]&masks[k]);}
  sort(order.begin(),order.end(),[&](int a,int b){return sig[a]<sig[b];});
  int count=0;nxt[order[0]]=0;
  for(unsigned j=1;j<N;++j){if(sig[order[j]]!=sig[order[j-1]])++count;nxt[order[j]]=count;}
  if(count+1==oldcount) return nxt;
  c=nxt;
 }
}
template<size_t N> unsigned image(unsigned m,const array<unsigned,N>&p){unsigned out=0;for(unsigned v=0;v<N;++v)if(m>>v&1)out|=1u<<p[v];return out;}
template<size_t N> bool preserves(const array<unsigned,N>&r,const array<unsigned,N>&p){for(unsigned v=0;v<N;++v)if(image(r[v],p)!=r[p[v]])return false;return true;}
template<size_t N> void validate_rows(const array<unsigned,N>&r){
 for(unsigned v=0;v<N;++v){need(r[v]<(1u<<N),"row width");need(!(r[v]>>v&1),"no loops");for(unsigned w=0;w<v;++w)need(((r[v]>>w&1)+(r[w]>>v&1))==1,"tournament");}
}
template<size_t N> void validate_actions(const array<unsigned,N>&r,const vector<array<unsigned,N>>&g){
 array<unsigned,N> id{};for(unsigned v=0;v<N;++v)id[v]=v;
 set<array<unsigned,N>> distinct(g.begin(),g.end());need(!g.empty()&&distinct.size()==g.size(),"nonempty distinct actions");need(distinct.count(id),"identity action");
 for(auto p:g){auto q=p;sort(q.begin(),q.end());need(q==id,"action permutation");need(preserves(r,p),"action preserves arcs");for(auto h:g){array<unsigned,N> c{};for(unsigned v=0;v<N;++v)c[v]=p[h[v]];need(distinct.count(c),"action group closed");}}
}
struct Parent{R8 rows;vector<P8> actions;};
struct Child{R9 rows;vector<P9> actions;bool residual;};
struct Canon{R8 rows{};P8 old_to_new{};bool initialized=false;};
void canonical_visit(const R8&r,const vector<vector<unsigned>>&cells,unsigned k,vector<unsigned>&order,Canon&out){
 if(k==cells.size()){
  ++canon_candidates;P8 p{};for(unsigned v=0;v<8;++v)p[order[v]]=v;R8 val{};for(unsigned v=0;v<8;++v)val[p[v]]=image(r[v],p);
  if(!out.initialized||val<out.rows){out.rows=val;out.old_to_new=p;out.initialized=true;}return;
 }
 auto part=cells[k];do{for(auto v:part)order.push_back(v);canonical_visit(r,cells,k+1,order,out);order.resize(order.size()-part.size());}while(next_permutation(part.begin(),part.end()));
}
Canon canonical_eight(const R8&r){
 // The signature enforces the order-eight-only source boundary.
 ++canon_calls;auto colors=refine(r);vector<vector<unsigned>> cells(*max_element(colors.begin(),colors.end())+1);
 for(unsigned v=0;v<8;++v)cells[colors[v]].push_back(v);
 Canon out;vector<unsigned> order;canonical_visit(r,cells,0,order,out);need(out.initialized,"eight canonical result");return out;
}
struct Pointed{pair<int,unsigned> key;unsigned mapped_subset;array<unsigned,8> inverse_child_labels;};
Pointed pointed_key(const R9&r,unsigned root,const map<R8,int>&index,const vector<Parent>&parents){
 ++reroots;R8 smaller{};array<unsigned,8> labels{};unsigned k=0;for(unsigned v=0;v<9;++v)if(v!=root)labels[k++]=v;
 unsigned incoming=0;for(unsigned v=0;v<8;++v){for(unsigned w=0;w<8;++w)if(r[labels[v]]>>labels[w]&1)smaller[v]|=1u<<w;if(r[labels[v]]>>root&1)incoming|=1u<<v;}
 Canon can=canonical_eight(smaller);auto it=index.find(can.rows);need(it!=index.end(),"deleted child canonical parent in complete order-eight atlas");
 Pointed out;out.key.first=it->second;out.mapped_subset=image(incoming,can.old_to_new);out.key.second=out.mapped_subset;
 for(const auto&p:parents[it->second].actions)out.key.second=min(out.key.second,image(out.mapped_subset,p));
 for(unsigned v=0;v<8;++v)out.inverse_child_labels[can.old_to_new[v]]=labels[v];
 return out;
}
unsigned error(const R9&r,const vector<P9>&g){auto c=refine(r);array<unsigned,9>size{};for(auto x:c)++size[x];array<bool,9>seen{};unsigned e=0;for(unsigned v=0;v<9;++v)if(!seen[v]){unsigned orbit=0;for(auto p:g)orbit|=1u<<p[v];need(orbit>>v&1,"orbit contains vertex");for(unsigned w=0;w<9;++w)if(orbit>>w&1){need(c[w]==c[v],"orbit preserves stable cell");seen[w]=true;}e+=size[c[v]]-bits(orbit);}return e;}
int main(int argc,char**argv){try{
 need(argc==4,"usage: build_nine parents-eight.txt parents-nine.txt receipt.json");ifstream input(argv[1]);need(bool(input),"input opened");int count=0;input>>count;need(bool(input)&&count>0,"valid nonempty input header");vector<Parent>parents(count);map<R8,int>index;uint64_t input_mass=0;
 for(int i=0;i<count;++i){auto&p=parents[i];for(auto&r:p.rows)input>>r;int ng=0;input>>ng;need(bool(input)&&ng>0&&ng<=40320,"valid positive parent action count");p.actions.resize(ng);for(auto&g:p.actions)for(auto&v:g)input>>v;need(bool(input),"complete input record");validate_rows(p.rows);validate_actions(p.rows,p.actions);need(40320u%p.actions.size()==0,"parent action order divides eight factorial");input_mass+=40320u/p.actions.size();need(index.emplace(p.rows,i).second,"distinct literal input parents");}
 string trailing;input>>trailing;need(trailing.empty(),"no trailing input");need(input_mass==(1ull<<28),"input labelled tournament mass");
 vector<Child>children;uint64_t states=0,with_singleton=0,without_singleton=0,singleton_children=0,residual_children=0,mass=0,total_error=0,upper9=0;
 map<unsigned,uint64_t>aut_hist,error_hist;map<unsigned,uint64_t>root_hist;
 for(int parent_id=0;parent_id<count;++parent_id){const auto&p=parents[parent_id];array<bool,256>seen{};
  for(unsigned mask=0;mask<256;++mask){if(seen[mask])continue;set<unsigned>orbit;for(auto g:p.actions)orbit.insert(image(mask,g));need(*orbit.begin()==mask,"minimum attachment representative");for(auto m:orbit){need(!seen[m],"disjoint subset orbit");seen[m]=true;}
   ++states;R9 r{};for(unsigned v=0;v<8;++v)r[v]=p.rows[v]|((mask>>v&1)<<8);r[8]=255u^mask;auto colors=refine(r);array<unsigned,9>sizes{};for(auto c:colors)++sizes[c];upper9+=sizes[colors[8]];++root_hist[sizes[colors[8]]];
   int selected=-1;for(unsigned v=0;v<9;++v)if(sizes[colors[v]]==1&&(selected<0||colors[v]<colors[selected]))selected=v;
   vector<P9>group;bool residual=selected<0;
   if(!residual){++with_singleton;if(selected!=8)continue;++singleton_children;
    for(auto h:p.actions)if(image(mask,h)==mask){P9 g{};for(unsigned v=0;v<8;++v)g[v]=h[v];g[8]=8;group.push_back(g);}
   }else{++without_singleton;vector<unsigned>partition;for(auto size:sizes)if(size)partition.push_back(size);sort(partition.begin(),partition.end());need(partition==vector<unsigned>{9}||partition==vector<unsigned>({3,3,3}),"singleton-free order-nine partition guard");
    array<Pointed,9>points;pair<int,unsigned>best={count,256};for(unsigned v=0;v<9;++v){points[v]=pointed_key(r,v,index,parents);best=min(best,points[v].key);}pair<int,unsigned>own={parent_id,mask};need(points[8].key==own,"own pointed key recovered");if(own!=best)continue;++residual_children;
    for(unsigned root=0;root<9;++root)if(points[root].key==own){for(auto h:p.actions)if(image(mask,h)==points[root].mapped_subset){P9 g{};for(unsigned v=0;v<8;++v)g[v]=points[root].inverse_child_labels[h[v]];g[8]=root;group.push_back(g);}}
   }
   validate_rows(r);validate_actions(r,group);need(362880u%group.size()==0,"child action order divides nine factorial");mass+=362880u/group.size();unsigned e=error(r,group);total_error+=e;++aut_hist[group.size()];++error_hist[e];children.push_back({r,group,residual});
  }
  need(all_of(seen.begin(),seen.end(),[](bool b){return b;}),"all parent attachments covered");
 }
 need(mass==(1ull<<36),"complete order-nine labelled orbit mass");need(total_error==758,"independent structural order-nine error matches constructed atlas");need(9*children.size()+total_error==upper9,"root trace agrees with selected class count");
 ofstream atlas(argv[2]);need(bool(atlas),"output atlas opened");atlas<<children.size()<<"\n";for(const auto&child:children){for(auto r:child.rows)atlas<<r<<" ";atlas<<child.actions.size()<<"\n";for(auto g:child.actions){for(unsigned v=0;v<9;++v){if(v)atlas<<" ";atlas<<g[v];}atlas<<"\n";}}atlas.close();need(bool(atlas),"atlas output complete");
 ofstream out(argv[3]);need(bool(out),"receipt opened");out<<"{\n\"status\":\"PASS\",\n\"parent_order\":8,\n\"constructed_order\":9,\n\"parent_classes\":"<<count<<",\n\"constructed_classes\":"<<children.size()<<",\n\"pointed_states\":"<<states<<",\n\"states_with_singleton\":"<<with_singleton<<",\n\"singleton_free_states\":"<<without_singleton<<",\n\"selected_singleton_classes\":"<<singleton_children<<",\n\"selected_singleton_free_classes\":"<<residual_children<<",\n\"labelled_orbit_mass\":"<<mass<<",\n\"refinement_error\":"<<total_error<<",\n\"upper_trace_order_nine\":"<<upper9<<",\n\"maximum_canonicalized_order\":8,\n\"order_nine_canonicalizations\":0,\n\"parent_canonicalizations\":"<<canon_calls<<",\n\"parent_labelling_candidates\":"<<canon_candidates<<",\n\"child_reroots\":"<<reroots<<",\n\"refinement_calls\":"<<refinements<<",\n\"refinement_rounds\":"<<rounds<<",\n\"checks\":"<<checks<<",\n\"automorphism_order_histogram\":[";
 bool first=true;for(auto x:aut_hist){if(!first)out<<",";first=false;out<<"["<<x.first<<","<<x.second<<"]";}out<<"],\n\"error_histogram\":[";first=true;for(auto x:error_hist){if(!first)out<<",";first=false;out<<"["<<x.first<<","<<x.second<<"]";}out<<"],\n\"stable_root_cell_histogram\":[";first=true;for(auto x:root_hist){if(!first)out<<",";first=false;out<<"["<<x.first<<","<<x.second<<"]";}out<<"]\n}\n";
 cout<<"classes="<<children.size()<<" singleton="<<singleton_children<<" residual="<<residual_children<<" residual_states="<<without_singleton<<" mass="<<mass<<" parent_canon="<<canon_calls<<"\n";
 }catch(exception&e){cerr<<e.what()<<"\n";return 1;}}
