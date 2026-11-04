#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <unordered_set>
#include <vector>
using Mask=std::uint32_t;using U=std::uint64_t;using Rows=std::array<unsigned,8>;
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
int pc(unsigned m){return __builtin_popcount(m);}
Rows read_graph(int n,Mask code){Rows a{};int k=0;for(int u=0;u<n;++u)for(int v=u+1;v<n;++v,++k)if(code&(1u<<k)){a[u]|=1u<<v;a[v]|=1u<<u;}return a;}
Mask code_at(int n,const Rows&a,const std::vector<int>&p){Mask c=0;int k=0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++k)if(a[p[i]]&(1u<<p[j]))c|=1u<<k;return c;}
bool possible(int n,const Rows&a){
 for(int u=0;u<n;++u)for(int v=u+1;v<n;++v)if(pc(a[u]&a[v])>2)return false;
 for(int center=0;center<n;++center){std::array<int,8> col;col.fill(-1);for(int root=0;root<n;++root)if((a[center]&(1u<<root))&&col[root]<0){std::vector<int> q{root};col[root]=0;for(size_t head=0;head<q.size();++head){int v=q[head];unsigned near=a[v]&a[center];while(near){int w=__builtin_ctz(near);near&=near-1;if(col[w]<0){col[w]=1-col[v];q.push_back(w);}else if(col[w]==col[v])return false;}}}}
 return true;
}
int main(int argc,char**argv){try{
 need(argc==2,"representative file");std::ifstream input(argv[1]);int count;input>>count;need(count==38,"38 input representatives");std::unordered_set<Mask> remaining;std::vector<int> p(8);std::iota(p.begin(),p.end(),0);U mass=0;
 for(int r=0;r<count;++r){Mask code;U size;input>>code>>size;need(bool(input),"input row");auto a=read_graph(8,code);need(possible(8,a),"representative passes independent filters");std::unordered_set<Mask> orbit;
  do{orbit.insert(code_at(8,a,p));}while(std::next_permutation(p.begin(),p.end()));need(orbit.size()==size,"independent orbit size");for(Mask m:orbit)need(remaining.insert(m).second,"input orbits disjoint");mass+=size;
 }
 U prefix_min2=0,prefix_filtered=0,children_tested=0,eligible=0;std::vector<int> id(8);std::iota(id.begin(),id.end(),0);
 for(Mask code=0;code<(1u<<21);++code){auto a=read_graph(7,code);bool okay=true;unsigned mandatory=0;for(int v=0;v<7;++v){int deg=pc(a[v]);okay &=deg>=2;if(deg==2)mandatory|=1u<<v;}if(!okay)continue;++prefix_min2;if(!possible(7,a))continue;++prefix_filtered;
  for(unsigned cut=0;cut<128;++cut){if(pc(cut)<3||(cut&mandatory)!=mandatory)continue;++children_tested;Rows b=a;b[7]=cut;for(int v=0;v<7;++v)if(cut&(1u<<v))b[v]|=128;if(!possible(8,b))continue;++eligible;need(remaining.erase(code_at(8,b,id))==1,"survivor exactly once in orbit union");}
 }
 need(remaining.empty()&&eligible==mass,"complete labelled survivor equality");std::cout<<"{\"universe\":268435456,\"designated_vertex7_deletion_parents\":2097152,\"parents_min_degree2\":"<<prefix_min2<<",\"parents_passing_hereditary_filters\":"<<prefix_filtered<<",\"admissible_children_checked\":"<<children_tested<<",\"eligible_graphs\":"<<eligible<<",\"representative_classes\":"<<count<<",\"orbit_union_equals_complete_survivor_set\":true,\"status\":\"FINITE-EXACT PASS\"}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
