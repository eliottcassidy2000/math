// Check a proposed reachable set by parent certificates and transition closure.
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>
using I=long long;using V=std::array<I,3>;using M=std::array<I,9>;
struct Hash{size_t operator()(const V&v)const{return v[0]+1000003*v[1]+1000033*v[2];}};
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
V step(V v,M p){V w{};for(int j=0;j<3;++j)for(int k=0;k<3;++k)w[j]+=v[k]*p[3*k+j];return w;}
int main(int argc,char**argv){
 need(argc==4,"matrices states bound");I bound=std::stoll(argv[3]);std::ifstream a(argv[1]),b(argv[2]);
 std::map<std::string,M> bank;std::map<M,bool> unique;std::string name;M p;
 while(a>>name){for(I&x:p)a>>x;bank[name]=p;unique[p]=true;}
 std::vector<V> list;std::unordered_set<V,Hash> states;V v;size_t parent;
 while(b>>v[0]>>v[1]>>v[2]>>parent>>name){
  if(list.empty())need(v==V{1,0,0}&&name=="_","root");
  else{need(parent<list.size()&&bank.count(name),"parent address");need(step(list[parent],bank.at(name))==v,"literal parent multiplication");}
  need(v[0]<=bound&&states.insert(v).second,"bound and uniqueness");list.push_back(v);
 }
 size_t transitions=0;
 for(const V&v:list)for(const auto&entry:unique){const M&p=entry.first;if(v[0]*p[0]>bound)break;
  if(v[0]*p[0]+v[1]*p[3]+v[2]*p[6]>bound)continue;
  need(states.count(step(v,p)),"transition closure");++transitions;
 }
 std::cout<<"closed states "<<states.size()<<" transitions "<<transitions<<" bound "<<bound<<'\n';
}
