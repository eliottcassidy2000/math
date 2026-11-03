#include <array>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
using I=long long;using V=std::array<I,3>;using M=std::array<I,9>;
struct Hash{size_t operator()(const V& v)const{return (v[0]*1000003+v[1])*1000033+v[2];}};
struct Node{V v;size_t parent;int letter;};
struct Atom{std::string name;M p;};
int main(int argc,char**argv){
 if(argc!=4)return 2;I bound=std::stoll(argv[1]);std::ifstream in(argv[2]);std::ofstream out(argv[3]);
 std::vector<Atom> atoms;Atom a;
 while(in>>a.name){for(I&x:a.p)in>>x;atoms.push_back(a);}
 std::sort(atoms.begin(),atoms.end(),[](const Atom&a,const Atom&b){return a.p<b.p||(a.p==b.p&&a.name<b.name);});
 atoms.erase(std::unique(atoms.begin(),atoms.end(),[](const Atom&a,const Atom&b){return a.p==b.p;}),atoms.end());
 std::vector<Node> nodes{{{1,0,0},0,-1}};std::unordered_map<V,size_t,Hash> seen;seen[{1,0,0}]=0;
 for(size_t z=0;z<nodes.size();++z){V v=nodes[z].v;
  for(int c=0;c<int(atoms.size());++c){const M&p=atoms[c].p;if(v[0]*p[0]>bound)break;
   I h=v[0]*p[0]+v[1]*p[3]+v[2]*p[6];if(h>bound)continue;
   V w={h,v[0]*p[1]+v[1]*p[4]+v[2]*p[7],v[0]*p[2]+v[1]*p[5]+v[2]*p[8]};
   if(seen.emplace(w,nodes.size()).second)nodes.push_back({w,z,c});
  }
 }
 for(size_t z=0;z<nodes.size();++z){const auto&n=nodes[z];
  out<<n.v[0]<<'\t'<<n.v[1]<<'\t'<<n.v[2]<<'\t'<<n.parent<<'\t';
  out<<(n.letter<0?"_":atoms[n.letter].name)<<'\n';
 }
 std::cout<<"matrices "<<atoms.size()<<" states "<<nodes.size()<<" bound "<<bound<<'\n';
}
