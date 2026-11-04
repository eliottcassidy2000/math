#include <array>
#include <vector>
#include <iostream>
#include <cstdint>
#include <functional>
#include <string>
using V=std::array<uint64_t,36>;
int index(int o,int a,int b,int c){return o*18+a*6+b*2+c;}
std::array<std::array<std::vector<int>,36>,4> tr;
V step(const V&v,int letter){V w{};for(int i=0;i<36;i++)if(v[i])for(int j:tr[letter][i])w[j]+=v[i];return w;}
uint64_t value(const V&v,int edge){uint64_t h=0;for(int o=0;o<2;o++)for(int a=0;a<2;a++)for(int b=0;b<2;b++)for(int c=0;c<2;c++)if(b || ((o==0)==bool(edge)))h+=v[index(o,a,b,c)];return h;}
int main(int argc,char**argv){int bound=argc>1?std::stoi(argv[1]):13,limit=argc>2?std::stoi(argv[2]):100;
 for(int l=0;l<4;l++)for(int o=0;o<2;o++)for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int c=0;c<2;c++){
  int u=(l>>o)&1,w=(l>>(1-o))&1,i=index(o,a,b,c);
  tr[l][i].push_back(index(o,u?2:1,b,c));
  if(b || u)tr[l][i].push_back(index(o,a,w?2:1,c));
  if(c || w)tr[l][i].push_back(index(o,a,b,1));
 }
 V init{};init[0]=init[18]=1;
 if(bound<0){int edge;std::string letters;while(std::cin>>edge>>letters){auto v=init;if(letters!="-")for(char c:letters)v=step(v,c-'0');std::cout<<value(v,edge)<<'\n';}return 0;}
 auto v1=step(init,1),v2=step(v1,1);std::string word;uint64_t nodes=0,pruned=0;int hits=0;
 std::function<void(const V&,const V&)> dfs=[&](const V&a,const V&b){
  if(hits>=limit)return;nodes++;bool good=false;
  for(int edge=0;edge<2;edge++){
   auto h=value(a,edge),k=value(b,edge);if(h<=4277&&k<=11735)good=true;
   if(h==4277&&k==11735){std::cout<<edge<<" "<<word<<"\n";std::cout.flush();hits++;}
  }
  if(!good){pruned++;return;}if((int)word.size()>=bound)return;
  for(int l=0;l<4;l++){word.push_back('0'+l);dfs(step(a,l),step(b,l));word.pop_back();}
 };
 dfs(v1,v2);std::cerr<<"nodes "<<nodes<<" pruned "<<pruned<<" hits "<<hits<<"\n";
}
