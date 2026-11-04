// Two feedback vertices over a transitive base. Each of three ordered gaps
// is empty, has a bad last boundary, or a good last boundary (last gap only 0/1).
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;using V=array<uint32_t,108>;
array<array<vector<int>,18>,4> maps;
array<array<vector<int>,2>,2> readout;
void init(){
 for(int c=0;c<4;c++)for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int d=0;d<2;d++){
  int id=(a*3+b)*2+d,w=c&1,u=c>>1;maps[c][id].push_back(((w?1:2)*3+b)*2+d);
  if(b||w)maps[c][id].push_back((a*3+(u?1:2))*2+d);
  if(d||u)maps[c][id].push_back((a*3+b)*2+1);
 }
 for(int order=0;order<2;order++)for(int edge=0;edge<2;edge++)for(int a:{0,2})for(int b:{0,2})for(int d=0;d<2;d++)
  if(b||edge==order)readout[order][edge].push_back((a*3+b)*2+d);
}
V step(const V&v,int c){V o{};for(int m=0;m<3;m++)for(int a=0;a<2;a++)for(int i=0;i<18;i++)if(auto x=v[m*36+a*18+i])
 for(int j:maps[a?((c&1)*2+(c>>1)):c][i])o[m*36+a*18+j]+=x;return o;}
array<uint64_t,3> H(const V&v,int edge){array<uint64_t,3>h{};for(int m=0;m<3;m++)for(int a=0;a<2;a++)for(int i:readout[a][edge])h[m]+=v[m*36+a*18+i];return h;}
uint64_t visits=0,pruned=0,hits=0,prefixes=0;int maxdepth=13,bulk,position=0;vector<int>word;
const array<array<uint64_t,3>,2>target={array<uint64_t,3>{4277,11735,26651},array<uint64_t,3>{4279,11737,26653}};
void dfs(const V&v,int depth,int lanes){visits++;int next=0;
 for(int edge=0;edge<2;edge++)if(lanes&(3<<(2*edge))){auto h=H(v,edge);
  for(int b=0;b<2;b++)if(lanes&(1<<(2*edge+b))){bool ok=true;for(int m=0;m<3;m++)if(h[m]>target[b][m])ok=false;if(ok){next|=1<<(2*edge+b);if(h==target[b]){hits++;cout<<"{\"offset\":"<<2+2*b<<",\"bulk\":"<<bulk<<",\"position\":"<<position<<",\"edge\":"<<edge<<",\"word\":[";for(int j=0;j<depth;j++){if(j)cout<<',';cout<<word[j];}cout<<"]}\n";cout.flush();}}}
 }
 if(!next){pruned++;return;}if(depth==maxdepth)return;
 for(int c=0;c<4;c++){word.push_back(c);dfs(step(v,c),depth+1,next);word.pop_back();}
}
void before(const V&v,int depth,int lanes){prefixes++;int next=0;for(int e=0;e<2;e++){auto h=H(v,e);for(int b=0;b<2;b++)if((lanes&(1<<(2*e+b)))&&h[0]<=target[b][0])next|=1<<(2*e+b);}if(!next)return;
 for(bulk=1;bulk<=3;bulk+=2){V u=v,z{};for(int m=0;m<3;m++){u=step(u,bulk);for(int i=0;i<36;i++)z[m*36+i]=u[m*36+i];}position=depth;dfs(z,depth,next);}
 if(depth==maxdepth)return;for(int c=0;c<4;c++){word.push_back(c);before(step(v,c),depth+1,next);word.pop_back();}
}
int main(int argc,char**argv){init();if(argc>1)maxdepth=stoi(argv[1]);
 if(argc>2){V v{};for(int m=0;m<3;m++)v[m*36]=v[m*36+18]=1;before(v,0,15);cerr<<"all_positions visits "<<visits<<" pruned "<<pruned<<" prefixes "<<prefixes<<" hits "<<hits<<'\n';return 0;}
 for(bulk=1;bulk<=3;bulk+=2){V v{};for(int m=0;m<3;m++){V u{};u[m*36]=u[m*36+18]=1;for(int j=0;j<=m;j++)u=step(u,bulk);for(int i=0;i<108;i++)v[i]+=u[i];}dfs(v,0,15);cerr<<"bulk "<<bulk<<" visits "<<visits<<" pruned "<<pruned<<" hits "<<hits<<'\n';}
 return 0;
}
