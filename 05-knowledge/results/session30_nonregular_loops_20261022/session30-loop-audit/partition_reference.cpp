#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;using I=long long;
I fact(int a){I f=1;for(int i=2;i<=a;i++)f*=i;return f;}
struct Audit{
 int L,n,scale;I S[4][4]{}; vector<int> runs,p,prev;vector<pair<int,int>> arcs;
 map<pair<int,int>,I> sums;unordered_map<string,I> cache;I connected=0,selected=0;
 I hom(vector<pair<int,int>> es,int V){
  int sign=1;for(auto &e:es)if(e.first>e.second){swap(e.first,e.second);sign=-sign;}sort(es.begin(),es.end());
  string key(1,char(V));for(auto[a,b]:es){key+=char(a);key+=char(b);}auto it=cache.find(key);if(it!=cache.end())return sign*it->second;
  int mult[12][12]={},deg[12]={};for(auto[a,b]:es){mult[a][b]++;mult[b][a]++;deg[a]++;deg[b]++;}
  vector<int> order;bool used[12]={};int first=max_element(deg,deg+V)-deg;order.push_back(first);used[first]=1;
  while((int)order.size()<V){int best=-1,bscore=-1;for(int v=0;v<V;v++)if(!used[v]){int score=0;for(int u:order)score+=mult[u][v];score=score*100+deg[v];if(score>bscore){best=v;bscore=score;}}order.push_back(best);used[best]=1;}
  int color[12]={};I total=0;
  function<void(int,I)> dfs=[&](int k,I val){if(k==V){total+=val;return;}int v=order[k];for(int c=0;c<n;c++){I next=val;for(int h=0;h<k;h++){int u=order[h],m=mult[u][v];if(!m)continue;I entry=u<v?S[color[u]][c]:S[c][color[u]];for(int j=0;j<m;j++)next*=entry;if(!next)break;}if(next){color[v]=c;dfs(k+1,next);}}};dfs(0,1);cache[key]=total;return sign*total;
 }
 void accept(int V){
  int deg[12]={},sz[12]={};vector<pair<int,int>> es;int parent[12];iota(parent,parent+V,0);function<int(int)> find=[&](int x){return parent[x]==x?x:parent[x]=find(parent[x]);};
  for(int a:p)sz[a]++;
  for(auto[a,b]:arcs){a=p[a];b=p[b];if(a==b)return;deg[a]++;deg[b]++;parent[find(a)]=find(b);es.emplace_back(a,b);}
  for(int j=1;j<V;j++)if(find(j)!=find(0))return;connected++;
  int ell=count(deg,deg+V,1),e=L-V,w=e+ell;if(w<1||w>2)return;selected++;
  I mu=1;for(int j=0;j<V;j++)mu*=((sz[j]-1)%2?-1:1)*fact(sz[j]-1);
  sums[{e,ell}]+=24/fact(runs.size())*mu*hom(es,V);
 }
 void partitions(int at,int maxlabel){if(at==(int)p.size()){accept(maxlabel+1);return;}for(int k=0;k<=maxlabel+1;k++){if(prev[at]>=0&&p[prev[at]]==k)continue;p[at]=k;partitions(at+1,max(maxlabel,k));}}
 void compositions(int remain){if(!remain){arcs.clear();prev.clear();int off=0;for(int a:runs){for(int j=0;j<=a;j++)prev.push_back(j?off+j-1:-1);for(int j=0;j<a;j++)arcs.emplace_back(off+j,off+j+1);off+=a+1;}p.assign(off,0);partitions(1,0);return;}for(int a=2;a<=remain;a+=2){runs.push_back(a);compositions(remain-a);runs.pop_back();}}
 void run(){compositions(L);cout<<"{\"L\":"<<L<<",\"n\":"<<n<<",\"scale\":"<<scale<<",\"connected\":"<<connected<<",\"selected\":"<<selected<<",\"cache\":"<<cache.size()<<",\"rows\":[";bool first=true;for(auto[key,num]:sums){if(!first)cout<<",";first=false;I den=24;for(int i=0;i<L;i++)den*=scale;I g=gcd(abs(num),den);cout<<"{\"excess\":"<<key.first<<",\"leaves\":"<<key.second<<",\"numerator\":"<<num/g<<",\"denominator\":"<<den/g<<"}";}cout<<"]}\n";}
};
int main(){for(int control=0;control<2;control++){Audit a;a.L=8;a.n=control?4:3;a.scale=control?4:7;if(!control){a.S[0][1]=2;a.S[0][2]=-1;a.S[1][2]=3;}else for(int i=0;i<4;i++)for(int j=i+1;j<4;j++)a.S[i][j]=1;for(int i=0;i<a.n;i++)for(int j=i+1;j<a.n;j++)a.S[j][i]=-a.S[i][j];a.run();}}
