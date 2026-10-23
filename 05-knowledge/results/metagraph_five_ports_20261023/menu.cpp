// Exact rooted transitive-five incidence menu on the complete inherited classes.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
struct T {uint64_t code; int aut; array<int,9> a{}; array<array<int,7>,9> f{};};
uint64_t checks=0;
void need(bool b,string s){checks++;if(!b)throw runtime_error(s);}
vector<int> subsets(int n,int k){vector<int>v;for(int s=0;s<(1<<n);s++)if(__builtin_popcount((unsigned)s)==k)v.push_back(s);return v;}
struct Menu{string name;vector<int>cols;};
vector<Menu>menus={{"v4",{0}},{"source5",{1}},{"source_sink5",{1,5}},{"containing5",{6}},{"roles5",{1,2,3,4,5}},{"v4_source5",{0,1}},{"v4_source_sink5",{0,1,5}},{"v4_containing5",{0,6}},{"v4_roles5",{0,1,2,3,4,5}}};
uint64_t fnv(uint64_t h,uint64_t z){for(int k=0;k<8;k++){h^=(z>>(8*k))&255;h*=1099511628211ULL;}return h;}
int main(int argc,char**argv){try{
 need(argc==4,"usage menu n input_code_aut.tsv output.json");int n=stoi(argv[1]);need(n>=5&&n<=9,"order");ifstream in(argv[2]);need(bool(in),"input");vector<T> ts;uint64_t code,mass=0;int aut;
 auto sub4=subsets(n,4),sub5=subsets(n,5);
 while(in>>code>>aut){T t;t.code=code;t.aut=aut;int bit=n*(n-1)/2-1;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,bit--){if(code>>bit&1)t.a[i]|=1<<j;else t.a[j]|=1<<i;}
  for(int s:sub4){array<int,4>d;int k=0;for(int i=0;i<n;i++)if(s>>i&1)d[k++]=__builtin_popcount((unsigned)(t.a[i]&s));sort(d.begin(),d.end());if(d==array<int,4>{1,1,2,2})for(int i=0;i<n;i++)if(s>>i&1)t.f[i][0]++;}
  for(int s:sub5){int seen=0;array<int,9>d{};for(int i=0;i<n;i++)if(s>>i&1){d[i]=__builtin_popcount((unsigned)(t.a[i]&s));seen|=1<<d[i];}if(seen==31)for(int i=0;i<n;i++)if(s>>i&1){t.f[i][5-d[i]]++;t.f[i][6]++;}}
  for(int i=0;i<n;i++)need(accumulate(t.f[i].begin()+1,t.f[i].begin()+6,0)==t.f[i][6],"five roles sum");
  uint64_t fact=1;for(int i=1;i<=n;i++)fact*=i;need(fact%aut==0,"aut divisibility");mass+=fact/aut;ts.push_back(t);
 }
 need(mass==(1ULL<<(n*(n-1)/2)),"complete orbit mass");vector<uint64_t> codes;for(auto&t:ts)codes.push_back(t.code);sort(codes.begin(),codes.end());need(adjacent_find(codes.begin(),codes.end())==codes.end(),"duplicate input");
 ofstream out(argv[3]);out<<"{\"n\":"<<n<<",\"classes\":"<<ts.size()<<",\"labelled_mass\":"<<mass<<",\"menus\":[";bool comma=false;
 for(auto&menu:menus){vector<array<int,9>>color(ts.size());map<vector<int>,int> bank;
  for(size_t j=0;j<ts.size();j++)for(int v=0;v<n;v++){vector<int>key;for(int c:menu.cols)key.push_back(ts[j].f[v][c]);auto [it,fresh]=bank.emplace(key,bank.size());color[j][v]=it->second;}
  for(int round=0;round<=2;round++){
   map<vector<int>,vector<int>> groups;for(size_t j=0;j<ts.size();j++){vector<int>k(color[j].begin(),color[j].begin()+n);sort(k.begin(),k.end());groups[k].push_back(j);}
   vector<vector<uint64_t>>bad;int largest=1,badclasses=0;uint64_t digest=1469598103934665603ULL;vector<vector<uint64_t>>all;
   for(auto&[key,ids]:groups){vector<uint64_t>x;for(int i:ids)x.push_back(ts[i].code);sort(x.begin(),x.end());all.push_back(x);if(x.size()>1){bad.push_back(x);largest=max(largest,(int)x.size());badclasses+=x.size();}}
   sort(all.begin(),all.end());for(auto&g:all){digest=fnv(digest,g.size());for(auto x:g)digest=fnv(digest,x);}sort(bad.begin(),bad.end());
   if(comma)out<<',';comma=true;out<<"{\"name\":\""<<menu.name<<"\",\"round\":"<<round<<",\"signatures\":"<<groups.size()<<",\"collision_groups\":"<<bad.size()<<",\"collision_classes\":"<<badclasses<<",\"largest_group\":"<<largest<<",\"partition_fnv64\":\""<<digest<<"\",\"collisions\":[";
   size_t limit=bad.size()<=500?bad.size():32;for(size_t i=0;i<limit;i++){if(i)out<<',';out<<'[';for(size_t j=0;j<bad[i].size();j++){if(j)out<<',';out<<bad[i][j];}out<<']';}out<<"],\"all_collisions_recorded\":"<<(limit==bad.size()?"true":"false")<<'}';
   cerr<<n<<' '<<menu.name<<' '<<round<<' '<<groups.size()<<" groups bad "<<bad.size()<<" classes "<<badclasses<<'\n';
   if(round<2){map<vector<int>,int> nextbank;vector<array<int,9>>next(ts.size());for(size_t j=0;j<ts.size();j++)for(int v=0;v<n;v++){vector<int>key{color[j][v]};vector<int>neighbors;for(int w=0;w<n;w++)if(ts[j].a[v]>>w&1)neighbors.push_back(color[j][w]);sort(neighbors.begin(),neighbors.end());key.insert(key.end(),neighbors.begin(),neighbors.end());auto[it,fresh]=nextbank.emplace(key,nextbank.size());next[j][v]=it->second;}color.swap(next);}
  }
 }
 out<<"],\"checks\":"<<checks<<"}\n";
 }catch(exception&e){cerr<<e.what()<<'\n';return 1;}return 0;}
