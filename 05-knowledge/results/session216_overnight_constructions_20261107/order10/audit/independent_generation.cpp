#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <unordered_set>
using namespace std;using Code=uint64_t;using Graph=vector<unsigned>;using Cells=vector<vector<int>>;
void must(bool b){if(!b)throw runtime_error("independent generation failure");}
Graph decode(Code c,int n){Graph a(n);int t=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,t++){if(c>>t&1)a[i]|=1u<<j;else a[j]|=1u<<i;}return a;}
Code bits(const Graph&a,const vector<int>&p){Code c=0;int t=0;for(int i=0;i<(int)p.size();i++)for(int j=i+1;j<(int)p.size();j++,t++)if(a[p[i]]>>p[j]&1)c|=Code(1)<<t;return c;}
// Independent orientation convention: incoming counts; initial degree partition;
// largest nonsingleton cell; distinguished vertex is LAST, not first.
Code canon(const Graph&a){int n=a.size();vector<unsigned> incoming(n);for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(a[j]>>i&1)incoming[i]|=1u<<j;map<int,vector<int>> degrees;for(int i=0;i<n;i++)degrees[__builtin_popcount(incoming[i])].push_back(i);Cells init;for(auto&kv:degrees)init.push_back(kv.second);Code best=~Code(0);
 auto recurse=[&](auto&& self,Cells cells)->void{
  while(true){vector<unsigned> masks;for(auto&cell:cells){unsigned m=0;for(int v:cell)m|=1u<<v;masks.push_back(m);}Cells refined;for(auto&cell:cells){map<vector<int>,vector<int>> split;for(int v:cell){vector<int>sig;for(auto m:masks)sig.push_back(__builtin_popcount(incoming[v]&m));split[sig].push_back(v);}for(auto&kv:split)refined.push_back(kv.second);}bool stable=refined.size()==cells.size();cells=move(refined);if(stable)break;}
  int chosen=-1;for(int k=0;k<(int)cells.size();k++)if(cells[k].size()>1&&(chosen<0||cells[k].size()>cells[chosen].size()))chosen=k;
  if(chosen<0){vector<int>p;for(auto&cell:cells)p.push_back(cell[0]);best=min(best,bits(a,p));return;}
  auto branch=cells[chosen];for(int v:branch){Cells next=cells;vector<int>other;for(int u:branch)if(u!=v)other.push_back(u);next[chosen]=other;next.insert(next.begin()+chosen+1,vector<int>{v});self(self,move(next));}
 };recurse(recurse,init);return best;}
vector<Graph> read(const string&name,int desired){ifstream f(name);must(bool(f));int count;f>>count;vector<Graph>out;for(int k=0;k<count;k++){int n;f>>n;Graph a(n);for(auto&v:a)f>>v;must(bool(f));for(int i=0;i<n;i++){must(!(a[i]>>n)&&!(a[i]>>i&1));for(int j=0;j<i;j++)must((a[i]>>j&1)+(a[j]>>i&1)==1);}if(n==desired)out.push_back(a);}return out;}
int main(int argc,char**argv){long long small=0;for(int n=1;n<=5;n++){map<Code,Code>forward,backward;for(Code c=0;c<(Code(1)<<(n*(n-1)/2));c++){auto a=decode(c,n);auto actual=canon(a);vector<int>p(n);for(int i=0;i<n;i++)p[i]=i;Code brute=~Code(0);do{brute=min(brute,bits(a,p));}while(next_permutation(p.begin(),p.end()));must(!forward.count(actual)||forward[actual]==brute);must(!backward.count(brute)||backward[brute]==actual);forward[actual]=brute;backward[brute]=actual;small++;}}
 auto eight=read(argv[1],8),nine=read(argv[2],9);unordered_set<Code> classes,seen;for(auto&a:nine)must(classes.insert(canon(a)).second);must(classes.size()==191536);unordered_set<Code>prior;for(auto&a:eight)must(prior.insert(canon(a)).second);must(prior.size()==6880);long long ext=0;
 for(auto a:eight)for(unsigned mask=0;mask<256;mask++){Graph g=a;g.push_back(mask);for(int i=0;i<8;i++)if(!(mask>>i&1))g[i]|=256;Code c=canon(g);must(classes.count(c));seen.insert(c);ext++;}
 must(seen==classes);cout<<"{\"small_bruteforce\":"<<small<<",\"independent_order8_classes\":"<<prior.size()<<",\"order9_classes\":"<<classes.size()<<",\"all_augmentations_checked\":"<<ext<<",\"method\":\"incoming refinement, largest-cell branching, last singleton\"}"<<endl;
}
