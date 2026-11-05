#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <chrono>
#include <stdexcept>
using namespace std;
struct G {vector<int> r;vector<vector<int>> cs;};
int n,full;vector<int> row;vector<int8_t> mem;vector<unsigned> stamp;unsigned epoch=0;
int countc(int mask,int last){
 if(mask==full)return !!(row[last]&1);
 int ix=mask*n+last;if(stamp[ix]==epoch)return mem[ix];stamp[ix]=epoch;
 int todo=row[last]&~mask,ans=0;
 while(todo){int b=todo&-todo;todo-=b;ans+=countc(mask|b,__builtin_ctz((unsigned)b));if(ans>=4){ans=4;break;}}
 mem[ix]=ans;return ans;
}
void collect(int mask,int last,vector<int>&w,vector<vector<int>>&cs){
 if(mask==full){if(row[last]&1)cs.push_back(w);return;}
 int todo=row[last]&~mask;
 while(todo){int b=todo&-todo;todo-=b;int v=__builtin_ctz((unsigned)b);if(countc(mask|b,v)){w.push_back(v);collect(mask|b,v,w,cs);w.pop_back();}}
}
string key(const vector<int>&r,const vector<vector<int>>&cs){
 int z=r.size();string best;vector<int>inv(z);
 for(auto &c:cs)for(int s=0;s<z;s++){
  for(int j=0;j<z;j++)inv[c[(s+j)%z]]=j;
  string k(2*z,'\0');
  for(int j=0;j<z;j++){int rr=0,u=c[(s+j)%z];for(int v=0;v<z;v++)if(r[u]>>v&1)rr|=1<<inv[v];k[2*j]=rr&255;k[2*j+1]=rr>>8;}
  if(best.empty()||k<best)best=k;
 }
 return best;
}
G decode(const string&k){G g;int z=k.size()/2;for(int i=0;i<z;i++)g.r.push_back((unsigned char)k[2*i]+((unsigned char)k[2*i+1]<<8));return g;}
void getcycles(G&g){for(int u=0;u<(int)g.r.size();u++)for(int v=0;v<(int)g.r.size();v++){if(u==v){if(g.r[u]>>u&1)throw runtime_error("loop");}else if(((g.r[u]>>v&1)+(g.r[v]>>u&1))!=1)throw runtime_error("not tournament");}row=g.r;n=row.size();full=(1<<n)-1;epoch++;int c=countc(1,0);if(c<1||c>3)throw runtime_error("bad retained count");vector<int>w={0};collect(1,0,w,g.cs);if((int)g.cs.size()!=c)throw runtime_error("cycle collection mismatch");}
int main(int argc,char**argv){if(argc!=2)throw runtime_error("usage: verifier graph-output-file");const int maxn=12;mem.resize((1<<maxn)*maxn);stamp.resize(mem.size());
 vector<G>gs;G base;base.r={2,4,9,3};getcycles(base);base=decode(key(base.r,base.cs));getcycles(base);gs.push_back(base);
 ofstream records(argv[1]);if(!records)throw runtime_error("cannot open graph records");long long candidates=0,filtered=0,kept=0;
 for(int order=4;order<=maxn;order++){
  long long hist[4]={},bad=0;for(auto &g:gs){hist[g.cs.size()]++;bool one=false;for(int rr:g.r){int d=__builtin_popcount((unsigned)rr);if(d==1||d==order-2)one=true;}if(g.cs.size()==3&&!one){bad++;cerr<<"C3 ATOM CANDIDATE n="<<order<<" rows";for(int rr:g.r)cerr<<" "<<rr;cerr<<"\n";}}
  vector<pair<string,const G*>> ordered;for(auto &g:gs)ordered.push_back({key(g.r,g.cs),&g});sort(ordered.begin(),ordered.end(),[](auto&a,auto&b){return a.first<b.first;});
  for(auto &kv:ordered){auto &g=*kv.second;records<<order<<","<<g.cs.size();for(int r:g.r)records<<","<<r;records<<"\n";}
  cout<<"{\"order\":"<<order<<",\"C1\":"<<hist[1]<<",\"C2\":"<<hist[2]<<",\"C3\":"<<hist[3]<<",\"no_degree_one_C3\":"<<bad<<",\"extension_candidates\":"<<candidates<<",\"insertion_filter_survivors\":"<<filtered<<",\"retained_labelled_extensions\":"<<kept<<"}"<<endl;
  if(bad)throw runtime_error("C3 atom candidate found");
  if(order==maxn)break;
  unordered_map<string,G>next;candidates=filtered=kept=0;
  n=order+1;full=(1<<n)-1;
  for(auto &g:gs){
   for(int mask=1;mask<(1<<order)-1;mask++){
    candidates++;int lower=0;
    for(auto &cy:g.cs){for(int j=0;j<order;j++)if(!(mask>>cy[j]&1)&&(mask>>cy[(j+1)%order]&1))lower++;if(lower>3)break;}
    if(lower>3)continue;filtered++;
    row=g.r;row.push_back(mask);for(int j=0;j<order;j++)if(!(mask>>j&1))row[j]|=1<<order;
    epoch++;int c=countc(1,0);if(c>3)continue;if(c<1)throw runtime_error("mixed extension not strong");kept++;
    vector<vector<int>>cs;vector<int>w={0};collect(1,0,w,cs);if((int)cs.size()!=c)throw runtime_error("wrong cycle list");
    string k=key(row,cs);if(next.find(k)==next.end()){G h=decode(k);next.emplace(k,move(h));}
   }
  }
  gs.clear();gs.reserve(next.size());for(auto &kv:next){G h=move(kv.second);getcycles(h);gs.push_back(move(h));}

 }
}
