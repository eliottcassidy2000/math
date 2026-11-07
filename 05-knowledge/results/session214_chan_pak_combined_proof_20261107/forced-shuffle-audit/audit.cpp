#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>
using namespace std;
using Word=vector<int>;
void need(bool b,const string&s){if(!b)throw runtime_error(s);}
struct Tour {int n; vector<int> a; vector<vector<Word>> words;
 Tour(int n_,uint64_t code):n(n_),a(n),words(1<<n){int k=0; for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++) if((code>>k)&1)a[i]|=1<<j;else a[j]|=1<<i; Word p;words[0].push_back(p); grow(p,0,-1);}
 void grow(Word&p,int mask,int last){int avail=((1<<n)-1)^mask;if(last>=0)avail&=a[last]; while(avail){int v=__builtin_ctz(avail);avail&=avail-1;p.push_back(v);words[mask|1<<v].push_back(p);grow(p,mask|1<<v,v);p.pop_back();}}
 bool edge(int u,int v)const{return(a[u]>>v)&1;}
 bool strong()const{for(int v=0;v<n;v++){int r=1<<v,old=0;while(r!=old){old=r;for(int i=0;i<n;i++)if((r>>i)&1)r|=a[i];}if(r!=(1<<n)-1)return false;}return true;}
};
Word merge(const Tour&t,const Word&a,const Word&b,int i,int j){Word c;int u=0,v=0;while(u<i&&v<j){if(t.edge(a[u],b[v]))c.push_back(a[u++]);else c.push_back(b[v++]);}while(u<i)c.push_back(a[u++]);while(v<j)c.push_back(b[v++]);return c;}
Word inject(const Tour&t,const Word&a,const Word&b){int i=-1;for(int z=0;z<(int)a.size();z++)for(int w:b)if(t.edge(a[z],w))i=z;need(i>=0,"no forward arc");int j=0;while(!t.edge(a[i],b[j]))j++;Word c=merge(t,a,b,i,j);c.push_back(a[i]);for(int z=j;z<(int)b.size();z++)c.push_back(b[z]);for(int z=i+1;z<(int)a.size();z++)c.push_back(a[z]);return c;}
uint64_t child_count(const vector<int>&a,int n,int mask,int last){if(mask==(1<<n)-1)return 1;int avail=((1<<n)-1)^mask;if(last>=0)avail&=a[last];uint64_t count=0;while(avail){int v=__builtin_ctz(avail);avail&=avail-1;count+=child_count(a,n,mask|1<<v,v);}return count;}
int main(){map<int,set<int>>core;map<array<int,3>,uint64_t>multiplicity;uint64_t all_pairs=0,all_cuts=0,all_children=0;cout<<"{\"status\":\"PASS\",\"by_order\":[";for(int n=1;n<=6;n++){uint64_t nt=1ULL<<(n*(n-1)/2),nc=0,np=0,neg=0,nchild=0,nstrong=0;for(uint64_t code=0;code<nt;code++){Tour t(n,code);int full=(1<<n)-1,h=t.words[full].size();need(h%2==1,"Redei parity");bool strong=t.strong();nstrong+=strong;if(strong&&n>=3)for(int v=0;v<n;v++){int h0=t.words[full^(1<<v)].size(),delta=h-h0;need(delta>0&&delta%2==0,"strong deletion increment");need(2*h>=3*h0+1,"positive deletion floor");if(delta<=8){core[delta/2].insert(h0);multiplicity[{n,delta/2,h0}]++;}}
 for(int am=1;am<full;am++){int bm=full^am;bool forward=false;for(int u=0;u<n;u++)if((am>>u)&1)if(t.a[u]&bm)forward=true;int q=0;for(const auto&p:t.words[full]){bool sawA=false,ordered=true;for(int v:p){if((am>>v)&1)sawA=true;else if(sawA)ordered=false;}q+=ordered;}if(!forward){need(q==h,"all backward cut");neg++;continue;}nc++;set<Word>images;for(const auto&a:t.words[am])for(const auto&b:t.words[bm]){Word c=inject(t,a,b);need((int)c.size()==n,"length");for(int z=1;z<n;z++)need(t.edge(c[z-1],c[z]),"invalid edge");Word aa,bb;bool sawA=false,ordered=true;for(int v:c)if((am>>v)&1){aa.push_back(v);sawA=true;}else{bb.push_back(v);if(sawA)ordered=false;}need(aa==a&&bb==b,"subsequence inverse");need(!ordered,"ordered image");need(images.insert(c).second,"collision");np++;}
 int P=t.words[am].size()*t.words[bm].size();need((int)images.size()==P,"product cardinality");need(h>=P+q&&2*q<=h-1,"shuffle floor");
 if(n<=5){vector<int>z=t.a;z.push_back(am);for(int b=0;b<n;b++)if((bm>>b)&1)z[b]|=1<<n;auto child=child_count(z,n+1,0,-1);need(child>=uint64_t(2*h-q),"slot floor");need(2*child>=uint64_t(3*h+1),"positive growth floor");nchild++;}}
 }
 if(n>1) { cout<<","; }
 cout<<"{\"n\":"<<n<<",\"labeled_tournaments\":"<<nt<<",\"strong_tournaments\":"<<nstrong<<",\"forward_cuts\":"<<nc<<",\"backward_controls\":"<<neg<<",\"injected_pairs\":"<<np<<",\"fresh_children\":"<<nchild<<"}";all_pairs+=np;all_cuts+=nc;all_children+=nchild;}
 cout<<"],\"core_parent_sets\":{";bool first=true;for(auto[k,s]:core){if(!first)cout<<",";first=false;cout<<"\""<<k<<"\":[";bool f=true;for(int x:s){if(!f)cout<<",";f=false;cout<<x;}cout<<"]";}cout<<"},\"all_order_parent_sets_after_external_factors\":{";first=true;for(int k=1;k<=4;k++){if(!first)cout<<",";first=false;set<int>s;for(int g:{1,3})if(k%g==0)for(int x:core[k/g])s.insert(g*x);cout<<"\""<<k<<"\":[";bool f=true;for(int x:s){if(!f)cout<<",";f=false;cout<<x;}cout<<"]";}cout<<"},\"strong_core_root_counts\":[";first=true;for(auto[key,value]:multiplicity){if(!first)cout<<",";first=false;cout<<"{\"order\":"<<key[0]<<",\"k\":"<<key[1]<<",\"parent\":"<<key[2]<<",\"labeled_roots\":"<<value<<"}";}cout<<"],\"totals\":{\"cuts\":"<<all_cuts<<",\"pairs\":"<<all_pairs<<",\"children\":"<<all_children<<"}}\n";
}
