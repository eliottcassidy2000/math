#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;using I=long long;using Rows=vector<int>;
void need(bool c,const char*m){if(!c)throw runtime_error(m);}
vector<I> starts(const Rows&a){int n=a.size(),M=1<<n;vector<I>dp(M*n),s(n);for(int i=0;i<n;i++)dp[(1<<i)*n+i]=1;for(int S=1;S<M;S++)for(int i=0;i<n;i++)if(S>>i&1){I z=dp[S*n+i];for(int j=0;j<n;j++)if(!(S>>j&1)&&(a[j]>>i&1))dp[(S|1<<j)*n+j]+=z;}for(int i=0;i<n;i++)s[i]=dp[(M-1)*n+i];return s;}
Rows sub(const Rows&a,const vector<int>&v){Rows b(v.size());for(int i=0;i<(int)v.size();i++)for(int j=0;j<(int)v.size();j++)if(a[v[i]]>>v[j]&1)b[i]|=1<<j;return b;}
vector<vector<int>> scc(const Rows&a,const vector<int>&v){Rows b=sub(a,v);int n=b.size();vector<int>reach(n);for(int i=0;i<n;i++){int S=1<<i,old=-1;while(S!=old){old=S;for(int j=0;j<n;j++)if(S>>j&1)S|=b[j];}reach[i]=S;}vector<vector<int>>cs;int used=0;for(int i=0;i<n;i++)if(!(used>>i&1)){vector<int>c;for(int j=0;j<n;j++)if((reach[i]>>j&1)&&(reach[j]>>i&1)){c.push_back(v[j]);used|=1<<j;}cs.push_back(c);}sort(cs.begin(),cs.end(),[&](auto&C,auto&D){return bool(a[C[0]]>>D[0]&1);});return cs;}
struct Word{const Rows&a;int K,W,maxlive=0,adds=0,forgets=0,steps=0;set<int>active,done;vector<int>order;Rows built;Word(const Rows&aa,int kk):a(aa),K(kk){auto m=[](int t){int z=0;while((z+1)*(z+1)<4*t)z++;return z;};W=0;for(int t=0;t<=K;t++)W=max(W,1+m(t)+m(K-t));}
 void add(int x){need(!done.count(x),"duplicate vertex");int row=0;for(int i=0;i<(int)order.size();i++){int y=order[i];if(a[x]>>y&1){need(active.count(y),"backward arc to forgotten vertex");row|=1<<i;}else built[i]|=1<<order.size();}built.push_back(row);order.push_back(x);done.insert(x);active.insert(x);maxlive=max(maxlive,(int)active.size());adds++;need(maxlive<=W,"register bound");}
 void forget(int x){need(active.erase(x)==1,"forget absent");forgets++;}
 void build(const vector<int>&v,const set<int>&keep){auto st=starts(sub(a,v));I budget=0;for(int i=0;i<(int)v.size();i++)if(keep.count(v[i]))budget+=st[i];need(budget<=K&&!keep.empty(),"marked budget");int root=*keep.begin();need(keep.size()*keep.size()<4*budget,"marked arity");if(v.size()==1){add(root);return;}vector<int>rest;for(int x:v)if(x!=root)rest.push_back(x);auto cs=scc(a,rest);auto C=cs[0];set<int>marks;for(int x:C)if(keep.count(x)||(a[root]>>x&1))marks.insert(x);need(!marks.empty(),"empty entry");build(C,marks);add(root);vector<int>drop;for(int x:active)if(!keep.count(x))drop.push_back(x);for(int x:drop)forget(x);for(int i=1;i<(int)cs.size();i++){need(cs[i].size()<=K,"large tail block");for(int x:cs[i])add(x);for(int x:cs[i])if(!keep.count(x))forget(x);}need(active==keep,"wrong final retained set");steps++;}
 void check(){vector<int>v(a.size());iota(v.begin(),v.end(),0);build(v,set<int>{0});need(done.size()==a.size(),"missing vertex");need(built==sub(a,order),"literal word mismatch");}
};
int main(int argc,char**argv){int limit=argc>2?stoi(argv[2]):1000000;ifstream f(argc>1?argv[1]:"session203-rooted-target-descent/run2001/ROOTED_BANK.tsv");need(bool(f),"bank not found");int n,s;I B,count=0,steps=0,adds=0,forgets=0;int maxlive=0,maxn=0,maxK=0;while(count<limit&&f>>n>>B>>s){Rows a(n);for(int&x:a)f>>x;auto actual=starts(a);need(actual[0]==s&&s<=26,"source start mismatch");vector<int>v(n);iota(v.begin(),v.end(),0);need(scc(a,v).size()==1,"source not strong");Word w(a,s);w.check();maxlive=max(maxlive,w.maxlive);maxn=max(maxn,n);maxK=max(maxK,s);steps+=w.steps;adds+=w.adds;forgets+=w.forgets;count++;}cout<<"{\n  \"status\": \"PASS\",\n  \"literal_sources\": "<<count<<",\n  \"largest_source_order\": "<<maxn<<",\n  \"largest_start_budget\": "<<maxK<<",\n  \"largest_observed_register_use\": "<<maxlive<<",\n  \"recursive_scc_steps\": "<<steps<<",\n  \"introductions\": "<<adds<<",\n  \"forget_operations\": "<<forgets<<",\n  \"scope\": \"Each supplied source independently start-counted, SCC-decomposed, register-compiled, and compared literally; bank completeness not used.\"\n}\n";}
