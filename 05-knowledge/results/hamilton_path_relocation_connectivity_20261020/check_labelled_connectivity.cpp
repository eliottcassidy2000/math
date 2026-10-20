// Independent direct-arc audit. Does not import candidate or prior checker code.
// Enumerates every labelled tournament and every Hamilton word for n <= 6.
#ifdef NDEBUG
#error "This verifier requires assertions; do not compile with -DNDEBUG."
#endif
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <vector>
using namespace std;
using Word = vector<int>;
bool A[8][8];
uint64_t enc(const Word &p) { uint64_t z=0; for (int x:p) z=(z<<3)|(x+1); return z; }
bool path(const Word &p) { for (size_t i=1;i<p.size();++i) if(!A[p[i-1]][p[i]]) return false; return true; }
Word erasev(Word p,int v) { auto it=find(p.begin(),p.end(),v); assert(it!=p.end());p.erase(it);return p; }
Word ins(Word p,int gap,int v) { p.insert(p.begin()+gap,v);return p; }
vector<Word> paths(const Word &vs) { vector<Word> out; Word p=vs; sort(p.begin(),p.end());do { if(path(p))out.push_back(p); }while(next_permutation(p.begin(),p.end()));return out; }
Word double_ins(const Word&r,int gu,int gv,int u,int v) {
 Word out;
 for(int g=0;g<=(int)r.size();++g) {
  if(g==gu && g==gv) { if(A[u][v]) {out.push_back(u);out.push_back(v);}else{out.push_back(v);out.push_back(u);} }
  else {if(g==gu)out.push_back(u);if(g==gv)out.push_back(v);}
  if(g<(int)r.size())out.push_back(r[g]);
 }
 return out;
}
int main(int argc,char**argv) {
 int maxn=argc>1?stoi(argv[1]):6;
 if(maxn<1 || maxn>6) { cerr<<"Supported exhaustive audit range is 1 <= n <= 6.\n"; return 2; }
 cout<<"n,tournaments,labelled_words,directed_relocation_edges,front_orbits,max_front_cycle,base_edge_lifts,shared_gap_lifts\n";
 for(int n=1;n<=maxn;++n) {
  uint64_t totals=0,edges=0,orbits=0,lifts=0,shared=0; int maxcyc=0;
  const int B=n*(n-1)/2, NT=1<<B;
  Word vs(n);iota(vs.begin(),vs.end(),0);
  for(int mask=0;mask<NT;++mask) {
   int b=0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j){A[i][j]=(mask>>b++)&1;A[j][i]=!A[i][j];}
   auto hp=paths(vs);assert(!hp.empty());totals+=hp.size();
   unordered_map<uint64_t,int> ids;for(int i=0;i<(int)hp.size();++i)ids.emplace(enc(hp[i]),i);
   vector<vector<int>> adj(hp.size());
   for(int i=0;i<(int)hp.size();++i) {
    const Word&p=hp[i]; unordered_set<int> nbr;
    for(int old=0;old<n;++old) {
     int u=p[old];Word r=erasev(p,u);
     for(int g=0;g<n;++g)if(g!=old) {
      Word q=ins(r,g,u);if(!path(q))continue;
      assert(path(r)); // Audit common-deletion lemma, not assumed in edge generation.
      assert(abs(g-old)>=2); // Adjacent swaps cannot be legal in a tournament.
      auto it=ids.find(enc(q));assert(it!=ids.end());nbr.insert(it->second);
     }
    }
    adj[i].assign(nbr.begin(),nbr.end());edges+=nbr.size();
   }
   vector<int> seen(hp.size(),0),todo{0};seen[0]=1;
   for(size_t k=0;k<todo.size();++k)for(int j:adj[todo[k]])if(!seen[j]){seen[j]=1;todo.push_back(j);}
   assert(todo.size()==hp.size());
   // Recover initial SCC from transitive closure, independently of path order.
   bool reach[8][8];for(int i=0;i<n;++i)for(int j=0;j<n;++j)reach[i][j]=(i==j)||A[i][j];
   for(int k=0;k<n;++k)for(int i=0;i<n;++i)for(int j=0;j<n;++j)reach[i][j]=reach[i][j]||(reach[i][k]&&reach[k][j]);
   unsigned C=0;for(int i=0;i<n;++i){bool initial=true;for(int j=0;j<n;++j)if(!reach[i][j])initial=false;if(initial)C|=1u<<i;}
   assert(C);int csize=__builtin_popcount(C);
   for(const Word&start:hp) {
    for(int i=0;i<n;++i)assert(bool(C&(1u<<start[i]))==(i<csize));
    if(csize==1)continue;
    ++orbits;unordered_map<uint64_t,int> visited;vector<int> first;Word p=start;
    while(!visited.count(enc(p))) {
     visited[enc(p)]=first.size();first.push_back(p[0]);int j=-1;
     for(int k=1;k<csize;++k)if(A[p[k]][p[0]])j=k;
     assert(j>=2);int u=p[0];Word r=erasev(p,u);Word q=ins(r,j,u);
     assert(path(q));assert(q[0]==p[1]);
     p=q;
    }
    int begin=visited.at(enc(p));unsigned recurrent=0;
    for(int i=begin;i<(int)first.size();++i)recurrent|=1u<<first[i];
    assert(recurrent==C);maxcyc=max(maxcyc,(int)first.size()-begin);
   }
   // Every deleted vertex, every base edge, and EVERY valid insertion gap of v in R.
   for(int v=0;v<n;++v) {
    auto base=paths(erasev(vs,v));
    for(const Word&p:base)for(int old=0;old<n-1;++old) {
     int u=p[old];Word r=erasev(p,u);
     for(int g=0;g<n-1;++g)if(g!=old) {
      Word q=ins(r,g,u);if(!path(q)||enc(q)<enc(p))continue;
      assert(path(r));
      int count=0;
      for(int gv=0;gv<=(int)r.size();++gv) {
       Word rv=ins(r,gv,v);if(!path(rv))continue;++count;
       Word pp=double_ins(r,old,gv,u,v),qq=double_ins(r,g,gv,u,v);
       assert(path(pp)&&path(qq)&&pp!=qq);
       assert(erasev(pp,v)==p&&erasev(qq,v)==q);
       assert(erasev(pp,u)==rv&&erasev(qq,u)==rv);
       ++lifts;if(gv==old||gv==g)++shared;
      }
      assert(count>0);
     }
    }
   }
  }
  uint64_t fact=1;for(int k=2;k<=n;++k)fact*=k;
  assert(totals==fact*(1ull<<(B-(n-1))));
  cout<<n<<','<<NT<<','<<totals<<','<<edges<<','<<orbits<<','<<maxcyc<<','<<lifts<<','<<shared<<'\n'<<flush;
 }
 cout<<"PASS: all labelled relocation graphs connected; all front recurrences and all edge lifts valid.\n";
}
