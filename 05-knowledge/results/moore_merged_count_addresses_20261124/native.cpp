#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <vector>
using namespace std;
using U=uint64_t;
using V=vector<int>;
struct TT {
 int n;
 V o;
 TT(int n,U x):n(n),o(n) {
  for(int i=0,k=0;i<n;i++)for(int j=i+1;j<n;j++,k++)
   if(x>>k&1)o[j]|=1<<i;else o[i]|=1<<j;
 }
 U read(const V& p) const {
  U m=0;for(int i=0,k=0;i<n;i++)for(int j=i+1;j<n;j++,k++)if((o[p[j]]>>p[i])&1)m|=U(1)<<k;
  return m;
 }
 U canon(vector<V> cells) const {
  while(true) {
   vector<V> more;
   for(auto c:cells) {
    map<V,V> split;
    for(int v:c) {
     V key;
     for(auto d:cells){int z=0;for(int j:d)z+=(o[v]>>j)&1;key.push_back(z);}
     split[key].push_back(v);
    }
    for(auto &t:split)more.push_back(t.second);
   }
   if(more.size()==cells.size())break;
   cells=move(more);
  }
  if((int)cells.size()==n){V p;for(auto& c:cells)p.push_back(c[0]);return read(p);}
  int j=0;while(cells[j].size()<2)++j;
  U best=~U(0);
  for(int v:cells[j]) {
   auto c=cells;auto &r=c[j];r.erase(find(r.begin(),r.end(),v));c.insert(c.begin()+j,V{v});
   best=min(best,canon(move(c)));
  }
  return best;
 }
 U key() const {V p(n);iota(p.begin(),p.end(),0);return canon({p});}
 int paths() const {
  int d[1024][10]={},ret=0;
  for(int m=1;m<(1<<n);m++)for(int b=m;b;b&=b-1){
   int v=__builtin_ctz((unsigned)b), z=m^(1<<v);if(!z)d[m][v]=1;
   for(int q=z&~o[v];q;q&=q-1)d[m][v]+=d[z][__builtin_ctz((unsigned)q)];
  }
  for(int i=0;i<n;i++)ret+=d[(1<<n)-1][i];return ret;
 }
};
int main(){
 for(int n:{5,10}) {
  int G[10][10]={};
  auto edge=[&](int a,int b){G[a][b]=G[b][a]=1;};
  for(int i=0;i<5;i++){edge(i,(i+1)%5);if(n==10){edge(5+i,5+(i+2)%5);edge(i,5+i);}}
  V p(n);iota(p.begin(),p.end(),0);
  set<U> marks[2], blueMarks[2];map<U,V> layout;U opt=0,all=(U(1)<<(n*(n-1)/2))-1;
  for(int i=0,k=0;i<n;i++)for(int j=i+1;j<n;j++,k++)if(j>i+1)opt|=U(1)<<k;
  do{
   int gp=1,cp=1;
   for(int i=1;i<n;i++){gp&=G[p[i-1]][p[i]];cp&=1-G[p[i-1]][p[i]];}
   if(!gp&&!cp)continue;
   U x=0,z=0;
   for(int i=0,k=0;i<n;i++)for(int j=i+1;j<n;j++,k++)if(j>i+1){
    if(G[p[i]][p[j]])x|=U(1)<<k;
    if(G[p[n-1-i]][p[n-1-j]])z|=U(1)<<k;
   }
   marks[cp].insert(x);layout.try_emplace(x,p);if(x==z)blueMarks[cp].insert(x);
  }while(next_permutation(p.begin(),p.end()));
  for(int s=0;s<2;s++){
   using N=pair<U,bool>;
   map<U,int> H;
   set<U> from,to;
   set<pair<U,U>> line;
   set<int> src,dst;
   set<U> self;map<U,U> sample;map<pair<U,U>, U> evidence;
   int loops=0, diff2=0;
   auto node=[&](U x){
    TT t(n,x);U a=t.key(),b=TT(n,x^all).key(),z=min(a,b);
    if(!H.count(z))H[z]=t.paths();
    if(a==b)self.insert(z);
    return z;
   };
   for(U x:marks[s]){
    U u=node(x),v=node(x^opt);from.insert(u);to.insert(v);line.insert(minmax(u,v));sample.try_emplace(u,x);sample.try_emplace(v,x^opt);evidence.try_emplace(minmax(u,v),x);
    src.insert(H[u]);dst.insert(H[v]);
   }
   for(auto [u,v]:line){loops+=u==v;diff2+=abs(H[u]-H[v])==2;}
   set<U> both;set_intersection(from.begin(),from.end(),to.begin(),to.end(),inserter(both,both.end()));
   int sca=0,scb=0;for(auto a:from)sca+=self.count(a);for(auto a:to)scb+=self.count(a);
   cout<<n<<" "<<(s?"complement-spine":"graph-spine")<<" marks/blue "<<marks[s].size()<<" "<<blueMarks[s].size()
   <<" source/mate/intersection "<<from.size()<<" "<<to.size()<<" "<<both.size()
   <<" sc "<<sca<<" "<<scb<<" edges/loops/+2 "<<line.size()<<" "<<loops<<" "<<diff2<<"\n";
   if(!src.empty())cout<<"source counts "<<src.size()<<" ["<<*src.begin()<<","<<*src.rbegin()<<"] mate "<<dst.size()<<" ["<<*dst.begin()<<","<<*dst.rbegin()<<"]\n";
   if(n<=5||s==0)for(auto [u,v]:line)cout<<"line "<<u<<":"<<H[u]<<" -- "<<v<<":"<<H[v]<<"\n";
   if(n==10 && s==1){
    for(auto [u,v]:line)if(abs(H[u]-H[v])==2){
     U x=evidence[{u,v}],t=node(x)==u?v:u;
     cout<<"seam "<<u<<":"<<H[u]<<" -- "<<v<<":"<<H[v]<<" native H "<<H[node(x)]<<" native mask "<<x<<" order";
     for(int i:layout[x])cout<<" "<<i;
     int physical=0;
     for(int a=0;a<n*(n-1)/2;a++) {
      U y=x^(U(1)<<a);
      physical+=min(TT(n,y).key(),TT(n,y^all).key())==t;
     }
     cout<<" physical one-arc realizations "<<physical<<"\n";
    }
    for(U u:self){
     TT t(n,sample[u]);set<U> tiles;V p;
     auto walk=[&](auto &me,int used)->void{
      if((int)p.size()==n){tiles.insert(t.read(p));return;}
      int out=(p.empty()?(1<<n)-1:t.o[p.back()])&~used;
      for(int z=out;z;z&=z-1){int j=__builtin_ctz((unsigned)z);p.push_back(j);me(me,used|(1<<j));p.pop_back();}
     };walk(walk,0);
     int blue=0;
     for(U y:tiles){
      U yy=0;V q(n);iota(q.rbegin(),q.rend(),0);yy=TT(n,y).read(q)^all;
      blue+=y==yy;
     }
     cout<<"self-converse node "<<u<<" H "<<H[u]<<" tiles "<<tiles.size()<<" blue "<<blue<<"\n";
    }
   }
  }
 }
}
