// Exact cycle alphabets; counts are rooted through an oriented edge.
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <algorithm>
#include <vector>
using U=std::uint64_t;
std::array<U,50> adj;
std::array<long long,51> cnt;
int root;
void chordless(int len,int cur,U used,U blocked){
  U free=adj[cur]&~(used|blocked);
  while(free){int j=__builtin_ctzll(free);free&=free-1;
    if(adj[root]>>j&1) {if(len>=3) ++cnt[len+1];}
    else chordless(len+1,j,used|(1ULL<<j),blocked|adj[cur]);
  }
}
void simple(int len,int cur,U used){
 if(len>=3&&(adj[root]>>cur&1))++cnt[len];
 U f=adj[cur]&~used;
 while(f){int x=__builtin_ctzll(f);f&=f-1;simple(len+1,x,used|(1ULL<<x));}
}
void through(int r,int s,bool induced){
 root=r;cnt.fill(0);U m=(1ULL<<r)|(1ULL<<s);
 if(!induced){simple(2,s,m);return;}
 U f=adj[s]&~m;
 while(f){int x=__builtin_ctzll(f);f&=f-1;
   if(adj[r]>>x&1)++cnt[3];
   else chordless(3,x,m|(1ULL<<x),adj[s]);
 }
}
bool witness(std::vector<int>& p,int n,U used){
 if((int)p.size()==n)return adj[p.back()]>>p[0]&1;
 U f=adj[p.back()]&~used;
 while(f){int x=__builtin_ctzll(f);f&=f-1;p.push_back(x);
  if(witness(p,n,used|(1ULL<<x)))return true;
  p.pop_back();
 }
 return false;
}
void add(int i,int j){adj[i]|=1ULL<<j;adj[j]|=1ULL<<i;}
void check_moore(int n,int d){
 for(int i=0;i<n;i++){
  assert(__builtin_popcountll(adj[i])==d);
  for(int j=0;j<i;j++)
   assert(__builtin_popcountll(adj[i]&adj[j])==!(adj[i]>>j&1));
 }
}
int main(){
 for(int i=0;i<5;i++){add(i,(i+1)%5);add(5+i,5+(i+2)%5);add(i,i+5);}
 check_moore(10,3);
 for(int i=0;i<10;i++)for(int j=0;j<i;j++)if(adj[i]>>j&1){
   through(i,j,false);
   for(int k=3;k<=10;k++)assert((cnt[k]>0)==(k==5||k==6||k==8||k==9));
   through(i,j,true);
   for(int k=3;k<=10;k++)assert((cnt[k]>0)==(k==5||k==6));
 }
 adj.fill(0);
 for(int a=0;a<5;a++)for(int x=0;x<5;x++){
  add(5*a+x,5*a+(x+1)%5);
  add(25+5*a+x,25+5*a+(x+2)%5);
  for(int b=0;b<5;b++)add(5*a+x,25+5*b+(x+a*b)%5);
 }
 check_moore(50,7);
 const std::vector<int> sizes={5,6,7,8,9,10,12,13,16,18};
 const std::vector<int> hits={36,180,720,1440,2160,5040,7920,9360,10080,12960};
 std::array<long long,51> want{};for(unsigned i=0;i<sizes.size();i++)want[sizes[i]]=hits[i];
 for(int i=0;i<50;i++)for(int j=0;j<i;j++)if(adj[i]>>j&1){
   through(i,j,true);assert(cnt==want);
 }
 // One explicitly marked induced 13-cycle separates every vertex after
 // two rounds of color refinement. This certifies pointwise stabilizer 1.
 std::vector<int> w={0,1,2,3,28,8,7,6,44,15,16,17,30};
 for(int i=0;i<13;i++)for(int j=0;j<i;j++)
   assert(bool(adj[w[i]]>>w[j]&1)==(i-j==1||i-j==12));
 std::array<int,50> color{};for(int i=0;i<13;i++)color[w[i]]=i+1;
 for(int t=0;t<2;t++){
   std::array<int,50> next{};std::map<std::vector<int>,int> keys;
   for(int i=0;i<50;i++){
    std::vector<int> c={color[i]};
    for(int j=0;j<50;j++)if(adj[i]>>j&1)c.push_back(color[j]);
    std::sort(c.begin()+1,c.end());
    if(!keys.count(c)){int pos=keys.size();keys[c]=pos;}
    next[i]=keys.at(c);
   }
   assert(keys.size()==(t==0?44:50));color=next;
 }
 assert(9360*175/13==126000);
 std::vector<int> p={0};assert(witness(p,21,1));
 std::cout<<"HS ordinary 21-cycle:";for(int x:p)std::cout<<" "<<x;
 std::cout<<"\nExact Petersen and HS induced spectra checked on every edge.\n";
}
