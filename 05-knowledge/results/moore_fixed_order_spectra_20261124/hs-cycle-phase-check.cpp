// Exact cycle alphabets; counts are rooted through an oriented edge.
#include <array>
#include <map>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <algorithm>
#include <vector>
using U=std::uint64_t;
std::array<U,50> adj;
void add(int i,int j){adj[i]|=1ULL<<j;adj[j]|=1ULL<<i;}
void check_moore(int n,int d){
 for(int i=0;i<n;i++){
  assert(__builtin_popcountll(adj[i])==d);
  for(int j=0;j<i;j++)
   assert(__builtin_popcountll(adj[i]&adj[j])==!(adj[i]>>j&1));
 }
}
// enumerate with minimum start, then choose one of two directions
#include <set>
#include <algorithm>
std::set<U> cycles;
void go(int r,int first,int cur,int len,U used,U forbidden) {
 U f=adj[cur]&~(forbidden|used)&~((1ULL<<r)-1);
 while(f){int j=__builtin_ctzll(f); f&=f-1;
  if(adj[r]>>j&1){if(len==12 && first<j)cycles.insert(used|(1ULL<<j));}
  else if(len<12)go(r,first,j,len+1,used|(1ULL<<j),forbidden|adj[cur]);}
}
int main(){
 adj.fill(0);
 for(int a=0;a<5;a++)for(int x=0;x<5;x++){
  add(5*a+x,5*a+(x+1)%5);
  add(25+5*a+x,25+5*a+(x+2)%5);
  for(int b=0;b<5;b++)add(5*a+x,25+5*b+(x+a*b)%5);
 }
 check_moore(50,7);
 for(int r=0;r<50;r++)for(int j=r+1;j<50;j++)if(adj[r]>>j&1)
  for(int k=r+1;k<50;k++)if(k!=j&&(adj[j]>>k&1)&&!(adj[r]>>k&1))go(r,j,k,3,(1ULL<<r)|(1ULL<<j)|(1ULL<<k),adj[j]);
 assert(cycles.size()==126000);
 // A_i,x -> A_i,x-i; B_j,y -> B_j+1,y.
 int move[50];for(int v=0;v<50;v++){
  int side=v/25,i=(v%25)/5,x=v%5;
  move[v]=side?25+5*((i+1)%5)+x:5*i+(x+5-i)%5;
 }
 for(int i=0;i<50;i++)for(int j=0;j<50;j++) assert((adj[i]>>j&1)==(adj[move[i]]>>move[j]&1));
 auto apply=[&](U m){U y=0;while(m){int i=__builtin_ctzll(m);m&=m-1;y|=1ULL<<move[i];}return y;};
 std::map<std::pair<int,int>,int> counts;int free=0;U first=0;
 for(U C:cycles){
  U D=apply(C),E=apply(D);int a=__builtin_popcountll(C&D),b=__builtin_popcountll(C&E);
  counts[{a==0,b==0}]++;
  if(a==0||b==0){free++;if(!first)first=C;}
 }
 const std::map<std::pair<int,int>,int> want={{{0,0},125900},{{0,1},50},{{1,0},50}};
 assert(counts==want && free==100);
 assert(first==1531568070816ULL);
 std::cout<<"126000 induced 13-cycles: 125900 no disjoint shifts, 50 each of the two allowed shift distances."<<std::endl;
}
