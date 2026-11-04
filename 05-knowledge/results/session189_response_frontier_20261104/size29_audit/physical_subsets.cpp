#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;
// Ordinary physical-subset endpoint recurrence, quotiented only by rotations
// of actual disjoint directed triangles. All four special vertices can instead
// remain separate singleton groups, so no special-child cover formula is used.
int orbit(int cycles,const vector<int>& out){
 int n=out.size(), singles=n-cycles, states=1<<(2*cycles+singles), stride=n+cycles;
 vector<uint8_t> dp((size_t)states*stride);vector<int> incoming(n),step(n);
 for(int i=0;i<n;i++){step[i]=1<<(i<cycles?2*i:cycles+i);for(int j=0;j<n;j++)if(out[j]>>i&1)incoming[i]|=1<<j;}
 for(int s=1;s<states;s++){
  auto *dest=&dp[(size_t)s*stride];
  for(int i=0;i<n;i++){
   int occupied=(s/step[i])&(i<cycles?3:1);if(!occupied)continue;
   int prev=s-step[i];const auto *src=&dp[(size_t)prev*stride];int external=prev==0;
   for(unsigned m=incoming[i];m;m&=m-1)external+=src[__builtin_ctz(m)];
   if(i>=cycles||occupied==1)dest[i]=external&7;
   else if(occupied==2){dest[i]=(2*external+src[i])&7;dest[n+i]=(external+src[i])&7;}
   else dest[i]=(3*(external+src[n+i]))&7;
  }
 }
 int answer=0;for(int i=0;i<n;i++)answer+=dp[(size_t)(states-1)*stride+i];return answer&7;
}
uint64_t literal(const vector<int>&out){
 int n=out.size(),states=1<<n;vector<uint64_t>dp((size_t)states*n);vector<int>incoming(n);
 for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(out[j]>>i&1)incoming[i]|=1<<j;
 for(int s=1;s<states;s++)for(unsigned todo=s;todo;todo&=todo-1){int i=__builtin_ctz(todo),prev=s^(1<<i);uint64_t x=prev==0;
  for(unsigned m=prev&incoming[i];m;m&=m-1) { x+=dp[(size_t)prev*n+__builtin_ctz(m)]; } dp[(size_t)s*n+i]=x;}
 uint64_t answer=0;for(int i=0;i<n;i++)answer+=dp[(size_t)(states-1)*n+i];return answer;
}
int main(){int cases;cin>>cases;for(int z=0;z<cases;z++){int id,cycles,n,mode;cin>>id>>cycles>>n>>mode;vector<int>r(n);for(int &x:r)cin>>x;cout<<id<<' '<<(mode?literal(r):orbit(cycles,r))<<'\n';cout.flush();if(z%100==0)cerr<<z+1<<'/'<<cases<<'\n';}}
