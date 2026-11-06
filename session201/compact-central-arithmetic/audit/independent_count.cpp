#include <cstdint>
#include <iostream>
#include <vector>
#include <string>
#include "independent_data.h"
static std::string decimal(__uint128_t x){if(!x)return "0";std::string s;while(x){s=char('0'+x%10)+s;x/=10;}return s;}
struct Counts{__uint128_t paths=0,cycles=0,nodes=0;};
void visit(int first,int last,uint32_t left,const uint32_t*rows,Counts&out){
 ++out.nodes;
 if(!left){++out.paths;if(first==0&&(rows[last]&1))++out.cycles;return;}
 uint32_t candidates=rows[last]&left;
 while(candidates){uint32_t bit=candidates&-candidates;candidates-=bit;visit(first,__builtin_ctz(bit),left^bit,rows,out);}
}
int main(){
 const int N=2134117;std::vector<int>g(N+1);g[1]=1;
 for(int i=2;i<=N;++i){if(!(i%2))g[i]=g[i/2];else{int q=(i+1)/2;g[i]=g[q]+q-q%2;}}
 std::vector<uint8_t>present(N+1,0);uint64_t tests=0;
 for(int f=0;f<NF;++f){
  // Both endpoint orders are deliberately included on the compact menu.
  for(int u=2;1LL*u*BF[f]<=2LL*N;++u){
   int64_t tau[]={g[u],u%2,u,1},row[4]={};
   for(int i=0;i<4;++i)for(int j=0;j<4;++j)row[j]+=tau[i]*AF[f][i][j];
   for(int v=2;1LL*u*v*BF[f]<=4LL*N;++v){
    int64_t n=row[0]*g[v]+row[1]*(v%2)+row[2]*v+row[3];++tests;
    if(n<3||!(n%2)){std::cerr<<"nonphysical output\n";return 2;}
    if(n<=N)present[n]=1;
   }
  }
 }
 std::cout<<"{\"cap\":"<<N<<",\"compact_forms\":12,\"directed_pair_tests\":"<<tests<<",\"gaps5\":[";
 bool first=true;for(int n=3;n<=N;n+=2)if(!present[n]){if(!first)std::cout<<",";first=false;std::cout<<n;}
 std::cout<<"],\"literal_graphs\":[";
 for(int k=0;k<NG;++k){
  int n=ORDERS[k];if(n>24)return 3;uint32_t full=(1u<<n)-1;
  for(int i=0;i<n;++i){if((ROWS[k][i]&~full)||(ROWS[k][i]&(1u<<i)))return 4;for(int j=0;j<i;++j)if(((ROWS[k][i]>>j)&1)+((ROWS[k][j]>>i)&1)!=1)return 5;}
  Counts c;for(int s=0;s<n;++s)visit(s,s,full^(1u<<s),ROWS[k],c);
  if(c.paths!=EXPECTED[k]||c.cycles!=1){std::cerr<<"literal count mismatch "<<EXPECTED[k]<<"\n";return 6;}
  if(k)std::cout<<",";
  std::cout<<"{\"target\":"<<EXPECTED[k]<<",\"vertices\":"<<n<<",\"paths\":"<<decimal(c.paths)<<",\"cycles\":"<<decimal(c.cycles)<<",\"partial_paths\":"<<decimal(c.nodes)<<"}";
 }
 std::cout<<"],\"universal_occupancy\":\"OPEN\"}\n";
}
