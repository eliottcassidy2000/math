#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <vector>
int main(int argc,char**argv){
 if(argc!=3)throw std::runtime_error("usage: orbits n file");
 int n=std::atoi(argv[1]); if(n<3||n>7)throw std::runtime_error("supported orders 3..7");
 int pos[7][7]={},m=0;std::vector<std::pair<int,int>> E;
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){pos[i][j]=m++;E.push_back({i,j});}
 int chunks=(m+7)/8;uint32_t size=1u<<m,unset=std::numeric_limits<uint32_t>::max();
 using Chunk=std::array<uint32_t,256>;
 std::vector<std::vector<Chunk>> tables(n-1,std::vector<Chunk>(chunks));
 std::vector<uint32_t> flip(n-1,0);
 for(int a=0;a<n-1;a++){
  int p[7];for(int i=0;i<n;i++)p[i]=i;std::swap(p[a],p[a+1]);
  for(int k=0;k<m;k++){
   int i=p[E[k].first],j=p[E[k].second];if(i>j){std::swap(i,j);flip[a]|=1u<<k;}
   int from=pos[i][j];for(int v=0;v<256;v++)if(v&(1u<<(from%8)))tables[a][from/8][v]|=1u<<k;
  }
 }
 std::vector<uint32_t> owner(size,unset),queue;uint32_t classes=0;
 for(uint32_t t=0;t<size;t++)if(owner[t]==unset){
  classes++;queue.clear();queue.push_back(t);owner[t]=t;
  for(size_t at=0;at<queue.size();at++)for(int a=0;a<n-1;a++){
   uint32_t u=queue[at],v=flip[a];for(int b=0;b<chunks;b++)v^=tables[a][b][(u>>(8*b))&255];
   if(owner[v]==unset){owner[v]=t;queue.push_back(v);}else if(owner[v]!=t)throw std::runtime_error("orbit collision");
  }
 }
 FILE*f=std::fopen(argv[2],"wb");if(!f)throw std::runtime_error("cannot open output");
 if(std::fwrite(owner.data(),sizeof(uint32_t),size,f)!=size)throw std::runtime_error("short write");std::fclose(f);
 std::printf("order=%d labelled=%u ordinary_classes=%u\n",n,size,classes);
}
