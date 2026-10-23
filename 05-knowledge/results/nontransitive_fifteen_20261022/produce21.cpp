#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <vector>
constexpr int N=21;
using Key=std::array<uint64_t,4>;
std::array<uint32_t,N> out;
std::array<int,N> path;
std::set<Key> keys;
std::vector<std::array<int,N>> words;
uint64_t nodes,completed;
int mode;
void need(bool x,const char*s){if(!x)throw std::runtime_error(s);}
int mul(int a,int b){
 if(!mode)return(a+b)%21;
 int x=a%7,y=a/7,u=b%7,v=b/7;
 return(x+(1<<y)*u)%7+7*((y+v)%3);
}
int inverse(int a){for(int b=0;b<N;++b)if(mul(a,b)==0&&mul(b,a)==0)return b;throw std::runtime_error("inverse");}
bool arc(int a,int b){return out[a]>>b&1;}
Key key(const std::array<int,N>&p){Key k{};int bit=0;
 for(int i=0;i<N;++i)for(int j=i+1;j<N;++j,++bit)if(arc(p[i],p[j]))k[bit/64]|=uint64_t(1)<<(bit%64);
 return k;
}
void record(){++completed;bool blue=true;
 for(int i=0;i<N;++i)for(int j=i+1;j<N;++j)if(arc(path[i],path[j])==arc(path[N-1-i],path[N-1-j]))blue=false;
 if(!blue&&keys.insert(key(path)).second)words.push_back(path);
}
void dfs(int a,uint32_t used,int depth){
 if(words.size()==7)return;
 if(++nodes>10000000)throw std::runtime_error("path discovery cap exhausted");
 if(depth==N){record();return;}
 uint32_t next=out[a]&~used;
 while(next&&words.size()<7){auto bit=next&-next;next-=bit;int b=std::countr_zero(bit);path[depth]=b;dfs(b,used|bit,depth+1);}
}
int main(int argc,char**argv){
 need(argc==2,"raw path");std::ofstream dst(argv[1],std::ios::binary);need(bool(dst),"output open");dst.write("VT21KEY1",8);
 uint64_t total=0,mx=0,mw=0;
 for(mode=0;mode<2;++mode){
  std::vector<int> reps;for(int a=1;a<N;++a)if(a<inverse(a))reps.push_back(a);need(reps.size()==10,"ten inverse pairs");
  for(int mask=0;mask<1024;++mask){
   std::set<int>S;for(int i=0;i<10;++i)S.insert(mask>>i&1?inverse(reps[i]):reps[i]);out.fill(0);
   for(int x=0;x<N;++x)for(int y=0;y<N;++y)if(x!=y&&S.count(mul(inverse(x),y)))out[x]|=1u<<y;
   keys.clear();words.clear();nodes=completed=0;path[0]=0;dfs(0,1,1);need(words.size()==7,"seven words");
   std::array<int,N> identity;for(int i=0;i<N;++i)identity[i]=i;Key src=key(identity);
   for(auto p:words){
    for(int x:p)dst.put(char('A'+x));Key target=src;
    for(int i=0;i<210;++i)target[i/64]^=uint64_t(1)<<(i%64);
    for(int j=0;j<N-1;++j){int a=p[j],b=p[j+1];if(a>b)std::swap(a,b);int k=a*(2*N-a-1)/2+b-a-1;target[k/64]^=uint64_t(1)<<(k%64);}
    for(int j=0;j<27;++j)dst.put(char(target[j/8]>>(8*(j%8))&255));
   }
   total+=nodes;mx=std::max(mx,nodes);mw=std::max(mw,completed);
  }
 }
 need(bool(dst),"output complete");std::cout<<"2048 frames,14336 words; nodes "<<total<<" maximum "<<mx<<" completed "<<mw<<'\n';
}
