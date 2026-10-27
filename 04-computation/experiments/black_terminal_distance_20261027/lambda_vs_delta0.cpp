#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <map>
int main(int argc,char**argv){int n=atoi(argv[1]);char line[256];long long tot=0,lt=0;std::map<std::pair<int,int>,long long>h;
 while(fgets(line,sizeof line,stdin)){if(line[0]!='0'&&line[0]!='1')continue;uint32_t out[16]={0};int k=0;
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){if(line[k]=='1')out[i]|=1u<<j;else out[j]|=1u<<i;k++;}
  int full=(1<<n)-1,best=99;for(int X=1;X<full;X++){int f=0;for(int x=0;x<n;x++)if(X>>x&1)f+=__builtin_popcount(out[x]&~X&full);if(f<best)best=f;}
  int d0=99;for(int v=0;v<n;v++){int s=__builtin_popcount(out[v]);d0=std::min(d0,std::min(s,n-1-s));}
  tot++;if(best<d0){lt++;h[{best,d0}]++;}}
 printf("n=%d strong classes %lld, with lambda < delta0: %lld",n,tot,lt);for(auto&p:h)printf(" (lambda %d, delta0 %d):%lld",p.first.first,p.first.second,p.second);printf("\n");}
