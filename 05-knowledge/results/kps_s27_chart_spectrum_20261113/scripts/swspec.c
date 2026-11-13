/* read tournaments (n rows) and print every H over all vertex switchings (switch set S containing no vertex 0) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static u64 endp[(1<<12)*12];
static char seen[1<<24];
int main(int argc,char**argv){
  long bound=atol(argv[1]); char line[4096];
  while(fgets(line,sizeof line,stdin)){
    char*p=line; int n=strtol(p,&p,10); unsigned r[16];
    for(int i=0;i<n;i++) r[i]=strtoul(p,&p,10);
    unsigned N=1u<<n;
    for(unsigned S=0;S<(1u<<(n-1));S++){ unsigned sw=S<<1; unsigned out[16];
      for(int i=0;i<n;i++){ unsigned o=r[i]; for(int j=0;j<n;j++) if(j!=i && (((sw>>i)^(sw>>j))&1)) o^=1u<<j; out[i]=o; }
      u64 H=0;
      for(unsigned T=1;T<N;T++) for(int v=0;v<n;v++) if(T>>v&1){ unsigned P=T&~(1u<<v); u64 a=0;
          if(!P) a=1; else for(int w=0;w<n;w++) if((P>>w&1)&&(out[w]>>v&1)) a+=endp[(size_t)P*n+w];
          endp[(size_t)T*n+v]=a; if(T==N-1) H+=a; }
      if(H<(u64)bound) seen[H]=1; }
  }
  for(long h=1;h<bound;h+=2) if(!seen[h]) printf("%ld ",h); printf("\n");
  return 0;
}
