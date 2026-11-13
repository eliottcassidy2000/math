/* rprof.c: full path-cover profile from adjacency rows.
   input lines: n r_0 ... r_{n-1};  output: n H F_1..F_n | r_0 ... r_{n-1}
   F_j = ordered covers by j nonempty vertex-disjoint paths (Q(t)=sum_j F_j binom(t,j)). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
typedef unsigned __int128 u128;
#define MAXN 16
static u64 Hs[1<<MAXN]; static u64 *endp; static u128 g[2][1<<MAXN];
static void pr128(u128 x){ char b[64]; int k=0; if(!x){ putchar('0'); return; } while(x){ b[k++]='0'+(int)(x%10); x/=10; } while(k) putchar(b[--k]); }
int main(void){
  char line[4096]; endp=malloc(sizeof(u64)*(1u<<MAXN)*MAXN);
  while(fgets(line,sizeof line,stdin)){
    char*p=line; int n=strtol(p,&p,10); if(n<=0) continue; unsigned out[MAXN];
    for(int i=0;i<n;i++) out[i]=strtoul(p,&p,10);
    unsigned N=1u<<n;
    Hs[0]=0;
    for(unsigned S=1;S<N;S++){ u64 tot=0;
      for(int v=0;v<n;v++) if(S>>v&1){ unsigned P=S&~(1u<<v); u64 a=0;
        if(!P) a=1; else for(int w=0;w<n;w++) if((P>>w&1)&&(out[w]>>v&1)) a+=endp[(size_t)P*n+w];
        endp[(size_t)S*n+v]=a; tot+=a; }
      Hs[S]=tot; }
    printf("%d %llu",n,(unsigned long long)Hs[N-1]);
    for(unsigned S=0;S<N;S++) g[0][S]=Hs[S];
    putchar(' '); pr128(g[0][N-1]);
    int cur=0;
    for(int j=2;j<=n;j++){ int nx=cur^1;
      for(unsigned S=0;S<N;S++){ u128 s=0; unsigned low=S&(~S+1u);
        /* ordered covers: last path B, earlier covers of A; require nothing special */
        for(unsigned B=S; B; B=(B-1)&S){ unsigned A=S&~B; if(!A) continue; s+=g[cur][A]*(u128)Hs[B]; }
        (void)low; g[nx][S]=s; }
      cur=nx; putchar(' '); pr128(g[cur][N-1]); }
    printf(" |"); for(int i=0;i<n;i++) printf(" %u",out[i]); putchar(10);
  }
  return 0;
}
