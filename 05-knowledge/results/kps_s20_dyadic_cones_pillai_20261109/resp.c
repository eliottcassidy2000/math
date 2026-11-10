/* resp.c: for each input tournament (line: n row_0 ... row_{n-1}, row_v = out-neighbour bitmask of v)
   print: H Gamma R_0 ... R_{n-1}   where R_z = sum_{S ni z} #(Ham paths of T[S] starting at z) * H(T[V\S])
   and Gamma = sum_S H(T[S]) H(T[V\S]).  All counts in unsigned 128-bit printed as decimal (fits for n<=20). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned __int128 u128;
static void pr(u128 x){char b[50];int i=49;b[i]=0;if(!x){putchar('0');return;}while(x){b[--i]='0'+(int)(x%10);x/=10;}fputs(b+i,stdout);}
int main(){
  int n; static unsigned out[32];
  while(scanf("%d",&n)==1){
    for(int i=0;i<n;i++) scanf("%u",&out[i]);
    unsigned N=1u<<n;
    u128 *Hs=calloc(N,sizeof(u128));
    u128 *end=calloc((size_t)N*n,sizeof(u128));   /* end[S*n+v]: paths covering S ending at v */
    u128 *st=calloc((size_t)N*n,sizeof(u128));    /* st[S*n+v]: paths covering S starting at v */
    Hs[0]=1;
    for(unsigned S=1;S<N;S++){
      u128 tot=0;
      for(int v=0;v<n;v++) if(S>>v&1){
        unsigned P=S&~(1u<<v);
        u128 a=0,b=0;
        if(!P){a=1;b=1;}
        else for(int w=0;w<n;w++) if(P>>w&1){
          if(out[w]>>v&1) a+=end[(size_t)P*n+w];
          if(out[v]>>w&1) b+=st[(size_t)P*n+w];
        }
        end[(size_t)S*n+v]=a; st[(size_t)S*n+v]=b; tot+=a;
      }
      Hs[S]=tot;
    }
    u128 G=0; for(unsigned S=0;S<N;S++) G+=Hs[S]*Hs[(N-1)&~S];
    pr(Hs[N-1]); putchar(' '); pr(G);
    for(int z=0;z<n;z++){
      u128 R=0;
      for(unsigned S=0;S<N;S++) if(S>>z&1) R+=st[(size_t)S*n+z]*Hs[(N-1)&~S];
      putchar(' '); pr(R);
    }
    printf(" |");
    for(int z=0;z<n;z++){ putchar(' '); pr(end[(size_t)(N-1)*n+z]); }
    putchar('\n'); fflush(stdout);
    free(Hs);free(end);free(st);
  }
  return 0;
}
