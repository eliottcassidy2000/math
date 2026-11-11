/* pairstate.c - literal state (r, M, h, E) of a pair (B, O) by subset dynamic programming (orders up to 20).
 * input (stdin), repeated:  n  out_0 ... out_{n-1}  Omask      (out_i = bit mask of the out-neighbours of vertex i)
 * output per pair:          n r M h E
 *   h = H(B); W = sum over ordered pairs e != f in O of #(P1 ending at e, P2 starting at f, both nonempty, covering B);
 *   r = h - W; M = sum_{o in O} #(P1 possibly empty, P2 starting at o); E = number of Hamiltonian paths of B ending in O. */
#include <stdio.h>
#include <stdlib.h>
typedef unsigned long long u64; typedef __int128 i128;
static void pr(i128 x){ if(x<0){putchar('-');x=-x;} char b[64]; int k=0; if(x==0)b[k++]='0'; while(x>0){b[k++]='0'+(int)(x%10);x/=10;} while(k--)putchar(b[k]); }
int main(void){ int n;
  while(scanf("%d",&n)==1){ unsigned out[24],inn[24]; for(int i=0;i<n;i++){ scanf("%u",&out[i]); inn[i]=0; }
    unsigned O; scanf("%u",&O); for(int i=0;i<n;i++) for(int j=0;j<n;j++) if(out[i]>>j&1) inn[j]|=1u<<i;
    size_t full=((size_t)1<<n)-1; u64 *en=calloc((full+1)*n,sizeof(u64)), *st=calloc((full+1)*n,sizeof(u64)), *hs=calloc(full+1,sizeof(u64));
    if(!en||!st||!hs){ fprintf(stderr,"memory\n"); return 1; }
    for(int v=0;v<n;v++){ en[((size_t)1<<v)*n+v]=1; st[((size_t)1<<v)*n+v]=1; }
    for(size_t s=1;s<=full;s++) for(int v=0;v<n;v++){ u64 c=en[s*n+v]; if(c){ unsigned nx=out[v]&~(unsigned)s; while(nx){int u=__builtin_ctz(nx);nx&=nx-1;en[(s|((size_t)1<<u))*n+u]+=c;} }
      c=st[s*n+v]; if(c){ unsigned nx=inn[v]&~(unsigned)s; while(nx){int u=__builtin_ctz(nx);nx&=nx-1;st[(s|((size_t)1<<u))*n+u]+=c;} } }
    hs[0]=1; for(size_t s=1;s<=full;s++){ u64 t=0; for(int v=0;v<n;v++) t+=en[s*n+v]; hs[s]=t; }
    i128 h=hs[full], W=0, M=0, E=0;
    for(int o=0;o<n;o++) if(O>>o&1){ E+=en[full*n+o];
      for(size_t X=1;X<=full;X++) if((X>>o&1)&&st[X*n+o]) M+=(i128)st[X*n+o]*(i128)hs[full^X]; }
    for(size_t X=1;X<full;X++){ size_t Y=full^X; i128 a=0,b=0; for(int e=0;e<n;e++) if((O>>e&1)){ if(X>>e&1) a+=en[X*n+e]; else b+=st[Y*n+e]; } W+=a*b; }
    printf("%d ",n); pr(h-W); putchar(' '); pr(M); putchar(' '); pr(h); putchar(' '); pr(E); putchar('\n'); fflush(stdout);
    free(en); free(st); free(hs); }
  return 0; }
