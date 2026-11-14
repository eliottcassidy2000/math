// Misses of root-change (O/Z/J) words with at most CMAX root changes (c only after b/c or at start), classified:
// prime / 7*prime / other.  Letters a=O, b=Z, c=J.  Per dyadic window counts.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64; typedef unsigned long long u64;
static i64 X; static u64 *bits; static int CMAX; static unsigned char *comp;
#define SETB(n) (bits[(n)>>7] |= 1ULL<<(((n)>>1)&63))
#define GETB(n) ((bits[(n)>>7]>>(((n)>>1)&63))&1ULL)
static void dfs(i64 H,i64 E,i64 L,i64 h,int last,int nc){
  SETB(H); i64 D=L-E; if(D>X-H) return;
  if(2*L<=X-H) dfs(H,0,2*L,h,1,nc);
  if(!(D==0 && L==h && E==h)) dfs(H+D,h,2*L-h,h,2,nc);
  if(nc<CMAX && last!=1) dfs(H+D,H-E,H+L,H,3,nc+1);
}
int main(int argc,char**argv){
  X=atoll(argv[1]); CMAX=atoi(argv[2]);
  bits=calloc((X>>7)+2,sizeof(u64)); comp=calloc(X/2+2,1);   // comp[n>>1]=1 if odd n composite
  for(i64 p=3;p*p<=X;p+=2) if(!comp[p>>1]) for(i64 q=p*p;q<=X;q+=2*p) comp[q>>1]=1;
  dfs(1,1,1,1,0,0);
  printf("CMAX=%d X=%lld: misses per window [2^k,2^(k+1)): all / prime / 7*prime\n",CMAX,X);
  for(int k=3;(1LL<<k)<X;k++){ i64 lo=1LL<<k,hi=1LL<<(k+1); if(hi>X) hi=X; i64 m=0,mp=0,m7=0; i64 lastp=0;
    for(i64 n=lo|1;n<hi;n+=2) if(!GETB(n)){ m++; if(!comp[n>>1] && n>1){mp++; lastp=n;} else if(n%7==0 && !comp[(n/7)>>1] && (n/7)>1) m7++; }
    printf("  k=%2d: %8lld  primes %7lld  7p %6lld  (last prime miss %lld)\n",k,m,mp,m7,lastp); }
  return 0; }
