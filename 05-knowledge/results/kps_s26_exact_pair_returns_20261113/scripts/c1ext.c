/* c1ext.c: extend uniquely Hamiltonian (C=1) strong tournaments of order n to order n+1.
   input lines: n r_0 r_1 ... r_{n-1}   (r_i = out-neighbour bitmask of vertex i)
   output lines (order n+1, canonical, deduplicated): n+1 r_0 ... r_n
   Completeness: every strong tournament of order >= 4 has a vertex v with T-v strong (Moon);
   by Session197 Lemma 1, T-v is then C=1, so every C=1 strong tournament of order n+1 >= 5 arises.
   Canonical form: relabel along the unique Hamiltonian cycle, minimise the pair-bit word over rotations. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
#define MAXN 16
static u64 cnt[(1<<MAXN)];      /* per-start DP, flattened [S*n+x] allocated dynamically */
static int n;
static unsigned rowv[MAXN];
/* hash set of canonical codes (pair words up to 120 bits -> two u64) */
typedef struct { u64 a,b; int used; } ent;
static ent *tab; static size_t tabsz, tabcnt;
static unsigned *store; static size_t nstore, capstore; /* rows of outputs */
static int insert(u64 a,u64 b){
  size_t h=(a*0x9E3779B97F4A7C15ULL ^ b*0xC2B2AE3D27D4EB4FULL)%tabsz;
  while(tab[h].used){ if(tab[h].a==a&&tab[h].b==b) return 0; h=(h+1)%tabsz; }
  tab[h].used=1; tab[h].a=a; tab[h].b=b; tabcnt++; return 1;
}
/* unique Hamiltonian cycle of tournament out[] of order m; returns 1 and fills cyc[] if exactly one */
static u64 *st; static size_t stsz;
static int hamcycle(const unsigned*out,int m,int*cyc){
  unsigned N=1u<<m, FULL=N-1;
  size_t need=(size_t)N*m; if(need>stsz){ free(st); st=malloc(need*sizeof(u64)); stsz=need; }
  memset(st,0,need*sizeof(u64));
  st[1u*m+0]=1;
  for(unsigned S=1;S<N;S++){ if(!(S&1)) continue;
    for(int v=0;v<m;v++){ u64 a=st[(size_t)S*m+v]; if(!a) continue;
      for(int w=0;w<m;w++) if(!(S>>w&1)&&(out[v]>>w&1)) st[(size_t)(S|1u<<w)*m+w]+=a; } }
  u64 C=0; int last=-1;
  for(int v=1;v<m;v++) if(out[v]&1u){ C+=st[(size_t)FULL*m+v]; if(st[(size_t)FULL*m+v]) last=v; }
  if(C!=1) return 0;
  /* backtrack */
  unsigned S=FULL; int v=last; int pos=m-1;
  while(1){ cyc[pos--]=v; if(v==0) break; unsigned P=S&~(1u<<v); int found=-1;
    for(int u=0;u<m;u++) if((P>>u&1)&&(out[u]>>v&1)&&st[(size_t)P*m+u]){ found=u; break; }
    S=P; v=found; if(v<0) return 0; }
  return 1;
}
static void canon(const unsigned*out,int m,const int*cyc,u64*ca,u64*cb,unsigned*rows){
  u64 bestA=~0ULL,bestB=~0ULL; int bests=0;
  for(int s=0;s<m;s++){
    /* new label of cyc[(s+i)%m] is i */
    u64 A=0,B=0; int k=0;
    for(int i=0;i<m;i++) for(int j=i+1;j<m;j++){
      int x=cyc[(s+i)%m], y=cyc[(s+j)%m];
      if(out[x]>>y&1){ if(k<64) A|=1ULL<<k; else B|=1ULL<<(k-64); }
      k++; }
    if(B<bestB||(B==bestB&&A<bestA)){ bestA=A; bestB=B; bests=s; }
  }
  *ca=bestA; *cb=bestB;
  for(int i=0;i<m;i++){ rows[i]=0; }
  for(int i=0;i<m;i++) for(int j=0;j<m;j++) if(i!=j){
    int x=cyc[(bests+i)%m], y=cyc[(bests+j)%m]; if(out[x]>>y&1) rows[i]|=1u<<j; }
}
int main(int argc,char**argv){
  char line[4096];
  tabsz=1<<22; tab=calloc(tabsz,sizeof(ent));
  capstore=1<<20; store=malloc(capstore*sizeof(unsigned));
  int outn=0; long nin=0;
  static u64 P[MAXN][MAXN];
  u64 *dp=NULL; size_t dpsz=0;
  while(fgets(line,sizeof line,stdin)){
    char*p=line; int m=strtol(p,&p,10); if(m<=0) continue;
    for(int i=0;i<m;i++) rowv[i]=strtoul(p,&p,10);
    n=m; nin++; outn=m+1;
    unsigned N=1u<<m;
    size_t need=(size_t)N*m; if(need>dpsz){ free(dp); dp=malloc(need*sizeof(u64)); dpsz=need; }
    /* P[y][x] = Hamiltonian paths of S from y to x */
    for(int y=0;y<m;y++){
      memset(dp,0,need*sizeof(u64)); dp[(size_t)(1u<<y)*m+y]=1;
      for(unsigned S=1;S<N;S++){ if(!(S>>y&1)) continue;
        for(int v=0;v<m;v++){ u64 a=dp[(size_t)S*m+v]; if(!a) continue;
          for(int w=0;w<m;w++) if(!(S>>w&1)&&(rowv[v]>>w&1)) dp[(size_t)(S|1u<<w)*m+w]+=a; } }
      for(int x=0;x<m;x++) P[y][x]=(x==y)?0:dp[(size_t)(N-1)*m+x];
    }
    for(unsigned O=1;O<N-1;O++){ /* new vertex m beats O, loses to the rest */
      u64 C=0;
      for(int x=0;x<m && C<2;x++) if(!(O>>x&1)) for(int y=0;y<m;y++) if(O>>y&1){ C+=P[y][x]; if(C>=2) break; }
      if(C!=1) continue;
      unsigned out[MAXN]; for(int i=0;i<m;i++) out[i]=rowv[i] | ((O>>i&1)?0u:(1u<<m));
      out[m]=O;
      int cyc[MAXN]; if(!hamcycle(out,m+1,cyc)){ fprintf(stderr,"cycle mismatch\n"); continue; }
      u64 a,b; unsigned rows[MAXN]; canon(out,m+1,cyc,&a,&b,rows);
      if(insert(a,b)){
        if(nstore+outn>capstore){ capstore*=2; store=realloc(store,capstore*sizeof(unsigned)); }
        for(int i=0;i<outn;i++) store[nstore++]=rows[i];
      }
    }
  }
  size_t cls=nstore/outn;
  for(size_t c=0;c<cls;c++){ printf("%d",outn); for(int i=0;i<outn;i++) printf(" %u",store[c*outn+i]); putchar(10); }
  fprintf(stderr,"input %ld classes of order %d -> %zu C=1 classes of order %d\n",nin,n,cls,outn);
  return 0;
}
