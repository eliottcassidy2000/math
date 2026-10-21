/* circ_counts.c -- exact Hamiltonian path / cycle counts of circulant tournaments.
 * mac-mini-2026-10-21-fanopsp.
 * For odd n, a circulant tournament Cay(Z_n, S) has S containing exactly one of
 * k, n-k for k = 1..m, m = (n-1)/2.  Classes are taken up to multipliers a in Z_n^*
 * (aS ~ S; a = -1 gives the converse, which has the same counts).
 * For each class we print
 *   S, mult (= #{a unit : aS = S}), H = n * (paths from 0), C (directed Ham cycles),
 *   rho0 (paths fixed by the anti-automorphism x -> -x composed with reversal),
 *   A = #{s in S : gcd(s,n) = 1} (rotation-invariant "arithmetic" Ham cycles).
 * Usage: ./circ_counts n [paley]   (paley: only the quadratic-residue set, n prime = 3 mod 4)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned __int128 u128;
static int n, m;
static int gcd(int a, int b){ while(b){int t=a%b;a=b;b=t;} return a; }
static void p128(u128 x){ char buf[64]; int i=63; buf[i]=0; if(x==0){printf("0");return;}
  while(x){ buf[--i]='0'+(int)(x%10); x/=10; } printf("%s",buf+i); }
static int setmask_of(const int *inS){ /* bit k-1 set iff n-k in S */
  int b=0; for(int k=1;k<=m;k++) if(inS[n-k]) b|=1<<(k-1); return b; }
int main(int argc,char**argv){
  n=atoi(argv[1]); m=(n-1)/2; int paley = (argc>2);
  int N1=n-1; /* vertices 1..n-1 indexed 0..n-2 */
  size_t full=((size_t)1<<N1);
  uint64_t *f = malloc(full*(size_t)N1*sizeof(uint64_t));
  if(!f){fprintf(stderr,"oom\n");return 1;}
  int units[64],nu=0; for(int a=1;a<n;a++) if(gcd(a,n)==1) units[nu++]=a;
  for(int b=0;b<(1<<m);b++){
    int inS[64]; memset(inS,0,sizeof inS);
    for(int k=1;k<=m;k++){ if(b>>(k-1)&1) inS[n-k]=1; else inS[k]=1; }
    if(paley){ /* keep only the QR set */
      int qr[64]; memset(qr,0,sizeof qr); for(int x=1;x<n;x++) qr[(x*x)%n]=1;
      int ok=1; for(int x=1;x<n;x++) if(inS[x]!=qr[x]) ok=0; if(!ok) continue; }
    else { /* canonical under multipliers */
      int canon=1, mult=0;
      for(int i=0;i<nu;i++){ int a=units[i]; int t[64]; memset(t,0,sizeof t);
        for(int x=1;x<n;x++) if(inS[x]) t[(a*x)%n]=1;
        int bb=setmask_of(t); if(bb<b){canon=0;break;} }
      if(!canon) continue; (void)mult; }
    int mult=0;
    for(int i=0;i<nu;i++){ int a=units[i]; int same=1;
      for(int x=1;x<n;x++) if(inS[x] && !inS[(a*x)%n]) {same=0;break;} if(same) mult++; }
    /* adjacency among 1..n-1 (index v-1) */
    uint32_t out[64]; for(int v=1;v<n;v++){ out[v-1]=0; for(int w=1;w<n;w++) if(w!=v && inS[((w-v)%n+n)%n]) out[v-1]|=1u<<(w-1); }
    memset(f,0,full*(size_t)N1*sizeof(uint64_t));
    for(int v=1;v<n;v++) if(inS[v]) f[((size_t)1<<(v-1))*N1+(v-1)]=1;
    for(size_t mask=1;mask<full;mask++){
      uint64_t *row=f+mask*N1;
      for(int v=0;v<N1;v++){ uint64_t c=row[v]; if(!c) continue;
        uint32_t cand=out[v] & ~(uint32_t)mask;
        while(cand){ int w=__builtin_ctz(cand); cand&=cand-1;
          f[(mask|((size_t)1<<w))*N1+w]+=c; } } }
    u128 P0=0, C=0; uint64_t *last=f+(full-1)*N1;
    for(int v=1;v<n;v++){ P0+=last[v-1]; if(inS[(n-v)%n]) C+=last[v-1]; } /* v->0 iff -v in S */
    /* rho0: paths 0,u1..um choosing one of each pair {x,-x}; DP over pair-sets */
    size_t pf=(size_t)1<<m; uint64_t *g=calloc(pf*(size_t)n,sizeof(uint64_t));
    for(int u=1;u<n;u++) if(inS[u]){ int pr=(u<=m)?u:n-u; g[((size_t)1<<(pr-1))*n+u]=1; }
    for(size_t pm=1;pm<pf;pm++) for(int v=1;v<n;v++){ uint64_t c=g[pm*n+v]; if(!c) continue;
      for(int w=1;w<n;w++){ int pr=(w<=m)?w:n-w; if(pm>>(pr-1)&1) continue; if(!inS[((w-v)%n+n)%n]) continue;
        g[(pm|((size_t)1<<(pr-1)))*n+w]+=c; } }
    u128 rho=0; for(int v=1;v<n;v++) rho+=g[(pf-1)*n+v]; free(g);
    int A=0; for(int s=1;s<n;s++) if(inS[s] && gcd(s,n)==1) A++;
    printf("n=%d S={",n); int first=1; for(int s=1;s<n;s++) if(inS[s]){ printf(first?"%d":",%d",s); first=0; }
    printf("} mult=%d H=",mult); p128(P0*(u128)n); printf(" H/n="); p128(P0); printf(" C="); p128(C);
    printf(" rho0="); p128(rho); printf(" A=%d\n",A); fflush(stdout);
  }
  free(f); return 0; }
