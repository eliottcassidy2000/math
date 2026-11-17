/* Wave 61. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long ull;
#define MAXN 13
static ull Fs[1 << MAXN], dpE[(1 << MAXN) * MAXN], dpS[(1 << MAXN) * MAXN];
static void regs(int n, const unsigned *out, ull *H, ull *r, ull *s) {
    int N = 1 << n;
    memset(dpE, 0, sizeof(ull) * N * n); memset(dpS, 0, sizeof(ull) * N * n);
    unsigned in[MAXN] = {0};
    for (int v = 0; v < n; v++) for (int u = 0; u < n; u++) if (out[v] >> u & 1) in[u] |= 1u << v;
    for (int v = 0; v < n; v++) { dpE[(1u << v) * n + v] = 1; dpS[(1u << v) * n + v] = 1; }
    for (int S = 1; S < N; S++) for (int v = 0; v < n; v++) {
        ull c = dpE[S * n + v]; if (c) { unsigned t = out[v] & ~(unsigned)S; while (t) { int u = __builtin_ctz(t); t &= t - 1; dpE[(S | 1u << u) * n + u] += c; } }
        ull d = dpS[S * n + v]; if (d) { unsigned t = in[v] & ~(unsigned)S; while (t) { int u = __builtin_ctz(t); t &= t - 1; dpS[(S | 1u << u) * n + u] += d; } } }
    Fs[0] = 1;
    for (int S = 1; S < N; S++) { ull a = 0; for (int v = 0; v < n; v++) a += dpE[S * n + v]; Fs[S] = a; }
    *H = Fs[N - 1];
    for (int v = 0; v < n; v++) {
        ull B = 0;
        for (int S = 1; S < N; S++) if (S >> v & 1) B += dpS[S * n + v] * Fs[(N - 1) & ~S];
        ull he = dpE[(N - 1) * n + v];
        r[v] = (B - he) / 2; s[v] = (B + he) / 2;
    }
}

/* Two-reading cone: near-sink marker X (beats only y) blown up into I7 or P7 gives a Z7-symmetric tournament with orbit count
   o_I = 25 a1 + 234 a2, o_P = 27 a1 + 234 a2, a1 = beta + H - e_y, a2 = beta - g (verified).  Bases: one-parent trees to NMAX. */
static unsigned char *hit; static unsigned long long LIM; static int NMAXg; static int par[MAXN];
static void eval(int n, const unsigned*out){ ull H, r[MAXN], s[MAXN]; regs(n,out,&H,r,s);
  for(int y=0;y<n;y++){ ull beta=r[y]+s[y], ey=s[y]-r[y], g=Fs[((1u<<n)-1)&~(1u<<y)];
    ull a1=beta+H-ey, a2=beta-g; ull oI=25*a1+234*a2, oP=27*a1+234*a2; if(oI<=LIM) hit[oI]|=1; if(oP<=LIM) hit[oP]|=2; } }
static void rec(int j){ int n=j; if(n>=1){ unsigned out[MAXN]={0}; for(int b=1;b<n;b++) for(int a=0;a<b;a++){ if(par[b]==a) out[b]|=1u<<a; else out[a]|=1u<<b; } eval(n,out); }
  if(j==NMAXg) return; for(int p=0;p<j;p++){ par[j]=p; rec(j+1); } }
static int isprime(ull q){ if(q<2) return 0; if(q%2==0) return q==2; for(ull d=3;d*d<=q;d+=2) if(q%d==0) return 0; return 1; }
int main(int argc,char**argv){ NMAXg=atoi(argv[1]); LIM=strtoull(argv[2],0,10); hit=calloc(LIM+1,1); par[0]=-1; rec(1);
  ull lastmissI=0; long odd=0,oddhit=0,pr=0,prhit=0; ull firstodd=0;
  /* report: coverage of odd numbers and primes, and the largest odd gap */
  ull X=LIM; for(ull q=1;q<=X;q+=2){ odd++; if(hit[q]) oddhit++; else lastmissI=q; if(isprime(q)){ pr++; if(hit[q]) prhit++; } }
  printf("one-parent bases to %d vertices: odd orbit counts <= %llu hit %ld of %ld; primes hit %ld of %ld; largest odd miss %llu\n",NMAXg,LIM,oddhit,odd,prhit,pr,lastmissI);
  /* first window of complete coverage: smallest K such that all odd numbers in [K, L] are hit, for L = LIM/4 */
  ull L=LIM/4, K=L; while(K>=1 && hit[K]) K-=2; printf("all odd numbers in (%llu, %llu] are orbit counts\n",K,L);
  for(ull lo=16; lo<LIM; lo*=2){ ull hi=lo*2>LIM?LIM:lo*2; long t=0,h=0; for(ull q=lo+1;q<hi;q+=2){ t++; if(hit[q]) h++; } printf("  window [%llu,%llu): odd hit %ld of %ld (%.1f%%)\n",lo,hi,h,t,100.0*h/t); }
  printf("first misses above 10000:"); int c=0; for(ull q=10001;q<=LIM&&c<15;q+=2) if(!hit[q]){ printf(" %llu",q); c++; } printf("\n"); return 0; }
