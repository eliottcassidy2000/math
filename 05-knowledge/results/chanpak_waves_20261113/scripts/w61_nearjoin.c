/* Wave 61: HYP-9453 near-join family over one-parent bases. */
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

/* Wave 61: near-join family N = h + 2^k B(T0,t0) (HYP-9453), k >= 1, over every rooted one-parent base T0 with at most NMAX vertices.
   Report odd N <= LIM not reached. */
static unsigned char *hit; static unsigned long long LIM; static int NMAXg; static int par[MAXN];
static void eval(int n,const unsigned*out){ ull H,r[MAXN],s[MAXN]; regs(n,out,&H,r,s);
  for(int t=0;t<n;t++){ ull B=r[t]+s[t]; for(ull p2=2; H+p2*B<=LIM; p2*=2) hit[H+p2*B]=1; } }
static void rec(int j){ int n=j; if(n>=1){ unsigned out[MAXN]={0}; for(int b=1;b<n;b++) for(int a=0;a<b;a++){ if(par[b]==a) out[b]|=1u<<a; else out[a]|=1u<<b; } eval(n,out); }
  if(j==NMAXg) return; for(int p=0;p<j;p++){ par[j]=p; rec(j+1); } }
int main(int argc,char**argv){ NMAXg=atoi(argv[1]); LIM=strtoull(argv[2],0,10); hit=calloc(LIM+1,1); rec(1);
  long miss=0; printf("near-joins over one-parent bases to %d vertices: odd N <= %llu not reached:",NMAXg,LIM);
  for(ull q=1;q<=LIM;q+=2) if(!hit[q]){ miss++; if(miss<=40) printf(" %llu",q); } printf("  (total %ld)\n",miss); return 0; }
