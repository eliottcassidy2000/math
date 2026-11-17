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

/* Wave 58 P115 test: one-parent tournaments (vertex j beats only par[j] < j among earlier vertices) up to NMAX vertices;
   counts h with a vertex t of ratio B(T,t) = 3h.  Then 7h = h + 2*3h is a k = 1 near-join over that base. */
static unsigned char has[1<<22]; static int NMAXg; static int par[MAXN];
static void rec(int j){ int n=j; /* evaluate current tree on n vertices */
  if(n>=2){ unsigned out[MAXN]={0}; for(int b=1;b<n;b++) for(int a=0;a<b;a++){ if(par[b]==a) out[b]|=1u<<a; else out[a]|=1u<<b; }
    ull H, r[MAXN], s[MAXN]; regs(n,out,&H,r,s); if(H<(1u<<22)) for(int v=0;v<n;v++) if(r[v]+s[v]==3*H){ has[H]=1; break; } }
  if(j==NMAXg) return; for(int p=0;p<j;p++){ par[j]=p; rec(j+1); } }
int main(int argc,char**argv){ NMAXg=atoi(argv[1]); rec(1); int c=0,cp=0; int firstmiss=0;
  for(int h=1;h<(1<<22);h++) if(has[h]){ c++; int pr=h>1; for(int d=2;d*d<=h;d++) if(h%d==0){pr=0;break;} if(pr) cp++; }
  printf("one-parent tournaments to %d vertices: %d counts with a ratio-3 root, %d of them prime\n",NMAXg,c,cp);
  printf("primes with a ratio-3 root:"); int k=0; for(int h=1;h<(1<<22)&&k<80;h++) if(has[h]){ int pr=h>1; for(int d=2;d*d<=h;d++) if(h%d==0){pr=0;break;} if(pr){ printf(" %d",h); k++; } } printf("\n");
  printf("primes 5..1000 WITHOUT a ratio-3 root:"); for(int h=5;h<=1000;h+=2){ int pr=1; for(int d=2;d*d<=h;d++) if(h%d==0){pr=0;break;} if(pr && h!=7 && !has[h]) printf(" %d",h); } printf("\n"); return 0; }
