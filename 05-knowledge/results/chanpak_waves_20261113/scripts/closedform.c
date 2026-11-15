// Union of two-tape carriers: spine [guard 1][x-tape][C][y-tape][guard 0] + one apex, over all strong marked cores C of
// orders 3..MAXO (gentourng classes, all masks). Counts via the exact cut transfer. Reports coverage of odd N <= X
// and, optionally, multiplicities in a window.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64;
#define MAXM 10
typedef struct { i64 a, q, s, e, t, f, n; } blk;
typedef struct { i64 h, N, S, d; } st;
static st stepb(st x, blk b) { st y; y.h = b.a * x.h; y.N = b.a * x.N + (b.n - b.a) * x.h + b.t * x.S + b.f * x.d; y.S = b.q * x.S + x.h * b.s; y.d = x.h * b.e; return y; }
static st step1(st x, int mu) { st y; y.h = x.h; y.N = x.N + (1 - mu) * (x.S + x.d); y.S = 2 * x.S + mu * x.h; y.d = mu * x.h; return y; }
static i64 Hsub[1 << MAXM], dpS[(1 << MAXM) * MAXM], dpR[(1 << MAXM) * MAXM];
static i64 hamcount(int n, const int *out) { int N = 1 << n; i64 *dp = calloc((size_t)N * n, sizeof(i64));
  for (int v = 0; v < n; v++) dp[(size_t)(1 << v) * n + v] = 1;
  for (int S = 1; S < N; S++) for (int v = 0; v < n; v++) { i64 x = dp[(size_t)S * n + v]; if (!x) continue; int nb = out[v] & ~S; while (nb) { int u = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | 1 << u) * n + u] += x; } }
  i64 s = 0; for (int v = 0; v < n; v++) s += dp[(size_t)(N - 1) * n + v]; free(dp); return s; }
static void vdata(int m, const int *out, i64 *Sv, i64 *ev, i64 *Tv, i64 *fv, i64 *a, i64 *q) {
  int N = 1 << m; memset(dpS, 0, sizeof(i64) * N * m); for (int v = 0; v < m; v++) dpS[(1 << v) * m + v] = 1;
  for (int S = 1; S < N; S++) for (int v = 0; v < m; v++) { i64 x = dpS[S * m + v]; if (!x) continue; int nb = out[v] & ~S; while (nb) { int u = __builtin_ctz(nb); nb &= nb - 1; dpS[(S | 1 << u) * m + u] += x; } }
  Hsub[0] = 1; for (int S = 1; S < N; S++) { i64 s = 0; for (int v = 0; v < m; v++) s += dpS[S * m + v]; Hsub[S] = s; }
  *a = Hsub[N - 1]; i64 qq = 0; for (int S = 0; S < N; S++) qq += Hsub[S] * Hsub[(N - 1) ^ S]; *q = qq;
  for (int v = 0; v < m; v++) Sv[v] = ev[v] = Tv[v] = fv[v] = 0;
  for (int S = 1; S < N; S++) { int sg = (__builtin_popcount(S) & 1) ? 1 : -1; for (int v = 0; v < m; v++) { i64 x = dpS[S * m + v]; if (!x) continue; i64 w = x * Hsub[(N - 1) ^ S]; Tv[v] += w; fv[v] += sg * w; } }
  int in[MAXM]; for (int v = 0; v < m; v++) in[v] = 0; for (int v = 0; v < m; v++) for (int u = 0; u < m; u++) if ((out[v] >> u) & 1) in[u] |= 1 << v;
  memset(dpR, 0, sizeof(i64) * N * m); for (int v = 0; v < m; v++) dpR[(1 << v) * m + v] = 1;
  for (int S = 1; S < N; S++) for (int v = 0; v < m; v++) { i64 x = dpR[S * m + v]; if (!x) continue; int nb = in[v] & ~S; while (nb) { int u = __builtin_ctz(nb); nb &= nb - 1; dpR[(S | 1 << u) * m + u] += x; } }
  for (int S = 1; S < N; S++) { int sg = (__builtin_popcount(S) & 1) ? 1 : -1; for (int v = 0; v < m; v++) { i64 x = dpR[S * m + v]; if (!x) continue; i64 w = x * Hsub[(N - 1) ^ S]; Sv[v] += w; ev[v] += sg * w; } } }
static i64 apexcount(int m, const int *out, int V) { int o2[MAXM + 1]; for (int v = 0; v < m; v++) o2[v] = out[v] | (((V >> v) & 1) ? 0 : (1 << m)); o2[m] = V; return hamcount(m + 1, o2); }
static void parse(const char *line, int n, int *out) { int k = 0; for (int i = 0; i < n; i++) out[i] = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') out[i] |= 1 << j; else out[j] |= 1 << i; k++; } }
static int strong(int n, const int *out) { int reach = 1, fr = 1; while (fr) { int nf = 0; for (int v = 0; v < n; v++) if ((fr >> v) & 1) nf |= out[v]; nf &= ~reach; reach |= nf; fr = nf; }
  if (reach != (1 << n) - 1) return 0; int in[16]; for (int v = 0; v < n; v++) in[v] = 0; for (int v = 0; v < n; v++) for (int u = 0; u < n; u++) if ((out[v] >> u) & 1) in[u] |= 1 << v;
  reach = 1; fr = 1; while (fr) { int nf = 0; for (int v = 0; v < n; v++) if ((fr >> v) & 1) nf |= in[v]; nf &= ~reach; reach |= nf; fr = nf; } return reach == (1 << n) - 1; }

static i64 pfun(i64 u){ if(u==0) return 0; if(u&1) return pfun(u>>1); i64 h=u>>1; return pfun(h)+h+(h&1); }
int main(void){ int bad=0, done=0, cores=0; i64 kmin=1LL<<60,kmax=-(1LL<<60);
  for(int m=3;m<=5;m++){ char cmd[64]; sprintf(cmd,"gentourng -c %d 2>/dev/null",m); FILE*fp=popen(cmd,"r"); char line[512];
    while(fgets(line,sizeof line,fp)){ if(line[0]=='>') continue; int out[MAXM]; parse(line,m,out); if(!strong(m,out)) continue;
      i64 Sv[MAXM],ev[MAXM],Tv[MAXM],fv[MAXM],a,q; vdata(m,out,Sv,ev,Tv,fv,&a,&q);
      for(int V=0;V<(1<<m);V++){ blk C; C.a=a; C.q=q; C.s=C.e=C.t=C.f=0;
        for(int v=0;v<m;v++){ if((V>>v)&1){C.s+=Sv[v];C.e+=ev[v];} else {C.t+=Tv[v];C.f+=fv[v];} } C.n=apexcount(m,out,V); cores++;
        i64 kappa=0; int first=1;
        for(i64 X=1;X<=70;X++) for(i64 Y=1;Y<=70;Y++){
          // build spine: guard 1, then bits of X after its leading 1 (MSB first), core, then y-bits, guard 0, where
          // Y = little-endian unmarked indicator including the final guard as its top bit
          st s={1,1,0,0}; s=step1(s,1); int lx=63-__builtin_clzll(X); for(int b=lx-1;b>=0;b--) s=step1(s,(X>>b)&1);
          s=stepb(s,C); int ly=63-__builtin_clzll(Y); for(int j=0;j<ly;j++){ int unmarked=(Y>>j)&1; s=step1(s,1-unmarked); } s=step1(s,0);
          i64 F=C.q*X*Y+C.t*X+C.s*Y+C.a*(pfun(X)+pfun(Y))+C.f*(X&1)+C.e*(Y&1);
          i64 k=s.N-F; if(first){kappa=k; first=0;} else if(k!=kappa) bad++; done++; }
        if(kappa<kmin) kmin=kappa; if(kappa>kmax) kmax=kappa; if(kappa!=C.n) { static int nd=0; nd++; if(nd<5) printf("kappa %lld != n %lld\n",kappa,C.n); }
        if(V==1&&m==3) printf("3-cycle, one marked vertex: (a,q,s,e,t,f,n)=(%lld,%lld,%lld,%lld,%lld,%lld,%lld), kappa=%lld\n",C.a,C.q,C.s,C.e,C.t,C.f,C.n,kappa);
      } }
    pclose(fp); }
  printf("closed form checked on %d (core,X,Y) triples over %d marked cores: %d inconsistencies; kappa range [%lld,%lld]\n",done,cores,bad,kmin,kmax); return 0; }
