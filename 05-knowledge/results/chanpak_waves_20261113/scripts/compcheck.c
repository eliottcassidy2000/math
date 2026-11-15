// Wave 26: passive-palette census. For every class R of order m (gentourng) and every mark set V
// (V = vertices beaten by the apex), print the block data of the exact cut transfer
//   a = H(R), q = P_R(2), (s,e) = total / odd-minus-even start weights from V,
//   (t,f) = total / odd-minus-even end weights into the complement of V, n = H(R + apex).
// Start weight of a nonempty path P = H(R - P). Mode "check" verifies the cut transfer on random
// two-block spines against direct Hamiltonian-path DP.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64;
#define MAXM 12

static i64 Hsub[1 << MAXM];
static i64 dpS[(1 << MAXM) * MAXM];

static i64 hamcount(int n, const int *out) {
  int N = 1 << n; i64 *dp = calloc((size_t)N * n, sizeof(i64));
  for (int v = 0; v < n; v++) dp[(size_t)(1 << v) * n + v] = 1;
  for (int S = 1; S < N; S++) for (int v = 0; v < n; v++) { i64 x = dp[(size_t)S * n + v]; if (!x) continue;
    int nb = out[v] & ~S; while (nb) { int u = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | 1 << u) * n + u] += x; } }
  i64 s = 0; for (int v = 0; v < n; v++) s += dp[(size_t)(N - 1) * n + v]; free(dp); return s;
}

// per-vertex start and end weights
static void vertex_data(int m, const int *out, i64 *Sv, i64 *ev, i64 *Tv, i64 *fv, i64 *a, i64 *q) {
  int N = 1 << m;
  // Hsub[S] = H(R[S])
  memset(dpS, 0, sizeof(i64) * N * m);
  for (int v = 0; v < m; v++) dpS[(1 << v) * m + v] = 1;
  for (int S = 1; S < N; S++) for (int v = 0; v < m; v++) { i64 x = dpS[S * m + v]; if (!x) continue;
    int nb = out[v] & ~S; while (nb) { int u = __builtin_ctz(nb); nb &= nb - 1; dpS[(S | 1 << u) * m + u] += x; } }
  Hsub[0] = 1; for (int S = 1; S < N; S++) { i64 s = 0; for (int v = 0; v < m; v++) s += dpS[S * m + v]; Hsub[S] = s; }
  *a = Hsub[N - 1];
  i64 qq = 0; for (int S = 0; S < N; S++) qq += Hsub[S] * Hsub[(N - 1) ^ S]; *q = qq;
  for (int v = 0; v < m; v++) { Sv[v] = ev[v] = Tv[v] = fv[v] = 0; }
  // ends: dpS[S*m+v] = # Ham paths of R[S] ending at v
  for (int S = 1; S < N; S++) { int pc = __builtin_popcount(S); int sg = (pc & 1) ? 1 : -1;
    for (int v = 0; v < m; v++) { i64 x = dpS[S * m + v]; if (!x) continue; i64 w = x * Hsub[(N - 1) ^ S];
      Tv[v] += w; fv[v] += sg * w; } }
  // starts: paths starting at v covering S = Ham paths of reversed R[S] ending at v
  static i64 dpR[(1 << MAXM) * MAXM];
  memset(dpR, 0, sizeof(i64) * N * m);
  int in[MAXM]; for (int v = 0; v < m; v++) in[v] = 0;
  for (int v = 0; v < m; v++) for (int u = 0; u < m; u++) if ((out[v] >> u) & 1) in[u] |= 1 << v;
  for (int v = 0; v < m; v++) dpR[(1 << v) * m + v] = 1;
  for (int S = 1; S < N; S++) for (int v = 0; v < m; v++) { i64 x = dpR[S * m + v]; if (!x) continue;
    int nb = in[v] & ~S; while (nb) { int u = __builtin_ctz(nb); nb &= nb - 1; dpR[(S | 1 << u) * m + u] += x; } }
  for (int S = 1; S < N; S++) { int pc = __builtin_popcount(S); int sg = (pc & 1) ? 1 : -1;
    for (int v = 0; v < m; v++) { i64 x = dpR[S * m + v]; if (!x) continue; i64 w = x * Hsub[(N - 1) ^ S];
      Sv[v] += w; ev[v] += sg * w; } }
}

static void parse(const char *line, int n, int *out) {
  int k = 0; for (int i = 0; i < n; i++) out[i] = 0;
  for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') out[i] |= 1 << j; else out[j] |= 1 << i; k++; }
}

// apex beats V: apex index m
static i64 apexcount(int m, const int *out, int V) {
  int o2[MAXM + 1]; for (int v = 0; v < m; v++) o2[v] = out[v] | (((V >> v) & 1) ? 0 : (1 << m)); o2[m] = V;
  return hamcount(m + 1, o2);
}


int main(void){
  srand(99); int bad=0;
  for(int trial=0; trial<400; trial++){
    int m1=1+rand()%4, m2=1+rand()%4, n=m1+m2; int out[MAXM]; for(int i=0;i<n;i++) out[i]=0;
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ int bi=i<m1, bj=j<m1; if(bi!=bj){ out[i]|=1<<j; continue; } if(rand()&1) out[i]|=1<<j; else out[j]|=1<<i; }
    int V=rand()%(1<<n);
    i64 S[MAXM],e[MAXM],T[MAXM],f[MAXM],a,q; vertex_data(n,out,S,e,T,f,&a,&q);
    i64 s=0,ee=0,t=0,ff=0; for(int v=0;v<n;v++){ if((V>>v)&1){s+=S[v];ee+=e[v];} else {t+=T[v];ff+=f[v];} }
    int o1[MAXM],o2[MAXM]; for(int i=0;i<m1;i++) o1[i]=out[i]&((1<<m1)-1); for(int i=0;i<m2;i++) o2[i]=(out[m1+i]>>m1)&((1<<m2)-1);
    int V1=V&((1<<m1)-1), V2=V>>m1;
    i64 S1[MAXM],e1[MAXM],T1[MAXM],f1[MAXM],a1,q1; vertex_data(m1,o1,S1,e1,T1,f1,&a1,&q1);
    i64 s1=0,ee1=0,t1=0,ff1=0; for(int v=0;v<m1;v++){ if((V1>>v)&1){s1+=S1[v];ee1+=e1[v];} else {t1+=T1[v];ff1+=f1[v];} }
    i64 S2[MAXM],e2[MAXM],T2[MAXM],f2[MAXM],a2,q2; vertex_data(m2,o2,S2,e2,T2,f2,&a2,&q2);
    i64 s2=0,ee2=0,t2=0,ff2=0; for(int v=0;v<m2;v++){ if((V2>>v)&1){s2+=S2[v];ee2+=e2[v];} else {t2+=T2[v];ff2+=f2[v];} }
    if(a!=a1*a2||q!=q1*q2||s!=q2*s1+a1*s2||ee!=a1*ee2||t!=q1*t2+a2*t1||ff!=a2*ff1) { bad++; if(bad<5) printf("fail trial %d\n",trial); }
    // mod-3 law per vertex
    for(int v=0;v<n;v++){ i64 lhs=((T[v]-((n%2)?1:-1)*S[v])%3+3)%3; if(lhs){ bad++; printf("mod3 fail\n"); break; } }
  }
  printf("composition rules (a,q,s,e,t,f) and T_v = (-1)^(m+1) S_v mod 3: %d failures in 400 random two-block spines\n",bad);
  return 0; }
