/* Independent audit counters (S19 audit). Written from scratch; shares no code with port_stats.c or packing_count.c.
 * Tournament on n <= 16 vertices given by out-neighbour bitmasks out[v] (bit w set iff v -> w).
 *  - cyc[S]   : number of directed Hamiltonian cycles of T[S] (|S| >= 3); start fixed at the HIGHEST vertex of S.
 *  - I_x[S]   : I(Omega(T[S]), x) = sum over sets of vertex-disjoint odd directed cycles inside S of x^(#cycles),
 *               recursion on the HIGHEST vertex of S.
 *  - ham(T)   : number of directed Hamiltonian paths (subset DP), and hend[v] = Hamiltonian paths ending at v.
 *  - port sums: E, O = sums of I_x[V - V(P)] over directed paths P starting at u of even / odd order (P = (u) included).
 */
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
typedef long long ll;

#define MAXN 16
static int TN;                 /* order */
static uint32_t TOUT[MAXN];    /* out-neighbour masks */
static ll *CYC;                /* size 2^n */

static void tc_alloc(int n) {
  static int cap = -1;
  if (cap < n) { free(CYC); CYC = (ll *)calloc((size_t)1 << n, sizeof(ll)); cap = n; }
}

/* cyc[S] for all S, highest vertex as anchor */
static void tc_cycles(void) {
  int n = TN; size_t N = (size_t)1 << n;
  memset(CYC, 0, N * sizeof(ll));
  ll *W = (ll *)calloc(((size_t)1 << (n > 0 ? n - 1 : 0)) * (size_t)n + 1, sizeof(ll));
  for (int t = 2; t < n; t++) {           /* anchor t, other vertices from {0..t-1} */
    size_t M = (size_t)1 << t;
    /* W[mask*n + v]: paths t -> ... -> v visiting exactly {t} u mask, v in mask */
    memset(W, 0, M * (size_t)n * sizeof(ll));
    for (int v = 0; v < t; v++) if (TOUT[t] >> v & 1) W[((size_t)1 << v) * n + v] = 1;
    for (size_t mask = 1; mask < M; mask++) {
      for (int v = 0; v < t; v++) {
        ll c = W[mask * n + v]; if (!c) continue;
        if (__builtin_popcountll(mask) >= 2 && (TOUT[v] >> t & 1)) CYC[mask | ((size_t)1 << t)] += c;
        uint32_t ext = TOUT[v] & (uint32_t)(M - 1) & ~(uint32_t)mask;
        while (ext) { int w = __builtin_ctz(ext); ext &= ext - 1; W[(mask | ((size_t)1 << w)) * n + w] += c; }
      }
    }
  }
  free(W);
}

/* I(Omega(T[S]), x) for all S into Ix (size 2^n) ; uses CYC */
static void tc_packing(ll x, ll *Ix) {
  int n = TN; size_t N = (size_t)1 << n;
  Ix[0] = 1;
  for (size_t S = 1; S < N; S++) {
    int t = 63 - __builtin_clzll((unsigned long long)S);
    size_t rest = S & ~((size_t)1 << t);
    ll val = Ix[rest];
    /* sub ranges over subsets of rest with even popcount >= 2 */
    for (size_t sub = rest; sub; sub = (sub - 1) & rest) {
      int pc = __builtin_popcountll(sub);
      if (pc & 1) continue;
      size_t U = sub | ((size_t)1 << t);
      ll c = CYC[U]; if (!c) continue;
      val += x * c * Ix[rest & ~sub];
    }
    Ix[S] = val;
  }
}

/* Hamiltonian paths: returns H(T); fills hend[v] = Hamiltonian paths ending at v */
static ll tc_ham(ll *hend) {
  int n = TN; size_t N = (size_t)1 << n;
  if (n == 0) return 1;
  ll *D = (ll *)calloc(N * (size_t)n, sizeof(ll));
  for (int v = 0; v < n; v++) D[((size_t)1 << v) * n + v] = 1;
  for (size_t mask = 1; mask < N; mask++)
    for (int v = 0; v < n; v++) {
      ll c = D[mask * n + v]; if (!c) continue;
      uint32_t ext = TOUT[v] & ~(uint32_t)mask;
      while (ext) { int w = __builtin_ctz(ext); ext &= ext - 1; D[(mask | ((size_t)1 << w)) * n + w] += c; }
    }
  ll H = 0;
  for (int v = 0; v < n; v++) { hend[v] = D[(N - 1) * n + v]; H += hend[v]; }
  free(D);
  return H;
}

/* port sums for port u using a table Ix of I(Omega(T[S]),x) */
static void tc_port(int u, const ll *Ix, ll *E, ll *O) {
  int n = TN; size_t N = (size_t)1 << n, FULL = N - 1;
  ll *D = (ll *)calloc(N * (size_t)n, sizeof(ll));
  D[((size_t)1 << u) * n + u] = 1;
  ll e = 0, o = 0;
  for (size_t mask = 1; mask < N; mask++) {
    if (!(mask >> u & 1)) continue;
    ll tot = 0;
    for (int v = 0; v < n; v++) {
      ll c = D[mask * n + v]; if (!c) continue;
      tot += c;
      uint32_t ext = TOUT[v] & ~(uint32_t)mask;
      while (ext) { int w = __builtin_ctz(ext); ext &= ext - 1; D[(mask | ((size_t)1 << w)) * n + w] += c; }
    }
    if (!tot) continue;
    ll wgt = tot * Ix[FULL & ~mask];
    if (__builtin_popcountll(mask) % 2 == 0) e += wgt; else o += wgt;
  }
  free(D);
  *E = e; *O = o;
}
