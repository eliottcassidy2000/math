/* AUDIT (independent): THM-5071 (1), (2), (3), (4a) for loopless digraphs.
 * For each digraph D (adjacency bits, no loops) compute
 *   H    = Hamiltonian paths (endpoint DP);
 *   GS1  = (1) organised by cycle vertex sets: sum over set partitions of V of prod_U w1(U), singletons 1, where
 *          w1(U) = #directed Ham cycles of D[U] + (-1)^(|U|-1) #directed Ham cycles of Dbar[U]  (|U| >= 2);
 *   PERM = (1) by explicit enumeration of all n! permutations (n <= 6 only), a second, independent route;
 *   GS2  = (2): explicit DFS enumeration of the directed cycles of D and of Dbar (length >= 2), each classified by
 *          its pair types: D-cycle with a B-pair -> 1, D-cycle with T-pairs only -> 2 (odd) / 0 (even),
 *          Dbar-cycle with a Z-pair -> (-1)^(len-1), Dbar-cycle with T-pairs only -> not admissible;
 *          then summed over vertex-disjoint packings (partition DP);
 *   PAR  = (3): number of packings of defect cycles (B-cycles of D, Z-cycles of Dbar), mod 2.
 * Checks H == GS1 == GS2 (== PERM when n <= 6) and H = PAR (mod 2).
 * Exhaustive: every loopless digraph of order <= 5 (4^10 = 1,048,576 at order 5). Random: orders 6, 7, 8.
 * (4a): TT_n plus sink back arcs to S, every S, n <= 10 (array bound): H = 1 + sum_{x in S} 2^(n-2-x).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef long long ll;
static int n, A[10][10];
static ll cD[1 << 10], cDb[1 << 10], w2[1 << 10], par[1 << 10], R[1 << 10];
static ll dp[1 << 10][10];

static ll ham(void) {
  int N = 1 << n;
  for (int S = 0; S < N; S++) for (int v = 0; v < n; v++) dp[S][v] = 0;
  for (int v = 0; v < n; v++) dp[1 << v][v] = 1;
  for (int S = 1; S < N; S++) for (int v = 0; v < n; v++) if (dp[S][v])
    for (int t = 0; t < n; t++) if (!(S >> t & 1) && A[v][t]) dp[S | 1 << t][t] += dp[S][v];
  ll h = 0;
  for (int v = 0; v < n; v++) h += dp[N - 1][v];
  return h;
}

/* directed Hamiltonian cycles of G[U] for all U (|U| >= 2), G = D (comp=0) or Dbar (comp=1) */
static void hamcycles(int comp, ll *out) {
  int N = 1 << n;
  memset(out, 0, sizeof(ll) * N);
  for (int s = 0; s < n; s++) {
    for (int S = 0; S < N; S++) for (int v = 0; v < n; v++) dp[S][v] = 0;
    dp[1 << s][s] = 1;
    for (int S = 1 << s; S < N; S++) {
      if (!(S >> s & 1) || (S & ((1 << s) - 1))) continue;
      for (int v = 0; v < n; v++) {
        ll c = dp[S][v];
        if (!c) continue;
        int back = comp ? (v != s && !A[v][s]) : A[v][s];
        if (__builtin_popcount(S) >= 2 && back) out[S] += c;
        for (int t = s + 1; t < n; t++) {
          int arc = comp ? !A[v][t] : A[v][t];
          if (!(S >> t & 1) && arc) dp[S | 1 << t][t] += c;
        }
      }
    }
  }
}

static ll partition_sum(ll *w) { /* sum over set partitions of prod w(U), singletons weight 1 */
  int N = 1 << n;
  R[0] = 1;
  for (int S = 1; S < N; S++) {
    int m = __builtin_ctz(S), rest = S & ~(1 << m);
    ll val = R[rest];
    for (int T = rest; T; T = (T - 1) & rest) { int U = T | (1 << m); if (w[U]) val += w[U] * R[S & ~U]; }
    R[S] = val;
  }
  return R[N - 1];
}

/* explicit DFS over cycles from the minimum vertex s, recording weights per vertex set */
static int cyc[10], clen;
static void dfs(int comp, int s, int v, int mask) {
  for (int t = s; t < n; t++) {
    int arc = comp ? (t != v && !A[v][t]) : A[v][t];
    if (!arc) continue;
    if (t == s) {
      if (clen < 2) continue;
      int hasB = 0, hasZ = 0;
      for (int i = 0; i < clen; i++) { int a = cyc[i], b = cyc[(i + 1) % clen];
        if (A[a][b] && A[b][a]) hasB = 1; if (!A[a][b] && !A[b][a]) hasZ = 1; }
      ll wt;
      if (!comp) wt = hasB ? 1 : (clen % 2 ? 2 : 0);
      else wt = hasZ ? ((clen - 1) % 2 ? -1 : 1) : 0;
      w2[mask] += wt;
      if ((!comp && hasB) || (comp && hasZ)) par[mask] += 1;
      continue;
    }
    if (mask >> t & 1) continue;
    cyc[clen++] = t;
    dfs(comp, s, t, mask | 1 << t);
    clen--;
  }
}

static int perm[10], used[10];
static ll permsum;
static void prec(int i) {
  if (i == n) {
    int seen[10] = {0}; ll w = 1;
    for (int s = 0; s < n && w; s++) if (!seen[s]) {
      int len = 0, inD = 1, inDb = 1, v = s;
      do { seen[v] = 1; int u = perm[v]; len++; if (u != v) { if (!A[v][u]) inD = 0; else inDb = 0; } v = u; } while (v != s);
      if (len > 1) w *= inD ? 1 : (inDb ? ((len - 1) % 2 ? -1 : 1) : 0);
    }
    permsum += w; return;
  }
  for (int v = 0; v < n; v++) if (!used[v]) { used[v] = 1; perm[i] = v; prec(i + 1); used[v] = 0; }
}

static unsigned long long st = 0x5DEECE66DULL;
static unsigned long long rnd(void) { unsigned long long z = (st += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; return z ^ (z >> 31); }

static ll checks, fails, pfails, permchecks;
static void check(void) {
  int N = 1 << n;
  ll H = ham();
  hamcycles(0, cD); hamcycles(1, cDb);
  static ll w1[1 << 10];
  for (int U = 0; U < N; U++) w1[U] = (__builtin_popcount(U) >= 2) ? cD[U] + (((__builtin_popcount(U) - 1) % 2) ? -cDb[U] : cDb[U]) : 0;
  ll g1 = partition_sum(w1);
  memset(w2, 0, sizeof(ll) * N); memset(par, 0, sizeof(ll) * N);
  for (int s = 0; s < n; s++) { cyc[0] = s; clen = 1; dfs(0, s, s, 1 << s); dfs(1, s, s, 1 << s); }
  ll g2 = partition_sum(w2);
  ll pp = partition_sum(par);
  int ok = (H == g1 && H == g2);
  if (n <= 6) { permsum = 0; memset(used, 0, sizeof used); prec(0); permchecks++; if (permsum != H) ok = 0; }
  checks++;
  if (!ok) { fails++; if (fails < 5) printf("FAIL n=%d H=%lld GS1=%lld GS2=%lld PERM=%lld\n", n, H, g1, g2, permsum); }
  if ((H - pp) % 2) { pfails++; if (pfails < 5) printf("PARITY FAIL n=%d H=%lld packings=%lld\n", n, H, pp); }
}

int main(void) {
  for (n = 1; n <= 5; n++) {
    int m = n * (n - 1) / 2; ll total = 1LL << (2 * m);
    for (ll code = 0; code < total; code++) {
      memset(A, 0, sizeof A); int k = 0;
      for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++, k++) { int s2 = (code >> (2 * k)) & 3; A[i][j] = s2 & 1; A[j][i] = s2 >> 1; }
      check();
    }
    printf("all loopless digraphs of order %d: cumulative checks %lld, formula failures %lld, parity failures %lld\n", n, checks, fails, pfails);
    fflush(stdout);
  }
  int reps[3] = {200000, 100000, 20000};
  for (n = 6; n <= 8; n++) {
    for (int r = 0; r < reps[n - 6]; r++) {
      memset(A, 0, sizeof A);
      unsigned long long b = rnd(); int mode = r % 4; /* 0 uniform, 1 semicomplete, 2 oriented, 3 near-tournament */
      for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
        int s2;
        if (mode == 0) s2 = rnd() & 3;
        else if (mode == 1) s2 = 1 + rnd() % 3;
        else if (mode == 2) s2 = rnd() % 3;
        else { s2 = 1 + (rnd() & 1); if (rnd() % 8 == 0) s2 = (rnd() & 1) ? 0 : 3; }
        A[i][j] = s2 & 1; A[j][i] = s2 >> 1; (void)b;
      }
      check();
    }
    printf("random digraphs of order %d (uniform / semicomplete / oriented / tournament+few defects): cumulative checks %lld,"
           " formula failures %lld, parity failures %lld\n", n, checks, fails, pfails);
    fflush(stdout);
  }
  printf("explicit-permutation route used on %lld of them (orders <= 6)\n", permchecks);
  /* (4a) sink back arcs */
  ll sk = 0, skf = 0;
  for (n = 2; n <= 10; n++) for (int S = 0; S < (1 << (n - 1)); S++) {
    memset(A, 0, sizeof A);
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) A[i][j] = 1;
    ll pred = 1;
    for (int x = 0; x < n - 1; x++) if (S >> x & 1) { A[n - 1][x] = 1; pred += 1LL << (n - 2 - x); }
    sk++; if (ham() != pred) skf++;
  }
  printf("(4a) sink back-arc construction, all S, orders 2..10: %lld digraphs, %lld failures\n", sk, skf);
  return 0;
}
