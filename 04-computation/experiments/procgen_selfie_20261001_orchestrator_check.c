/* Orchestrator's independent audit engine for the selfie lane (2026-10-01).
 * Written from the definitions; the lane's code was not read.
 * Tournament encoding (gentourng and labeled modes): pairs (i,j), i<j, in row-major order;
 * bit/char 1 means i->j.
 * Driver: procgen_selfie_20261001_orchestrator_check.py (compiles this file into a temp dir).
 * Modes:
 *   sa labeled N        all 2^C(N,2) labeled tournaments, N<=7
 *   sa classes          gentourng lines on stdin
 *   sa one STRING       exact arc counts of one tournament (N<=20)
 *   sa par2 STRING      arc-count parities by bitset DP (N<=24)
 *   sa paley q mode     QR_q (mode 0) or QR_q minus vertex 0 (mode 1), q prime = 3 mod 4; exact if q<=19 else mod 2
 *   sa circ N           all circulant tournaments on Z_N (N odd)
 *   sa z3z3             all Cayley tournaments on Z_3 x Z_3
 *   sa beta STRING      HP-blocking number beta, hall bound, sigma bound
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef unsigned long long u64;
static int N;
static uint32_t OUT[32], IN[32];
static u64 *F, *B, *G;
static u64 C[32][32];
static long fails = 0;

#define FAIL(...) do { fails++; if (fails < 20) { printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } } while (0)

static void alloc(int n) { size_t s = (size_t)(1u << n) * n; F = malloc(s * 8); B = malloc(s * 8); G = malloc(s * 8);
  if (!F || !B || !G) { fprintf(stderr, "oom\n"); exit(1); } }

static void set_from_string(const char *s) {
  int len = strlen(s); N = 0; while (N * (N - 1) / 2 < len) N++;
  if (N * (N - 1) / 2 != len) { fprintf(stderr, "bad length %d\n", len); exit(1); }
  memset(OUT, 0, sizeof OUT); memset(IN, 0, sizeof IN);
  int k = 0;
  for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++, k++) {
    if (s[k] == '1') { OUT[i] |= 1u << j; IN[j] |= 1u << i; } else { OUT[j] |= 1u << i; IN[i] |= 1u << j; }
  }
}
static void set_from_code(u64 code) {
  memset(OUT, 0, sizeof OUT); memset(IN, 0, sizeof IN);
  int k = 0;
  for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++, k++) {
    if (code >> k & 1) { OUT[i] |= 1u << j; IN[j] |= 1u << i; } else { OUT[j] |= 1u << i; IN[i] |= 1u << j; }
  }
}
/* F[m][v] = #directed paths with vertex set m ending at v; B[m][v] = ... starting at v */
static void dp(void) {
  int full = 1 << N; size_t s = (size_t)full * N;
  memset(F, 0, s * 8); memset(B, 0, s * 8);
  for (int v = 0; v < N; v++) { F[(size_t)(1 << v) * N + v] = 1; B[(size_t)(1 << v) * N + v] = 1; }
  for (int m = 1; m < full; m++) {
    for (int v = 0; v < N; v++) {
      if (!(m >> v & 1)) continue;
      u64 fv = F[(size_t)m * N + v], bv = B[(size_t)m * N + v];
      if (fv) { uint32_t t = OUT[v] & ~(uint32_t)m; while (t) { int w = __builtin_ctz(t); t &= t - 1; F[(size_t)(m | 1 << w) * N + w] += fv; } }
      if (bv) { uint32_t t = IN[v] & ~(uint32_t)m; while (t) { int w = __builtin_ctz(t); t &= t - 1; B[(size_t)(m | 1 << w) * N + w] += bv; } }
    }
  }
}
static u64 Hsub(uint32_t U) { if (!U) return 1; u64 h = 0; for (int v = 0; v < N; v++) if (U >> v & 1) h += F[(size_t)U * N + v]; return h; }
/* arc counts inside the subtournament on S: C[u][v] = #HPs of T[S] through u->v */
static void arcs_sub(uint32_t S) {
  for (int u = 0; u < N; u++) for (int v = 0; v < N; v++) C[u][v] = 0;
  for (uint32_t A = (0 - S) & S; A && A != S; A = (A - S) & S) {
    uint32_t R = S ^ A;
    for (int u = 0; u < N; u++) if (A >> u & 1) {
      u64 fu = F[(size_t)A * N + u]; if (!fu) continue;
      uint32_t t = OUT[u] & R;
      while (t) { int v = __builtin_ctz(t); t &= t - 1; C[u][v] += fu * B[(size_t)R * N + v]; }
    }
  }
}
/* directed Hamiltonian cycles */
static u64 hamcycles(void) {
  int full = 1 << N; size_t s = (size_t)full * N; memset(G, 0, s * 8);
  G[1 * N + 0] = 1;
  for (int m = 1; m < full; m += 2) for (int v = 0; v < N; v++) {
    u64 g = G[(size_t)m * N + v]; if (!g) continue;
    uint32_t t = OUT[v] & ~(uint32_t)m; while (t) { int w = __builtin_ctz(t); t &= t - 1; G[(size_t)(m | 1 << w) * N + w] += g; }
  }
  u64 hc = 0; for (int v = 1; v < N; v++) if (OUT[v] & 1u) hc += G[(size_t)(full - 1) * N + v];
  return hc;
}
/* strong components in order C_1 => C_2 => ... ; returns k, comps[] */
static int components(uint32_t *comps) {
  uint32_t reach[32];
  for (int v = 0; v < N; v++) { uint32_t r = 1u << v, fr = r; while (fr) { uint32_t nf = 0; uint32_t t = fr; while (t) { int w = __builtin_ctz(t); t &= t - 1; nf |= OUT[w]; } nf &= ~r; r |= nf; fr = nf; } reach[v] = r; }
  uint32_t seen = 0; int k = 0;
  /* order by size of reach set, decreasing */
  for (;;) {
    int best = -1; for (int v = 0; v < N; v++) if (!(seen >> v & 1)) { if (best < 0 || __builtin_popcount(reach[v]) > __builtin_popcount(reach[best])) best = v; }
    if (best < 0) break;
    uint32_t c = 0; for (int w = 0; w < N; w++) if ((reach[best] >> w & 1) && (reach[w] >> best & 1)) c |= 1u << w;
    comps[k++] = c; seen |= c;
  }
  return k;
}
/* z via formula, recursively on components */
static long zformula(void) {
  uint32_t comps[32]; int k = components(comps);
  long z = 0;
  for (int i = 0; i < k; i++) for (int j = i + 2; j < k; j++) z += (long)__builtin_popcount(comps[i]) * __builtin_popcount(comps[j]);
  for (int i = 0; i < k; i++) if (__builtin_popcount(comps[i]) >= 3) {
    arcs_sub(comps[i]);
    for (int u = 0; u < N; u++) if (comps[i] >> u & 1) { uint32_t t = OUT[u] & comps[i]; while (t) { int v = __builtin_ctz(t); t &= t - 1; if (C[u][v] == 0) z++; } }
  }
  return z;
}

typedef struct { long cls; u64 lab; } cnt;
static void addc(cnt *c, u64 w) { c->cls++; c->lab += w; }

/* ------------------------------------------------------------------ labeled mode */
static void labeled(int n) {
  N = n; alloc(N);
  int E = N * (N - 1) / 2; u64 total = 1ull << E; uint32_t full = (1u << N) - 1;
  u64 *Harr = malloc(total * 8);
  u64 *cut = malloc((size_t)(1u << N) * 8);
  for (uint32_t L = 0; L <= full; L++) { u64 m = 0; int k = 0; for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++, k++) if (((L >> i) ^ (L >> j)) & 1) m |= 1ull << k; cut[L] = m; }
  u64 cover = 0, allodd = 0, alleven = 0, bad = 0, zfail = 0, oddpar = 0, x0fail = 0, redeifail = 0;
  int minodd = 99, maxodd = -1;
  u64 nfact = 1; for (int i = 2; i <= N; i++) nfact *= i;
  for (u64 code = 0; code < total; code++) {
    set_from_code(code); dp(); arcs_sub(full);
    u64 H = Hsub(full); Harr[code] = H;
    int nodd = 0, ncov = 1, onall = 0; long z = 0;
    for (int u = 0; u < N; u++) { uint32_t t = OUT[u]; while (t) { int v = __builtin_ctz(t); t &= t - 1; if (C[u][v] & 1) nodd++; if (!C[u][v]) { ncov = 0; z++; } if (C[u][v] == H) onall = 1; } }
    if (!(H & 1)) FAIL("Redei H even code=%llu", code);
    if ((nodd & 1) != ((N - 1) & 1)) oddpar++;
    if (nodd < minodd) minodd = nodd; if (nodd > maxodd) maxodd = nodd;
    if (nodd == E) allodd++; if (nodd == 0) alleven++;
    if (ncov) cover++; if (ncov && onall) bad++;
    if (zformula() != z) zfail++;
    /* x=0 fixed-point OCF: sum_U (-1)^(N-|U|) H(T[U]) */
    long long s0 = 0; for (uint32_t U = 0; U <= full; U++) { long long h = (long long)Hsub(U); s0 += ((N - __builtin_popcount(U)) & 1) ? -h : h; }
    long long expect;
    if (N == 4) expect = 0;
    else if (N == 3 || N == 5 || N == 7) expect = 2 * (long long)hamcycles();
    else { /* N = 6: 4 * #partitions into two cyclic triples */
      long long cnt3 = 0;
      for (uint32_t S = 0; S <= full; S++) if (__builtin_popcount(S) == 3 && (S & 1)) {
        uint32_t R = full ^ S; int cyc1 = 0, cyc2 = 0;
        int a[3], b[3], ia = 0, ib = 0; for (int v = 0; v < N; v++) { if (S >> v & 1) a[ia++] = v; else b[ib++] = v; }
        cyc1 = (__builtin_popcount(OUT[a[0]] & S) == 1 && __builtin_popcount(OUT[a[1]] & S) == 1);
        cyc2 = (__builtin_popcount(OUT[b[0]] & R) == 1 && __builtin_popcount(OUT[b[1]] & R) == 1);
        if (cyc1 && cyc2) cnt3++;
      }
      expect = 4 * cnt3;
    }
    if (s0 != expect) x0fail++;
    /* selfie Redei (N <= 6): Z_L = sum_{M subset L} 2^|M| H(T - M) odd and = H + 2|L| mod 4 */
    if (N <= 6) for (uint32_t L = 0; L <= full; L++) {
      u64 Z = 0; for (uint32_t M = L;; M = (M - 1) & L) { Z += (1ull << __builtin_popcount(M)) * Hsub(full ^ M); if (!M) break; }
      if ((Z & 3) != ((H + 2 * __builtin_popcount(L)) & 3)) redeifail++;
    }
  }
  printf("labeled N=%d: total=%llu cover=%llu allodd=%llu alleven=%llu minodd=%d maxodd=%d bad(cover&on-all)=%llu zformula_fail=%llu oddcount_parity_fail=%llu x0_fail=%llu selfieRedei_fail=%llu%s\n",
         N, total, cover, allodd, alleven, minodd, maxodd, bad, zfail, oddpar, x0fail, redeifail, N <= 6 ? "" : "(not run)");
  /* switching sums */
  u64 sumfail = 0, altfail = 0;
  for (u64 code = 0; code < total; code++) {
    u64 s = 0; long long a = 0;
    for (uint32_t L = 0; L <= full; L++) { u64 h = Harr[code ^ cut[L]]; s += h; a += (__builtin_popcount(L) & 1) ? -(long long)h : (long long)h; }
    if (s != 2 * nfact) sumfail++; if (a != 0) altfail++;
  }
  printf("labeled N=%d: sum_L H(switch_L T) != 2N!: %llu ; sum_L (-1)^|L| H != 0: %llu\n", N, sumfail, altfail);
  /* Walsh degree on tiling representatives (base path i+1 -> i, i.e. pair (i,i+1) bit 0) */
  u64 basemask = 0; { int k = 0; for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++, k++) if (j == i + 1) basemask |= 1ull << k; }
  long long *W = malloc((size_t)(1u << N) * 8);
  u64 reps = 0, degfail = 0, topatt = 0, consth = 0; int topdeg = 2 * (N / 3);
  for (u64 code = 0; code < total; code++) {
    if (code & basemask) continue;
    reps++;
    for (uint32_t L = 0; L <= full; L++) W[L] = (long long)Harr[code ^ cut[L]];
    int isconst = 1; for (uint32_t L = 1; L <= full; L++) if (W[L] != W[0]) isconst = 0;
    if (isconst) { consth++; printf("  constant-H switching class N=%d code=%llu H=%lld\n", N, code, W[0]); }
    for (int h = 1; h <= (int)full; h <<= 1) for (uint32_t i = 0; i <= full; i += 2 * h) for (uint32_t j = i; j < i + h; j++) { long long x = W[j], y = W[j + h]; W[j] = x + y; W[j + h] = x - y; }
    int maxd = 0;
    for (uint32_t A = 0; A <= full; A++) if (W[A]) { int d = __builtin_popcount(A); if ((d & 1) || 3 * d > 2 * N) degfail++; if (d > maxd) maxd = d; }
    if (maxd == topdeg) topatt++;
    if (W[0] != (long long)nfact * 2) FAIL("W0");
  }
  printf("labeled N=%d: tiling reps=%llu walsh_degree_fail=%llu top_degree_%d_attained=%llu constant_H_classes=%llu\n", N, reps, degfail, topdeg, topatt, consth);
  free(W); free(Harr); free(cut);
}

/* ------------------------------------------------------------------ classes mode */
static int betaval(int maxk);
static long hallval(long *sigma);
static int shave_ok(void);
static long noshave = 0;
static void classes(int dobeta, int doshave) {
  char line[256]; int n0 = -1;
  long ncls = 0, cover = 0, allodd = 0, alleven = 0, strongdead = 0, zfail = 0, bh_fail = 0, hall_lt_sigma = 0;
  int minodd = 999, maxodd = -1; int maxbeta = -1;
  while (fgets(line, sizeof line, stdin)) {
    if (line[0] != '0' && line[0] != '1') continue;
    line[strcspn(line, "\r\n")] = 0;
    set_from_string(line);
    if (n0 < 0) { n0 = N; alloc(N); }
    uint32_t full = (1u << N) - 1;
    dp(); arcs_sub(full); u64 H = Hsub(full);
    int nodd = 0, ncov = 1; long z = 0; int E = N * (N - 1) / 2;
    for (int u = 0; u < N; u++) { uint32_t t = OUT[u]; while (t) { int v = __builtin_ctz(t); t &= t - 1; if (C[u][v] & 1) nodd++; if (!C[u][v]) { ncov = 0; z++; } } }
    ncls++;
    if (!(H & 1)) FAIL("Redei");
    if (nodd < minodd) minodd = nodd; if (nodd > maxodd) maxodd = nodd;
    if (nodd == E) { allodd++; printf("  ALLODD %s H=%llu\n", line, H); }
    if (nodd == 0) alleven++;
    if (ncov) cover++;
    uint32_t comps[32]; int k = components(comps);
    if (k == 1 && z > 0) strongdead++;
    if (zformula() != z) zfail++;
    if (dobeta) {
      long sigma; long hb = hallval(&sigma); int b = betaval(N);
      if (b != hb) { bh_fail++; printf("  beta!=hall %s beta=%d hall=%ld sigma=%ld\n", line, b, hb, sigma); }
      if (hb < sigma) hall_lt_sigma++;
      if (b > maxbeta) maxbeta = b;
      int reg = (N & 1); for (int v = 0; v < N; v++) if (__builtin_popcount(OUT[v]) != (N - 1) / 2) reg = 0;
      if (reg) printf("  REGULAR %s beta=%d\n", line, b);
      if (doshave) { int r = shave_ok(); if (!r) { noshave++; printf("  NO_REDEI_SHAVING %s H=%llu\n", line, H); } }
    }
  }
  printf("classes N=%d: n=%ld cover=%ld allodd=%ld alleven=%ld minodd=%d maxodd=%d strong_with_dead_arc=%ld zformula_fail=%ld", n0, ncls, cover, allodd, alleven, minodd, maxodd, strongdead, zfail);
  if (dobeta) printf(" beta!=hall=%ld hall<sigma=%ld max_beta=%d", bh_fail, hall_lt_sigma, maxbeta);
  if (doshave) printf(" no_redei_shaving=%ld", noshave);
  printf("\n");
}

/* ------------------------------------------------------------------ beta / hall */
static uint32_t DOUT[32];
static int has_hp(void) { /* reach[m] = set of v such that some path on m ends at v, arcs DOUT */
  static uint32_t R[1 << 16]; int full = 1 << N;
  memset(R, 0, sizeof(uint32_t) * full);
  for (int v = 0; v < N; v++) R[1 << v] = 1u << v;
  for (int m = 1; m < full; m++) { uint32_t r = R[m]; if (!r) continue; uint32_t ext = 0; uint32_t t = r; while (t) { int v = __builtin_ctz(t); t &= t - 1; ext |= DOUT[v]; } ext &= ~(uint32_t)m;
    while (ext) { int w = __builtin_ctz(ext); ext &= ext - 1; uint32_t t2 = r; int ok = 0; while (t2) { int v = __builtin_ctz(t2); t2 &= t2 - 1; if (DOUT[v] >> w & 1) { ok = 1; break; } } if (ok) R[m | 1 << w] |= 1u << w; } }
  return R[full - 1] != 0;
}
static int arcu[64], arcv[64], narcs;
static int rec(int start, int left) {
  if (left == 0) return !has_hp();
  for (int i = start; i < narcs; i++) { DOUT[arcu[i]] &= ~(1u << arcv[i]); int r = rec(i + 1, left - 1); DOUT[arcu[i]] |= 1u << arcv[i]; if (r) return 1; }
  return 0;
}
static int betaval(int maxk) {
  narcs = 0; for (int u = 0; u < N; u++) { uint32_t t = OUT[u]; while (t) { int v = __builtin_ctz(t); t &= t - 1; arcu[narcs] = u; arcv[narcs] = v; narcs++; } }
  for (int v = 0; v < N; v++) DOUT[v] = OUT[v];
  for (int k = 0; k <= maxk; k++) if (rec(0, k)) return k;
  return 99;
}
static long hallval(long *sigma) {
  uint32_t full = (1u << N) - 1; long best = 1 << 30; long sg = 1 << 30;
  for (uint32_t X = 1; X <= full; X++) {
    int x = __builtin_popcount(X); if (x < 2) continue;
    for (int dir = 0; dir < 2; dir++) {
      int w[32]; long tot = 0;
      for (int a = 0; a < N; a++) { w[a] = __builtin_popcount((dir ? IN[a] : OUT[a]) & X); tot += w[a]; }
      /* remove the x-2 largest */
      int ww[32]; memcpy(ww, w, sizeof w);
      for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) if (ww[j] > ww[i]) { int t = ww[i]; ww[i] = ww[j]; ww[j] = t; }
      long val = tot; for (int i = 0; i < x - 2; i++) val -= ww[i];
      if (val < best) best = val;
      if (x == 2 && val < sg) sg = val;
    }
  }
  *sigma = sg; return best;
}

/* ------------------------------------------------------------------ Redei shaving (N <= 7)
 * state = set of remaining arcs; a move deletes an arc e with c(e) even, so H stays odd (H(D - e) = H(D) - c(e));
 * success when H = 1. DFS with a visited bitmap over arc subsets. */
static int sarcu[32], sarcv[32], snarcs; static unsigned char *visited;
static uint32_t SOUT0[32], SIN0[32];
static void load_state(uint32_t st) {
  memset(OUT, 0, sizeof OUT); memset(IN, 0, sizeof IN);
  for (int i = 0; i < snarcs; i++) if (st >> i & 1) { OUT[sarcu[i]] |= 1u << sarcv[i]; IN[sarcv[i]] |= 1u << sarcu[i]; }
}
static int shave_dfs(uint32_t st) {
  if (visited[st >> 3] >> (st & 7) & 1) return 0;
  visited[st >> 3] |= 1 << (st & 7);
  load_state(st); dp(); uint32_t full = (1u << N) - 1; u64 H = Hsub(full);
  if (!(H & 1)) { fails++; return 0; }
  if (H == 1) return 1;
  arcs_sub(full);
  u64 cs[32]; for (int i = 0; i < snarcs; i++) cs[i] = (st >> i & 1) ? C[sarcu[i]][sarcv[i]] : 1;
  /* try arcs with the largest even c first */
  for (int pass = 0; pass < snarcs; pass++) {
    int best = -1; for (int i = 0; i < snarcs; i++) if ((st >> i & 1) && !(cs[i] & 1) && cs[i] != (u64)-1 && (best < 0 || cs[i] > cs[best])) best = i;
    if (best < 0) break;
    cs[best] = (u64)-1;
    if (shave_dfs(st & ~(1u << best))) return 1;
  }
  return 0;
}
static int shave_ok(void) {
  snarcs = 0; for (int u = 0; u < N; u++) { uint32_t t = OUT[u]; while (t) { int v = __builtin_ctz(t); t &= t - 1; sarcu[snarcs] = u; sarcv[snarcs] = v; snarcs++; } }
  memcpy(SOUT0, OUT, sizeof OUT); memcpy(SIN0, IN, sizeof IN);
  size_t nb = ((size_t)1 << snarcs) / 8 + 1; visited = calloc(nb, 1);
  int r = shave_dfs((uint32_t)(((u64)1 << snarcs) - 1));
  free(visited); memcpy(OUT, SOUT0, sizeof OUT); memcpy(IN, SIN0, sizeof IN);
  return r;
}

/* ------------------------------------------------------------------ parity DP (bitset), N <= 24 */
static void par2(void) {
  uint32_t full = (1u << N) - 1; size_t M = (size_t)1 << N;
  uint32_t *f = calloc(M, 4), *b = calloc(M, 4);
  if (!f || !b) { fprintf(stderr, "oom\n"); exit(1); }
  for (int v = 0; v < N; v++) { f[1u << v] = 1u << v; b[1u << v] = 1u << v; }
  for (uint32_t m = 1; m <= full; m++) {
    uint32_t fm = f[m], bm = b[m];
    for (int w = 0; w < N; w++) if (!(m >> w & 1)) {
      if (__builtin_popcount(fm & IN[w]) & 1) f[m | 1u << w] ^= 1u << w;
      if (__builtin_popcount(bm & OUT[w]) & 1) b[m | 1u << w] ^= 1u << w;
    }
  }
  uint32_t P[32]; memset(P, 0, sizeof P); /* P[u] bit v = parity of c(u->v) */
  for (uint32_t A = 1; A < full; A++) {
    uint32_t U = f[A]; if (!U) continue; uint32_t Rb = b[full ^ A]; if (!Rb) continue;
    while (U) { int u = __builtin_ctz(U); U &= U - 1; P[u] ^= Rb & OUT[u]; }
  }
  int nodd = 0; for (int u = 0; u < N; u++) nodd += __builtin_popcount(P[u] & OUT[u]);
  printf("par2 N=%d: H mod 2 = %d, #odd arcs = %d of %d\n", N, __builtin_popcount(f[full]) & 1, nodd, N * (N - 1) / 2);
  free(f); free(b);
}
static void one(void) {
  alloc(N); uint32_t full = (1u << N) - 1; dp(); arcs_sub(full); u64 H = Hsub(full);
  int nodd = 0; u64 mn = ~0ull, mx = 0;
  for (int u = 0; u < N; u++) { uint32_t t = OUT[u]; while (t) { int v = __builtin_ctz(t); t &= t - 1; if (C[u][v] & 1) nodd++; if (C[u][v] < mn) mn = C[u][v]; if (C[u][v] > mx) mx = C[u][v]; } }
  /* values multiset */
  printf("one N=%d: H=%llu #odd arcs=%d of %d, min c=%llu max c=%llu; c values:", N, H, nodd, N * (N - 1) / 2, mn, mx);
  u64 vals[64]; int nv = 0; for (int u = 0; u < N; u++) { uint32_t t = OUT[u]; while (t) { int v = __builtin_ctz(t); t &= t - 1; int f2 = 0; for (int i = 0; i < nv; i++) if (vals[i] == C[u][v]) f2 = 1; if (!f2 && nv < 64) vals[nv++] = C[u][v]; } }
  for (int i = 0; i < nv; i++) { int cnt = 0; for (int u = 0; u < N; u++) { uint32_t t = OUT[u]; while (t) { int v = __builtin_ctz(t); t &= t - 1; if (C[u][v] == vals[i]) cnt++; } } printf(" %llux%d", vals[i], cnt); }
  printf("\n");
}
static int isqr(int x, int q) { for (int y = 1; y < q; y++) if ((y * y) % q == x) return 1; return 0; }
static void paley(int q, int minus) {
  int n = minus ? q - 1 : q; N = n; memset(OUT, 0, sizeof OUT); memset(IN, 0, sizeof IN);
  int lab[32], k = 0; for (int x = 0; x < q; x++) if (!(minus && x == 0)) lab[k++] = x;
  for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j && isqr(((lab[j] - lab[i]) % q + q) % q, q)) { OUT[i] |= 1u << j; IN[j] |= 1u << i; }
  printf("QR_%d%s: ", q, minus ? " minus 0" : "");
  if (n <= 18) one(); else par2();
}
static void circ(int n) {
  N = n; alloc(N); int h = (n - 1) / 2; long tot = 0, alleven = 0;
  for (int S = 0; S < (1 << h); S++) {
    /* connection set: for d = 1..h, include d if bit set, else n-d */
    uint32_t conn = 0; for (int d = 1; d <= h; d++) conn |= 1u << ((S >> (d - 1) & 1) ? d : n - d);
    memset(OUT, 0, sizeof OUT); memset(IN, 0, sizeof IN);
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j && (conn >> (((j - i) % n + n) % n) & 1)) { OUT[i] |= 1u << j; IN[j] |= 1u << i; }
    uint32_t full = (1u << N) - 1; dp(); arcs_sub(full);
    int nodd = 0; for (int u = 0; u < N; u++) { uint32_t t = OUT[u]; while (t) { int v = __builtin_ctz(t); t &= t - 1; if (C[u][v] & 1) nodd++; } }
    tot++; if (nodd == 0) alleven++;
  }
  printf("circulants Z_%d: %ld tournaments, all-even: %ld\n", n, tot, alleven);
}
static void z3z3(void) {
  N = 9; alloc(N); long tot = 0, alleven = 0;
  /* nonzero elements of Z3xZ3 in 4 antipodal pairs: (1,0),(0,1),(1,1),(1,2) and negatives */
  int reps[4][2] = {{1, 0}, {0, 1}, {1, 1}, {1, 2}};
  for (int S = 0; S < 16; S++) {
    int conn[3][3] = {{0}};
    for (int i = 0; i < 4; i++) { int a = reps[i][0], b = reps[i][1]; if (S >> i & 1) conn[a][b] = 1; else conn[(3 - a) % 3][(3 - b) % 3] = 1; }
    memset(OUT, 0, sizeof OUT); memset(IN, 0, sizeof IN);
    for (int x = 0; x < 9; x++) for (int y = 0; y < 9; y++) if (x != y) { int da = ((y / 3 - x / 3) % 3 + 3) % 3, db = ((y % 3 - x % 3) % 3 + 3) % 3; if (conn[da][db]) { OUT[x] |= 1u << y; IN[y] |= 1u << x; } }
    uint32_t full = (1u << N) - 1; dp(); arcs_sub(full);
    int nodd = 0; for (int u = 0; u < N; u++) { uint32_t t = OUT[u]; while (t) { int v = __builtin_ctz(t); t &= t - 1; if (C[u][v] & 1) nodd++; } }
    tot++; if (nodd == 0) alleven++;
  }
  printf("Cayley Z3xZ3: %ld tournaments, all-even: %ld\n", tot, alleven);
}

int main(int argc, char **argv) {
  if (argc < 2) return 1;
  if (!strcmp(argv[1], "labeled")) labeled(atoi(argv[2]));
  else if (!strcmp(argv[1], "classes")) classes(argc > 2 && !strcmp(argv[2], "beta"), argc > 3 && !strcmp(argv[3], "shave"));
  else if (!strcmp(argv[1], "one")) { set_from_string(argv[2]); one(); }
  else if (!strcmp(argv[1], "par2")) { set_from_string(argv[2]); par2(); }
  else if (!strcmp(argv[1], "paley")) paley(atoi(argv[2]), atoi(argv[3]));
  else if (!strcmp(argv[1], "circ")) circ(atoi(argv[2]));
  else if (!strcmp(argv[1], "z3z3")) z3z3();
  else if (!strcmp(argv[1], "beta")) { set_from_string(argv[2]); long sg; long hb = hallval(&sg); int b = betaval(N); printf("beta %s: beta=%d hall=%ld sigma=%ld\n", argv[2], b, hb, sg); }
  printf("fails=%ld\n", fails);
  return fails ? 2 : 0;
}
