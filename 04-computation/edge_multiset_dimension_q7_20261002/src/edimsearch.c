// edimsearch.c -- exhaustive search for edge-multiset resolving k-subsets of Q_d (2 <= d <= 7) up to
// Aut(Q_d), in the MINIMUM-IMBALANCE normal form ("normal form B", generalised to any d).
// (scratch lane edim2/q7, 2026-10-01)
//
// NORMAL FORM (proof in q7_report.md).  beta_i(S) = #{s in S: s_i = 0} - #{s in S: s_i = 1}.  Every
// k-set S is Aut(Q_d)-equivalent to  S = A x {0}  u  B x {1}  (split coordinate = d-1) with
//   |A| = a >= b = |B|,  a - b = min_i |beta_i(S)|, hence |beta_i(S)| >= a - b for i = 0..d-2,
//   A = the canonical representative of its Aut(Q_{d-1})-orbit (read from the input file),
//   B = any b-subset of Q_{d-1} satisfying the imbalance constraint.
// For one value of a and a range of representatives A the program enumerates ALL b-subsets B (DFS in
// increasing vertex order; a subtree is skipped only if NO completion can satisfy the final imbalance
// constraint) and tests every resulting k-set S.  The leaf count is the number of pairs (A, B) in the
// normal-form domain; it is checked against an independent dynamic-programming count.
//
// EDGE KEYS.  Edge e = {u, u + e_i} (u_i = 0), vertex s: d(e,s) = popcount((u^s) & ~(1<<i)) (projection
// lemma).  key(e) = sum_{s in S} C(d(e,s)) with C(t) = 32^t for t <= d-2 and C(d-1) = 0: the histogram
// levels 0..d-2 in 5-bit fields; level d-1 is implied because every histogram sums to k.  For k <= 31
// the key is an injective function of the histogram, so edges collide iff their keys are equal.
//
// LEAF TEST.  The last B-level is handled in bulk.  At a node where one element v remains, the admissible
// v form a 64-bit mask Cm.  Edge pair (e,f) kills v iff key(e) + C(d(e,v)) = key(f) + C(d(f,v)), i.e.
// C(d(f,v)) - C(d(e,v)) = de := key(e) - key(f).  Because (x,y) -> C(y) - C(x) is injective on x != y
// (base-32 digits), the killed set is {v: d(e,v) = d(f,v)} if de = 0, the class {v: d(e,v)=x, d(f,v)=y}
// if de = C(y) - C(x), and empty otherwise (exact).  Two engines:
//   engine 0 (default, "bulk"): a hot list of edge pairs is evaluated in branchless blocks of 8 (kill
//            masks via the difference decoder) until every candidate is killed; the completing pair is
//            moved to the front.
//   engine 1 ("per-leaf"): for each candidate v separately, the hot list is scanned with the direct test
//            key(e) + C(d(e,v)) == key(f) + C(d(f,v)) until the first hit (move-to-front).
// In both engines a candidate that survives the hot list gets a FULL TEST of all edge keys (hash table);
// the first colliding pair found is inserted at the front of the hot list (and asserted to kill v);
// a candidate without collision is a resolving set; it is re-verified from scratch and printed.  Hence
// every admissible leaf is either killed by an exactly computed collision or fully tested.
// Option -v (verification): every leaf is additionally tested from scratch (keys recomputed from the
// definition d(e,s) = min over endpoints) and the program aborts if any killed leaf has no collision,
// if any kill mask differs from a brute-force recomputation, or if incremental keys differ from scratch.
//
// usage: edimsearch d k a repfile first last [-e engine] [-H hotcap] [-v]
//        edimsearch d 0 -x        (check mode: read vertex lists from stdin, print RESOLVING/NOT + defect)
// repfile: binary little-endian uint64 masks (bit v = vertex v of Q_{d-1}); reps with index in [first,last).
// Output: "RESOLVING d= k= set=v1,v2,..." per resolving leaf; final "SUMMARY ..." line.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#define MAXE 448
#define HOTMAX 4096
static int D, N, N1, E, K, A_, B_, DELTA;
static uint32_t C[128][MAXE] __attribute__((aligned(64)));   // C[s][e] = C(d(e,s))
static uint32_t ZERO[MAXE] __attribute__((aligned(64)));
static uint64_t DM[MAXE][8];            // DM[e][x] = { w in 0..N1-1 : d(e, N1 + w) = x };  DM[e][7] = 0
static int EU[MAXE], EI[MAXE];
static uint64_t ZI[7], OI[7], FULL1;   // ZI[i] = {w : w_i = 0}, OI[i] = {w : w_i = 1}  (layer-1 indices)
static uint32_t Cval[8];
typedef struct { int32_t k; uint8_t x, y, pad0, pad1; } dent_t;
static dent_t dtab[256]; static uint32_t dmagic;               // perfect hash of {0} u {C(y)-C(x)}
typedef struct { uint16_t e, f; uint32_t pad; uint64_t eq; } hot_t;
static hot_t hot[HOTMAX + 16]; static int nhot = 0, hotcap = 1024, engine = 0, verify = 0;

static uint32_t Kst[8][MAXE] __attribute__((aligned(64)));
static int beta[8], chosen[8];
static uint64_t curA;
static long long st_leaves = 0, st_nodes = 0, st_full = 0, st_found = 0, st_scan = 0, st_vleaves = 0, st_reps = 0, st_repsfeas = 0;

static inline int dist_e(int e, int s) { return __builtin_popcount((EU[e] ^ s) & ~(1 << EI[e])); }
static void fatal(const char *m) { fprintf(stderr, "FATAL: %s\n", m); printf("FATAL: %s\n", m); fflush(stdout); exit(9); }

static void setup(void) {
  N = 1 << D; N1 = N >> 1; E = 0;
  for (int u = 0; u < N; u++) for (int i = 0; i < D; i++) if (!(u >> i & 1)) { EU[E] = u; EI[E] = i; E++; }
  if (E != D * N1) fatal("edge count");
  for (int t = 0; t < D; t++) Cval[t] = (t <= D - 2) ? (1u << (5 * t)) : 0u;
  for (int s = 0; s < N; s++) for (int e = 0; e < E; e++) C[s][e] = Cval[dist_e(e, s)];
  memset(ZERO, 0, sizeof ZERO);
  for (int e = 0; e < E; e++) { for (int x = 0; x < 8; x++) DM[e][x] = 0;
    for (int w = 0; w < N1; w++) DM[e][dist_e(e, N1 + w)] |= 1ULL << w; }
  FULL1 = (N1 == 64) ? ~0ULL : ((1ULL << N1) - 1);
  for (int i = 0; i < D - 1; i++) { ZI[i] = OI[i] = 0; for (int w = 0; w < N1; w++) { if (w >> i & 1) OI[i] |= 1ULL << w; else ZI[i] |= 1ULL << w; } }
  int32_t vals[64]; int8_t vx[64], vy[64]; int nv = 0;
  vals[nv] = 0; vx[nv] = 7; vy[nv] = 7; nv++;
  for (int x = 0; x < D; x++) for (int y = 0; y < D; y++) if (x != y) { vals[nv] = (int32_t)(Cval[y] - Cval[x]); vx[nv] = x; vy[nv] = y; nv++; }
  for (int a = 0; a < nv; a++) for (int b = a + 1; b < nv; b++) if (vals[a] == vals[b]) fatal("difference map not injective");
  uint64_t z = 0x243F6A8885A308D3ULL; int ok = 0;
  for (int tries = 0; tries < 1000000 && !ok; tries++) {
    z ^= z << 13; z ^= z >> 7; z ^= z << 17; dmagic = (uint32_t)z | 1u;
    uint8_t used[256] = {0}; ok = 1;
    for (int a = 0; a < nv && ok; a++) { uint32_t h = ((uint32_t)vals[a] * dmagic) >> 24; if (used[h]) ok = 0; used[h] = 1; }
  }
  if (!ok) fatal("no perfect hash");
  for (int h = 0; h < 256; h++) { dtab[h].k = INT32_MIN; dtab[h].x = 7; dtab[h].y = 7; }
  for (int a = 0; a < nv; a++) { uint32_t h = ((uint32_t)vals[a] * dmagic) >> 24; dtab[h].k = vals[a]; dtab[h].x = vx[a]; dtab[h].y = vy[a]; }
  for (int t = 0; t < HOTMAX + 16; t++) { hot[t].e = hot[t].f = 0; hot[t].eq = 0; }   // null pairs kill nothing
}
static inline uint64_t eqmask(int e, int f) { uint64_t m = 0; for (int x = 0; x < D; x++) m |= DM[e][x] & DM[f][x]; return m; }
// kill set of pair p at the node whose keys are Kb + cu (branchless; null pair e = f = 0, eq = 0 gives 0)
static inline uint64_t killmask(const uint32_t *Kb, const uint32_t *cu, const hot_t *p) {
  int32_t de = (int32_t)((Kb[p->e] + cu[p->e]) - (Kb[p->f] + cu[p->f]));
  const dent_t *t = &dtab[((uint32_t)de * dmagic) >> 24];
  uint64_t sel = -(uint64_t)(t->k == de), z = -(uint64_t)(de == 0);
  uint64_t m = DM[p->e][t->x] & DM[p->f][t->y] & sel;
  return (m & ~z) | (p->eq & z);
}
// reference kill set (slow, direct), used by -v
static uint64_t killmask_ref(const uint32_t *Kb, const uint32_t *cu, int e, int f, uint64_t Cm) {
  uint64_t m = 0;
  for (uint64_t c = Cm; c; c &= c - 1) { int v = __builtin_ctzll(c);
    if (Kb[e] + cu[e] + C[N1 + v][e] == Kb[f] + cu[f] + C[N1 + v][f]) m |= 1ULL << v; }
  return m;
}
// full test: returns 1 and a colliding pair (pe, pf) if two edge keys are equal, else 0
static uint32_t htk[2048], hts[2048]; static uint16_t hte[2048]; static uint32_t stampgen = 0;
static int fullcheck(const uint32_t *Kb, const uint32_t *cu, const uint32_t *cv, int *pe, int *pf) {
  if (++stampgen == 0) { memset(hts, 0, sizeof hts); stampgen = 1; }
  for (int e = 0; e < E; e++) {
    uint32_t k = Kb[e] + cu[e] + cv[e], h = (k * 2654435761u) >> 21;
    while (hts[h] == stampgen) { if (htk[h] == k) { *pe = hte[h]; *pf = e; return 1; } h = (h + 1) & 2047; }
    hts[h] = stampgen; htk[h] = k; hte[h] = (uint16_t)e;
  }
  return 0;
}
// keys of a set computed from scratch, straight from the definition d(e,s) = min(d(u,s), d(v,s))
static void keys_from_scratch(const int *S, int k, uint32_t *out) {
  for (int e = 0; e < E; e++) { uint32_t x = 0; int u = EU[e], v = EU[e] | (1 << EI[e]);
    for (int j = 0; j < k; j++) { int du = __builtin_popcount(u ^ S[j]), dv = __builtin_popcount(v ^ S[j]); x += Cval[du < dv ? du : dv]; }
    out[e] = x; }
}
static int current_set(int u, int v, int *S) {   // A u {N1+chosen[0..B_-3]} u {N1+u} u {N1+v}  (u,v < 0: absent)
  int k = 0;
  for (uint64_t m = curA; m; m &= m - 1) S[k++] = __builtin_ctzll(m);
  for (int j = 0; j < B_ - 2; j++) S[k++] = N1 + chosen[j];
  if (u >= 0) S[k++] = N1 + u;
  if (v >= 0) S[k++] = N1 + v;
  return k;
}
static void report(int u, int v) {
  int S[64]; int k = current_set(u, v, S);
  if (k != K) fatal("report size");
  uint32_t Kv[MAXE]; keys_from_scratch(S, k, Kv); int pe, pf;   // re-verify from scratch before printing
  if (fullcheck(Kv, ZERO, ZERO, &pe, &pf)) fatal("reported set is not resolving when recomputed from scratch");
  for (int i = 1; i < k; i++) for (int j = i; j > 0 && S[j - 1] > S[j]; j--) { int t = S[j]; S[j] = S[j - 1]; S[j - 1] = t; }
  printf("RESOLVING d=%d k=%d set=", D, K);
  for (int i = 0; i < k; i++) printf(i ? ",%d" : "%d", S[i]);
  printf("\n"); fflush(stdout); st_found++;
}
static void hot_insert(int e, int f) {   // at the front; hot[nhot..nhot+7] remain null pairs
  if (nhot < hotcap) nhot++;
  memmove(hot + 1, hot, sizeof(hot_t) * (nhot - 1));
  hot[0].e = (uint16_t)e; hot[0].f = (uint16_t)f; hot[0].eq = eqmask(e, f);
  for (int t = 0; t < 8; t++) { hot[nhot + t].e = 0; hot[nhot + t].f = 0; hot[nhot + t].eq = 0; }
}
static inline void mtf(int i) { if (i > 0) { hot_t t = hot[i]; memmove(hot + 1, hot, sizeof(hot_t) * i); hot[0] = t; } }

// node with one remaining element; keys of the current (k-1)-set are Kb + cu; Cm = admissible last elements
static void last_level(const uint32_t *Kb, const uint32_t *cu, uint64_t Cm, int u) {
  st_leaves += __builtin_popcountll(Cm); st_nodes++;
  uint64_t truth = 0, rem = Cm;
  if (verify) {   // exact status of every leaf, independent of the kill logic
    int S[64]; int k = current_set(u, -1, S); uint32_t Kv[MAXE]; keys_from_scratch(S, k, Kv);
    for (int e = 0; e < E; e++) if (Kv[e] != Kb[e] + cu[e]) fatal("incremental keys differ from scratch keys");
    for (uint64_t m = Cm; m; m &= m - 1) { int v = __builtin_ctzll(m), pe, pf;
      int S2[64]; int k2 = current_set(u, v, S2); uint32_t K2[MAXE]; keys_from_scratch(S2, k2, K2);
      if (fullcheck(K2, ZERO, ZERO, &pe, &pf)) truth |= 1ULL << v; }
    st_vleaves += __builtin_popcountll(Cm);
    for (int j = 0; j < nhot; j++) if ((killmask(Kb, cu, &hot[j]) & Cm) != killmask_ref(Kb, cu, hot[j].e, hot[j].f, Cm)) fatal("killmask != reference");
  }
  if (engine == 0) {   // bulk, branchless blocks of 8
    int jj;
    for (jj = 0; jj < nhot; jj += 8) {
      uint64_t km[8], all = 0;
      for (int t = 0; t < 8; t++) { km[t] = killmask(Kb, cu, &hot[jj + t]); all |= km[t]; }
      if (!(rem & ~all)) {
        int t = 0; for (; t < 8; t++) { rem &= ~km[t]; if (!rem) break; }
        st_scan += jj + t + 1; mtf(jj + t);
        break;
      }
      rem &= ~all;
    }
    if (jj >= nhot) st_scan += nhot;
  } else {             // per-leaf: first hit with the direct key comparison, move-to-front
    uint64_t r2 = 0;
    for (uint64_t m = Cm; m; m &= m - 1) {
      int v = __builtin_ctzll(m), hit = -1; const uint32_t *cv = C[N1 + v];
      for (int i = 0; i < nhot; i++) { int e = hot[i].e, f = hot[i].f; if (Kb[e] + cu[e] + cv[e] == Kb[f] + cu[f] + cv[f]) { hit = i; break; } }
      st_scan += hit < 0 ? nhot : hit + 1;
      if (hit < 0) r2 |= 1ULL << v; else mtf(hit);
    }
    rem = r2;
  }
  if (verify && ((Cm & ~rem) & ~truth)) fatal("hot list killed a leaf without collision");
  while (rem) {   // full tests of the survivors
    int v = __builtin_ctzll(rem), pe, pf; st_full++;
    if (!fullcheck(Kb, cu, C[N1 + v], &pe, &pf)) {
      if (verify && (truth >> v & 1)) fatal("verify mismatch (full test says resolving)");
      report(u, v); rem &= rem - 1; continue;
    }
    if (verify && !(truth >> v & 1)) fatal("verify mismatch (full test says collision)");
    hot_insert(pe, pf);
    uint64_t km = killmask(Kb, cu, &hot[0]);
    if (!(km >> v & 1)) fatal("kill logic inconsistent with full test");
    if (verify && (km & rem & ~truth)) fatal("new pair killed a leaf without collision");
    rem &= ~km;
  }
}
// imbalance (coordinates 0..D-2 of A plus the chosen B elements)
static inline int feasible(const int *bt, int r) {   // some completion with r more elements can reach |beta_i| >= DELTA
  for (int i = 0; i < D - 1; i++) if (bt[i] + r < DELTA && bt[i] - r > -DELTA) return 0;
  return 1;
}
static inline uint64_t allowed_last(const int *bt) {   // last elements w giving |beta_i + sigma_i(w)| >= DELTA for all i
  if (DELTA <= 1) return FULL1;                       // vacuous: beta_i(S) has the parity of k, and DELTA = k mod 2
  uint64_t m = FULL1;
  for (int i = 0; i < D - 1; i++) {
    uint64_t mi = 0;
    if (bt[i] + 1 >= DELTA || bt[i] + 1 <= -DELTA) mi |= ZI[i];
    if (bt[i] - 1 >= DELTA || bt[i] - 1 <= -DELTA) mi |= OI[i];
    m &= mi;
  }
  return m;
}
static inline void addsig(int *bt, int w, int sgn) { for (int i = 0; i < D - 1; i++) bt[i] += ((w >> i) & 1) ? -sgn : sgn; }
static inline uint64_t above(int u) { return (u + 1 >= 64) ? 0 : (FULL1 & ~((2ULL << u) - 1)); }

static void dfs(int depth, int start, const uint32_t *Kb) {
  int r = B_ - depth;
  if (r == 1) {   // only reached when B_ == 1
    uint64_t Cm = (start >= N1) ? 0 : (FULL1 & ~((1ULL << start) - 1));
    Cm &= allowed_last(beta);
    if (Cm) last_level(Kb, ZERO, Cm, -1);
    return;
  }
  if (r == 2) {
    for (int u = start; u <= N1 - 2; u++) {
      addsig(beta, u, 1);
      uint64_t Cm = above(u) & allowed_last(beta);
      if (Cm) { chosen[depth] = u; last_level(Kb, C[N1 + u], Cm, u); }
      addsig(beta, u, -1);
    }
    return;
  }
  uint32_t *Kn = Kst[depth + 1];
  for (int w = start; w <= N1 - r; w++) {
    addsig(beta, w, 1);
    if (feasible(beta, r - 1)) {
      const uint32_t *cw = C[N1 + w];
      for (int e = 0; e < E; e++) Kn[e] = Kb[e] + cw[e];
      chosen[depth] = w;
      dfs(depth + 1, w + 1, Kn);
    }
    addsig(beta, w, -1);
  }
}
static int check_mode(void) {
  char line[8192];
  while (fgets(line, sizeof line, stdin)) {
    int S[256], k = 0; char *p = line;
    while (*p) { while (*p && (*p < '0' || *p > '9')) p++; if (!*p) break; S[k++] = (int)strtol(p, &p, 10); }
    if (k == 0) continue;
    if (k > 31) { printf("SKIP k=%d>31\n", k); continue; }
    for (int j = 0; j < k; j++) if (S[j] < 0 || S[j] >= N) fatal("vertex out of range");
    uint32_t Kv[MAXE]; keys_from_scratch(S, k, Kv);
    uint32_t t[MAXE]; memcpy(t, Kv, sizeof(uint32_t) * E);
    for (int i = 1; i < E; i++) { uint32_t x = t[i]; int j = i; while (j > 0 && t[j - 1] > x) { t[j] = t[j - 1]; j--; } t[j] = x; }
    int distinct = 1; for (int i = 1; i < E; i++) distinct += t[i] != t[i - 1];
    int pe, pf; int coll = fullcheck(Kv, ZERO, ZERO, &pe, &pf);
    if ((coll == 0) != (distinct == E)) fatal("check mode: hash test and sort test disagree");
    printf("%s k=%d distinct=%d defect=%d\n", distinct == E ? "RESOLVING" : "NOT", k, distinct, E - distinct);
  }
  return 0;
}
int main(int argc, char **argv) {
  if (argc >= 4 && !strcmp(argv[3], "-x")) { D = atoi(argv[1]); if (D < 2 || D > 7) return 1; setup(); return check_mode(); }
  if (argc < 7) { fprintf(stderr, "usage: edimsearch d k a repfile first last [-e engine] [-H hotcap] [-v]\n"); return 1; }
  D = atoi(argv[1]); K = atoi(argv[2]); A_ = atoi(argv[3]); B_ = K - A_; DELTA = A_ - B_;
  const char *rf = argv[4]; long long first = atoll(argv[5]), last = atoll(argv[6]);
  for (int i = 7; i < argc; i++) {
    if (!strcmp(argv[i], "-v")) verify = 1;
    else if (!strcmp(argv[i], "-H") && i + 1 < argc) hotcap = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-e") && i + 1 < argc) engine = atoi(argv[++i]);
    else { fprintf(stderr, "bad option %s\n", argv[i]); return 1; }
  }
  if (D < 2 || D > 7 || K < 1 || K > 31 || DELTA < 0 || B_ < 0 || B_ > 7 || hotcap < 1 || hotcap > HOTMAX || engine < 0 || engine > 1) { fprintf(stderr, "bad parameters\n"); return 1; }
  setup();
  FILE *f = fopen(rf, "rb"); if (!f) { perror(rf); return 3; }
  fseek(f, 0, SEEK_END); long long nrep = ftell(f) / 8;
  if (last > nrep) last = nrep;
  if (first < 0) first = 0;
  uint64_t *reps = malloc(8 * (size_t)(last > first ? last - first : 1));
  if (last > first) { fseek(f, first * 8, SEEK_SET); if (fread(reps, 8, (size_t)(last - first), f) != (size_t)(last - first)) { perror("read"); return 3; } }
  fclose(f);
  struct timespec t0, t1; clock_gettime(CLOCK_MONOTONIC, &t0);
  for (long long r = first; r < last; r++) {
    uint64_t M = reps[r - first];
    if (__builtin_popcountll(M) != A_ || (N1 < 64 && (M >> N1))) fatal("bad representative");
    st_reps++; curA = M;
    for (int i = 0; i < D - 1; i++) beta[i] = 0;
    for (uint64_t m = M; m; m &= m - 1) addsig(beta, __builtin_ctzll(m), 1);
    if (!feasible(beta, B_)) continue;
    st_repsfeas++;
    uint32_t *K0 = Kst[0]; memset(K0, 0, sizeof(uint32_t) * E);
    for (uint64_t m = M; m; m &= m - 1) { const uint32_t *cs = C[__builtin_ctzll(m)]; for (int e = 0; e < E; e++) K0[e] += cs[e]; }
    if (B_ == 0) {   // S = A x {0}: exactly one leaf (feasible(beta, 0) <=> all |beta_i| >= DELTA)
      st_leaves++; int pe, pf;
      if (!fullcheck(K0, ZERO, ZERO, &pe, &pf)) report(-1, -1);
      continue;
    }
    dfs(0, 0, K0);
  }
  clock_gettime(CLOCK_MONOTONIC, &t1);
  double sec = (t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec);
  printf("SUMMARY d=%d k=%d a=%d b=%d first=%lld last=%lld reps=%lld feasible_reps=%lld leaves=%lld nodes=%lld fullchecks=%lld found=%lld scan=%lld vleaves=%lld engine=%d hotcap=%d verify=%d sec=%.3f ns_per_leaf=%.3f\n",
         D, K, A_, B_, first, last, st_reps, st_repsfeas, st_leaves, st_nodes, st_full, st_found, st_scan, st_vleaves, engine, hotcap, verify, sec, st_leaves ? 1e9 * sec / st_leaves : 0.0);
  return 0;
}
