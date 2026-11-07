/* aud_enum.c - independent auditor's enumerator for THM-4922 / HYP-9512 (first creators of H values).
 *
 * Written from scratch by the auditor; it shares no code with classgen.c / order10.c.
 * Method (different from the packet's global hash set of minimal codes):
 *   - canonical labelling by individualization-refinement, with an initial colouring by (score, cyclic triangles
 *     through v), refinement by IN-neighbour counts per cell, individualization in the SMALLEST non-singleton cell,
 *     and the lexicographically MAXIMAL leaf code (first pair = most significant bit);
 *   - all leaves are explored (no pruning), so the leaves with the best code are exactly an Aut(T)-coset; this gives
 *     |Aut(T)| and the Aut(T)-orbits (union-find over the automorphisms best_leaf -> leaf);
 *   - generation by canonical augmentation (McKay): an extension R + v of an order-(n-1) representative R is accepted
 *     iff the new vertex lies in the Aut-orbit of the vertex in the last canonical position; then each class arises
 *     from exactly one representative, and duplicates can only occur among the extensions of the same R (removed
 *     locally). No global hash set is used.
 * Checks built in: sum over classes of n!/|Aut| = 2^C(n,2) (orbit-stabilizer), and (mode "brute") a comparison of the
 * canonical form with a brute-force canonical form (max over all n! labellings) for n <= 7.
 *
 * usage:  aud_enum NMAX OUTDIR            (NMAX <= 10)
 *         aud_enum brute 7                (self-test of the canonical form against brute force, orders 2..7)
 * Output, n <= 9: OUTDIR/aud_cls_<n>.txt, one line per class:
 *         code H c3 strong SC aut hc | H(T-v_0) ... H(T-v_{n-1})
 * Output, n = 10: OUTDIR/aud_o10_byH.txt (per H: all classes, SC classes; first-creator fiber, SC first creators,
 *         fully critical, non-strong first creators, min k, min c3 and count at min c3, min hc and count at min hc)
 *         and a summary on stdout.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define MAXN 10
#define HCAP 32768
typedef unsigned long long u64;
typedef unsigned short m16;

/* ------------------------------------------------------------------ canonical labelling */
typedef struct {
  int n;
  m16 out[MAXN];
  int have;
  u64 best;
  int lab0[MAXN];
  long nbest;
  int uf[MAXN];
} cctx;

static int uf_find(int *uf, int x) { while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; } return x; }
static void uf_union(int *uf, int a, int b) { a = uf_find(uf, a); b = uf_find(uf, b); if (a != b) uf[a] = b; }

static u64 leaf_code(const cctx *C, const int *lab) {
  u64 code = 0; int n = C->n;
  for (int i = 0; i < n; i++) {
    m16 oi = C->out[lab[i]];
    for (int j = i + 1; j < n; j++) code = (code << 1) | (u64)((oi >> lab[j]) & 1);
  }
  return code;
}

/* map arbitrary ordered colour values to 0..k-1 preserving order; returns k */
static int normalize(int *col, int n) {
  int vals[MAXN], nv = 0;
  for (int v = 0; v < n; v++) {
    int x = col[v], seen = 0;
    for (int i = 0; i < nv; i++) if (vals[i] == x) { seen = 1; break; }
    if (!seen) vals[nv++] = x;
  }
  for (int i = 1; i < nv; i++) { int x = vals[i], j = i - 1; while (j >= 0 && vals[j] > x) { vals[j + 1] = vals[j]; j--; } vals[j + 1] = x; }
  for (int v = 0; v < n; v++) for (int i = 0; i < nv; i++) if (vals[i] == col[v]) { col[v] = i; break; }
  return nv;
}

static int sigcmp(const int *a, const int *b, int len) {
  for (int i = 0; i < len; i++) if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
  return 0;
}

/* equitable refinement by in-neighbour counts; col is normalized on return; returns number of cells */
static int refine(const cctx *C, int *col) {
  int n = C->n;
  int k = normalize(col, n);
  for (;;) {
    int sig[MAXN][MAXN + 1];
    for (int v = 0; v < n; v++) { sig[v][0] = col[v]; for (int c = 0; c < k; c++) sig[v][1 + c] = 0; }
    for (int u = 0; u < n; u++) {
      m16 o = C->out[u];
      while (o) { int v = __builtin_ctz(o); o &= o - 1; sig[v][1 + col[u]]++; }   /* u -> v : u is an in-neighbour of v */
    }
    int ord[MAXN];
    for (int v = 0; v < n; v++) ord[v] = v;
    for (int i = 1; i < n; i++) { int x = ord[i], j = i - 1; while (j >= 0 && sigcmp(sig[ord[j]], sig[x], k + 1) > 0) { ord[j + 1] = ord[j]; j--; } ord[j + 1] = x; }
    int nc[MAXN], c = 0; nc[ord[0]] = 0;
    for (int i = 1; i < n; i++) { if (sigcmp(sig[ord[i]], sig[ord[i - 1]], k + 1) != 0) c++; nc[ord[i]] = c; }
    for (int v = 0; v < n; v++) col[v] = nc[v];
    if (c + 1 == k) return k;
    k = c + 1;
  }
}

static void search(cctx *C, int *col, int k) {
  int n = C->n;
  if (k == n) {
    int lab[MAXN];
    for (int v = 0; v < n; v++) lab[col[v]] = v;
    u64 code = leaf_code(C, lab);
    if (!C->have || code > C->best) {
      C->have = 1; C->best = code; memcpy(C->lab0, lab, sizeof(int) * n); C->nbest = 1;
      for (int i = 0; i < n; i++) C->uf[i] = i;
    } else if (code == C->best) {
      C->nbest++;
      for (int i = 0; i < n; i++) uf_union(C->uf, C->lab0[i], lab[i]);
    }
    return;
  }
  int size[MAXN] = {0};
  for (int v = 0; v < n; v++) size[col[v]]++;
  int t = -1;
  for (int c = 0; c < k; c++) if (size[c] > 1 && (t < 0 || size[c] < size[t])) t = c;
  for (int v = 0; v < n; v++) if (col[v] == t) {
    int c2[MAXN];
    for (int u = 0; u < n; u++) c2[u] = 2 * col[u] + ((col[u] == t && u != v) ? 1 : 0);
    int k2 = refine(C, c2);
    search(C, c2, k2);
  }
}

static void canon(cctx *C, const m16 *out, int n) {
  C->n = n; memcpy(C->out, out, sizeof(m16) * n); C->have = 0; C->nbest = 0;
  int col[MAXN];
  for (int v = 0; v < n; v++) {
    int s = __builtin_popcount(out[v]);
    m16 in = 0; for (int u = 0; u < n; u++) if ((out[u] >> v) & 1) in |= (m16)(1u << u);
    int t3 = 0; m16 o = out[v];
    while (o) { int u = __builtin_ctz(o); o &= o - 1; t3 += __builtin_popcount(out[u] & in); }  /* v->u->w->v */
    col[v] = s * 100 + t3;
  }
  int k = refine(C, col);
  search(C, col, k);
}

static void decode(u64 code, int n, m16 *out) {
  memset(out, 0, sizeof(m16) * MAXN);
  int b = n * (n - 1) / 2 - 1;
  for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
    if ((code >> b) & 1) out[i] |= (m16)(1u << j); else out[j] |= (m16)(1u << i);
    b--;
  }
}

static void converse(const m16 *out, int n, m16 *op) {
  m16 full = (m16)((1u << n) - 1);
  for (int v = 0; v < n; v++) op[v] = (m16)(full & ~out[v] & ~(1u << v));
}

/* ------------------------------------------------------------------ invariants */
typedef struct { u64 f[1 << MAXN][MAXN]; } dpbuf;

/* H(T) and dh[x] = H(T - x) */
static u64 hpaths(dpbuf *D, const m16 *out, int n, u64 *dh) {
  int F = 1 << n;
  for (int S = 0; S < F; S++) for (int v = 0; v < n; v++) D->f[S][v] = 0;
  for (int v = 0; v < n; v++) D->f[1 << v][v] = 1;
  for (int S = 1; S < F; S++) for (int v = 0; v < n; v++) {
    u64 x = D->f[S][v]; if (!x) continue;
    m16 nx = (m16)(out[v] & ~S);
    while (nx) { int u = __builtin_ctz(nx); nx &= nx - 1; D->f[S | (1 << u)][u] += x; }
  }
  u64 h = 0; for (int v = 0; v < n; v++) h += D->f[F - 1][v];
  if (dh) for (int x = 0; x < n; x++) {
    int S = (F - 1) & ~(1 << x); u64 s = 0;
    if (n == 1) { dh[x] = 1; continue; }        /* empty tournament: one (empty) path by convention; unused */
    for (int v = 0; v < n; v++) s += D->f[S][v];
    dh[x] = s;
  }
  return h;
}

/* number of Hamiltonian cycles (directed), n >= 3; 0 for n < 3 */
static u64 hcycles(dpbuf *D, const m16 *out, int n) {
  if (n < 3) return 0;
  int F = 1 << n;
  for (int S = 0; S < F; S++) for (int v = 0; v < n; v++) D->f[S][v] = 0;
  D->f[1][0] = 1;
  for (int S = 1; S < F; S += 2) for (int v = 0; v < n; v++) {
    u64 x = D->f[S][v]; if (!x) continue;
    m16 nx = (m16)(out[v] & ~S);
    while (nx) { int u = __builtin_ctz(nx); nx &= nx - 1; D->f[S | (1 << u)][u] += x; }
  }
  u64 c = 0; for (int v = 1; v < n; v++) if (out[v] & 1) c += D->f[F - 1][v];
  return c;
}

static int c3count(const m16 *out, int n) {
  int s = n * (n - 1) * (n - 2) / 6;
  for (int v = 0; v < n; v++) { int d = __builtin_popcount(out[v]); s -= d * (d - 1) / 2; }
  return s;
}

static int is_strong(const m16 *out, int n) {
  if (n <= 1) return 1;
  m16 full = (m16)((1u << n) - 1), r = 1, prev = 0;
  while (r != prev) { prev = r; for (int v = 0; v < n; v++) if ((r >> v) & 1) r |= out[v]; }
  if (r != full) return 0;
  m16 in[MAXN] = {0};
  for (int v = 0; v < n; v++) for (int u = 0; u < n; u++) if ((out[v] >> u) & 1) in[u] |= (m16)(1u << v);
  r = 1; prev = 0;
  while (r != prev) { prev = r; for (int v = 0; v < n; v++) if ((r >> v) & 1) r |= in[v]; }
  return r == full;
}

/* ------------------------------------------------------------------ generation by canonical augmentation */
static int cmp_u64(const void *a, const void *b) { u64 x = *(const u64 *)a, y = *(const u64 *)b; return x < y ? -1 : x > y; }

/* all classes of order n generated from representative code rc of order n-1; returns count, codes in res */
static int extend_rep(u64 rc, int n, u64 *res) {
  m16 base[MAXN]; decode(rc, n - 1, base);
  cctx C; int cnt = 0;
  for (int P = 0; P < (1 << (n - 1)); P++) {
    m16 t[MAXN];
    for (int v = 0; v < n - 1; v++) t[v] = base[v];
    t[n - 1] = 0;
    for (int v = 0; v < n - 1; v++) { if ((P >> v) & 1) t[n - 1] |= (m16)(1u << v); else t[v] |= (m16)(1u << (n - 1)); }
    canon(&C, t, n);
    if (uf_find(C.uf, n - 1) != uf_find(C.uf, C.lab0[n - 1])) continue;   /* new vertex not in the canonical orbit */
    res[cnt++] = C.best;
  }
  qsort(res, cnt, sizeof(u64), cmp_u64);
  int m = 0;
  for (int i = 0; i < cnt; i++) if (i == 0 || res[i] != res[i - 1]) res[m++] = res[i];
  return m;
}

static u64 factorial(int n) { u64 f = 1; for (int i = 2; i <= n; i++) f *= i; return f; }

/* ------------------------------------------------------------------ brute-force self-test */
static u64 brute_code(const m16 *out, int n, long *aut) {
  int perm[MAXN]; for (int i = 0; i < n; i++) perm[i] = i;
  cctx C; C.n = n; memcpy(C.out, out, sizeof(m16) * n);
  u64 best = 0; long cnt = 0; int first = 1;
  /* Heap's algorithm */
  int c[MAXN] = {0};
  for (;;) {
    u64 code = leaf_code(&C, perm);
    if (first || code > best) { best = code; cnt = 1; first = 0; } else if (code == best) cnt++;
    int i = 1;
    while (i < n && c[i] >= i) { c[i] = 0; i++; }
    if (i >= n) break;
    if (i % 2 == 0) { int tmp = perm[0]; perm[0] = perm[i]; perm[i] = tmp; } else { int tmp = perm[c[i]]; perm[c[i]] = perm[i]; perm[i] = tmp; }
    c[i]++;
  }
  if (aut) *aut = cnt;
  return best;
}

static int brute_selftest(int nmax) {
  /* order-by-order: reps from canonical augmentation; for every extension compare canonical forms */
  u64 *reps = malloc(sizeof(u64)); size_t nreps = 1; reps[0] = 0;
  int ok = 1;
  for (int n = 2; n <= nmax; n++) {
    size_t cap = nreps << (n - 1);
    u64 *mine = malloc(cap * sizeof(u64)), *br = malloc(cap * sizeof(u64)); size_t ne = 0;
    long badaut = 0;
    for (size_t r = 0; r < nreps; r++) {
      m16 base[MAXN]; decode(reps[r], n - 1, base);
      for (int P = 0; P < (1 << (n - 1)); P++) {
        m16 t[MAXN];
        for (int v = 0; v < n - 1; v++) t[v] = base[v];
        t[n - 1] = 0;
        for (int v = 0; v < n - 1; v++) { if ((P >> v) & 1) t[n - 1] |= (m16)(1u << v); else t[v] |= (m16)(1u << (n - 1)); }
        cctx C; canon(&C, t, n); long ab;
        u64 b = brute_code(t, n, &ab);
        if (ab != C.nbest) badaut++;
        /* the canonical code must be an isomorphic copy: brute code of decode(C.best) must equal b */
        m16 d[MAXN]; decode(C.best, n, d);
        if (brute_code(d, n, NULL) != b) { ok = 0; }
        mine[ne] = C.best; br[ne] = b; ne++;
      }
    }
    /* bijection between distinct mine-codes and distinct brute codes: count distinct pairs */
    u64 *pairs = malloc(ne * 2 * sizeof(u64));
    qsort(mine, ne, sizeof(u64), cmp_u64); qsort(br, ne, sizeof(u64), cmp_u64);
    size_t dm = 0, db = 0;
    for (size_t i = 0; i < ne; i++) { if (i == 0 || mine[i] != mine[i - 1]) dm++; if (i == 0 || br[i] != br[i - 1]) db++; }
    free(pairs);
    /* classes by canonical augmentation */
    u64 *nr = malloc(cap * sizeof(u64)); size_t nn = 0; u64 tmp[1 << (MAXN - 1)];
    for (size_t r = 0; r < nreps; r++) { int m = extend_rep(reps[r], n, tmp); for (int i = 0; i < m; i++) nr[nn++] = tmp[i]; }
    qsort(nr, nn, sizeof(u64), cmp_u64);
    size_t dup = 0; for (size_t i = 1; i < nn; i++) if (nr[i] == nr[i - 1]) dup++;
    printf("brute n=%d: extensions %zu, distinct canonical codes (mine) %zu, distinct brute codes %zu, augmentation classes %zu (duplicates %zu), |Aut| mismatches %ld\n",
           n, ne, dm, db, nn, dup, badaut);
    if (dm != db || dm != nn || dup || badaut) ok = 0;
    free(mine); free(br); free(reps); reps = nr; nreps = nn;
  }
  printf("brute self-test %s\n", ok ? "PASSED" : "FAILED");
  free(reps);
  return ok ? 0 : 1;
}

/* ------------------------------------------------------------------ main */
static unsigned char inS[MAXN + 1][HCAP];   /* inS[m][h] = h attained at order m */
static int leastord[HCAP];

int main(int argc, char **argv) {
  if (argc >= 2 && !strcmp(argv[1], "brute")) return brute_selftest(argc >= 3 ? atoi(argv[2]) : 7);
  if (argc < 3) { fprintf(stderr, "usage: aud_enum NMAX OUTDIR | aud_enum brute 7\n"); return 1; }
  int NMAX = atoi(argv[1]); const char *dir = argv[2];
  if (NMAX > MAXN) NMAX = MAXN;
  int nthr = 1;
#ifdef _OPENMP
  nthr = omp_get_max_threads();
#endif
  dpbuf **DP = malloc(nthr * sizeof(dpbuf *));
  for (int i = 0; i < nthr; i++) DP[i] = malloc(sizeof(dpbuf));
  u64 *reps = malloc(sizeof(u64)); size_t nreps = 1; reps[0] = 0;   /* order 1 */
  for (int n = 1; n <= NMAX; n++) {
    time_t t0 = time(NULL);
    if (n == 10) break;   /* order 10 handled below without storing classes */
    if (n > 1) {
      /* generate order n */
      int *cnt = malloc(nreps * sizeof(int));
      u64 **out = malloc(nreps * sizeof(u64 *));
      #pragma omp parallel for schedule(dynamic, 16)
      for (long r = 0; r < (long)nreps; r++) {
        u64 tmp[1 << (MAXN - 1)];
        int m = extend_rep(reps[r], n, tmp);
        out[r] = malloc((m ? m : 1) * sizeof(u64)); memcpy(out[r], tmp, m * sizeof(u64)); cnt[r] = m;
      }
      size_t tot = 0; for (size_t r = 0; r < nreps; r++) tot += cnt[r];
      u64 *nr = malloc(tot * sizeof(u64)); size_t k = 0;
      for (size_t r = 0; r < nreps; r++) { memcpy(nr + k, out[r], cnt[r] * sizeof(u64)); k += cnt[r]; free(out[r]); }
      free(out); free(cnt); free(reps); reps = nr; nreps = tot;
    }
    /* per-class records */
    char fn[1024]; snprintf(fn, sizeof fn, "%s/aud_cls_%d.txt", dir, n);
    FILE *fo = fopen(fn, "wb");
    if (!fo) { fprintf(stderr, "cannot open %s\n", fn); return 1; }
    char **lines = malloc(nreps * sizeof(char *));
    u64 orbsum = 0; long nsc = 0;
    #pragma omp parallel for schedule(dynamic, 64) reduction(+:orbsum, nsc)
    for (long r = 0; r < (long)nreps; r++) {
      int tid = 0;
#ifdef _OPENMP
      tid = omp_get_thread_num();
#endif
      m16 t[MAXN]; decode(reps[r], n, t);
      u64 dh[MAXN]; u64 h = hpaths(DP[tid], t, n, dh);
      u64 hc = hcycles(DP[tid], t, n);
      cctx C; canon(&C, t, n);
      long aut = C.nbest;
      m16 op[MAXN]; converse(t, n, op); cctx C2; canon(&C2, op, n);
      int sc = (C2.best == C.best);
      if (C.best != reps[r]) { fprintf(stderr, "non-idempotent canonical form!\n"); exit(2); }
      orbsum += factorial(n) / (u64)aut; nsc += sc;
      char buf[512]; int p = snprintf(buf, sizeof buf, "%llu %llu %d %d %d %ld %llu |", reps[r], h, c3count(t, n), is_strong(t, n), sc, aut, hc);
      for (int x = 0; x < n; x++) p += snprintf(buf + p, sizeof buf - p, " %llu", n > 1 ? dh[x] : 0ULL);
      buf[p++] = '\n'; buf[p] = 0;
      lines[r] = strdup(buf);
    }
    for (size_t r = 0; r < nreps; r++) { fputs(lines[r], fo); free(lines[r]); }
    free(lines); fclose(fo);
    for (size_t r = 0; r < nreps; r++) { m16 t[MAXN]; decode(reps[r], n, t); u64 h = hpaths(DP[0], t, n, NULL); if (h >= HCAP) { fprintf(stderr, "H cap\n"); return 1; } inS[n][h] = 1; }
    int nb = n * (n - 1) / 2;
    printf("order %d: classes %zu, SC classes %ld, sum n!/|Aut| = %llu (2^%d = %llu) %s  [%lds]\n", n, nreps, nsc, orbsum, nb, 1ULL << nb,
           orbsum == (1ULL << nb) ? "OK" : "MISMATCH", (long)(time(NULL) - t0));
    fflush(stdout);
  }
  for (int h = 0; h < HCAP; h++) { leastord[h] = 0; for (int m = 1; m <= 9; m++) if (inS[m][h]) { leastord[h] = m; break; } }
  if (NMAX < 10) return 0;

  /* ------------------------------------------------ order 10: stream, aggregate, do not store */
  {
    const int n = 10; time_t t0 = time(NULL);
    typedef struct {
      long all[HCAP], allsc[HCAP], fib[HCAP], fibsc[HCAP], full[HCAP], nonstrong[HCAP], kpos[HCAP];
      int kmin[HCAP], minc3[HCAP], minhc[HCAP];
      long cntc3[HCAP], cnthc[HCAP];
      long kdist[MAXN + 1], hypo, ncls, nfc, scmod4[4], hypo_list_n;
      u64 orbsum;
    } stats;
    stats **ST = malloc(nthr * sizeof(stats *));
    for (int i = 0; i < nthr; i++) {
      ST[i] = calloc(1, sizeof(stats));
      for (int h = 0; h < HCAP; h++) { ST[i]->kmin[h] = 99; ST[i]->minc3[h] = 1 << 30; ST[i]->minhc[h] = 1 << 30; }
    }
    long done = 0;
    #pragma omp parallel for schedule(dynamic, 32)
    for (long r = 0; r < (long)nreps; r++) {
      int tid = 0;
#ifdef _OPENMP
      tid = omp_get_thread_num();
#endif
      stats *S = ST[tid];
      u64 tmp[1 << (MAXN - 1)];
      int m = extend_rep(reps[r], n, tmp);
      for (int i = 0; i < m; i++) {
        m16 t[MAXN]; decode(tmp[i], n, t);
        u64 dh[MAXN]; u64 h = hpaths(DP[tid], t, n, dh);
        if (h >= HCAP) { fprintf(stderr, "H cap\n"); exit(1); }
        int hc = (int)hcycles(DP[tid], t, n);
        cctx C; canon(&C, t, n);
        m16 op[MAXN]; converse(t, n, op); cctx C2; canon(&C2, op, n);
        int sc = (C2.best == C.best);
        int c3 = c3count(t, n);
        int k = 0; for (int x = 0; x < n; x++) if (leastord[dh[x]] == 9) k++;
        S->ncls++; S->orbsum += factorial(n) / (u64)C.nbest;
        S->all[h]++; if (sc) { S->allsc[h]++; S->scmod4[h & 3]++; }
        if (c3 < S->minc3[h]) { S->minc3[h] = c3; S->cntc3[h] = 0; }
        if (c3 == S->minc3[h]) S->cntc3[h]++;
        if (hc < S->minhc[h]) { S->minhc[h] = hc; S->cnthc[h] = 0; }
        if (hc == S->minhc[h]) S->cnthc[h]++;
        if (leastord[h] == 0) {               /* H not attained at any order <= 9: first creator at order 10 */
          S->nfc++; S->fib[h]++; S->kdist[k]++;
          if (sc) S->fibsc[h]++;
          if (k == n) S->full[h]++;
          if (k > 0) S->kpos[h]++;
          if (k < S->kmin[h]) S->kmin[h] = k;
          if (!is_strong(t, n)) S->nonstrong[h]++;
        } else if (k == n) S->hypo++;
      }
      #pragma omp atomic
      done++;
      if ((done & 32767) == 0) { fprintf(stderr, "order 10: %ld / %zu reps [%lds]\n", done, nreps, (long)(time(NULL) - t0)); }
    }
    /* merge */
    stats *A = ST[0];
    for (int i = 1; i < nthr; i++) {
      stats *B = ST[i];
      for (int h = 0; h < HCAP; h++) {
        A->all[h] += B->all[h]; A->allsc[h] += B->allsc[h]; A->fib[h] += B->fib[h]; A->fibsc[h] += B->fibsc[h];
        A->full[h] += B->full[h]; A->nonstrong[h] += B->nonstrong[h]; A->kpos[h] += B->kpos[h];
        if (B->kmin[h] < A->kmin[h]) A->kmin[h] = B->kmin[h];
        if (B->minc3[h] < A->minc3[h]) { A->minc3[h] = B->minc3[h]; A->cntc3[h] = B->cntc3[h]; }
        else if (B->minc3[h] == A->minc3[h]) A->cntc3[h] += B->cntc3[h];
        if (B->minhc[h] < A->minhc[h]) { A->minhc[h] = B->minhc[h]; A->cnthc[h] = B->cnthc[h]; }
        else if (B->minhc[h] == A->minhc[h]) A->cnthc[h] += B->cnthc[h];
      }
      for (int k = 0; k <= MAXN; k++) A->kdist[k] += B->kdist[k];
      A->hypo += B->hypo; A->ncls += B->ncls; A->nfc += B->nfc; A->orbsum += B->orbsum;
      for (int j = 0; j < 4; j++) A->scmod4[j] += B->scmod4[j];
    }
    char fn[1024]; snprintf(fn, sizeof fn, "%s/aud_o10_byH.txt", dir);
    FILE *fo = fopen(fn, "wb");
    fprintf(fo, "# h all_classes SC_classes fc_fiber fc_SC fc_pairs fc_fully_critical fc_nonstrong fc_kmin minc3 count_at_minc3 minhc count_at_minhc\n");
    long vals = 0, uniq = 0, sconly = 0, nsonly = 0, mixed = 0, noparent = 0, nonstrong = 0, maxfib = 0, maxfib_h = 0;
    for (int h = 1; h < HCAP; h++) if (A->all[h]) {
      long s = A->fibsc[h], p = (A->fib[h] - s) / 2;
      if (A->fib[h]) {
        vals++;
        if ((s == 1 && p == 0) || (s == 0 && p == 1)) uniq++;
        if (p == 0) sconly++; else if (s == 0) nsonly++; else mixed++;
        if (!A->kpos[h]) noparent++;
        nonstrong += A->nonstrong[h];
        if (A->fib[h] > maxfib) { maxfib = A->fib[h]; maxfib_h = h; }
        if ((A->fib[h] - s) % 2) printf("ODD NON-SC COUNT at h=%d\n", h);
      }
      fprintf(fo, "%d %ld %ld %ld %ld %ld %ld %ld %d %d %ld %d %ld\n", h, A->all[h], A->allsc[h], A->fib[h], s, A->fib[h] ? p : 0,
              A->full[h], A->nonstrong[h], A->fib[h] ? A->kmin[h] : -1, A->minc3[h], A->cntc3[h], A->minhc[h], A->cnthc[h]);
    }
    fclose(fo);
    printf("order 10: classes %ld, sum 10!/|Aut| = %llu (2^45 = %llu) %s, first creators %ld, hypo analogues %ld  [%lds]\n",
           A->ncls, A->orbsum, 1ULL << 45, A->orbsum == (1ULL << 45) ? "OK" : "MISMATCH", A->nfc, A->hypo, (long)(time(NULL) - t0));
    printf("order 10 k distribution (first creators):"); for (int k = 0; k <= n; k++) printf(" %d:%ld", k, A->kdist[k]); printf("\n");
    printf("order 10 SC classes by H mod 4: 1:%ld 3:%ld (total %ld)\n", A->scmod4[1], A->scmod4[3], A->scmod4[1] + A->scmod4[3]);
    printf("order 10 new values %ld; dichotomy %ld; SC-only %ld NS-only %ld mixed %ld; values without lineage %ld; non-strong first creators %ld; max fiber %ld at h=%ld\n",
           vals, uniq, sconly, nsonly, mixed, noparent, nonstrong, maxfib, maxfib_h);
  }
  return 0;
}
