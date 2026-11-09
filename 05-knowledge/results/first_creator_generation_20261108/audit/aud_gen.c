/* aud_gen.c - independent auditor's census for the packet first_creator_generation_20261108 (PREREGISTRATION.md,
 * PREREGISTRATION-2.md, PROOF.md sections 3.2-3.3).  Auditor: independent audit of kind-pasteur-2026-11-08-S17.
 *
 * Class enumeration: canonical augmentation with the canonical labelling of the PREVIOUS independent audit
 * (05-knowledge/results/first_creators_genealogy_20261107/audit/aud_enum.c; the block between the two "copied"
 * markers below is copied verbatim from it: cctx, refine, search, canon, decode, converse).  Everything else
 * (Hamiltonian-path DP, arc-strong connectivity, module test, categories, attachment form, statistics) is written
 * from scratch for this audit and shares no code with the packet's ext.c / gen.c / lam.c / dump10.c / merge.py.
 *
 * Definitions used (taken from the packet's text, not its code):
 *   S(m)      = set of H over all m-tournaments (computed here from the classes of order m);
 *   FC        = H(T) not in S(n-1);
 *   EXT       = not FC, and some T - v has H(T - v) not in S(n-2) (a first creator of order n-1);
 *   SUM       = otherwise, T not strong;
 *   SUB       = otherwise, T has a module M with 2 <= |M| <= n-1 (every outside vertex beats all of M or none of M);
 *   REST      = otherwise;
 *   pendant   = some vertex of out-degree 1 or in-degree 1;
 *   lambda(T) = min over nonempty proper X of #arcs from V\X into X  (brute force over all 2^n - 2 sets X);
 *   f_k(n)    = min H over n-tournaments with lambda >= k;  mindeg = min over v of min(outdeg, indeg).
 * Passes: A = orders 1..9 (classes stored), C = order 10 (streamed, OpenMP), B = orders 3..9 extended by every
 * pattern through the attachment form (fertility, value genealogy, the "reached only from first creators" share).
 * usage: aud_gen OUTDIR [nthreads]
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
typedef unsigned long long u64;
typedef unsigned short m16;

/* ======================== copied verbatim from aud_enum.c (previous independent audit) ======================== */
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
/* ======================== end of copied block; everything below is the auditor's own code ===================== */

#define HCAP 100000

typedef struct { u64 f[1 << MAXN][MAXN]; } dpw;

/* H(T) and dh[x] = H(T - x): subset DP over paths counted by their last vertex */
static u64 ham_all(dpw *D, const m16 *out, int n, u64 *dh) {
  int F = 1 << n;
  memset(D->f, 0, sizeof(u64) * MAXN * (size_t)F);
  for (int v = 0; v < n; v++) D->f[1 << v][v] = 1;
  for (int S = 1; S < F; S++) {
    for (int v = 0; v < n; v++) {
      u64 x = D->f[S][v];
      if (!x) continue;
      m16 nx = (m16)(out[v] & ~S);
      while (nx) { int u = __builtin_ctz(nx); nx &= nx - 1; D->f[S | (1 << u)][u] += x; }
    }
  }
  u64 h = 0;
  for (int v = 0; v < n; v++) h += D->f[F - 1][v];
  if (dh) for (int x = 0; x < n; x++) {
    if (n == 1) { dh[x] = 1; continue; }
    int S = (F - 1) & ~(1 << x); u64 s = 0;
    for (int v = 0; v < n; v++) s += D->f[S][v];
    dh[x] = s;
  }
  return h;
}

static void in_masks(const m16 *out, int n, m16 *in) {
  for (int v = 0; v < n; v++) in[v] = 0;
  for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) if ((out[u] >> v) & 1) in[v] |= (m16)(1u << u);
}

/* arc-strong connectivity by brute force over all nonempty proper X: arcs from V\X into X = sum_{x in X} |in(x) \ X| */
static int lambda_arc(const m16 *out, int n) {
  if (n <= 1) return 0;
  m16 in[MAXN]; in_masks(out, n, in);
  int full = (1 << n) - 1, best = 1 << 30;
  for (int X = 1; X < full; X++) {
    int b = 0; m16 nX = (m16)(full & ~X);
    for (int x = 0; x < n && b < best; x++) if ((X >> x) & 1) b += __builtin_popcount(in[x] & nX);
    if (b < best) best = b;
  }
  return best;
}

/* strong connectivity of the subtournament induced by mask (forward and backward closure from one vertex) */
static int strong_mask(const m16 *out, int n, int mask) {
  int c = __builtin_popcount(mask);
  if (c <= 1) return 1;
  m16 in[MAXN]; in_masks(out, n, in);
  int s = __builtin_ctz(mask);
  int r = 1 << s, prev = 0;
  while (r != prev) { prev = r; for (int v = 0; v < n; v++) if ((r >> v) & 1) r |= out[v] & mask; }
  if (r != mask) return 0;
  r = 1 << s; prev = 0;
  while (r != prev) { prev = r; for (int v = 0; v < n; v++) if ((r >> v) & 1) r |= in[v] & mask; }
  return r == mask;
}

/* a nontrivial module (2 <= |M| <= n-1), brute force over all subsets; returns its mask or 0 */
static int find_module(const m16 *out, int n) {
  int full = (1 << n) - 1;
  for (int M = 1; M < full; M++) {
    if (__builtin_popcount(M) < 2) continue;
    int ok = 1;
    for (int z = 0; z < n && ok; z++) if (!((M >> z) & 1)) { int o = out[z] & M; if (o != 0 && o != M) ok = 0; }
    if (ok) return M;
  }
  return 0;
}

static u64 factorial(int n) { u64 f = 1; for (int i = 2; i <= n; i++) f *= i; return f; }
static int cmp_pair(const void *a, const void *b) {
  u64 x = ((const u64 *)a)[0], y = ((const u64 *)b)[0]; return x < y ? -1 : x > y;
}

/* canonical augmentation (same acceptance rule as aud_enum.c's extend_rep), also returning |Aut| of each child */
static int extend_rep2(u64 rc, int n, u64 *code, long *aut) {
  m16 base[MAXN]; decode(rc, n - 1, base);
  cctx C; int cnt = 0;
  u64 pairs[2 * (1 << (MAXN - 1))];
  for (int P = 0; P < (1 << (n - 1)); P++) {
    m16 t[MAXN];
    for (int v = 0; v < n - 1; v++) t[v] = base[v];
    t[n - 1] = 0;
    for (int v = 0; v < n - 1; v++) { if ((P >> v) & 1) t[n - 1] |= (m16)(1u << v); else t[v] |= (m16)(1u << (n - 1)); }
    canon(&C, t, n);
    if (uf_find(C.uf, n - 1) != uf_find(C.uf, C.lab0[n - 1])) continue;
    pairs[2 * cnt] = C.best; pairs[2 * cnt + 1] = (u64)C.nbest; cnt++;
  }
  qsort(pairs, cnt, 2 * sizeof(u64), cmp_pair);
  int m = 0;
  for (int i = 0; i < cnt; i++) if (i == 0 || pairs[2 * i] != pairs[2 * (i - 1)]) { code[m] = pairs[2 * i]; aut[m] = (long)pairs[2 * i + 1]; m++; }
  return m;
}

/* ------------------------------------------------------------------ statistics */
typedef struct {
  long ncls, nfc, next, nsum, nsub, nrest, nstrong;
  long rest_pend, rest_graft, rest_neither, rest_nsdel, rest_no_nsdel, rest_nonpend_nsdel;
  long rest_lam[MAXN + 1], rest_no_nsdel_lam[MAXN + 1];
  long lam_cnt[MAXN + 1], lam_fc[MAXN + 1], lam_nfc[MAXN + 1];
  u64 lam_minH[MAXN + 1], lam_maxH[MAXN + 1], lam_nfc_minH[MAXN + 1], lam_nfc_maxH[MAXN + 1];
  long lam_minH_cnt[MAXN + 1]; u64 lam_minH_code[MAXN + 1];
  long deg_cnt[MAXN + 1], deg_nfc[MAXN + 1], deg_fc[MAXN + 1];
  long lam_deg[MAXN + 1][MAXN + 1], lam_deg_nfc[MAXN + 1][MAXN + 1];
  long kdist[MAXN + 2], hypo;
  u64 orbsum, maxH;
  long sum_nonstrong_fc;
} stats;

static void stats_init(stats *S) {
  memset(S, 0, sizeof(stats));
  for (int l = 0; l <= MAXN; l++) { S->lam_minH[l] = S->lam_nfc_minH[l] = ~0ULL; }
}

static void stats_merge(stats *A, const stats *B) {
  long *a = (long *)A; const long *b = (const long *)B;
  /* add every long counter field by field (explicitly, to keep min/max fields separate) */
  (void)a; (void)b;
  A->ncls += B->ncls; A->nfc += B->nfc; A->next += B->next; A->nsum += B->nsum; A->nsub += B->nsub; A->nrest += B->nrest;
  A->nstrong += B->nstrong;
  A->rest_pend += B->rest_pend; A->rest_graft += B->rest_graft; A->rest_neither += B->rest_neither;
  A->rest_nsdel += B->rest_nsdel; A->rest_no_nsdel += B->rest_no_nsdel; A->rest_nonpend_nsdel += B->rest_nonpend_nsdel;
  for (int l = 0; l <= MAXN; l++) {
    A->rest_lam[l] += B->rest_lam[l]; A->rest_no_nsdel_lam[l] += B->rest_no_nsdel_lam[l];
    A->lam_cnt[l] += B->lam_cnt[l]; A->lam_fc[l] += B->lam_fc[l]; A->lam_nfc[l] += B->lam_nfc[l];
    if (B->lam_minH[l] < A->lam_minH[l]) { A->lam_minH[l] = B->lam_minH[l]; A->lam_minH_cnt[l] = B->lam_minH_cnt[l]; A->lam_minH_code[l] = B->lam_minH_code[l]; }
    else if (B->lam_minH[l] == A->lam_minH[l] && B->lam_minH[l] != ~0ULL) A->lam_minH_cnt[l] += B->lam_minH_cnt[l];
    if (B->lam_maxH[l] > A->lam_maxH[l]) A->lam_maxH[l] = B->lam_maxH[l];
    if (B->lam_nfc_minH[l] < A->lam_nfc_minH[l]) A->lam_nfc_minH[l] = B->lam_nfc_minH[l];
    if (B->lam_nfc_maxH[l] > A->lam_nfc_maxH[l]) A->lam_nfc_maxH[l] = B->lam_nfc_maxH[l];
    A->deg_cnt[l] += B->deg_cnt[l]; A->deg_nfc[l] += B->deg_nfc[l]; A->deg_fc[l] += B->deg_fc[l];
    for (int d = 0; d <= MAXN; d++) { A->lam_deg[l][d] += B->lam_deg[l][d]; A->lam_deg_nfc[l][d] += B->lam_deg_nfc[l][d]; }
  }
  for (int k = 0; k <= MAXN + 1; k++) A->kdist[k] += B->kdist[k];
  A->hypo += B->hypo; A->orbsum += B->orbsum; if (B->maxH > A->maxH) A->maxH = B->maxH;
  A->sum_nonstrong_fc += B->sum_nonstrong_fc;
}

static unsigned char inS[MAXN + 1][HCAP];      /* inS[m][h]: h attained at order m (m <= 10) */
static unsigned char vNewSeen[MAXN + 1][HCAP];  /* some first creator of order m has H = h */
static unsigned char vReached[MAXN + 1][HCAP];  /* some first creator of order m with H = h has a first-creator deletion */
static unsigned char vNotFull[MAXN + 1][HCAP];  /* some first creator of order m with H = h has a non-first-creator deletion */
static FILE *frest;                              /* REST non-pendant classes */

typedef struct { unsigned char reached[HCAP], notfull[HCAP], newseen[HCAP], inS[HCAP]; } vbuf;

static void analyze(dpw *D, const m16 *t, int n, u64 code, long aut, stats *S, vbuf *V) {
  u64 dh[MAXN];
  u64 h = ham_all(D, t, n, dh);
  if (h >= HCAP) { fprintf(stderr, "H cap\n"); exit(3); }
  S->ncls++; S->orbsum += factorial(n) / (u64)aut;
  if (h > S->maxH) S->maxH = h;
  V->inS[h] = 1;
  int fc = !inS[n - 1][h];
  int k = 0;
  if (n >= 3) for (int x = 0; x < n; x++) if (!inS[n - 2][dh[x]]) k++;
  int lam = lambda_arc(t, n);
  int mindeg = n;
  for (int v = 0; v < n; v++) { int od = __builtin_popcount(t[v]), id = n - 1 - od; if (od < mindeg) mindeg = od; if (id < mindeg) mindeg = id; }
  int strong = lam >= 1;
  if (strong) S->nstrong++;
  S->lam_cnt[lam]++;
  if (h < S->lam_minH[lam]) { S->lam_minH[lam] = h; S->lam_minH_cnt[lam] = 1; S->lam_minH_code[lam] = code; }
  else if (h == S->lam_minH[lam]) S->lam_minH_cnt[lam]++;
  if (h > S->lam_maxH[lam]) S->lam_maxH[lam] = h;
  S->deg_cnt[mindeg]++;
  S->lam_deg[lam][mindeg]++;
  if (fc) {
    S->nfc++; S->lam_fc[lam]++; S->deg_fc[mindeg]++; S->kdist[k]++;
    if (!strong) S->sum_nonstrong_fc++;
    V->newseen[h] = 1;
    if (k >= 1) V->reached[h] = 1;
    if (k < n) V->notfull[h] = 1;
    return;
  }
  S->lam_nfc[lam]++; S->deg_nfc[mindeg]++; S->lam_deg_nfc[lam][mindeg]++;
  if (h < S->lam_nfc_minH[lam]) S->lam_nfc_minH[lam] = h;
  if (h > S->lam_nfc_maxH[lam]) S->lam_nfc_maxH[lam] = h;
  if (n >= 3 && k == n) S->hypo++;
  if (k >= 1) { S->next++; return; }
  if (!strong) { S->nsum++; return; }
  if (find_module(t, n)) { S->nsub++; return; }
  S->nrest++;
  S->rest_lam[lam]++;
  int pend = 0;
  for (int v = 0; v < n; v++) { int od = __builtin_popcount(t[v]); if (od == 1 || n - 1 - od == 1) pend = 1; }
  int nsdel = 0;
  for (int x = 0; x < n && !nsdel; x++) if (!strong_mask(t, n, ((1 << n) - 1) & ~(1 << x))) nsdel = 1;
  if (nsdel) S->rest_nsdel++; else { S->rest_no_nsdel++; S->rest_no_nsdel_lam[lam]++; }
  if (pend) S->rest_pend++;
  else {
    if (lam == 1) S->rest_graft++; else S->rest_neither++;
    if (nsdel) S->rest_nonpend_nsdel++;
    /* record: order, code, H, lambda, mindeg, nonstrong-deletion flag, SC, |Aut| */
    m16 op[MAXN]; converse(t, n, op);
    cctx C1, C2; canon(&C1, t, n); canon(&C2, op, n);
    #ifdef _OPENMP
    #pragma omp critical
    #endif
    {
      fprintf(frest, "%d %llu %llu %d %d %d %d %ld\n", n, code, h, lam, mindeg, nsdel, C1.best == C2.best, C1.nbest);
    }
  }
}

static void print_stats(int n, const stats *S, FILE *fo) {
  fprintf(fo, "ORDER %d: classes %ld (strong %ld), sum n!/|Aut| = %llu vs 2^%d = %llu %s; max H %llu\n", n, S->ncls, S->nstrong,
          S->orbsum, n * (n - 1) / 2, 1ULL << (n * (n - 1) / 2), S->orbsum == (1ULL << (n * (n - 1) / 2)) ? "OK" : "MISMATCH", S->maxH);
  long nonfc = S->ncls - S->nfc;
  fprintf(fo, "  FC %ld; non-FC %ld = EXT %ld + SUM %ld + SUB %ld + REST %ld (sum %ld)\n", S->nfc, nonfc, S->next, S->nsum, S->nsub,
          S->nrest, S->next + S->nsum + S->nsub + S->nrest);
  if (nonfc) fprintf(fo, "  shares of non-FC: EXT %.4f%%, REST %.4f%%\n", 100.0 * S->next / nonfc, 100.0 * S->nrest / nonfc);
  fprintf(fo, "  REST: pendant %ld, non-pendant lambda=1 %ld, non-pendant lambda>=2 %ld; with a non-strong deletion %ld, without %ld;"
          " non-pendant with a non-strong deletion %ld\n", S->rest_pend, S->rest_graft, S->rest_neither, S->rest_nsdel, S->rest_no_nsdel,
          S->rest_nonpend_nsdel);
  if (S->nrest) fprintf(fo, "  REST pendant share %.4f%%\n", 100.0 * S->rest_pend / S->nrest);
  fprintf(fo, "  REST lambda histogram:"); for (int l = 0; l <= MAXN; l++) if (S->rest_lam[l]) fprintf(fo, " %d:%ld", l, S->rest_lam[l]); fprintf(fo, "\n");
  fprintf(fo, "  REST without non-strong deletion, lambda histogram:"); for (int l = 0; l <= MAXN; l++) if (S->rest_no_nsdel_lam[l]) fprintf(fo, " %d:%ld", l, S->rest_no_nsdel_lam[l]); fprintf(fo, "\n");
  for (int l = 0; l <= MAXN; l++) if (S->lam_cnt[l])
    fprintf(fo, "  lambda=%d: classes %ld, FC %ld, non-FC %ld; H in [%llu, %llu] (min attained by %ld classes, e.g. code %llu); non-FC H in [%lld, %lld]\n",
            l, S->lam_cnt[l], S->lam_fc[l], S->lam_nfc[l], S->lam_minH[l], S->lam_maxH[l], S->lam_minH_cnt[l], S->lam_minH_code[l],
            S->lam_nfc[l] ? (long long)S->lam_nfc_minH[l] : -1LL, S->lam_nfc[l] ? (long long)S->lam_nfc_maxH[l] : -1LL);
  fprintf(fo, "  f_k(%d) = min H over lambda >= k:", n);
  for (int k = 1; k <= MAXN; k++) { u64 m = ~0ULL; for (int l = k; l <= MAXN; l++) if (S->lam_cnt[l] && S->lam_minH[l] < m) m = S->lam_minH[l]; if (m != ~0ULL) fprintf(fo, " f_%d=%llu", k, m); }
  fprintf(fo, "\n");
  for (int d = 0; d <= MAXN; d++) if (S->deg_cnt[d]) fprintf(fo, "  mindeg=%d: classes %ld, FC %ld, non-FC %ld\n", d, S->deg_cnt[d], S->deg_fc[d], S->deg_nfc[d]);
  long mism = 0; for (int l = 0; l <= MAXN; l++) for (int d = 0; d <= MAXN; d++) if ((l >= 3) != (d >= 3)) mism += S->lam_deg[l][d];
  fprintf(fo, "  (lambda>=3) != (mindeg>=3): %ld classes\n", mism);
  fprintf(fo, "  lambda x mindeg (classes / non-FC):");
  for (int l = 0; l <= MAXN; l++) for (int d = 0; d <= MAXN; d++) if (S->lam_deg[l][d]) fprintf(fo, " (%d,%d):%ld/%ld", l, d, S->lam_deg[l][d], S->lam_deg_nfc[l][d]);
  fprintf(fo, "\n");
  fprintf(fo, "  FC by number k of first-creator deletions:"); for (int k = 0; k <= n; k++) if (S->kdist[k]) fprintf(fo, " %d:%ld", k, S->kdist[k]); fprintf(fo, "\n");
  fprintf(fo, "  non-FC with every deletion a first creator (hypo analogue): %ld; non-strong FC: %ld\n", S->hypo, S->sum_nonstrong_fc);
  fflush(fo);
}

static void merge_values(int n, vbuf *V) {
  for (int h = 0; h < HCAP; h++) {
    if (V->inS[h]) inS[n][h] = 1;
    if (V->newseen[h]) vNewSeen[n][h] = 1;
    if (V->reached[h]) vReached[n][h] = 1;
    if (V->notfull[h]) vNotFull[n][h] = 1;
  }
}

static void value_report(int n, FILE *fo) {
  long nS = 0, nnew = 0, unreached = 0, onlyfc = 0; long mx = 0;
  fprintf(fo, "  VALUES order %d:", n);
  for (int h = 0; h < HCAP; h++) if (inS[n][h]) { nS++; mx = h; }
  for (int h = 0; h < HCAP; h++) if (vNewSeen[n][h]) { nnew++; if (!vReached[n][h]) unreached++; if (!vNotFull[n][h]) onlyfc++; }
  fprintf(fo, " |S(%d)| = %ld (max %ld), new values %ld, new values with no first creator having a first-creator deletion %ld",
          n, nS, mx, nnew, unreached);
  if (unreached && unreached < 20) { fprintf(fo, " ["); for (int h = 0; h < HCAP; h++) if (vNewSeen[n][h] && !vReached[n][h]) fprintf(fo, " %d", h); fprintf(fo, " ]"); }
  fprintf(fo, ", new values all of whose first creators have ONLY first-creator deletions %ld (%.4f)\n", onlyfc, nnew ? (double)onlyfc / nnew : 0.0);
  fflush(fo);
}

/* ------------------------------------------------------------------ attachment form (pass B) */
typedef struct { long long fe[1 << (MAXN - 1)][MAXN], gs[1 << (MAXN - 1)][MAXN]; } fgbuf;

/* for T of order n <= 9: H(T + v_x) for all x via H + sum x_a L_a - sum_{a<b} x_a x_b M_ab */
static void ext_values(fgbuf *G, const m16 *t, int n, long long *hx) {
  int F = 1 << n;
  for (int S = 0; S < F; S++) for (int a = 0; a < n; a++) { G->fe[S][a] = 0; G->gs[S][a] = 0; }
  for (int a = 0; a < n; a++) { G->fe[1 << a][a] = 1; G->gs[1 << a][a] = 1; }
  m16 in[MAXN]; in_masks(t, n, in);
  for (int S = 1; S < F; S++) for (int a = 0; a < n; a++) if ((S >> a) & 1) {
    long long x = G->fe[S][a];
    if (x) { m16 o = (m16)(t[a] & ~S); while (o) { int b = __builtin_ctz(o); o &= o - 1; G->fe[S | (1 << b)][b] += x; } }
    long long y = G->gs[S][a];
    if (y) { m16 o = (m16)(in[a] & ~S); while (o) { int b = __builtin_ctz(o); o &= o - 1; G->gs[S | (1 << b)][b] += y; } }
  }
  long long H = 0; for (int a = 0; a < n; a++) H += G->fe[F - 1][a];
  long long W[MAXN][MAXN]; memset(W, 0, sizeof W);
  for (int S = 1; S < F - 1; S++) { int C = (F - 1) & ~S;
    for (int a = 0; a < n; a++) if (((S >> a) & 1) && G->fe[S][a]) for (int b = 0; b < n; b++) if ((C >> b) & 1) W[a][b] += G->fe[S][a] * G->gs[C][b]; }
  long long L[MAXN], M[MAXN][MAXN];
  for (int a = 0; a < n; a++) { L[a] = G->fe[F - 1][a] - G->gs[F - 1][a]; for (int b = 0; b < n; b++) if (b != a) L[a] += W[a][b]; }
  for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) M[a][b] = W[a][b] + W[b][a];
  for (int x = 0; x < F; x++) {
    long long v = H;
    for (int a = 0; a < n; a++) if ((x >> a) & 1) { v += L[a]; for (int b = a + 1; b < n; b++) if ((x >> b) & 1) v -= M[a][b]; }
    hx[x] = v;
  }
}

/* ------------------------------------------------------------------ main */
int main(int argc, char **argv) {
  if (argc < 2) { fprintf(stderr, "usage: aud_gen OUTDIR [nthreads]\n"); return 1; }
  const char *dir = argv[1];
  int nthr = 1;
#ifdef _OPENMP
  if (argc >= 3) omp_set_num_threads(atoi(argv[2]));
  nthr = omp_get_max_threads();
#endif
  char fn[1024];
  snprintf(fn, sizeof fn, "%s/aud_gen.out", dir); FILE *fo = fopen(fn, "w");
  snprintf(fn, sizeof fn, "%s/aud_rest_nonpendant.txt", dir); frest = fopen(fn, "w");
  fprintf(frest, "# order code H lambda mindeg has_nonstrong_deletion self_converse |Aut|  (REST classes with no vertex of in- or out-degree 1)\n");
  if (!fo || !frest) { fprintf(stderr, "cannot open outputs\n"); return 1; }
  inS[0][1] = 1;   /* empty tournament: one (empty) path; only used as S(0) for order-1/2 bookkeeping */
  dpw **D = malloc(sizeof(dpw *) * nthr); for (int i = 0; i < nthr; i++) D[i] = malloc(sizeof(dpw));
  vbuf **VB = malloc(sizeof(vbuf *) * nthr); for (int i = 0; i < nthr; i++) VB[i] = calloc(1, sizeof(vbuf));
  stats **ST = malloc(sizeof(stats *) * nthr); for (int i = 0; i < nthr; i++) ST[i] = malloc(sizeof(stats));
  u64 *reps[MAXN + 1]; size_t nreps[MAXN + 1];
  reps[1] = malloc(sizeof(u64)); reps[1][0] = 0; nreps[1] = 1;
  /* ---------------- pass A: orders 1..9 */
  for (int n = 1; n <= 9; n++) {
    time_t t0 = time(NULL);
    if (n > 1) {
      size_t cap = nreps[n - 1] << (n - 1);
      reps[n] = malloc(sizeof(u64) * cap); size_t m = 0;
      u64 code[1 << (MAXN - 1)]; long aut[1 << (MAXN - 1)];
      for (size_t r = 0; r < nreps[n - 1]; r++) { int c = extend_rep2(reps[n - 1][r], n, code, aut); for (int i = 0; i < c; i++) reps[n][m++] = code[i]; }
      nreps[n] = m;
    }
    stats S; stats_init(&S); vbuf *V = VB[0]; memset(V, 0, sizeof(vbuf));
    for (size_t r = 0; r < nreps[n]; r++) {
      m16 t[MAXN]; decode(reps[n][r], n, t);
      cctx C; canon(&C, t, n);
      if (C.best != reps[n][r]) { fprintf(stderr, "non-idempotent canonical code\n"); return 2; }
      analyze(D[0], t, n, reps[n][r], C.nbest, &S, V);
    }
    merge_values(n, V);
    print_stats(n, &S, fo); value_report(n, fo);
    fprintf(fo, "  [order %d: %lds]\n", n, (long)(time(NULL) - t0));
    snprintf(fn, sizeof fn, "%s/aud_S%d.txt", dir, n); FILE *fs = fopen(fn, "w");
    for (int h = 0; h < HCAP; h++) if (inS[n][h]) fprintf(fs, "%d\n", h);
    fclose(fs);
    fflush(frest);
  }
  /* ---------------- pass C: order 10, streamed */
  {
    const int n = 10; time_t t0 = time(NULL);
    for (int i = 0; i < nthr; i++) { stats_init(ST[i]); memset(VB[i], 0, sizeof(vbuf)); }
    snprintf(fn, sizeof fn, "%s/aud_o10_codes.bin", argc >= 4 ? argv[3] : dir);
    FILE *fb = (argc >= 4) ? fopen(fn, "wb") : NULL;
    long done = 0;
    #pragma omp parallel for schedule(dynamic, 16)
    for (long r = 0; r < (long)nreps[9]; r++) {
      int tid = 0;
#ifdef _OPENMP
      tid = omp_get_thread_num();
#endif
      u64 code[1 << (MAXN - 1)]; long aut[1 << (MAXN - 1)];
      int c = extend_rep2(reps[9][r], n, code, aut);
      for (int i = 0; i < c; i++) { m16 t[MAXN]; decode(code[i], n, t); analyze(D[tid], t, n, code[i], aut[i], ST[tid], VB[tid]); }
      if (fb) {
        #pragma omp critical
        { fwrite(code, sizeof(u64), c, fb); }
      }
      #pragma omp atomic
      done++;
      if ((done & 16383) == 0) { fprintf(stderr, "order 10: %ld / %zu [%lds]\n", done, nreps[9], (long)(time(NULL) - t0)); }
    }
    if (fb) fclose(fb);
    for (int i = 1; i < nthr; i++) stats_merge(ST[0], ST[i]);
    for (int i = 0; i < nthr; i++) merge_values(n, VB[i]);
    print_stats(n, ST[0], fo); value_report(n, fo);
    fprintf(fo, "  [order 10: %lds, %d threads]\n", (long)(time(NULL) - t0), nthr);
    snprintf(fn, sizeof fn, "%s/aud_S10.txt", dir); FILE *fs = fopen(fn, "w");
    for (int h = 0; h < HCAP; h++) if (inS[n][h]) fprintf(fs, "%d\n", h);
    fclose(fs);
  }
  /* ---------------- pass B: orders 3..9 extended by every pattern (attachment form) */
  {
    fgbuf *G = malloc(sizeof(fgbuf));
    static long long A[HCAP], Fp[HCAP];
    static long long hx[1 << (MAXN - 1)];
    for (int n = 3; n <= 9; n++) {
      memset(A, 0, sizeof A); memset(Fp, 0, sizeof Fp);
      long nfc = 0, sterile = 0, below_half = 0; long minfert = 1 << 30; u64 minfert_code = 0; long long minfert_h = 0;
      fprintf(fo, "EXTENSION order %d -> %d:\n", n, n + 1);
      for (size_t r = 0; r < nreps[n]; r++) {
        m16 t[MAXN]; decode(reps[n][r], n, t);
        u64 h = ham_all(D[0], t, n, NULL);
        int fc = !inS[n - 1][h];
        ext_values(G, t, n, hx);
        long fert = 0;
        for (int x = 0; x < (1 << n); x++) {
          long long v = hx[x];
          if (v < 0 || v >= HCAP) { fprintf(stderr, "ext value out of range\n"); return 4; }
          A[v]++;
          if (fc) { Fp[v]++; if (!inS[n][v]) fert++; }
        }
        if (fc) {
          nfc++;
          if (fert == 0) sterile++;
          if (2 * fert < (1L << n)) { below_half++; fprintf(fo, "  FC below half: code %llu H %llu fertile %ld of %d\n", reps[n][r], h, fert, 1 << n); }
          if (fert < minfert) { minfert = fert; minfert_code = reps[n][r]; minfert_h = (long long)h; }
        }
      }
      long got = 0, extra = 0, missing = 0, newv = 0, unreached = 0, onlyfc = 0;
      for (int v = 0; v < HCAP; v++) {
        if (A[v]) { got++; if (!inS[n + 1][v]) extra++; }
        if (inS[n + 1][v] && !A[v]) missing++;
        if (inS[n + 1][v] && !inS[n][v]) { newv++; if (!Fp[v]) unreached++; if (Fp[v] == A[v]) onlyfc++; }
      }
      fprintf(fo, "  values of all extensions %ld (missing from S(%d): %ld, extra: %ld); new values at %d: %ld; not reached from an FC: %ld",
              got, n + 1, missing, extra, n + 1, newv, unreached);
      if (unreached && unreached < 20) { fprintf(fo, " ["); for (int v = 0; v < HCAP; v++) if (inS[n + 1][v] && !inS[n][v] && !Fp[v]) fprintf(fo, " %d", v); fprintf(fo, " ]"); }
      fprintf(fo, "; reached ONLY from FC parents: %ld (%.4f)\n", onlyfc, newv ? (double)onlyfc / newv : 0.0);
      fprintf(fo, "  FC %ld: sterile %ld, below half %ld, min fertile %ld of %d (code %llu, H %lld)\n", nfc, sterile, below_half, minfert,
              1 << n, minfert_code, minfert_h);
      fflush(fo);
    }
  }
  fclose(fo); fclose(frest);
  return 0;
}
