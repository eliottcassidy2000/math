/* aud_o11_cls.c (derived from the auditor's aud_o11.c; same engine) - class-level genealogy at order 11.
 * An order-11 first creator T' with NO first-creator deletion has T = T' - v a non-first-creator of order 10 for every v,
 * so extending every NON-first-creator class of order 10 by every pattern finds all of them (no order-11 canonical form
 * is needed).  For each such T and x with H(T + v_x) not in S(10), test every deletion T' - u (u in V) by the deletion forms;
 * T' - v = T itself is a non-first-creator.  usage: aud_o11_cls CODES.bin S9.txt S10.txt OUT.txt [nthreads]
 *
 * original header of aud_o11.c follows:
 * aud_o11.c - independent auditor's exhaustive order-10 -> order-11 extension run (phase 2 of the audit of
 * 05-knowledge/results/first_creator_generation_20261108).  Written from scratch; shares no code with ext.c / lam11.c.
 *
 * Input: the auditor's OWN order-10 classes (aud_gen.c: canonical augmentation; codes in aud_enum.c's bit order, the
 * pair (0,1) being the most significant of the 45 bits, bit set = i -> j).
 * For every class T and every pattern x (bit a = [a -> v]):
 *   H(T + v_x) by the attachment form H + sum x_a L_a - sum_{a<b} x_a x_b M_ab (validated by aud_attach.c, and here
 *   against a direct subset DP of the 11-vertex tournament on a sample of classes);
 *   for first creators also every deletion H(T + v_x - u), u in V, by the deletion forms (the f/g tables of T restricted
 *   to V \ u; also validated against the direct DP on the sample).
 * Reports: R1 (every value new at 11 has a first-creator parent), R2 (hypo pairs), R3/R4 (fertility), R5 (new values
 * whose every (class, pattern) realization has a first-creator parent), realization multiplicities, S(11) vs klein;
 * and the exhaustive lambda >= 3 search: parents with min in/out-degree >= 2, patterns giving min degree >= 3.
 * usage: aud_o11 CODES.bin S9.txt S10.txt least_order_11.txt OUTDIR [nthreads] [sample_stride]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef N
#define N 10
#endif
#define HC 100000
typedef unsigned long long u64;
typedef long long i64;
typedef unsigned short m16;

static unsigned char inS9[HC], inS10[HC], inS11k[HC];

static void decode_msb(u64 code, int n, m16 *out) {
  memset(out, 0, sizeof(m16) * 16);
  int b = n * (n - 1) / 2 - 1;
  for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
    if ((code >> b) & 1) out[i] |= (m16)(1u << j); else out[j] |= (m16)(1u << i);
    b--;
  }
}

typedef struct {
  i64 fe[1 << N][N], gs[1 << N][N];
  i64 W[N][N], L[N], M[N][N], H;
  i64 Lu[N][N], Mu[N][N][N], Hu[N];
  i64 f11[1 << (N + 1)][N + 1];
} work;

typedef struct {
  u64 A[HC], F[HC];
  long ncls, nfc, sterile, below512, hypo, fert_hist[1025];
  long minfert; u64 minfert_codes[16]; int nminfert; i64 minfert_H;
  long long lam3_parents, lam3_ext; i64 lam3_min; u64 lam3_min_code; int lam3_min_x;   /* long long: 'long' is 32-bit on Windows */
  long sample_cls, sample_fc, sample_checks, sample_fail, sample_fert_mismatch, sample_hypo_bf;
  long hypo_rec;
} tstats;

static FILE *flam, *fhypo;

static void dp_tables(work *w, const m16 *t, const m16 *in) {
  const int F = 1 << N;
  memset(w->fe, 0, sizeof(w->fe)); memset(w->gs, 0, sizeof(w->gs));
  for (int a = 0; a < N; a++) { w->fe[1 << a][a] = 1; w->gs[1 << a][a] = 1; }
  for (int S = 1; S < F; S++) {
    int s = S;
    while (s) {
      int a = __builtin_ctz(s); s &= s - 1;
      i64 x = w->fe[S][a];
      if (x) { int o = t[a] & ~S; while (o) { int b = __builtin_ctz(o); o &= o - 1; w->fe[S | (1 << b)][b] += x; } }
      i64 y = w->gs[S][a];
      if (y) { int o = in[a] & ~S; while (o) { int b = __builtin_ctz(o); o &= o - 1; w->gs[S | (1 << b)][b] += y; } }
    }
  }
}

/* form over the vertex set U (U = V for T, U = V \ u for a deletion): returns H(T[U]); fills W (|V| x |V| indices) */
static i64 form_on(work *w, int U, i64 Wm[N][N], i64 *Lv) {
  for (int a = 0; a < N; a++) for (int b = 0; b < N; b++) Wm[a][b] = 0;
  for (int S = (U - 1) & U; S > 0; S = (S - 1) & U) {
    int C = U ^ S;
    int s = S;
    while (s) {
      int a = __builtin_ctz(s); s &= s - 1;
      i64 x = w->fe[S][a]; if (!x) continue;
      int c = C;
      while (c) { int b = __builtin_ctz(c); c &= c - 1; Wm[a][b] += x * w->gs[C][b]; }
    }
  }
  i64 H = 0;
  for (int a = 0; a < N; a++) if ((U >> a) & 1) {
    H += w->fe[U][a];
    Lv[a] = w->fe[U][a] - w->gs[U][a];
    for (int b = 0; b < N; b++) if (b != a && ((U >> b) & 1)) Lv[a] += Wm[a][b];
  }
  return H;
}

static inline i64 eval(i64 H, const i64 *L, i64 M[N][N], int x) {
  i64 v = H;
  while (x) {
    int a = __builtin_ctz(x); x &= x - 1;
    v += L[a];
    int r = x;
    while (r) { int b = __builtin_ctz(r); r &= r - 1; v -= M[a][b]; }
  }
  return v;
}

/* direct DP on the 11-vertex tournament T + v_x: H and every deletion */
static i64 direct11(work *w, const m16 *t, int x, i64 *hdel) {
  const int n = N + 1, F = 1 << n;
  m16 o[N + 1];
  for (int a = 0; a < N; a++) o[a] = t[a];
  o[N] = 0;
  for (int a = 0; a < N; a++) { if ((x >> a) & 1) o[a] |= (m16)(1u << N); else o[N] |= (m16)(1u << a); }
  memset(w->f11, 0, sizeof(w->f11));
  for (int a = 0; a < n; a++) w->f11[1 << a][a] = 1;
  for (int S = 1; S < F; S++) for (int a = 0; a < n; a++) {
    i64 y = w->f11[S][a]; if (!y) continue;
    int q = o[a] & ~S;
    while (q) { int b = __builtin_ctz(q); q &= q - 1; w->f11[S | (1 << b)][b] += y; }
  }
  i64 h = 0; for (int a = 0; a < n; a++) h += w->f11[F - 1][a];
  for (int u = 0; u < n; u++) { int S = (F - 1) ^ (1 << u); i64 s = 0; for (int a = 0; a < n; a++) s += w->f11[S][a]; hdel[u] = s; }
  return h;
}

static int lambda11(const m16 *o, int n) {
  m16 in[16] = {0};
  for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) if ((o[u] >> v) & 1) in[v] |= (m16)(1u << u);
  int full = (1 << n) - 1, best = 1 << 30;
  for (int X = 1; X < full; X++) {
    int b = 0; int nX = full & ~X;
    for (int a = 0; a < n && b < best; a++) if ((X >> a) & 1) b += __builtin_popcount(in[a] & nX);
    if (b < best) best = b;
  }
  return best;
}

static void loadset(const char *fn, unsigned char *s) {
  FILE *f = fopen(fn, "r"); if (!f) { perror(fn); exit(1); }
  long v; while (fscanf(f, "%ld", &v) == 1) if (v >= 0 && v < HC) s[v] = 1;
  fclose(f);
}

static int cmp_u64(const void *a, const void *b) { u64 x = *(const u64 *)a, y = *(const u64 *)b; return x < y ? -1 : x > y; }


int main(int argc, char **argv) {
  if (argc < 5) { fprintf(stderr, "usage: aud_o11_cls CODES.bin S9.txt S10.txt OUT.txt [nthreads]\n"); return 1; }
  loadset(argv[2], inS9); loadset(argv[3], inS10);
#ifdef _OPENMP
  if (argc >= 6) omp_set_num_threads(atoi(argv[5]));
  int nthr = omp_get_max_threads();
#else
  int nthr = 1;
#endif
  FILE *fc = fopen(argv[1], "rb"); if (!fc) { perror(argv[1]); return 1; }
  fseek(fc, 0, SEEK_END); long nb = ftell(fc); fseek(fc, 0, SEEK_SET);
  long ncls = nb / 8; u64 *codes = malloc(nb);
  if (fread(codes, 8, ncls, fc) != (size_t)ncls) { fprintf(stderr, "read error\n"); return 1; }
  fclose(fc);
  work **WK = malloc(sizeof(work *) * nthr); for (int i = 0; i < nthr; i++) WK[i] = malloc(sizeof(work));
  long long nonfc = 0, fcext = 0, orphan_pairs = 0; time_t t0 = time(NULL);
  static unsigned char orphan_val[HC];
  FILE *fo = fopen(argv[4], "w");
  #pragma omp parallel for schedule(dynamic, 64) reduction(+:nonfc, fcext, orphan_pairs)
  for (long r = 0; r < ncls; r++) {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    work *w = WK[tid];
    m16 t[16], in[16] = {0};
    decode_msb(codes[r], N, t);
    for (int u = 0; u < N; u++) for (int v = 0; v < N; v++) if ((t[u] >> v) & 1) in[v] |= (m16)(1u << u);
    dp_tables(w, t, in);
    const int FULL = (1 << N) - 1;
    w->H = form_on(w, FULL, w->W, w->L);
    if (!inS9[w->H]) continue;                 /* only non-first-creator parents */
    nonfc++;
    for (int a = 0; a < N; a++) for (int b = 0; b < N; b++) w->M[a][b] = w->W[a][b] + w->W[b][a];
    for (int u = 0; u < N; u++) {
      i64 Wm[N][N];
      w->Hu[u] = form_on(w, FULL ^ (1 << u), Wm, w->Lu[u]);
      for (int a = 0; a < N; a++) for (int b = 0; b < N; b++) w->Mu[u][a][b] = Wm[a][b] + Wm[b][a];
      w->Lu[u][u] = 0;
    }
    for (int x = 0; x <= FULL; x++) {
      i64 hv = eval(w->H, w->L, w->M, x);
      if (inS10[hv]) continue;                  /* T' not a first creator */
      fcext++;
      int hasfc = 0;
      for (int u = 0; u < N && !hasfc; u++) { i64 hu = eval(w->Hu[u], w->Lu[u], w->Mu[u], x & ~(1 << u)); if (!inS9[hu]) hasfc = 1; }
      if (!hasfc) {
        orphan_pairs++;
        #pragma omp critical
        { orphan_val[hv] = 1; fprintf(fo, "ORPHAN parent %llu x %d H %lld\n", codes[r], x, hv); }
      }
    }
  }
  long nv = 0; for (int h = 0; h < HC; h++) if (orphan_val[h]) nv++;
  fprintf(fo, "class-level genealogy at order %d: non-first-creator parents %lld; (parent, pattern) pairs giving an order-%d first creator %lld; "
              "of these with no first-creator deletion (class-level orphans) %lld, values %ld; %lds, %d threads\n",
          N + 1, nonfc, N + 1, fcext, orphan_pairs, nv, (long)(time(NULL) - t0), nthr);
  fclose(fo);
  return 0;
}
