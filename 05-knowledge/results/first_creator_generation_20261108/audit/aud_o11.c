/* aud_o11.c - independent auditor's exhaustive order-10 -> order-11 extension run (phase 2 of the audit of
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

#define N 10
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
  if (argc < 6) { fprintf(stderr, "usage: aud_o11 CODES.bin S9.txt S10.txt least_order_11.txt OUTDIR [nthreads] [stride]\n"); return 1; }
  loadset(argv[2], inS9); loadset(argv[3], inS10);
  { FILE *f = fopen(argv[4], "r"); char line[8192]; if (!f) { perror(argv[4]); return 1; }
    while (fgets(line, sizeof line, f)) { if (line[0] == '#') continue; long h; char nn[32];
      if (sscanf(line, "%ld %31s", &h, nn) == 2 && nn[0] != '>') { int k = atoi(nn); if (k <= 11 && h < HC) inS11k[h] = 1; } }
    fclose(f); }
  const char *dir = argv[5];
#ifdef _OPENMP
  if (argc >= 7) omp_set_num_threads(atoi(argv[6]));
  int nthr = omp_get_max_threads();
#else
  int nthr = 1;
#endif
  long stride = argc >= 8 ? atol(argv[7]) : 2500;
  FILE *fc = fopen(argv[1], "rb"); if (!fc) { perror(argv[1]); return 1; }
  fseek(fc, 0, SEEK_END); long nb = ftell(fc); fseek(fc, 0, SEEK_SET);
  long ncls = nb / 8; u64 *codes = malloc(nb);
  if (fread(codes, 8, ncls, fc) != (size_t)ncls) { fprintf(stderr, "read error\n"); return 1; }
  fclose(fc);
  qsort(codes, ncls, sizeof(u64), cmp_u64);   /* deterministic order (the file was written by threads) */
  char fn[1024];
  snprintf(fn, sizeof fn, "%s/aud_o11_lambda3.txt", dir); flam = fopen(fn, "w");
  fprintf(flam, "# auditor's exhaustive lambda>=3 search at order 11, H <= 15745: parent code (aud_enum bit order) x H lambda(direct) H(direct) status\n");
  snprintf(fn, sizeof fn, "%s/aud_o11_hypo.txt", dir); fhypo = fopen(fn, "w");
  fprintf(fhypo, "# hypo pairs found (parent code, x, H)\n");
  work **WK = malloc(sizeof(work *) * nthr); tstats **TS = malloc(sizeof(tstats *) * nthr);
  for (int i = 0; i < nthr; i++) { WK[i] = malloc(sizeof(work)); TS[i] = calloc(1, sizeof(tstats)); TS[i]->minfert = 1 << 30; TS[i]->lam3_min = 1LL << 60; }
  time_t t0 = time(NULL); long done = 0;
  #pragma omp parallel for schedule(dynamic, 64)
  for (long r = 0; r < ncls; r++) {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    work *w = WK[tid]; tstats *st = TS[tid];
    m16 t[16], in[16] = {0};
    decode_msb(codes[r], N, t);
    for (int u = 0; u < N; u++) for (int v = 0; v < N; v++) if ((t[u] >> v) & 1) in[v] |= (m16)(1u << u);
    int mindeg = N, req1 = 0, req0 = 0;
    for (int a = 0; a < N; a++) {
      int od = __builtin_popcount(t[a]), id = N - 1 - od;
      if (od < mindeg) mindeg = od; if (id < mindeg) mindeg = id;
      if (od == 2) req1 |= 1 << a;    /* needs a -> v to reach out-degree 3 */
      if (id == 2) req0 |= 1 << a;    /* needs v -> a to reach in-degree 3 */
    }
    dp_tables(w, t, in);
    const int FULL = (1 << N) - 1;
    w->H = form_on(w, FULL, w->W, w->L);
    for (int a = 0; a < N; a++) for (int b = 0; b < N; b++) w->M[a][b] = w->W[a][b] + w->W[b][a];
    int fcq = !inS9[w->H];
    st->ncls++;
    if (fcq) {
      st->nfc++;
      for (int u = 0; u < N; u++) {
        i64 Wm[N][N];
        w->Hu[u] = form_on(w, FULL ^ (1 << u), Wm, w->Lu[u]);
        for (int a = 0; a < N; a++) for (int b = 0; b < N; b++) w->Mu[u][a][b] = Wm[a][b] + Wm[b][a];
        w->Lu[u][u] = 0;
      }
    }
    int sample = (r % stride) == 0;
    long fert = 0, fert_bf = 0, hypo_bf = 0;
    for (int x = 0; x <= FULL; x++) {
      i64 hv = eval(w->H, w->L, w->M, x);
      if (hv < 0 || hv >= HC) { fprintf(stderr, "value out of range %lld\n", hv); exit(5); }
      st->A[hv]++;
      int ishypo = 0;
      if (fcq) {
        st->F[hv]++;
        if (!inS10[hv]) fert++;
        else if (x != 0 && x != FULL) {
          ishypo = 1;
          for (int u = 0; u < N && ishypo; u++) {
            i64 hu = eval(w->Hu[u], w->Lu[u], w->Mu[u], x & ~(1 << u));
            if (inS9[hu]) ishypo = 0;
          }
          if (ishypo) {
            st->hypo++;
            #pragma omp critical
            { fprintf(fhypo, "%llu %d %lld\n", codes[r], x, hv); }
          }
        }
      }
      if (mindeg >= 2) {
        int pc = __builtin_popcount(x);
        if ((x & req1) == req1 && (x & req0) == 0 && pc >= 3 && pc <= N - 3) {
          st->lam3_ext++;
          if (hv < st->lam3_min) { st->lam3_min = hv; st->lam3_min_code = codes[r]; st->lam3_min_x = x; }
          if (hv <= 15745) {
            m16 o[16]; for (int a = 0; a < N; a++) o[a] = t[a]; o[N] = 0;
            for (int a = 0; a < N; a++) { if ((x >> a) & 1) o[a] |= (m16)(1u << N); else o[N] |= (m16)(1u << a); }
            int lam = lambda11(o, N + 1);
            i64 hd[N + 1]; i64 hdir = direct11(w, t, x, hd);
            #pragma omp critical
            { fprintf(flam, "%llu %d %lld %d %lld %s\n", codes[r], x, hv, lam, hdir, inS10[hv] ? "inS10" : (inS11k[hv] ? "new" : "NOT-IN-S11?")); }
          }
        }
      }
      if (sample) {
        i64 hd[N + 1]; i64 hdir = direct11(w, t, x, hd);
        st->sample_checks++;
        if (hdir != hv) st->sample_fail++;
        if (fcq) {
          if (!inS10[hdir]) fert_bf++;
          /* every deletion u in V by the deletion form must match the direct DP */
          for (int u = 0; u < N; u++) { i64 hu = eval(w->Hu[u], w->Lu[u], w->Mu[u], x & ~(1 << u)); if (hu != hd[u]) st->sample_fail++; }
          if (hd[N] != w->H) st->sample_fail++;
          if (inS10[hdir]) { int all = 1; for (int u = 0; u <= N; u++) if (inS9[hd[u]]) all = 0; if (all) hypo_bf++; }
          if (ishypo != (inS10[hdir] && x != 0 && x != FULL && ({ int all = 1; for (int u = 0; u <= N; u++) if (inS9[hd[u]]) all = 0; all; }))) st->sample_fail++;
        }
      }
    }
    if (sample) { st->sample_cls++; if (fcq) { st->sample_fc++; if (fert_bf != fert) st->sample_fert_mismatch++; st->sample_hypo_bf += hypo_bf; } }
    if (fcq) {
      st->fert_hist[fert]++;
      if (fert == 0) st->sterile++;
      if (fert < 512) st->below512++;
      if (fert < st->minfert) { st->minfert = fert; st->nminfert = 0; st->minfert_H = w->H; }
      if (fert == st->minfert && st->nminfert < 16) st->minfert_codes[st->nminfert++] = codes[r];
    }
    if (mindeg >= 2) st->lam3_parents++;
    #pragma omp atomic
    done++;
    if ((done % 500000) == 0) { fprintf(stderr, "  %ld / %ld classes [%lds]\n", done, ncls, (long)(time(NULL) - t0)); }
  }
  /* merge */
  tstats *S = TS[0];
  for (int i = 1; i < nthr; i++) {
    tstats *B = TS[i];
    for (int h = 0; h < HC; h++) { S->A[h] += B->A[h]; S->F[h] += B->F[h]; }
    S->ncls += B->ncls; S->nfc += B->nfc; S->sterile += B->sterile; S->below512 += B->below512; S->hypo += B->hypo;
    for (int k = 0; k <= 1024; k++) S->fert_hist[k] += B->fert_hist[k];
    if (B->minfert < S->minfert) { S->minfert = B->minfert; S->nminfert = 0; S->minfert_H = B->minfert_H; }
    if (B->minfert == S->minfert) for (int j = 0; j < B->nminfert && S->nminfert < 16; j++) S->minfert_codes[S->nminfert++] = B->minfert_codes[j];
    S->lam3_parents += B->lam3_parents; S->lam3_ext += B->lam3_ext;
    if (B->lam3_min < S->lam3_min) { S->lam3_min = B->lam3_min; S->lam3_min_code = B->lam3_min_code; S->lam3_min_x = B->lam3_min_x; }
    S->sample_cls += B->sample_cls; S->sample_fc += B->sample_fc; S->sample_checks += B->sample_checks; S->sample_fail += B->sample_fail;
    S->sample_fert_mismatch += B->sample_fert_mismatch; S->sample_hypo_bf += B->sample_hypo_bf;
  }
  snprintf(fn, sizeof fn, "%s/aud_o11.out", dir); FILE *fo = fopen(fn, "w");
  fprintf(fo, "auditor's order 10 -> 11 run: classes %ld (first creators %ld), extensions %lld, %d threads, %lds\n", S->ncls, S->nfc,
          (long long)S->ncls * 1024LL, nthr, (long)(time(NULL) - t0));
  long s11 = 0, s11k = 0, miss = 0, extra = 0, newv = 0, unreached = 0, onlyfc = 0;
  for (int h = 0; h < HC; h++) {
    if (S->A[h]) s11++;
    if (inS11k[h]) s11k++;
    if (inS11k[h] && !S->A[h]) miss++;
    if (S->A[h] && !inS11k[h]) extra++;
  }
  fprintf(fo, "S(11) from the extensions: %ld values; klein's S(11): %ld; missing %ld, extra %ld\n", s11, s11k, miss, extra);
  /* new values at 11 */
  static u64 fl[HC]; long nf = 0;
  for (int h = 0; h < HC; h++) if (S->A[h] && !inS10[h]) {
    newv++;
    if (!S->F[h]) unreached++;
    if (S->F[h] == S->A[h]) onlyfc++;
    fl[nf++] = S->F[h];
  }
  qsort(fl, nf, sizeof(u64), cmp_u64);
  fprintf(fo, "R1: new values at 11: %ld; with no first-creator parent: %ld\n", newv, unreached);
  fprintf(fo, "R5: new values all of whose (class, pattern) realizations have a first-creator parent: %ld (%.4f)\n", onlyfc, (double)onlyfc / newv);
  fprintf(fo, "first-creator realizations per new value: min %llu, median %llu (lower) / %llu (upper), max %llu\n", fl[0], fl[(nf - 1) / 2], fl[nf / 2], fl[nf - 1]);
  fprintf(fo, "new values with at most 9 first-creator realizations (h: F, A):");
  for (int h = 0; h < HC; h++) if (S->A[h] && !inS10[h] && S->F[h] <= 9) fprintf(fo, " %d:%llu,%llu", h, S->F[h], S->A[h]);
  fprintf(fo, "\n");
  fprintf(fo, "R2: hypo pairs %ld\n", S->hypo);
  fprintf(fo, "R3: sterile first creators %ld\n", S->sterile);
  fprintf(fo, "R4: first creators with fewer than 512 fertile patterns %ld; minimum fertile %ld (H %lld), attained by:", S->below512, S->minfert, S->minfert_H);
  for (int j = 0; j < S->nminfert; j++) fprintf(fo, " %llu", S->minfert_codes[j]);
  fprintf(fo, "\nfertile-count distribution (lowest values):");
  { int shown = 0; for (int k = 0; k <= 1024 && shown < 12; k++) if (S->fert_hist[k]) { fprintf(fo, " %d:%ld", k, S->fert_hist[k]); shown++; } }
  fprintf(fo, "\nlambda>=3 search: parents with min degree >= 2: %lld; degree-valid patterns %lld; min H over all of them %lld (parent %llu, x %d)\n",
          S->lam3_parents, S->lam3_ext, S->lam3_min, S->lam3_min_code, S->lam3_min_x);
  fprintf(fo, "sample (every %ld-th class): classes %ld (first creators %ld), direct-DP pattern checks %ld, mismatches %ld, fertile-count mismatches %ld, brute-force hypo pairs %ld\n",
          stride, S->sample_cls, S->sample_fc, S->sample_checks, S->sample_fail, S->sample_fert_mismatch, S->sample_hypo_bf);
  fclose(fo); fclose(flam); fclose(fhypo);
  return 0;
}
