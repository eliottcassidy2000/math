/* aud_o11_check.c - independent auditor's checks of the order-11 witness files of first_creator_generation_20261108.
 * The block between the "copied" markers is the canonical labelling of the PREVIOUS independent audit
 * (first_creators_genealogy_20261107/audit/aud_enum.c), compiled here with MAXN 11; all other code is the auditor's.
 * Modes (all read the packet's files with the PACKET's bit order: bit b <-> pair (i<j) in lexicographic order,
 * b = 0 for (0,1); bit set = i -> j):
 *   wit  W.txt S9 S10 KLEIN      every line "h code x": H(T) not in S(9), H(T + v_x) = h by direct DP, h new at 11,
 *                                 and the h's are exactly klein's new values at 11, each once
 *   lam3 L.txt MINE.txt S9 S10   every pair of lambda3_witnesses_order11.txt (and the auditor's own list): lambda and H
 *                                 directly, canonical classes with explicit isomorphism certificates, |Aut|, SC,
 *                                 FC / EXT / module, C3[R5,R5,1]
 *   fert S10 code...             fertility of order-10 classes (packet bit order) by direct DP over all 1024 patterns,
 *                                 and their canonical codes in the auditor's (aud_enum) bit order
 *   root63 S7 S8                 order-9 first creators with a deletion of H = 63 also have another first-creator deletion
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define MAXN 11
typedef unsigned long long u64;
typedef unsigned short m16;
typedef long long i64;
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
/* ======================== end of copied block; the rest is the auditor's own code ===================== */
#define HC 100000
static unsigned char S7[HC], S8[HC], S9[HC], S10[HC], S11k[HC], new11[HC];

static void loadset(const char *fn, unsigned char *s) {
  FILE *f = fopen(fn, "r"); if (!f) { perror(fn); exit(1); }
  long v; while (fscanf(f, "%ld", &v) == 1) if (v >= 0 && v < HC) s[v] = 1;
  fclose(f);
}
static void decode_lsb(u64 c, int n, m16 *out) {   /* the packet's bit order */
  memset(out, 0, sizeof(m16) * MAXN); int b = 0;
  for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if ((c >> b) & 1) out[i] |= (m16)(1u << j); else out[j] |= (m16)(1u << i); b++; }
}
static void extend(const m16 *t, int n, int x, m16 *o) {   /* o = t + v_x, v = vertex n, bit a of x = [a -> v] */
  for (int a = 0; a < n; a++) o[a] = t[a];
  o[n] = 0;
  for (int a = 0; a < n; a++) { if ((x >> a) & 1) o[a] |= (m16)(1u << n); else o[n] |= (m16)(1u << a); }
}
static i64 dpf[1 << MAXN][MAXN];
static i64 hpath(const m16 *o, int n, i64 *dh) {   /* direct subset DP; dh[u] = H(T - u) */
  int F = 1 << n;
  memset(dpf, 0, sizeof(i64) * MAXN * (size_t)F);
  for (int a = 0; a < n; a++) dpf[1 << a][a] = 1;
  for (int S = 1; S < F; S++) for (int a = 0; a < n; a++) { i64 y = dpf[S][a]; if (!y) continue; int q = o[a] & ~S; while (q) { int b = __builtin_ctz(q); q &= q - 1; dpf[S | (1 << b)][b] += y; } }
  i64 h = 0; for (int a = 0; a < n; a++) h += dpf[F - 1][a];
  if (dh) for (int u = 0; u < n; u++) { int S = (F - 1) ^ (1 << u); i64 s = 0; for (int a = 0; a < n; a++) s += dpf[S][a]; dh[u] = s; }
  return h;
}
static int lambda_bf(const m16 *o, int n) {
  m16 in[MAXN] = {0};
  for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) if ((o[u] >> v) & 1) in[v] |= (m16)(1u << u);
  int full = (1 << n) - 1, best = 1 << 30;
  for (int X = 1; X < full; X++) { int b = 0, nX = full & ~X; for (int a = 0; a < n; a++) if ((X >> a) & 1) b += __builtin_popcount(in[a] & nX); if (b < best) best = b; }
  return best;
}
static int module_bf(const m16 *o, int n) {
  int full = (1 << n) - 1;
  for (int M = 1; M < full; M++) { if (__builtin_popcount(M) < 2) continue; int ok = 1;
    for (int z = 0; z < n && ok; z++) if (!((M >> z) & 1)) { int q = o[z] & M; if (q && q != M) ok = 0; }
    if (ok) return M; }
  return 0;
}
static int mindeg(const m16 *o, int n) { int m = n; for (int a = 0; a < n; a++) { int od = __builtin_popcount(o[a]); if (od < m) m = od; if (n - 1 - od < m) m = n - 1 - od; } return m; }

/* canonical code + |Aut| + self-converse + canonical labelling */
typedef struct { u64 code; long aut; int sc; int lab[MAXN]; } cinfo;
static cinfo cinf(const m16 *o, int n) {
  cinfo r; cctx C; canon(&C, o, n); r.code = C.best; r.aut = C.nbest; memcpy(r.lab, C.lab0, sizeof(int) * n);
  m16 op[MAXN]; converse(o, n, op); cctx C2; canon(&C2, op, n); r.sc = (C2.best == C.best);
  return r;
}
/* explicit isomorphism check: phi = lab2^{-1} o lab1 maps T1 to T2 iff codes equal; verify arc by arc */
static int iso_cert(const m16 *o1, const cinfo *c1, const m16 *o2, const cinfo *c2, int n) {
  int phi[MAXN];   /* lab[i] = vertex at canonical position i */
  for (int i = 0; i < n; i++) phi[c1->lab[i]] = c2->lab[i];
  for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) if (a != b) if (((o1[a] >> b) & 1) != ((o2[phi[a]] >> phi[b]) & 1)) return 0;
  return 1;
}

static int mode_wit(int argc, char **argv) {
  loadset(argv[3], S9); loadset(argv[4], S10);
  FILE *k = fopen(argv[5], "r"); char line[8192]; long nnew = 0;
  while (fgets(line, sizeof line, k)) { if (line[0] == '#') continue; long h; char nn[32]; if (sscanf(line, "%ld %31s", &h, nn) == 2 && nn[0] != '>' && atoi(nn) <= 11 && h < HC) { S11k[h] = 1; if (atoi(nn) == 11) { new11[h] = 1; nnew++; } } }
  fclose(k);
  long mism_new = 0; for (int h = 0; h < HC; h++) if (new11[h] != (S11k[h] && !S10[h])) mism_new++;
  FILE *f = fopen(argv[2], "r"); static int seen[HC];
  long lines = 0, badfc = 0, badh = 0, notnew = 0, dup = 0;
  while (fgets(line, sizeof line, f)) {
    if (line[0] == '#') continue;
    long h; u64 code; int x; if (sscanf(line, "%ld %llu %d", &h, &code, &x) != 3) continue;
    lines++;
    m16 t[MAXN], o[MAXN]; decode_lsb(code, 10, t);
    i64 H = hpath(t, 10, NULL);
    if (H >= HC || S9[H]) badfc++;
    extend(t, 10, x, o);
    i64 h2 = hpath(o, 11, NULL);
    if (h2 != h) badh++;
    if (h < 0 || h >= HC || !new11[h]) notnew++; else { if (seen[h]) dup++; seen[h]++; }
  }
  fclose(f);
  long missing = 0; for (int h = 0; h < HC; h++) if (new11[h] && !seen[h]) missing++;
  printf("witness check: klein new values at 11 (n(h) = 11): %ld (agree with S(11) minus the auditor's S(10): %s)\n", nnew, mism_new ? "NO" : "yes");
  printf("  lines %ld; parent not a first creator (H in S(9)): %ld; H(T + v_x) != h by direct DP: %ld; h not new at 11: %ld; duplicate h: %ld; new values without a line: %ld\n",
         lines, badfc, badh, notnew, dup, missing);
  return 0;
}

#define MAXW 512
static int mode_lam3(int argc, char **argv) {
  loadset(argv[4], S9); loadset(argv[5], S10);
  m16 T[MAXW][MAXN]; cinfo CI[MAXW]; i64 HH[MAXW]; int LAM[MAXW], nw = 0, src[MAXW];
  for (int file = 0; file < 2; file++) {
    FILE *f = fopen(argv[2 + file], "r"); if (!f) { perror(argv[2 + file]); return 1; }
    char line[8192]; long bad = 0, n0 = nw;
    while (fgets(line, sizeof line, f)) {
      if (line[0] == '#') continue;
      u64 code; int x; long h; char st[64];
      m16 t[MAXN];
      if (file == 0) { if (sscanf(line, "W %llu %d %ld %63s", &code, &x, &h, st) != 4) continue; decode_lsb(code, 10, t); }
      else { int lam; long hd; if (sscanf(line, "%llu %d %ld %d %ld %63s", &code, &x, &h, &lam, &hd, st) != 6) continue; decode(code, 10, t); }
      extend(t, 10, x, T[nw]);
      HH[nw] = hpath(T[nw], 11, NULL); LAM[nw] = lambda_bf(T[nw], 11);
      if (HH[nw] != h) bad++;
      if (LAM[nw] < 3) bad++;
      int instat = (strcmp(st, "inS10") == 0);
      if (instat != S10[HH[nw]]) bad++;
      CI[nw] = cinf(T[nw], 11); src[nw] = file; nw++;
    }
    fclose(f);
    printf("lam3 file %d (%s): pairs %ld, direct-check failures (H, lambda >= 3, status) %ld\n", file, file ? "auditor's own search" : "packet", nw - n0, bad);
  }
  /* classes */
  u64 cls[MAXW]; int ncls = 0, rep[MAXW], cnt[2][MAXW];
  memset(cnt, 0, sizeof cnt);
  long certfail = 0;
  for (int i = 0; i < nw; i++) {
    int k = -1; for (int j = 0; j < ncls; j++) if (cls[j] == CI[i].code) k = j;
    if (k < 0) { k = ncls++; cls[k] = CI[i].code; rep[k] = i; }
    else if (!iso_cert(T[i], &CI[i], T[rep[k]], &CI[rep[k]], 11)) certfail++;
    cnt[src[i]][k]++;
  }
  printf("canonical classes: %d; isomorphism certificates failed: %ld\n", ncls, certfail);
  printf("  class: H status |Aut| SC lambda mindeg module(mask) EXT(#first-creator deletions) deletionH | pairs(packet, auditor)\n");
  for (int k = 0; k < ncls; k++) {
    int i = rep[k]; i64 dh[MAXN]; hpath(T[i], 11, dh);
    int ext = 0; for (int u = 0; u < 11; u++) if (dh[u] < HC && !S9[dh[u]]) ext++;
    int mod = module_bf(T[i], 11);
    printf("  %llu: H %lld %s |Aut| %ld SC %d lambda %d mindeg %d module %d (size %d) EXT %d dels", cls[k], HH[i], S10[HH[i]] ? "inS10(non-FC)" : "new(FC)",
           CI[i].aut, CI[i].sc, LAM[i], mindeg(T[i], 11), mod, __builtin_popcount(mod), ext);
    for (int u = 0; u < 11; u++) printf(" %lld", dh[u]);
    printf(" | %d, %d\n", cnt[0][k], cnt[1][k]);
  }
  /* C3[R5, R5, 1]: A = 0..4, B = 5..9 regular (i -> i+1, i+2 mod 5), A => B => v => A, v = 10 */
  m16 c[MAXN] = {0};
  for (int blk = 0; blk < 2; blk++) for (int i = 0; i < 5; i++) for (int d = 1; d <= 2; d++) c[5 * blk + i] |= (m16)(1u << (5 * blk + (i + d) % 5));
  for (int i = 0; i < 5; i++) for (int j = 5; j < 10; j++) c[i] |= (m16)(1u << j);
  for (int j = 5; j < 10; j++) c[j] |= (m16)(1u << 10);
  for (int i = 0; i < 5; i++) c[10] |= (m16)(1u << i);
  cinfo cc = cinf(c, 11); i64 hc = hpath(c, 11, NULL);
  int match = -1; for (int k = 0; k < ncls; k++) if (cls[k] == cc.code) match = k;
  printf("C3[R5,R5,1]: H %lld, lambda %d, mindeg %d, |Aut| %ld, SC %d, canonical code %llu, equals class #%d%s\n", hc, lambda_bf(c, 11), mindeg(c, 11), cc.aut, cc.sc, cc.code,
         match, match >= 0 && iso_cert(c, &cc, T[rep[match]], &CI[rep[match]], 11) ? " (isomorphism certified)" : "");
  return 0;
}

static int mode_fert(int argc, char **argv) {
  loadset(argv[2], S10);
  for (int i = 3; i < argc; i++) {
    u64 code = strtoull(argv[i], 0, 10); m16 t[MAXN], o[MAXN]; decode_lsb(code, 10, t);
    i64 H = hpath(t, 10, NULL); long fert = 0;
    for (int x = 0; x < 1024; x++) { extend(t, 10, x, o); i64 h = hpath(o, 11, NULL); if (!S10[h]) fert++; }
    cinfo ci = cinf(t, 10);
    m16 op[MAXN]; converse(t, 10, op); cinfo cco = cinf(op, 10);
    printf("order-10 class (packet code %llu): H %lld, fertile patterns by direct DP %ld of 1024; auditor canonical code %llu (converse %llu), |Aut| %ld, SC %d\n",
           code, H, fert, ci.code, cco.code, ci.aut, ci.sc);
  }
  return 0;
}

/* canonical augmentation (same acceptance rule as aud_enum.c) to regenerate orders <= 9 */
static int cmpu(const void *a, const void *b) { u64 x = *(const u64 *)a, y = *(const u64 *)b; return x < y ? -1 : x > y; }
static int ext_rep(u64 rc, int n, u64 *res) {
  m16 base[MAXN]; decode(rc, n - 1, base); cctx C; int cnt = 0;
  for (int P = 0; P < (1 << (n - 1)); P++) {
    m16 t[MAXN]; for (int v = 0; v < n - 1; v++) t[v] = base[v]; t[n - 1] = 0;
    for (int v = 0; v < n - 1; v++) { if ((P >> v) & 1) t[n - 1] |= (m16)(1u << v); else t[v] |= (m16)(1u << (n - 1)); }
    canon(&C, t, n);
    if (uf_find(C.uf, n - 1) != uf_find(C.uf, C.lab0[n - 1])) continue;
    res[cnt++] = C.best;
  }
  qsort(res, cnt, sizeof(u64), cmpu); int m = 0;
  for (int i = 0; i < cnt; i++) if (i == 0 || res[i] != res[i - 1]) res[m++] = res[i];
  return m;
}
static int mode_root63(int argc, char **argv) {
  loadset(argv[2], S7); loadset(argv[3], S8);
  u64 *reps = malloc(8); size_t nr = 1; reps[0] = 0;
  for (int n = 2; n <= 9; n++) {
    u64 *nx = malloc(8 * (nr << (n - 1))); size_t m = 0; u64 tmp[1 << 10];
    for (size_t r = 0; r < nr; r++) { int c = ext_rep(reps[r], n, tmp); for (int i = 0; i < c; i++) nx[m++] = tmp[i]; }
    free(reps); reps = nx; nr = m;
  }
  long nfc = 0, with63 = 0, with63_other = 0, with63_only = 0;
  for (size_t r = 0; r < nr; r++) {
    m16 t[MAXN]; decode(reps[r], 9, t); i64 dh[MAXN]; i64 h = hpath(t, 9, dh);
    if (S8[h]) continue;
    nfc++;
    int has63 = 0, other = 0;
    for (int u = 0; u < 9; u++) { if (dh[u] == 63) has63 = 1; else if (!S7[dh[u]]) other = 1; }
    if (has63) { with63++; if (other) with63_other++; else with63_only++; }
  }
  printf("order 9: classes %zu, first creators %ld; with a deletion of H = 63 (a creator of 63): %ld; of these with another first-creator deletion %ld, without %ld\n",
         nr, nfc, with63, with63_other, with63_only);
  return 0;
}

/* lineage: rebuild the packet's lineage_le10.txt from the auditor's own classes (orders 4-9 regenerated here by canonical
   augmentation, order 10 read from the auditor's aud_o10_codes.bin in aud_enum bit order).
   usage: aud_o11_check lineage DIR_WITH_aud_S*.txt O10CODES.bin OUT.txt */
#define LW 64                 /* parent values are < 3358 < 64*64 */
static unsigned char Sx[11][HC];
typedef struct { long ncls[HC]; long npar[HC]; u64 par[HC][LW]; } linbuf;
static void lin_class(const m16 *t, int n, linbuf *L) {
  i64 dh[MAXN]; i64 h = hpath(t, n, dh);
  if (Sx[n - 1][h]) return;                      /* not a first creator */
  L->ncls[h]++;
  int any = 0;
  for (int u = 0; u < n; u++) if (!Sx[n - 2][dh[u]]) { any = 1; L->par[h][dh[u] >> 6] |= 1ULL << (dh[u] & 63); }
  if (any) L->npar[h]++;
}
static int mode_lineage(int argc, char **argv) {
  char fn[2048];
  for (int m = 1; m <= 10; m++) { snprintf(fn, sizeof fn, "%s/aud_S%d.txt", argv[2], m); loadset(fn, Sx[m]); }
  Sx[0][1] = 1;
  linbuf *L = calloc(1, sizeof(linbuf));
  u64 *reps = malloc(8); size_t nr = 1; reps[0] = 0;
  for (int n = 2; n <= 9; n++) {
    u64 *nx = malloc(8 * (nr << (n - 1))); size_t m = 0; u64 tmp[1 << 10];
    for (size_t r = 0; r < nr; r++) { int c = ext_rep(reps[r], n, tmp); for (int i = 0; i < c; i++) nx[m++] = tmp[i]; }
    free(reps); reps = nx; nr = m;
    if (n >= 4) for (size_t r = 0; r < nr; r++) { m16 t[MAXN]; decode(reps[r], n, t); lin_class(t, n, L); }
  }
  FILE *fc = fopen(argv[3], "rb"); if (!fc) { perror(argv[3]); return 1; }
  u64 code; long n10 = 0;
  while (fread(&code, 8, 1, fc) == 1) { m16 t[MAXN]; decode(code, 10, t); lin_class(t, 10, L); n10++; }
  fclose(fc);
  FILE *fo = fopen(argv[4], "w");
  fprintf(fo, "# auditor's rebuild of lineage_le10.txt: h n(h) #first-creator classes #with a first-creator parent parent values\n");
  long lines = 0;
  for (int n = 4; n <= 10; n++) for (int h = 0; h < HC; h++) if (L->ncls[h] && Sx[n][h] && !Sx[n - 1][h]) {
    fprintf(fo, "%d %d %ld %ld ", h, n, L->ncls[h], L->npar[h]);
    if (!L->npar[h]) fprintf(fo, "ORPHAN");
    int first = 1;
    for (int p = 0; p < 64 * LW; p++) if ((L->par[h][p >> 6] >> (p & 63)) & 1) { fprintf(fo, first ? "%d" : ",%d", p); first = 0; }
    fprintf(fo, "\n"); lines++;
  }
  fclose(fo);
  printf("lineage rebuilt: %ld order-10 classes read, %ld value lines written\n", n10, lines);
  return 0;
}

/* dump: the auditor's class codes (aud_enum bit order) of one order n <= 9, as raw u64, for positive controls */
static int mode_dump(int argc, char **argv) {
  int nt = atoi(argv[2]);
  u64 *reps = malloc(8); size_t nr = 1; reps[0] = 0;
  for (int n = 2; n <= nt; n++) {
    u64 *nx = malloc(8 * (nr << (n - 1))); size_t m = 0; u64 tmp[1 << 10];
    for (size_t r = 0; r < nr; r++) { int c = ext_rep(reps[r], n, tmp); for (int i = 0; i < c; i++) nx[m++] = tmp[i]; }
    free(reps); reps = nx; nr = m;
  }
  FILE *f = fopen(argv[3], "wb"); fwrite(reps, 8, nr, f); fclose(f);
  printf("order %d: %zu classes written\n", nt, nr);
  return 0;
}

int main(int argc, char **argv) {
  if (argc >= 4 && !strcmp(argv[1], "dump")) return mode_dump(argc, argv);
  if (argc >= 5 && !strcmp(argv[1], "lineage")) return mode_lineage(argc, argv);
  if (argc >= 6 && !strcmp(argv[1], "wit")) return mode_wit(argc, argv);
  if (argc >= 6 && !strcmp(argv[1], "lam3")) return mode_lam3(argc, argv);
  if (argc >= 4 && !strcmp(argv[1], "fert")) return mode_fert(argc, argv);
  if (argc >= 4 && !strcmp(argv[1], "root63")) return mode_root63(argc, argv);
  fprintf(stderr, "usage: see header\n"); return 1;
}
