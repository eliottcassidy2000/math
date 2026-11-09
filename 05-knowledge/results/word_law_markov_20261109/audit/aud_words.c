/* Independent audit code for THM-5020 (opus-2026-11-09-S88), written by the auditor.
 *
 * Convention (THM-4951): vertices 0..n-1, base path i -> i+1, tile (i, j) with j >= i + 2,
 * reversed in t iff j -> i; u = t + 1 reverses every tile.  A word w of length L = n - 2;
 * A = {a : w_a = 1}, B = {a : w_a = 0}; P(X) = max 3-separated subset of X, computed here by
 * the DP f(i) = max(f(i-1), [i-1 in X](1 + f(i-3))) (NOT by the greedy rule).
 * R(w) reverses exactly the tiles inside a window [p, q + 2] of some 1-run [p, q].
 *
 * Modes:
 *   validate N      : all lines {t, t+1} for n = N (N <= 8): census by unordered nu-pair and colour,
 *                     to compare with THM-4951 (4) (validates the nu solver).
 *   words NMIN NMAX : every word, n = NMIN..NMAX.  Checks
 *       (J)   P(A) + P(B) >= ceil((L + 2J)/3), and the deletion step for every single junction;
 *       (2)   nu(R(w) + 1) = P(B);
 *       (3)   nu(R(w)) = P(A) when w has no 0-junction (no factor 11011);
 *       proof lemmas: every 3-cycle of u has the orientation a->b->c->a and one of the types T1-T4
 *       exactly as stated, its certificate x(D) is a 0-position with x(D) + 1 in D, and disjoint
 *       3-cycles have certificates >= 3 apart; every 3-cycle of t is a window triangle or a
 *       shared-vertex triangle as stated, and (no 0-junction) meets F = {g_i + 2} (greedy packing g of A).
 *       Also records how often R(w) fails on words WITH a 0-junction.
 *   tight NMIN NMAX : every tight colouring (pruned DFS), n = NMIN..NMAX: J <= 1, the THM-5020 (4)
 *       realisation (R(w) if no 11011, else R(complement)) is a tight line with an endpoint of triple
 *       word w, and it is blue iff w is a palindrome.  Counts tight colourings, palindromic ones,
 *       and those with a 0-junction / 1-junction / both; records whether R(w) itself fails.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 40
typedef uint64_t u64;

static int n;
static u64 outm[MAXN];
static unsigned char rv[MAXN][MAXN]; /* rv[i][j] = 1 iff tile (i,j) reversed in t */

static void build(int comp) {
  for (int i = 0; i < n; i++) outm[i] = 0;
  for (int i = 0; i + 1 < n; i++) outm[i] |= 1ULL << (i + 1);
  for (int i = 0; i < n; i++)
    for (int j = i + 2; j < n; j++) {
      int r = rv[i][j] ^ comp;
      if (r) outm[j] |= 1ULL << i; else outm[i] |= 1ULL << j;
    }
}
static int arc(int x, int y) { return (outm[x] >> y) & 1; }

/* ---------- exact nu: max number of vertex-disjoint 3-cycles ---------- */
#define MAXT 20000
static u64 tri_by[MAXN][MAXT / 4];
static int ntri_by[MAXN];
static int best;
static long long nodes;
static void rec(u64 avail, int cur) {
  nodes++;
  if (cur > best) best = cur;
  if (!avail) return;
  if (cur + __builtin_popcountll(avail) / 3 <= best) return;
  int v = __builtin_ctzll(avail);
  for (int t = 0; t < ntri_by[v]; t++) {
    u64 T = tri_by[v][t];
    if ((T & avail) == T) rec(avail & ~T, cur + 1);
  }
  rec(avail & ~(1ULL << v), cur);
}
static int nu_cur(void) {
  for (int v = 0; v < n; v++) ntri_by[v] = 0;
  for (int a = 0; a < n; a++)
    for (int b = a + 1; b < n; b++)
      for (int c = b + 1; c < n; c++) {
        int cyc = (arc(a, b) && arc(b, c) && arc(c, a)) || (arc(a, c) && arc(c, b) && arc(b, a));
        if (cyc) tri_by[a][ntri_by[a]++] = (1ULL << a) | (1ULL << b) | (1ULL << c);
      }
  best = 0;
  rec((n >= 64) ? ~0ULL : ((1ULL << n) - 1), 0);
  return best;
}

/* ---------- words ---------- */
static int L;
static int w[MAXN];

static int Pcol(int col) { /* DP over positions */
  int f[MAXN + 3];
  f[0] = f[1] = f[2] = 0; /* f[i+3] = packing within positions < i */
  for (int i = 0; i < L; i++) {
    int a = f[i + 2];                    /* f(i)   : positions < i   -> index i+3-1 */
    int b = (w[i] == col) ? 1 + f[i] : 0; /* positions < i-2 -> index i */
    f[i + 3] = a > b ? a : b;
  }
  return f[L + 2];
}
static int Pcol_of(const int *x, int len, int col) {
  int f[MAXN + 3];
  f[0] = f[1] = f[2] = 0;
  for (int i = 0; i < len; i++) {
    int a = f[i + 2];
    int b = (x[i] == col) ? 1 + f[i] : 0;
    f[i + 3] = a > b ? a : b;
  }
  return f[len + 2];
}
/* junction centred at z of type col = letter at z (0-junction: w = 11011 around z) */
static int is_junction(const int *x, int len, int z) {
  if (z < 2 || z + 2 >= len) return 0;
  int c = x[z];
  return x[z - 2] != c && x[z - 1] != c && x[z + 1] != c && x[z + 2] != c;
}
static int count_j(const int *x, int len, int col) {
  int s = 0;
  for (int z = 2; z + 2 < len; z++) if (x[z] == col && is_junction(x, len, z)) s++;
  return s;
}

/* run structure of 1-runs: run_start[a], run_end[a] for a in a 1-run */
static void set_runwindow(const int *x) {
  for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) rv[i][j] = 0;
  int a = 0;
  while (a < L) {
    if (x[a] != 1) { a++; continue; }
    int p = a;
    while (a < L && x[a] == 1) a++;
    int q = a - 1;
    for (int i = p; i <= q + 2; i++)
      for (int j = i + 2; j <= q + 2; j++) rv[i][j] = 1;
  }
}

static long long fail_class = 0;
static long long n_tri_u = 0, n_pairs_u = 0, n_tri_t = 0;
static int mutate = 0; /* 1: deliberately wrong T2 certificate (x = p - 2) to test that the checker bites */
static void report(const char *what, int nn) {
  fail_class++;
  if (fail_class <= 20) {
    printf("  FAIL %s n=%d w=", what, nn);
    for (int i = 0; i < L; i++) putchar('0' + w[i]);
    putchar('\n');
  }
}

/* is position p the start / end of a 1-run */
static int is_start(int p) { return p >= 0 && p < L && w[p] == 1 && (p == 0 || w[p - 1] == 0); }
static int is_end(int q) { return q >= 0 && q < L && w[q] == 1 && (q == L - 1 || w[q + 1] == 0); }
static int run_start_of(int a) { while (a > 0 && w[a - 1] == 1) a--; return a; }
static int run_end_of(int a) { while (a < L - 1 && w[a + 1] == 1) a++; return a; }

/* certificate lemma on u = R(w) + 1; returns number of violations */
static int check_u_lemma(void) {
  static int ta[MAXT], tb[MAXT], tc[MAXT], tx[MAXT];
  int nt = 0, bad = 0;
  for (int a = 0; a < n; a++)
    for (int b = a + 1; b < n; b++)
      for (int c = b + 1; c < n; c++) {
        int o1 = arc(a, b) && arc(b, c) && arc(c, a);
        int o2 = arc(a, c) && arc(c, b) && arc(b, a);
        if (!o1 && !o2) continue;
        if (o2) { bad++; continue; } /* a->c->b->a must be impossible in u */
        int x = -100;
        if (b == a + 1 && c == b + 1) { /* T1 */
          if (a < L && w[a] == 0) x = a; else bad++;
        } else if (b == a + 1 && c >= b + 2) { /* T2: b = p start of 1-run, c in window, c >= p+2 */
          int p = b;
          if (is_start(p) && p >= 1 && c <= run_end_of(p) + 2) x = mutate ? p - 2 : p - 1; else bad++;
        } else if (b >= a + 2 && c == b + 1) { /* T3: b = q+2, q end of 1-run [p,q], a in [p,q] */
          int q = b - 2;
          if (is_end(q) && a >= run_start_of(q) && a <= q && q + 3 <= n - 1) x = q + 1; else bad++;
        } else { /* T4: b shared by windows [p'', b] and [b, e'] */
          int ok = is_end(b - 2) && b - 1 >= 0 && w[b - 1] == 0 && is_start(b) &&
                   a >= run_start_of(b - 2) && a <= b - 2 && c >= b + 2 && c <= run_end_of(b) + 2;
          if (ok) x = b - 1; else bad++;
        }
        if (x == -100) continue;
        if (!(x >= 0 && x < L && w[x] == 0)) bad++;
        if (!(x + 1 == a || x + 1 == b || x + 1 == c)) bad++;
        ta[nt] = a; tb[nt] = b; tc[nt] = c; tx[nt] = x; nt++; n_tri_u++;
      }
  for (int i = 0; i < nt; i++) {
    u64 Di = (1ULL << ta[i]) | (1ULL << tb[i]) | (1ULL << tc[i]);
    for (int j = i + 1; j < nt; j++) {
      u64 Dj = (1ULL << ta[j]) | (1ULL << tb[j]) | (1ULL << tc[j]);
      if (Di & Dj) continue;
      n_pairs_u++;
      int d = tx[i] - tx[j];
      if (d < 0) d = -d;
      if (d < 3) bad++;
    }
  }
  return bad;
}

/* t = R(w): triangle classification, and (no 0-junction) F = {g_i + 2} meets all 3-cycles */
static int check_t_lemma(int J0) {
  int g[MAXN], ng = 0, last = -10;
  for (int a = 0; a < L; a++) if (w[a] == 1 && a - last >= 3) { g[ng++] = a; last = a; }
  u64 F = 0;
  for (int i = 0; i < ng; i++) F |= 1ULL << (g[i] + 2);
  int bad = 0;
  for (int a = 0; a < n; a++)
    for (int b = a + 1; b < n; b++)
      for (int c = b + 1; c < n; c++) {
        int o1 = arc(a, b) && arc(b, c) && arc(c, a);
        int o2 = arc(a, c) && arc(c, b) && arc(b, a);
        if (!o1 && !o2) continue;
        n_tri_t++;
        if (o1) {
          if (!(b == a + 1 && c == a + 2 && a < L && w[a] == 1)) bad++;
        } else {
          int ok = is_end(b - 2) && w[b - 1] == 0 && is_start(b) && a >= run_start_of(b - 2) && a <= b - 2 &&
                   c >= b + 2 && c <= run_end_of(b) + 2;
          if (!ok) bad++;
        }
        if (J0 == 0) {
          u64 D = (1ULL << a) | (1ULL << b) | (1ULL << c);
          if (!(D & F)) bad++;
        }
      }
  if (J0 == 0 && ng != Pcol(1)) bad++;
  return bad;
}

static int is_pal(void) { for (int i = 0; i < L; i++) if (w[i] != w[L - 1 - i]) return 0; return 1; }
static int sigma_fixed(void) {
  for (int i = 0; i < n; i++) for (int j = i + 2; j < n; j++) if (rv[i][j] != rv[n - 1 - j][n - 1 - i]) return 0;
  return 1;
}

static void mode_words(int nmin, int nmax) {
  for (n = nmin; n <= nmax; n++) {
    L = n - 2;
    int k = n / 3;
    long long cnt = 0, bJ = 0, bDel = 0, bDelJ = 0, b2 = 0, b3 = 0, bU = 0, bT = 0, nJ0 = 0, failJ0 = 0, tightJ0 = 0,
              tightJ0fail = 0;
    int ex_print = 0;
    for (u64 code = 0; code < (1ULL << L); code++) {
      for (int i = 0; i < L; i++) w[i] = (code >> (L - 1 - i)) & 1;
      cnt++;
      int pa = Pcol(1), pb = Pcol(0);
      int J0 = count_j(w, L, 0), J1 = count_j(w, L, 1), J = J0 + J1;
      if (3 * (pa + pb) < L + 2 * J) { bJ++; report("junction-bound", n); }
      /* deletion step for each junction z */
      for (int z = 2; z + 2 < L; z++) {
        if (!is_junction(w, L, z)) continue;
        int x[MAXN], m = 0;
        for (int i = 0; i < L; i++) if (i != z) x[m++] = w[i];
        int qa = Pcol_of(x, m, 1), qb = Pcol_of(x, m, 0);
        if (pa + pb < qa + qb + 1) { bDel++; report("deletion", n); }
        int Jx = count_j(x, m, 0) + count_j(x, m, 1);
        if (Jx != J - 1) { bDelJ++; report("deletion-junction-count", n); }
      }
      set_runwindow(w);
      build(0);
      int nt = nu_cur();
      int tb = check_t_lemma(J0);
      build(1);
      int nuu = nu_cur();
      int ub = check_u_lemma();
      if (nuu != pb) { b2++; report("(2)", n); }
      if (ub) { bU++; report("u-lemma", n); }
      if (tb) { bT++; report("t-lemma", n); }
      if (J0 == 0 && nt != pa) { b3++; report("(3)", n); }
      if (J0 > 0) {
        nJ0++;
        if (nt != pa) failJ0++;
        if (pa + pb == k) { tightJ0++; if (nt != pa) tightJ0fail++; }
        if (nt != pa && ex_print < 2 && n <= 9) {
          printf("  n=%d word with 0-junction where R(w) fails: ", n);
          for (int i = 0; i < L; i++) putchar('0' + w[i]);
          printf(" nu(R(w)) = %d, P(A) = %d\n", nt, pa);
          ex_print++;
        }
      }
    }
    printf("n=%d L=%d words=%lld: junction-bound fails %lld, deletion-step fails %lld, deletion junction-count "
           "mismatches %lld, (2) fails %lld, (3) fails %lld, u-lemma violations %lld, t-lemma violations %lld | "
           "words with a 0-junction %lld, of which R(w) has nu != P(A): %lld; tight ones %lld, failing %lld\n",
           n, L, cnt, bJ, bDel, bDelJ, b2, b3, bU, bT, nJ0, failJ0, tightJ0, tightJ0fail);
    fflush(stdout);
  }
}

/* ---------- tight colourings ---------- */
static long long t_cnt, t_pal, t_j0, t_j1, t_both, t_fail, t_rw_fail, t_j0_x0, t_pal_j;
static int K;
static void tight_leaf(void) {
  int pa = Pcol(1), pb = Pcol(0);
  if (pa + pb != K) return;
  t_cnt++;
  int pal = is_pal();
  t_pal += pal;
  int J0 = count_j(w, L, 0), J1 = count_j(w, L, 1);
  if (J0) t_j0++;
  if (J1) t_j1++;
  if (J0 && J1) t_both++;
  if (J0 && w[0] == 1) t_j0_x0++;
  if (pal && (J0 + J1)) t_pal_j++;
  if (J0 + J1 > 1) { t_fail++; report("tight with J > 1", n); return; }
  /* realisation */
  int wb[MAXN];
  for (int i = 0; i < L; i++) wb[i] = 1 - w[i];
  const int *base = (J0 == 0) ? w : wb;
  set_runwindow(base);
  build(0);
  int a = nu_cur();
  build(1);
  int b = nu_cur();
  if (a + b != K) { t_fail++; report("realisation not tight", n); }
  /* endpoint words: t has triple word = base, t+1 has complement of base */
  int ok_word = 1;
  for (int i = 0; i < L; i++) {
    int tw = rv[i][i + 2];
    int want = (J0 == 0) ? w[i] : 1 - w[i]; /* endpoint t = R(base) has word base; t+1 has word w when base = wb */
    if (tw != want) ok_word = 0;
  }
  if (!ok_word) { t_fail++; report("word not realised", n); }
  if (J0 == 0) { if (a != pa || b != pb) { t_fail++; report("pair mismatch", n); } }
  else { if (a != pb || b != pa) { t_fail++; report("pair mismatch (complement)", n); } }
  if (sigma_fixed() != pal) { t_fail++; report("colour", n); }
  /* does R(w) itself work? */
  if (J0) {
    set_runwindow(w);
    build(0);
    int a2 = nu_cur();
    if (a2 != pa) t_rw_fail++;
  }
}
static void tight_dfs(int pos) {
  if (pos >= 1) {
    /* prune: packings of the prefix are lower bounds for the full word */
    int save = L;
    L = pos;
    int s = Pcol(1) + Pcol(0);
    L = save;
    if (s > K) return;
  }
  if (pos == L) { tight_leaf(); return; }
  w[pos] = 0; tight_dfs(pos + 1);
  w[pos] = 1; tight_dfs(pos + 1);
}
static void mode_tight(int nmin, int nmax) {
  for (n = nmin; n <= nmax; n++) {
    L = n - 2;
    K = n / 3;
    t_cnt = t_pal = t_j0 = t_j1 = t_both = t_fail = t_rw_fail = t_j0_x0 = t_pal_j = 0;
    tight_dfs(0);
    printf("n=%d k=%d: tight colourings %lld (palindromic %lld); with 0-junction %lld (x_0 = 1: %lld), with "
           "1-junction %lld, both %lld; palindromic with a junction %lld; realisation failures %lld; R(w) alone "
           "fails on %lld of the 0-junction words\n",
           n, K, t_cnt, t_pal, t_j0, t_j0_x0, t_j1, t_both, t_pal_j, t_fail, t_rw_fail);
    fflush(stdout);
  }
}

/* ---------- validation: all lines at small n ---------- */
static void mode_validate(int N) {
  n = N;
  int m = 0, ti[64], tj[64];
  for (int i = 0; i < n; i++) for (int j = i + 2; j < n; j++) { ti[m] = i; tj[m] = j; m++; }
  int k = n / 3;
  static long long black[12][12], blue[12][12];
  memset(black, 0, sizeof black);
  memset(blue, 0, sizeof blue);
  /* lines {t, t+1}: take t with tile 0 = (0,2) not reversed... any representative: t with bit (m-1) = 0 */
  for (u64 code = 0; code < (1ULL << (m - 1)); code++) {
    for (int e = 0; e < m; e++) rv[ti[e]][tj[e]] = (code >> e) & 1;
    build(0);
    int a = nu_cur();
    build(1);
    int b = nu_cur();
    int lo = a < b ? a : b, hi = a < b ? b : a;
    if (sigma_fixed()) blue[lo][hi]++; else black[lo][hi]++;
  }
  printf("all lines at n = %d (black + blue) by unordered nu-pair:\n", n);
  for (int i = 0; i <= k; i++)
    for (int j = i; j <= k; j++)
      if (black[i][j] + blue[i][j]) printf("  (%d,%d): %lld + %lld\n", i, j, black[i][j], blue[i][j]);
}

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage\n"); return 1; }
  if (getenv("AUDIT_MUTATE")) mutate = 1;
  if (!strcmp(argv[1], "validate")) mode_validate(atoi(argv[2]));
  else if (!strcmp(argv[1], "words")) mode_words(atoi(argv[2]), atoi(argv[3]));
  else if (!strcmp(argv[1], "tight")) mode_tight(atoi(argv[2]), atoi(argv[3]));
  printf("total reported failures: %lld (3-cycles of u classified %lld, disjoint pairs of them checked %lld, "
         "3-cycles of t classified %lld)%s\n",
         fail_class, n_tri_u, n_pairs_u, n_tri_t, mutate ? " [MUTATION TEST: wrong T2 certificate]" : "");
  return 0;
}
