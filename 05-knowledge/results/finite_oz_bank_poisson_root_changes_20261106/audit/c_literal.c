/* c_literal.c - auditor's literal-tournament check of the source moves and of the O/Z tail formula.
 * Construction (Session 198 bilinear-descent section 2, Session 206 ROOT_OBSTRUCTION section 1):
 *   rooted edge t: vertices 0 (root c) and 1, arc 0 -> 1.
 *   O: append a universal sink w (every existing vertex beats w); keep the root.
 *   Z: append w with w -> root and x -> w for every other existing x; keep the root.
 *   J: append w as for Z, then make w the root.
 * For each tournament we count by subset DP (no formulas):
 *   D = #Hamiltonian paths, E = #Hamiltonian paths ending at the root, d = #Hamiltonian paths of V - root,
 *   U = #ordered covers (P1, P2) by two vertex-disjoint directed paths, P2 nonempty and starting at the root,
 *       P1 possibly empty (H(empty) = 1),
 *   S = #Hamiltonian paths starting at the root (must be 1: unit start).
 * Mode W L: every word over {O,Z,J} of length <= L: compare (D,E,U,d) with the move arithmetic, and S == 1.
 *          Also JZ = OJ on the literal responses for every prefix of length <= L-2.
 * Mode T L K: every prefix word of length <= L and every O/Z tail of length <= K (parameter m < 2^K, letters
 *          chronological = binary digits low to high, O = 0, Z = 1): compare literal D with
 *          h + R m - d Q(m) - e (m mod 2), (h,R,e,d) = (D,U,E,d) of the prefix, Q by its definition.
 * Mode S file K: for each seed word in file (one word per line, letters OZJ, '-' for the empty word): same as T.
 * Mode R K: root symmetry: literal H of the converse tournament with reversed chain = H, and T-pairing.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef long long ll; typedef unsigned long long u64;
#define MAXN 19
static int n, root; static uint32_t beats[MAXN]; /* beats[i] = set of j with i -> j */
static u64 *f;   /* f[S*n+v]: #paths covering S ending at v (any start) */
static u64 *g;   /* g[S*n+v]: #paths covering S starting at root, ending at v */
static u64 *Hs;  /* H(S) */
static void init_t(void){ n = 2; root = 0; memset(beats, 0, sizeof beats); beats[0] = 1u << 1; }
static void move(char c){
  if (n >= MAXN){ fprintf(stderr, "too many vertices\n"); exit(2); }
  int w = n++; beats[w] = 0;
  for (int x = 0; x < w; x++) beats[x] |= 1u << w;            /* everyone beats w ... */
  if (c == 'Z' || c == 'J'){ beats[root] &= ~(1u << w); beats[w] |= 1u << root; } /* ... except the root */
  if (c == 'J') root = w;
}
static void count(ll *D, ll *E, ll *U, ll *d, ll *S){
  u64 full = (1ULL << n) - 1, NS = 1ULL << n;
  memset(f, 0, NS * n * sizeof(u64)); memset(g, 0, NS * n * sizeof(u64));
  for (int v = 0; v < n; v++) f[(1ULL << v) * n + v] = 1;
  g[(1ULL << root) * n + root] = 1;
  for (u64 Sm = 1; Sm < NS; Sm++){
    for (int v = 0; v < n; v++){
      if (!((Sm >> v) & 1)) continue;
      u64 fv = f[Sm * n + v], gv = g[Sm * n + v];
      if (!fv && !gv) continue;
      uint32_t out = beats[v] & ~(uint32_t)Sm;
      while (out){ int w = __builtin_ctz(out); out &= out - 1;
        f[(Sm | (1ULL << w)) * n + w] += fv; g[(Sm | (1ULL << w)) * n + w] += gv; }
    }
  }
  Hs[0] = 1;
  for (u64 Sm = 1; Sm < NS; Sm++){ u64 s = 0; for (int v = 0; v < n; v++) s += f[Sm * n + v]; Hs[Sm] = s; }
  *D = (ll)Hs[full]; *E = (ll)f[full * n + root]; *d = (ll)Hs[full & ~(1ULL << root)];
  u64 Ssum = 0; for (int v = 0; v < n; v++) Ssum += g[full * n + v]; *S = (ll)Ssum;
  u64 Uc = 0;
  for (u64 Sm = 1; Sm < NS; Sm++){ if (!((Sm >> root) & 1)) continue;
    u64 G = 0; for (int v = 0; v < n; v++) G += g[Sm * n + v];
    if (G) Uc += G * Hs[full & ~Sm]; }
  *U = (ll)Uc;
}
static ll Qdef(ll m){ /* Session 205 pair-sum definition, independent of the recurrence */
  int L = 0; while ((m >> L)) L++;
  ll q = 0;
  for (int i = 0; i < L; i++) for (int j = 0; j < i; j++) if (((m >> i) & 1) && ((m >> j) & 1)) q += 1LL << (i - 1 - j);
  for (int i = 1; i < L; i++) if (((m >> i) & 1) && ((m >> (i - 1)) & 1)) q += 1;
  return q;
}
static void arith(char c, ll *D, ll *E, ll *U, ll *d){
  ll a = *D, b = *E, u = *U, e = *d;
  if (c == 'O'){ *D = a; *E = 0; *U = 2 * u; *d = e; }
  else if (c == 'Z'){ *D = a - b + u; *E = e; *U = 2 * u - e; *d = e; }
  else { *D = a - b + u; *E = a - b; *U = a + u; *d = a; }
}
static long long nchecked = 0, nbad = 0, nunit_bad = 0, njz = 0, njz_bad = 0, ntail = 0, ntail_bad = 0;
static char word[64];
static void build(const char *w){ init_t(); for (const char *p = w; *p; p++) move(*p); }
static void check_word(const char *w){
  ll D, E, U, d, S; build(w); count(&D, &E, &U, &d, &S);
  ll aD = 1, aE = 0, aU = 2, ad = 1; for (const char *p = w; *p; p++) arith(*p, &aD, &aE, &aU, &ad);
  nchecked++;
  if (D != aD || E != aE || U != aU || d != ad){ nbad++; if (nbad < 10) printf("MISMATCH word %s literal (%lld,%lld,%lld,%lld) arith (%lld,%lld,%lld,%lld)\n", w[0] ? w : "-", D, E, U, d, aD, aE, aU, ad); }
  if (S != 1){ nunit_bad++; if (nunit_bad < 10) printf("NOT UNIT-START word %s S=%lld\n", w, S); }
}
static void rec_words(int len, int L){
  word[len] = 0; check_word(word);
  if (len >= 2 && word[len-2] == 'J' && word[len-1] == 'Z'){
    /* compare with the same prefix followed by O J */
    char w2[64]; memcpy(w2, word, len); w2[len-2] = 'O'; w2[len-1] = 'J'; w2[len] = 0;
    ll D1, E1, U1, d1, S1, D2, E2, U2, d2, S2; build(word); count(&D1, &E1, &U1, &d1, &S1); build(w2); count(&D2, &E2, &U2, &d2, &S2);
    njz++; if (D1 != D2 || E1 != E2 || U1 != U2 || d1 != d2) njz_bad++;
  }
  if (len == L) return;
  const char *A = "OZJ"; for (int i = 0; i < 3; i++){ word[len] = A[i]; rec_words(len + 1, L); }
}
static void check_tails(const char *w, int K){
  ll D, E, U, d, S; build(w); count(&D, &E, &U, &d, &S);
  ll h = D, R = U, e = E, dd = d;
  for (ll m = 0; m < (1LL << K); m++){
    char w2[96]; int L = (int)strlen(w); memcpy(w2, w, L);
    int nb = 0; while ((m >> nb)) nb++;
    for (int i = 0; i < nb; i++) w2[L + i] = ((m >> i) & 1) ? 'Z' : 'O';
    w2[L + nb] = 0;
    ll D2, E2, U2, d2, S2; build(w2); count(&D2, &E2, &U2, &d2, &S2);
    ll pred = h + R * m - dd * Qdef(m) - e * (m & 1);
    ntail++; if (D2 != pred){ ntail_bad++; if (ntail_bad < 10) printf("TAIL MISMATCH prefix %s m=%lld literal %lld formula %lld\n", w[0] ? w : "-", m, D2, pred); }
    if (S2 != 1) nunit_bad++;
  }
}
static void rec_prefix(int len, int L, int K){
  word[len] = 0; check_tails(word, K);
  if (len == L) return;
  const char *A = "OZJ"; for (int i = 0; i < 3; i++){ word[len] = A[i]; rec_prefix(len + 1, L, K); }
}
int main(int argc, char **argv){
  size_t NS = 1ULL << MAXN;
  f = malloc(NS * MAXN * sizeof(u64)); g = malloc(NS * MAXN * sizeof(u64)); Hs = malloc(NS * sizeof(u64));
  if (!f || !g || !Hs){ fprintf(stderr, "oom\n"); return 1; }
  if (argv[1][0] == 'W'){ int L = atoi(argv[2]); rec_words(0, L);
    printf("mode W L=%d: words checked %lld, response mismatches %lld, not unit-start %lld; JZ vs OJ pairs %lld, mismatches %lld\n", L, nchecked, nbad, nunit_bad, njz, njz_bad); }
  else if (argv[1][0] == 'T'){ int L = atoi(argv[2]), K = atoi(argv[3]); rec_prefix(0, L, K);
    printf("mode T L=%d K=%d: tail tournaments %lld, formula mismatches %lld, not unit-start %lld\n", L, K, ntail, ntail_bad, nunit_bad); }
  else if (argv[1][0] == 'S'){ FILE *fp = fopen(argv[2], "r"); int K = atoi(argv[3]); char buf[128];
    while (fgets(buf, sizeof buf, fp)){ buf[strcspn(buf, "\r\n ")] = 0; if (!buf[0]) continue; const char *w = strcmp(buf, "-") ? buf : "";
      ll D, E, U, d, S; build(w); count(&D, &E, &U, &d, &S);
      long long before = ntail_bad; check_tails(w, K);
      printf("seed word %-14s (h,R,e,d)=(%lld,%lld,%lld,%lld) n=%d tails m<2^%d: %s\n", buf, D, U, E, d, n, K, ntail_bad == before ? "all match" : "MISMATCH"); }
    printf("mode S: tail tournaments %lld, formula mismatches %lld, not unit-start %lld\n", ntail, ntail_bad, nunit_bad); }
  return 0;
}
