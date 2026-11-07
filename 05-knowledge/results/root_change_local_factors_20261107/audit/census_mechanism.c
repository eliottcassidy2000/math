/* census_mechanism.c - auditor's attribution of the zero-class excess of M_2 at census scale (window [LO,HI)).
 * Same enumeration as census_lf_audit.c (letter-by-letter DFS over prefixes w with <= 2 J, counted words wZ and, if
 * #J(w) <= 1, wJ, both of value F(t.w); F strictly increases along letters, so pruning at F > X is exact).
 * Every counted word is put in one category, for each prime p in {3,5,7,11,13,17,23,31} (first match wins):
 *   PREF : the word begins with the root rectangle Z^a O^b J (a, b >= 1) and p | 2^a + 1. Then t.Z^a O^b J =
 *          (h, g, h, g) with g = 2^a + 1, h = g(2^b + 1): every coordinate is 0 mod g, so EVERY continuation has
 *          value divisible by 2^a + 1 (asserted below for every such prefix);
 *   MULT : the prefix ends with a multiplier block after its last J (THM-4916 (1c) and its OJ twin), so F(t.w) is
 *          (2^k+1) times an integer: tail after the last J is Z O^k (k >= 1)  [value (2^k+1)(D+2U) of the pre-J state],
 *          or O^k (k >= 1) with E = 0 in the pre-J state                        [value (2^k+1)(D+U)],
 *          or, with no J, the root prefix Z O^k (k >= 1)                       [value 3(2^k+1), t.Z = t.J];
 *   SEED : not MULT, and the last root change is seed-confined mod p: p | gcd(D,U) of the pre-J state (THM-4916 (1a));
 *   REST : everything else.
 * For each p and category the program prints the class-0 mean, the mean over nonzero classes, and the share of the total
 * zero-class excess (mean(N = 0 mod p) - mean(N != 0 mod p)) carried by the category; rho_p with and without MULT.
 * usage: census_mechanism X LO HI nthreads
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
typedef long long ll; typedef unsigned long long ull;
#define NP 8
static const int PR[NP] = {3, 5, 7, 11, 13, 17, 23, 31};
static ll X, LO, HI;
typedef struct { ull c[NP][4][32]; ull words, prefviol, prefchk; } acc_t;
static acc_t *ACC;
enum { EMPTY = 0, ZOS = 1, OS = 2, OTHER = 3 };
enum { R_EMPTY = 0, R_Z = 1, R_ZO = 2, R_OTHER = 3 };
typedef struct { ll D, E, U, d; int j; int conf; int shape, k, preE0, isroot; int rs, ra, pa; } node;

static inline void record(acc_t *A, ll N, const node *w, int asJ){
  A->words++;
  if (w->pa > 0){ A->prefchk++; if (N % ((1LL << w->pa) + 1)) A->prefviol++; }
  if (N < LO || N >= HI) return;
  int mult, conf;
  if (asJ){ /* word wJ: value F(t.w); its category is decided by the prefix w */
    mult = 0; conf = w->conf;
  } else { mult = 0; conf = w->conf; }
  if (w->shape == ZOS && w->k >= 1) mult = 1;
  else if (w->shape == OS && w->k >= 1 && !w->isroot && w->preE0) mult = 1;
  for (int i = 0; i < NP; i++){
    int pref = w->pa > 0 && ((1LL << w->pa) + 1) % PR[i] == 0;
    int cat = pref ? 3 : mult ? 0 : ((conf >> i) & 1) ? 1 : 2;
    A->c[i][cat][N % PR[i]]++;
  }
}
static void visit(acc_t *A, const node *w){
  ll D = w->D, E = w->E, U = w->U, d = w->d, N = D - E + U;
  record(A, N, w, 0);
  if (w->j <= 1) record(A, N, w, 1);
  node c;
  if (D + 2*U <= X){ c = *w; c.E = 0; c.U = 2*U;
    c.shape = w->shape == EMPTY ? OS : w->shape == OS ? OS : w->shape == ZOS ? ZOS : OTHER; c.k = (w->shape == EMPTY) ? 1 : w->k + 1;
    if (w->j == 0){ c.rs = (w->rs == R_Z || w->rs == R_ZO) ? R_ZO : R_OTHER; }
    visit(A, &c); }
  if (D - E + 3*U - 2*d <= X){ c = *w; c.D = D - E + U; c.E = d; c.U = 2*U - d;
    c.shape = w->shape == EMPTY ? ZOS : OTHER; c.k = 0;
    if (w->j == 0){ c.rs = w->rs == R_EMPTY ? R_Z : w->rs == R_Z ? R_Z : R_OTHER; c.ra = w->rs == R_Z ? w->ra + 1 : 1; }
    visit(A, &c); }
  if (w->j < 2 && D + 2*U <= X){
    c.D = D - E + U; c.E = D - E; c.U = D + U; c.d = D; c.j = w->j + 1; c.isroot = 0;
    c.conf = 0; for (int i = 0; i < NP; i++) if (D % PR[i] == 0 && U % PR[i] == 0) c.conf |= 1 << i;
    c.shape = EMPTY; c.k = 0; c.preE0 = (E == 0);
    c.rs = R_OTHER; c.ra = w->ra; c.pa = (w->j == 0 && w->rs == R_ZO) ? w->ra : w->pa;
    visit(A, &c); }
}
static node *T; static size_t nT = 0, capT = 0;
static void collect(acc_t *A, const node *w, int depth){
  if (depth == 12){ if (nT == capT){ capT = capT ? 2*capT : 4096; T = realloc(T, capT*sizeof(node)); } T[nT++] = *w; return; }
  ll D = w->D, E = w->E, U = w->U, d = w->d, N = D - E + U;
  record(A, N, w, 0); if (w->j <= 1) record(A, N, w, 1);
  node c;
  if (D + 2*U <= X){ c = *w; c.E = 0; c.U = 2*U; c.shape = w->shape == EMPTY ? OS : w->shape == OS ? OS : w->shape == ZOS ? ZOS : OTHER; c.k = (w->shape == EMPTY) ? 1 : w->k + 1; if (w->j == 0){ c.rs = (w->rs == R_Z || w->rs == R_ZO) ? R_ZO : R_OTHER; } collect(A, &c, depth + 1); }
  if (D - E + 3*U - 2*d <= X){ c = *w; c.D = D - E + U; c.E = d; c.U = 2*U - d; c.shape = w->shape == EMPTY ? ZOS : OTHER; c.k = 0; if (w->j == 0){ c.rs = w->rs == R_EMPTY ? R_Z : w->rs == R_Z ? R_Z : R_OTHER; c.ra = w->rs == R_Z ? w->ra + 1 : 1; } collect(A, &c, depth + 1); }
  if (w->j < 2 && D + 2*U <= X){ c.D = D - E + U; c.E = D - E; c.U = D + U; c.d = D; c.j = w->j + 1; c.isroot = 0;
    c.conf = 0; for (int i = 0; i < NP; i++) if (D % PR[i] == 0 && U % PR[i] == 0) c.conf |= 1 << i;
    c.shape = EMPTY; c.k = 0; c.preE0 = (E == 0); c.rs = R_OTHER; c.ra = w->ra; c.pa = (w->j == 0 && w->rs == R_ZO) ? w->ra : w->pa; collect(A, &c, depth + 1); }
}
int main(int argc, char **argv){
  if (argc < 5){ fprintf(stderr, "usage: census_mechanism X LO HI nthreads\n"); return 1; }
  X = atoll(argv[1]); LO = atoll(argv[2]); HI = atoll(argv[3]); int nth = atoi(argv[4]); omp_set_num_threads(nth);
  ACC = calloc(nth + 1, sizeof(acc_t));
  node root = { 1, 0, 2, 1, 0, 0, EMPTY, 0, 1, 1, R_EMPTY, 0, 0 };
  collect(&ACC[nth], &root, 0);
  #pragma omp parallel for schedule(dynamic, 64)
  for (ll i = 0; i < (ll)nT; i++) visit(&ACC[omp_get_thread_num()], &T[i]);
  acc_t *S = &ACC[nth];
  for (int t = 0; t < nth; t++){ S->words += ACC[t].words; S->prefviol += ACC[t].prefviol; S->prefchk += ACC[t].prefchk; for (int i = 0; i < NP; i++) for (int a = 0; a < 4; a++) for (int r = 0; r < 32; r++) S->c[i][a][r] += ACC[t].c[i][a][r]; }
  printf("census_mechanism X=%lld window [%lld,%lld) counted words %llu; counted words (all values <= X) with a root-rectangle prefix Z^aO^bJ: %llu, of which NOT divisible by 2^a+1: %llu\n", X, LO, HI, S->words, S->prefchk, S->prefviol);
  const char *nm[4] = {"MULT", "SEED", "REST", "PREF"};
  for (int i = 0; i < NP; i++){
    int p = PR[i]; ll n0 = 0, nz = 0; for (ll v = LO | 1; v < HI; v += 2){ if (v % p) nz++; else n0++; }
    double m0[4], mz[4], tot0 = 0, totz = 0;
    for (int a = 0; a < 4; a++){ ull s0 = S->c[i][a][0], sz = 0; for (int r = 1; r < p; r++) sz += S->c[i][a][r];
      m0[a] = (double)s0 / n0; mz[a] = (double)sz / nz; tot0 += m0[a]; totz += mz[a]; }
    double ex = tot0 - totz;
    printf("p=%2d rho=%.4f | without MULT %.4f | without PREF %.4f | without PREF,MULT %.4f | REST only %.4f | excess %.3f:", p, tot0/totz,
           (tot0 - m0[0])/(totz - mz[0]), (tot0 - m0[3])/(totz - mz[3]), (tot0 - m0[0] - m0[3])/(totz - mz[0] - mz[3]), m0[2]/mz[2], ex);
    for (int a = 3; a >= 0; a--) printf(" %s mean0 %.3f meanNZ %.3f share %.1f%%;", nm[a], m0[a], mz[a], 100*(m0[a]-mz[a])/ex);
    printf("\n");
  }
  return 0;
}
