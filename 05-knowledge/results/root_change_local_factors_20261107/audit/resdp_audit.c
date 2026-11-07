/* resdp_audit.c - independent auditor's residue automaton for the letter-count model of HYP-9506 (A).
 *
 * Model (stated in the packet, re-implemented here from the definition): count every prefix w with exactly n letters
 * O/Z and at most two letters J (J anywhere: first, last, repeated), with weight 2 if #J(w) <= 1 (it is counted as wZ
 * and wJ in M_2) and weight 1 if #J(w) = 2 (only wZ), and record F(t.w) = D - E + U mod q. L_n(r) = q * share(r).
 *
 * Implementation: three level vectors A_j (exactly j letters J, any last letter) over the q^4 residue states with
 * the recurrence A_0(n) = P A_0(n-1), A_1(n) = J A_0(n) + P A_1(n-1), A_2(n) = J A_1(n) + P A_2(n-1), where P = O + Z.
 * In mode "exact" the counts are exact 64-bit integers (no normalisation, small n only); otherwise long double with
 * renormalisation (shares are scale-free).
 *
 * usage: resdp_audit q nmax list|all [exact]
 *   prints, for each n in the comma list (or every n <= nmax with "all"), the model shares by class and the level-1
 *   and level-2 shares separately.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
typedef unsigned long long ull;
typedef long double LD;

static int q; static int S;
static int *tO, *tZ, *tJ, *Fc;

static int md(ll x){ x %= q; if (x < 0) x += q; return (int)x; }
static int enc(int D, int E, int U, int d){ return ((D*q + E)*q + U)*q + d; }

int main(int argc, char **argv){
  if (argc < 4){ fprintf(stderr, "usage: resdp_audit q nmax list|all [exact]\n"); return 1; }
  q = atoi(argv[1]); int nmax = atoi(argv[2]); int all = !strcmp(argv[3], "all"); int exact = argc > 4 && !strcmp(argv[4], "exact");
  char *want = calloc(nmax + 1, 1);
  if (!all){ char *s = strdup(argv[3]); for (char *tok = strtok(s, ","); tok; tok = strtok(NULL, ",")){ int n = atoi(tok); if (n >= 0 && n <= nmax) want[n] = 1; } }
  S = q*q*q*q;
  tO = malloc(S*sizeof(int)); tZ = malloc(S*sizeof(int)); tJ = malloc(S*sizeof(int)); Fc = malloc(S*sizeof(int));
  for (int D = 0; D < q; D++) for (int E = 0; E < q; E++) for (int U = 0; U < q; U++) for (int d = 0; d < q; d++){
    int s = enc(D,E,U,d);
    tO[s] = enc(D, 0, md(2LL*U), d);
    tZ[s] = enc(md((ll)D - E + U), d, md(2LL*U - d), d);
    tJ[s] = enc(md((ll)D - E + U), md((ll)D - E), md((ll)D + U), D);
    Fc[s] = md((ll)D - E + U);
  }
  int t0 = enc(1 % q, 0, 2 % q, 1 % q);
  if (exact){
    ull *A[3], *B[3]; for (int j = 0; j < 3; j++){ A[j] = calloc(S, sizeof(ull)); B[j] = calloc(S, sizeof(ull)); }
    A[0][t0] = 1;
    for (int n = 0; n <= nmax; n++){
      for (int s = 0; s < S; s++) if (A[0][s]) A[1][tJ[s]] += A[0][s];
      for (int s = 0; s < S; s++) if (A[1][s]) A[2][tJ[s]] += A[1][s];    /* includes the J just added: JJ allowed */
      if (all || want[n]){
        ull *cls = calloc(q, sizeof(ull)), tot = 0;
        for (int j = 0; j < 3; j++){ ull w = j < 2 ? 2 : 1; for (int s = 0; s < S; s++) if (A[j][s]){ cls[Fc[s]] += w*A[j][s]; tot += w*A[j][s]; } }
        printf("q=%d n=%d exact weighted counts (total %llu):", q, n, tot); for (int r = 0; r < q; r++) printf(" %llu", cls[r]); printf("\n"); free(cls);
      }
      if (n == nmax) break;
      for (int j = 0; j < 3; j++){ memset(B[j], 0, S*sizeof(ull)); for (int s = 0; s < S; s++) if (A[j][s]){ B[j][tO[s]] += A[j][s]; B[j][tZ[s]] += A[j][s]; } ull *x = A[j]; A[j] = B[j]; B[j] = x; }
    }
    return 0;
  }
  LD *A[3], *B[3]; for (int j = 0; j < 3; j++){ A[j] = calloc(S, sizeof(LD)); B[j] = calloc(S, sizeof(LD)); }
  A[0][t0] = 1;
  /* level weights relative to the stored (renormalised) vectors: all three share one scale factor */
  for (int n = 0; n <= nmax; n++){
    for (int s = 0; s < S; s++) if (A[0][s] != 0) A[1][tJ[s]] += A[0][s];
    for (int s = 0; s < S; s++) if (A[1][s] != 0) A[2][tJ[s]] += A[1][s];
    if (all || want[n]){
      LD *cls = calloc(q, sizeof(LD)), *c1 = calloc(q, sizeof(LD)), *c2 = calloc(q, sizeof(LD)); LD tot = 0, t1 = 0, t2 = 0;
      for (int j = 0; j < 3; j++){ LD w = j < 2 ? 2 : 1; for (int s = 0; s < S; s++) if (A[j][s] != 0){ cls[Fc[s]] += w*A[j][s]; tot += w*A[j][s];
        if (j == 1){ c1[Fc[s]] += A[j][s]; t1 += A[j][s]; } if (j == 2){ c2[Fc[s]] += A[j][s]; t2 += A[j][s]; } } }
      printf("q=%d n=%d model L:", q, n); for (int r = 0; r < q; r++) printf(" %.10Lf", q*cls[r]/tot);
      printf(" | level1 L:"); for (int r = 0; r < q; r++) printf(" %.10Lf", q*c1[r]/t1);
      printf(" | level2 L:"); for (int r = 0; r < q; r++) printf(" %.10Lf", q*c2[r]/t2);
      printf("\n"); fflush(stdout); free(cls); free(c1); free(c2);
    }
    if (n == nmax) break;
    LD tot = 0;
    for (int j = 0; j < 3; j++){ memset(B[j], 0, S*sizeof(LD)); for (int s = 0; s < S; s++) if (A[j][s] != 0){ B[j][tO[s]] += A[j][s]; B[j][tZ[s]] += A[j][s]; } }
    for (int j = 0; j < 3; j++) for (int s = 0; s < S; s++) tot += B[j][s];
    for (int j = 0; j < 3; j++){ for (int s = 0; s < S; s++) B[j][s] /= tot; LD *x = A[j]; A[j] = B[j]; B[j] = x; }
  }
  return 0;
}
