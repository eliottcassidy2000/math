/* cesaro_lazy.c - auditor's floating-point long-block limit of the letter-count model, by a different method from the
 * packet's exact_limit.py: the Cesaro projector of the O/Z block chain P is computed as lim_k ((I+P)/2)^k (the lazy
 * chain is aperiodic and has the same closed classes and stationary laws, so its powers converge to the Cesaro
 * projector of P even when P has periodic classes). Outputs the F-class laws of Pi J Pi delta_t (one J) and
 * Pi J Pi J Pi delta_t (two J), the closed forms for prime q, and the convergence gap between K/2 and K lazy steps.
 * usage: cesaro_lazy q K
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef long long ll;
static int q, S; static int *tO, *tZ, *tJ, *Fc;
static int md(ll x){ x %= q; if (x < 0) x += q; return (int)x; }
static int enc(int D, int E, int U, int d){ return ((D*q + E)*q + U)*q + d; }
static double *tmp;
static void Pi(double *v, int K){
  for (int k = 0; k < K; k++){
    for (int s = 0; s < S; s++) tmp[s] = 0.5*v[s];
    for (int s = 0; s < S; s++) if (v[s] != 0){ tmp[tO[s]] += 0.25*v[s]; tmp[tZ[s]] += 0.25*v[s]; }
    memcpy(v, tmp, S*sizeof(double));
  }
}
static void Jop(double *v){ memset(tmp, 0, S*sizeof(double)); for (int s = 0; s < S; s++) if (v[s] != 0) tmp[tJ[s]] += v[s]; memcpy(v, tmp, S*sizeof(double)); }
static void law(double *v, double *L){ double tot = 0; for (int r = 0; r < q; r++) L[r] = 0; for (int s = 0; s < S; s++){ L[Fc[s]] += v[s]; tot += v[s]; } for (int r = 0; r < q; r++) L[r] *= q/tot; }
static void limits(int K, double *L1, double *L2){
  double *v = calloc(S, sizeof(double)); v[enc(1 % q, 0, 2 % q, 1 % q)] = 1;
  Pi(v, K); Jop(v); Pi(v, K); law(v, L1); Jop(v); Pi(v, K); law(v, L2); free(v);
}
int main(int argc, char **argv){
  q = atoi(argv[1]); int K = atoi(argv[2]); S = q*q*q*q;
  tO = malloc(S*sizeof(int)); tZ = malloc(S*sizeof(int)); tJ = malloc(S*sizeof(int)); Fc = malloc(S*sizeof(int)); tmp = malloc(S*sizeof(double));
  for (int D = 0; D < q; D++) for (int E = 0; E < q; E++) for (int U = 0; U < q; U++) for (int d = 0; d < q; d++){
    int s = enc(D,E,U,d); tO[s] = enc(D, 0, md(2LL*U), d); tZ[s] = enc(md((ll)D-E+U), d, md(2LL*U-d), d);
    tJ[s] = enc(md((ll)D-E+U), md((ll)D-E), md((ll)D+U), D); Fc[s] = md((ll)D-E+U); }
  double *a1 = calloc(q, sizeof(double)), *a2 = calloc(q, sizeof(double)), *b1 = calloc(q, sizeof(double)), *b2 = calloc(q, sizeof(double));
  limits(K/2, a1, a2); limits(K, b1, b2);
  double gap = 0; for (int r = 0; r < q; r++){ gap = fmax(gap, fabs(a1[r]-b1[r])); gap = fmax(gap, fabs(a2[r]-b2[r])); }
  printf("q=%d lazy steps %d (gap to %d steps %.2e)\n  two-J L:", q, K, K/2, gap); for (int r = 0; r < q; r++) printf(" %.12f", b2[r]);
  printf("\n  one-J L:"); for (int r = 0; r < q; r++) printf(" %.12f", b1[r]); printf("\n");
  int prime = q > 2; for (int k = 2; k < q; k++) if (q % k == 0) prime = 0;
  if (prime){ double p = q;
    double c2 = 1 + (3*p-2)*(p-1)/(2*p*p*p), c2n = 1 - (3*p-2)/(2*p*p*p), c1 = 1 + (3*p-4)/(4*p*p), c1m = 1 - 1/(p*p) + 1/(4*p), c1n = 1 - 1/(p*p);
    double e2 = fabs(b2[0]-c2), e1 = fabs(b1[0]-c1) ; for (int r = 1; r < q; r++){ e2 = fmax(e2, fabs(b2[r]-c2n)); e1 = fmax(e1, fabs(b1[r] - (r == q-1 ? c1m : c1n))); }
    printf("  closed forms: two-J L(0) = %.12f (max |diff| over classes %.2e), one-J L(0) = %.12f (max |diff| %.2e); (L(0)-1)p = %.6f\n", c2, e2, c1, e1, (b2[0]-1)*p);
  }
  return 0;
}
