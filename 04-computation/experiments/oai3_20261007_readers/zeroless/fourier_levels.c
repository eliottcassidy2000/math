/* fourier_levels.c -- level sums of the Fourier transform of the zeroless-digit measure mu on Z_2.
   For m >= 1 and t odd mod 2^m:  P_m(zeta^t) = prod_{i<m} F(10^i t / 2^m),  F(th) = sum_{a=1}^9 e(a th).
   Reindexing u = 5^(m-1) t:  P_m = prod_{l=0}^{m-1} F(omega_l(u)),  omega_l(u) = (5^-l u mod 2^(l+1)) / 2^(l+1),
   which depends only on u mod 2^(l+1): one binary DFS gives all levels.
   Outputs  G_m = sum_t P_m(zeta^t)  (an integer: the trace of the cyclotomic unit P_m(zeta_{2^m}); G_m = 2^(m-1) Delta_{m-1})
   and      H_m = sum_t |P_m(zeta^t)|   (L1 Fourier mass at level m; sum_m H_m/9^m < inf implies a continuous density). */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>
static int M;
static double complex Gs[64]; static double Hs[64], Mx[64];
static uint64_t inv5pow[64];
static inline double complex Fval(double th) {
  double s1 = sin(M_PI * th);
  double r = sin(9.0 * M_PI * th) / s1;
  return r * cexp(2.0 * M_PI * I * 5.0 * th);
}
static void dfs(int l, uint64_t u, double complex prod) {
  /* u known mod 2^(l+1) */
  uint64_t mod = (l + 1 >= 64) ? ~0ULL : ((1ULL << (l + 1)) - 1);
  uint64_t num = (inv5pow[l] * u) & mod;
  double th = (double)num / (double)(1ULL << (l + 1));
  double complex p = prod * Fval(th);
  Gs[l + 1] += p; double a = cabs(p); Hs[l + 1] += a; if (a > Mx[l + 1]) Mx[l + 1] = a;
  if (l + 1 >= M) return;
  dfs(l + 1, u, p);
  dfs(l + 1, u | (1ULL << (l + 1)), p);
}
int main(int argc, char **argv) {
  M = atoi(argv[1]);
  uint64_t inv = 1, five = 5; for (int i = 0; i < 7; i++) inv *= 2 - five * inv;  /* 5^-1 mod 2^64 */
  inv5pow[0] = 1; for (int l = 1; l < 64; l++) inv5pow[l] = inv5pow[l-1] * inv;
  dfs(0, 1, 1.0);
  for (int m = 1; m <= M; m++)
    printf("m=%d G_m=%.6f (imag %.2e) H_m=%.6e maxterm=%.6e  H_m/9^m=%.6e  H_m^(1/m)=%.5f  (H_m/2^(m-1))^(1/m)=%.5f\n", m, creal(Gs[m]), cimag(Gs[m]), Hs[m], Mx[m],
           Hs[m] / pow(9.0, m), pow(Hs[m], 1.0 / m), pow(Hs[m] / pow(2.0, m - 1), 1.0 / m));
  return 0;
}
