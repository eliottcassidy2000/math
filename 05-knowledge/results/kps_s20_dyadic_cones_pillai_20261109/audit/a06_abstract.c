/* a06_abstract.c (S20 audit): own abstract search for THM-5109 (D)/(E2)/(E3).
   States (h, e, R = 2^j): h odd, 1 <= h <= 2^j, 0 <= e <= h (optionally e = R mod 2).
   Cone step: transitive gadget TT_m, entry p >= 2, exit q (any q in 1..m), raw closed forms (Theorem B, verified
   against brute force in a01):
     q >= p      : R' = 2^(m-q) (h + 2^(p-1) R)
     2 <= q < p  : R' = 2^(m-q) h + 3^(p-q-1) 2^(m-p+q) R
     q = 1 < p   : R' = 2^(m-1) h - 2^(m-p) e + 3^(p-2) 2^(m-p+1) R
   Being a power of two does not depend on m (common factor 2^(m-p) when m >= p), so we test m = p; e' depends
   on m only through [m = p]:  e' = 2^(q-1) R if m = p (q < p), 0 if m > p.
   For every state with a dyadic step we check the (E2) classification and test the outputs (both m = p and
   m > p variants) for a second dyadic step with p' <= PMAX (any q', including q' >= p').
   q >= 2 steps do not depend on e and are evaluated once per (h, j).
   usage: a06_abstract JMAX PMAX PARITY(0/1) */
#include <stdio.h>
#include <stdlib.h>
typedef __int128 i128;
static i128 P3[200], P2[200];
static int ispow(i128 x){ return x > 0 && (x & (x - 1)) == 0; }
static int lg(i128 x){ int k = 0; while (x > 1) { x >>= 1; k++; } return k; }
static int PMAX;
static long long unexpected = 0;
/* does (h,e,2^j) admit a dyadic cone step (p >= 2, any exit)? */
static int has_step(i128 h, i128 e, int j){
  /* overflow guard: 2^(PMAX) h and 3^(PMAX-2) 2^(j+2) must stay below 2^126 */
  if (j > 30 || h > P2[64] || PMAX > 60) { printf("OVERFLOW RISK j=%d\n", j); exit(4); }
  i128 R = P2[j];
  for (int p = 2; p <= PMAX; p++) {
    for (int q = 1; q <= p + 1; q++) {   /* q = p, p+1 (with m = p+1) test the 'exit at/after entry' rows */
      int m = (q <= p) ? p : p + 1;
      i128 Rp;
      if (q >= p) Rp = P2[m - q] * (h + P2[p - 1] * R);
      else if (q >= 2) Rp = P2[m - q] * h + P3[p - q - 1] * P2[m - p + q] * R;
      else Rp = P2[m - 1] * h - P2[m - p] * e + P3[p - 2] * P2[m - p + 1] * R;
      if (ispow(Rp)) { if (q >= p) unexpected++; return 1; }
    }
  }
  return 0;
}
static int classified(long long h, long long e, int j){
  if (h == 1 || h == 5 || h == 13) return 1;
  long long dp[3] = {1, 5, 13}; int ps[3] = {3, 5, 7};
  for (int k = 0; k < 3; k++) if (e == (1LL << (ps[k] - 1)) * h - dp[k] * (1LL << (j + 1))) return 1;
  return 0;
}
static long long chains = 0;
static void test_outputs(long long h, long long e, int j, int p, int q, i128 Rp, i128 R){
  i128 hp = (i128)h + P2[p - 1] * R; int jp = lg(Rp);
  i128 oe[2] = { P2[q - 1] * R, 0 }; int oj[2] = { jp, jp + 1 };
  for (int k = 0; k < 2; k++) {
    if (has_step(hp, oe[k], oj[k])) {
      chains++;
      if (chains <= 10) printf("CHAIN: (%lld,%lld,2^%d) -(p=%d,q=%d,%s)-> (%lld,%lld,2^%d) -> dyadic\n", h, e, j, p, q,
                               k ? "m>p" : "m=p", (long long)hp, (long long)oe[k], oj[k]);
    }
  }
}
int main(int argc, char **argv){
  int JMAX = atoi(argv[1]); PMAX = atoi(argv[2]); int PAR = atoi(argv[3]);
  P3[0] = 1; P2[0] = 1; for (int i = 1; i < 200; i++) { P3[i] = P3[i - 1] * 3; P2[i] = P2[i - 1] * 2; }
  long long states = 0, withstep = 0, viol = 0, a1 = 0, a513 = 0, other = 0, qge2states = 0;
  for (int j = 0; j <= JMAX; j++) {
    long long Rl = 1LL << j; i128 R = Rl;
    for (long long h = 1; h <= Rl; h += 2) {
      /* q >= 2 dyadic steps of (h, *, 2^j): list them once */
      int nq = 0; int qp[64], qq[64]; i128 qR[64];
      for (int p = 3; p <= PMAX; p++) for (int q = 2; q < p; q++) {
        i128 Rp = P2[p - q] * (i128)h + P3[p - q - 1] * P2[q] * R;
        if (ispow(Rp) && nq < 64) { qp[nq] = p; qq[nq] = q; qR[nq] = Rp; nq++; }
      }
      for (long long e = 0; e <= h; e++) {
        if (PAR && ((e - Rl) & 1)) continue;
        states++;
        int found = 0;
        for (int k = 0; k < nq; k++) { found = 1; test_outputs(h, e, j, qp[k], qq[k], qR[k], R); }
        if (nq) qge2states++;
        for (int p = 2; p <= PMAX; p++) {
          i128 Rp = P2[p - 1] * (i128)h - (i128)e + P3[p - 2] * P2[1] * R;     /* q = 1, m = p */
          if (!ispow(Rp)) continue;
          found = 1; test_outputs(h, e, j, p, 1, Rp, R);
        }
        if (found) {
          withstep++;
          if (!classified(h, e, j)) { viol++; if (viol <= 10) printf("UNCLASSIFIED step state (%lld,%lld,2^%d)\n", h, e, j); }
          if (h == 1) a1++; else if (h == 5 || h == 13) a513++; else other++;
        }
      }
    }
  }
  printf("JMAX=%d PMAX=%d parity=%d: states %lld, with a dyadic step %lld (h=1: %lld, h in {5,13}: %lld, other h: %lld; "
         "states with a q>=2 step: %lld), classification violations %lld, two-step chains %lld, dyadic exits q>=p: %lld\n",
         JMAX, PMAX, PAR, states, withstep, a1, a513, other, qge2states, viol, chains, unexpected);
  return 0;
}
