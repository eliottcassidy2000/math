/* S19 audit: second, structurally different coverage computation (no state closure).
 * Every fixed-port tail value over a core (h0,R0,D0,hu) is h(m) = h0 + (x/2)[R0 m + D0 [m odd] - hu Q(m)], m >= 0,
 * (THM-5085 (2), verified by check_transfer.py), with Q(0)=0, Q(2m)=Q(m), Q(2m+1)=Q(m)+m+(m mod 2).
 * Enumerate m directly, with an a-priori bound from 0 <= Q(m) <= m-1 (m >= 1):
 *  - non-sink port (R0 > hu): h(m) >= h0 + (x/2)[(R0-hu)m + hu - |D0|], so m <= (LIM - h0 + (x/2)|D0|)/((x/2)(R0-hu)) + 1;
 *  - sink port (R0 == hu, D0 == -hu): h(2k+1) = h(k), and h(2k) >= h0 + (x/2) hu (k+1), so only m = 0 and
 *    m = 2k with k <= (LIM - h0)/((x/2) hu) are needed.
 * Usage: cover2 CORES X LIM [HITS_FROM_cover TO_COMPARE]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r"); ll X = atoll(argv[2]), LIM = atoll(argv[3]);
  unsigned char *hit = calloc(LIM + 1, 1);
  ll QMAX = 2 * LIM + 100; int *Q = malloc((QMAX + 1) * sizeof(int)); Q[0] = 0; /* Q(m) <= m - 1 fits in 32 bits */
  for (ll m = 1; m <= QMAX; m++) { ll k = m >> 1; Q[m] = (int)((m & 1) ? Q[k] + k + (k & 1) : Q[k]); }
  ll h0, R0, D0, hu, nc = 0, evals = 0;
  while (fscanf(f, "%lld %lld %lld %lld", &h0, &R0, &D0, &hu) == 4) {
    nc++;
    if (h0 <= LIM) hit[h0] = 1;
    if (R0 > hu) {
      ll den = X * (R0 - hu); /* (x/2)(R0-hu) * 2 */
      ll mb = (2 * (LIM - h0) + X * llabs(D0)) / den + 2;
      if (mb > QMAX) { fprintf(stderr, "QMAX too small\n"); return 2; }
      for (ll m = 1; m <= mb; m++) {
        ll br = R0 * m + ((m & 1) ? D0 : 0) - hu * (ll)Q[m];
        if (br % 2) { fprintf(stderr, "odd bracket\n"); return 3; }
        ll h = h0 + X * (br / 2); evals++;
        if (h >= 0 && h <= LIM) hit[h] = 1;
      }
    } else if (R0 == hu && D0 == -hu) {
      ll kb = (2 * (LIM - h0)) / (X * hu) + 2;
      for (ll k = 1; k <= kb; k++) {
        ll m = 2 * k; if (m > QMAX) { fprintf(stderr, "QMAX too small (sink)\n"); return 2; }
        ll br = R0 * m - hu * (ll)Q[m];
        ll h = h0 + X * (br / 2); evals++;
        if (h >= 0 && h <= LIM) hit[h] = 1;
      }
    } else { fprintf(stderr, "unexpected core with R0 <= hu: %lld %lld %lld %lld\n", h0, R0, D0, hu); return 4; }
  }
  ll nmiss = 0;
  printf("closed-form enumeration: x=%lld LIM=%lld cores=%lld evaluations=%lld misses:", X, LIM, nc, evals);
  for (ll v = 1; v <= LIM; v++) { if (v % X != 1 % X || (X == 2 && v % 2 == 0)) continue; if (!hit[v]) { nmiss++; if (nmiss <= 400) printf(" %lld", v); } }
  printf("\ncount %lld\n", nmiss);
  if (argc > 4) {
    FILE *g = fopen(argv[4], "rb"); unsigned char *h2 = malloc(LIM + 1);
    if (fread(h2, 1, LIM + 1, g) != (size_t)(LIM + 1)) { fprintf(stderr, "short read\n"); return 5; }
    ll diff = 0; for (ll v = 0; v <= LIM; v++) if (!!h2[v] != !!hit[v]) diff++;
    printf("hit sets compared with the state-closure run: %lld differing values in [0, %lld]\n", diff, LIM);
  }
  return 0;
}
