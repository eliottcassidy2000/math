/* h_parity.c - auditor's numerical check of the parity-split capacity (THM-4888 (C) second half; Session 206
 * weighted-parameter-capacity section 1):  #{m even : H_s(m) <= X} ~ (X/d) ln(R/(R-d/2)),
 * #{m odd : H_s(m) <= X} ~ (X/d) ln((R-d/2)/(R-d)).  Full linear scan of m, Q array by the recurrence.
 * usage: h_parity X  (seeds hard-coded below)  prints counts, predictions and (count - prediction)/sqrt(X). */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
typedef long long ll;
int main(int argc, char **argv){
  ll X = atoll(argv[1]);
  ll S[][4] = {{1,2,0,1},{3,3,1,1},{5,6,2,3},{9,11,3,5},{11,12,4,5},{15,19,5,9},{81,99,27,45},{51,67,17,33},{159,160,78,79},{25,30,10,15}};
  int ns = sizeof S / sizeof S[0];
  ll mmax = X + 2; uint32_t *Q = malloc((size_t)(mmax + 1) * 4); Q[0] = 0;
  for (ll m = 1; m <= mmax; m++){ ll k = m >> 1; Q[m] = (m & 1) ? (uint32_t)(Q[k] + k + (k & 1)) : Q[k]; }
  for (int s = 0; s < ns; s++){
    ll h = S[s][0], R = S[s][1], e = S[s][2], d = S[s][3]; ll ce = 0, co = 0;
    ll lim = (X - h) / (R - d);
    for (ll m = 0; m <= lim; m++){ ll H = h + R * m - d * (ll)Q[m] - e * (m & 1); if (H <= X){ if (m & 1) co++; else ce++; } }
    double pe = (double)X / d * log((double)R / (R - d / 2.0)), po = (double)X / d * log((R - d / 2.0) / (R - d));
    printf("seed (%lld,%lld,%lld,%lld) X=%lld even %lld pred %.1f err/sqrtX %+.3f | odd %lld pred %.1f err/sqrtX %+.3f | ratio even/odd obs %.5f pred %.5f\n",
      h, R, e, d, X, ce, pe, (ce - pe) / sqrt((double)X), co, po, (co - po) / sqrt((double)X), (double)ce / co, pe / po);
  }
  return 0;
}
