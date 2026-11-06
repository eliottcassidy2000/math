/* greedyF2.c (opus-2026-11-06-S82, exploratory): greedy minimal sub-bank Y of the order-<=7 rooted triples such that
 *   F2(Y) = { y_h g(v) - y_e e(v) + y_R (1+v) : y in Y, v >= 0 }
 * covers every odd N in [1, 2^LOGX] except 7 and 21 (F2 alone, no F1).  Prints the picks in order with their gains.
 * Usage: greedyF2 LOGX bankfile [maxpicks]
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>
typedef long long ll;
static uint32_t *gtab; static ll X; static uint8_t *cov;
int main(int argc, char **argv) {
    int LOGX = atoi(argv[1]); X = 1LL << LOGX; int maxp = argc > 3 ? atoi(argv[3]) : 400;
    FILE *f = fopen(argv[2], "r"); int nb = 0; static ll yh[4096], ye[4096], yR[4096];
    while (fscanf(f, "%lld %lld %lld", &yh[nb], &ye[nb], &yR[nb]) == 3) nb++;
    fclose(f);
    gtab = malloc(sizeof(uint32_t) * (size_t)(X + 2)); gtab[0] = 1;
    for (ll v = 1; v <= X + 1; v++) { if (v & 1) gtab[v] = gtab[v >> 1]; else { ll u = v >> 1; gtab[v] = (uint32_t)((ll)gtab[u] + u + (u & 1)); } }
    cov = calloc((size_t)(X / 2 + 1), 1);
    cov[7 >> 1] = 1; cov[21 >> 1] = 1;                         /* never attainable: excluded from the universe */
    ll uncovered = X / 2 - 2; int *used = calloc(nb, sizeof(int)); ll *gain = malloc(sizeof(ll) * nb);
    printf("# greedyF2 LOGX=%d bank=%d: universe = odd N in [1, 2^%d] minus {7, 21}: %lld\n", LOGX, nb, LOGX, uncovered);
    for (int pick = 0; pick < maxp && uncovered > 0; pick++) {
        #pragma omp parallel for schedule(dynamic, 4)
        for (int i = 0; i < nb; i++) {
            gain[i] = 0; if (used[i]) continue;
            ll cnt = 0;
            for (ll v = 0; (1 + v) * yR[i] <= X; v++) {
                ll ev = (v == 0 || !(v & 1)); ll N = (ll)gtab[v] * yh[i] - ev * ye[i] + (1 + v) * yR[i];
                if (N >= 1 && N <= X && (N & 1) && !cov[N >> 1]) cnt++;   /* counts multiplicity within one family too */
            }
            gain[i] = cnt;
        }
        int best = -1; for (int i = 0; i < nb; i++) if (!used[i] && (best < 0 || gain[i] > gain[best])) best = i;
        if (best < 0 || gain[best] == 0) break;
        used[best] = 1; ll newc = 0;
        for (ll v = 0; (1 + v) * yR[best] <= X; v++) {
            ll ev = (v == 0 || !(v & 1)); ll N = (ll)gtab[v] * yh[best] - ev * ye[best] + (1 + v) * yR[best];
            if (N >= 1 && N <= X && (N & 1) && !cov[N >> 1]) { cov[N >> 1] = 1; newc++; }
        }
        uncovered -= newc;
        printf("pick %3d: y=(%lld,%lld,%lld) new %lld remaining %lld\n", pick + 1, yh[best], ye[best], yR[best], newc, uncovered);
        fflush(stdout);
    }
    if (uncovered > 0) { printf("UNCOVERED (first 40):"); ll s = 0; for (ll N = 1; N <= X && s < 40; N += 2) if (!cov[N >> 1]) { printf(" %lld", N); s++; } printf("\n"); }
    return 0;
}
