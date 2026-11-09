/* S19 audit, Q4: core covers of the targets in [LO, LIM] (odd for x = 2, = 1 mod 3 for x = 3).
 * Usage: q4 CORES X LO LIM S19_GREEDY_OUT SPARSE_OUT KMAX
 *  1. per-core value sets (closed-form enumeration, as in cover2.c) restricted to the targets;
 *  2. checks that the core states listed in S19's greedy output are all in CORES and that together they cover every target;
 *  3. an independent greedy cover (ties broken by the order of CORES);
 *  4. a lower bound: a set of targets whose covering-core sets are pairwise disjoint (any cover needs a distinct core
 *     for each of them), built greedily from the targets with the fewest covering cores, then re-verified;
 *  5. writes, for every target covered by at most KMAX cores, the list of its covering cores (for an LP / MILP bound).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef long long ll;
static ll X, LO, LIM, NT;
static ll tid(ll v) { if (v < LO || v > LIM) return -1; if (X == 2) { if (v % 2 == 0) return -1; return (v - LO) / 2; } if (v % X != 1) return -1; return (v - LO) / X; }
static ll tval(ll t) { return LO + (X == 2 ? 2 : X) * t; }
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r"); X = atoll(argv[2]); LO = atoll(argv[3]); LIM = atoll(argv[4]);
  const char *s19 = argv[5], *sparse = argv[6]; int KMAX = atoi(argv[7]);
  if (X == 2 && LO % 2 == 0) LO++;
  while (X != 2 && LO % X != 1) LO++;
  NT = (LIM - LO) / (X == 2 ? 2 : X) + 1;
  ll W = (NT + 63) / 64;
  ll (*core)[4] = malloc(sizeof(*core) * 10000); int nc = 0;
  while (fscanf(f, "%lld %lld %lld %lld", &core[nc][0], &core[nc][1], &core[nc][2], &core[nc][3]) == 4) nc++;
  uint64_t **set = malloc(sizeof(uint64_t *) * nc);
  ll QMAX = 2 * LIM + 100; int *Q = malloc((QMAX + 1) * sizeof(int)); Q[0] = 0;
  for (ll m = 1; m <= QMAX; m++) { ll k = m >> 1; Q[m] = (int)((m & 1) ? Q[k] + k + (k & 1) : Q[k]); }
  for (int c = 0; c < nc; c++) {
    set[c] = calloc(W, 8);
    ll h0 = core[c][0], R0 = core[c][1], D0 = core[c][2], hu = core[c][3];
    ll t = tid(h0); if (t >= 0) set[c][t >> 6] |= 1ULL << (t & 63);
    if (h0 > LIM) continue;
    if (R0 > hu) {
      ll mb = (2 * (LIM - h0) + X * llabs(D0)) / (X * (R0 - hu)) + 2;
      for (ll m = 1; m <= mb; m++) { ll br = R0 * m + ((m & 1) ? D0 : 0) - hu * (ll)Q[m]; ll h = h0 + X * (br / 2); t = tid(h); if (t >= 0) set[c][t >> 6] |= 1ULL << (t & 63); }
    } else {
      ll kb = (2 * (LIM - h0)) / (X * hu) + 2;
      for (ll k = 1; k <= kb; k++) { ll m = 2 * k; ll br = R0 * m - hu * (ll)Q[m]; ll h = h0 + X * (br / 2); t = tid(h); if (t >= 0) set[c][t >> 6] |= 1ULL << (t & 63); }
    }
  }
  /* union and per-target cover counts */
  int *cnt = calloc(NT, sizeof(int));
  for (int c = 0; c < nc; c++) for (ll t = 0; t < NT; t++) if (set[c][t >> 6] >> (t & 63) & 1) cnt[t]++;
  ll unc = 0; for (ll t = 0; t < NT; t++) if (!cnt[t]) unc++;
  printf("x=%lld targets in [%lld,%lld]: %lld; core states %d; targets covered by no core: %lld\n", X, LO, LIM, NT, nc, unc);
  int mincnt = 1 << 30; for (ll t = 0; t < NT; t++) if (cnt[t] < mincnt) mincnt = cnt[t];
  ll nmin = 0; for (ll t = 0; t < NT; t++) if (cnt[t] == mincnt) nmin++;
  printf("minimum number of covering cores of a target: %d (attained by %lld targets)\n", mincnt, nmin);
  /* 2. S19's greedy list */
  {
    FILE *g = fopen(s19, "r"); char line[512]; int k = 0, missing = 0; uint64_t *u = calloc(W, 8);
    while (fgets(line, sizeof line, g)) {
      ll a, b, d, e; char *p = strstr(line, "core (h0="); if (!p) continue;
      if (sscanf(p, "core (h0=%lld R0=%lld D0=%lld hu=%lld)", &a, &b, &d, &e) != 4) continue;
      int found = -1; for (int c = 0; c < nc; c++) if (core[c][0] == a && core[c][1] == b && core[c][2] == d && core[c][3] == e) { found = c; break; }
      k++; if (found < 0) { missing++; continue; }
      for (ll i = 0; i < W; i++) u[i] |= set[found][i];
    }
    ll left = 0; for (ll t = 0; t < NT; t++) if (!(u[t >> 6] >> (t & 63) & 1)) left++;
    printf("S19 greedy list: %d core states, %d not among the independently generated states; targets left uncovered: %lld\n", k, missing, left);
  }
  /* 3. independent greedy */
  {
    uint64_t *rem = calloc(W, 8); for (ll t = 0; t < NT; t++) rem[t >> 6] |= 1ULL << (t & 63);
    ll left = NT; int rounds = 0; char *used = calloc(nc, 1);
    while (left > 0) {
      int best = -1; ll bc = 0;
      for (int c = 0; c < nc; c++) { if (used[c]) continue; ll s = 0; for (ll i = 0; i < W; i++) s += __builtin_popcountll(set[c][i] & rem[i]); if (s > bc) { bc = s; best = c; } }
      if (best < 0) break; used[best] = 1; rounds++;
      for (ll i = 0; i < W; i++) rem[i] &= ~set[best][i]; left -= bc;
    }
    printf("independent greedy cover: %d core states (uncovered left %lld)\n", rounds, left);
  }
  /* 4. disjoint covering sets lower bound: targets sorted by cnt ascending, take if disjoint from used cores */
  {
    ll *ord = malloc(NT * sizeof(ll)); ll m = 0;
    for (int c0 = 1; c0 <= nc; c0++) for (ll t = 0; t < NT; t++) if (cnt[t] == c0) ord[m++] = t;
    char *used = calloc(nc, 1); ll *chosen = malloc(NT * sizeof(ll)); ll nch = 0;
    for (ll i = 0; i < m; i++) {
      ll t = ord[i]; int ok = 1;
      for (int c = 0; c < nc && ok; c++) if ((set[c][t >> 6] >> (t & 63) & 1) && used[c]) ok = 0;
      if (!ok) continue;
      for (int c = 0; c < nc; c++) if (set[c][t >> 6] >> (t & 63) & 1) used[c] = 1;
      chosen[nch++] = t;
    }
    /* re-verify pairwise disjointness explicitly */
    ll viol = 0;
    for (ll a = 0; a < nch; a++) for (ll b = a + 1; b < nch; b++)
      for (int c = 0; c < nc; c++) if ((set[c][chosen[a] >> 6] >> (chosen[a] & 63) & 1) && (set[c][chosen[b] >> 6] >> (chosen[b] & 63) & 1)) { viol++; break; }
    printf("disjoint-covering-set lower bound: %lld targets (pairwise intersections found: %lld):", nch, viol);
    for (ll a = 0; a < nch; a++) printf(" %lld(%d)", tval(chosen[a]), cnt[chosen[a]]);
    printf("\n");
  }
  /* 5. sparse rows for the LP / MILP bound */
  {
    FILE *g = fopen(sparse, "w"); ll rows = 0;
    for (ll t = 0; t < NT; t++) if (cnt[t] <= KMAX) {
      fprintf(g, "%lld", tval(t));
      for (int c = 0; c < nc; c++) if (set[c][t >> 6] >> (t & 63) & 1) fprintf(g, " %d", c);
      fprintf(g, "\n"); rows++;
    }
    fclose(g);
    printf("wrote %lld targets with at most %d covering cores to the sparse file\n", rows, KMAX);
  }
  return 0;
}
