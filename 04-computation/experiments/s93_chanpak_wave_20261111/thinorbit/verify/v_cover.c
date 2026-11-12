/* Verifier (lane thinorbit, S93): independent orbit code, written without reference to orbit_cover.c's design.
 *
 * usage: v_cover MATFILE SUBSET MODE LIM [K]
 *   SUBSET: all | inv (det != 0) | uni (|det| = 1) | sing (det = 0)
 *   MODE exact : exact closure of ALL orbit states (h,E,O) of (1,0,1) with h <= LIM.  No saturation: E and O are stored
 *                exactly ((64 - bits(h))/2 bits each); any overflow is counted and makes the run invalid.
 *                Prints the cumulative state counts at powers of two, the list of odd v <= LIM that are not first
 *                coordinates, and the statistics of r(v) = #states with first coordinate v on dyadic windows.
 *   MODE capped: forward search keeping at most K states per value of h (first come, exact states, no saturation,
 *                64-bit E and O), processed in increasing h.  Every value marked is a TRUE orbit value (the stored
 *                states are exact images of exact states); the miss list is only an upper bound for the true misses.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef long long ll;
static ll M[4096][9];
static int nm;
static ll det3(const ll *m) {
  return m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]);
}
static uint64_t *tab, tsize, tcount;
static uint64_t mix(uint64_t x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x; }
static void regrow(void) {
  uint64_t os = tsize, *ot = tab;
  tsize = os ? os * 2 : (1ULL << 20);
  tab = calloc(tsize, 8);
  if (!tab) { fprintf(stderr, "oom\n"); exit(1); }
  for (uint64_t i = 0; i < os; i++) if (ot[i]) { uint64_t j = mix(ot[i]) & (tsize - 1); while (tab[j]) j = (j + 1) & (tsize - 1); tab[j] = ot[i]; }
  free(ot);
}
static int add(uint64_t key) {
  if (2 * (tcount + 1) > tsize) regrow();
  uint64_t j = mix(key) & (tsize - 1);
  while (tab[j]) { if (tab[j] == key) return 0; j = (j + 1) & (tsize - 1); }
  tab[j] = key; tcount++; return 1;
}
int main(int argc, char **argv) {
  if (argc < 5) { fprintf(stderr, "usage\n"); return 1; }
  FILE *f = fopen(argv[1], "r");
  if (!f) { fprintf(stderr, "cannot open\n"); return 1; }
  const char *sub = argv[2];
  ll LIM = atoll(argv[4]);
  int K = argc > 5 ? atoi(argv[5]) : 16;
  ll t[9]; int nread = 0;
  while (fscanf(f, "%lld %lld %lld %lld %lld %lld %lld %lld %lld", t, t + 1, t + 2, t + 3, t + 4, t + 5, t + 6, t + 7, t + 8) == 9) {
    nread++;
    ll d = det3(t);
    int keep = !strcmp(sub, "all") || (!strcmp(sub, "inv") && d != 0) || (!strcmp(sub, "uni") && (d == 1 || d == -1)) || (!strcmp(sub, "sing") && d == 0);
    if (keep) { memcpy(M[nm], t, sizeof(t)); nm++; }
    int c; while ((c = fgetc(f)) != '\n' && c != EOF) {}
  }
  fclose(f);
  int minq = 1 << 30; for (int k = 0; k < nm; k++) { if (M[k][1] < minq) minq = (int)M[k][1]; for (int j = 0; j < 9; j++) if (M[k][j] < 0) { fprintf(stderr, "negative entry\n"); return 1; } if (M[k][0] < 1) { fprintf(stderr, "eta < 1\n"); return 1; } }
  printf("file %s: %d matrices read, subset %s: %d generators; min coefficient of E in the h-row: %d; LIM %lld\n", argv[1], nread, sub, nm, minq, LIM);
  unsigned char *hit = calloc(LIM + 2, 1);
  if (!strcmp(argv[3], "exact")) {
    int bh = 1; while ((1LL << bh) <= LIM) bh++;
    int BE = (64 - bh) / 2; if (BE > 31) BE = 31;
    const int SH2 = 2 * BE, SH1 = BE;
    uint64_t qcap = 1 << 20, qlen = 0, *queue = malloc(qcap * 8), overflow = 0;
    const uint64_t m23 = (1ULL << BE) - 1;
    regrow();
    uint64_t s0 = (1ULL << SH2) | 1ULL;
    add(s0); queue[qlen++] = s0;
    ll maxE = 0, maxO = 0;
    for (uint64_t qi = 0; qi < qlen; qi++) {
      uint64_t s = queue[qi];
      ll h = (ll)(s >> SH2), e = (ll)((s >> SH1) & m23), o = (ll)(s & m23);
      for (int k = 0; k < nm; k++) {
        const ll *m = M[k];
        ll h2 = m[0] * h + m[1] * e + m[2] * o;
        if (h2 > LIM) continue;
        ll e2 = m[3] * h + m[4] * e + m[5] * o, o2 = m[6] * h + m[7] * e + m[8] * o;
        if (e2 > (ll)m23 || o2 > (ll)m23) { overflow++; continue; }
        if (e2 > maxE) maxE = e2;
        if (o2 > maxO) maxO = o2;
        uint64_t key = ((uint64_t)h2 << SH2) | ((uint64_t)e2 << SH1) | (uint64_t)o2;
        if (add(key)) { if (qlen == qcap) { qcap *= 2; queue = realloc(queue, qcap * 8); if (!queue) { fprintf(stderr, "oom\n"); return 1; } } queue[qlen++] = key; }
      }
    }
    printf("exact closure: %llu distinct states with h <= %lld; overflow events %llu (must be 0); max E %lld, max O %lld\n",
           (unsigned long long)qlen, LIM, (unsigned long long)overflow, maxE, maxO);
    uint32_t *r = calloc(LIM + 2, 4);
    for (uint64_t qi = 0; qi < qlen; qi++) r[queue[qi] >> SH2]++;
    uint64_t cum = 0; ll v = 1;
    printf("cumulative #states with h <= 2^k:");
    for (int k = 0; (1LL << k) <= LIM; k++) { for (; v <= (1LL << k); v++) cum += r[v]; printf(" %d:%llu", k, (unsigned long long)cum); }
    printf("\n");
    ll nmiss = 0; printf("odd v <= %lld with r(v) = 0:", LIM);
    for (ll x = 1; x <= LIM; x += 2) if (!r[x]) { nmiss++; if (nmiss <= 150) printf(" %lld", x); }
    printf("\n#misses %lld; even v with r > 0: ", nmiss);
    { ll ev = 0; for (ll x = 2; x <= LIM; x += 2) if (r[x]) ev++; printf("%lld\n", ev); }
    printf("window k | #odd | #r=0 | min r (v) | mean r | max r | CS share (sum r)^2/(#odd sum r^2)\n");
    for (int k = 1; (1LL << k) <= LIM; k++) {
      ll lo = 1LL << k, hi = (2LL << k) - 1; if (hi > LIM) hi = LIM;
      ll n = 0, z = 0, mn = 1LL << 60, mnv = 0, mx = 0; double s1 = 0, s2 = 0;
      for (ll x = lo + 1; x <= hi; x += 2) { n++; ll c = r[x]; if (!c) z++; if (c < mn) { mn = c; mnv = x; } if (c > mx) mx = c; s1 += c; s2 += (double)c * c; }
      if (n) printf("  %2d | %lld | %lld | %lld (%lld) | %.2f | %lld | %.4f\n", k, n, z, mn, mnv, s1 / n, mx, s2 > 0 ? s1 * s1 / (n * s2) : 0.0);
    }
  } else {
    ll nv = LIM / 2 + 1;
    int *cnt = calloc(nv, sizeof(int));
    ll *SE = malloc(sizeof(ll) * nv * K), *SO = malloc(sizeof(ll) * nv * K);
    if (!cnt || !SE || !SO) { fprintf(stderr, "oom\n"); return 1; }
    cnt[0] = 1; SE[0] = 0; SO[0] = 1; hit[1] = 1;
    uint64_t stored = 1, apps = 0;
    int prune = minq >= 2;
    for (ll v = 1; v <= LIM; v += 2) {
      ll iv = v / 2;
      for (int i = 0; i < cnt[iv]; i++) {
        ll e = SE[iv * K + i], o = SO[iv * K + i];
        for (int k = 0; k < nm; k++) {
          const ll *m = M[k];
          ll h2 = m[0] * v + m[1] * e + m[2] * o;
          if (h2 > LIM) continue;
          apps++;
          if (!(h2 & 1)) { fprintf(stderr, "even h reached\n"); return 1; }
          ll e2 = m[3] * v + m[4] * e + m[5] * o, o2 = m[6] * v + m[7] * e + m[8] * o;
          hit[h2] = 1;
          if (prune && h2 + 2 * e2 > LIM) continue;  /* every child has h >= h2 + 2 e2 > LIM */
          ll j = h2 / 2; int c = cnt[j];
          if (c >= K) continue;
          int dup = 0; for (int a = 0; a < c; a++) if (SE[j * K + a] == e2 && SO[j * K + a] == o2) { dup = 1; break; }
          if (dup) continue;
          SE[j * K + c] = e2; SO[j * K + c] = o2; cnt[j] = c + 1; stored++;
        }
      }
    }
    printf("capped search K = %d: %llu states stored, %llu images with h <= LIM\n", K, (unsigned long long)stored, (unsigned long long)apps);
    ll nmiss = 0, last = -1; printf("odd v <= %lld not marked:", LIM);
    for (ll x = 1; x <= LIM; x += 2) if (!hit[x]) { nmiss++; last = x; if (nmiss <= 150) printf(" %lld", x); }
    printf("\n#not marked %lld, largest %lld\n", nmiss, last);
  }
  return 0;
}
