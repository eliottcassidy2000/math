/* mirror_sample.c - opus-2026-11-09-S87 (PREREGISTRATION-4, HYP-9540 (D)).
 * Blue tilings at odd n = 2c + 1 (path 0 -> ... -> n-1; tile (i, j), j >= i + 2, reversed iff j -> i;
 * sigma(i, j) = (n-1-j, n-1-i)).  Coordinates: x = the c mirror tiles (i, n-1-i), y = one bit per 2-element sigma-orbit
 * of tiles (N orbits).  By THM-4922 (4) (H - 1)/2 = A(y) + x.B(y) (mod 2), with N(u) the number of mirror-free
 * directed paths c -> ... -> u (each mirror pair {i, n-1-i} met at most once), B_i = N(i) + N(n-1-i),
 * A = sum_{i < c} N(i).  N(u) mod 2 is computed by a DP over the sets of used mirror pairs.
 * For each y: y is in D iff B(y) = B(y + 1) (y + 1 = all y-bits flipped); then e(y) = A(y) + A(y + 1) + 1.
 * The (1,3)-share of blue lines is 1/2 + (|D_1| - |D_0|) / 2^(N + 1)  (THM-4951 (5d)).
 * usage: mirror_sample n MODE SAMPLES SEED SLICE NSLICES
 *        MODE 0: exhaustive over the y with index = SLICE (mod NSLICES); MODE 1: SAMPLES uniform random y (seeded).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
typedef unsigned long long u64;
#define MAXN 17
static int n, c, m, ti[200], tj[200], np, oa[100], ob[100];
static uint32_t cur[1 << 8];
static void adjacency(u64 y, uint32_t *out){
  memset(out, 0, sizeof(uint32_t) * MAXN);
  for (int i = 0; i + 1 < n; i++) out[i] |= 1u << (i + 1);
  /* all non-path tiles forward first, then reverse the y-orbits that are set; mirror tiles stay forward (x = 0) */
  for (int b = 0; b < m; b++) out[ti[b]] |= 1u << tj[b];
  for (int q = 0; q < np; q++) if (y >> q & 1){ int bs[2] = {oa[q], ob[q]};
    for (int z = 0; z < 2; z++){ int b = bs[z]; out[ti[b]] &= ~(1u << tj[b]); out[tj[b]] |= 1u << ti[b]; } }
}
static void parities(const uint32_t *out, int *Np){        /* Np[u] = N(u) mod 2 */
  int P = 1 << c;
  memset(cur, 0, sizeof(uint32_t) * P);
  uint32_t o = out[c];
  while (o){ int w = __builtin_ctz(o); o &= o - 1; int p = w < c ? w : n - 1 - w; cur[1 << p] ^= 1u << w; }
  for (int mask = 1; mask < P; mask++){ uint32_t s = cur[mask];
    while (s){ int u = __builtin_ctz(s); s &= s - 1; uint32_t nb = out[u] & ~(1u << c);
      while (nb){ int w = __builtin_ctz(nb); nb &= nb - 1; int p = w < c ? w : n - 1 - w; if (mask >> p & 1) continue;
        cur[mask | (1 << p)] ^= 1u << w; } } }
  for (int u = 0; u < n; u++) Np[u] = 0;
  for (int mask = 1; mask < P; mask++){ uint32_t s = cur[mask]; while (s){ int u = __builtin_ctz(s); s &= s - 1; Np[u] ^= 1; } }
}
static void AB(u64 y, int *A, uint32_t *B){
  uint32_t out[MAXN]; int Np[MAXN]; adjacency(y, out); parities(out, Np);
  *A = 0; *B = 0; for (int i = 0; i < c; i++){ *A ^= Np[i]; if (Np[i] ^ Np[n - 1 - i]) *B |= 1u << i; }
}
static u64 rs;
static u64 rnd(void){ rs ^= rs >> 12; rs ^= rs << 25; rs ^= rs >> 27; return rs * 0x2545F4914F6CDD1DULL; }
int main(int argc, char **argv){
  n = atoi(argv[1]); int mode = atoi(argv[2]); long long samples = atoll(argv[3]); rs = 0x9E3779B97F4A7C15ULL ^ (u64)atoll(argv[4]);
  int slice = atoi(argv[5]), ns = atoi(argv[6]); c = n / 2;
  int idx[MAXN][MAXN]; m = 0;
  for (int i = 0; i < n; i++) for (int j = i + 2; j < n; j++){ ti[m] = i; tj[m] = j; idx[i][j] = m; m++; }
  np = 0; for (int b = 0; b < m; b++){ int s = idx[n - 1 - tj[b]][n - 1 - ti[b]]; if (b < s){ oa[np] = b; ob[np] = s; np++; } }
  u64 full = np == 64 ? ~0ULL : ((1ULL << np) - 1);
  long long cnt = 0, d = 0, d1 = 0;
  u64 next = (u64)slice;
  rnd(); rnd();
  for (;;){
    u64 y;
    if (mode == 0){ if (next > full) break; y = next; next += (u64)ns; }
    else { if (cnt >= samples) break; y = rnd() & full; }
    int A, Ap; uint32_t B, Bp; AB(y, &A, &B); AB(y ^ full, &Ap, &Bp);
    cnt++;
    if (B == Bp){ d++; if ((A ^ Ap ^ 1) == 1) d1++; }
  }
  double z = (double)(d1 - (d - d1)) / cnt, ez2 = (double)d / cnt;
  double se = sqrt((ez2 - z * z) / cnt) / 2;
  printf("n=%d c=%d N=%d mode=%d slice %d/%d: y %lld, |D| %lld, |D_1| %lld, |D_0| %lld; share = 1/2 + %.6f (se %.6f)\n",
         n, c, np, mode, slice, ns, cnt, d, d1, d - d1, z / 2, se);
  return 0;
}
