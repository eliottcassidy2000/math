/* aud_witness16.c - audit K2 (2026-11-08): audit K1's aud_witness.c with MAXN = 16 (orders up to 16, for sc_walk.out). Written from scratch.
 * Re-checks every line "W h n strong rev i-j ..." of the packet's sc_witness.out:
 *   build TT_n (i -> j for i < j), reverse each listed arc i-j and its mirror (n-1-j, n-1-i);
 *   check (1) the result is a tournament, (2) mu(v) = n-1-v is an anti-automorphism (u -> v iff mu(v) -> mu(u)),
 *   which proves self-converseness, (3) H (pull subset DP) equals h, (4) Landau strongness equals the flag.
 * Additionally, for witnesses with n <= 9 (and a sample above), self-converseness is confirmed by an independent
 * backtracking anti-isomorphism search that does not use mu.
 * usage: aud_witness < sc_witness.out
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef uint32_t u32;
#define MAXN 16
static int n; static uint16_t outm[MAXN], inm[MAXN]; static u32 F[(1 << MAXN) * MAXN];
static u32 hpaths(void){
  int FN = 1 << n;
  for (int v = 0; v < n; v++) F[(1u << v) * n + v] = 1;
  for (int S = 3; S < FN; S++){ if (!(S & (S - 1))) continue;
    for (int v = 0; v < n; v++) if (S >> v & 1){ int R = S & ~(1 << v); uint16_t m = inm[v] & (uint16_t)R; u32 s = 0;
      while (m){ int u = __builtin_ctz(m); m &= m - 1; s += F[(u32)R * n + u]; } F[(u32)S * n + v] = s; } }
  u32 h = 0; for (int v = 0; v < n; v++) h += F[(u32)(FN - 1) * n + v]; return h;
}
static int landau_strong(void){
  if (n == 1) return 1;
  int sc[MAXN]; for (int v = 0; v < n; v++) sc[v] = __builtin_popcount(outm[v]);
  for (int i = 1; i < n; i++){ int x = sc[i], j = i - 1; while (j >= 0 && sc[j] > x){ sc[j + 1] = sc[j]; j--; } sc[j + 1] = x; }
  int s = 0; for (int k = 1; k < n; k++){ s += sc[k - 1]; if (s == k * (k - 1) / 2) return 0; }
  return 1;
}
static void flip(int i, int j){ /* make j -> i where i -> j was */
  outm[i] &= ~(1u << j); inm[j] &= ~(1u << i); outm[j] |= 1u << i; inm[i] |= 1u << j;
}
static int p[MAXN], used[MAXN], found;
static void bt(int k){
  if (found) return;
  if (k == n){ found = 1; return; }
  int dk = __builtin_popcount(outm[k]);
  for (int w = 0; w < n && !found; w++) if (!used[w]){
    if (__builtin_popcount(outm[w]) != n - 1 - dk) continue;
    int ok = 1;
    for (int j = 0; j < k && ok; j++) if (((outm[j] >> k) & 1) != ((outm[w] >> p[j]) & 1)) ok = 0;
    if (!ok) continue;
    p[k] = w; used[w] = 1; bt(k + 1); used[w] = 0;
  }
}
int main(void){
  char line[512]; long nw = 0, bad = 0, btchecked = 0;
  while (fgets(line, sizeof line, stdin)){
    if (line[0] != 'W') continue;
    long h; int nn, st; char *q = line + 1; int off;
    if (sscanf(q, "%ld %d %d rev%n", &h, &nn, &st, &off) != 3){ fprintf(stderr, "parse error: %s", line); bad++; continue; }
    n = nn; q += off;
    memset(outm, 0, sizeof outm); memset(inm, 0, sizeof inm);
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++){ outm[i] |= 1u << j; inm[j] |= 1u << i; }
    int i, j, o; int err = 0;
    while (sscanf(q, " %d-%d%n", &i, &j, &o) == 2){
      q += o;
      if (!(i < j && j < n)){ err = 1; break; }
      int mi = n - 1 - j, mj = n - 1 - i;
      if (outm[i] >> j & 1) flip(i, j); else err = 1;               /* each listed arc reversed once */
      if (!(mi == i && mj == j)){ if (outm[mi] >> mj & 1) flip(mi, mj); else err = 1; }
    }
    /* tournament + mu anti-invariance */
    for (int u = 0; u < n && !err; u++) for (int v = 0; v < n; v++) if (u != v){
      int uv = outm[u] >> v & 1, vu = outm[v] >> u & 1; if (uv + vu != 1){ err = 2; break; }
      if (uv != ((outm[n - 1 - v] >> (n - 1 - u)) & 1)){ err = 3; break; } }
    u32 hh = err ? 0 : hpaths(); int ss = err ? -1 : landau_strong();
    if (!err && (n <= 9 || nw % 25 == 0)){ memset(used, 0, sizeof used); found = 0; bt(0); btchecked++; if (!found) err = 4; }
    if (err || hh != (u32)h || ss != st){ bad++; fprintf(stderr, "BAD (err %d, H %u, strong %d): %s", err, hh, ss, line); }
    nw++;
  }
  printf("witness lines checked: %ld; failures: %ld; independent backtracking self-converse confirmations: %ld\n", nw, bad, btchecked);
  return bad != 0;
}
