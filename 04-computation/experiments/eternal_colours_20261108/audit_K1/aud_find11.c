/* aud_find11.c - audit K1 (2026-11-08). Explicit self-converse 11-tournaments for given H values (here 1551 and 1627,
 * found in the auditor's complete order-11 census, out/o11_full.txt). Since H(T - 10) <= H(T), it scans the reduced
 * order-10 codes of aud_sc.c with H < max target, and their 32 extensions by the fixed vertex 10 (10 -> v iff
 * sigma(v) -> 10). For the first tournament found with each target H it re-checks anti-invariance under
 * sigma = (0 1)(2 3)...(8 9), recomputes H, and prints the out-neighbour lists and the upper-triangle code (bit b for the
 * b-th pair i<j, set iff i -> j; THM-4922's census convention).
 * usage: aud_find11 h1 h2 ...
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef uint32_t u32; typedef uint64_t u64;
#define MAXN 12
static int n, sig[MAXN], oi[128], oj[128], nother, pi_[16], pj_[16], npair;
static uint16_t outm[MAXN], inm[MAXN]; static u32 F[(1 << 11) * 11];
static void setarc(int a, int b){ outm[a] |= (uint16_t)(1u << b); inm[b] |= (uint16_t)(1u << a); }
static void orbits10(void){
  n = 10; for (int v = 0; v < 10; v++) sig[v] = v ^ 1;
  static int vis[MAXN][MAXN]; memset(vis, 0, sizeof vis); nother = npair = 0;
  for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++){
    if (vis[u][v]) continue;
    int su = sig[u], sv = sig[v]; int a = su < sv ? su : sv, b = su < sv ? sv : su;
    vis[u][v] = 1; vis[a][b] = 1;
    if (a == u && b == v){ pi_[npair] = u; pj_[npair] = v; npair++; } else { oi[nother] = u; oj[nother] = v; nother++; }
  }
}
static void build10(u64 c){
  u64 code = c | (((1ULL << npair) - 1) << nother);
  memset(outm, 0, sizeof outm); memset(inm, 0, sizeof inm);
  for (int k = 0; k < nother; k++){ int a = oi[k], b = oj[k]; if (!((code >> k) & 1)){ int t = a; a = b; b = t; } setarc(a, b); setarc(sig[b], sig[a]); }
  for (int k = 0; k < npair; k++){ int a = pi_[k], b = pj_[k]; if (!((code >> (nother + k)) & 1)){ int t = a; a = b; b = t; } setarc(a, b); }
}
static u32 hpaths(int nn){
  int FN = 1 << nn;
  for (int v = 0; v < nn; v++) F[(1u << v) * nn + v] = 1;
  for (int S = 3; S < FN; S++){ if (!(S & (S - 1))) continue;
    for (int v = 0; v < nn; v++) if (S >> v & 1){ int R = S & ~(1 << v); uint16_t m = inm[v] & (uint16_t)R; u32 s = 0;
      while (m){ int u = __builtin_ctz(m); m &= m - 1; s += F[(u32)R * nn + u]; } F[(u32)S * nn + v] = s; } }
  u32 h = 0; for (int v = 0; v < nn; v++) h += F[(u32)(FN - 1) * nn + v]; return h;
}
int main(int argc, char **argv){
  int nt = argc - 1; long tg[16]; int done[16] = {0}; long mx = 0;
  for (int i = 0; i < nt; i++){ tg[i] = atol(argv[i + 1]); if (tg[i] > mx) mx = tg[i]; }
  orbits10();
  int left = nt;
  for (u64 c = 0; c < (1ULL << nother) && left; c++){
    orbits10(); build10(c);
    if (hpaths(10) > (u32)mx) continue;
    uint16_t o10[MAXN], i10[MAXN]; memcpy(o10, outm, sizeof o10); memcpy(i10, inm, sizeof i10);
    for (int e = 0; e < 32 && left; e++){
      memcpy(outm, o10, sizeof o10); memcpy(inm, i10, sizeof i10); outm[10] = inm[10] = 0;
      for (int k = 0; k < 5; k++){            /* orbit {2k, 10} ~ {2k+1, 10}: 10 -> 2k iff 2k+1 -> 10 */
        if (e >> k & 1){ setarc(10, 2 * k); setarc(2 * k + 1, 10); } else { setarc(2 * k, 10); setarc(10, 2 * k + 1); }
      }
      u32 h = hpaths(11);
      for (int i = 0; i < nt; i++) if (!done[i] && h == (u32)tg[i]){
        int sg[11]; for (int v = 0; v < 11; v++) sg[v] = v < 10 ? (v ^ 1) : 10;
        int ok = 1;
        for (int u = 0; u < 11; u++) for (int v = 0; v < 11; v++) if (u != v){
          int uv = outm[u] >> v & 1, vu = outm[v] >> u & 1;
          if (uv + vu != 1 || uv != ((outm[sg[v]] >> sg[u]) & 1)) ok = 0; }
        u64 code = 0; int b = 0; for (int x = 0; x < 11; x++) for (int y = x + 1; y < 11; y++){ if (outm[x] >> y & 1) code |= 1ULL << b; b++; }
        printf("H = %u, order 11, anti-invariant under (0 1)(2 3)(4 5)(6 7)(8 9): %s; upper-triangle code %llu\n", h, ok ? "yes" : "NO", (unsigned long long)code);
        for (int v = 0; v < 11; v++){ printf("  %2d ->", v); for (int u = 0; u < 11; u++) if (outm[v] >> u & 1) printf(" %d", u); printf("\n"); }
        done[i] = 1; left--;
      }
    }
  }
  return left != 0;
}
