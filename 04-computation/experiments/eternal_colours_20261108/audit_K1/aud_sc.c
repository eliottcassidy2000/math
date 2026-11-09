/* aud_sc.c - audit K1 of THM-4930 (2026-11-08). Written from scratch; shares no code with the packet.
 *
 * Census of the Hamiltonian-path counts of self-converse tournaments of order n via the anti-invariant model.
 * Convention (deliberately different from the packet's (i, i+5) pairing):
 *   sigma(v) = v ^ 1 for v < 2m (pairs {0,1}, {2,3}, ...), and sigma(n-1) = n-1 when n = 2m+1 is odd.
 * A labelled tournament T is sigma-anti-invariant iff  u -> v  <=>  sigma(v) -> sigma(u)  for all u != v.
 * The arc orbits are computed generically: the unordered pair {u,v} is sent to {sigma u, sigma v}; one bit per orbit
 * orients a representative u -> v, and the image arc sigma(v) -> sigma(u) is then forced.
 * Reduction (mode "red"): conjugating by the transposition (2k 2k+1) commutes with sigma and reverses the arc inside
 * pair k, so every isomorphism class is met with all within-pair arcs oriented 2k -> 2k+1. Mode "full" enumerates
 * every orbit bit.
 * H by a PULL subset DP (F[S][v] = number of paths covering S ending at v = sum over in-neighbours u of v in S-v of
 * F[S-v][u]); strongness by Landau's criterion on the sorted score sequence (strong iff every proper prefix sum
 * exceeds C(k,2)). With CHECKANTI every generated tournament is also tested for tournament-ness and anti-invariance.
 *
 * usage: aud_sc n red|full [C0 C1]   -> lines "h count strongcount" (count = labelled tournaments in the range)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef uint32_t u32; typedef uint64_t u64;
#define MAXN 12
#define HCAP (1 << 18)
static int n, sig[MAXN];
static int oi[128], oj[128], nother; /* non-pair orbit representatives */
static int pi_[16], pj_[16], npair;   /* within-pair orbits {2k, 2k+1} */
static uint16_t outm[MAXN], inm[MAXN];
static u32 *F;
static void setarc(int a, int b){ outm[a] |= (uint16_t)(1u << b); inm[b] |= (uint16_t)(1u << a); }
static void build(u64 code){
  memset(outm, 0, sizeof outm); memset(inm, 0, sizeof inm);
  for (int k = 0; k < nother; k++){
    int a = oi[k], b = oj[k]; if (!((code >> k) & 1)){ int t = a; a = b; b = t; }
    setarc(a, b); setarc(sig[b], sig[a]);
  }
  for (int k = 0; k < npair; k++){
    int a = pi_[k], b = pj_[k]; if (!((code >> (nother + k)) & 1)){ int t = a; a = b; b = t; }
    setarc(a, b);
  }
}
static int checkanti(void){
  for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) if (u != v){
    int uv = outm[u] >> v & 1, vu = outm[v] >> u & 1;
    if (uv + vu != 1) return 0;
    int img = outm[sig[v]] >> sig[u] & 1;
    if (uv != img) return 0;
  }
  return 1;
}
static u32 hpaths(void){
  int FN = 1 << n;
  for (int v = 0; v < n; v++) F[(1u << v) * n + v] = 1;
  for (int S = 3; S < FN; S++){
    if (!(S & (S - 1))) continue;
    for (int v = 0; v < n; v++) if (S >> v & 1){
      int R = S & ~(1 << v); uint16_t m = inm[v] & (uint16_t)R; u32 s = 0;
      while (m){ int u = __builtin_ctz(m); m &= m - 1; s += F[(u32)R * n + u]; }
      F[(u32)S * n + v] = s;
    }
  }
  u32 h = 0; for (int v = 0; v < n; v++) h += F[(u32)(FN - 1) * n + v]; return h;
}
static int landau_strong(void){
  if (n == 1) return 1;
  int sc[MAXN]; for (int v = 0; v < n; v++) sc[v] = __builtin_popcount(outm[v]);
  for (int i = 1; i < n; i++){ int x = sc[i], j = i - 1; while (j >= 0 && sc[j] > x){ sc[j + 1] = sc[j]; j--; } sc[j + 1] = x; }
  int s = 0; for (int k = 1; k < n; k++){ s += sc[k - 1]; if (s == k * (k - 1) / 2) return 0; }
  return 1;
}
int main(int argc, char **argv){
  if (argc < 3){ fprintf(stderr, "usage: aud_sc n red|full [C0 C1]\n"); return 1; }
  n = atoi(argv[1]); int red = !strcmp(argv[2], "red");
  int m = n / 2;
  for (int v = 0; v < n; v++) sig[v] = (v < 2 * m) ? (v ^ 1) : v;
  static int vis[MAXN][MAXN];
  nother = 0; npair = 0;
  for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++){
    if (vis[u][v]) continue;
    int su = sig[u], sv = sig[v]; int a = su < sv ? su : sv, b = su < sv ? sv : su;
    vis[u][v] = 1; vis[a][b] = 1;
    if (a == u && b == v){ /* orbit of size 1: must be a sigma 2-cycle */
      if (!(sig[u] == v)){ fprintf(stderr, "unexpected fixed pair\n"); return 1; }
      pi_[npair] = u; pj_[npair] = v; npair++;
    } else { oi[nother] = u; oj[nother] = v; nother++; }
  }
  int nb = nother + npair;
  fprintf(stderr, "n=%d: %d orbits (%d within-pair), mode %s\n", n, nb, npair, red ? "red" : "full");
  u64 total = red ? (1ULL << nother) : (1ULL << nb);
  u64 C0 = argc > 3 ? strtoull(argv[3], 0, 10) : 0, C1 = argc > 4 ? strtoull(argv[4], 0, 10) : total;
  F = calloc((size_t)(1u << n) * n, sizeof(u32));
  u64 *cnt = calloc(HCAP, sizeof(u64)), *scnt = calloc(HCAP, sizeof(u64));
  u64 pairmask = ((1ULL << npair) - 1) << nother;
  for (u64 c = C0; c < C1; c++){
    u64 code = red ? (c | pairmask) : c;
    build(code);
#ifdef CHECKANTI
    if (!checkanti()){ fprintf(stderr, "ANTI-INVARIANCE FAILURE at code %llu\n", (unsigned long long)code); return 2; }
#endif
    u32 h = hpaths();
    if (h >= HCAP){ fprintf(stderr, "H overflow\n"); return 3; }
    cnt[h]++; if (landau_strong()) scnt[h]++;
  }
  for (u32 h = 0; h < HCAP; h++) if (cnt[h]) printf("%u %llu %llu\n", h, (unsigned long long)cnt[h], (unsigned long long)scnt[h]);
  return 0;
}
