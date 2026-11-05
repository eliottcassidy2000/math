/*
 * join_witness_check_20261104.c  (opus-2026-11-04-S78)
 * End-to-end check of the two witness families of THM-4875 (D) on actual tournaments, by direct Hamiltonian-path counts.
 *   F1: J = T(w1) => T(w2)^op with x2 -> x1 reversed:  H(J) = (1+u)(1+v) + g(u) g(v) - e(u) e(v)   (v = val w1, u = val w2)
 *   F2: J = T(w) => S^op with s -> x reversed, S any tournament with root s:  H(J) = g(v) y_h - e(v) y_e + (1+v) y_R,
 *       y = (h, e_s, R_s)(S).
 * T(w): chain 1 -> ... -> n plus apex x with x -> i iff bit i of w (most significant first) is 1.
 * Usage: jwc TRIALS seed
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef long long ll;
#define MX 20
static ll gfun(ll v) { if (v == 0) return 1; if (v & 1) return gfun(v >> 1); ll u = v >> 1; return gfun(u) + u + (u & 1); }
static ll *dpb, *stb, *Hsb;
static void triple(int n, const uint32_t *out, int r, ll *h, ll *e, ll *R) {
    int full = (1 << n) - 1; size_t sz = (size_t)(full + 1) * n;
    memset(dpb, 0, sz * sizeof(ll)); memset(stb, 0, sz * sizeof(ll));
    for (int x = 0; x < n; x++) { dpb[(size_t)(1 << x) * n + x] = 1; stb[(size_t)(1 << x) * n + x] = 1; }
    uint32_t inn[MX]; for (int v = 0; v < n; v++) { inn[v] = 0; for (int w = 0; w < n; w++) if (out[w] >> v & 1) inn[v] |= 1u << w; }
    for (int S = 1; S <= full; S++) for (int x = 0; x < n; x++) {
        ll c = dpb[(size_t)S * n + x]; if (c) { uint32_t nb = out[x] & ~(uint32_t)S; while (nb) { int y = __builtin_ctz(nb); nb &= nb - 1; dpb[(size_t)(S | (1 << y)) * n + y] += c; } }
        ll d = stb[(size_t)S * n + x]; if (d) { uint32_t nb = inn[x] & ~(uint32_t)S; while (nb) { int y = __builtin_ctz(nb); nb &= nb - 1; stb[(size_t)(S | (1 << y)) * n + y] += d; } } }
    Hsb[0] = 1; for (int S = 1; S <= full; S++) { ll t = 0; for (int x = 0; x < n; x++) t += dpb[(size_t)S * n + x]; Hsb[S] = t; }
    *h = Hsb[full]; *e = dpb[(size_t)full * n + r]; ll t = 0;
    for (int S = 0; S <= full; S++) if (!(S >> r & 1)) t += Hsb[S] * stb[(size_t)(full & ~S) * n + r];
    *R = t;
}
static ll hcount(int n, const uint32_t *out) { ll h, e, R; triple(n, out, 0, &h, &e, &R); return h; }
/* one-apex T(w) on vertices base..base+L (apex = base + L) inside out[] */
static void oneapex(uint32_t *out, int base, int L, ll w) {
    for (int i = 0; i < L; i++) for (int j = i + 1; j < L; j++) out[base + i] |= 1u << (base + j);
    int x = base + L;
    for (int i = 0; i < L; i++) { int bit = (w >> (L - 1 - i)) & 1; if (bit) out[x] |= 1u << (base + i); else out[base + i] |= 1u << x; }
}
static uint64_t rs = 0x9E3779B97F4A7C15ULL;
static uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
int main(int argc, char **argv) {
    long T = argc > 1 ? atol(argv[1]) : 20000; if (argc > 2) rs ^= (uint64_t)atoll(argv[2]) * 0x2545F4914F6CDD1DULL;
    size_t big = ((size_t)1 << 16) * 16;
    dpb = malloc(big * sizeof(ll)); stb = malloc(big * sizeof(ll)); Hsb = malloc(sizeof(ll) << 16);
    long f1 = 0, f2 = 0, n1 = 0, n2 = 0; ll maxH = 0;
    for (long t = 0; t < T; t++) {
        /* F2 */
        { int L = (int)(rnd() % 10), m = 1 + (int)(rnd() % 6); ll w = L ? (ll)(rnd() % (1ULL << L)) : 0;
          uint32_t S[MX]; memset(S, 0, sizeof S);
          for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { if (rnd() & 1) S[i] |= 1u << j; else S[j] |= 1u << i; }
          int s = (int)(rnd() % m); ll yh, ye, yR; triple(m, S, s, &yh, &ye, &yR);
          int n = L + 1 + m; uint32_t o[MX]; memset(o, 0, sizeof o);
          oneapex(o, 0, L, w); int x = L;
          for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) if (S[i] >> j & 1) o[L + 1 + j] |= 1u << (L + 1 + i);   /* S^op */
          for (int a = 0; a <= L; a++) for (int b = 0; b < m; b++) { if (a == x && b == s) o[L + 1 + b] |= 1u << a; else o[a] |= 1u << (L + 1 + b); }
          ll H = hcount(n, o), v = w, ev = (v == 0 || !(v & 1));
          ll pred = gfun(v) * yh - ev * ye + (1 + v) * yR; n2++; if (H != pred) f2++; if (H > maxH) maxH = H; }
        /* F1 */
        { int L1 = (int)(rnd() % 7), L2 = (int)(rnd() % 7); ll w1 = L1 ? (ll)(rnd() % (1ULL << L1)) : 0, w2 = L2 ? (ll)(rnd() % (1ULL << L2)) : 0;
          int n = L1 + 1 + L2 + 1; uint32_t o[MX], B[MX]; memset(o, 0, sizeof o); memset(B, 0, sizeof B);
          oneapex(o, 0, L1, w1); int x1 = L1;
          oneapex(B, 0, L2, w2); int x2 = L2;
          for (int i = 0; i <= L2; i++) for (int j = 0; j <= L2; j++) if (B[i] >> j & 1) o[L1 + 1 + j] |= 1u << (L1 + 1 + i);   /* T(w2)^op */
          for (int a = 0; a <= L1; a++) for (int b = 0; b <= L2; b++) { if (a == x1 && b == x2) o[L1 + 1 + b] |= 1u << a; else o[a] |= 1u << (L1 + 1 + b); }
          ll H = hcount(n, o), v = w1, u = w2, ev = (v == 0 || !(v & 1)), eu = (u == 0 || !(u & 1));
          ll pred = (1 + u) * (1 + v) + gfun(u) * gfun(v) - eu * ev; n1++; if (H != pred) f1++; if (H > maxH) maxH = H; }
    }
    printf("# join_witness_check_20261104: F1 trials %ld failures %ld | F2 trials %ld failures %ld | largest H %lld\n", n1, f1, n2, f2, maxH);
    return 0;
}
