/*
 * procgen_seven2_20260926_restrict.c -- Min's energy game with RESTRICTED signs, incremental (warm-started)
 * (lane "seven2", session collatz-procgen-20260922, 2026-09-26).  Shared library driven by ctypes from
 * procgen_seven2_20260926_lib.py.  Exploratory search engine: whatever it reports is re-checked by the
 * caller (the final sign strategy is evaluated by procgen_seven2_20260926_rhomax.c and its certificate by
 * procgen_seven_20260926_verify.c).
 *
 * Pair form (THM-4486, Prop. P of the seven note).  Level k, odd q, H = 2^(k-1); pairs P in Z/H, lifts
 * x_b = P + b H.  Each odd node x carries a mask of allowed signs (bit0: '+', bit1: '-').  Threshold
 * F = fn/fd, e(P) = fd - fn (P odd), -fn (P even).  Least fixed point of
 *     G(P) = max(0, e(P) + max_b min_{s allowed at x_b} G(T_s(x_b) mod H))    (T_s(x) = x/2 or (qx+s)/2)
 * with values above cap set to TOP.  If no pair is TOP, the argmin signs form a strategy within the
 * masks with rho_max <= F (Lemma G1 / Prop. P upper); if some pair is TOP (and the cap is adequate), no
 * such strategy exists -- a cap can only make the answer "infeasible" wrongly, never "feasible" wrongly.
 * Restricting masks can only raise the least fixed point, so after a restriction the iteration is
 * warm-started from the previous fixed point (still a lower bound of the new one).
 */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef int64_t i64;
#define TOPV 0x7fffffff
static int K; static i64 N, H, Q, QINV;
static int32_t *G = 0, *SAVE = 0;
static uint8_t *MASK = 0;          /* per odd node index i = (x-1)/2 */
static uint64_t *DIRTY = 0;
static i64 NW;

int r_init(int k, i64 q) {
    free(G); free(SAVE); free(MASK); free(DIRTY);
    K = k; N = (i64)1 << k; H = N >> 1; Q = q;
    uint64_t x = (uint64_t)q; for (int i = 0; i < 6; i++) x = x * (2ULL - (uint64_t)q * x);
    QINV = (i64)(x & (uint64_t)(N - 1));
    G = (int32_t *)calloc((size_t)H, 4); SAVE = (int32_t *)malloc((size_t)H * 4);
    MASK = (uint8_t *)malloc((size_t)H); memset(MASK, 3, (size_t)H);
    NW = (H + 63) / 64;
    DIRTY = (uint64_t *)calloc((size_t)NW, 8);
    return G && SAVE && MASK && DIRTY ? 0 : -1;
}
void r_reset_values(void) { memset(G, 0, (size_t)H * 4); for (i64 w = 0; w < NW; w++) DIRTY[w] = ~0ULL; if (H % 64) DIRTY[NW - 1] = (1ULL << (H % 64)) - 1; }
void r_save(void) { memcpy(SAVE, G, (size_t)H * 4); }
void r_restore(void) { memcpy(G, SAVE, (size_t)H * 4); memset(DIRTY, 0, (size_t)NW * 8); }
void r_get_mask(uint8_t *out) { memcpy(out, MASK, (size_t)H); }
void r_set_mask(const uint8_t *m) { memcpy(MASK, m, (size_t)H); }
void r_get_values(int32_t *out) { memcpy(out, G, (size_t)H * 4); }
static inline void mark(i64 P) { DIRTY[P >> 6] |= 1ULL << (P & 63); }
/* set the mask of the odd nodes x = c + j*M (M = 2^d, c odd) to m; marks their pairs dirty */
void r_mask_class(i64 c, int d, int m) {
    i64 M = (i64)1 << d;
    for (i64 x = c; x < N; x += M) { MASK[(x - 1) >> 1] = (uint8_t)m; mark(x & (H - 1)); }
}
static inline i64 val_node(i64 x) {          /* min over allowed signs of G(target) */
    if (!(x & 1)) return G[(x >> 1) & (H - 1)];
    int m = MASK[(x - 1) >> 1];
    i64 t = Q * x, best = TOPV;
    if (m & 1) { i64 v = G[((t + 1) >> 1) & (H - 1)]; if (v < best) best = v; }
    if (m & 2) { i64 v = G[((t - 1) >> 1) & (H - 1)]; if (v < best) best = v; }
    return best;
}
/* run to the fixed point; returns the number of TOP pairs */
i64 r_solve(i64 fn, i64 fd, i64 cap) {
    i64 eo = fd - fn, ee = -fn;
    int any = 1;
    while (any) {
        any = 0;
        for (i64 w = 0; w < NW; w++) {
            uint64_t bits = DIRTY[w]; DIRTY[w] = 0;
            while (bits) {
                int bit = __builtin_ctzll(bits); bits &= bits - 1;
                i64 P = w * 64 + bit;
                if (G[P] == TOPV) continue;
                i64 a = val_node(P), b = val_node(P + H);
                i64 m = a > b ? a : b;
                i64 v = (m == TOPV) ? TOPV : m + ((P & 1) ? eo : ee);
                if (v < 0) v = 0;
                if (v > cap) v = TOPV;
                if (v > G[P]) {
                    G[P] = (int32_t)v; any = 1;
                    mark((2 * P) & (H - 1));
                    mark((QINV * ((2 * P - 1) & (H - 1))) & (H - 1));
                    mark((QINV * ((2 * P + 1) & (H - 1))) & (H - 1));
                }
            }
        }
    }
    i64 c = 0;
    for (i64 P = 0; P < H; P++) if (G[P] == TOPV) c++;
    return c;
}
/* extract an argmin strategy (1 = sign '-') within the masks, ties to '+' if allowed */
void r_strategy(uint8_t *minus) {
    for (i64 i = 0; i < H; i++) {
        i64 x = 2 * i + 1, t = Q * x; int m = MASK[i];
        i64 vp = (m & 1) ? G[((t + 1) >> 1) & (H - 1)] : TOPV, vm = (m & 2) ? G[((t - 1) >> 1) & (H - 1)] : TOPV;
        minus[i] = vm < vp ? 1 : 0;
    }
}
