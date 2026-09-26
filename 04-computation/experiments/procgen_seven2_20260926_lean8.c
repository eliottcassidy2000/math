/*
 * procgen_seven2_20260926_lean8.c -- one-sided energy fixed points of the min-max cycle-density game (THM-4486) with
 * ONE BYTE per negation-orbit representative, for levels k = 30, 31 under a 700 MB memory cap
 * (lane "seven2", session collatz-procgen-20260922, 2026-09-26).  Search engine only: every certificate it writes is
 * re-checked by procgen_seven2_20260926_verify8.c (an independent exact checker).
 *
 * This is the "lower"/"upper" mode of the seven lane's solver (procgen_seven_20260926_game.c) with uint8 values:
 *   level k, odd q, N = 2^k, H = 2^(k-1); pairs P in Z/H with lifts x_b = P + b H; threshold F = fn/fd;
 *   e(P) = fd - fn (P odd), -fn (P even); options of x: x/2 mod H (x even), (q x +- 1)/2 mod H (x odd).
 *   MAX (lower):  g(P) = max(0, -e(P) + min_b max_{Q option of x_b} g(Q))
 *   MIN (upper):  G(P) = max(0,  e(P) + max_b min_{Q option of x_b} G(Q))
 * least fixed points from 0, values above cap (<= 254) are TOP (255).  Both operators commute with the negation
 * x -> -x (it maps the lifts of P onto those of -P and the option set of x onto minus that of -x), so the least fixed
 * points are negation-invariant and are computed on the representatives r = min(P, H-P), 0 <= r <= H/2
 * (2^(k-2) + 1 bytes, plus a dirty bitset of 2^(k-2)/8 bytes).
 * Caps only turn finite values into TOP: a small cap can only make a test fail, never produce a false certificate.
 *
 * Offset option (lower mode, argument 'shift'): the Max potential satisfies g(P) >= fn v(P), v(P) = the number of
 * forced halvings from P (the 2-adic valuation of P, and k-1 for P = 0), so the engine stores h = g - fn v(P) in one
 * byte instead of g (values of h above the cap are TOP, which can only shrink W).  The file is then prefix.lo_h8.
 * Output files (the seven lane's compact format, on all H pairs):
 *   lower: prefix.lo_taub (bit P = lift tau(P), argmin with ties to 0), prefix.lo_g8 (g(P), 255 = outside W)
 *   upper: prefix.up_sigb (bit i = 1 iff the odd node 2i+1 takes '-', argmin with ties to '+'), prefix.up_psi8
 * usage: procgen_seven2_20260926_lean8 lower|upper q k fn fd prefix [cap]
 * prints "LOWER-WRITTEN" / "UPPER-WRITTEN" on success, "LOWER-FAILED" / "UPPER-FAILED" otherwise.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

typedef int64_t i64;
#define TOP8 255
static i64 N, H, HH, Q, QINV, NR;
static uint8_t *V;
static uint64_t *DIRTY;
static int MODE, SHIFT = 0, KK;
static i64 FN, FD, CAP;
static inline i64 vv(i64 P) { if (P == 0) return KK - 1; return __builtin_ctzll((unsigned long long)P); }
/* the stored value of a pair P as a potential (TOP stays TOP) */
static inline i64 gval(i64 P) { i64 h = V[P <= HH ? P : H - P]; if (h == TOP8) return 1 << 30; return SHIFT ? h + FN * vv(P) : h; }

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static inline i64 rep(i64 P) { return P <= HH ? P : H - P; }
static inline int options(i64 P, int b, i64 *Qo) {
    i64 x = P + (b ? H : 0);
    if (!(P & 1)) { Qo[0] = (x >> 1) & (H - 1); return 1; }
    i64 t = Q * x;
    Qo[0] = ((t + 1) >> 1) & (H - 1);
    Qo[1] = ((t - 1) >> 1) & (H - 1);
    return 2;
}
static inline int eval(i64 P) {
    i64 e = (P & 1) ? (FD - FN) : -FN;
    i64 Qo[2];
    i64 v;
    if (MODE == 1) {
        i64 best = 1 << 30;
        for (int b = 0; b < 2; b++) {
            int no = options(P, b, Qo); i64 w = -1;
            for (int j = 0; j < no; j++) { i64 x = gval(Qo[j]); if (x > w) w = x; }
            if (w < best) best = w;
        }
        if (best >= (1 << 30)) return TOP8;
        v = best - e; if (v < 0) v = 0;
        if (SHIFT) v -= FN * vv(P);
        if (v < 0) v = 0;          /* cannot happen for the least fixed point (g >= fn v); harmless: raising h keeps the test sound */
    } else {
        i64 best = -1;
        for (int b = 0; b < 2; b++) {
            int no = options(P, b, Qo); i64 w = 1 << 30;
            for (int j = 0; j < no; j++) { i64 x = gval(Qo[j]); if (x < w) w = x; }
            if (w > best) best = w;
        }
        if (best >= (1 << 30)) return TOP8;
        v = best + e; if (v < 0) v = 0;
    }
    return v > CAP ? TOP8 : (int)v;
}

int main(int argc, char **argv) {
    if (argc < 7) { fprintf(stderr, "usage: %s lower|upper q k fn fd prefix [cap [shift]]\n", argv[0]); return 2; }
    MODE = !strcmp(argv[1], "lower") ? 1 : 0;
    Q = atoll(argv[2]); int k = atoi(argv[3]); FN = atoll(argv[4]); FD = atoll(argv[5]);
    const char *prefix = argv[6];
    if (!(Q & 1) || k < 4 || k > 32) { fprintf(stderr, "bad q/k\n"); return 2; }
    N = (i64)1 << k; H = N >> 1; HH = H >> 1; NR = HH;
    { uint64_t x = (uint64_t)Q; for (int i = 0; i < 6; i++) x = x * (2ULL - (uint64_t)Q * x); QINV = (i64)(x & (uint64_t)(N - 1)); }
    CAP = argc > 7 ? atoll(argv[7]) : 254;
    SHIFT = (argc > 8 && MODE == 1) ? atoi(argv[8]) : 0;
    KK = k;
    if (CAP > 254) CAP = 254;
    double t0 = now();
    V = (uint8_t *)calloc((size_t)(NR + 1), 1);
    i64 nw = (NR + 1 + 63) / 64;
    DIRTY = (uint64_t *)malloc(sizeof(uint64_t) * (size_t)nw);
    if (!V || !DIRTY) { fprintf(stderr, "out of memory\n"); return 2; }
    for (i64 w = 0; w < nw; w++) DIRTY[w] = ~0ULL;
    if ((NR + 1) % 64) DIRTY[nw - 1] = (1ULL << ((NR + 1) % 64)) - 1;
    i64 work = 0; int sweeps = 0, any = 1, lost = 0;
    while (any && !lost) {
        any = 0; sweeps++;
        for (i64 w = 0; w < nw && !lost; w++) {
            while (DIRTY[w]) {
                int bit = __builtin_ctzll(DIRTY[w]);
                DIRTY[w] &= DIRTY[w] - 1;
                i64 r = w * 64 + bit;
                work++;
                if (V[r] == TOP8) continue;
                int v = eval(r);
                if (v > V[r]) {
                    V[r] = (uint8_t)v; any = 1;
                    if (v == TOP8 && MODE == 0) { lost = 1; break; }   /* MIN: one TOP pair refutes the test */
                    i64 pr[3];
                    pr[0] = rep((2 * r) & (H - 1));
                    pr[1] = rep((QINV * ((2 * r - 1) & (H - 1))) & (H - 1));
                    pr[2] = rep((QINV * ((2 * r + 1) & (H - 1))) & (H - 1));
                    for (int j = 0; j < 3; j++) if (V[pr[j]] != TOP8) DIRTY[pr[j] >> 6] |= 1ULL << (pr[j] & 63);
                }
            }
        }
    }
    i64 ntop = 0, vmax = 0;
    for (i64 r = 0; r <= NR; r++) { if (V[r] == TOP8) ntop++; else if (V[r] > vmax) vmax = V[r]; }
    printf("  %s F=%lld/%lld cap %lld: work %lld (%d sweeps), TOP reps %lld of %lld, max value %lld, %.1f s\n",
           MODE ? "LOWER" : "UPPER", (long long)FN, (long long)FD, (long long)CAP, (long long)work, sweeps,
           (long long)ntop, (long long)(NR + 1), (long long)vmax, now() - t0);
    if ((MODE == 1 && ntop == NR + 1) || (MODE == 0 && (ntop > 0 || lost))) { printf("%s-FAILED\n", MODE ? "LOWER" : "UPPER"); return 1; }
    const i64 CH = 1 << 22;
    uint8_t *b8 = (uint8_t *)malloc((size_t)CH);
    char fnm[1024]; FILE *f;
    snprintf(fnm, sizeof fnm, "%s.%s", prefix, MODE ? "lo_taub" : "up_sigb"); f = fopen(fnm, "wb");
    for (i64 lo = 0; lo < H; lo += 8 * CH) {
        i64 n = H - lo < 8 * CH ? H - lo : 8 * CH;
        memset(b8, 0, (size_t)((n + 7) / 8));
        for (i64 i = 0; i < n; i++) {
            i64 P = lo + i; int bitv = 0;
            if (MODE == 1) {
                if (V[rep(P)] == TOP8) continue;
                i64 Qo[2], bv[2];
                for (int b = 0; b < 2; b++) {
                    int no = options(P, b, Qo); i64 v = -1;
                    for (int j = 0; j < no; j++) { i64 x = gval(Qo[j]); if (x > v) v = x; }
                    bv[b] = v;
                }
                bitv = bv[1] < bv[0];
            } else {
                i64 x = 2 * P + 1, t = Q * x;
                i64 Qp = ((t + 1) >> 1) & (H - 1), Qm = ((t - 1) >> 1) & (H - 1);
                bitv = gval(Qm) < gval(Qp);
            }
            if (bitv) b8[i >> 3] |= (uint8_t)(1u << (i & 7));
        }
        fwrite(b8, 1, (size_t)((n + 7) / 8), f);
    }
    fclose(f);
    snprintf(fnm, sizeof fnm, "%s.%s", prefix, MODE ? (SHIFT ? "lo_h8" : "lo_g8") : "up_psi8"); f = fopen(fnm, "wb");
    for (i64 lo = 0; lo < H; lo += CH) {
        i64 n = H - lo < CH ? H - lo : CH;
        for (i64 i = 0; i < n; i++) b8[i] = V[rep(lo + i)];
        fwrite(b8, 1, (size_t)n, f);
    }
    fclose(f);
    printf("%s-WRITTEN %lld %lld time=%.2f\n", MODE ? "LOWER" : "UPPER", (long long)FN, (long long)FD, now() - t0);
    return 0;
}
