/*
 * procgen_seven2_20260926_rhomax.c -- exact rho_max(sigma) of ONE sign strategy of q n +- 1
 * (lane "seven2", session collatz-procgen-20260922, 2026-09-26).
 *
 * Setting (THM-4486).  Level k, odd q, N = 2^k, H = 2^(k-1).  sigma assigns a sign to every odd residue
 * mod N; T_sigma(x) = x/2 (x even), (q x + sigma(x))/2 (x odd).  G_sigma has nodes Z/N and edges from x to
 * both lifts of the pair T_sigma(x) mod H.  rho_max(sigma) = largest odd density of a cycle of G_sigma.
 *
 * PAIR GRAPH.  A cycle of G_sigma is a sequence of pairs P_i in Z/H and lift bits b_i with
 * P_(i+1) = T_sigma(P_i + b_i H) mod H; its odd density is the fraction of odd P_i.  So rho_max(sigma) is
 * the maximum mean cycle of the graph on Z/H with out-edges P -> succ_b(P) = T_sigma(P + b H) mod H
 * (b = 0, 1) and weight [P odd].
 *
 * METHOD (Dinkelbach iteration with exact integers).  For a threshold F = fn/fd put e(P) = fd - fn
 * (P odd), -fn (P even), and iterate f(P) <- max(0, e(P) + max_b f(succ_b(P))) from f = 0 (Gauss-Seidel
 * sweeps over a dirty bitset), recording at every strict increase the maximizing lift pi(P).
 *   Invariant: f(P) <= e(P) + f(succ_pi(P)(P)) for every P with f(P) > 0.
 *   (a) If a sweep changes nothing, f is a potential: f(succ_b(P)) + e(P) <= f(P) for all P, b; summing
 *       around any cycle gives e-weight <= 0, i.e. rho_max(sigma) <= F  (upper certificate; written
 *       to disk and re-checked by procgen_seven_20260926_verify.c).
 *   (b) Every cycle of the pointer graph P -> succ_pi(P)(P) (on {f > 0}) has e-weight > 0, i.e. odd
 *       density > F: take the last pointer set on the cycle, at P_j; just before that update
 *       f(P_j) < e(P_j) + f(P_(j+1)) and the other edges satisfy the invariant, so the telescoping sum of
 *       e(P_i) + f(P_(i+1)) - f(P_i) >= 0 around the cycle is > 0.  Such a cycle is an explicit cycle of
 *       G_sigma (nodes P_i + pi(P_i) H), printed, and re-checked exactly by the caller.
 *   After each sweep the pointer graph is searched for a cycle; if one is found, F becomes its density
 *   (strictly larger) and the iteration restarts.  The final F is rho_max(sigma) exactly: a witness
 *   cycle of density F and a potential proving <= F.
 *
 * usage: procgen_seven2_20260926_rhomax q k sigfile [F0n F0d [certprefix]]
 *   sigfile: H bits, LSB first; bit i = 1 iff the sign of the odd node 2i+1 is '-'.
 *   F0 = F0n/F0d: starting threshold (default 0/1); if the iteration at F0 converges at once the
 *   program reports "RHOMAX_LE F0n F0d" (rho_max <= F0, no witness).
 * prints "RHOMAX a p" (reduced), "WITNESS len odd" and "CYCLE x_0 x_1 ... " (nodes of G_sigma), and,
 * if certprefix is given, writes certprefix.up_sig (uint8[H], 1 = '-') and certprefix.up_psi (int32[H])
 * in the full format of procgen_seven_20260926_verify.c (an upper-only certificate at F = a/p).
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

typedef int64_t i64;
static i64 Q, N, H, QINV;
static int K;
static uint8_t *SIG;        /* bits */
static int32_t *F32;        /* potential */
static uint8_t *PTR;        /* pointer bits (lift) */
static uint64_t *DIRTY;
static uint8_t *COL;

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static inline int sigbit(i64 x) { i64 i = (x - 1) >> 1; return (SIG[i >> 3] >> (i & 7)) & 1; }
static inline i64 succ(i64 P, int b) {
    i64 x = P + (b ? H : 0);
    if (!(x & 1)) return (x >> 1) & (H - 1);
    i64 s = sigbit(x) ? -1 : 1;
    return ((Q * x + s) >> 1) & (H - 1);
}
static inline int getptr(i64 P) { return (PTR[P >> 3] >> (P & 7)) & 1; }
static inline void setptr(i64 P, int b) { if (b) PTR[P >> 3] |= (uint8_t)(1u << (P & 7)); else PTR[P >> 3] &= (uint8_t)~(1u << (P & 7)); }
static i64 gcdll(i64 a, i64 b) { while (b) { i64 t = a % b; a = b; b = t; } return a < 0 ? -a : a; }

/* one run at threshold fn/fd; returns 0 = converged (certificate in F32), 1 = pointer cycle found
   (returned in cyc[], length *clen), 2 = overflow/failure */
static int run(i64 fn, i64 fd, i64 **cycp, i64 *clen, i64 maxsweeps) {
    i64 nw = (H + 63) / 64;
    for (i64 P = 0; P < H; P++) F32[P] = 0;
    memset(PTR, 0, (size_t)((H + 7) / 8));
    for (i64 w = 0; w < nw; w++) DIRTY[w] = ~0ULL;
    if (H % 64) DIRTY[nw - 1] = (1ULL << (H % 64)) - 1;
    i64 eo = fd - fn, ee = -fn;
    for (i64 sweep = 0; sweep < maxsweeps; sweep++) {
        int any = 0;
        for (i64 w = 0; w < nw; w++) {
            uint64_t bits = DIRTY[w];          /* snapshot: every pair is evaluated at most once per sweep */
            DIRTY[w] = 0;
            while (bits) {
                int bit = __builtin_ctzll(bits);
                bits &= bits - 1;
                i64 P = w * 64 + bit;
                i64 e = (P & 1) ? eo : ee;
                i64 v0 = F32[succ(P, 0)], v1 = F32[succ(P, 1)];
                int b = v1 > v0;
                i64 v = e + (b ? v1 : v0);
                if (v > F32[P]) {
                    if (v > (1LL << 30)) return 2;
                    F32[P] = (int32_t)v; setptr(P, b); any = 1;
                    i64 pr[3];
                    pr[0] = (2 * P) & (H - 1);
                    pr[1] = (QINV * ((2 * P - 1) & (H - 1))) & (H - 1);
                    pr[2] = (QINV * ((2 * P + 1) & (H - 1))) & (H - 1);
                    for (int j = 0; j < 3; j++) DIRTY[pr[j] >> 6] |= 1ULL << (pr[j] & 63);
                }
            }
        }
        if (!any) return 0;
        /* pointer-cycle search on {f > 0}: colors 0 = unvisited, 1 = on current walk, 2 = done */
        memset(COL, 0, (size_t)H);
        for (i64 s = 0; s < H; s++) {
            if (COL[s] || F32[s] <= 0) continue;
            i64 u = s;
            while (F32[u] > 0 && COL[u] == 0) { COL[u] = 1; u = succ(u, getptr(u)); }
            if (F32[u] > 0 && COL[u] == 1) {            /* cycle through u */
                i64 len = 0, v = u;
                do { len++; v = succ(v, getptr(v)); } while (v != u);
                i64 *cyc = (i64 *)malloc(sizeof(i64) * len);
                v = u;
                for (i64 i = 0; i < len; i++) { cyc[i] = v + (getptr(v) ? H : 0); v = succ(v, getptr(v)); }
                *cycp = cyc; *clen = len;
                return 1;
            }
            u = s;                                      /* mark the walk done */
            while (F32[u] > 0 && COL[u] == 1) { COL[u] = 2; u = succ(u, getptr(u)); }
        }
    }
    return 2;
}

/* BAD mode (exploratory search objective, not a proof step): least fixed point of
   f(P) = max(0, e(P) + max_b f(succ_b(P))) at F = fn/fd with values capped at cap (TOP above); prints the number of
   TOP pairs (pairs from which a path of weight > cap starts, in particular every pair that reaches a cycle of density
   > F).  Sweeps without the pointer machinery. */
static i64 bad_count(i64 fn, i64 fd, i64 cap) {
    i64 nw = (H + 63) / 64;
    for (i64 P = 0; P < H; P++) F32[P] = 0;
    for (i64 w = 0; w < nw; w++) DIRTY[w] = ~0ULL;
    if (H % 64) DIRTY[nw - 1] = (1ULL << (H % 64)) - 1;
    i64 eo = fd - fn, ee = -fn;
    const int32_t TOPV = 0x7fffffff;
    int any = 1;
    while (any) {
        any = 0;
        for (i64 w = 0; w < nw; w++) {
            uint64_t bits = DIRTY[w]; DIRTY[w] = 0;
            while (bits) {
                int bit = __builtin_ctzll(bits); bits &= bits - 1;
                i64 P = w * 64 + bit;
                if (F32[P] == TOPV) continue;
                i64 e = (P & 1) ? eo : ee;
                i64 v0 = F32[succ(P, 0)], v1 = F32[succ(P, 1)];
                i64 m = v0 > v1 ? v0 : v1;
                i64 v = (m == TOPV) ? TOPV : e + m;
                if (v > cap) v = TOPV;
                if (v > F32[P]) {
                    F32[P] = (int32_t)v; any = 1;
                    i64 pr[3];
                    pr[0] = (2 * P) & (H - 1);
                    pr[1] = (QINV * ((2 * P - 1) & (H - 1))) & (H - 1);
                    pr[2] = (QINV * ((2 * P + 1) & (H - 1))) & (H - 1);
                    for (int j = 0; j < 3; j++) DIRTY[pr[j] >> 6] |= 1ULL << (pr[j] & 63);
                }
            }
        }
    }
    i64 c = 0;
    for (i64 P = 0; P < H; P++) if (F32[P] == TOPV) c++;
    return c;
}

int main(int argc, char **argv) {
    if (argc >= 2 && !strcmp(argv[1], "bad")) {
        if (argc < 8) { fprintf(stderr, "usage: %s bad q k sigfile fn fd cap\n", argv[0]); return 2; }
        Q = atoll(argv[2]); K = atoi(argv[3]); N = 1LL << K; H = N >> 1;
        { uint64_t x = (uint64_t)Q; for (int i = 0; i < 6; i++) x = x * (2ULL - (uint64_t)Q * x); QINV = (i64)(x & (uint64_t)(N - 1)); }
        SIG = (uint8_t *)malloc((size_t)((H + 7) / 8));
        FILE *f = fopen(argv[4], "rb");
        if (!f || fread(SIG, 1, (size_t)((H + 7) / 8), f) != (size_t)((H + 7) / 8)) { fprintf(stderr, "cannot read sigfile\n"); return 2; }
        fclose(f);
        F32 = (int32_t *)malloc(sizeof(int32_t) * (size_t)H);
        DIRTY = (uint64_t *)malloc(sizeof(uint64_t) * (size_t)((H + 63) / 64));
        printf("BAD %lld\n", (long long)bad_count(atoll(argv[5]), atoll(argv[6]), atoll(argv[7])));
        return 0;
    }
    if (argc < 4) { fprintf(stderr, "usage: %s q k sigfile [F0n F0d [certprefix]]\n", argv[0]); return 2; }
    Q = atoll(argv[1]); K = atoi(argv[2]);
    if (!(Q & 1) || K < 3 || K > 30) { fprintf(stderr, "bad q/k\n"); return 2; }
    N = 1LL << K; H = N >> 1;
    { uint64_t x = (uint64_t)Q; for (int i = 0; i < 6; i++) x = x * (2ULL - (uint64_t)Q * x); QINV = (i64)(x & (uint64_t)(N - 1)); }
    SIG = (uint8_t *)malloc((size_t)((H + 7) / 8));
    FILE *f = fopen(argv[3], "rb");
    if (!f || fread(SIG, 1, (size_t)((H + 7) / 8), f) != (size_t)((H + 7) / 8)) { fprintf(stderr, "cannot read sigfile\n"); return 2; }
    fclose(f);
    i64 fn = argc > 5 ? atoll(argv[4]) : 0, fd = argc > 5 ? atoll(argv[5]) : 1;
    const char *cert = argc > 6 ? argv[6] : 0;
    F32 = (int32_t *)malloc(sizeof(int32_t) * (size_t)H);
    PTR = (uint8_t *)malloc((size_t)((H + 7) / 8));
    DIRTY = (uint64_t *)malloc(sizeof(uint64_t) * (size_t)((H + 63) / 64));
    COL = (uint8_t *)malloc((size_t)H);
    if (!F32 || !PTR || !DIRTY || !COL) { fprintf(stderr, "out of memory\n"); return 2; }
    double t0 = now();
    i64 wlen = 0, wodd = 0; i64 *wit = 0;
    int iters = 0;
    while (1) {
        i64 *cyc = 0, clen = 0;
        int r = run(fn, fd, &cyc, &clen, 1000000);
        iters++;
        if (r == 2) { printf("FAILED at F=%lld/%lld\n", (long long)fn, (long long)fd); return 1; }
        if (r == 0) break;
        i64 odd = 0;
        for (i64 i = 0; i < clen; i++) odd += cyc[i] & 1;
        /* the pointer cycle has density > F (proved); check it */
        if (!(odd * fd > fn * clen)) { printf("INTERNAL: pointer cycle not denser than F\n"); return 1; }
        i64 g = gcdll(odd, clen);
        fn = odd / g; fd = clen / g;
        if (wit) free(wit);
        wit = cyc; wlen = clen; wodd = odd;
        fprintf(stderr, "  iteration %d: witness cycle of length %lld, density %lld/%lld (%.1f s)\n", iters, (long long)clen,
                (long long)fn, (long long)fd, now() - t0);
    }
    if (!wit) { printf("RHOMAX_LE %lld %lld\n", (long long)fn, (long long)fd); }
    else {
        printf("RHOMAX %lld %lld\n", (long long)fn, (long long)fd);
        printf("WITNESS %lld %lld\n", (long long)wlen, (long long)wodd);
        printf("CYCLE");
        for (i64 i = 0; i < wlen && i < 100000; i++) printf(" %lld", (long long)wit[i]);
        printf("\n");
    }
    i64 fmax = 0;
    for (i64 P = 0; P < H; P++) if (F32[P] > fmax) fmax = F32[P];
    printf("STATS k=%d q=%lld iterations=%d max_psi=%lld time=%.2f\n", K, (long long)Q, iters, (long long)fmax, now() - t0);
    if (cert) {
        char fnm[1024];
        snprintf(fnm, sizeof fnm, "%s.up_sig", cert); f = fopen(fnm, "wb");
        for (i64 i = 0; i < H; i++) { uint8_t b = (SIG[i >> 3] >> (i & 7)) & 1; fwrite(&b, 1, 1, f); }
        fclose(f);
        snprintf(fnm, sizeof fnm, "%s.up_psi", cert); f = fopen(fnm, "wb");
        fwrite(F32, sizeof(int32_t), (size_t)H, f);
        fclose(f);
    }
    return 0;
}
