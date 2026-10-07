/* qsurv.c -- survival of "no merge with any partner" for several THM-4581 pair chains driven by the SAME bits.
   lane core, golden session 2026-10-07.

   Driver v (Haar 2-adic) with fair parity bits beta_t (Terras: fresh fair coins).  Partner i is u_i = 3^{k_i} v + e_i,
   tracked exactly as (k, N) with e = N / 3^max(0,-k)  (THM-4581 table; same as coal_check.step).
   q(T) = P(no partner chain absorbed at (0,0) by Terras time T).  Partners whose (k,N) coincide have merged with each
   other (they are the same orbit from then on) and are deduplicated.

   modes:
     trans R          : partners v+1, ..., v+R      (chains (0, r)); this is q_R of the brief
     mers  DMAX       : S19 Mersenne any-lag family: driver p = 1 mod 16 (forced bits 1,0,1,0 then fair),
                        partners q_D = 3^{-D} p + (3^{-D} - 1), odd D <= DMAX  (chains (-D, 1 - 3^D))
     single k0 N0     : one chain from (k0, N0)
   usage: qsurv MODE PARAM PATHS TMAX SEED [extra]
   output: T  survivors  q(T)  sqrt(T) q(T)  on a log grid; plus censoring counts (bignum capacity). */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#define CAP 96               /* limbs: 6144 bits; |k| up to ~3800 */
typedef unsigned __int128 u128;
typedef struct { int n; uint64_t d[CAP]; } big;   /* two's complement, n active limbs (sign = top bit of d[n-1]) */

static big P3[4200];         /* 3^j, j < 4200 (unsigned magnitudes, stored as nonnegative two's complement) */
static int P3MAX = 0;

static inline int sgn_top(const big *x) { return (int)(x->d[x->n - 1] >> 63); }
static inline void normalize(big *x) {
    while (x->n > 1) {
        uint64_t top = x->d[x->n - 1], nxt = x->d[x->n - 2];
        if ((top == 0 && !(nxt >> 63)) || (top == ~0ULL && (nxt >> 63))) x->n--; else break;
    }
}
static inline void extend(big *x, int m) {           /* sign-extend to m limbs */
    uint64_t f = sgn_top(x) ? ~0ULL : 0;
    while (x->n < m) x->d[x->n++] = f;
}
static inline void set_small(big *x, int64_t v) { x->n = 1; x->d[0] = (uint64_t)v; }
static inline int is_zero(const big *x) { return x->n == 1 && x->d[0] == 0; }
static inline int is_odd(const big *x) { return (int)(x->d[0] & 1); }
static int overflowed = 0;
static long nwit = 0; static int witprinted = 0; static char wbits[4096];
static inline void half(big *x) {                    /* exact arithmetic shift right by 1 */
    for (int i = 0; i < x->n - 1; i++) x->d[i] = (x->d[i] >> 1) | (x->d[i + 1] << 63);
    x->d[x->n - 1] = (uint64_t)((int64_t)x->d[x->n - 1] >> 1);
    normalize(x);
}
static inline void mul3add(big *x, int64_t c) {      /* x = 3x + c */
    if (x->n + 1 > CAP) { overflowed = 1; return; }
    extend(x, x->n + 1);
    u128 carry = 0;
    for (int i = 0; i < x->n; i++) { u128 v = (u128)x->d[i] * 3 + carry; x->d[i] = (uint64_t)v; carry = v >> 64; }
    /* add signed c */
    uint64_t cc = (uint64_t)c, fill = c < 0 ? ~0ULL : 0; u128 cr = 0;
    for (int i = 0; i < x->n; i++) {
        u128 v = (u128)x->d[i] + (i == 0 ? cc : fill) + cr; x->d[i] = (uint64_t)v; cr = v >> 64;
    }
    normalize(x);
}
static inline void addpow3(big *x, int j, int sign) { /* x += sign * 3^j */
    const big *p = &P3[j];
    int m = (x->n > p->n ? x->n : p->n) + 1;
    if (m > CAP) { overflowed = 1; return; }
    extend(x, m);
    if (sign > 0) {
        u128 cr = 0;
        for (int i = 0; i < m; i++) { u128 v = (u128)x->d[i] + (i < p->n ? p->d[i] : 0) + cr; x->d[i] = (uint64_t)v; cr = v >> 64; }
    } else {
        u128 br = 0; /* x - p */
        for (int i = 0; i < m; i++) {
            uint64_t pi = (i < p->n ? p->d[i] : 0);
            u128 v = (u128)x->d[i] - pi - br; x->d[i] = (uint64_t)v; br = (v >> 64) ? 1 : 0;
        }
    }
    normalize(x);
}
static int big_eq(const big *a, const big *b) {
    if (a->n != b->n) return 0;
    return memcmp(a->d, b->d, sizeof(uint64_t) * a->n) == 0;
}

typedef struct { int k; big N; int alive; } chain;

/* one Terras step; returns 1 if absorbed */
static inline int cstep(chain *c, int beta) {
    int k = c->k; big *N = &c->N;
    if (k >= 0) {
        if (!is_odd(N)) {
            if (!beta) half(N);
            else { mul3add(N, 1); addpow3(N, k, -1); half(N); }
        } else {
            if (!beta) { mul3add(N, 1); half(N); c->k = k + 1; }
            else if (k >= 1) { addpow3(N, k - 1, -1); half(N); c->k = k - 1; }
            else { mul3add(N, -1); half(N); c->k = -1; }
        }
    } else {
        int a = -k;
        if (a + 1 >= P3MAX) { overflowed = 1; return 0; }
        if (!is_odd(N)) {
            if (!beta) half(N);
            else { mul3add(N, -1); addpow3(N, a, +1); half(N); }
        } else {
            if (!beta) { addpow3(N, a - 1, +1); half(N); c->k = k + 1; }
            else { mul3add(N, -1); half(N); c->k = k - 1; }
        }
    }
    if (c->k >= P3MAX - 1) overflowed = 1;
    return c->k == 0 && is_zero(N);
}

/* xoshiro256** */
static uint64_t s[4];
static inline uint64_t rotl(const uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static uint64_t next64(void) {
    const uint64_t result = rotl(s[1] * 5, 7) * 9; const uint64_t t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45); return result;
}
static uint64_t splitmix(uint64_t *x) { uint64_t z = (*x += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; return z ^ (z >> 31); }

static uint64_t hashc(const chain *c) {
    uint64_t h = 1469598103934665603ULL ^ (uint64_t)(int64_t)c->k;
    for (int i = 0; i < c->N.n; i++) { h ^= c->N.d[i]; h *= 1099511628211ULL; h ^= h >> 29; }
    return h ^ (uint64_t)c->N.n * 0x9e3779b97f4a7c15ULL;
}
typedef struct { uint64_t h; int i; } hi;
static int cmphi(const void *a, const void *b) { uint64_t x = ((const hi *)a)->h, y = ((const hi *)b)->h; return x < y ? -1 : x > y; }

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: qsurv MODE PARAM PATHS TMAX SEED [N0]\n"); return 1; }
    const char *mode = argv[1]; int param = atoi(argv[2]); long paths = atol(argv[3]); long TMAX = atol(argv[4]);
    uint64_t seed = strtoull(argv[5], 0, 10);
    /* powers of 3 */
    P3MAX = 4200;
    set_small(&P3[0], 1);
    for (int j = 1; j < P3MAX; j++) { P3[j] = P3[j - 1]; mul3add(&P3[j], 0); if (overflowed) { P3MAX = j; overflowed = 0; break; } }
    int nch0 = 0; chain *init = NULL; int forced = 0; int forced_bits[4] = {1, 0, 1, 0};
    if (!strcmp(mode, "trans")) {
        nch0 = param; init = calloc(nch0, sizeof(chain));
        for (int r = 1; r <= param; r++) { init[r - 1].k = 0; set_small(&init[r - 1].N, r); }
    } else if (!strcmp(mode, "mers")) {
        for (int D = 1; D <= param; D += 2) nch0++;
        init = calloc(nch0, sizeof(chain)); int i = 0;
        for (int D = 1; D <= param; D += 2, i++) {           /* N = 1 - 3^D */
            init[i].k = -D; set_small(&init[i].N, 1); addpow3(&init[i].N, D, -1);
        }
        forced = 4;
    } else if (!strcmp(mode, "single")) {
        nch0 = 1; init = calloc(1, sizeof(chain)); init[0].k = param; set_small(&init[0].N, atoll(argv[6]));
    } else { fprintf(stderr, "bad mode\n"); return 1; }
    /* log grid */
    int NG = 0; long grid[400];
    for (double x = 1.0; ; x *= pow(10.0, 0.1)) { long g = (long)floor(x + 0.5); if (g > TMAX) break; if (NG == 0 || g != grid[NG - 1]) grid[NG++] = g; }
    if (grid[NG - 1] != TMAX) grid[NG++] = TMAX;
    long *surv = calloc(NG, sizeof(long)); long censored = 0; double total_steps = 0;
    chain *ch = malloc(sizeof(chain) * nch0); hi *H = malloc(sizeof(hi) * nch0);
    uint64_t sm = seed * 0x2545F4914F6CDD1DULL + 12345;
    for (long path = 0; path < paths; path++) {
        for (int j = 0; j < 4; j++) s[j] = splitmix(&sm);
        memcpy(ch, init, sizeof(chain) * nch0);
        int m = nch0; long t; long tau = -1; overflowed = 0;
        uint64_t bits = 0; int nb = 0;
        for (t = 0; t < TMAX; t++) {
            int beta;
            if (t < forced) beta = forced_bits[t];
            else { if (nb == 0) { bits = next64(); nb = 64; } beta = bits & 1; bits >>= 1; nb--; }
            if (t < 4096) wbits[t] = (char)beta;
            int absorbed = 0;
            for (int i = 0; i < m; i++) if (cstep(&ch[i], beta)) { absorbed = 1; }
            if (overflowed) break;
            if (absorbed) { tau = t + 1; break; }
            if (m == 1 && nch0 > 1) break;   /* all partners coalesced: witness, stop path */
            if (m > 1 && ((t & 7) == 7)) {                 /* lazy dedup */
                for (int i = 0; i < m; i++) { H[i].h = hashc(&ch[i]); H[i].i = i; }
                qsort(H, m, sizeof(hi), cmphi);
                char *kill = calloc(m, 1); int nk = 0;
                for (int i = 1; i < m; i++) if (H[i].h == H[i - 1].h) {
                    int a = H[i - 1].i, b = H[i].i;
                    if (ch[a].k == ch[b].k && big_eq(&ch[a].N, &ch[b].N)) { kill[b] = 1; nk++; H[i].i = a; }
                }
                if (nk) { int w = 0; for (int i = 0; i < m; i++) if (!kill[i]) ch[w++] = ch[i]; m = w; }
                free(kill);
            }
        }
        if (m == 1 && tau < 0 && !overflowed && nch0 > 1) {
            nwit++;
            if (!witprinted) { witprinted = 1;
                printf("# first witness: path %ld, all %d partners coalesced by Terras time %ld, driver unmerged; first %ld driver bits: ", path, nch0, t, t + 1);
                for (long j = 0; j <= t && j < 4096; j++) putchar('0' + wbits[j]); printf("\n"); }
        }
        total_steps += (double)t;
        if (overflowed) { censored++; /* count as surviving only up to t (conservative: drop from grid beyond t) */
            for (int g = 0; g < NG; g++) if (grid[g] <= t) surv[g]++;
            continue; }
        for (int g = 0; g < NG; g++) if (tau < 0 || tau > grid[g]) surv[g]++;
    }
    printf("# witnesses (all partners coalesced before any merge with driver): %ld of %ld paths\n", nwit, paths);
    printf("# qsurv mode=%s param=%d paths=%ld TMAX=%ld seed=%llu chains=%d censored=%ld mean_steps=%.1f\n",
           mode, param, paths, TMAX, (unsigned long long)seed, nch0, censored, total_steps / paths);
    printf("# T survivors q sqrtT_q\n");
    for (int g = 0; g < NG; g++) {
        double q = (double)surv[g] / paths;
        printf("%ld %ld %.6e %.4f\n", grid[g], surv[g], q, sqrt((double)grid[g]) * q);
    }
    return 0;
}
