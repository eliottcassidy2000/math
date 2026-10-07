/* A2 audit: per-class flags for K: descent, maximal-branch, join (c(w) not max in its (a, c mod 3^a) group).
   Counts: descent-unc, descent+join unc (Roosendaal), maximal unc, maximal+join unc, join-only (join but not branch). */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef unsigned __int128 u128;
static uint64_t P3[64];
static int K, S_, A_;
static uint8_t *dflag, *mflag;
static uint64_t *cval; static uint8_t *aval;
static int bw(int b, int E, uint64_t rr) {
    if (E + 1 < S_ - A_) {
        uint64_t m = P3[A_ - b]; uint64_t r2 = (uint64_t)(((u128)rr * 2) % m);
        if ((u128)P3[A_ - b] << (b + E + 1) < ((u128)1 << S_)) return 1;
        if (bw(b, E + 1, r2)) return 1;
    }
    if (b < A_ && rr % 3 == 2) {
        uint64_t m = P3[A_ - b]; uint64_t t = (uint64_t)(((u128)rr * 2 + m - 1) % m); uint64_t r2 = (t / 3) % P3[A_ - b - 1];
        if ((u128)P3[A_ - b - 1] << (b + 1 + E) < ((u128)1 << S_)) return 1;
        if (bw(b + 1, E, r2)) return 1;
    }
    return 0;
}
/* rho = n mod 2^s; c = exact constant: x_s = (3^a n + c)/2^s */
static void fw(int s, int a, uint64_t r, uint64_t rho, uint64_t c, int dc, int mc) {
    if (!dc && s > 0 && (u128)P3[a] < ((u128)1 << s)) dc = 1;
    if (dc) mc = 1;
    if (!mc && s > 0) { S_ = s; A_ = a; if (bw(0, 0, r)) mc = 1; }
    if (s == K) { dflag[rho] = dc; mflag[rho] = mc; cval[rho] = c; aval[rho] = a; return; }
    /* parity of x_s at n = rho */
    u128 x = ((u128)P3[a] * rho + c) >> s;     /* exact since divisible */
    int par = (int)(x & 1);
    for (int beta = 0; beta <= 1; beta++) {
        uint64_t rho2 = (par == beta) ? rho : rho + (1ULL << s);
        if (beta == 0) {
            uint64_t m = P3[a]; uint64_t r2 = (uint64_t)(((u128)r * ((m + 1) / 2)) % m);
            fw(s + 1, a, r2, rho2, c, dc, mc);
        } else {
            uint64_t m = P3[a + 1]; uint64_t r2 = (uint64_t)(((u128)(3 * (u128)r + 1) * ((m + 1) / 2)) % m);
            /* x_{s+1} = (3 x_s + 1)/2 = (3^{a+1} n + 3c + 2^s)/2^{s+1} */
            fw(s + 1, a + 1, r2, rho2, 3 * c + (1ULL << s), dc, mc);
        }
    }
}
typedef struct { uint64_t key, c; uint32_t y; } rec;
static int cmpr(const void *p, const void *q) { const rec *x = p, *y = q; if (x->key != y->key) return x->key < y->key ? -1 : 1; return x->c < y->c ? -1 : x->c > y->c; }
int main(int argc, char **argv) {
    K = atoi(argv[1]);
    P3[0] = 1; for (int i = 1; i < 40; i++) P3[i] = P3[i - 1] * 3;
    uint64_t n = 1ULL << K;
    dflag = calloc(n, 1); mflag = calloc(n, 1); cval = calloc(n, 8); aval = calloc(n, 1);
    fw(0, 0, 0, 0, 0, 0, 0);
    /* verify c against direct computation for a few residues */
    for (uint64_t y = 0; y < n; y += n / 64 + 1) {
        u128 x = y; int a = 0; for (int i = 0; i < K; i++) { if (x & 1) { x = (3 * x + 1) / 2; a++; } else x /= 2; }
        u128 cc = (x << K) - (u128)P3[a] * y;
        if (cc != cval[y] || a != aval[y]) { printf("MISMATCH at %llu\n", (unsigned long long)y); return 1; }
    }
    rec *R = malloc(n * sizeof(rec));
    for (uint64_t y = 0; y < n; y++) { R[y].key = (uint64_t)aval[y] * P3[K] + cval[y] % P3[aval[y]]; R[y].c = cval[y]; R[y].y = (uint32_t)y; }
    qsort(R, n, sizeof(rec), cmpr);
    uint8_t *jflag = calloc(n, 1);
    for (uint64_t i = 0; i < n; i++) { if (i + 1 < n && R[i + 1].key == R[i].key) jflag[R[i].y] = 1; }  /* not the max of its group */
    uint64_t du = 0, dju = 0, mu = 0, mju = 0, jonly = 0;
    for (uint64_t y = 0; y < n; y++) {
        if (!dflag[y]) du++;
        if (!dflag[y] && !jflag[y]) dju++;
        if (!mflag[y]) mu++;
        if (!mflag[y] && !jflag[y]) mju++;
        if (!mflag[y] && jflag[y]) jonly++;
    }
    printf("K=%d descent_unc=%llu descent+joins_unc=%llu maximal_unc=%llu maximal+joins_unc=%llu join_only(not branch)=%llu\n", K,
           (unsigned long long)du, (unsigned long long)dju, (unsigned long long)mu, (unsigned long long)mju, (unsigned long long)jonly);
    /* list maximal-certified but descent+join-uncertified classes at small K */
    if (K <= 12) {
        printf("classes certified by branches but not by descent/joins at K=%d:", K);
        for (uint64_t y = 0; y < n; y++) if (mflag[y] && !dflag[y] && !jflag[y]) printf(" %llu", (unsigned long long)y);
        printf("\n");
    }
    return 0;
}
