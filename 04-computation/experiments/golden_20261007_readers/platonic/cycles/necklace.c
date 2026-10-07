/* necklace.c -- integer cycles of the Terras map T(x) = x/2 | (3x+1)/2 on Z via Lyndon words.
 * For each L, enumerate binary Lyndon words a_1..a_L (FKM algorithm) and maintain
 *   d <- 3d + 2^(t-1) at a 1 in position t (1-indexed), j = number of ones,
 * so T^L(x) = (3^j x + d)/2^L for the 2-adic x whose parity vector begins a_1..a_L.
 * The word is the parity word (from the element x) of an integer cycle iff (2^L - 3^j) | d,
 * and then x = d/(2^L - 3^j).  Each primitive cycle <-> exactly one Lyndon word.
 * Pruning: for L >= 2 an integer cycle needs j >= L/2 (product identity), so a subtree is cut when
 *   2*(j + #remaining positions) < L.
 * Leaves are counted per j and compared with the exact Lyndon count
 *   (1/L) sum_{e | gcd(L,j)} mu(e) C(L/e, j/e).
 * d is kept as unsigned __int128 (exact; d < 3^L).  Parallel over a frontier at depth T0.
 * usage: necklace L1 L2 NTHREADS
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <stdatomic.h>
typedef unsigned __int128 u128;
static int L, T0;
static u128 mod_j[200]; static int pos_j[200];
typedef struct { int t, p, j; u128 d; unsigned char a[200]; } state_t;
static state_t *front; static int nfront, capfront;
static _Atomic int nextitem;
static pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;
static uint64_t leaves_j[200];
static void u128str(u128 v, char *buf) { char t[64]; int n = 0; if (!v) { strcpy(buf, "0"); return; }
    while (v) { t[n++] = '0' + (int)(v % 10); v /= 10; } for (int i = 0; i < n; i++) buf[i] = t[n - 1 - i]; buf[n] = 0; }
static void found(const unsigned char *a, u128 d, int j) {
    u128 x = d / mod_j[j]; char buf[64]; u128str(x, buf);
    pthread_mutex_lock(&mu);
    printf("FOUND L=%d k=%d x=%s%s word=", L, j, pos_j[j] ? "" : "-", buf);
    for (int i = 1; i <= L; i++) putchar('0' + a[i]);
    putchar('\n'); fflush(stdout);
    pthread_mutex_unlock(&mu);
}
typedef struct { unsigned char a[200]; uint64_t cnt[200]; } ctx_t;
static void gen(ctx_t *c, int t, int p, u128 d, int j) {
    if (t > L) {
        if (p == L) { c->cnt[j]++; if (L <= 40 ? ((uint64_t)d % (uint64_t)mod_j[j] == 0) : (d % mod_j[j] == 0)) found(c->a, d, j); }
        return;
    }
    if (L >= 2 && 2 * (j + (L - t + 1)) < L) return;
    int b = c->a[t - p];
    c->a[t] = (unsigned char)b;
    gen(c, t + 1, p, b ? 3 * d + ((u128)1 << (t - 1)) : d, j + b);
    if (b == 0) { c->a[t] = 1; gen(c, t + 1, t, 3 * d + ((u128)1 << (t - 1)), j + 1); }
}
static void genfront(unsigned char *a, int t, int p, u128 d, int j) {
    if (t > L || t == T0) {
        if (nfront == capfront) { capfront = capfront ? 2 * capfront : 1024; front = realloc(front, capfront * sizeof(state_t)); }
        state_t *s = &front[nfront++]; s->t = t; s->p = p; s->j = j; s->d = d; memcpy(s->a, a, 200);
        return;
    }
    if (L >= 2 && 2 * (j + (L - t + 1)) < L) return;
    int b = a[t - p]; a[t] = (unsigned char)b;
    genfront(a, t + 1, p, b ? 3 * d + ((u128)1 << (t - 1)) : d, j + b);
    if (b == 0) { a[t] = 1; genfront(a, t + 1, t, 3 * d + ((u128)1 << (t - 1)), j + 1); }
}
static void *worker(void *arg) {
    ctx_t *c = calloc(1, sizeof(ctx_t));
    for (;;) {
        int i = atomic_fetch_add(&nextitem, 1); if (i >= nfront) break;
        state_t *s = &front[i]; memcpy(c->a, s->a, 200);
        gen(c, s->t, s->p, s->d, s->j);
    }
    pthread_mutex_lock(&mu); for (int j = 0; j <= L; j++) leaves_j[j] += c->cnt[j]; pthread_mutex_unlock(&mu);
    free(c); return 0;
}
static int mobius(int n) { int r = 1; for (int q = 2; q * q <= n; q++) if (n % q == 0) { n /= q; if (n % q == 0) return 0; r = -r; } if (n > 1) r = -r; return r; }
static long double binom(int n, int k) { long double r = 1; for (int i = 1; i <= k; i++) r = r * (n - k + i) / i; return r; }
static int gcd(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }
int main(int argc, char **argv) {
    int L1 = atoi(argv[1]), L2 = atoi(argv[2]), nt = atoi(argv[3]);
    for (L = L1; L <= L2; L++) {
        u128 two = (u128)1 << L, p3 = 1;
        for (int j = 0; j <= L; j++) { if (two > p3) { mod_j[j] = two - p3; pos_j[j] = 1; } else { mod_j[j] = p3 - two; pos_j[j] = 0; } p3 *= 3; }
        memset(leaves_j, 0, sizeof leaves_j);
        T0 = L > 22 ? 18 : 2; nfront = 0; atomic_store(&nextitem, 0);
        unsigned char a[200]; memset(a, 0, 200);
        genfront(a, 1, 1, 0, 0);
        pthread_t th[64]; for (int i = 0; i < nt; i++) pthread_create(&th[i], 0, worker, 0);
        for (int i = 0; i < nt; i++) pthread_join(th[i], 0);
        /* check leaf counts against the Lyndon formula for j >= L/2 (all j when L == 1) */
        int ok = 1; long double tot = 0;
        for (int j = 0; j <= L; j++) {
            if (L >= 2 && 2 * j < L) continue;
            long double s = 0; int g = gcd(L, j); if (j == 0) g = L;
            for (int e = 1; e <= g; e++) if (g % e == 0) s += mobius(e) * binom(L / e, j / e);
            s /= L; tot += s;
            if ((long double)leaves_j[j] != s) { ok = 0; printf("  MISMATCH L=%d j=%d leaves=%llu formula=%.0Lf\n", L, j, (unsigned long long)leaves_j[j], s); }
        }
        printf("L=%d lyndon_leaves(j>=L/2)=%.0Lf count_check=%s frontier=%d\n", L, tot, ok ? "OK" : "FAIL", nfront); fflush(stdout);
    }
    return 0;
}
