/* Exhaustive (FINITE-EXACT) certified merge shares for the T-clock relation chain.
   State (j, N): x = 3^j y + c, N = c * 3^max(0,-j).  Coins p = y's T-parities (i.i.d. fair for Haar y).
   For each depth s <= K, count coin words of length s whose chain first reaches (0,0) at step s.
   share(K') = sum_{s <= K'} count(s) 2^{-s} = Haar measure of {y : y and x merge within K' T-steps}.
   Usage: ./tclock_certified K j0 N0 [forced prefix as a string of 0/1]
   Exact integer arithmetic in __int128 (overflow-checked).                                         */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef __int128 i128;
static i128 P3[80];
static unsigned long long cnt[80];
static int K;
static int overflow = 0;
static const i128 LIM = ((i128)1) << 120;

static inline int step(int j, i128 N, int p, int *jo, i128 *No) {
    int e = (int)(N & 1);
    if (N < 0) e = (int)((-N) & 1);
    i128 r;
    if (p == 0 && e == 0) { *jo = j; r = N / 2; }
    else if (p == 0 && e == 1) {
        if (j >= 0) { *jo = j + 1; r = (3 * N + 1) / 2; }
        else { int d = -j; *jo = j + 1; r = (N + P3[d - 1]) / 2; }
    } else if (p == 1 && e == 0) {
        if (j >= 0) { *jo = j; r = (3 * N + 1 - P3[j]) / 2; }
        else { int d = -j; *jo = j; r = (3 * N + P3[d] - 1) / 2; }
    } else {
        if (j >= 1) { *jo = j - 1; r = (N - P3[j - 1]) / 2; }
        else { *jo = j - 1; r = (3 * N - 1) / 2; }
    }
    if (r > LIM || r < -LIM) overflow = 1;
    *No = r;
    return 0;
}

static void dfs(int s, int j, i128 N) {
    if (j == 0 && N == 0) { cnt[s]++; return; }
    if (s == K) return;
    int j1; i128 N1;
    step(j, N, 0, &j1, &N1); dfs(s + 1, j1, N1);
    step(j, N, 1, &j1, &N1); dfs(s + 1, j1, N1);
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s K j0 N0 [prefix]\n", argv[0]); return 1; }
    K = atoi(argv[1]);
    int j = atoi(argv[2]);
    i128 N = (i128)atoll(argv[3]);
    const char *pre = argc > 4 ? argv[4] : "";
    P3[0] = 1; for (int i = 1; i < 80; i++) P3[i] = 3 * P3[i - 1];
    int s = 0;
    for (const char *q = pre; *q; q++) {
        if (j == 0 && N == 0) break;
        int j1; i128 N1; step(j, N, *q - '0', &j1, &N1); j = j1; N = N1; s++;
    }
    if (j == 0 && N == 0) { printf("merged inside the forced prefix\n"); return 0; }
    /* depth is counted from the start (forced steps included); coins after the prefix are fair */
    int s0 = s;
    K = K;  /* total depth */
    /* run DFS with fair coins from depth s0 to K */
    {
        /* reuse dfs with offset: shift counts */
        dfs(s0, j, N);
    }
    if (overflow) { printf("OVERFLOW\n"); return 2; }
    double share = 0.0;
    printf("# start (j0,N0)=(%s,%s) prefix '%s'  K=%d\n# s  count(s)  share(<=s) = sum count 2^-(s - s0)\n",
           argv[2], argv[3], pre, K);
    for (int t = s0; t <= K; t++) {
        share += (double)cnt[t] * __builtin_ldexp(1.0, -(t - s0));
        printf("%d %llu %.10f\n", t, cnt[t], share);
    }
    return 0;
}
