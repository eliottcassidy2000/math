/* Audit2 (D, E): exact Hamiltonian path counts, written independently of the session's hk64.c / hk_big.c / ie_ham64.c.
   Modes:
     hcount hk            : read upper-triangle codes ('1' at (i,j), i<j, iff i -> j) on stdin, print H per line
                            (Held-Karp over subsets, uint64; n <= 22)
     hcount iecirc n s1 s2 ...: circulant Cay(Z_n, {s_i}); H by inclusion-exclusion over vertex subsets
                            H = sum_X (-1)^(n-|X|) 1^T A_X^(n-1) 1, one subset per rotation orbit (orbit size = period),
                            __int128 arithmetic (32-bit masks for n <= 31, 64-bit masks for 32 <= n <= 40)
     hcount ie            : read codes on stdin, inclusion-exclusion without symmetry (n <= 26)
     hcount fneg n s1 ... : signed count F(T,-1) = sum over orderings of (-1)^(#forward steps) for a circulant, by
                            Held-Karp with weights S_vw = +-1 over subsets containing vertex 0, times n (wrapping uint64
                            arithmetic is exact mod 2^64; |F| < 2^63 required)
   Output of H as decimal. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef __int128 i128;
static int n;
static int A[32][32];
static int B[64][64];   /* adjacency for the 64-bit-mask mode */

static void print_i128(i128 x) {
    char buf[64]; int k = 0, neg = 0;
    if (x < 0) { neg = 1; x = -x; }
    if (x == 0) buf[k++] = '0';
    while (x > 0) { buf[k++] = '0' + (int)(x % 10); x /= 10; }
    if (neg) putchar('-');
    while (k > 0) putchar(buf[--k]);
}

static int read_code(char *s) {
    int m = (int)strlen(s);
    while (m > 0 && (s[m - 1] == '\n' || s[m - 1] == '\r' || s[m - 1] == ' ')) s[--m] = 0;
    if (m == 0) return 0;
    int nn = 1; while (nn * (nn - 1) / 2 < m) nn++;
    if (nn * (nn - 1) / 2 != m) { fprintf(stderr, "bad code length %d\n", m); exit(1); }
    n = nn; int k = 0;
    memset(A, 0, sizeof A);
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (s[k] == '1') A[i][j] = 1; else A[j][i] = 1; k++; }
    return 1;
}

static uint64_t held_karp(void) {
    size_t M = (size_t)1 << n;
    uint64_t *dp = calloc(M * n, sizeof(uint64_t));
    if (!dp) { fprintf(stderr, "alloc\n"); exit(1); }
    for (int v = 0; v < n; v++) dp[((size_t)1 << v) * n + v] = 1;
    for (size_t m = 1; m < M; m++)
        for (int v = 0; v < n; v++) {
            uint64_t c = dp[m * n + v];
            if (!c) continue;
            for (int w = 0; w < n; w++) if (A[v][w] && !((m >> w) & 1)) dp[(m | ((size_t)1 << w)) * n + w] += c;
        }
    uint64_t H = 0;
    for (int v = 0; v < n; v++) H += dp[(M - 1) * n + v];
    free(dp);
    return H;
}

static i128 walks(uint32_t X) {
    /* number of walks with n vertices inside T[X] */
    i128 u[32], w[32];
    int vs[32], k = 0;
    for (int v = 0; v < n; v++) if ((X >> v) & 1) vs[k++] = v;
    for (int i = 0; i < k; i++) u[i] = 1;
    for (int step = 0; step < n - 1; step++) {
        for (int j = 0; j < k; j++) {
            i128 s = 0; int b = vs[j];
            for (int i = 0; i < k; i++) if (A[vs[i]][b]) s += u[i];
            w[j] = s;
        }
        memcpy(u, w, sizeof(i128) * k);
    }
    i128 t = 0;
    for (int i = 0; i < k; i++) t += u[i];
    return t;
}

static i128 ie_all(void) {
    i128 H = 0;
    uint32_t M = (uint32_t)1 << n;
    for (uint32_t X = 1; X < M; X++) {
        int sz = __builtin_popcount(X);
        i128 W = walks(X);
        if ((n - sz) & 1) H -= W; else H += W;
    }
    return H;
}

static i128 ie_circ(void) {
    i128 H = 0;
    uint32_t full = (n == 32) ? 0xffffffffu : (((uint32_t)1 << n) - 1);
    uint32_t M = (uint32_t)1 << n;
    for (uint32_t X = 1; X < M; X++) {
        /* canonical: X minimal among rotations; period = smallest r with rot^r(X) = X */
        uint32_t y = X; int minimal = 1, period = n;
        for (int r = 1; r < n; r++) {
            y = ((y << 1) | (y >> (n - 1))) & full;
            if (y < X) { minimal = 0; break; }
            if (y == X && period == n) period = r;
        }
        if (!minimal) continue;
        int sz = __builtin_popcount(X);
        i128 W = walks(X) * period;                 /* orbit size = period */
        if ((n - sz) & 1) H -= W; else H += W;
    }
    return H;
}


static i128 walks64(uint64_t X) {
    i128 u[64], w[64];
    int vs[64], k = 0;
    for (int v = 0; v < n; v++) if ((X >> v) & 1) vs[k++] = v;
    for (int i = 0; i < k; i++) u[i] = 1;
    for (int step = 0; step < n - 1; step++) {
        for (int j = 0; j < k; j++) {
            i128 s = 0; int b = vs[j];
            for (int i = 0; i < k; i++) if (B[vs[i]][b]) s += u[i];
            w[j] = s;
        }
        memcpy(u, w, sizeof(i128) * k);
    }
    i128 t = 0;
    for (int i = 0; i < k; i++) t += u[i];
    return t;
}

static i128 ie_circ64(void) {
    /* same as ie_circ with 64-bit masks (n <= 40) */
    i128 H = 0;
    uint64_t full = (((uint64_t)1) << n) - 1;
    uint64_t M = ((uint64_t)1) << n;
    for (uint64_t X = 1; X < M; X++) {
        uint64_t y = X; int minimal = 1, period = n;
        for (int r = 1; r < n; r++) {
            y = ((y << 1) | (y >> (n - 1))) & full;
            if (y < X) { minimal = 0; break; }
            if (y == X && period == n) period = r;
        }
        if (!minimal) continue;
        int sz = __builtin_popcountll(X);
        i128 W = walks64(X) * period;
        if ((n - sz) & 1) H -= W; else H += W;
    }
    return H;
}


static i128 ie_circ64_part(int k, int K) {
    /* ie_circ64 restricted to masks X with X % K == k (to split the work across processes) */
    i128 H = 0;
    uint64_t full = (((uint64_t)1) << n) - 1;
    uint64_t M = ((uint64_t)1) << n;
    uint64_t X0 = (uint64_t)(k == 0 ? K : k);
    for (uint64_t X = X0; X < M; X += (uint64_t)K) {
        uint64_t y = X; int minimal = 1, period = n;
        for (int r = 1; r < n; r++) {
            y = ((y << 1) | (y >> (n - 1))) & full;
            if (y < X) { minimal = 0; break; }
            if (y == X && period == n) period = r;
        }
        if (!minimal) continue;
        int sz = __builtin_popcountll(X);
        i128 W = walks64(X) * period;
        if ((n - sz) & 1) H -= W; else H += W;
    }
    return H;
}

static void set_circ(int argc, char **argv, int first) {
    n = atoi(argv[first]);
    memset(A, 0, sizeof A); memset(B, 0, sizeof B);
    for (int i = first + 1; i < argc; i++) {
        int s = atoi(argv[i]) % n;
        for (int u = 0; u < n; u++) { B[u][(u + s) % n] = 1; if (n <= 32) A[u][(u + s) % n] = 1; }
    }
    for (int u = 0; u < n; u++) for (int v = 0; v < n; v++)
        if (u != v && B[u][v] == B[v][u]) { fprintf(stderr, "not a tournament\n"); exit(1); }
}

static void fneg(void) {
    /* signed paths starting at 0: weights S = A - A^T; masks containing 0 */
    size_t M = (size_t)1 << (n - 1);              /* bits for vertices 1..n-1 */
    uint64_t *dp = calloc(M * n, sizeof(uint64_t));
    if (!dp) { fprintf(stderr, "alloc\n"); exit(1); }
    dp[0 * n + 0] = 1;
    for (size_t m = 0; m < M; m++)
        for (int v = 0; v < n; v++) {
            uint64_t c = dp[m * n + v];
            if (!c) continue;
            for (int w = 1; w < n; w++) {
                if ((m >> (w - 1)) & 1) continue;
                size_t m2 = m | ((size_t)1 << (w - 1));
                if (A[v][w]) dp[m2 * n + w] += c; else dp[m2 * n + w] -= c;
            }
        }
    uint64_t tot = 0;
    for (int v = 0; v < n; v++) tot += dp[(M - 1) * n + v];
    free(dp);
    /* F(T,-1) = sum (-1)^(#forward) = (-1)^(n-1) * sum (-1)^(#backward) = (-1)^(n-1) * n * tot */
    int64_t s = (int64_t)tot;
    i128 F = (i128)s * n;
    if ((n - 1) & 1) F = -F;
    printf("n=%d F(T,-1)=", n); print_i128(F); printf("\n");
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage\n"); return 1; }
    char line[4096];
    if (!strcmp(argv[1], "hk")) {
        while (fgets(line, sizeof line, stdin)) if (read_code(line)) printf("%llu\n", (unsigned long long)held_karp());
    } else if (!strcmp(argv[1], "ie")) {
        while (fgets(line, sizeof line, stdin)) if (read_code(line)) { print_i128(ie_all()); printf("\n"); fflush(stdout); }
    } else if (!strcmp(argv[1], "iecirc")) {
        set_circ(argc, argv, 2); print_i128(n <= 31 ? ie_circ() : ie_circ64()); printf("\n");
    } else if (!strcmp(argv[1], "iecircpart")) {        /* hcount iecircpart k K n s1 s2 ...: partial sum over X = k mod K */
        int k = atoi(argv[2]), K = atoi(argv[3]);
        set_circ(argc, argv, 4); print_i128(ie_circ64_part(k, K)); printf("\n");
    } else if (!strcmp(argv[1], "iecirc64")) {          /* 64-bit path forced (self-test) */
        set_circ(argc, argv, 2); print_i128(ie_circ64()); printf("\n");
    } else if (!strcmp(argv[1], "fneg")) {
        set_circ(argc, argv, 2); fneg();
    }
    return 0;
}
