/* AUDIT a2 (independent auditor, 2026-10-21).  s(m) = # strong labelled tournaments, reduced mod M, computed as
 * minus the Hurwitz inverse of T(m) = 2^C(m,2):  b_0 = 1, b_m = -sum_{k=1..m} C(m,k) T(k) b_{m-k}, s(m) = -b_m.
 * Binomials by an in-place Pascal row (no division, so any modulus M works).  All arithmetic in unsigned __int128.
 * Modes:
 *   period PMAX K      : every odd prime p <= PMAX: check s(m+p-1) = eps s(m) and D(m+p-1) = eps D(m) (mod p)
 *                        for 1 <= m <= K*p, eps := 2^((p-1)/2) mod p (cross-checked against p mod 8).
 *   direct NMIN NMAX   : every composite n in [NMIN,NMAX]: print "n pass s(n)%n c(n)%n".
 *   single n [M]       : s(n) mod M (default M = n) and c(n) mod M.
 *   wief PMAX          : odd primes p <= PMAX with D(p) = s(p) - 2^((p-1)/2) = 0 (mod p^2).
 *   seq N M            : print s(m) mod M for m = 1..N.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
typedef unsigned __int128 u128;
static u64 mulm(u64 a, u64 b, u64 m) { return (u64)((u128)a * b % m); }
static u64 powm(u64 b, u64 e, u64 m) { u64 r = 1 % m; b %= m; while (e) { if (e & 1) r = mulm(r, b, m); b = mulm(b, b, m); e >>= 1; } return r; }
static int isprime(u64 n) { if (n < 2) return 0; for (u64 d = 2; d * d <= n; d++) if (n % d == 0) return 0; return 1; }
/* fills s[1..N] with s(m) mod M */
static void strong_mod(int N, u64 M, u64 *s) {
    u64 *T = malloc((N + 1) * sizeof(u64)), *bb = malloc((N + 1) * sizeof(u64)), *row = calloc(N + 2, sizeof(u64));
    T[0] = 1 % M;
    for (int k = 1; k <= N; k++) T[k] = mulm(T[k - 1], powm(2, k - 1, M), M);
    bb[0] = 1 % M; row[0] = 1 % M;
    int small = (M < (1ULL << 31));
    for (int m = 1; m <= N; m++) {
        for (int k = m; k >= 1; k--) { row[k] = row[k] + row[k - 1]; if (row[k] >= M) row[k] -= M; } /* row = C(m,.) */
        u64 acc = 0;
        if (small) {
            for (int k = 1; k <= m; k++) { acc += (row[k] * T[k] % M) * bb[m - k] % M; if (acc >= M) acc -= M; }
        } else {
            for (int k = 1; k <= m; k++) { acc += mulm(mulm(row[k], T[k], M), bb[m - k], M); if (acc >= M) acc -= M; }
        }
        bb[m] = (M - acc) % M;
        s[m] = acc % M; /* s(m) = -b_m = acc */
    }
    free(T); free(bb); free(row);
}
static u64 circ_mod(u64 n, u64 M) { if (n == 1) return 1 % M; return (n & 1) ? powm(2, (n - 1) / 2, M) : 0; }
int main(int argc, char **argv) {
    if (argc < 2) return 1;
    if (!strcmp(argv[1], "period")) {
        int PMAX = atoi(argv[2]), K = atoi(argv[3]); long fails = 0, checks = 0;
        for (int p = 3; p <= PMAX; p += 2) {
            if (!isprime(p)) continue;
            int N = (K + 1) * p + 2; u64 *s = malloc((N + 1) * sizeof(u64)); strong_mod(N, p, s);
            u64 e = powm(2, (p - 1) / 2, p); int eps = (e == 1) ? 1 : -1;
            int eps8 = (p % 8 == 1 || p % 8 == 7) ? 1 : -1;
            if (e != 1 && e != (u64)p - 1) printf("EULER?? p=%d\n", p);
            if (eps != eps8) printf("EPS MISMATCH p=%d\n", p);
            for (int m = 1; m <= K * p; m++) {
                u64 lhs = s[m + p - 1], rhs = (eps == 1) ? s[m] : (p - s[m]) % p;
                u64 Dl = (s[m + p - 1] + p - circ_mod(m + p - 1, p)) % p, Dm = (s[m] + p - circ_mod(m, p)) % p;
                u64 Dr = (eps == 1) ? Dm : (p - Dm) % p;
                checks++;
                if (lhs != rhs || Dl != Dr) { fails++; printf("FAIL p=%d m=%d\n", p, m); }
            }
            /* the m = 0 case is NOT claimed: report whether s(p-1) = 0 mod p */
            printf("p=%d eps=%+d  s(p-1) mod p = %llu\n", p, eps, s[p - 1]);
            free(s);
        }
        printf("period checks %ld failures %ld\n", checks, fails);
    } else if (!strcmp(argv[1], "direct")) {
        int a = atoi(argv[2]), bnd = atoi(argv[3]);
        for (int n = a; n <= bnd; n++) {
            if (n < 4 || isprime(n)) continue;
            u64 *s = malloc((n + 1) * sizeof(u64)); strong_mod(n, n, s);
            u64 cc = circ_mod(n, n);
            printf("%d %d %llu %llu\n", n, s[n] == cc, s[n], cc); fflush(stdout); free(s);
        }
    } else if (!strcmp(argv[1], "single")) {
        int n = atoi(argv[2]); u64 M = (argc > 3) ? strtoull(argv[3], 0, 10) : (u64)n;
        u64 *s = malloc((n + 1) * sizeof(u64)); strong_mod(n, M, s);
        printf("n=%d M=%llu s(n) mod M = %llu  c(n) mod M = %llu  pass=%d\n", n, M, s[n], circ_mod(n, M), s[n] == circ_mod(n, M));
        free(s);
    } else if (!strcmp(argv[1], "wief")) {
        int PMAX = atoi(argv[2]); int cnt = 0;
        printf("odd primes p <= %d with D(p) = 0 mod p^2:", PMAX);
        for (int p = 3; p <= PMAX; p += 2) {
            if (!isprime(p)) continue; cnt++;
            u64 M = (u64)p * p; u64 *s = malloc((p + 1) * sizeof(u64)); strong_mod(p, M, s);
            if (s[p] == powm(2, (p - 1) / 2, M)) { printf(" %d", p); fflush(stdout); }
            free(s);
        }
        printf("\n(%d odd primes tested)\n", cnt);
    } else if (!strcmp(argv[1], "seq")) {
        int N = atoi(argv[2]); u64 M = strtoull(argv[3], 0, 10);
        u64 *s = malloc((N + 1) * sizeof(u64)); strong_mod(N, M, s);
        for (int m = 1; m <= N; m++) printf("%d %llu\n", m, s[m]);
    }
    return 0;
}
