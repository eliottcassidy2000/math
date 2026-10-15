/* Audit B (2026-10-14): independent Carmichael census for THM-4651 (4).
 *
 * Method (different from the session's vectorised Korselt test over a smallest-prime-factor sieve):
 *   every Carmichael number is an odd base-2 Fermat pseudoprime, so
 *   1. byte sieve of Eratosthenes for primality up to N + 2;
 *   2. for every odd composite n < N, test 2^(n-1) == 1 (mod n)  (64-bit arithmetic, n < 2^32);
 *   3. factor each base-2 pseudoprime by trial division and apply Korselt's criterion
 *      (squarefree, p - 1 | n - 1 for every p | n).
 * For each Carmichael number print its factorisation, omega, the Giuga score
 * G = #{p | n : p | n/p - 1}, and whether n - 2, n + 2 are prime.
 *
 * Usage: carmichael_psp2 N
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef uint64_t u64;
static u64 powmod(u64 b, u64 e, u64 m) {
    u64 r = 1 % m; b %= m;
    while (e) { if (e & 1) r = (r * b) % m; b = (b * b) % m; e >>= 1; }
    return r;
}

int main(int argc, char **argv) {
    u64 N = argc > 1 ? strtoull(argv[1], 0, 10) : 100000000ULL;
    u64 M = N + 3;
    unsigned char *comp = calloc(M + 1, 1);
    comp[0] = comp[1] = 1;
    for (u64 i = 2; i * i <= M; i++) if (!comp[i]) for (u64 j = i * i; j <= M; j += i) comp[j] = 1;
    u64 smallp[20000]; int nsp = 0;
    for (u64 i = 2; i * i <= M + 100000 && nsp < 20000; i++) if (i <= M && !comp[i]) smallp[nsp++] = i;
    u64 npsp = 0, ncarm = 0, ntw = 0;
    long dist[12][12] = {{0}};
    for (u64 n = 3; n < N; n += 2) {
        if (!comp[n]) continue;
        if (powmod(2, n - 1, n) != 1) continue;
        npsp++;
        /* factor */
        u64 m = n, f[40]; int nf = 0, sqfree = 1;
        for (int i = 0; i < nsp && smallp[i] * smallp[i] <= m; i++) {
            u64 p = smallp[i];
            if (m % p == 0) { f[nf++] = p; m /= p; if (m % p == 0) { sqfree = 0; break; } }
        }
        if (!sqfree) continue;
        if (m > 1) f[nf++] = m;
        int kor = 1;
        for (int i = 0; i < nf; i++) if ((n - 1) % (f[i] - 1) != 0) { kor = 0; break; }
        if (!kor || nf < 2) continue;
        ncarm++;
        int G = 0;
        char gl[256]; int gpos = 0; gl[0] = 0;
        for (int i = 0; i < nf; i++) if (((n / f[i]) - 1) % f[i] == 0) { G++; gpos += snprintf(gl + gpos, sizeof gl - gpos, "%llu,", (unsigned long long)f[i]); }
        dist[G][nf]++;
        int lo = !comp[n - 2], hi = !comp[n + 2];
        if (lo || hi) ntw++;
        printf("C %llu omega %d G %d giuga_at [%s] factors", (unsigned long long)n, nf, G, gl);
        for (int i = 0; i < nf; i++) printf(" %llu", (unsigned long long)f[i]);
        printf(" n-2prime %d n+2prime %d\n", lo, hi);
    }
    printf("N = %llu: base-2 pseudoprimes (odd composite) %llu; Carmichael %llu; with n-2 or n+2 prime %llu\n",
           (unsigned long long)N, (unsigned long long)npsp, (unsigned long long)ncarm, (unsigned long long)ntw);
    printf("distribution (G, omega):");
    for (int g = 0; g < 12; g++) for (int w = 0; w < 12; w++) if (dist[g][w]) printf(" (%d,%d): %ld", g, w, dist[g][w]);
    printf("\n");
    return 0;
}
