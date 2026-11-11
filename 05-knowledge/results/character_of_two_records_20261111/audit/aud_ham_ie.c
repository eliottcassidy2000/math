/* AUDIT program (independent, second method): number of Hamiltonian paths modulo 2^64 by inclusion-exclusion
 *     H = sum over vertex sets S of (-1)^(n - |S|) * #{walks on n vertices that stay inside S},
 * no subset dynamic programme, no shared code with aud_ham.c.  n <= 26 (time 2^n * n^3 / 4).
 * usage: aud_ham_ie < input     (input: n, then n rows of 0/1)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

int main(void) {
    int n, i, j;
    char row[64];
    uint32_t in[32];
    if (scanf("%d", &n) != 1 || n < 1 || n > 26) return 1;
    for (i = 0; i < n; i++) in[i] = 0;
    for (i = 0; i < n; i++) {
        if (scanf("%63s", row) != 1) return 1;
        for (j = 0; j < n; j++) if (row[j] == '1') in[j] |= 1u << i;   /* i -> j */
    }
    uint64_t total = 0;
    uint32_t S, full = (n == 32) ? 0xffffffffu : ((1u << n) - 1);
    for (S = 1; S <= full; S++) {
        int members[32], k = 0;
        uint32_t t;
        for (t = S; t; t &= t - 1) members[k++] = __builtin_ctz(t);
        uint64_t c[32], d[32];
        for (i = 0; i < k; i++) c[members[i]] = 1;
        int step;
        for (step = 1; step < n; step++) {
            for (i = 0; i < k; i++) {
                int w = members[i];
                uint64_t s = 0;
                for (t = in[w] & S; t; t &= t - 1) s += c[__builtin_ctz(t)];
                d[w] = s;
            }
            for (i = 0; i < k; i++) c[members[i]] = d[members[i]];
        }
        uint64_t walks = 0;
        for (i = 0; i < k; i++) walks += c[members[i]];
        if ((n - k) & 1) total -= walks; else total += walks;
        if (S == full) break;
    }
    printf("%llu\n", (unsigned long long)total);
    return 0;
}
