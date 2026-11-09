/* ham_ie.c - audit P1 (independent): Hamiltonian path counts by inclusion-exclusion over walk counts.
 *
 * Algorithm (different from the session's subset DP): for a digraph on V (|V| = n),
 *     H(T)     = sum_{X subset V} (-1)^(n - |X|)     W_{n-1}(X),
 *     H(T - v) = sum_{X subset V \ v} (-1)^(n-1-|X|) W_{n-2}(X),
 * where W_k(X) is the number of walks with k arcs inside X. All arithmetic is mod 2^64 (unsigned wraparound);
 * the true counts are < 2^64, so the residues are the exact values. Memory is O(n).
 * Input lines: name n out_0 ... out_{n-1}. Output: name n H and H(T - v) for every v.
 * Compile: cc -O2 -o ham_ie ham_ie.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;

int main(void) {
    char name[64]; int n;
    while (scanf("%63s %d", name, &n) == 2) {
        uint32_t out[32], in[32];
        memset(in, 0, sizeof in);
        for (int i = 0; i < n; i++) { unsigned long long m; if (scanf("%llu", &m) != 1) return 1; out[i] = (uint32_t)m; }
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (out[i] >> j & 1) in[j] |= 1u << i;
        uint32_t full = (1u << n) - 1;
        u64 H = 0, Hd[32]; memset(Hd, 0, sizeof Hd);
        u64 w[32], nw[32];
        int mem[32];
        for (uint32_t X = 1; X <= full; X++) {
            int k = 0; uint32_t s = X;
            while (s) { mem[k++] = __builtin_ctz(s); s &= s - 1; }
            for (int i = 0; i < k; i++) w[mem[i]] = 1;
            u64 W2 = 0, W1 = 0;
            if (n - 2 == 0) { for (int i = 0; i < k; i++) W2 += w[mem[i]]; }
            for (int step = 1; step <= n - 1; step++) {
                for (int i = 0; i < k; i++) {
                    int y = mem[i]; uint32_t t = in[y] & X; u64 acc = 0;
                    while (t) { int x = __builtin_ctz(t); t &= t - 1; acc += w[x]; }
                    nw[y] = acc;
                }
                for (int i = 0; i < k; i++) w[mem[i]] = nw[mem[i]];
                if (step == n - 2) for (int i = 0; i < k; i++) W2 += w[mem[i]];
                if (step == n - 1) for (int i = 0; i < k; i++) W1 += w[mem[i]];
            }
            if ((n - k) & 1) H -= W1; else H += W1;
            /* deletion of v not in X: sign (-1)^(n-1-k) */
            uint32_t c = ~X & full;
            while (c) { int v = __builtin_ctz(c); c &= c - 1; if ((n - 1 - k) & 1) Hd[v] -= W2; else Hd[v] += W2; }
        }
        printf("%s n=%d H=%llu del:", name, n, H);
        for (int v = 0; v < n; v++) printf(" %llu", Hd[v]);
        printf("\n");
        fflush(stdout);
    }
    return 0;
}
