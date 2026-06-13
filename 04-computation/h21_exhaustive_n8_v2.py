"""
h21_exhaustive_n8_v2.py — Exhaustive H=21 check at n=8, optimized

Strategy:
1. Skip tournaments with source/sink (score 0 or 7) → induction from n=7
2. Skip tournaments with any vertex not in a 3-cycle → vertex-free induction
3. Skip tournaments with t3 > 10 → alpha_1 > 10 → H > 21
4. Only compute H via Held-Karp for the remaining (rare) tournaments

Author: opus-2026-03-07-S41
"""

import os, sys, tempfile, subprocess

C_CODE = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 8
#define NBITS (N*(N-1)/2)
#define TOTAL (1U << NBITS)

static inline int arc(unsigned int T, int i, int j) {
    if (i < j) {
        int pos = i*(2*N-i-1)/2 + (j-i-1);
        return (T >> pos) & 1;
    } else {
        int pos = j*(2*N-j-1)/2 + (i-j-1);
        return 1 - ((T >> pos) & 1);
    }
}

/* Compute scores and check for source/sink */
static inline int has_source_sink(unsigned int T) {
    for (int i = 0; i < N; i++) {
        int s = 0;
        for (int j = 0; j < N; j++) {
            if (i != j && arc(T, i, j)) s++;
        }
        if (s == 0 || s == N-1) return 1;
    }
    return 0;
}

/* Count 3-cycles */
static inline int count_3cycles(unsigned int T) {
    int t3 = 0;
    for (int a = 0; a < N; a++)
        for (int b = a+1; b < N; b++)
            for (int c = b+1; c < N; c++) {
                int ab = arc(T,a,b), bc = arc(T,b,c), ca = arc(T,c,a);
                if ((ab && bc && ca) || ((1-ab) && (1-bc) && (1-ca)))
                    t3++;
            }
    return t3;
}

/* Check if all vertices are in at least one 3-cycle */
static inline int all_in_3cycle(unsigned int T) {
    for (int v = 0; v < N; v++) {
        int found = 0;
        for (int a = 0; a < N && !found; a++) {
            if (a == v) continue;
            for (int b = a+1; b < N && !found; b++) {
                if (b == v) continue;
                int va = arc(T,v,a), ab = arc(T,a,b), bv = arc(T,b,v);
                if (va && ab && bv) { found = 1; break; }
                int vb = arc(T,v,b), ba = arc(T,b,a), av = arc(T,a,v);
                if (vb && ba && av) { found = 1; break; }
            }
        }
        if (!found) return 0;
    }
    return 1;
}

/* Held-Karp for Hamiltonian path count */
long long hamiltonian_paths(unsigned int T) {
    static long long dp[1<<N][N];
    memset(dp, 0, sizeof(dp));
    for (int v = 0; v < N; v++)
        dp[1<<v][v] = 1;
    for (unsigned int mask = 1; mask < (1U<<N); mask++) {
        for (int v = 0; v < N; v++) {
            if (!(mask & (1<<v))) continue;
            if (dp[mask][v] == 0) continue;
            for (int u = 0; u < N; u++) {
                if (mask & (1<<u)) continue;
                if (arc(T, v, u))
                    dp[mask|(1<<u)][u] += dp[mask][v];
            }
        }
    }
    long long total = 0;
    unsigned int full = (1U<<N)-1;
    for (int v = 0; v < N; v++) total += dp[full][v];
    return total;
}

int main() {
    long long h21_count = 0;
    long long skipped_src_sink = 0;
    long long skipped_free_vertex = 0;
    long long skipped_high_t3 = 0;
    long long checked_hk = 0;

    for (unsigned int T = 0; T < TOTAL; T++) {
        /* Filter 1: source/sink */
        if (has_source_sink(T)) {
            skipped_src_sink++;
            continue;
        }

        /* Filter 2: count 3-cycles */
        int t3 = count_3cycles(T);
        if (t3 > 10) {
            skipped_high_t3++;
            continue;
        }

        /* Filter 3: all vertices in 3-cycle */
        if (!all_in_3cycle(T)) {
            skipped_free_vertex++;
            continue;
        }

        /* Compute H via Held-Karp */
        checked_hk++;
        long long H = hamiltonian_paths(T);

        if (H == 21) {
            h21_count++;
            printf("H=21 FOUND! T=%u H=%lld t3=%d\n", T, H, t3);
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    if (i == j) printf(".");
                    else printf("%d", arc(T,i,j));
                }
                printf("\n");
            }
            fflush(stdout);
        }

        if ((T+1) % 50000000 == 0) {
            fprintf(stderr, "Progress: %u/%u (%.1f%%). src_sink=%lld free=%lld high_t3=%lld checked=%lld h21=%lld\n",
                    T+1, TOTAL, 100.0*(T+1)/TOTAL,
                    skipped_src_sink, skipped_free_vertex, skipped_high_t3, checked_hk, h21_count);
        }
    }

    printf("\n=== EXHAUSTIVE n=8 RESULTS ===\n");
    printf("Total tournaments: %u\n", TOTAL);
    printf("Skipped (source/sink): %lld\n", skipped_src_sink);
    printf("Skipped (free vertex, t3<=10): %lld\n", skipped_free_vertex);
    printf("Skipped (t3 > 10): %lld\n", skipped_high_t3);
    printf("Checked via Held-Karp: %lld\n", checked_hk);
    printf("H=21 found: %lld\n", h21_count);

    return 0;
}
"""

def main():
    tmpdir = tempfile.mkdtemp()
    src = os.path.join(tmpdir, "h21_n8_v2.c")
    exe = os.path.join(tmpdir, "h21_n8_v2")

    with open(src, 'w') as f:
        f.write(C_CODE)

    subprocess.run(["gcc", "-O3", "-o", exe, src], check=True)
    print("Compiled. Running EXHAUSTIVE n=8 check with optimizations...")
    sys.stdout.flush()

    proc = subprocess.run([exe], capture_output=True, text=True, timeout=600)
    print(proc.stdout)
    if proc.stderr:
        for l in proc.stderr.strip().split('\n')[-5:]:
            print(l, file=sys.stderr)

if __name__ == "__main__":
    main()
