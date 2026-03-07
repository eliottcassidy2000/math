"""
h21_exhaustive_n8.py — Exhaustive H=21 check at n=8

Checks all 2^28 = 268,435,456 tournaments on 8 vertices.
Uses C with Held-Karp for H(T) computation.

Key optimization: for each tournament, first check if it has a vertex
with no 3-cycle (such vertex is in no cycle, so we can reduce to n=7).
Then check H via Held-Karp only for tournaments where every vertex
is in some 3-cycle.

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

/* Check if vertex v is in any 3-cycle */
static inline int in_any_3cycle(unsigned int T, int v) {
    for (int a = 0; a < N; a++) {
        if (a == v) continue;
        for (int b = a+1; b < N; b++) {
            if (b == v) continue;
            /* Check {v,a,b} for 3-cycle */
            int va = arc(T,v,a), ab = arc(T,a,b), bv = arc(T,b,v);
            if (va && ab && bv) return 1; /* v->a->b->v */
            int vb = arc(T,v,b), ba = arc(T,b,a), av = arc(T,a,v);
            if (vb && ba && av) return 1; /* v->b->a->v */
        }
    }
    return 0;
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
    long long all_in_3cycle_count = 0;
    long long has_free_vertex_count = 0;

    for (unsigned int T = 0; T < TOTAL; T++) {
        /* Check if every vertex is in at least one 3-cycle */
        int all_in = 1;
        for (int v = 0; v < N; v++) {
            if (!in_any_3cycle(T, v)) {
                all_in = 0;
                break;
            }
        }

        if (!all_in) {
            has_free_vertex_count++;
            /* A vertex with no 3-cycle is in no cycle at all.
               By induction (verified at n=7), H(T-v) != 21.
               So H(T) = H(T-v) != 21. Skip. */
            continue;
        }

        all_in_3cycle_count++;

        /* Compute H via Held-Karp */
        long long H = hamiltonian_paths(T);

        if (H == 21) {
            h21_count++;
            printf("H=21 FOUND! T=%u\n", T);
            /* Print adjacency matrix */
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    if (i == j) printf(".");
                    else printf("%d", arc(T,i,j));
                }
                printf("\n");
            }
            fflush(stdout);
        }

        if ((T+1) % 10000000 == 0) {
            fprintf(stderr, "Progress: %u/%u (%.1f%%). free=%lld all_3c=%lld h21=%lld\n",
                    T+1, TOTAL, 100.0*(T+1)/TOTAL,
                    has_free_vertex_count, all_in_3cycle_count, h21_count);
        }
    }

    printf("\n=== EXHAUSTIVE n=8 RESULTS ===\n");
    printf("Total tournaments: %u\n", TOTAL);
    printf("Has free vertex (skipped): %lld\n", has_free_vertex_count);
    printf("All vertices in 3-cycle (checked): %lld\n", all_in_3cycle_count);
    printf("H=21 found: %lld\n", h21_count);

    return 0;
}
"""

def main():
    tmpdir = tempfile.mkdtemp()
    src = os.path.join(tmpdir, "h21_n8.c")
    exe = os.path.join(tmpdir, "h21_n8")

    with open(src, 'w') as f:
        f.write(C_CODE)

    subprocess.run(["gcc", "-O3", "-o", exe, src], check=True)
    print("Compiled. Running EXHAUSTIVE n=8 check (268M tournaments)...")
    print("This will take several minutes...")
    sys.stdout.flush()

    proc = subprocess.Popen([exe], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # Stream output
    import select
    while True:
        reads = [proc.stdout, proc.stderr]
        readable, _, _ = select.select(reads, [], [], 1)
        for stream in readable:
            line = stream.readline()
            if line:
                if stream == proc.stderr:
                    print(line.strip(), file=sys.stderr)
                else:
                    print(line.strip())
        if proc.poll() is not None:
            # Read remaining
            for line in proc.stdout:
                print(line.strip())
            for line in proc.stderr:
                print(line.strip(), file=sys.stderr)
            break

if __name__ == "__main__":
    main()
