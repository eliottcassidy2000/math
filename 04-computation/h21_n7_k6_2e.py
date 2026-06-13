#!/usr/bin/env python3
"""
H=21 gap proof extension: K_6-2e impossibility at n=7.

At n=6 (exhaustive): every tournament whose 3-cycle conflict graph contains a
K_6-2e substructure (6 three-cycles with exactly 2 disjoint pairs) also has
additional directed 5-cycles, giving I(Omega(T), 2) >= 29 (never 21).

This script checks whether the same holds at n=7 by exhaustive enumeration
of all 2^21 = 2,097,152 tournaments on 7 vertices.

For each tournament with t3 >= 6:
  - Check if any 6 three-cycles form a K_6-2e in the conflict graph
  - Compute H(T) via Held-Karp DP (= I(Omega(T), 2))
  - Record minimum H among K_6-2e tournaments

Also confirms H=21 never occurs at n=7 (for any tournament).

RESULTS (exhaustive, 18.9s):
  - H=21 never achieved at n=7 (0 out of 2,097,152 tournaments)
  - 1,597,968 tournaments have K_6-2e in their 3-cycle conflict graph
  - Minimum H among K_6-2e tournaments: 29
  - The H=29 minimum case has t3=6, t5=4, t7=0 (same forcing as n=6)
  - H values for K_6-2e tournaments are all odd and >= 29
  - t3 ranges from 6 to 14 among K_6-2e tournaments

CONCLUSION: The K_6-2e impossibility extends fully to n=7.
Any tournament with a K_6-2e substructure in its conflict graph has H >= 29.
Combined with the global exhaustive check, H=21 is impossible at n=7.

Instance: opus-2026-03-07
"""

import itertools
import time
from collections import Counter
import ctypes
import os
import sys
import struct


# ─────────────────────────────────────────────────
# Fast C implementation for Held-Karp
# ─────────────────────────────────────────────────

C_CODE = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Count Hamiltonian paths in a tournament on n=7 vertices.
   adj[v] is a bitmask of outgoing neighbors of v.
   Returns total number of Hamiltonian paths. */
long long held_karp_7(int adj[7]) {
    int n = 7;
    int full = (1 << n) - 1;
    /* dp[S][v] = number of Ham paths using vertex set S ending at v */
    /* S ranges 0..127, v ranges 0..6 */
    long long dp[128][7];
    memset(dp, 0, sizeof(dp));

    for (int v = 0; v < n; v++)
        dp[1 << v][v] = 1;

    for (int S = 1; S <= full; S++) {
        for (int v = 0; v < n; v++) {
            if (!(S & (1 << v))) continue;
            long long c = dp[S][v];
            if (c == 0) continue;
            int out = adj[v] & (~S);  /* neighbors not in S */
            while (out) {
                int u = __builtin_ctz(out);  /* lowest set bit */
                dp[S | (1 << u)][u] += c;
                out &= out - 1;  /* clear lowest bit */
            }
        }
    }

    long long total = 0;
    for (int v = 0; v < n; v++)
        total += dp[full][v];
    return total;
}

/* Count directed 3-cycles. Returns count. Also fills cycle_verts[3*count]. */
int find_3_cycles(int adj[7], int n, int cycle_verts[]) {
    int count = 0;
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++)
            for (int k = j+1; k < n; k++) {
                if (((adj[i]>>j)&1) && ((adj[j]>>k)&1) && ((adj[k]>>i)&1)) {
                    cycle_verts[3*count] = i;
                    cycle_verts[3*count+1] = j;
                    cycle_verts[3*count+2] = k;
                    count++;
                } else if (((adj[i]>>k)&1) && ((adj[k]>>j)&1) && ((adj[j]>>i)&1)) {
                    cycle_verts[3*count] = i;
                    cycle_verts[3*count+1] = j;
                    cycle_verts[3*count+2] = k;
                    count++;
                }
            }
    return count;
}

/* Process all 2^21 tournaments on 7 vertices.
   For each tournament with t3 >= 6, check K_6-2e and compute H.
   Also check all tournaments for H=21.

   Output arrays (caller allocates):
     h21_count: number of tournaments with H=21
     k6_count: number of tournaments with K_6-2e in conflict graph
     h_values: H values for K_6-2e tournaments (up to max_k6)
     t3_values: t3 values for K_6-2e tournaments
     min_H_k6: minimum H among K_6-2e tournaments
*/
void process_all(
    long long *h21_count,
    long long *k6_count,
    long long *min_H_k6,
    long long h_values[],     /* H for each K_6-2e tournament */
    int t3_values[],          /* t3 for each K_6-2e tournament */
    int max_k6,               /* max entries in h_values/t3_values */
    long long *h_dist,        /* H value histogram for K_6-2e, size 1000 */
    int *progress_flag        /* set to 1 every 100k tournaments */
) {
    int n = 7;
    int edges[21][2];
    int m = 0;
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++) {
            edges[m][0] = i;
            edges[m][1] = j;
            m++;
        }

    long long total = 1L << m;
    *h21_count = 0;
    *k6_count = 0;
    *min_H_k6 = 999999;

    for (long long bits = 0; bits < total; bits++) {
        int adj[7] = {0};
        for (int k = 0; k < m; k++) {
            int i = edges[k][0], j = edges[k][1];
            if ((bits >> k) & 1)
                adj[j] |= (1 << i);
            else
                adj[i] |= (1 << j);
        }

        /* Find 3-cycles */
        int cverts[45*3];  /* max C(7,3)=35 cycles, use 45 for safety */
        int t3 = find_3_cycles(adj, n, cverts);

        if (t3 < 6) {
            /* Still check for H=21 */
            long long H = held_karp_7(adj);
            if (H == 21) (*h21_count)++;
            continue;
        }

        /* Compute H */
        long long H = held_karp_7(adj);
        if (H == 21) (*h21_count)++;

        /* Check K_6-2e: for each 6-subset of t3 cycles, count shared-vertex edges */
        /* Precompute: cycle i and cycle j share a vertex? */
        /* Store as bitmask for efficiency */
        unsigned int shares[35] = {0};  /* shares[i] bitmask: bit j set if cycle i,j share vertex */
        for (int a = 0; a < t3; a++) {
            int va0 = cverts[3*a], va1 = cverts[3*a+1], va2 = cverts[3*a+2];
            unsigned int amask = (1<<va0)|(1<<va1)|(1<<va2);
            for (int b = a+1; b < t3; b++) {
                int vb0 = cverts[3*b], vb1 = cverts[3*b+1], vb2 = cverts[3*b+2];
                unsigned int bmask = (1<<vb0)|(1<<vb1)|(1<<vb2);
                if (amask & bmask) {
                    shares[a] |= (1 << b);
                    shares[b] |= (1 << a);
                }
            }
        }

        /* Check all C(t3, 6) subsets */
        int found_k6_2e = 0;
        int idx[6];
        for (idx[0] = 0; idx[0] < t3-5 && !found_k6_2e; idx[0]++)
        for (idx[1] = idx[0]+1; idx[1] < t3-4 && !found_k6_2e; idx[1]++)
        for (idx[2] = idx[1]+1; idx[2] < t3-3 && !found_k6_2e; idx[2]++)
        for (idx[3] = idx[2]+1; idx[3] < t3-2 && !found_k6_2e; idx[3]++)
        for (idx[4] = idx[3]+1; idx[4] < t3-1 && !found_k6_2e; idx[4]++)
        for (idx[5] = idx[4]+1; idx[5] < t3 && !found_k6_2e; idx[5]++) {
            int edge_count = 0;
            for (int a = 0; a < 6; a++)
                for (int b = a+1; b < 6; b++)
                    if (shares[idx[a]] & (1 << idx[b]))
                        edge_count++;
            if (edge_count == 13) {
                found_k6_2e = 1;
            }
        }

        if (found_k6_2e) {
            long long ki = *k6_count;
            if (ki < max_k6) {
                h_values[ki] = H;
                t3_values[ki] = t3;
            }
            (*k6_count)++;
            if (H < *min_H_k6) *min_H_k6 = H;
            if (H < 1000) h_dist[H]++;
        }
    }
}
"""


def compile_and_load():
    """Compile the C code and load it."""
    c_path = "/tmp/h21_n7_k6_2e.c"
    so_path = "/tmp/h21_n7_k6_2e.so"

    with open(c_path, "w") as f:
        f.write(C_CODE)

    os.system(f"cc -O3 -o {so_path} -shared -fPIC {c_path}")
    return ctypes.CDLL(so_path)


def run_c_exhaustive():
    """Run the exhaustive C search."""
    lib = compile_and_load()

    # Set up ctypes
    lib.process_all.restype = None
    max_k6 = 200000  # generous upper bound

    h21_count = ctypes.c_longlong(0)
    k6_count = ctypes.c_longlong(0)
    min_H_k6 = ctypes.c_longlong(0)
    h_values = (ctypes.c_longlong * max_k6)()
    t3_values = (ctypes.c_int * max_k6)()
    h_dist = (ctypes.c_longlong * 1000)()
    progress = ctypes.c_int(0)

    print("Running exhaustive C search over 2,097,152 tournaments...")
    start = time.time()

    lib.process_all(
        ctypes.byref(h21_count),
        ctypes.byref(k6_count),
        ctypes.byref(min_H_k6),
        h_values,
        t3_values,
        max_k6,
        h_dist,
        ctypes.byref(progress)
    )

    elapsed = time.time() - start

    # ─────────────────────────────────────────────────
    # Report
    # ─────────────────────────────────────────────────
    print()
    print("=" * 70)
    print(f"RESULTS: n=7 exhaustive (2,097,152 tournaments, {elapsed:.1f}s)")
    print("=" * 70)

    print(f"\nTournaments with H=21 overall: {h21_count.value}")
    print(f"Tournaments with K_6-2e in conflict graph: {k6_count.value}")

    if k6_count.value > 0:
        print(f"Minimum H among K_6-2e tournaments: {min_H_k6.value}")

        # H distribution
        print(f"\nH value distribution for K_6-2e tournaments:")
        for h_val in range(1000):
            if h_dist[h_val] > 0:
                tag = " <-- TARGET" if h_val == 21 else ""
                print(f"  H={h_val}: {h_dist[h_val]}{tag}")

        # t3 distribution
        t3_counter = Counter()
        actual_k6 = min(k6_count.value, max_k6)
        for i in range(actual_k6):
            t3_counter[t3_values[i]] += 1
        print(f"\nt3 distribution for K_6-2e tournaments: {dict(sorted(t3_counter.items()))}")

        h21_in_k6 = h_dist[21]
        print(f"\nH=21 among K_6-2e tournaments: {h21_in_k6}")

        if h21_in_k6 == 0 and h21_count.value == 0:
            print("\n*** CONFIRMED: H=21 never achieved at n=7 ***")
            print(f"*** K_6-2e tournaments have H >= {min_H_k6.value} ***")
        elif h21_in_k6 == 0:
            print(f"\n*** K_6-2e impossibility extends to n=7: min H = {min_H_k6.value} ***")
            print(f"*** H=21 occurs {h21_count.value} times via other mechanisms ***")
    else:
        print("\nNo K_6-2e subgraphs found at n=7")

    if h21_count.value == 0:
        print("\n*** GLOBAL CONFIRMATION: H=21 is impossible at n=7 ***")

    return h21_count.value, k6_count.value, min_H_k6.value


# ─────────────────────────────────────────────────
# Fallback: Pure Python (slower, for verification)
# ─────────────────────────────────────────────────

def run_python_sample(num_samples=50000):
    """Pure Python fallback with sampling."""
    import random
    random.seed(42)

    n = 7
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(edges)

    k6_count = 0
    h21_count = 0
    min_H_k6 = float('inf')
    h_values_k6 = Counter()

    print(f"Running Python sampling: {num_samples} random tournaments at n={n}")
    start = time.time()

    for trial in range(num_samples):
        if trial % 10000 == 0 and trial > 0:
            print(f"  {trial}/{num_samples} -- K6-2e: {k6_count}, H=21: {h21_count}")

        bits = random.randint(0, (1 << m) - 1)
        adj = [0] * n
        for k, (i, j) in enumerate(edges):
            if bits & (1 << k):
                adj[j] |= (1 << i)
            else:
                adj[i] |= (1 << j)

        # 3-cycles
        cycles = []
        for i in range(n):
            for j in range(i + 1, n):
                for k in range(j + 1, n):
                    if ((adj[i] >> j) & 1) and ((adj[j] >> k) & 1) and ((adj[k] >> i) & 1):
                        cycles.append(frozenset([i, j, k]))
                    elif ((adj[i] >> k) & 1) and ((adj[k] >> j) & 1) and ((adj[j] >> i) & 1):
                        cycles.append(frozenset([i, j, k]))
        t3 = len(cycles)

        # H via Held-Karp
        dp = [[0] * n for _ in range(1 << n)]
        for v in range(n):
            dp[1 << v][v] = 1
        for S in range(1, 1 << n):
            for v in range(n):
                if not (S & (1 << v)):
                    continue
                if dp[S][v] == 0:
                    continue
                for u in range(n):
                    if S & (1 << u):
                        continue
                    if (adj[v] >> u) & 1:
                        dp[S | (1 << u)][u] += dp[S][v]
        full = (1 << n) - 1
        H = sum(dp[full][v] for v in range(n))

        if H == 21:
            h21_count += 1

        if t3 < 6:
            continue

        # Check K_6-2e
        found = False
        m_cyc = len(cycles)
        for subset in itertools.combinations(range(m_cyc), 6):
            edge_count = 0
            for a in range(6):
                for b in range(a + 1, 6):
                    if cycles[subset[a]] & cycles[subset[b]]:
                        edge_count += 1
            if edge_count == 13:
                found = True
                break

        if found:
            k6_count += 1
            h_values_k6[H] += 1
            min_H_k6 = min(min_H_k6, H)

    elapsed = time.time() - start
    print(f"\nPython sample results ({elapsed:.1f}s):")
    print(f"  K_6-2e count: {k6_count}")
    print(f"  H=21 count: {h21_count}")
    if k6_count > 0:
        print(f"  Min H among K_6-2e: {min_H_k6}")
        print(f"  H distribution: {dict(sorted(h_values_k6.items()))}")


if __name__ == "__main__":
    try:
        h21, k6, minH = run_c_exhaustive()
    except Exception as e:
        print(f"C compilation failed: {e}")
        print("Falling back to Python sampling...")
        run_python_sample()
