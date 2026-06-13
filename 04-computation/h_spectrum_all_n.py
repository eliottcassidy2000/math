#!/usr/bin/env python3
"""
H-spectrum gaps at n=3,4,5,6,7 (exhaustive).

Compute all achievable H values and identify permanent gaps.

Instance: opus-2026-03-07-S40
"""

import ctypes
import os
import time

C_CODE = r"""
#include <string.h>
#include <stdlib.h>

/* Generic Held-Karp for n <= 10 */
long long held_karp(int *adj, int n) {
    int full = (1 << n) - 1;
    int states = 1 << n;
    long long *dp = (long long *)calloc(states * n, sizeof(long long));

    for (int v = 0; v < n; v++)
        dp[(1 << v) * n + v] = 1;

    for (int S = 1; S <= full; S++) {
        for (int v = 0; v < n; v++) {
            if (!(S & (1 << v))) continue;
            long long c = dp[S * n + v];
            if (c == 0) continue;
            for (int u = 0; u < n; u++) {
                if (S & (1 << u)) continue;
                if (adj[v * n + u])
                    dp[(S | (1 << u)) * n + u] += c;
            }
        }
    }

    long long total = 0;
    for (int v = 0; v < n; v++)
        total += dp[full * n + v];
    free(dp);
    return total;
}

/* Compute H-spectrum for given n */
void h_spectrum_n(int n, long long *hist, int hist_size) {
    int m = n * (n - 1) / 2;
    int edges[100][2];
    int ei = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) {
            edges[ei][0] = i;
            edges[ei][1] = j;
            ei++;
        }

    long long total = 1L << m;
    for (long long bits = 0; bits < total; bits++) {
        int adj[100] = {0};
        for (int k = 0; k < m; k++) {
            int i = edges[k][0], j = edges[k][1];
            if ((bits >> k) & 1)
                adj[j * n + i] = 1;
            else
                adj[i * n + j] = 1;
        }
        long long H = held_karp(adj, n);
        if (H < hist_size)
            hist[H]++;
    }
}
"""

def main():
    c_path = "/tmp/h_spectrum_all.c"
    so_path = "/tmp/h_spectrum_all.so"
    with open(c_path, "w") as f:
        f.write(C_CODE)
    os.system(f"cc -O3 -o {so_path} -shared -fPIC {c_path}")
    lib = ctypes.CDLL(so_path)

    lib.h_spectrum_n.restype = None

    hist_size = 2000
    all_achieved = {}

    for n in range(3, 8):
        hist = (ctypes.c_longlong * hist_size)()
        print(f"n={n} ({2**(n*(n-1)//2)} tournaments)...", end=" ", flush=True)
        start = time.time()
        lib.h_spectrum_n(n, hist, hist_size)
        elapsed = time.time() - start

        achieved = set()
        for h in range(hist_size):
            if hist[h] > 0:
                achieved.add(h)
        all_achieved[n] = achieved

        max_h = max(achieved)
        odd_gaps = [h for h in range(1, max_h + 1, 2) if h not in achieved]
        print(f"done ({elapsed:.1f}s). {len(achieved)} distinct H values, range [1,{max_h}]")
        print(f"  Odd gaps in [1..{max_h}]: {odd_gaps}")

    # Find permanent gaps (missing at ALL n where they could appear)
    print("\n" + "=" * 60)
    print("PERMANENT GAP ANALYSIS")
    print("=" * 60)

    # Values that are gaps at every n checked
    max_check = 50  # check gaps up to H=50
    for h in range(1, max_check + 1, 2):
        statuses = []
        for n in range(3, 8):
            max_h_n = max(all_achieved[n])
            if h > max_h_n:
                statuses.append(f"n={n}:N/A")
            elif h in all_achieved[n]:
                statuses.append(f"n={n}:YES")
            else:
                statuses.append(f"n={n}:GAP")

        is_permanent_gap = all(
            h > max(all_achieved[n]) or h not in all_achieved[n]
            for n in range(3, 8)
        )
        # Only show if it's relevant (achievable at some n or a gap at all n)
        relevant_ns = [n for n in range(3, 8) if h <= max(all_achieved[n])]
        if relevant_ns:
            tag = " *** PERMANENT GAP ***" if is_permanent_gap else ""
            gap_at = [n for n in relevant_ns if h not in all_achieved[n]]
            achieved_at = [n for n in relevant_ns if h in all_achieved[n]]
            if gap_at:
                print(f"  H={h}: gap at n={gap_at}, achieved at n={achieved_at}{tag}")


if __name__ == "__main__":
    main()
