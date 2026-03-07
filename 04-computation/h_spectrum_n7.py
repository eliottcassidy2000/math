#!/usr/bin/env python3
"""
Complete H-spectrum at n=7 (exhaustive).

Uses C via ctypes for speed. Records all H values that occur,
identifies gaps, and checks which odd values are missing.

Instance: opus-2026-03-07-S40
"""

import ctypes
import os
import time
from collections import Counter

C_CODE = r"""
#include <string.h>

long long held_karp_7(int adj[7]) {
    int n = 7;
    int full = (1 << n) - 1;
    long long dp[128][7];
    memset(dp, 0, sizeof(dp));
    for (int v = 0; v < n; v++)
        dp[1 << v][v] = 1;
    for (int S = 1; S <= full; S++) {
        for (int v = 0; v < n; v++) {
            if (!(S & (1 << v))) continue;
            long long c = dp[S][v];
            if (c == 0) continue;
            int out = adj[v] & (~S);
            while (out) {
                int u = __builtin_ctz(out);
                dp[S | (1 << u)][u] += c;
                out &= out - 1;
            }
        }
    }
    long long total = 0;
    for (int v = 0; v < n; v++)
        total += dp[full][v];
    return total;
}

/* Compute H for all 2^21 tournaments, store histogram */
void h_spectrum(long long *hist, int hist_size) {
    int n = 7, m = 21;
    int edges[21][2];
    int ei = 0;
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++) {
            edges[ei][0] = i;
            edges[ei][1] = j;
            ei++;
        }

    long long total = 1L << m;
    for (long long bits = 0; bits < total; bits++) {
        int adj[7] = {0};
        for (int k = 0; k < m; k++) {
            int i = edges[k][0], j = edges[k][1];
            if ((bits >> k) & 1)
                adj[j] |= (1 << i);
            else
                adj[i] |= (1 << j);
        }
        long long H = held_karp_7(adj);
        if (H < hist_size)
            hist[H]++;
    }
}
"""

def main():
    c_path = "/tmp/h_spectrum_n7.c"
    so_path = "/tmp/h_spectrum_n7.so"
    with open(c_path, "w") as f:
        f.write(C_CODE)
    os.system(f"cc -O3 -o {so_path} -shared -fPIC {c_path}")
    lib = ctypes.CDLL(so_path)

    hist_size = 2000
    hist = (ctypes.c_longlong * hist_size)()

    print("Computing full H-spectrum at n=7 (2,097,152 tournaments)...")
    start = time.time()
    lib.h_spectrum(hist, hist_size)
    elapsed = time.time() - start
    print(f"Done in {elapsed:.1f}s\n")

    print("=" * 60)
    print("H-SPECTRUM AT n=7")
    print("=" * 60)

    all_h = []
    for h in range(hist_size):
        if hist[h] > 0:
            all_h.append((h, hist[h]))

    total = sum(c for _, c in all_h)
    print(f"Total tournaments: {total}")
    print(f"Distinct H values: {len(all_h)}")
    print(f"H range: {all_h[0][0]} to {all_h[-1][0]}")

    # All achieved values
    achieved = set(h for h, _ in all_h)

    # Print all values
    print(f"\nAll H values (count):")
    for h, c in all_h:
        print(f"  H={h}: {c}")

    # Find gaps in odd numbers
    max_h = all_h[-1][0]
    odd_achieved = sorted(h for h in achieved if h % 2 == 1)
    odd_missing = [h for h in range(1, max_h + 1, 2) if h not in achieved]

    print(f"\nOdd H values achieved: {len(odd_achieved)}")
    print(f"Odd H values missing (gaps) in [1..{max_h}]: {odd_missing}")

    # Check specific values
    for target in [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29]:
        status = "ACHIEVED" if target in achieved else "GAP"
        count = hist[target] if target in achieved else 0
        print(f"  H={target}: {status} ({count} tournaments)")

    # Check if H is always odd
    even_h = [h for h in achieved if h % 2 == 0]
    if even_h:
        print(f"\nEven H values found: {even_h}")
    else:
        print(f"\nConfirmed: H is always odd at n=7")


if __name__ == "__main__":
    main()
