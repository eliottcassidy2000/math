#!/usr/bin/env python3
"""
Sample n=8 tournaments for H=21 and H=63 (and full low-range gaps).

Uses C Held-Karp for speed. 500k random tournaments.

Instance: opus-2026-03-07-S40
"""

import ctypes
import os
import time
import random
from collections import Counter

C_CODE = r"""
#include <string.h>

long long held_karp_8(int adj[8]) {
    int n = 8;
    int full = (1 << n) - 1;
    long long dp[256][8];
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

/* Process one tournament given edge bits. Return H value. */
long long process_one(long long bits) {
    int n = 8, m = 28;
    int edges[28][2];
    int ei = 0;
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++) {
            edges[ei][0] = i;
            edges[ei][1] = j;
            ei++;
        }

    int adj[8] = {0};
    for (int k = 0; k < m; k++) {
        int i = edges[k][0], j = edges[k][1];
        if ((bits >> k) & 1)
            adj[j] |= (1 << i);
        else
            adj[i] |= (1 << j);
    }
    return held_karp_8(adj);
}
"""

def main():
    c_path = "/tmp/h21_n8_sample.c"
    so_path = "/tmp/h21_n8_sample.so"
    with open(c_path, "w") as f:
        f.write(C_CODE)
    os.system(f"cc -O3 -o {so_path} -shared -fPIC {c_path}")
    lib = ctypes.CDLL(so_path)
    lib.process_one.restype = ctypes.c_longlong

    random.seed(42)
    num_samples = 500000
    m = 28

    hist = Counter()
    print(f"Sampling {num_samples} random tournaments at n=8...")
    start = time.time()

    for trial in range(num_samples):
        if trial % 100000 == 0 and trial > 0:
            print(f"  {trial}/{num_samples}")
        bits = random.randint(0, (1 << m) - 1)
        H = lib.process_one(bits)
        hist[H] += 1

    elapsed = time.time() - start
    print(f"Done in {elapsed:.1f}s\n")

    # Focus on low-range gaps
    print("=" * 60)
    print(f"H-SPECTRUM SAMPLE AT n=8 ({num_samples} tournaments)")
    print("=" * 60)

    all_h = sorted(hist.keys())
    print(f"H range: {all_h[0]} to {all_h[-1]}")
    print(f"Distinct H values seen: {len(all_h)}")

    # Show all values up to H=100
    print(f"\nH values up to 100:")
    for h in range(1, 101, 2):
        if h in hist:
            print(f"  H={h}: {hist[h]} times")
        else:
            print(f"  H={h}: NOT SEEN (potential gap)")

    # Key targets
    print(f"\nKey values:")
    for target in [7, 21, 63]:
        if target in hist:
            print(f"  H={target}: ACHIEVED ({hist[target]} times)")
        else:
            print(f"  H={target}: NOT SEEN in {num_samples} samples")

    # Gaps in odd values up to 200
    seen = set(hist.keys())
    max_seen = max(seen)
    odd_missing = [h for h in range(1, min(201, max_seen), 2) if h not in seen]
    print(f"\nOdd values NOT seen in [1..200]: {odd_missing}")


if __name__ == "__main__":
    main()
