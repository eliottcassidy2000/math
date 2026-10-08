#!/usr/bin/env python3
"""Odd-step counts and Terras lengths of M_K = 2^K - 1 to first 1 for K in [K0, K1] (appends to a file)."""
import sys
K0, K1, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
with open(out, 'a') as f:
    for K in range(K0, K1 + 1):
        n = (1 << K) - 1; odd = 0; ter = 0
        while n != 1:
            if n & 1:
                n = (3*n + 1) >> 1; odd += 1; ter += 1
            else:
                z = (n & -n).bit_length() - 1
                n >>= z; ter += z
        f.write(f"{K} {odd} {ter}\n")
        if K % 200 == 0: f.flush()
