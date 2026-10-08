#!/usr/bin/env python3
"""Odd-step counts sigma(M_K) and Terras lengths of M_K = 2^K - 1 to first 1, K = 2..KMAX (mac-mini-2026-10-07 continuation).
Writes K, odd steps, Terras steps per line. Used for the orphan law (equal-time Mersenne partners)."""
import sys
KMAX = int(sys.argv[1]); out = sys.argv[2]
with open(out, 'w') as f:
    for K in range(2, KMAX + 1):
        n = (1 << K) - 1; odd = 0; ter = 0
        while n != 1:
            if n & 1:
                n = (3*n + 1) >> 1; odd += 1; ter += 1
            else:
                z = (n & -n).bit_length() - 1
                n >>= z; ter += z
        f.write(f"{K} {odd} {ter}\n")
        if K % 500 == 0: f.flush()
