#!/usr/bin/env python3
"""Orphan law for general residual sources n = 2^K t - 1 (first reset 2), t of B random bits (mac-mini-2026-10-07 cont.).
A source is a deletion orphan if no child h_D = (n+1)/2^D - 1, 1 <= D <= K-1, has the same odd-step count to 1 and Terras
length sigma(h_D) = sigma(n) - D (equal-time merge before 1). Reports orphan fraction vs B for fixed K."""
import random, sys, math
def v2(x): return (x & -x).bit_length() - 1
def ot(n):
    o = 0; t = 0
    while n != 1:
        if n & 1: n = (3*n + 1) >> 1; o += 1; t += 1
        else:
            z = v2(n); n >>= z; t += z
    return o, t
rnd = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 7)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
for K in (9, 17):
    for B in (25, 50, 100, 200, 400, 800):
        orph = 0; d1 = 0; tot = 0
        while tot < N:
            t = rnd.getrandbits(B) | (1 << (B-1)) | 1
            if v2(3**K * t - 1) != 1: continue
            n = (1 << K) * t - 1
            on, tn = ot(n)
            found = False; f1 = False
            for D in range(1, K):
                h = (n + 1) // (1 << D) - 1
                oh, th = ot(h)
                if oh == on and tn - th == D:
                    found = True
                    if D <= 2: f1 = True
                    break
            tot += 1; orph += (not found); d1 += f1
        print(f"K={K:2d} B={B:4d}: orphan fraction {orph/N:.3f} (x sqrt(K+B) = {orph/N*math.sqrt(K+B):.2f});  reset pair D<=2 merges {d1/N:.3f}", flush=True)
