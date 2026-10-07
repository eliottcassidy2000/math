#!/usr/bin/env python3
"""Chain-free check: Terras orbits of y and u = 3^k y + e (exact 2-adic arithmetic mod 2^K,
y uniform), count equal-time merges T^n(u) == T^n(y) within n <= Nmax.
For k < 0 or non-integer e we work mod 2^K with 3^-1 as a 2-adic unit."""
import random, math, sys
K = 9000; Nmax = 8000          # Nmax < K: each Terras step consumes one bit of precision
def run(k, e_num, e_den_pow3, trials, seed):
    rng = random.Random(seed); merged = []; mod0 = 1 << K
    inv3 = pow(3, -1, mod0)
    for _ in range(trials):
        y = rng.getrandbits(K); m = mod0
        slope = pow(3, k, m) if k >= 0 else pow(inv3, -k, m)
        u = (slope * y + e_num * pow(inv3, e_den_pow3, m)) % m
        v = y; hit = None
        for n in range(1, Nmax + 1):
            v = v // 2 if v % 2 == 0 else (3 * v + 1) // 2
            u = u // 2 if u % 2 == 0 else (3 * u + 1) // 2
            m >>= 1; v %= m; u %= m
            if u == v and m > (1 << 64):      # equal mod 2^(K-n) with >64 bits of margin
                hit = n; break
        merged.append(hit)
    return merged
cases = [("y vs y+1", 0, 1, 0), ("y vs 3y", 1, 0, 0), ("y vs y/3", -1, 0, 0), ("y vs 9y+5", 2, 5, 0),
         ("y vs y-7", 0, -7, 0), ("y vs 27y+1/3", 3, 1, 1)]
trials = int(sys.argv[1]) if len(sys.argv) > 1 else 120
for name, k, en, ed in cases:
    res = run(k, en, ed, trials, 11 + k)
    hits = sorted(h for h in res if h is not None)
    frac = len(hits) / trials
    med = hits[len(hits)//2] if hits else None
    print(f"{name:14s}: merged within {Nmax} steps: {len(hits)}/{trials} = {frac:.3f}; median merge time {med}; "
          f"chain-predicted unmerged ~ c/sqrt({Nmax}) = c*{1/math.sqrt(Nmax):.4f}")
    sys.stdout.flush()
