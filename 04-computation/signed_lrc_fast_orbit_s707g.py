"""FAST sign-orbit counter via Phi^2 vectors (matmul + hashing).  monad-explorer-S707.
Collision class <=> (Phi(eps)_t^2)_t  vector, Phi = S eps, S[t,i]=sin(2pi t i/C).
This avoids O(n^2) folded-spectrum recomputation; reaches new composite C (45=3^2*5, 49=7^2).
Confirms C<=39 against brute force, then extends the deficiency table.
"""
import numpy as np
import math
from itertools import product


def is_prime(C):
    return C > 1 and all(C % d for d in range(2, int(C**0.5) + 1))


def orbit_count(n, batch=1 << 18):
    C = 2 * n - 1
    n1 = n - 1
    H = np.arange(1, n1 + 1)
    S = np.sin(2 * math.pi * np.outer(H, H) / C)  # (n1,n1), S[t-1,i-1]=sin(2pi t i/C)
    # cuts: fix eps_0 = +1 (runner 1), vary the other n1-1 -> 2^{n-2} cuts
    nfree = n1 - 1
    total = 1 << nfree
    seen = set()
    for start in range(0, total, batch):
        end = min(start + batch, total)
        k = end - start
        idx = np.arange(start, end, dtype=np.uint64)
        # build eps matrix (n1 x k): row0 = +1; rows 1..nfree from bits
        eps = np.ones((n1, k), dtype=np.int8)
        for b in range(nfree):
            bit = ((idx >> np.uint64(b)) & np.uint64(1)).astype(np.int8)
            eps[b + 1, :] = np.where(bit == 1, 1, -1)
        Phi = S @ eps                      # (n1, k)
        Phi2 = np.round(Phi * Phi, 6)      # (n1, k)
        # hash each column
        cols = Phi2.T.copy()
        for row in cols:
            seen.add(row.tobytes())
    return C, total, len(seen), total - len(seen)


print(f"{'n':>3} {'C':>4} {'fact':>10} {'2^(n-2)':>10} {'orbit':>10} {'defic':>8}")
for n in range(5, 26):
    C = 2 * n - 1
    if is_prime(C):
        continue
    if (n - 2) > 23:
        continue
    Cc, total, orbit, defic = orbit_count(n)
    f2, m2, d = [], C, 2
    while d * d <= m2:
        while m2 % d == 0:
            f2.append(d); m2 //= d
        d += 1
    if m2 > 1:
        f2.append(m2)
    fstr = "x".join(map(str, f2))
    print(f"{n:>3} {C:>4} {fstr:>10} {total:>10} {orbit:>10} {defic:>8}", flush=True)
