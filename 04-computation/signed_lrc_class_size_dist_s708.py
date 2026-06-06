"""Collision-class SIZE DISTRIBUTION for signed-LRC AP_n sign-orbit.  monad-explorer-S708.

Goal: crack the prime-power/mixed deficiency count law (THM-417 handoff).
Instead of only counting orbits, record the SIZE of each collision class.
Two cuts eps, eps' collide iff (Phi(eps)_t^2)_t == (Phi(eps')_t^2)_t  (folded clock-multiset),
Phi = S eps, S[t,i]=sin(2pi t i/C),  t,i in 1..n-1.  Fix eps_0=+1 (global-flip quotient).

deficiency = total_cuts - n_classes = sum_classes (size-1).
If silent moves form a group G_eps acting freely on each class, every class size is a power of 2
and deficiency = sum_eps (1 - 1/|G_eps|).  We tabulate the size histogram to expose the structure.
"""
import numpy as np
import math
from collections import Counter


def is_prime(C):
    return C > 1 and all(C % d for d in range(2, int(C**0.5) + 1))


def factor_str(C):
    f, m, d = [], C, 2
    while d * d <= m:
        while m % d == 0:
            f.append(d); m //= d
        d += 1
    if m > 1:
        f.append(m)
    return "x".join(map(str, f))


def class_sizes(n, batch=1 << 18):
    C = 2 * n - 1
    n1 = n - 1
    H = np.arange(1, n1 + 1)
    S = np.sin(2 * math.pi * np.outer(H, H) / C)  # S[t-1,i-1]
    nfree = n1 - 1
    total = 1 << nfree
    counts = Counter()
    for start in range(0, total, batch):
        end = min(start + batch, total)
        k = end - start
        idx = np.arange(start, end, dtype=np.uint64)
        eps = np.ones((n1, k), dtype=np.int8)
        for b in range(nfree):
            bit = ((idx >> np.uint64(b)) & np.uint64(1)).astype(np.int8)
            eps[b + 1, :] = np.where(bit == 1, 1, -1)
        Phi = S @ eps
        Phi2 = np.round(Phi * Phi, 6)
        for row in Phi2.T.copy():
            counts[row.tobytes()] += 1
    # histogram of class sizes
    hist = Counter(counts.values())
    n_classes = len(counts)
    defic = total - n_classes
    return C, total, n_classes, defic, hist


print(f"{'n':>3} {'C':>4} {'fact':>9} {'total':>9} {'classes':>9} {'defic':>7}  size-histogram {{size:count}}")
for n in range(5, 24):
    C = 2 * n - 1
    if is_prime(C):
        continue
    if (n - 2) > 21:
        continue
    Cc, total, ncl, defic, hist = class_sizes(n)
    hs = {int(s): int(c) for s, c in sorted(hist.items())}
    # sanity: sum size*count == total ; sum (size-1)*count == defic
    chk_tot = sum(s * c for s, c in hist.items())
    chk_def = sum((s - 1) * c for s, c in hist.items())
    flag = "" if (chk_tot == total and chk_def == defic) else "  <<CHECK FAIL>>"
    print(f"{n:>3} {C:>4} {factor_str(C):>9} {total:>9} {ncl:>9} {defic:>7}  {hs}{flag}", flush=True)
