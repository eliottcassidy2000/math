#!/usr/bin/env python3
"""Audit A, item 3: return drift E|e_ret|^theta vs rho|e|^theta + C; skeleton visits per level (claimed 2);
steps per visit at level h (claimed <= m_h + 2); C_exc. Excursions capped at CAP steps (counted separately)."""
import random, math, sys
from collections import defaultdict
P3 = [3 ** i for i in range(5000)]
def v2(x): return (x & -x).bit_length() - 1
def step(k, N, beta):
    sig = N & 1
    if k >= 1:
        if not sig: return (k, N >> 1) if not beta else (k, (3 * N + 1 - P3[k]) >> 1)
        return (k + 1, (3 * N + 1) >> 1) if not beta else (k - 1, (N - P3[k - 1]) >> 1)
    if k == 0:
        if not sig: return (0, N >> 1) if not beta else (0, (3 * N) >> 1)
        return (1, (3 * N + 1) >> 1) if not beta else (-1, (3 * N - 1) >> 1)
    h = -k
    if not sig: return (k, N >> 1) if not beta else (k, (3 * N + P3[h] - 1) >> 1)
    return (k + 1, (N + P3[h - 1]) >> 1) if not beta else (k - 1, (3 * N - 1) >> 1)

def excursion(e, rng, CAP, visits=None, steps_at=None):
    """from (0, e) run until the next arrival at k = 0 (after departing); return (e_ret, steps) or (None, steps)"""
    k, N = 0, e
    n = 0
    while k == 0 and (N & 1) == 0:
        if N == 0: return 0, 0
        k, N = step(k, N, rng.getrandbits(1)); n += 1
    # departure flip
    k, N = step(k, N, rng.getrandbits(1)); n += 1
    lastk = None
    while k != 0:
        h = abs(k)
        if visits is not None:
            if steps_at is not None: steps_at[h] += 1
            if (N & 1) and False: pass
        sig = N & 1
        if visits is not None and lastk != k:
            visits[h] += 1   # arrival at level h (a new visit)
        lastk = k
        k, N = step(k, N, rng.getrandbits(1)); n += 1
        if n > CAP: return None, n
    return N, n

if __name__ == '__main__':
    th = 0.5
    s = 2 ** th * (1 - math.sqrt(1 - 0.75 ** th)); rho = s * 2 ** -th
    mh = lambda h: v2(3 ** h - 1)
    Cexc = sum(2 * (mh(h) + 2) * 2 ** -th * s ** (h - 1) for h in range(1, 3000))
    C = s * 6 ** -th + Cexc
    print(f"theta=1/2: s={s:.6f} rho={rho:.6f}  C_exc={Cexc:.3f}  C={C:.3f}  C/(1-rho)={C/(1-rho):.2f}")
    rng = random.Random(99)
    CAP = 200000
    # (1) visits per level and steps per visit, from many excursions started at (0, odd e)
    visits = defaultdict(int); steps_at = defaultdict(int); nexc = 0; capped = 0
    for i in range(20000):
        e = rng.choice([1, -1, 3, 5, -7, 101, 2 ** 40 + 1, -(3 ** 30)])
        r, n = excursion(e, rng, CAP, visits, steps_at)
        nexc += 1; capped += r is None
    print(f"[visits] {nexc} excursions ({capped} capped at {CAP} steps); E[#visits to level h] (claim 2), steps/visit (claim <= m_h+2):")
    for h in (1, 2, 3, 4, 5, 6, 8, 12, 16, 24, 32):
        vis = visits[h] / nexc; spv = steps_at[h] / max(1, visits[h])
        print(f"   h={h:3d}: visits {vis:.3f}   steps/visit {spv:.3f}   m_h+2 = {mh(h)+2}")
    # (2) return drift for several starting e
    print("[drift] E|e_ret|^1/2 / |e|^1/2 against rho = %.4f, and E|e_ret|^1/2 - rho|e|^1/2 against C = %.2f" % (rho, C))
    for e in (1, -1, 3, 7, 21, 1001, -1001, 2 ** 20 + 1, 10 ** 12 + 39, -(10 ** 12 + 39), 3 ** 40, 2 ** 200 + 1):
        tot = 0.0; n = 0; cap = 0; M = 4000
        for i in range(M):
            r, st = excursion(e, rng, CAP)
            if r is None: cap += 1; continue
            tot += abs(r) ** th; n += 1
        m = tot / n
        print(f"   e0={e if abs(e) < 10**13 else ('2^200+1' if e == 2**200+1 else e):>16}: ratio {m / abs(e) ** th:.4f}   excess {m - rho * abs(e) ** th:8.3f}   (n={n}, capped={cap})")
