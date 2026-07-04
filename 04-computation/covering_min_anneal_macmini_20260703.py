#!/usr/bin/env python3
"""
STRONG covering-min search: numpy-fast approx-M annealing + EXACT rational M on finalists.
(mac-mini-2026-07-03-S30) Court case (convergent-not-covering-min) open item 3a: covering-min trajectory n>=10.
MUST reproduce n=7,8,9 winners (2/13,2/15,4/33), then extend to n=14. Quantify margin M/(1/n) (uniform looseness).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import numpy as np
import random

def gcd_all(xs): return reduce(gcd, xs)
def is_covering(sp, n): return all(any(v % q == 0 for v in sp) for q in range(2, n+1))

_TGRID = None
def approxM(sp):
    """fast float view M over a fine grid (guides annealing; not exact)."""
    global _TGRID
    if _TGRID is None:
        _TGRID = np.arange(1, 6000)/6000.0
    v = np.array(sp, dtype=np.float64)
    ph = np.outer(v, _TGRID) % 1.0          # (k, G)
    d = np.minimum(ph, 1.0-ph)              # dist to nearest int
    mn = d.min(axis=0)                       # min over speeds, per t
    return mn.max()

def nd(x):
    x = x % 1
    return x if x <= 1-x else 1-x
def M_exact(sp):
    """exact rational view M over the complete breakpoint set."""
    vs = sorted(set(sp)); Q = set()
    for i in range(len(vs)):
        Q.add(2*vs[i])
        for j in range(i+1, len(vs)):
            Q.add(vs[i]+vs[j]); Q.add(vs[j]-vs[i])
    best = F(0); bt = None
    for q in Q:
        if q < 2: continue
        for a in range(1, q):
            t = F(a, q); m = min(nd(v*t) for v in sp)
            if m > best: best, bt = m, t
    return best, bt

def repair_cover(sp, n, hi):
    sp = set(x for x in sp if x >= 1)
    for q in range(2, n+1):
        if not any(v % q == 0 for v in sp):
            for mult in range(1, hi//q+1):
                if q*mult not in sp: sp.add(q*mult); break
    return sorted(sp)

def anneal_one(seed, n, hi, iters, rng):
    k = n-1
    cur = sorted(set(seed))
    if len(cur) != k or not is_covering(cur, n) or gcd_all(cur) != 1: return (1.0, None)
    curM = approxM(cur); best = (curM, cur[:]); T = 0.4
    for it in range(iters):
        T *= 0.9994
        cand = cur[:]; idx = rng.randrange(k); old = cand[idx]
        new = old + rng.choice([-4,-3,-2,-1,1,2,3,4]) if rng.random() < 0.6 else rng.randint(1, hi)
        if new < 1 or new in cand: continue
        cand[idx] = new; cand = sorted(set(cand))
        if len(cand) != k or gcd_all(cand) != 1 or not is_covering(cand, n): continue
        candM = approxM(cand); dM = candM - curM
        if candM < curM or rng.random() < pow(2.718, -dM/max(T*0.05,1e-4)):
            cur, curM = cand, candM
            if curM < best[0]: best = (curM, cur[:])
    return best

if __name__ == "__main__":
    print("STRONG covering-min: numpy annealing + exact-M finalists. Validate n=7,8,9, extend to 14.")
    print("="*94)
    known = {7: F(2,13), 8: F(2,15), 9: F(4,33)}
    print(f"{'n':>3} {'1/n':>9} {'n/Phi6':>9} {'best M (exact)':>16} {'vs known':>16} {'M/(1/n)':>9} {'best set':>36}")
    for n in range(7, 15):
        rng = random.Random(7000+n); k = n-1
        phi6 = n*n-n+1; apex = (n-1)*n
        hi = max(3*n, 36) if n <= 9 else max(5*n, 70)
        seeds = [list(range(1, n-1))+[apex], list(range(1, n)), list(range(2, n+1))]
        for _ in range(60):
            base = repair_cover(sorted(set(rng.sample(range(1, hi), min(k, hi-1)))), n, hi)
            if len(base) >= k: seeds.append(base[:k] if len(base) > k else base)
        iters = 4000 if n <= 10 else 6000
        # collect top approx candidates, then exact-eval them
        results = []
        for sd in seeds:
            m, s = anneal_one(sd, n, hi, iters, rng)
            if s is not None: results.append((m, tuple(s)))
        results = sorted(set(results))[:25]  # exact-eval the 25 best-by-approx
        bestE = (F(1), None)
        for _, s in results:
            Me, _ = M_exact(list(s))
            if Me < bestE[0]: bestE = (Me, list(s))
        M, sset = bestE
        kv = known.get(n)
        vs = (f"{str(kv)}: " + ("MATCH" if M == kv else ("BEATS!" if M < kv else "above"))) if kv else "-"
        print(f"{n:>3} {float(F(1,n)):>9.5f} {float(F(n,phi6)):>9.5f} {str(M)+f'={float(M):.5f}':>16} "
              f"{vs:>16} {float(M)/(1.0/n):>9.4f} {str(sset):>36}")
    print("\n=> validation: n=7,8,9 -> 2/13,2/15,4/33 (or lower). n=14: is best M < 14/183=0.07650?")
    print("   M/(1/n) column = the uniform-looseness margin: if bounded ~1.07-1.10, primitive covering is LOOSE.")
