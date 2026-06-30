#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Covering-min as an ILP over the danger-circulant (set-cover on Z/m). mac-mini-2026-06-30-S51.

For witness modulus m and level j: a primitive covering (n-1)-set S has max_{i in Z/m} min_{v in S}||v i/m||
<= j/m  iff  for every i in Z/m, some v in S has dist(v i mod m, 0) <= j  (S's level-j danger arcs COVER Z/m).
This is a SET COVER ILP (HiGHS via scipy.milp). Binary-search the smallest feasible j -> mu(m)=j/m, a LOWER
bound on the covering-min M*; the returned S is then VERIFIED by its true M(S) (all breakpoints). If
M(S)=mu(m), then M* = mu(m) exactly. Min over m pins the covering-min, esp. for ODD n (the frontier).

The set-cover LP-dual is a PACKING = a fractional INDEPENDENT SET in the danger circulant (the lonely-point
packing); its independence polynomial is the OCF-type object -> the chromatic<->OCF bridge, computational.
"""
from __future__ import annotations
import functools, math
from fractions import Fraction as F
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
print = functools.partial(print, flush=True)


def Mexact(S):
    S = sorted(set(S)); cand = set()
    for i in range(len(S)):
        for j in range(len(S)):
            for d in (S[i]-S[j], S[i]+S[j]):
                if d > 0:
                    for k in range(1, d):
                        cand.add(F(k, d))
    best = F(0); at = None
    for t in cand:
        g = min(min((v*t) % 1, 1-((v*t) % 1)) for v in S)
        if g > best:
            best = g; at = t
    return best, at


def cdist(a, m):
    a %= m
    return min(a, m-a)


def feasible(n, m, j, V):
    """ILP: exists primitive-or-not covering (n-1)-set, speeds in 1..V, level-j danger arcs cover Z/m?
    Returns the set S (list) or None."""
    speeds = list(range(1, V+1))
    P = len(speeds)
    rows = []
    lb = []
    ub = []
    # size: sum x = n-1
    rows.append([1.0]*P); lb.append(n-1); ub.append(n-1)
    # covering: for q=2..n, sum_{q|v} x >= 1
    for q in range(2, n+1):
        rows.append([1.0 if v % q == 0 else 0.0 for v in speeds]); lb.append(1); ub.append(P)
    # danger-cover: for i=1..m-1, sum_{v: cdist(v i, m) <= j} x >= 1
    for i in range(1, m):
        rows.append([1.0 if cdist(v*i, m) <= j else 0.0 for v in speeds]); lb.append(1); ub.append(P)
    A = np.array(rows)
    cons = LinearConstraint(A, lb, ub)
    res = milp(c=np.zeros(P), constraints=cons, integrality=np.ones(P),
               bounds=Bounds(0, 1), options={"time_limit": 30})
    if res.success and res.x is not None:
        return [speeds[i] for i in range(P) if res.x[i] > 0.5]
    return None


def covering_min_ip(n, m_range, V):
    best = (F(1), None, None)
    for m in m_range:
        # binary-search smallest feasible j in [1, m//2]
        lo, hi, found = 1, m//2, None
        # linear scan up (small j first) -- the min is small
        for j in range(1, m//2 + 1):
            S = feasible(n, m, j, V)
            if S is not None:
                found = (j, S); break
        if found is None:
            continue
        j, S = found
        mu = F(j, m)
        trueM, at = Mexact(S)
        tight = (trueM == mu)
        if trueM < best[0]:
            best = (trueM, m, S)
        print(f"    m={m:>3}: mu(m)=j/m={j}/{m}={float(mu):.5f}  S={S}  trueM={trueM}={float(trueM):.5f}"
              f"  tight:{tight}  prim:{math.gcd(*S) if len(S)>1 else S[0]}=={'1' if math.gcd(*S)==1 else 'NO'}")
    return best


def main():
    print("COVERING-MIN via ILP (set-cover on the danger circulant Z/m)\n")
    for n, m_range, V in [(9, range(15, 36), 40), (11, range(19, 40), 44), (13, range(23, 46), 52)]:
        print(f"n={n}: 1/n={1/n:.5f}, scanning m in [{m_range.start},{m_range.stop-1}], speeds<= {V}:")
        bM, bm, bS = covering_min_ip(n, m_range, V)
        if bS:
            print(f"  => best (verified) covering set: M={bM}={float(bM):.5f} at m={bm}, S={bS}, "
                  f"margin={float(bM)-1/n:.5f}, primitive={math.gcd(*bS)==1}\n")
        else:
            print("  => none found\n")


if __name__ == "__main__":
    main()
