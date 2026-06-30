#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Covering-min as a CORRECT set-cover ILP. mac-mini-2026-06-30-S51 (v2, fixes the coarse-grid bug).

Universe of witnesses = ALL breakpoints tau=k/d, 1<=k<d<=2V (the max-min of ||v t|| is always attained at
such a tau, d a pairwise sum/diff of speeds <= 2V). M(S) <= r  iff  for every tau in the universe, some v in S
has ||v tau|| <= r (S's radius-r danger arcs cover the universe). Binary-search the smallest feasible r over
the finite set of attained distances => the EXACT covering-min. ILP per r via scipy.milp (HiGHS).

The set-cover (cover the witness universe by danger arcs) has LP-dual = a PACKING of pairwise-incompatible
lonely witnesses = an INDEPENDENT SET in the 'danger conflict graph'; its independence polynomial is the
OCF-type object (I(Omega,2)) -> the chromatic<->OCF bridge, computational.
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


def norm(x):  # ||x|| distance to nearest integer, x a Fraction
    f = x - int(x)
    if f < 0:
        f += 1
    return min(f, 1-f)


def covering_min_ip(n, V):
    speeds = list(range(1, V+1))
    P = len(speeds)
    # witness universe: reduced k/d, 1<=k<d<=2V
    univ = sorted({F(k, d) for d in range(2, 2*V+1) for k in range(1, d) if math.gcd(k, d) == 1})
    U = len(univ)
    # dist[p][u] = ||speeds[p]*univ[u]||
    dist = [[norm(speeds[p]*univ[u]) for u in range(U)] for p in range(P)]
    # candidate r values (sorted distinct attained distances), only those that could be the covering-min
    rvals = sorted({dist[p][u] for p in range(P) for u in range(U)})

    # fixed constraints: size + divisibility covering + PRIMITIVITY (gcd=1: for each prime p, some speed coprime to p)
    base_rows = [[1.0]*P]; base_lb = [n-1]; base_ub = [n-1]
    for q in range(2, n+1):
        base_rows.append([1.0 if v % q == 0 else 0.0 for v in speeds]); base_lb.append(1); base_ub.append(P)
    primes = [p for p in range(2, V+1) if all(p % d for d in range(2, int(p**0.5)+1))]
    for p in primes:  # at least one selected speed NOT divisible by p  => gcd(S)=1
        base_rows.append([1.0 if v % p != 0 else 0.0 for v in speeds]); base_lb.append(1); base_ub.append(P)

    def feasible(r):
        rows = [row[:] for row in base_rows]; lb = base_lb[:]; ub = base_ub[:]
        for u in range(U):
            rows.append([1.0 if dist[p][u] <= r else 0.0 for p in range(P)]); lb.append(1); ub.append(P)
        A = np.array(rows)
        res = milp(c=np.zeros(P), constraints=LinearConstraint(A, lb, ub),
                   integrality=np.ones(P), bounds=Bounds(0, 1), options={"time_limit": 25})
        if res.success and res.x is not None:
            return [speeds[i] for i in range(P) if res.x[i] > 0.5]
        return None

    # binary search smallest feasible r
    lo, hi = 0, len(rvals)-1
    best = None
    while lo <= hi:
        mid = (lo+hi)//2
        S = feasible(rvals[mid])
        if S is not None:
            best = (rvals[mid], S); hi = mid-1
        else:
            lo = mid+1
    return best


def main():
    print("COVERING-MIN via CORRECT set-cover ILP (universe = all breakpoints k/d, d<=2V)\n")
    for n, V in [(7, 22), (8, 24), (9, 38), (11, 46), (13, 56)]:
        res = covering_min_ip(n, V)
        if res is None:
            print(f"n={n}: infeasible/none (V={V})"); continue
        r, S = res
        trueM, at = Mexact(S)
        g = 0
        for v in S:
            g = math.gcd(g, v)
        print(f"n={n:>2} (V={V}): covering-min r={r}={float(r):.5f}  S={S}  trueM={trueM}={float(trueM):.5f}  "
              f"witness t*={at}  margin={float(trueM)-1/n:.5f}  primitive={g==1}  (1/n={1/n:.5f})")


if __name__ == "__main__":
    main()
