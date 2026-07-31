#!/usr/bin/env python3
"""rho(m) = minimal within-block deadline ratio for the dyadic block [m,2m).

System (*):  sum_{k,i} p_{k,i} C(L_k, j-1-i) = C(m,j)/2,  1 <= j <= m-1,
             0 <= p_{k,i} <= C(a_k,i)  integers,   L_k = m-1-k-a_k.

Solved with HiGHS (scipy.optimize.milp); every FEASIBLE answer is then
re-verified in exact integer arithmetic, and every answer is cross-checked
against an exact DFS whenever the DFS search space is small enough.
"""
import numpy as np
from fractions import Fraction
from itertools import product
from math import comb
from scipy.optimize import milp, LinearConstraint, Bounds


def build(m, a):
    L = [m - 1 - k - a[k] for k in range(m)]
    assert all(x >= 0 for x in L)
    varlist = [(k, i) for k in range(m) for i in range(a[k] + 1)]
    A = np.zeros((m - 1, len(varlist)))
    for c, (k, i) in enumerate(varlist):
        for j in range(1, m):
            t = j - 1 - i
            if 0 <= t <= L[k]:
                A[j - 1, c] = comb(L[k], t)
    b = np.array([comb(m, j) // 2 for j in range(1, m)], dtype=float)
    ub = np.array([comb(a[k], i) for (k, i) in varlist], dtype=float)
    return varlist, A, b, ub, L


def feasible_ilp(m, a):
    varlist, A, b, ub, L = build(m, a)
    res = milp(c=np.zeros(len(varlist)),
               constraints=LinearConstraint(A, b, b),
               bounds=Bounds(np.zeros(len(varlist)), ub),
               integrality=np.ones(len(varlist)))
    if res.x is None:
        return False, None
    x = [int(round(v)) for v in res.x]
    # exact re-verification
    for j in range(1, m):
        s = sum(x[c] * (comb(L[k], j - 1 - i) if 0 <= j - 1 - i <= L[k] else 0)
                for c, (k, i) in enumerate(varlist))
        if s != comb(m, j) // 2:
            return False, None
    for c, (k, i) in enumerate(varlist):
        if not (0 <= x[c] <= comb(a[k], i)):
            return False, None
    return True, dict(zip(varlist, x))


def feasible_exact(m, a, cap=3_000_000):
    """Exact DFS: variables grouped by the highest equation they touch."""
    L = [m - 1 - k - a[k] for k in range(m)]
    rhs = [0] + [comb(m, j) // 2 for j in range(1, m)]
    groups = {}
    for k in range(m):
        for i in range(a[k] + 1):
            groups.setdefault(m - k - a[k] + i, []).append((k, i))
    order = sorted(groups, reverse=True)
    nodes = 0
    val = {}

    def contrib(j):
        s = 0
        for (k, i), v in val.items():
            t = j - 1 - i
            if v and 0 <= t <= L[k]:
                s += v * comb(L[k], t)
        return s

    def rec(idx):
        nonlocal nodes
        if idx == len(order):
            return all(contrib(j) == rhs[j] for j in range(1, m))
        j = order[idx]
        vs = groups[j]
        for combo in product(*[range(comb(a[k], i) + 1) for (k, i) in vs]):
            nodes += 1
            if nodes > cap:
                raise RuntimeError("cap")
            for (kv, c) in zip(vs, combo):
                val[kv] = c
            if 1 <= j <= m - 1 and contrib(j) != rhs[j]:
                continue
            if any(contrib(jj) > rhs[jj] for jj in range(1, m)):
                continue
            if rec(idx + 1):
                return True
        for kv in vs:
            val.pop(kv, None)
        return False

    return rec(0), nodes


def rho(m, cross_check=True):
    cands = sorted({Fraction(m + k + 1 + aa, m + k)
                    for k in range(m) for aa in range(m - k)})

    def prof(C):
        a = []
        for k in range(m):
            n = m + k
            cap = Fraction(C) * n - n - 1
            if cap < 0:
                return None
            a.append(min(m - 1 - k, int(cap)))
        return a

    lo, hi, best = 0, len(cands) - 1, None
    while lo <= hi:
        mid = (lo + hi) // 2
        a = prof(cands[mid])
        ok = False
        if a is not None:
            ok, w = feasible_ilp(m, a)
            if cross_check:
                try:
                    ok2, nodes = feasible_exact(m, a)
                    assert ok == ok2, ("ILP/DFS disagree", m, a, ok, ok2)
                except RuntimeError:
                    pass
        if ok:
            best = (cands[mid], a, w)
            hi = mid - 1
        else:
            lo = mid + 1
    return best


if __name__ == "__main__":
    print("m=4 sanity (vs exhaustive brute force):")
    for a, want in (((1, 1, 0, 0), True), ((2, 0, 0, 0), True),
                    ((1, 0, 0, 0), False), ((0, 0, 0, 0), False)):
        got, w = feasible_ilp(4, list(a))
        print(f"   a={a}  feasible={got}  expected={want}  "
              f"{'OK' if got == want else 'MISMATCH'}")
    print()
    for r in range(2, 7):
        m = 2 ** r
        res = rho(m, cross_check=(m <= 8))
        if res is None:
            print(f"m={m:4d}: no feasible profile")
            continue
        C, a, w = res
        T = [m + k + 1 + a[k] for k in range(m)]
        print(f"m={m:4d}:  rho(m) = {C} = {float(C):.5f}")
        print(f"        a = {a}")
        print(f"        T = {T}")
        print(f"        argmax ratio at n = {m + max(range(m), key=lambda k: Fraction(T[k], m+k))}")
