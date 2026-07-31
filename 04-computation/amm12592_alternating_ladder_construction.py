#!/usr/bin/env python3
"""AMM 12592: the ladder over the FULL atom family, via the modulation P_k.

DOMINATION IS A MODULATION.  E_k is dominated by (1+u)^{a_k} exactly when

        [u^i] E_k = C(a_k, i) * P_k(a_k - i),      |P_k| <= 1.

So the legal atoms are the binomial profile modulated by any function bounded
by 1.  Two extremes:  P_k = +-1 gives E_k = +-(1+u)^{a_k}; P_k(y) = +-(-1)^y
gives E_k = +-(u-1)^{a_k}, the ALTERNATING atom.  In w = u/(1+u) coordinates
these are Lam_k = +-1 and Lam_k = +-(2w-1)^{a_k}.

WHY ALTERNATION IS THE RIGHT ATOM.  Lam_k(1) = P_k(0) is forced to +-1, and
Lam'_k(1) = -a_k (P_k(1) - P_k(0)), so |Lam'_k(1)| <= 2 a_k with equality
exactly at the alternating atom -- twice the reach of the unit atoms
u^c(1+u)^{a-c}, and it is |Delta^j P(0)| <= 2^j that saturates (ARCH).

THE LADDER.  Taking P_k integer-valued in {-1,0,1} makes every coefficient
C(a_k,i)P_k(a_k-i) an integer automatically.  Expanding the system in powers
of (1-w) turns it into, for each d,

        sum_{k=0}^{d} C(a_k, d-k) * Delta^{d-k} P_k(0) = eps(-1)^d C(m-1,d),

and since Delta^j P(0) = sum_i (-1)^{j-i} C(j,i) P(i), the variable
P_k(d-k) is NEW at level d and enters with coefficient C(a_k, d-k).  So each
level is a SIGNED SUBSET-SUM with coins C(a_k,d-k), k = 0..d, over
{-1,0,+1} -- a clean triangular search.
"""
import sys
from fractions import Fraction
from math import comb

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(1000000)


def profile(m, C):
    a = []
    for k in range(m):
        n = m + k
        cap = Fraction(C) * n - n - 1
        if cap < 0:
            return None
        a.append(min(m - 1 - k, int(cap)))
    return a


def solve_level(coins, target, cap_nodes=200000):
    """x_k in {-1,0,1} with sum coins[k]*x_k == target.  Greedy by size, then
    bounded DFS.  Returns the vector or None."""
    idx = sorted(range(len(coins)), key=lambda i: -coins[i])
    n = len(idx)
    suffix = [0] * (n + 1)
    for i in range(n - 1, -1, -1):
        suffix[i] = suffix[i + 1] + coins[idx[i]]
    x = [0] * len(coins)
    nodes = 0

    def rec(i, rem):
        nonlocal nodes
        if abs(rem) > suffix[i]:
            return False
        if i == n:
            return rem == 0
        nodes += 1
        if nodes > cap_nodes:
            raise RuntimeError("nodes")
        c = coins[idx[i]]
        order = (1, 0, -1) if rem > 0 else ((-1, 0, 1) if rem < 0 else (0, 1, -1))
        for v in order:
            x[idx[i]] = v
            if rec(i + 1, rem - c * v):
                return True
        x[idx[i]] = 0
        return False

    try:
        return x if rec(0, target) else None
    except RuntimeError:
        return None


def build(m, a, eps=1):
    """Returns P (list of lists) or None."""
    P = [[None] * (a[k] + 1) for k in range(m)]
    for d in range(m):
        ks = [k for k in range(min(d, m - 1) + 1) if d - k <= a[k]]
        if not ks:
            continue
        # known part: contributions of already-fixed P_k(i), i < d-k
        known = 0
        for k in ks:
            j = d - k
            s = 0
            for i in range(j):                     # i < j are already fixed
                s += (-1) ** (j - i) * comb(j, i) * P[k][i]
            known += comb(a[k], j) * s
        target = eps * (-1) ** d * comb(m - 1, d) - known
        coins = [comb(a[k], d - k) for k in ks]
        sol = solve_level(coins, target)
        if sol is None:
            return None, d
        for k, v in zip(ks, sol):
            P[k][d - k] = v
    for k in range(m):
        for i in range(a[k] + 1):
            if P[k][i] is None:
                P[k][i] = 0
    return P, None


def verify(m, a, P, eps):
    """Rebuild E_k, check the box and the identity sum_k E_k (1+u)^{L_k}."""
    L = [m - 1 - k - a[k] for k in range(m)]
    acc = [0] * (2 * m + 2)
    for k in range(m):
        for i in range(a[k] + 1):
            e = comb(a[k], i) * P[k][a[k] - i]
            if abs(e) > comb(a[k], i):
                return "BOX VIOLATED"
            if e:
                for t in range(L[k] + 1):
                    acc[i + t] += e * comb(L[k], t)
    tgt = [0] * (2 * m + 2)
    tgt[m - 1] = eps
    return "OK" if acc == tgt else "IDENTITY FAILS"


if __name__ == "__main__":
    print("ladder over the full (modulation) atom family")
    print("  m   best verified slope C           profile a_0")
    CANDS = [Fraction(p, q) for p, q in
             ((8, 5), (13, 8), (5, 3), (12, 7), (7, 4), (11, 6), (15, 8), (2, 1))]
    for m in (4, 8, 16, 32, 64):
        best = None
        for C in CANDS:
            a = profile(m, C)
            if a is None:
                continue
            for eps in (1, -1):
                P, fail = build(m, a, eps)
                if P is not None and verify(m, a, P, eps) == "OK":
                    best = (C, a[0])
                    break
            if best:
                break
        print(f"{m:4d}   " + (f"{str(best[0]):8s} = {float(best[0]):.4f}"
                              f"          a_0={best[1]}" if best else "none of the tested slopes"))
        sys.stdout.flush()
