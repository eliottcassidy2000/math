"""
Exact Lonely Runner max-min collar M(S) computation.
Self-contained, from the definition. Uses fractions.Fraction.

M(S) = max_{t in (0,1)} min_{v in S} ||v t||,  ||x|| = dist(x, nearest int).

Method (crossing-candidates): M(S) is attained at t = m/d where d in {v_i+-v_j},
m=1..d-1. Evaluate min_v ||v t|| at each candidate, take max.

We also include candidates t = m/(2 v_i) and t = m/v_i (endpoints of each runner's
danger structure) to be safe; the binding-pair theorem guarantees the optimum is at
t = m/(v_i +- v_j) for some pair (including i=j giving 2 v_i, and the trivial v_i),
so the candidate set {m/d : d in D} with D = all |v_i +- v_j| and 2 v_i and v_i
is complete.
"""
from fractions import Fraction
from math import gcd
from functools import reduce


def frac_norm(t):
    # ||t|| = distance to nearest integer, t a Fraction
    f = t - int(t)            # fractional part in (-1,1)
    if f < 0:
        f += 1                # now in [0,1)
    return f if f <= Fraction(1, 2) else 1 - f


def min_collar_at(t, S):
    return min(frac_norm(v * t) for v in S)


def candidate_denoms(S):
    D = set()
    n = len(S)
    Sl = list(S)
    for i in range(n):
        vi = Sl[i]
        D.add(2 * vi)
        D.add(vi)
        for j in range(i + 1, n):
            vj = Sl[j]
            D.add(vi + vj)
            d = abs(vi - vj)
            if d > 0:
                D.add(d)
    return D


def M_exact(S):
    """Return M(S) as a Fraction (exact)."""
    S = sorted(set(S))
    D = candidate_denoms(S)
    best = Fraction(0)
    best_t = None
    for d in D:
        if d <= 1:
            continue
        for m in range(1, d):
            t = Fraction(m, d)
            val = min_collar_at(t, S)
            if val > best:
                best = val
                best_t = t
    return best, best_t


def lcm(a, b):
    return a * b // gcd(a, b)


def lcm_list(xs):
    return reduce(lcm, xs, 1)


if __name__ == "__main__":
    # sanity: AP {1,2,3} k=3, floor 1/4
    for k in range(2, 8):
        S = list(range(1, k + 1))
        M, t = M_exact(S)
        print(f"AP k={k}: M={M} (={float(M):.6f}), floor 1/(k+1)={Fraction(1,k+1)}, t*={t}")
