"""
Reproduce the primorial-family claim and explore.

Claim: a family achieving max-min level a/(a(k+1)-1) where g*k^2 ~ 1/a,
when k-1 is a primorial (k-1 in {6,30,210,...}), with a ~ omega(k-1)
(number of distinct small prime factors of k-1).

We need a CONSTRUCTION. The natural LRC near-extremal set: take the AP {1,..,k}
but it gives exactly the floor (tight). Non-tight near-extremal sets are obtained
by perturbing. Let's first just measure level a := M/(1 - (k) * M)??
Define level a from M = a/(a(k+1)-1)  =>  solve for a:
   M*(a(k+1)-1) = a
   M a (k+1) - M = a
   a (M(k+1) - 1) = M
   a = M / (M(k+1) - 1)
So given M and k, a = M / (M*(k+1) - 1). For M slightly above 1/(k+1), denominator
M(k+1)-1 is small positive, a is large. g = M - 1/(k+1), g*k^2 ~ 1/a.

So "level a large" == "M close to floor". a >= sqrt(k) means M - 1/(k+1) <= ~1/(sqrt(k)*k^2).
"""
import sys
from fractions import Fraction
from math import gcd
sys.path.insert(0, ".")
from lrc_maxmin import M_exact


def level_a(M, k):
    den = M * (k + 1) - 1
    if den <= 0:
        return None
    return M / den


def report(S, k=None):
    S = sorted(set(S))
    if k is None:
        k = len(S)
    M, t = M_exact(S)
    floor = Fraction(1, k + 1)
    g = M - floor
    a = level_a(M, k)
    gk2 = g * k * k
    print(f"k={k} |S|={len(S)} M={M} (~{float(M):.8f}) floor={floor} g={g} g*k^2={float(gk2):.6f} a={float(a) if a else None}")
    return M, g, gk2, a


if __name__ == "__main__":
    # AP gives exactly the floor (tight) -> excluded from sigma_2.
    # Try perturbations: replace one element of {1..k}.
    pass
