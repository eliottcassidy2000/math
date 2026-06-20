"""
Convergence probe for I_B(u,v) -> Phi_2(B).

We scan growing coprime pairs (u,v) that are "generic":
  - gcd(u,v)=1, both coprime to 7
  - u,v large compared to max(B)
  - |u-v| large, and v not a small rational multiple of u
We report I_B(u,v) - Phi_2 and look at whether it shrinks as scale grows.
"""

from fractions import Fraction as F
from math import floor, gcd

SECTORS = 7


def measure_by_misscount(E):
    pts = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        denom = SECTORS * e
        for k in range(0, denom + 1):
            pts.add(F(k, denom))
    pts = sorted(pts)
    out = {}
    for i in range(len(pts) - 1):
        a, b = pts[i], pts[i + 1]
        if b == a:
            continue
        mid = (a + b) / 2
        hit = set()
        for e in E:
            t = e * mid
            t = t - floor(t)
            hit.add(floor(t * SECTORS))
        miss = SECTORS - len(hit)
        out[miss] = out.get(miss, F(0)) + (b - a)
    return out


def p0(E):
    return measure_by_misscount(E).get(0, F(0))


def I_B(B, u, v):
    Bt = tuple(B)
    return p0(Bt + (u, v)) - p0(Bt + (u,)) - p0(Bt + (v,)) + p0(Bt)


def Phi2(B):
    d = measure_by_misscount(B)
    p1 = d.get(1, F(0))
    p2 = d.get(2, F(0))
    return (2 * p2 - p1) / 49


def generic_pairs():
    """A sequence of growing coprime pairs designed to avoid small relations.
    Use Fibonacci-like ratio (golden) scaled up, forcing coprime & coprime-to-7."""
    pairs = []
    # golden ratio convergents (Fibonacci) are the LEAST resonant ratios
    fib = [13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181, 6765,
           10946, 17711, 28657, 46368, 75025]
    for i in range(len(fib) - 1):
        u, v = fib[i], fib[i + 1]
        if gcd(u, v) == 1 and u % 7 != 0 and v % 7 != 0:
            pairs.append((u, v))
    return pairs


if __name__ == "__main__":
    for B in [(0, 1, 2, 3, 4, 5, 6, 7), (0, 1, 2, 4, 8)]:
        phi = Phi2(B)
        print("=" * 72)
        print(f"B = {B}    Phi_2 = {phi} ~ {float(phi):.12f}")
        print("=" * 72)
        print(f"{'(u,v)':>16} {'I_B(u,v)':>20} {'I_B - Phi_2':>16}")
        for (u, v) in generic_pairs():
            iv = I_B(B, u, v)
            diff = iv - phi
            print(f"{str((u,v)):>16} {float(iv):>20.12f} {float(diff):>16.2e}")
