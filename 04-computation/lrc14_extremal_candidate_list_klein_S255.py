# -*- coding: utf-8 -*-
# klein-2026-07-11-S255: THE EXTREMAL-CANDIDATE EVALUATOR (guard vs MISTAKE-127/138).
#
# The extremal families for LRC(14) functionals are ALGEBRAICALLY COHERENT (APs, mod-p residue
# grids, GW/doubling points) -- measure-zero-thin, aligned, INVISIBLE to any diameter box or local
# hill-climb (reflection: the-extremals-are-algebraically-special-invisible-to-local-search).
# BEFORE claiming any extremum, evaluate the functional on THIS explicit candidate list.
#
# Usage: import and call evaluate(func, k) where func(E)->value (E a tuple of speeds).
# Reports the value on every algebraic candidate + flags the argmin/argmax.

from math import gcd
from fractions import Fraction as F

def lcm(a, b): return a // gcd(a, b) * b

def candidates(k):
    """The algebraic extremal candidates for a k-core (the families sampling misses)."""
    C = {}
    C['consec (AP step1)'] = tuple(range(1, k + 1))
    for step in (2, 3, 7, 14):                      # AP grids incl. the mod-7 / mod-14 resonances
        C[f'AP step{step} off1'] = tuple(1 + step * i for i in range(k))
    for r in range(1, 7):                            # all six mod-7 residue grids (BUNCH poles)
        C[f'mod7 res{r}'] = tuple(r + 7 * i for i in range(k))
    C['mod14 res1'] = tuple(1 + 14 * i for i in range(k))
    # GW / doubling flavor for k=9 (accelerate the top element)
    if k >= 3:
        base = list(range(1, k))
        C['doubling top'] = tuple(base + [2 * (k - 1)])
    # tight-AP flavor {1..k} already = consec; dilations are invariant so omitted
    return C

def Ndist(E):
    nz = [abs(e) for e in E if e]; has0 = len(nz) < len(E)
    L = 1
    for e in nz: L = lcm(L, e)
    D = 7 * L; pts = set([0, D])
    for e in nz:
        st = L // e; pts.update(range(0, D + 1, st))
    pts = sorted(pts); pn = [0] * 8; b0 = 1 if has0 else 0
    for t1, t2 in zip(pts, pts[1:]):
        s = t1 + t2; hit = b0
        for e in nz: hit |= 1 << ((7 * e * s // (2 * D)) % 7)
        pn[7 - bin(hit).count("1")] += t2 - t1
    return [F(x, D) for x in pn]

# standard functionals (as functions of the core)
def J(E):
    p = Ndist(E); return sum(F(j * (7 - j)) * p[j] for j in range(7))
def POS(E):
    p = Ndist(E); T = [sum(p[j:]) for j in range(8)]; return 6*T[1]+4*T[2]+2*T[3]
def BUNCH(E):
    p = Ndist(E); T = [sum(p[j:]) for j in range(8)]; return 2*T[5]+4*T[6]
def nu(E):
    p = Ndist(E); return 1 - p[0]

def evaluate(func, k, name=""):
    C = candidates(k)
    rows = sorted(((func(E), lbl, E) for lbl, E in C.items()), key=lambda r: float(r[0]))
    print(f"  {name} on the {len(C)} algebraic candidates (k={k}):")
    for v, lbl, E in rows:
        print(f"    {lbl:18s} = {str(v):>16} ~ {float(v):.5f}")
    print(f"    ARGMIN: {rows[0][1]}   ARGMAX: {rows[-1][1]}")
    return rows

if __name__ == "__main__":
    print("VERIFY the evaluator reproduces the known extremals (k=9):")
    print("J (base functional; min should be consec):")
    evaluate(J, 9, "J")
    print("POS (min should be consec):")
    evaluate(POS, 9, "POS")
    print("BUNCH (max should be mod-7 res1 = the pole box search missed):")
    evaluate(BUNCH, 9, "BUNCH")
    print("nu = covering good-set measure (min should be consec):")
    evaluate(nu, 9, "nu")
