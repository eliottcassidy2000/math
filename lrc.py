from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
import random

def frac_norm(x):
    # ||x|| distance to nearest integer, x a Fraction
    r = x - (x.numerator // x.denominator)  # in [0,1)
    # r in [0,1)
    return r if r <= F(1,2) else 1 - r

def candidates(vs):
    """Return set of candidate t in [0,1) where D(t) max can occur."""
    C = set()
    vs = list(vs)
    # tent peaks: t=(2k+1)/(2 v_i)
    for v in vs:
        d = 2*v
        for k in range(0, d):
            t = F(2*k+1, d)
            if 0 <= t < 1:
                C.add(t)
    # tent crossings: t = k/(v_i+v_j) and k/(v_i-v_j)
    for i in range(len(vs)):
        for j in range(len(vs)):
            if i == j:
                continue
            for s in (vs[i]+vs[j], abs(vs[i]-vs[j])):
                if s == 0:
                    continue
                for k in range(0, s):
                    C.add(F(k, s))
    return C

def D(t, vs):
    return min(frac_norm(v*t) for v in vs)

def M_exact(vs):
    """Exact max over t of D(t). Returns (M, argmax_t)."""
    C = candidates(vs)
    best = None
    bestt = None
    for t in C:
        d = D(t, vs)
        if best is None or d > best:
            best = d
            bestt = t
    return best, bestt

def verify_with_walls(vs, claimed_M):
    """Independent check: D is piecewise linear; all breakpoints are tent
    vertices (t=k/(2 v_i)) plus all pairwise crossings. The max over [0,1)
    is attained at a breakpoint. Compute max over the full breakpoint set,
    and also scan cell midpoints as a sanity check that no interior point
    exceeds the breakpoint max."""
    bp = set()
    bp.add(F(0)); bp.add(F(1))
    for v in vs:
        d = 2*v
        for k in range(0, d+1):
            bp.add(F(k, d))
    # crossings included as candidates already; add them explicitly
    bp |= candidates(vs)
    bp = sorted(b for b in bp if 0 <= b <= 1)
    mx = F(0)
    for w in bp:
        if 0 <= w < 1:
            d = D(w, vs)
            if d > mx:
                mx = d
    # midpoint scan: interior of a linear cell can only be <= endpoints
    for a, b in zip(bp, bp[1:]):
        mid = (a+b)/2
        d = D(mid, vs)
        if d > mx:
            mx = d
    return mx

def primitive(vs):
    return reduce(gcd, vs) == 1
