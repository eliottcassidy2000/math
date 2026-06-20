#!/usr/bin/env python3
"""
INDEPENDENT verification (written from scratch) of the fully-decorrelated
boundary value claim:

  P_r(B) = sum_{t=0}^6 prof_t(B) * c_t(r)
  prof_t(B) = meas{ x in [0,1) : B misses exactly t of the 7 sectors }
  c_t(r)    = sum_{i=0}^t (-1)^i C(t,i) (1 - i/7)^r

Claim to test (try to REFUTE):
  For dangerous rows k=8,9,10 with caps cap_8=2243/5880, cap_9=1979/4004,
  cap_10=4/7, we have
       max over B subset {0..14}, |B|=k-r, 0 in B   of  P_r(B)  <=  cap_k
  with margin >= ~0.17.
  Specifically (k,r)=(9,2) expect max ~0.246 (consec_7),
                 (k,r)=(10,2) expect max ~0.399 (consec_8).

Also verify P_r(B) is the actual limit of p0(B u F) for a dissociated far
set F (e.g. F=[23,41] for r=2).

Everything exact (Fraction) where possible. Independent re-derivations:
 - prof computed two independent ways (interval scan AND rational-point sampling
   on a fine common grid) to cross-check.
 - p0 (measure of full-coverage set) computed independently from prof of B u F.
"""
from fractions import Fraction as F
from math import comb
from itertools import combinations
import sys
sys.stdout.reconfigure(line_buffering=True)

DEN = 7  # number of sectors

# ----------------------------------------------------------------------
# Sector of e*x for x in [0,1): floor( (e*x mod 1) * 7 ).
# For x a rational p/q this is exact.
# ----------------------------------------------------------------------
def sector(e, x):
    """sector index 0..6 of e*x, x a Fraction in [0,1)."""
    frac = (e * x) % 1          # Fraction in [0,1)
    return int(frac * DEN)      # floor since frac*7 >= 0; int() truncates toward 0

def breakpoints(E):
    """All x in [0,1] where any sector(e*.) can jump: e*x = m/7 => x = m/(7e)."""
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for m in range(0, 7 * e + 1):
            bps.add(F(m, 7 * e))
    return sorted(bps)

def distinct_sectors_at(E, x):
    return {sector(e, x) for e in E}

def miss_profile(B):
    """
    prof[t] = meas{ x in [0,1) : number of sectors B FAILS to hit == t }, t=0..6.
    A constant-sector cell is [x0,x1); evaluate at midpoint.
    """
    B = sorted(set(B))
    bps = breakpoints(B)
    prof = [F(0)] * 7
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        hit = len(distinct_sectors_at(B, xm))
        t = 7 - hit
        assert 0 <= t <= 6, (B, t, hit)
        prof[t] += (x1 - x0)
    return prof

def p0_fullcover(E):
    """meas{ x : the set E hits ALL 7 sectors } = prof_0 of E."""
    return miss_profile(E)[0]

# Independent c_t(r) as exact Fraction (inclusion-exclusion over which of the
# t sectors stay uncovered by r iid uniform-on-7 runners).
def c_t(t, r):
    s = F(0)
    for i in range(t + 1):
        s += F((-1) ** i) * comb(t, i) * F(7 - i, 7) ** r
    return s

def P_r(B, r):
    prof = miss_profile(B)
    return sum(prof[t] * c_t(t, r) for t in range(7))

# ----------------------------------------------------------------------
# SECOND, fully independent prof computation by exact rational sampling:
# Build common refinement grid from all breakpoints, then sum cell lengths
# but determine sector-count by the ALGEBRAICALLY-different route of building
# a fine uniform grid of denominator L = lcm(7*e). This double-checks the
# interval-scan above.
# ----------------------------------------------------------------------
from math import lcm
def miss_profile_gridcheck(B):
    B = sorted(set(B))
    nz = [e for e in B if e != 0]
    if not nz:
        # all hit sector 0 only => misses 6 sectors everywhere
        prof = [F(0)] * 7
        prof[6] = F(1)
        return prof
    L = 1
    for e in nz:
        L = lcm(L, 7 * e)
    # cells are [j/L,(j+1)/L); sector const on each since L multiple of 7e
    prof = [F(0)] * 7
    cell = F(1, L)
    for j in range(L):
        xm = F(2 * j + 1, 2 * L)
        hit = len(distinct_sectors_at(B, xm))
        prof[7 - hit] += cell
    return prof

# ----------------------------------------------------------------------
caps = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(4, 7)}

def maximize(k, r):
    """max P_r(B) over B subset {0..14}, |B|=k-r, 0 in B."""
    size = k - r
    universe = list(range(0, 15))  # {0..14}
    rest = [e for e in universe if e != 0]
    best = None
    best_B = None
    # choose size-1 elements from {1..14}, always include 0
    for combo in combinations(rest, size - 1):
        B = (0,) + combo
        val = P_r(B, r)
        if best is None or val > best:
            best = val
            best_B = B
    return best, best_B

print("=" * 78)
print("INDEPENDENT cross-check of prof_t on a few B (interval scan vs grid)")
print("=" * 78)
for B in [(0,1,2,3,4,5,6), (0,1,2,3,4,5,6,7), (0,4,6,8,10,12,14)]:
    p1 = miss_profile(B)
    p2 = miss_profile_gridcheck(B)
    ok = (p1 == p2)
    print(f"B={B}")
    print(f"   prof(scan) = {[str(x) for x in p1]}")
    print(f"   prof(grid) = {[str(x) for x in p2]}")
    print(f"   total={sum(p1)}  MATCH={ok}")
    assert ok, "prof methods disagree!"

print("\n" + "=" * 78)
print("c_t(r) exact values (r=2):")
print("=" * 78)
for t in range(7):
    print(f"   c_{t}(2) = {c_t(t,2)} = {float(c_t(t,2)):.6f}")

print("\n" + "=" * 78)
print("MAXIMIZE P_r(B) over B subset {0..14}, |B|=k-r, 0 in B")
print("=" * 78)
for (k, r) in [(8, 2), (9, 2), (10, 2)]:
    best, bestB = maximize(k, r)
    cap = caps[k]
    margin = cap - best
    status = "OK (P_r < cap)" if best < cap else "OVER CAP!"
    print(f"\n(k={k}, r={r}): size |B|={k-r}")
    print(f"   argmax B          = {bestB}")
    print(f"   max P_r(B)        = {best} = {float(best):.6f}")
    print(f"   cap_{k}            = {cap} = {float(cap):.6f}")
    print(f"   margin = cap-max  = {margin} = {float(margin):.6f}   [{status}]")

print("\n" + "=" * 78)
print("LIMIT CHECK: p0(B u F) -> P_r(B) for dissociated F=[23,41], r=2")
print("=" * 78)
Fset = [23, 41]
for B in [(0,1,2,3,4,5,6,7), (0,1,2,3,4,5,6,7,8)]:
    Pr = P_r(B, 2)
    p0 = p0_fullcover(list(B) + Fset)
    print(f"\nB={B}")
    print(f"   P_r(B)          (boundary value)     = {Pr} = {float(Pr):.6f}")
    print(f"   p0(B u [23,41]) (actual, exact)      = {p0} = {float(p0):.6f}")
    print(f"   |diff|                               = {float(abs(Pr - p0)):.6e}")

print("\n" + "=" * 78)
print("LIMIT CHECK with FARTHER dissociated set F to show convergence")
print("=" * 78)
B = (0,1,2,3,4,5,6,7)
Pr = P_r(B, 2)
print(f"B={B},  P_r(B) target = {float(Pr):.8f}")
for Fset in [[23,41],[101,211],[1009,2003],[100003,200003]]:
    p0 = p0_fullcover(list(B) + Fset)
    print(f"   F={Fset}: p0 = {float(p0):.8f}  |p0-P_r|={float(abs(p0-Pr)):.3e}")

print("\nDONE.")
