#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S17 -- HYP-4452: THE DENSITY FLOOR via multi-lens synthesis.

The sole remaining (G) obligation (S16b): safe(2/25) > 0 for non-AP, uniform in
height, where safe(rho) = Leb{t : ||v_i t|| >= rho for all i}.  REFRAME:
safe > 0 <=> the 12 danger arcs (radius rho) do NOT cover the circle <=> a
COVERING/TILING statement -- the AP is the unique tiler (safe=0); non-AP leaves
a gap.

This probes, exactly (rational arithmetic on arc endpoints):
 (1) safe(2/25) for the AP (=0) and non-AP families (>0?) -- the floor object;
 (2) the MINIMUM of safe over non-AP -- is it bounded below (the floor)?
 (3) THE SCALE-GAP THRESHOLD R_0: families with a ratio-R gap -> decorrelation
     -> M >= 2/25 (safe > 0).  R_0 => height bound max < R_0^11 => finite check;
 (4) n-DEPENDENCE (S15): safe structure at n=7 (floor FAILS, gap members) vs 13.
"""
import itertools, random
from fractions import Fraction as F
from math import gcd

def safe_measure(W, rho):
    """EXACT Leb{t in [0,1): ||v t|| >= rho for all v in W}, via arc endpoints.
    danger for comb v: t within rho/v of m/v.  safe = complement."""
    # collect danger intervals (as fractions) mod 1
    ivals = []
    for v in W:
        r = F(rho) / v
        for m in range(v):
            c = F(m, v)
            ivals.append((c - r, c + r))
    # normalize endpoints into [0,1) breakpoints
    pts = set([F(0), F(1)])
    for lo, hi in ivals:
        for x in (lo, hi):
            xr = x - int(x)
            if xr < 0:
                xr += 1
            pts.add(xr)
    spts = sorted(pts)
    def in_danger(t):
        for v in W:
            y = v * t
            y = y - int(y)
            d = min(y, 1 - y)
            if d < rho:
                return True
        return False
    safe = F(0)
    for i in range(len(spts) - 1):
        a, b = spts[i], spts[i + 1]
        mid = (a + b) / 2
        if not in_danger(mid):
            safe += (b - a)
    return safe

def primitive(W):
    g = 0
    for v in W:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in W)) if g > 1 else tuple(sorted(W))

RHO = F(2, 25)

def part1_safe():
    print("=" * 78)
    print("PART 1: safe(2/25) for AP vs non-AP (the density-floor object)")
    print("=" * 78)
    AP = tuple(range(1, 13))
    print(f"  AP {{1..12}}: safe(2/25) = {safe_measure(AP, RHO)} "
          f"(expect 0: M(AP)=1/13 < 2/25, arcs TILE)")
    random.seed(17)
    mins = []
    zero_nonAP = []
    for _ in range(3000):
        W = primitive(tuple(sorted(random.sample(range(1, 40), 12))))
        if len(set(W)) != 12:
            continue
        s = safe_measure(W, RHO)
        mins.append((s, W))
        if s == 0:
            zero_nonAP.append(W)
    mins.sort()
    print(f"  non-AP families sampled: {len(mins)}")
    print(f"  MIN safe(2/25) over non-AP: {mins[0][0]} = {float(mins[0][0]):.5f} "
          f"at {list(mins[0][1])}")
    print(f"  non-AP families with safe = 0 (would be gap/tight): {len(zero_nonAP)}")
    for W in zero_nonAP[:5]:
        print(f"    safe=0: {list(W)}")
    print(f"  smallest 5 safe values (the floor's approach):")
    for s, W in mins[:5]:
        print(f"    safe={s}={float(s):.5f}  {list(W)}")

def part3_scalegap():
    print()
    print("=" * 78)
    print("PART 3: THE SCALE-GAP THRESHOLD R_0 (the height bound)")
    print("(a family with a ratio-R scale gap decorrelates => M >= 2/25 => safe>0)")
    print("=" * 78)
    import sys; sys.path.insert(0, '04-computation')
    # small part A (<=6 elts) + large part B (>= R*maxA), check safe>0 for growing R
    random.seed(171)
    print(f"  {'ratio R':>8} {'#tested':>8} {'safe=0 (no decorr)':>20} {'min safe':>12}")
    for Rband in [(2, 3), (3, 5), (5, 8), (8, 13), (13, 21), (21, 40)]:
        tested = 0; zeros = 0; minsafe = None
        for _ in range(400):
            na = random.randint(4, 8)
            A = random.sample(range(1, 13), na)
            maxA = max(A)
            nb = 12 - na
            lo = Rband[0] * maxA; hi = Rband[1] * maxA
            if hi - lo < nb:
                continue
            B = random.sample(range(lo, hi + 1), nb)
            W = primitive(tuple(sorted(set(A + B))))
            if len(set(W)) != 12:
                continue
            tested += 1
            s = safe_measure(W, RHO)
            if s == 0:
                zeros += 1
            if minsafe is None or s < minsafe:
                minsafe = s
        print(f"  {str(Rband):>8} {tested:>8} {zeros:>20} "
              f"{(float(minsafe) if minsafe is not None else 0):>12.5f}")
    print("  => R_0 = the ratio where safe=0 STOPS occurring (decorrelation guaranteed).")
    print("     max element < R_0^11 => (G) is a FINITE CHECK (with the S16 q<=2max lever).")

def part4_ndep():
    print()
    print("=" * 78)
    print("PART 4: n-DEPENDENCE -- safe at n=7 (floor FAILS) vs the structure")
    print("=" * 78)
    # n=7 gap member has safe(2/13) = 0 (M=5/33 < 2/13)
    m7 = (1, 5, 6, 11, 16, 17)
    print(f"  n=7 gap member {list(m7)}: safe(2/13) = {safe_measure(m7, F(2,13))} "
          f"(=0: it's a gap member, floor FAILS at n=7)")
    # a non-gap n=7 family for contrast
    print(f"  n=7 AP {{1..6}}: safe(2/13) = {safe_measure(tuple(range(1,7)), F(2,13))} "
          f"(=0: M(AP6)=1/7 < 2/13, the tiler)")
    print(f"  n=7 loose (1,2,4,8,9,11): safe(2/13) = "
          f"{safe_measure((1,2,4,8,9,11), F(2,13))}")
    print("  => at n=7 both the AP AND the gap member have safe=0; the floor is")
    print("     the statement that at n=13, ONLY the AP does (the two walls, S15).")

if __name__ == "__main__":
    part1_safe()
    part3_scalegap()
    part4_ndep()
