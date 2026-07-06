#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S7 -- HYP-4312 Parts A2 + B.

PART A2: THE DELICATE CASE of the Q50 census bound.  By the uniform cell lemma
(HYP-4252), a family with M = 2/25 has maximizer q* = 25k (k = q*/25).  Because
a 2/25-attainer HAS a 2/25-point, kps's cluster-gcd ladder (which requires a
no-2/25-point family) does NOT bound k -- so a priori q* could be 50, 75, ....
A 2/25-attainer with q* > 50 would appear to opus-S98's q <= 50 census as NOT
clearing (its only 2/25-witness lives beyond the census cap) => MISCLASSIFIED as
a gap family => the census bound FAILS.

So the census bound q <= 50 is SOUND iff every primitive 2/25-attainer has
q* <= 50 (i.e. k <= 2).  This searches hard for large-q* 2/25-attainers.

PART B: the covering-floor discrepancy.  My S5 adversarial search found
consec[1..11] TILES at radius 2/25 (phi_worst -> 0); kps-S20d reports the
distinct-freq floor >= 0.06 through l=14.  EXACT rational resolution: does the
union of radius-2/25 combs of frequencies [1..L] cover [0,1) at some phase?
"""
import itertools, random
from fractions import Fraction as F
from math import gcd

def dist_int(x, q):
    r = x % q
    return min(r, q - r)

def M_at(W, a, q):
    """margin of W at t = a/q, as a Fraction."""
    return min(F(dist_int(v * a, q), q) for v in W)

def exact_M_and_all_maximizers(W):
    """M(W) and ALL reduced maximizer points a/s (for q* analysis)."""
    dens = set()
    for v, w in itertools.combinations(W, 2):
        dens.add(v + w)
        if v != w:
            dens.add(abs(v - w))
    for v in W:
        dens.add(2 * v)
    best = F(0)
    pts = []
    seen = set()
    for s in dens:
        if s == 0:
            continue
        for j in range(1, s):
            t = F(j, s)
            if t in seen:
                continue
            seen.add(t)
            mv = min(abs(v * t - round(v * t)) for v in W)
            if mv > best:
                best, pts = mv, [t]
            elif mv == best:
                pts.append(t)
    return best, pts

def primitive(W):
    g = 0
    for v in W:
        g = gcd(g, v)
    return tuple(v // g for v in W) if g > 1 else tuple(W)

def part_A2():
    print("=" * 78)
    print("PART A2: hunting large-q* 2/25-attainers (the census bound's failure mode)")
    print("=" * 78)
    RHO = F(2, 25)
    random.seed(77)
    found = {}   # q* -> list of families
    tested = 0
    # (i) the known attainer + its residue-structured relatives (k=1, q*=25)
    known = [tuple([1,2,3,5,7,8,9,10,11,12,17,19])]
    # (ii) 2-LIFT candidates targeting k=2 (q*=50): residues mod 50 with the
    #      block-lift structure doubled.  Base {1..12}\{4,6} then lift 4,6 by
    #      multiples of 25 (the cell-50 grid): 4->4+25=29, 6->6+25=31 (k=1 on 25)
    #      or engineer mod-50 structure.
    cand = []
    base10 = [1,2,3,5,7,8,9,10,11,12]
    # single/double lifts by multiples that could push q* up
    for L1 in range(13, 60):
        for L2 in range(L1+1, 70):
            cand.append(tuple(sorted(base10 + [L1, L2])))
    # also: dilations/near-attainer perturbations
    for _ in range(4000):
        W = sorted(random.sample(range(1, 55), 12))
        cand.append(primitive(W))
    # residue-mod-25 engineered (unit pairs summing to 25 = the cell-25 shape)
    for _ in range(3000):
        # pick a unit pair (a,25-a) then fill with clearers of level 2 mod 25
        a = random.choice([1,2,3,4,6,7,8,9,11,12])
        pair = [a, 25 - a]
        m = random.choice([x for x in range(1,25) if gcd(x,25)==1])
        pool = [v for v in range(1, 55) if dist_int(v*m,25) >= 2 and v not in pair]
        if len(pool) >= 10:
            W = tuple(sorted(pair + random.sample(pool, 10)))
            cand.append(primitive(W))
    seen = set()
    for W in known + cand:
        if len(set(W)) != 12 or W in seen:
            continue
        seen.add(W)
        tested += 1
        M, pts = exact_M_and_all_maximizers(W)
        if M == RHO:
            qstar = min(t.denominator for t in pts)  # smallest maximizer denom
            found.setdefault(qstar, []).append(W)
    print(f"  tested {tested} distinct families")
    if found:
        print(f"  2/25-EXACT families found, by maximizer q*:")
        for qs in sorted(found):
            k = qs // 25 if qs % 25 == 0 else F(qs, 25)
            print(f"    q* = {qs} (k={k}): {len(found[qs])} families; e.g. {found[qs][0]}")
        maxq = max(found)
        print(f"  => MAX maximizer q* among 2/25-attainers = {maxq} "
              f"({'<= 50: opus Q50 census SOUND for the boundary case' if maxq<=50 else '*** > 50: census cap must rise to '+str(maxq)})")
    else:
        print("  no 2/25-exact families found in this search "
              "(the attainer species is sparse; try the known one directly)")
        M, pts = exact_M_and_all_maximizers(known[0])
        print(f"  known attainer M={M}, maximizer denoms {sorted(set(t.denominator for t in pts))}")

def covers_exactly(freqs, phases, rho):
    """EXACT: does the union of combs {f s + phi : dist < rho} cover [0,1)?
    Check every gap between consecutive danger-arc endpoints.  All rational."""
    # danger arcs for comb (f, phi): dist(f s + phi) < rho
    #   <=> f s + phi in (m - rho, m + rho) for some integer m
    #   <=> s in ((m - rho - phi)/f, (m + rho - phi)/f)
    intervals = []
    for f, phi in zip(freqs, phases):
        # m ranges so that the interval meets [0,1): f*s+phi in (m-rho, m+rho)
        # s in [0,1) => f*s+phi in [phi, f+phi)
        mlo = -2
        mhi = f + 2
        for m in range(int(mlo), int(mhi) + 1):
            lo = (F(m) - rho - phi) / f
            hi = (F(m) + rho - phi) / f
            intervals.append((lo, hi))
    # normalize onto [0,1): collect breakpoints, test midpoint of each gap
    pts = set()
    for lo, hi in intervals:
        for x in (lo, hi):
            xr = x - int(x)          # frac part in [0,1)
            if xr < 0: xr += 1
            pts.add(xr)
    pts.add(F(0)); pts.add(F(1))
    spts = sorted(pts)
    def covered(s):
        for f, phi in zip(freqs, phases):
            y = f * s + phi
            yr = y - round(y)
            if abs(yr) < rho:
                return True
        return False
    # test the midpoint of every elementary interval
    for i in range(len(spts) - 1):
        mid = (spts[i] + spts[i+1]) / 2
        if not covered(mid):
            return False, mid
    return True, None

def part_B():
    print()
    print("=" * 78)
    print("PART B: does distinct-freq {1..L} cover the circle at radius 2/25? (EXACT)")
    print("=" * 78)
    rho = F(2, 25)
    random.seed(5)
    for L in range(7, 13):
        freqs = list(range(1, L + 1))
        # search for a covering rational phase (denominators up to D)
        covered_found = None
        # try structured phases first: k/(2L) offsets, then random rationals
        Dtry = 2 * max(freqs) * 5
        trials = 0
        best_gap = None
        for _ in range(3000):
            phases = [F(random.randint(0, Dtry - 1), Dtry) for _ in freqs]
            cov, gap = covers_exactly(freqs, phases, rho)
            trials += 1
            if cov:
                covered_found = phases
                break
        status = "COVERS (a rational phase tiles!)" if covered_found else f"no cover found in {trials} rational phases"
        print(f"  L={L} freqs {freqs}: {status}")
        if covered_found:
            print(f"    covering phase (x{Dtry}): {[p*Dtry for p in covered_found]}")

if __name__ == "__main__":
    part_A2()
    part_B()
