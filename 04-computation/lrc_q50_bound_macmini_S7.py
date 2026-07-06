#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S7 -- HYP-4312 Part A: HARDENING THE Q50 CENSUS BOUND.

opus-S98's census reduction (LRCRayTransport.margin_of_residue_witness) makes
the spectral gap a FINITE census IF every gap-clearing 12-family has a
margin-2/25 witness at bounded denominator q <= Q_max.  opus supports Q_max = 50
by 300/300 sampling.  This script hardens that bound with EXACT rational
arithmetic and locates its structural reason.

THE THEORETICAL CLAIM (to test): the spectrum just above 2/25 consists of
SMALL-DENOMINATOR Farey fractions (2/23, 1/11, 3/34, 2/21, ...), so a family
with M(v) >= 2/25 clears at a small-q point.  The one delicate case is
M(v) = 2/25 EXACTLY: the safe set is a single maximizer point at q* -- is q*
bounded?

MEASURES (exact, via dist_int):
 (1) MIN WITNESS DENOM: for M(v) >= 2/25, the least q with some a/q having
     margin >= 2/25.  Distribution over many families; max = the empirical Q_max.
 (2) SPECTRUM DENOMINATORS in [2/25, 1/6]: enumerate achieved M-values, confirm
     reduced denominators are small (the Farey reason).
 (3) THE M=2/25-EXACT BOUNDARY: primitive families with M exactly 2/25 -- their
     maximizer denominator q* (the delicate case for the census).
 (4) HIGH-HEIGHT via the residue bridge: a height-1e6 family's witness = its
     residue family's witness (opus-S98) -- confirm min witness denom is a
     residue-class property (bounded independent of height).
"""
import itertools, random
from fractions import Fraction as F
from math import gcd

def dist_int(x, q):
    r = x % q
    return min(r, q - r)

def exact_M_argmax(W):
    """M(W) = max_t min_i ||v_i t|| and one maximizer a/s, via the critical grid
    (sums/diffs/half-integers -- THM-592 grid attainment)."""
    dens = set()
    for v, w in itertools.combinations(W, 2):
        dens.add(v + w)
        if v != w:
            dens.add(abs(v - w))
    for v in W:
        dens.add(2 * v)
    best = F(0); best_t = F(0)
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
                best, best_t = mv, t
    return best, best_t

RHO = F(2, 25)

def min_witness_denom(W, Qcap=400):
    """least q with some a in [1,q-1] s.t. margin(W, a/q) >= 2/25 (exact).
    margin >= 2/25 <=> dist_int(v_i a, q) >= 2q/25 for all i."""
    for q in range(1, Qcap + 1):
        thr = F(2, 25) * q      # need dist_int(v_i a, q) >= thr (as rationals)
        for a in range(1, q):
            ok = True
            for v in W:
                if dist_int(v * a, q) < thr:
                    ok = False
                    break
            if ok:
                return q
    return None  # no witness up to Qcap (M < 2/25 or witness beyond cap)

def primitive(W):
    g = 0
    for v in W:
        g = gcd(g, v)
    return [v // g for v in W] if g > 1 else list(W)

def make_families(seed=7):
    random.seed(seed)
    fams = []
    # the exact 2/25 attainer + relatives
    fams.append(("attainer 2/25", [1,2,3,5,7,8,9,10,11,12,17,19]))
    # deep well (n=13 extremizer analog), lifts
    fams.append(("deep-well {1..11,168}", list(range(1,12))+[168]))
    fams.append(("14r-ladder r=7", [1,2,3,4,5,6,8,9,10,11,12,98]))
    # double-lift species (S4): near 2/17
    base = [1,2,3,5,7,8,9,10,11,12]
    for s in (1,2,3,17):
        fams.append((f"dbl-lift s={s}", base+[4+13*s,6+13*s]))
    # random primitive families of various heights
    for _ in range(400):
        n = 12
        W = sorted(random.sample(range(1, 80), n))
        fams.append(("rand", primitive(W)))
    # near-gap engineered: families with M just above 2/25 (small-denom Farey)
    # 2/23-attainers: lift structure targeting q=23
    for _ in range(60):
        W = sorted(random.sample(range(1,50), 12))
        fams.append(("rand-lowH", primitive(W)))
    # high-height families (test the bridge): add a big multiple
    for _ in range(60):
        W = sorted(random.sample(range(1,40), 11))
        W.append(random.randint(1000, 5000))
        fams.append(("high-height", primitive(sorted(set(W)))))
    return [(nm, W) for nm, W in fams if len(set(W)) == 12]

def part_A():
    print("=" * 78)
    print("PART A: the Q50 witness-denominator bound (exact)")
    print("=" * 78)
    fams = make_families()
    clearing = []      # (name, W, M, minq)
    gap_hits = []
    spectrum_in_band = {}   # M-value -> count, for M in [2/25, 1/6]
    maxq_exact = []    # M == 2/25 exactly: maximizer denom
    for nm, W in fams:
        M, t = exact_M_argmax(W)
        if F(1,13) < M < F(2,25):
            gap_hits.append((nm, W, M, t))
            continue
        if M >= F(2,25):
            mq = min_witness_denom(W)
            clearing.append((nm, W, M, mq))
            if F(2,25) <= M <= F(1,6):
                spectrum_in_band[M] = spectrum_in_band.get(M, 0) + 1
            if M == F(2,25):
                maxq_exact.append((nm, W, t.denominator))
    # report
    print(f"  families: {len(fams)};  gap hits (M in (1/13,2/25)): {len(gap_hits)}")
    if gap_hits:
        for nm,W,M,t in gap_hits[:5]:
            print(f"    *** GAP: {nm} M={M} W={W}")
    valid = [c for c in clearing if c[3] is not None]
    if valid:
        qs = [c[3] for c in valid]
        mx = max(qs)
        argmax = [c for c in valid if c[3] == mx][0]
        print(f"  clearing families (M >= 2/25): {len(valid)}")
        print(f"  MIN-WITNESS-DENOM distribution: min={min(qs)}, median={sorted(qs)[len(qs)//2]}, MAX={mx}")
        print(f"    argmax (largest min-witness q): {argmax[0]} M={argmax[2]} -> q={mx}")
        over50 = [c for c in valid if c[3] > 50]
        print(f"  families needing q > 50: {len(over50)}  "
              f"({'*** Q50 BOUND VIOLATED' if over50 else 'Q50 bound holds on this sample'})")
        for c in over50[:6]:
            print(f"    q={c[3]}: {c[0]} M={c[2]} W={c[1]}")
    none_q = [c for c in clearing if c[3] is None]
    if none_q:
        print(f"  families with no witness up to Qcap=400: {len(none_q)} (M>=2/25 but witness beyond cap?)")
        for c in none_q[:4]:
            print(f"    {c[0]} M={c[2]} W={c[1]}")
    print()
    print("  SPECTRUM VALUES in [2/25, 1/6] (the band just above the gap):")
    for M in sorted(spectrum_in_band):
        print(f"    {M} = {float(M):.5f}  (denom {M.denominator})  x{spectrum_in_band[M]}")
    dens = [M.denominator for M in spectrum_in_band]
    if dens:
        print(f"  => band denominators: max = {max(dens)} "
              f"({'all <= 50: Farey-small, the bound reason' if max(dens)<=50 else 'some large!'})")
    print()
    print("  M = 2/25 EXACTLY (the delicate boundary case) maximizer denominators:")
    for nm,W,q in maxq_exact[:10]:
        print(f"    {nm}: maximizer q* = {q}")
    if maxq_exact:
        print(f"  => max maximizer q* among 2/25-exact = {max(q for _,_,q in maxq_exact)}")

if __name__ == "__main__":
    part_A()
