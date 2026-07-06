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
    # (1) MIN-WITNESS-DENOM for ALL families (fast, height-independent) -- the
    #     bound question.  None => no witness up to Qcap (M < 2/25).
    clearing = []      # (name, W, minq)
    nowit = []
    Qcap = 130
    for nm, W in fams:
        mq = min_witness_denom(W, Qcap=Qcap)
        if mq is None:
            nowit.append((nm, W))
        else:
            clearing.append((nm, W, mq))
    qs = [c[2] for c in clearing]
    mx = max(qs); argmax = [c for c in clearing if c[2] == mx][0]
    print(f"  families: {len(fams)};  with a 2/25-witness (M>=2/25): {len(clearing)};  "
          f"no witness up to q={Qcap} (M<2/25): {len(nowit)}")
    print(f"  MIN-WITNESS-DENOM: min={min(qs)}, median={sorted(qs)[len(qs)//2]}, MAX={mx}")
    print(f"    largest min-witness q: {argmax[0]} q={mx}  W={argmax[1]}")
    over50 = [c for c in clearing if c[2] > 50]
    print(f"  families needing q > 50: {len(over50)}  "
          f"({'*** Q50 BOUND VIOLATED (Q_max=%d)'%mx if over50 else 'Q50 bound holds on this sample'})")
    for c in over50[:8]:
        print(f"    q={c[2]}: {c[0]}  W={c[1]}")
    # (2) SPECTRUM in [2/25, 1/6] + (3) 2/25-exact maximizer -- exact M only on
    #     SMALL-height families (speeds <= 80: exact_M grid is cheap there).
    print()
    print("  SPECTRUM in [2/25, 1/6] and 2/25-exact maximizers (small-height subsample):")
    spectrum_in_band = {}
    maxq_exact = []
    small = [(nm, W) for nm, W in fams if max(W) <= 80]
    for nm, W in small:
        M, t = exact_M_argmax(W)
        if F(1,13) < M < F(2,25):
            print(f"    *** GAP HIT: {nm} M={M} W={W}")
        if F(2,25) <= M <= F(1,6):
            spectrum_in_band[M] = spectrum_in_band.get(M, 0) + 1
        if M == F(2,25):
            maxq_exact.append((nm, t.denominator))
    for M in sorted(spectrum_in_band):
        print(f"    {M} = {float(M):.5f}  (denom {M.denominator})  x{spectrum_in_band[M]}")
    dens = [M.denominator for M in spectrum_in_band]
    if dens:
        print(f"  => band denominators: max = {max(dens)} "
              f"({'all <= 50 (Farey-small: the bound reason)' if max(dens)<=50 else 'some > 50!'})")
    if maxq_exact:
        print(f"  M=2/25-exact maximizer q*: max = {max(q for _,q in maxq_exact)} "
              f"(over {len(maxq_exact)} families)")

if __name__ == "__main__":
    part_A()
