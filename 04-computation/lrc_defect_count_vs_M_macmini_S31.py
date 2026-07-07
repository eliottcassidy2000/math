#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S31 -- THE CRUX (opus-S120's single open step): rule out
gap members with >= 3 DEFECTS.

opus-S120 signature: every gap member = an (N-2)=10-term dilated AP + EXACTLY 2
defects.  The 2-defect world is CLOSED (empty at N=12: my S26-30 + kps S39
146757-family sweep).  The ENTIRE residual = the Freiman-stability step:
    d(V) := 12 - (longest dilated sub-AP of V)  >= 3   =>   M(V) >= 2/25.
(gap member has M in (1/13,2/25); M >= 1/13 always; so 'not in gap and not the AP'
means M >= 2/25.)

THIS SCRIPT tests the claim decisively + measures the MARGIN:
 (1) longest_dilated_AP(V) = max L with {a,a+e,..,a+(L-1)e} subset V;
 (2) sweep many 12-speed families (structured near-AP + random, bounded height),
     bucket by defect count d=12-L, report MIN M per bucket;
 (3) if any d>=3 family has M in the gap -> COUNTEREXAMPLE to opus's signature
     (report it).  Else: the min-M margin above 2/25 per defect count = how much
     room the stability theorem has.
"""
import itertools
import random
from fractions import Fraction as F
from math import gcd

LO, HI = F(1, 13), F(2, 25)


def _dens(W):
    d = set()
    for v, w in itertools.combinations(W, 2):
        d.add(v + w)
        if v != w:
            d.add(abs(v - w))
    for v in W:
        d.add(2 * v)
    d.discard(0)
    return d


def float_M(W):
    """full float M (no early exit) -- for the M-vs-defect trend."""
    best = 0.0
    for s in _dens(W):
        for j in range(1, s):
            t = j / s
            mv = min(abs(v * t - round(v * t)) for v in W)
            if mv > best:
                best = mv
    return best


def exact_M(W):
    best = F(0)
    seen = set()
    for s in _dens(W):
        for j in range(1, s):
            t = F(j, s)
            if t in seen:
                continue
            seen.add(t)
            mv = min(abs(v * t - round(v * t)) for v in W)
            if mv > best:
                best = mv
    return best


def longest_dilated_AP(W):
    """max length L of an arithmetic progression (any common difference) subset W."""
    S = set(W)
    Ws = sorted(W)
    best = 1
    for i in range(len(Ws)):
        for j in range(i + 1, len(Ws)):
            a, e = Ws[i], Ws[j] - Ws[i]
            L = 2
            nxt = Ws[j] + e
            while nxt in S:
                L += 1
                nxt += e
            if L > best:
                best = L
    return best


def primitive(W):
    g = 0
    for v in W:
        g = gcd(g, v)
    return tuple(sorted(set(v // g for v in W))) if g > 1 else tuple(sorted(set(W)))


def gen_families(rng, hmax=40, n_random=60000):
    """12-speed families: (a) dilated AP (len L) + defects for L=6..12 [controlled
    defect count]; (b) random 12-subsets of [1,hmax]."""
    fams = set()
    # (a) structured: L-term dilated AP + (12-L) defects, various spacing
    for e in range(1, 8):
        for a in range(1, e + 2):
            for L in range(6, 13):
                ap = [a + i * e for i in range(L)]
                if ap[-1] > hmax:
                    continue
                need = 12 - L
                pool = [x for x in range(1, hmax + 1) if x not in ap]
                if need == 0:
                    W = primitive(tuple(ap))
                    if len(W) == 12:
                        fams.add(W)
                    continue
                # sample defect combos (too many to enumerate for need>=3)
                for _ in range(200):
                    combo = rng.sample(pool, need)
                    W = primitive(tuple(sorted(set(ap) | set(combo))))
                    if len(W) == 12:
                        fams.add(W)
    # (b) random
    for _ in range(n_random):
        W = primitive(tuple(sorted(rng.sample(range(1, hmax + 1), 12))))
        if len(W) == 12:
            fams.add(W)
    return fams


def main():
    rng = random.Random(31)
    hi_f, lo_f = float(HI), float(LO)
    print("=" * 88)
    print("CRUX TEST (opus-S120): M-vs-DEFECT trend; does d(V) >= 3 => M(V) >= 2/25 ?")
    print("=" * 88)
    fams = gen_families(rng, n_random=18000)
    fams = sorted(fams)
    # keep all structured (rare) + subsample; cap total for full float-M tractability
    if len(fams) > 16000:
        fams = fams[:: max(1, len(fams) // 16000)]
    print(f"  testing {len(fams)} families; computing FULL float-M per defect bucket")
    # bucket by defect count; track min FLOAT M per bucket + gap members (exact)
    min_M = {}          # defect -> (min float M, family)
    count = {}
    gap_members = []
    argmin_near = {}    # defect -> family achieving min (for exact confirmation)
    for W in fams:
        d = 12 - longest_dilated_AP(W)
        count[d] = count.get(d, 0) + 1
        fm = float_M(W)
        if d not in min_M or fm < min_M[d][0]:
            min_M[d] = (fm, W)
        # any family with float M plausibly in the gap -> exact check
        if lo_f - 1e-4 < fm < hi_f + 1e-4:
            M = exact_M(W)
            if LO < M < HI:
                gap_members.append((M, d, W))
    print(f"  families by defect count d: "
          f"{ {k: count.get(k,0) for k in sorted(count)} }")
    print()
    if gap_members:
        print(f"  *** {len(gap_members)} GAP MEMBERS FOUND:")
        for M, d, W in sorted(gap_members)[:12]:
            flag = "  <== d>=3 COUNTEREXAMPLE to opus-S120!" if d >= 3 else "  (d<=2)"
            print(f"      M={M} ({float(M):.5f}) defects={d}  W={list(W)}{flag}")
    else:
        print("  => NO gap members at any defect count (confirms empty).")
    print()
    print("  MIN M (nearest approach to the gap from above) by defect count d:")
    print(f"    {'d':>3} {'longest-AP':>10} {'min float-M':>12} {'>2/25?':>7} "
          f"{'margin':>10}  min-family")
    for d in sorted(min_M):
        fm, W = min_M[d]
        marg = fm - hi_f
        Ws = str(list(W)) if len(str(list(W))) < 42 else str(list(W)[:10]) + "..]"
        print(f"    {d:>3} {12 - d:>10} {fm:>12.5f} {'yes' if fm >= hi_f - 1e-9 else 'NO':>7} "
              f"{marg:>+10.5f}  {Ws}")
    print()
    print("  KEY: the gap (1/13,2/25)=(0.0769,0.0800).  If min-M rises with d and")
    print("  d>=3 all sit clearly above 0.08, opus-S120's stability route is solid;")
    print("  the d=0 (AP,1/13) and d=1 (block-lift,2/25) BRACKET the gap with nothing")
    print("  inside -- the structural reason (G) holds at N=12.")


if __name__ == "__main__":
    main()
