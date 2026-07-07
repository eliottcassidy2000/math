#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S30 -- the OPEN part: do the ORDER-3+ (dilated-AP) species
land in the N=12 gap?  Integrating kps-S39 + opus-S119 + my THM-632.

THE TENSION opus-S119 flagged: 'N even => first gap empty' CANNOT be universal --
N=6 is EVEN yet NONEMPTY via {1,5,6,11,16,17}=5/33 (order 3).  My THM-632/parity
covers the order-2 (canonical mediant, spacing-1 bordered AP).  kps-S39: the
order-3 members are DILATED APs (spacing d = numerator) + boundary defects.  This
species is what my spacing-1 sweep MISSED and where N=6 gets its member.

QUESTION: does any dilated-AP (spacing d>=2) + defect family at N=12 land in the
gap (1/13, 2/25)?  If NOT (and my construction REPRODUCES the N=6 member =
validation), that closes the order-3+ species at N=12 empirically and localizes
WHY N=12 differs from N=6.

For each candidate not in the gap: record the OVERSHOOT (M and its denominator,
parity, and whether an even speed drives it) -- the opus-S119/HYP-4592 binder gate.
"""
import itertools
from fractions import Fraction as F
from math import gcd

GAP = {6: (F(1, 7), F(2, 13)), 12: (F(1, 13), F(2, 25))}


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


def classify_fast(W, hi_f, slack):
    """early-exit: 'above' if some t beats hi+slack; else true float max."""
    cutoff = hi_f + slack
    best = 0.0
    for s in _dens(W):
        for j in range(1, s):
            t = j / s
            mv = min(abs(v * t - round(v * t)) for v in W)
            if mv > best:
                best = mv
                if best > cutoff:
                    return 'above', best
    return 'scan', best


def primitive(W):
    g = 0
    for v in W:
        g = gcd(g, v)
    return tuple(sorted(set(v // g for v in W))) if g > 1 else tuple(sorted(set(W)))


def dilated_ap_candidates(N, dmax=8, hmax=None):
    """Dilated-AP core {a, a+d, ..., a+(L-1)d} (spacing d) + boundary defects,
    assembled to N distinct primitive speeds.  Models kps-S39's order-3 species:
    N=6 target {1,5,6,11,16,17} = core {1,6,11,16} (d=5) + {5=6-1, 17=16+1}."""
    if hmax is None:
        hmax = 4 * N
    fams = set()
    for d in range(2, dmax + 1):
        for L in range(max(2, N - 4), N + 1):
            for a in range(1, d + 2):
                core = [a + i * d for i in range(L)]
                if core[-1] > hmax or len(set(core)) != L:
                    continue
                need = N - L
                # defect pool: inner borders (core[i] +- 1), outer borders,
                # and the "fill" points between consecutive core elements
                pool = set()
                for c in core:
                    pool.update([c - 1, c + 1])
                pool.update([core[0] - d + 1, core[-1] + 1, core[-1] + d - 1])
                pool = sorted(x for x in pool if x >= 1 and x not in core)
                if need == 0:
                    W = primitive(tuple(core))
                    if len(W) == N:
                        fams.add(W)
                elif need > 0 and len(pool) >= need:
                    for combo in itertools.combinations(pool, need):
                        W = primitive(tuple(sorted(set(core) | set(combo))))
                        if len(W) == N:
                            fams.add(W)
    return fams


def analyze(N, dmax=8, sample=None):
    lo, hi = GAP[N]
    lo_f, hi_f = float(lo), float(hi)
    fams = sorted(dilated_ap_candidates(N, dmax=dmax))
    if sample and len(fams) > sample:
        stride = len(fams) // sample
        fams = fams[::stride]
    in_gap = []
    above = 0
    below = 0
    overshoot_q = {}   # denominator -> count (for near-gap escapes)
    even_driven = 0
    near_examples = []
    for W in fams:
        tag, fm = classify_fast(W, hi_f, 0.004)
        if tag == 'above':
            above += 1
            continue
        if fm < lo_f - 0.01:
            below += 1
            continue
        M = exact_M(W)
        if lo < M < hi:
            in_gap.append((M, W))
        elif M >= hi:
            above += 1
            q = M.denominator
            overshoot_q[q] = overshoot_q.get(q, 0) + 1
            # is the overshoot driven by an even speed? (M's maximizer resonates a 2w)
            if any(v % 2 == 0 for v in W):
                even_driven += 1
            if len(near_examples) < 8 and M < lo + (hi - lo) * 6:
                near_examples.append((M, q, W))
        else:
            below += 1
    return fams, in_gap, above, below, overshoot_q, even_driven, near_examples


def main():
    print("=" * 86)
    print("VALIDATION at N=6: does the dilated-AP construction REPRODUCE {1,5,6,11,16,17}=5/33?")
    print("=" * 86)
    fams6, in_gap6, above6, below6, oq6, ed6, ne6 = analyze(6, dmax=6)
    print(f"  {len(fams6)} dilated-AP candidates at N=6; IN-GAP: {len(in_gap6)}")
    target6 = (1, 5, 6, 11, 16, 17)
    found = [ (M,W) for M,W in in_gap6 if W == target6 ]
    print(f"  target {list(target6)} in candidate set: {target6 in set(fams6)}")
    for M, W in sorted(in_gap6)[:8]:
        star = "  <== the known 5/33 member" if W == target6 else ""
        print(f"    IN GAP: M={M} ({float(M):.4f})  W={list(W)}{star}")
    if not in_gap6:
        print("    (none found -- construction still misses the interior-defect member; note it)")

    print()
    print("=" * 86)
    print("THE OPEN TEST at N=12: any dilated-AP (spacing d>=2) + defect in the gap (1/13,2/25)?")
    print("=" * 86)
    fams12, in_gap12, above12, below12, oq12, ed12, ne12 = analyze(12, dmax=8, sample=14000)
    print(f"  tested {len(fams12)} dilated-AP candidates at N=12")
    print(f"    IN GAP:  {len(in_gap12)}")
    print(f"    above:   {above12}")
    print(f"    below:   {below12}")
    if in_gap12:
        print("  *** GAP MEMBERS FOUND -- would REFUTE (G) at N=12! investigate:")
        for M, W in sorted(in_gap12)[:12]:
            print(f"      M={M} ({float(M):.5f})  W={list(W)}")
    else:
        print("  => NO dilated-AP order-3+ member in the gap at N=12 (empirical).")
    print()
    print(f"  near-gap escape denominators (M just above 2/25): top by count")
    for q, c in sorted(oq12.items(), key=lambda x: -x[1])[:10]:
        print(f"    q={q} ({'odd' if q % 2 else 'even'}):  {c} families")
    print(f"  fraction of above-gap escapes with an even speed present: {ed12}/{above12}")
    print()
    print("  near-gap examples (M in the first sixth above the gap):")
    for M, q, W in sorted(ne12)[:8]:
        print(f"    M={M} ({float(M):.5f}) q={q}  W={list(W)[:9]}..")
    print()
    print("  INTERPRETATION: N=6 nonempty (order-3 dilated AP) vs N=12 -- if N=12's")
    print("  dilated-AP species is ALSO gap-empty, the difference from N=6 is n-specific")
    print("  arithmetic (the spacing-d resonances at N=12 all overshoot), localizing")
    print("  opus-S119's residual to a finite dilated-AP + overshoot check.")


if __name__ == "__main__":
    main()
