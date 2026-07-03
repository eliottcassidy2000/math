#!/usr/bin/env python3
"""
WIDE RESIDUAL coverage diagnostic (mac-mini-2026-07-03-S26): does {census + THM-608(speed-tight) +
scale_separation_phase(phase-tight/Eisenstein)} cover ALL compressed covering families, or is there a
generic/sweep gap?

For a compressed covering gcd=1 family (slow base <=22 + far cluster ~N, possibly WIDE), classify its
minimal rational witness a/q:
 - CENSUS: q small (<= 45) -- pinned/loose, bounded-magnitude reachable.
 - PHASE-TIGHT at t0: at SOME good rational t0 (small q), the far cluster {(N+c_i)t0 mod 1} has phase-span
   < 1/14 (fits one danger-arc) => scale_separation_phase peels it (opus). [includes the 13-comb Eisenstein.]
 - SPEED-TIGHT: far speed-spread D small (THM-608 condition ii) -- near-equal.
 - GAP (sweep): none of the above at any small q -- would need the ergodic sweep / large-denominator witness.
"""
from math import gcd
from functools import reduce
import random

def gcd_all(xs): return reduce(gcd, xs)
def nd(x): x = x % 1; return min(x, 1-x)
def is_covering(sp): return all(any(v % q == 0 for v in sp) for q in range(2,15))

def min_witness_q(speeds, qmax=200):
    for q in range(2, qmax+1):
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            if all(nd(v*a/q) >= 1/14 for v in speeds): return q, a
    return None, None

def phase_span_far(far, a, q):
    """span of far phases {(N+c_i)*a/q mod 1} on the circle (min covering arc)."""
    ph = sorted((v*a/q) % 1 for v in far)
    if not ph: return 0.0
    gaps = [ph[(i+1)%len(ph)] - ph[i] for i in range(len(ph)-1)] + [1 - ph[-1] + ph[0]]
    return 1 - max(gaps)   # span = 1 - largest gap (points fit in an arc of this length)

def classify(speeds):
    near = sorted(v for v in speeds if v <= 22)
    far = sorted(v for v in speeds if v > 22)
    D = (far[-1]-far[0]) if far else 0
    q, a = min_witness_q(speeds, qmax=200)
    if q is None: return "NO-WITNESS<=200", D, None
    # census: small q
    if q <= 45:
        # phase-tight at this witness?
        pspan = phase_span_far(far, a, q)
        tag = "CENSUS(q<=45)"
        if pspan < 1/14: tag += "+phase-tight"
        return tag, D, q
    # q > 45: check phase-tight at some small good t0
    for qq in range(15, 46):
        for aa in range(1, qq):
            if gcd(aa,qq)!=1: continue
            if all(nd(v*aa/qq) >= 1/14 for v in near):   # base safe
                if far and phase_span_far(far, aa, qq) < 1/14:
                    return "PHASE-TIGHT(scale_sep_phase)", D, q
    if D <= 12:
        return "SPEED-TIGHT(THM-608)", D, q
    return "GAP(sweep?)", D, q

if __name__ == "__main__":
    rng = random.Random(26)
    print("Wide-residual coverage: classify compressed covering gcd=1 families (slow base + far cluster).")
    print("=" * 92)
    counts = {}; gaps = []
    for _ in range(4000):
        # compressed: slow base (<=6, <=22) + far cluster ~N (7+), possibly WIDE
        nnear = rng.randint(0,6); nfar = 13 - nnear
        near = rng.sample(range(1,23), nnear) if nnear else []
        N = rng.choice([50, 200, 1000, 5000])
        width = rng.choice([6, 30, N//4, N//2, N])   # near-equal .. wide
        far = sorted({N + rng.randint(0, max(1,width)) for _ in range(nfar)})
        if len(far) < nfar: continue
        speeds = sorted(set(near + far))
        if len(speeds) != 13 or gcd_all(speeds) != 1 or not is_covering(speeds): continue
        if not any(v > 22 for v in speeds): continue
        tag, D, q = classify(speeds)
        counts[tag] = counts.get(tag, 0) + 1
        if tag.startswith("GAP"):
            gaps.append((speeds, D, q))
    print("route coverage counts:")
    for k in sorted(counts, key=lambda k: -counts[k]):
        print(f"   {k:>34}: {counts[k]}")
    print(f"\nGAP (sweep?) families: {len(gaps)}")
    for sp, D, q in gaps[:6]:
        print(f"   far-span D={D}, min-witness q={q}: {sp}")
    print("\n=> if GAP=0: census + THM-608 + scale_separation_phase COVER all compressed families (wide residual")
    print("   closed by the three routes). if GAP>0: the generic/sweep case is the genuine remaining crux.")
