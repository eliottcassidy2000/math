#!/usr/bin/env python3
"""
WIDE RESIDUAL gap hunt (mac-mini-2026-07-03-S26): find compressed covering gcd=1 families that are
NOT covered by {census(q<=45), THM-608(speed-tight), scale_separation_phase(far phase-tight at witness)}.

Correct test: at a family's MINIMAL witness a/q (all runners safe), the FAR cluster phases sit in the safe
region [1/14,13/14]. scale_separation_phase can peel the far ONLY if the far phase-span < 1/14 (fits one
danger arc, placeable by the base's slack). THM-608 needs far SPEED-spread D small. Census needs q<=45.
GAP = large witness q (>45) AND far phase-span >= 1/14 at every good small t0 -- lonely only by the full
(possibly large-denominator) witness, no cluster peel, census unreachable.

Hunt over: tight locus (deep well + dilations gcd=1), aligned band-blockers, wide clusters, GW-like, random.
"""
from math import gcd
from functools import reduce
import random

def gcd_all(xs): return reduce(gcd, xs)
def nd(x): x = x % 1; return min(x, 1-x)
def is_covering(sp): return all(any(v % q == 0 for v in sp) for q in range(2,15))

def min_witness(speeds, qmax=400):
    for q in range(2, qmax+1):
        for a in range(1, q):
            if gcd(a,q)!=1: continue
            if all(nd(v*a/q) >= 1/14 for v in speeds): return q, a
    return None, None

def far_phase_span(far, a, q):
    ph = sorted((v*a/q) % 1 for v in far)
    if len(ph) < 2: return 0.0
    gaps = [ph[i+1]-ph[i] for i in range(len(ph)-1)] + [1 - ph[-1] + ph[0]]
    return 1 - max(gaps)

def peelable_phase_tight(near, far, qmax=45):
    """exists a small good t0 (base=near lonely) at which far phase-span < 1/14?"""
    for q in range(15, qmax+1):
        for a in range(1, q):
            if gcd(a,q)!=1: continue
            if near and not all(nd(v*a/q) >= 1/14 for v in near): continue
            if far_phase_span(far, a, q) < 1/14: return True
    return False

def classify(speeds):
    near = [v for v in speeds if v <= 22]; far = [v for v in speeds if v > 22]
    D = max(far)-min(far) if far else 0
    q, a = min_witness(speeds, qmax=400)
    if q is None: return "NO-WITNESS", D, None
    if q <= 45: return "census", D, q
    if D <= 12: return "THM-608(speed-tight)", D, q
    if peelable_phase_tight(near, far): return "scale_sep_phase(phase-tight)", D, q
    return "GAP", D, q

if __name__ == "__main__":
    rng = random.Random(260)
    fams = []
    # tight / near-tight
    fams.append(("deep well {1..12,182}", list(range(1,13))+[182]))
    for d in [2,3,5,7]:
        fams.append((f"{d}*deepwell", [d*x for x in list(range(1,13))+[182]]))
    # aligned band-blockers (near-equal far, q~45-53)
    for N in [1000, 5000, 30000]:
        band=list(range(15,60)); rng.shuffle(band)
        far=sorted({q*round(N/q) for q in band[:9]}); far=[f for f in far if f>22]
        sp=far[:]
        for q in [8,9,5,7,11,13,2,3,4,6]:
            if len(sp)>=13: break
            if not any(s%q==0 for s in sp): sp.append(q)
        while len(sp)<13: sp.append(rng.randint(2,22))
        fams.append((f"aligned~{N}", sorted(set(sp))[:13]))
    # engineered WIDE clusters (far spanning a wide range)
    for _ in range(2000):
        nnear=rng.randint(0,6); nfar=13-nnear
        near=rng.sample(range(1,23),nnear) if nnear else []
        N=rng.choice([100,500,2000]); D=rng.choice([N//2, N, 2*N])
        far=sorted({rng.randint(23, 23+N+D) for _ in range(nfar)})
        if len(far)<nfar: continue
        fams.append(("wide", sorted(set(near+far))))

    print("Wide-residual gap hunt: classify hard compressed covering gcd=1 families.")
    print("=" * 90)
    counts={}; gaps=[]
    for name, sp in fams:
        sp=sorted(set(int(x) for x in sp))
        if len(sp)!=13 or gcd_all(sp)!=1 or not is_covering(sp) or not any(v>22 for v in sp): continue
        tag, D, q = classify(sp)
        counts[tag]=counts.get(tag,0)+1
        if tag=="GAP": gaps.append((name, sp, D, q))
    for k in sorted(counts, key=lambda k:-counts[k]):
        print(f"   {k:>32}: {counts[k]}")
    print(f"\nGAP families (large witness q>45 AND far NOT phase-tight): {len(gaps)}")
    for name, sp, D, q in gaps[:8]:
        far=[v for v in sp if v>22]
        print(f"   [{name}] far-span D={D}, witness q={q}: far={far[:6]}{'...' if len(far)>6 else ''}")
    print("\n=> GAP=0: the three routes cover the compressed case (wide residual closed).")
    print("   GAP>0: those families are the genuine remaining crux (lonely only by a large-q witness, no peel).")
