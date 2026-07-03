#!/usr/bin/env python3
"""
DECISIVE gap hunt (mac-mini-2026-07-03-S26): is every compressed covering gcd=1 family either
SMALL-q (census-able, general p/q) or SMALL-D (far speed-spread small => THM-608 speed-tight)?
If so, {census + THM-608} cover the compressed case and there is NO wide-large-q gap.
The genuine gap = a family with BOTH large min-witness q AND large far-span D (and not phase-tight).

Correction context: opus-S52 claimed the deep well {1..12,182} needs q=183 (census fails). FALSE: it is
lonely at 3/40 (general p/q, min-dist 3/40 > 1/14). The census uses general p/q (kps lonely14_of_ratio),
not just the covering sieve p=1. So the Eisenstein 183 is not needed for the deep well.
"""
from math import gcd
from functools import reduce
import random

def gcd_all(xs): return reduce(gcd, xs)
def nd(x): x = x % 1; return min(x, 1-x)
def is_covering(sp): return all(any(v % q == 0 for v in sp) for q in range(2,15))
def min_witness_q(sp, qmax=600):
    for q in range(2, qmax+1):
        for a in range(1, q):
            if gcd(a,q)!=1: continue
            if all(nd(v*a/q) >= 1/14 for v in sp): return q
    return None
def far_span(sp):
    far=[v for v in sp if v>22]; return (max(far)-min(far)) if far else 0

if __name__ == "__main__":
    rng = random.Random(2626)
    print("Joint (min-witness q, far-span D) over compressed covering gcd=1 families. Hunt large-q AND large-D.")
    print("=" * 92)
    worst_q = (0, None); worst_qD = (0, None)   # max q, and max q among wide-D (D>40) families
    bothbig = []   # q>45 AND D>40
    ntest = 0
    for trial in range(120000):
        nnear = rng.randint(0,6); nfar = 13-nnear
        near = rng.sample(range(1,23), nnear) if nnear else []
        N = rng.choice([50,200,1000,5000,20000])
        style = rng.choice(["near-equal","wide","aligned","2cluster"])
        if style=="near-equal":
            far = sorted({N+rng.randint(0, rng.choice([6,20,50])) for _ in range(nfar)})
        elif style=="wide":
            far = sorted({rng.randint(23, N+1) for _ in range(nfar)})
        elif style=="aligned":
            band=list(range(15,70)); rng.shuffle(band)
            far=sorted({q*round(N/q) for q in band[:nfar]}); far=[f for f in far if f>22]
        else:
            a1=rng.randint(23,N//2 or 24); a2=rng.randint(N//2 or 24, N+1)
            far=sorted({rng.choice([a1,a2])+rng.randint(0,15) for _ in range(nfar)})
        sp = sorted(set(near + far))
        if len(sp)!=13 or gcd_all(sp)!=1 or not is_covering(sp) or not any(v>22 for v in sp): continue
        ntest += 1
        q = min_witness_q(sp, qmax=600)
        if q is None:
            print("  NO WITNESS <=600:", sp); continue
        D = far_span(sp)
        if q > worst_q[0]: worst_q = (q, (sp, D))
        if D > 40 and q > worst_qD[0]: worst_qD = (q, (sp, D))
        if q > 45 and D > 40:
            bothbig.append((q, D, sp))
    print(f"tested {ntest} compressed covering gcd=1 families")
    print(f"MAX min-witness q = {worst_q[0]}  (D={worst_q[1][1]}, far-span)")
    print(f"MAX q among WIDE (D>40) families = {worst_qD[0]}")
    print(f"families with BOTH q>45 AND D>40 (the potential gap): {len(bothbig)}")
    for q, D, sp in sorted(bothbig, reverse=True)[:8]:
        print(f"   q={q}, D={D}: {[v for v in sp if v>22][:6]}...")
    print("\n=> if all large-q families have SMALL D (worst_qD small): large-q => near-equal => THM-608;")
    print("   {census(small q) + THM-608(small D)} COVER the compressed case, wide residual CLOSED, no")
    print("   Eisenstein/scale_sep_phase gap. if bothbig>0: those are the genuine remaining crux.")
