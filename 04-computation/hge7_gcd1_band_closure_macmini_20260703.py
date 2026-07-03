#!/usr/bin/env python3
"""
DECISIVE test (mac-mini-2026-07-03-S21, HYP-3877):  does  gcd=1 + covering + >=7 far  =>  band-closable?
Loneliness is DILATION-INVARIANT: M(d*S)=M(S), lonely time scales by 1/d. So WLOG gcd(speeds)=1.
Under gcd=1 the tight dilated-AP {d,2d,..,13d} REDUCES to the window family {1,..,13} (farCount=0, handled
by hwindow). Claim: every gcd=1 covering family with >=7 far (>22) is 1/14-lonely at some a/q, q in {15..33}.

Adversarial hunt for a band FAILURE: random + near-AP (gcd=1) + GW-like + composite-far + min-M families.
Any covering gcd=1 hge7 family needing q>33 (or q>600) is a genuine counterexample to the band closure.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random

def gcd_all(xs): return reduce(gcd, xs)

def danger_residues(q):
    return {r for r in range(q) if min(r, q - r) * 14 < q}

def lonely_at_q(speeds, q):
    dang = danger_residues(q)
    for a in range(1, q):
        if gcd(a, q) == 1 and all((v * a) % q not in dang for v in speeds):
            return a
    return None

def smallest_band_q(speeds, qmax=600):
    for q in range(15, qmax + 1):
        if lonely_at_q(speeds, q) is not None:
            return q
    return None

def is_covering(speeds):
    return all(any(v % q == 0 for v in speeds) for q in range(2, 15))

def far_count(speeds): return sum(1 for v in speeds if abs(v) > 22)

if __name__ == "__main__":
    rng = random.Random(31337)
    print("DECISIVE: gcd=1 + covering + >=7 far  =>  band {15..33} closable?  hunting failures.")
    print("=" * 88)
    worst = (0, None); nfail = 0; ntot = 0
    over33 = []

    # (A) broad random gcd=1 covering hge7
    for _ in range(60000):
        nfar = rng.choice([7,8,9,10,11,12,13])
        nnear = 13 - nfar
        near = rng.sample(range(1, 23), nnear) if nnear else []
        kind = rng.choice(["near-equal","spread","2cluster","composite"])
        if kind == "near-equal":
            w1 = rng.randint(23, 8000); far = sorted({w1 + rng.randint(0, rng.choice([6,50,w1//2 or 1])) for _ in range(nfar)})
        elif kind == "spread":
            far = sorted(rng.sample(range(23, 9000), nfar))
        elif kind == "2cluster":
            a1 = rng.randint(23, 4000); a2 = rng.randint(23, 8000)
            far = sorted({rng.choice([a1,a2]) + rng.randint(0,20) for _ in range(nfar)})
        else:  # composite far (highly divisible in band)
            base = rng.choice([2520, 27720, 5040, 840, 1680, 9240])
            far = sorted({base + rng.randint(-3,3) + 24*rng.randint(0,3) for _ in range(nfar)} | {base})
            far = [f for f in far if f > 22]
        far = [f for f in far if f > 22]
        speeds = near + far
        if len(set(speeds)) != 13: continue
        if gcd_all(speeds) != 1: continue
        if far_count(speeds) < 7: continue
        if not is_covering(speeds): continue
        ntot += 1
        q = smallest_band_q(speeds, qmax=600)
        if q is None:
            nfail += 1
            if nfail <= 8: print("  BAND FAILURE (q>600):", speeds)
        else:
            if q > worst[0]: worst = (q, speeds)
            if q > 33: over33.append((q, speeds))

    # (B) targeted near-AP gcd=1 covering hge7 (APs that ARE covering with gcd 1)
    ap_hits = 0
    for start in range(1, 40):
        for step in range(1, 40):
            speeds = [start + step*i for i in range(13)]
            if len(set(speeds)) != 13: continue
            if gcd_all(speeds) != 1: continue
            if far_count(speeds) < 7: continue
            if not is_covering(speeds): continue
            ap_hits += 1; ntot += 1
            q = smallest_band_q(speeds, qmax=600)
            if q is None:
                nfail += 1; print("  AP BAND FAILURE:", speeds)
            else:
                if q > worst[0]: worst = (q, speeds)
                if q > 33: over33.append((q, speeds))

    print(f"\ntested {ntot} gcd=1 covering hge7 families ({ap_hits} covering gcd=1 APs among them)")
    print(f"band FAILURES (no q<=600): {nfail}")
    print(f"WORST band-q = {worst[0]}")
    print(f"   worst family: {worst[1]}")
    print(f"families needing q in (33,600]: {len(over33)}")
    for q, sp in sorted(over33, reverse=True)[:10]:
        print(f"   q={q}: {sp}  (gcd={gcd_all(sp)}, far={far_count(sp)})")
    verdict = "CLOSES" if nfail == 0 else "has FAILURES"
    print(f"\n=> gcd=1 + covering + >=7 far  {verdict}  by band {{15..{worst[0]}}}  (0 failures needed for closure)")
