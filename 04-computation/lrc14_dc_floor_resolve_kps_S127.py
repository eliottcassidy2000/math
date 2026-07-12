# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont52: SHARPEN THE CRUX -- resolve the exact DC (covering) M-floor.
# The crux is now DC => M > 1/14 (klein-S266: >=1/13; my cont.41/51: extremal at 1/12=0.0833).
# OPEN: is the floor EXACTLY 1/12, or is there a DC family with M in [1/13, 1/12) (=> floor lower)?
# This searches hard for sub-1/12 DC families, seeded by (a) the tight non-DC families AP {1..13} and
# Goddyn-Wong {1..11,13,24} perturbed toward DC, (b) the known 1/12 extremal, (c) 2-block/near-AP structures,
# and characterizes the extremal(s). Also checks: is GW covering (DC)? what forces covering => M>=1/13?
from math import gcd
from functools import reduce
from fractions import Fraction as F
import random, itertools

def is_DC(v):  return all(any(x % d == 0 for x in v) for d in range(2, 15))
def prim(v):   return reduce(gcd, v) == 1
def lrun(v):
    v = sorted(set(v)); b = m = 1
    for i in range(1, len(v)):
        if v[i] == v[i-1] + 1: m += 1; b = max(b, m)
        else: m = 1
    return b
def valid(v): return len(set(v)) == 13 and prim(v) and is_DC(v)
def norm(x): r = x - int(x); r = r + 1 if r < 0 else r; return min(r, 1 - r)
def Mfloat(v):
    qcap = 3 * max(v) + 2; best = 0.0
    for q in range(2, qcap):
        for p in range(1, q):
            if gcd(p, q) == 1:
                m = min(norm(vi * p / q) for vi in v)
                if m > best: best = m
    return best
def Mexact(v):
    qcap = 3 * max(v) + 2; best = F(0)
    for q in range(2, qcap):
        for p in range(1, q):
            if gcd(p, q) == 1:
                m = min(norm(F(vi * p, q)) for vi in v)
                if m > best: best = m
    return best

def main():
    F12, F13, F14 = F(1,12), F(1,13), F(1,14)
    print(f"targets: 1/12={float(F12):.5f}  1/13={float(F13):.5f}  1/14={float(F14):.5f}")

    # (A) the two tight families -- are they covering (DC)?
    AP = list(range(1,14)); GW = [1,2,3,4,5,6,7,8,9,10,11,13,24]
    for name,v in [("AP {1..13}",AP),("Goddyn-Wong {1..11,13,24}",GW)]:
        miss=[d for d in range(2,15) if not any(x%d==0 for x in v)]
        print(f"  {name:26s} M={Mexact(v)}  DC?{not miss}  missing mults of {miss}")
    print("  => both tight families are NON-covering (missing mult of 14) => THM-366 (t=1/14). The crux excludes them.")

    # (B) exhaustive-ish over 2-block + near-{known extremal} DC families for M < 1/12
    print("\n(B) hunting DC families with M < 1/12 (below the known extremal):")
    known = [1,2,3,4,10,11,12,13,14,15,16,17,18]
    found_below = []
    cand = set()
    # perturb known extremal + AP + GW by swaps, keep DC primitive
    seeds = [known, AP, GW, [1,2,3,10,11,12,13,14,15,16,17,18,19], [2,3,4,5,6,7,8,9,10,11,12,13,14]]
    rng = random.Random(52)
    for s in seeds:
        for _ in range(4000):
            w = s[:]
            for _ in range(rng.randint(1,3)):
                w[rng.randrange(13)] = rng.randint(1, 60)
            w = sorted(set(w))
            if len(w)==13 and valid(w): cand.add(tuple(w))
    # also all 2-block families {a..a+r} u {b..b+s} small
    for a in range(1,6):
        for r in range(2,7):
            for b in range(a+r+1, 30):
                for s2 in range(2,10):
                    blk=list(range(a,a+r))+list(range(b,b+s2))
                    if len(blk)==13:
                        blk=sorted(set(blk))
                        if len(blk)==13 and valid(blk): cand.add(tuple(blk))
    print(f"    candidate DC families: {len(cand)}")
    minM=F(1); minfam=None
    for t in cand:
        v=list(t); mf=Mfloat(v)
        if mf < float(F12)+1e-9:
            me=Mexact(v)
            if me < F12: found_below.append((me,v))
            if me < minM: minM,minfam=me,v
    print(f"    min exact M over candidates = {minM}={float(minM):.5f} at {minfam}")
    print(f"    DC families with M STRICTLY < 1/12: {len(found_below)}")
    for me,v in sorted(found_below)[:6]:
        print(f"      M={me}={float(me):.5f}  in [1/13,1/12)? {F13<=me<F12}  {v}  (2-block? lrun={lrun(v)})")

    # (C) verdict
    print("\n(C) VERDICT on the DC floor:")
    if not found_below:
        print(f"    NO DC family below 1/12 found => DC floor = 1/12 (klein's >=1/13 is a conservative provable bound;")
        print(f"    the exact extremal is 1/12 at {known}).")
    else:
        fl = min(me for me,_ in found_below)
        print(f"    DC families BELOW 1/12 exist; sharpest found M={fl}={float(fl):.5f} => the DC floor is LOWER than 1/12,")
        print(f"    approaching klein's 1/13 -- the crux margin is smaller than cont.41/51 claimed.")

if __name__ == "__main__":
    main()
