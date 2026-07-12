# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont50: probing the DC M-FLOOR -- is inf M over primitive divisor-complete families
# bounded away from 1/14, or do near-tight DC families approach it?
#
# CONTEXT: post-MISTAKE-140 (boxeph), M is dilation-invariant, so the DC class stratifies by STRUCTURE not
# diameter. The near-dilate slice has M in [1/13, 1/11] (boxeph H*). My cont.41 exhaustive census (Vmax<=22)
# found M-floor 1/12 at {1,2,3,4,10..18}. THE OPEN QUESTION: is the DC M-floor bounded away from 1/14 (=> LRC14
# for DC holds WITH MARGIN), or can DC families be made near-tight (M -> 1/14)? The tight point {1..13} (M=1/14)
# is NON-DC (no mult of 14); DC excludes it -- but how close can DC get? This searches hard, WITH the
# MISTAKE-140-mandated near-dilate + near-{1..13} seed battery (not just random/narrow-band).
from math import gcd
from functools import reduce
from fractions import Fraction as F
import random

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
def Mexact(v, qcap=None):
    if qcap is None: qcap = 3 * max(v) + 2
    best = F(0)
    for q in range(2, qcap):
        # only q that can beat current best: need q with a p giving min-frac > best; scan all p coprime
        for p in range(1, q):
            if gcd(p, q) == 1:
                m = min(norm(F(vi * p, q)) for vi in v)
                if m > best: best = m
    return best

def main():
    T = F(1, 14)
    print(f"target tight bound 1/14 = {float(T):.5f};  1/13={1/13:.5f}  1/12={1/12:.5f}  1/11={1/11:.5f}")

    # (1) KNOWN low-M DC families (exact M)
    known = {
        "cont.41 extremal {1,2,3,4,10..18}": [1,2,3,4,10,11,12,13,14,15,16,17,18],
        "boxeph H* + delta (near-dilate core)": [1,2,3,4,8,9,10,11,12,13,14,16,15],  # H*(12) + 15 as the 13th
        "shift-AP {2..14}": [2,3,4,5,6,7,8,9,10,11,12,13,14],
    }
    print("\n(1) known / seed DC families -- exact M:")
    floor = F(1); fam0 = None
    for name, v in known.items():
        if valid(v):
            m = Mexact(v); print(f"    {name:42s} M = {m} = {float(m):.5f}")
            if m < floor: floor, fam0 = m, v
        else:
            print(f"    {name:42s} NOT valid DC/primitive (skip)")

    # (2) near-{1..13} DC families: perturb the tight point to make it DC, see how low M goes
    print("\n(2) near-tight DC families (perturbations of {1..13} that ARE divisor-complete):")
    base = list(range(1, 14))
    cand = []
    # replace one element with a value that adds the missing mult-of-14 while staying primitive+DC
    for drop in range(1, 14):
        for add in range(14, 60):
            v = sorted(set([x for x in base if x != drop] + [add]))
            if len(v) == 13 and valid(v):
                cand.append(v)
    # also two-swaps
    for d1 in range(1,14):
        for a1 in [14,28,42]:
            for d2 in range(1,14):
                for a2 in range(15,40):
                    v = sorted(set([x for x in base if x not in (d1,d2)] + [a1,a2]))
                    if len(v)==13 and valid(v): cand.append(v)
    seen=set(); uniq=[]
    for v in cand:
        t=tuple(v)
        if t not in seen: seen.add(t); uniq.append(v)
    lo = F(1); loex=None
    for v in uniq:
        m = Mexact(v)
        if m < lo: lo, loex = m, v
    print(f"    searched {len(uniq)} near-{{1..13}} DC families; MIN M = {lo} = {float(lo):.5f} at {loex}")
    if lo < floor: floor, fam0 = lo, loex

    # (3) adversarial hill-climb minimizing M (bounded diameter, exact M), near-tight seeds
    print("\n(3) adversarial hill-climb minimizing M (Vmax<=45, exact M):")
    rng = random.Random(50)
    best = floor
    seeds = [fam0] + uniq[:20] + [[1,2,3,4,10,11,12,13,14,15,16,17,18]]
    for s in range(80):
        cur = seeds[s % len(seeds)][:] if s < len(seeds) else None
        if cur is None:
            for _ in range(300):
                w = sorted(rng.sample(range(1,46),13))
                if valid(w): cur = w; break
        if cur is None or not valid(cur): continue
        curM = Mexact(cur)
        for _ in range(150):
            w = cur[:]; i = rng.randrange(13); w[i] = rng.randint(1, 45)
            w = sorted(set(w))
            if len(w)!=13 or not valid(w): continue
            m = Mexact(w)
            if m <= curM: cur, curM = w, m
        if curM < best: best, fam0 = curM, cur
    print(f"    hill-climb MIN M = {best} = {float(best):.5f} at {fam0}")

    print(f"\n=> DC M-FLOOR (this search) = {best} = {float(best):.5f}")
    print(f"   margin over 1/14: {float(best) - float(T):.5f}  ({'BOUNDED AWAY' if best > T else 'AT/BELOW'} 1/14)")
    print(f"   => if the floor holds at ~{float(best):.4f} >> 1/14={float(T):.4f}, LRC(14)-for-DC has a REAL margin;")
    print(f"      the hard core is NOT near-tight -- consistent with cont.41 (1/12) + boxeph near-dilate [1/13,1/11].")

if __name__ == "__main__":
    main()
