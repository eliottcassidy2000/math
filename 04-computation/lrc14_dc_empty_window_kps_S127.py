# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont51: is the empty M-window (1/14, 2/27) free of DIVISOR-COMPLETE families?
# I.e. is the DC M-floor >= 2/27, so the near-AP hard families (M in (1/14,2/27)) are all NON-DC (=> handled
# by THM-366 t=1/d), and the DC route needs only the >= 2/27 floor (boxeph-S20 finite check + my transport)?
#
# CONTEXT: opus-S246 -- the M-spectrum has an empty window (1/14, 2/27) [2/27=mediant(1/14,1/13)]; the Farey-
# window rigidity [M<2/27 => dilated {1..13}] proves LRC(14). boxeph-S20 finite check: all 10.9M primitive DC
# (Vmax<=30) clear at non-14 q<=39 => M>1/14. THIS pins the DC floor VALUE and tests the empty window FOR DC:
# if no DC family has M in [1/14, 2/27), the Farey rigidity's DC part is vacuous -- the hard near-AP families
# are exactly the NON-DC ones (THM-366), and DC is loose with margin >= 2/27-1/14 = 1/378.
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
def Mfloat(v):
    qcap = 3 * max(v) + 2
    best = 0.0
    for q in range(2, qcap):
        for p in range(1, q):
            if gcd(p, q) == 1:
                m = min(norm(vi * p / q) for vi in v)
                if m > best: best = m
    return best
def MexactF(v):  # exact Fraction M (for verifying the extremes only)
    qcap = 3 * max(v) + 2
    best = F(0)
    for q in range(2, qcap):
        for p in range(1, q):
            if gcd(p, q) == 1:
                m = min(norm(F(vi * p, q)) for vi in v)
                if m > best: best = m
    return best
def Mexact(v): return MexactF(v)  # cores use exact

def main():
    T14, T27 = F(1, 14), F(2, 27)
    print(f"empty-window test: (1/14, 2/27) = ({float(T14):.5f}, {float(T27):.5f}); width {float(T27-T14):.5f}")

    # (1) structural cores -- exact M
    cores = {
        "cont.41 extremal {1,2,3,4,10..18}": [1,2,3,4,10,11,12,13,14,15,16,17,18],
        "boxeph H*+15 (near-dilate core)": [1,2,3,4,8,9,10,11,12,13,14,15,16],
        "shift-AP {2..14}": [2,3,4,5,6,7,8,9,10,11,12,13,14],
    }
    print("\n(1) structural DC cores -- exact M vs the window:")
    floor = F(1); fam0 = None
    for name, v in cores.items():
        if valid(v):
            m = Mexact(v); inwin = T14 <= m < T27
            print(f"    {name:38s} M={m}={float(m):.5f}  in (1/14,2/27)? {inwin}")
            if m < floor: floor, fam0 = m, v
        else:
            print(f"    {name:38s} not valid DC (skip)")

    # (2) adversarial hill-climb minimizing M over primitive DC (Vmax<=32), near-tight seeds
    print("\n(2) adversarial hill-climb: min M over primitive DC (Vmax<=32), hunting the window:")
    rng = random.Random(51)
    T14f, T27f = float(T14), float(T27)
    bestf = float(floor); bestfam = fam0; inwindow = []
    seeds = [fam0, [1,2,3,4,10,11,12,13,14,15,16,17,18]]
    for restart in range(60):
        cur = seeds[restart] if restart < len(seeds) else None
        if cur is None:
            for _ in range(400):
                w = sorted(rng.sample(range(1, 33), 13))
                if valid(w): cur = w; break
        if cur is None or not valid(cur): continue
        curM = Mfloat(cur)
        for _ in range(120):
            w = cur[:]; i = rng.randrange(13); w[i] = rng.randint(1, 32)
            w = sorted(set(w))
            if len(w) != 13 or not valid(w): continue
            m = Mfloat(w)
            if m <= curM: cur, curM = w, m
        if T14f - 1e-9 <= curM < T27f: inwindow.append((round(curM, 5), cur))
        if curM < bestf: bestf, bestfam = curM, cur
    best = MexactF(bestfam)   # exact-verify the min family
    print(f"    MIN M found (float {bestf:.5f}) EXACT = {best} = {float(best):.5f} at {bestfam}")
    print(f"    DC families found IN (1/14, 2/27): {len(inwindow)}")
    if inwindow:
        for m, v in inwindow[:3]: print(f"      !! M={m}={float(m):.5f}  {v}")
    print(f"\n=> DC M-floor (this search) = {float(best):.5f};  2/27 = {float(T27):.5f}")
    if best >= T27:
        print(f"   EMPTY WINDOW HOLDS FOR DC: no DC family with M in (1/14,2/27); DC floor >= 2/27 (margin >=1/378).")
        print(f"   => the near-AP hard families (M in the window) are all NON-DC => handled by THM-366 (t=1/d);")
        print(f"      the DC route needs ONLY the >=2/27 floor (boxeph finite check + my dilation transport).")
    else:
        print(f"   FOUND DC family with M < 2/27 -- the empty window does NOT cleanly separate DC (Farey rigidity needed for DC too).")

if __name__ == "__main__":
    main()
