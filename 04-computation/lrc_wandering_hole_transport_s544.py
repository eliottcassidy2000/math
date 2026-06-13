#!/usr/bin/env python3
"""
lrc_wandering_hole_transport_s544.py    oracle-2026-06-01-S544

Creative reframe: "global spread guarantees local emptiness."

Among the n points (observer at 0 + n-1 runners at v_i t), the gaps sum to 1, so by
PIGEONHOLE the LARGEST gap is always >= 1/n. So a region of "local emptiness" (a
hole of width >= 1/n) ALWAYS EXISTS, at every instant -- guaranteed by the global
spread. LRC is therefore NOT about *creating* emptiness; it is about TRANSPORTING
the guaranteed hole to the observer.

We track the wandering hole and separate the two halves:
  G(t) = the largest gap among all n points        (>= 1/n ALWAYS -- the guarantee)
  O(t) = the observer's collar = min_i ||v_i t||    (>= 1/n = LRC, the transport)
The hole exists for free; LRC = a gap of width >= 2/n sits AROUND the observer at
some t (so O(t) >= 1/n). The apex/source-sink (S530) is the hole; the fat collar
(S541) is O. The hard part is the hole's TRANSPORT to x=0, not its existence.
"""
from fractions import Fraction
from functools import reduce
from math import gcd

def frac(x): return x - (x.numerator // x.denominator)
def d0(x):
    f = frac(x); return min(f, 1 - f)

def sample_times(speeds, n):
    W = set([Fraction(0), Fraction(1)])
    # collision walls (two points coincide) AND danger walls (runner at 1/n from obs)
    runs = speeds[1:]
    for i in range(len(runs)):
        for j in range(len(runs)):
            d = abs(runs[i] - runs[j])
            if d == 0: continue
            for k in range(0, 2 * d): W.add(frac(Fraction(k, 2 * d)))
        s = runs[i]
        if s == 0: continue
        for m in range(0, abs(s) + 1):
            for sg in (1, -1):
                t = Fraction(n * m + sg, n * abs(s))
                if 0 <= t <= 1: W.add(t)
    Wl = sorted(w for w in W if 0 <= w < 1)
    mids = [(a + b) / 2 for a, b in zip(Wl, Wl[1:] + [Fraction(1)])]
    return Wl + mids

def gaps(positions):
    p = sorted(positions)
    g = [p[(i + 1) % len(p)] - p[i] for i in range(len(p) - 1)]
    g.append(1 - p[-1] + p[0])
    return g

def analyze(speeds, n):
    thr = Fraction(1, n)
    minG = Fraction(2); maxO = Fraction(0); twofold = 0; total = 0
    for t in sample_times(speeds, n):
        pos = [frac(Fraction(s) * t) for s in speeds]   # observer at 0 included (speed 0)
        G = max(gaps(pos))
        O = min(d0(Fraction(s) * t) for s in speeds[1:])  # observer collar
        minG = min(minG, G); maxO = max(maxO, O)
        if G >= 2 * thr: twofold += 1
        total += 1
    return dict(minG=minG, maxO=maxO, twofold_frac=Fraction(twofold, total), thr=thr)

def primitive(s): return reduce(gcd, [x for x in s if x]) == 1

def main():
    print("The wandering hole: global spread guarantees local emptiness (oracle-S544)\n")
    print("  guarantee: min_t (largest gap) >= 1/n always.  transport: max_t (observer collar) >= 1/n = LRC.")
    print(f"  {'system':<18}{'1/n':>7}{'min_t G (hole)':>16}{'>=1/n?':>8}{'max_t O (collar)':>18}{'LRC':>6}{'frac G>=2/n':>13}")
    cases = [("AP n=5 (tight)",(0,1,2,3,4),5),("generic n=5",(0,1,3,5,7),5),
             ("AP n=7 (tight)",(0,1,2,3,4,5,6),7),("generic n=7",(0,1,2,4,8,9,11),7),
             ("AP n=14 (tight)",tuple(range(14)),14),("generic n=14",(0,1,2,3,4,5,6,7,8,9,10,11,12,20),14)]
    for name, sp, n in cases:
        if not primitive(sp): continue
        r = analyze(sp, n)
        print(f"  {name:<18}{str(r['thr']):>7}{str(r['minG']):>16}{str(r['minG']>=r['thr']):>8}"
              f"{str(r['maxO']):>18}{str(r['maxO']>=r['thr']):>6}{float(r['twofold_frac']):>13.3f}")
    print()
    print("READING: min_t G >= 1/n for EVERY system (the hole/local-emptiness is GUARANTEED")
    print("at all times by global spread -- pigeonhole). LRC = max_t O >= 1/n (the guaranteed")
    print("hole reaches the observer). For the tight AP, min_t G = 1/n AND max_t O = 1/n: the")
    print("hole is exactly the minimum size and only barely reaches the observer (boundary).")
    print("Generic: G and O comfortably exceed 1/n. So the conjecture is NOT 'does emptiness")
    print("exist' (it always does) but 'is the guaranteed defect transportable to x=0' -- and")
    print("the AP is the system where the transport is critically tight.")

if __name__ == "__main__":
    main()
