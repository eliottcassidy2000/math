#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD 3 -- SHARPEN ANGLE 2: the POSITIVE resonance error is what must be bounded for the
WIDE upper bound, and it shrinks with k while the margin grows. Test scaling k=8..12.

For the wide upper bound p0(E) <= cap_k we have p0 = p0_decorr + err, p0_decorr <= Q(k-1),
margin = cap_k - Q(k-1). The bound CLOSES iff  err <= margin  for every wide E.
Only the POSITIVE part of err matters. This script:
  (1) over consec base + commensurable far pair, finds max POSITIVE err per k;
  (2) over MANY bounded bases + commensurable far pair, finds the global worst positive err per k;
  (3) reports margin and the ratio err+ / margin (the "lossiness" needed).
If err+ / margin stays small and ideally shrinks, the loose E_2/Koksma bound is ample.
"""
import sys, random
from fractions import Fraction as F
from functools import reduce
from math import gcd, lcm
sys.stdout.reconfigure(line_buffering=True)

caps = {8: F(2243,5880), 9: F(1979,4004), 10: F(55,91), 11: F(66,91), 12: F(6,7)}
# Q(k-1) values from canon decorrelated audit:
Qkm1 = {8: F(289,1470), 9: F(621,1715), 10: F(1229,2744), 11: F(65599,123480), 12: F(14873,24696)}

def p0(E):
    E = sorted(set(int(e) for e in E))
    bps = {F(0), F(1)}
    for e in E:
        ae = abs(e)
        if ae == 0: continue
        for m in range(7 * ae + 1): bps.add(F(m, 7 * ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    tot = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        mid = (lo + hi) / 2
        s = set()
        for e in E:
            v = F(e) * mid; v = v - F(v.numerator // v.denominator)
            s.add((v.numerator * 7) // v.denominator)
        if len(s & set(range(1, 7))) == 6: tot += hi - lo
    return tot

def p0_decorr(B, nfar=2, samples=8000, seed=0):
    Bs = sorted(set(int(e) for e in B))
    bps = {F(0), F(1)}
    for e in Bs:
        ae = abs(e)
        if ae == 0: continue
        for m in range(7 * ae + 1): bps.add(F(m, 7 * ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    rng = random.Random(seed); tot = 0.0
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        mid = (lo + hi) / 2
        s = set()
        for e in Bs:
            v = F(e) * mid; v = v - F(v.numerator // v.denominator)
            s.add((v.numerator * 7) // v.denominator)
        need = set(range(1, 7)) - s
        if not need:
            tot += float(hi - lo); continue
        hit = 0
        for _ in range(samples):
            got = set(rng.randrange(7) for _ in range(nfar))
            if need <= got: hit += 1
        tot += float(hi - lo) * hit / samples
    return tot

def ratios():
    R = []
    for q in range(1, 7):
        for p in range(q + 1, 3 * q + 1):
            if gcd(p, q) == 1 and 1 < p / q <= 2.2: R.append((p, q))
    return sorted(set(R), key=lambda t: t[0]/t[1])

def main():
    C = 300; Rs = ratios()
    print("k | margin=cap-Q | consec-base max(err+) | global max(err+) over bases | err+/margin")
    for k in (8, 9, 10, 11, 12):
        margin = float(caps[k] - Qkm1[k])
        consec = list(range(k - 2))
        bdec_c = p0_decorr(consec, 2, 9000, seed=k)
        # consec base err+
        cmax = -9
        for (p, q) in Rs:
            E = sorted(set(consec + [C*q, C*p]))
            if len(E) != k or reduce(gcd, E) != 1: continue
            err = float(p0(E)) - bdec_c
            cmax = max(cmax, err)
        # global over a few bounded bases (consec, even-AP, random spreads)
        bases = [consec, [0] + [2*i for i in range(1, k-2)]]
        rng = random.Random(100 + k)
        for _ in range(12):
            bases.append([0] + sorted(rng.sample(range(1, 14), k - 3)))
        gmax = -9
        for B in bases:
            B = sorted(set(B))
            if len(B) != k - 2: continue
            bdec = p0_decorr(B, 2, 4000, seed=k*7)
            for (p, q) in Rs:
                E = sorted(set(B + [C*q, C*p]))
                if len(E) != k or reduce(gcd, E) != 1: continue
                err = float(p0(E)) - bdec
                gmax = max(gmax, err)
        print(f"{k} | {margin:.4f} | {cmax:+.4f} | {gmax:+.4f} | {gmax/margin:.3f}  "
              f"(=> need {1/(gmax/margin):.0f}x-lossy bound)")
    print("\n=> if err+/margin shrinks (or stays <<1), the loose E_2/Koksma resonance bound closes WIDE.")

if __name__ == "__main__":
    main()
