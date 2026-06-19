#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Does the LRC spectral gap dip below Theta(1/k^2) UNBOUNDEDLY?  i.e. does a_max(k)
grow, where sigma_2(k) = a_max/(a_max(k+1)-1) (the Stern-Brocot mediant just above
the floor 1/(k+1))?   g(k)*k^2 -> 1/a_max.   kind-pasteur-2026-06-19-S9.

Method: FILTERED exhaustive scan for configs with M(S) < mediant = 2/(2k+1).
SOUND FILTER: M(S) < mediant  <=>  NO time t has min_v||vt|| >= mediant.
So probe small-denominator times; if ANY gives min >= mediant, M >= mediant => SKIP
(never discards a true below-mediant config). Survivors get an EXACT maxmin.
This makes a much larger speed box tractable (below-mediant configs are very rare).

For every below-mediant config we record M=p/q, a (if e=p(k+1)-q==1), maxS, the
binding pair, and t*.  a_max(k) = max a seen.  If a_max grows with k => gap dips.
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None

def fn(x):  # ||x||
    r = x - (x.numerator // x.denominator)
    return min(r, 1 - r)

def maxmin_exact(S):
    S = sorted(set(S)); cand = set()
    for i in range(len(S)):
        for j in range(i, len(S)):
            for d in (S[i] + S[j], S[i] - S[j]):
                if d == 0: continue
                d = abs(d)
                for m in range(1, d): cand.add(Fraction(m, d))
    best = Fraction(0); bt = Fraction(0)
    for t in cand:
        mv = min(fn(v * t) for v in S)
        if mv > best: best, bt = mv, t
    return best, bt

def primitive(S): return reduce(gcd, S) == 1

def build_probe_times(k, Qmax):
    """reduced fractions m/q, 1<=m<q, q<=Qmax (the small-denominator candidate optima)."""
    pts = []
    for q in range(2, Qmax + 1):
        for m in range(1, q):
            if gcd(m, q) == 1:
                pts.append(Fraction(m, q))
    return pts

def scan_below_mediant(k, box, Qprobe=None):
    """filtered exhaustive scan over 1=s_1<...<s_k<=box, gcd 1. Returns list of
       (M, S, tstar) with M < mediant, plus the minimal M (=sigma_2 if any below)."""
    mediant = Fraction(2, 2 * k + 1)
    floor = Fraction(1, k + 1)
    if Qprobe is None: Qprobe = 3 * k + 3
    probes = build_probe_times(k, Qprobe)
    survivors = []
    n_scanned = 0
    n_passed = 0
    for tail in itertools.combinations(range(2, box + 1), k - 1):
        S = (1,) + tail
        if not primitive(S): continue
        n_scanned += 1
        # FILTER: any probe time with min >= mediant => M>=mediant => skip
        skip = False
        for t in probes:
            mv = min(fn(v * t) for v in S)
            if mv >= mediant:
                skip = True; break
        if skip: continue
        n_passed += 1
        M, ts = maxmin_exact(S)
        if floor < M < mediant:
            survivors.append((M, S, ts))
    survivors.sort(key=lambda x: x[0])
    return survivors, n_scanned, n_passed

def express(M, k):
    p, q = M.numerator, M.denominator
    e = p * (k + 1) - q
    return (p if e == 1 else None), e, p, q

if __name__ == "__main__":
    # box must exceed q/2 ~ a*(k+1)/2 to catch level a. For a up to ~5: box ~ 2.5k.
    config = {6: 22, 7: 24, 8: 26, 9: 28, 10: 30, 11: 30, 12: 30}
    print("=== filtered below-mediant scan: a_max(k) growth test ===")
    print(f"  (catches level a if box > (a(k+1)-1)/2)\n")
    summary = []
    for k in sorted(config):
        box = config[k]
        survivors, ns, npass = scan_below_mediant(k, box)
        floor = Fraction(1, k + 1)
        max_box_a = (2 * box) // (k + 1)   # largest a catchable: q=a(k+1)-1<=2box
        if not survivors:
            print(f"  k={k:2d} box={box} (catch a<= {max_box_a}): NO below-mediant config. "
                  f"scanned={ns} passed_filter={npass}")
            summary.append((k, None, None, None, max_box_a))
            continue
        sigma2 = survivors[0][0]
        a, e, p, q = express(sigma2, k)
        amax = max((express(M, k)[0] or 0) for M, _, _ in survivors)
        g = sigma2 - floor
        # distinct a-values seen among e==1 survivors
        avals = sorted(set(express(M, k)[0] for M, _, _ in survivors if express(M, k)[0]))
        print(f"  k={k:2d} box={box} (catch a<= {max_box_a}): sigma_2={sigma2} (a={a},e={e}) "
              f"g*k^2={float(g)*k*k:.4f}  #below={len(survivors)} mediant-a seen={avals}  "
              f"witness={survivors[0][1]}")
        summary.append((k, sigma2, a, avals, max_box_a))
    print("\n=== a_max(k) and the catchable ceiling ===")
    for k, s2, a, avals, mba in summary:
        print(f"  k={k:2d}: sigma_2-a={a}  mediant-a values seen={avals}  (could detect up to a={mba})")
    print("\n  If 'a seen' stays in {2,3} while box could catch a>=4,5 => a_max BOUNDED => g=Theta(1/k^2).")
