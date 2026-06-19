#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC max-min spectrum: the SECOND point sigma_2(k) and the spectral gap
    g(k) = sigma_2(k) - 1/(k+1).
kind-pasteur-2026-06-19-S9.

Setup (user's k = #speeds; codex/oracle's n = #runners = k+1):
  S = set of k distinct gcd-1 positive integer speeds, observer at 0.
  M(S) = max_tau min_{v in S} ||v tau||   (the max-min collar / gap).
  LRC(k+1): M(S) >= 1/(k+1) for every such S.  Floor 1/(k+1) hit by AP {1..k}.

Claims under test:
  (A) M_k = {1,...,k-1, 2k} (AP with apex doubled) gives M = 2/(2k+1) = mediant(1/(k+1),1/k).
  (B) sigma_2(k) = min over non-tight S of M(S).  Is 1/(k+1) < sigma_2(k) <= 2/(2k+1)?
      mediant-tight (=) for k<=5, strict (<) for k>=6 ?
  (C) g(k) = sigma_2(k) - 1/(k+1).  Upper bound g(k) <= 1/((2k+1)(k+1)) = Theta(1/k^2).
      LIVE QUESTION: lower bound. Does g(k)*k^2 stay bounded below, or dip to 0?

EXACT max-min via crossing candidates: f(tau)=min_v||v tau|| is piecewise linear;
its local maxima sit where two sawtooths tie: v_i tau == +- v_j tau (mod 1),
i.e. tau = m/(v_i +- v_j).  Enumerate all such rationals in (0,1), evaluate f
exactly (Fraction), take the max.  Provably exact for finite speed sets.
"""
import sys, itertools
from fractions import Fraction
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None

def frac_norm(x):
    """||x|| = distance from x (a Fraction) to nearest integer, as a Fraction in [0,1/2]."""
    r = x - (x.numerator // x.denominator)   # in [0,1)
    return min(r, 1 - r)

def maxmin(S):
    """Exact M(S) = max_tau min_{v in S} ||v tau||.  Returns (Fraction M, Fraction tau*)."""
    S = sorted(set(S))
    cand = set()
    for i in range(len(S)):
        for j in range(i, len(S)):
            for d in (S[i] + S[j], S[i] - S[j]):
                if d == 0:
                    continue
                d = abs(d)
                for m in range(1, d):
                    cand.add(Fraction(m, d))
    best = Fraction(0); bestt = Fraction(0)
    for t in cand:
        mv = min(frac_norm(v * t) for v in S)
        if mv > best:
            best = mv; bestt = t
    return best, bestt

def is_primitive(S):
    from functools import reduce
    return reduce(gcd, S) == 1

def spectrum_scan(k, max_speed, want_floor=None):
    """Scan all k-subsets of {1..max_speed} (gcd-1), return sorted list of (M, S)
       grouped, plus the floor and the second-smallest distinct value sigma_2."""
    floor = Fraction(1, k + 1)
    vals = {}   # Fraction M -> example S (smallest sum/lex)
    count_below_mediant = 0
    mediant = Fraction(2, 2 * k + 1)
    below_examples = []
    for S in itertools.combinations(range(1, max_speed + 1), k):
        if S[0] != 1:            # WLOG 1 in S (dilation/normalization) -- speeds gcd 1, and
            continue             # any primitive set is a dilation-free rep; require min speed 1
        if not is_primitive(S):
            continue
        M, t = maxmin(S)
        if M not in vals:
            vals[M] = S
        if floor < M < mediant:
            count_below_mediant += 1
            if len(below_examples) < 12:
                below_examples.append((M, S))
    sorted_vals = sorted(vals.keys())
    sigma1 = sorted_vals[0]
    sigma2 = sorted_vals[1] if len(sorted_vals) > 1 else None
    return dict(floor=floor, mediant=mediant, sigma1=sigma1, sigma2=sigma2,
                sigma2_witness=vals.get(sigma2), n_below_mediant=count_below_mediant,
                below_examples=below_examples, all_vals=sorted_vals[:8],
                vals_map=vals)

if __name__ == "__main__":
    print("=== (A) verify M_k = {1,...,k-1,2k} gives 2/(2k+1) ===")
    for k in range(2, 16):
        Mk = list(range(1, k)) + [2 * k]
        M, t = maxmin(Mk)
        target = Fraction(2, 2 * k + 1)
        floor = Fraction(1, k + 1)
        ok = "OK" if M == target else f"!! got {M}"
        print(f"  k={k:2d}  M_k={Mk if k<=8 else '{1..%d,%d}'%(k-1,2*k)}  M={M}  target 2/(2k+1)={target}  [{ok}]  (floor 1/(k+1)={floor})")

    print("\n=== verify AP {1..k} is tight (M = 1/(k+1)) ===")
    for k in range(2, 16):
        AP = list(range(1, k + 1))
        M, t = maxmin(AP)
        floor = Fraction(1, k + 1)
        print(f"  k={k:2d}  M(AP)={M}  floor={floor}  [{'OK' if M==floor else '!!'}]  at tau*={t}")

    print("\n=== (B,C) spectrum scan: sigma_2(k), gap g(k), g(k)*k^2 ===")
    # box sizes: keep combinations tractable. max_speed chosen so C(max_speed-1, k-1) reasonable.
    boxes = {2: 26, 3: 24, 4: 22, 5: 20, 6: 19, 7: 18, 8: 17}
    print(f"  {'k':>2} {'floor':>8} {'sigma_2':>10} {'mediant':>9} {'g(k)':>10} {'g*k^2':>8} {'#in(fl,med)':>11}  sigma_2 witness")
    rows = []
    for k in range(2, 9):
        ms = boxes[k]
        r = spectrum_scan(k, ms)
        g = r['sigma2'] - r['floor'] if r['sigma2'] else None
        gk2 = float(g) * k * k if g else None
        rows.append((k, r, g, gk2))
        wit = r['sigma2_witness']
        print(f"  {k:>2} {str(r['floor']):>8} {str(r['sigma2']):>10} {str(r['mediant']):>9} "
              f"{str(g):>10} {gk2:>8.4f} {r['n_below_mediant']:>11}  {wit}")
    print("\n  (sigma_2 < mediant  =>  strict; sigma_2 == mediant => mediant-tight)")
    print("\n=== below-mediant examples (the open-gap configs) ===")
    for k in range(6, 9):
        r = spectrum_scan(k, boxes[k])
        print(f"  k={k}: floor={r['floor']} mediant={r['mediant']}")
        for M, S in r['below_examples'][:6]:
            print(f"      M={M} ({float(M):.5f})  S={S}")
