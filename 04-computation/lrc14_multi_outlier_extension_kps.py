#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
HYP-3950 EXTENSION: close the r >= 5 outlier counts + the packed-cluster regime (r >= 7) + the
exact pair-overlap law. Companion to lrc14_multi_outlier_flat_ratio_peel_kps.py.

 (1) LEDGER m = 1..6 (min lonely measure over bounded m-cores) -> union floors r = 5, 6.
 (2) r >= 7 (union floor vacuous, 1 - r/7 <= 0): the outliers are all > V_bound; if they are SPREAD the
     multiplicative peel applies; if PACKED (a cluster {W..W+11} at large W) the lonely measure is
     governed by the 2-torus offset structure: meas -> E_s[meas(∩_j (safe0 - j*s))] ~ (6/7)^r-ish >> 1/36.
     We probe packed clusters directly (W = 200, 1000, 5000) and mixed bounded+packed configs.
 (3) EXACT PAIR-OVERLAP LAW: meas(D_w ∩ D_{w'}) depends only on the reduced ratio p/q; = 1/(7k) for
     ratio k (exact multiples), -> 1/49 for coprime (independence). Verify + tabulate f(p,q)*49.
"""
import sys, itertools, random
from fractions import Fraction as Fr
from math import gcd, floor
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
BAND = 1.0/14.0; T36 = 1.0/36.0

def clip(arcs, v, band=BAND):
    out = []
    for (lo, hi) in arcs:
        k0 = int(floor(lo*v)); k1 = int(floor(hi*v))
        for k in range(k0, k1+1):
            a = (k + band)/v; b = (k + 1 - band)/v
            l = lo if lo > a else a
            h = hi if hi < b else b
            if l < h - 1e-15: out.append((l, h))
    return out
def lonely_arcs(S):
    arcs = [(0.0, 1.0)]
    for v in sorted(S): arcs = clip(arcs, v)
    return arcs
def measN(S):
    a = lonely_arcs(S); return sum(h-l for l, h in a), len(a)

print("="*104)
print(" (1) LEDGER m=1..6 + union floors r=5,6  [floor: min_m*(1-r/7) - (N/7)*r/W  for outliers > W]")
print("="*104)
rng = random.Random(29)
led = {}
for m in range(1, 7):
    pool = set()
    for V in range(m, 16):   # V up to 15: must include drops of {1..13}/{1..14} ({1,11,12,13} = cap_9 argmin!)
        from math import comb
        d = V - m
        if comb(V, d) <= 30000:
            for drop in itertools.combinations(range(1, V+1), d):
                C = tuple(sorted(set(range(1, V+1)) - set(drop)))
                g = 0
                for c in C: g = gcd(g, c)
                pool.add(tuple(c//g for c in C))
    for _ in range(3000):
        C = tuple(sorted(rng.sample(range(1, 26), m)))
        g = 0
        for c in C: g = gcd(g, c)
        pool.add(tuple(c//g for c in C))
    best = min(((measN(C)), C) for C in pool)
    (mu, N), C = best
    led[m] = (mu, N, C)
    r = 11 - m
    uf = mu*(1.0 - r/7.0)
    print(f"  m={m}: min meas={mu:.6f} N={N:2d} argmin={C}   union floor (1-{r}/7)*min = {uf:+.6f}"
          f"   {'clears 1/36' if uf >= T36 else 'FAILS -> peel/cluster regime'}")
print("  (m=1: exactly 6/7 for any speed, scale-invariance)")
print("\n  LEDGER = CAP ATLAS cross-check (min_m should equal canon cap_{13-m}):")
from fractions import Fraction as F2
for m, capname, capval in [(2, 'cap_11', F2(66,91)), (3, 'cap_10', F2(55,91)),
                           (4, 'cap_9', F2(45,91)-F2(1,4004)), (5, 'cap_8', F2(2243,5880))]:
    mu = led[m][0]
    print(f"    m={m}: ledger {mu:.6f}  vs {capname} = {float(capval):.6f}   match: {abs(mu-float(capval))<1e-9}")
print("  => the small-m lonely ledger IS the sector route's cap atlas: min_m = cap_{13-m} (same extremizers).")
print("  Two extremal PHASES: {1} u top (Dirichlet-type) for m<=4; odd-unit-heavy near-AP for m>=5")
print("  (the phase transition = THM-576's cap-pattern break at j=5).")

print("\n" + "="*104)
print(" (2) r>=7: PACKED far clusters (union floor vacuous). meas of {W..W+r-1} alone and with bounded part")
print("="*104)
for W in (200, 1000, 5000):
    for r in (7, 9, 11):
        S = tuple(range(W, W+r))
        mu, N = measN(S)
        print(f"  cluster {{{W}..{W+r-1}}} (r={r}): meas={mu:.6f}  [(6/7)^{r}={((6/7)**r):.6f}]  ratio={mu/((6/7)**r):.3f}"
              f"   >= 1/36? {mu >= T36}")
# mixed: small bounded part + packed far cluster (the adversary's best hope at r>=7)
print("\n  mixed bounded + packed cluster (11 speeds total):")
for B in [(1,2,3,4), (1,2,3), (1,2), (1,)]:
    r = 11 - len(B)
    for W in (200, 1000):
        S = tuple(B) + tuple(range(W, W+r))
        mu, N = measN(S)
        print(f"    B={B} + {{{W}..{W+r-1}}}: meas={mu:.6f}   >= 1/36? {mu >= T36}   (x{mu/T36:.2f})")
# adversarial: packed cluster placed at a resonant scale of the bounded part + dilated clusters
print("\n  adversarial variants (resonant scale, dilated cluster, two clusters):")
for S in [ (1,2,3,4)+tuple(range(84, 91)),           # cluster at 84 = 12*7 resonance
           (1,2,3,4)+tuple(84*k for k in range(1, 8)),   # AP cluster step 84
           (1,2,3)+tuple(50*k for k in range(1, 9)),     # AP cluster step 50
           (1,2,3,4)+(50,51,100,101,200,201,400),        # paired doublets
           (1,2,3,4,5,7,8)+(19,20,39,40) ]:              # near-min 7-core + two consecutive doublets
    mu, N = measN(S)
    print(f"    S={S}: meas={mu:.6f}  >= 1/36? {mu >= T36}  (x{mu/T36:.2f})")

print("\n" + "="*104)
print(" (3) EXACT PAIR-OVERLAP LAW: meas(D_p*g ∩ D_q*g) = F(p,q) (ratio-only); F(k,1)=1/(7k); coprime->1/49")
print("="*104)
full = [(0.0, 1.0)]
def overlap(w1, w2):
    S1 = clip(full, w1); S12 = clip(S1, w2)
    return 1.0 - 2.0*(6.0/7.0) + sum(h-l for l, h in S12)
print("  F(k,1)*7k (should be 1.0 for exact multiples):")
for k in range(2, 11):
    ov = overlap(30, 30*k)
    print(f"    k={k:2d}: overlap={ov:.6f}  x 7k = {ov*7*k:.4f}")
print("\n  ratio-only test: F(p,q) at different scales g (rows should match):")
for (p, q) in [(3,2), (5,3), (7,5), (4,3), (5,2)]:
    vals = []
    for g in (12, 24, 60):
        vals.append(overlap(q*g, p*g))
    print(f"    p/q={p}/{q}: " + "  ".join(f"g={g}: {v:.6f}" for g, v in zip((12,24,60), vals)) +
          f"   x49 = {vals[-1]*49:.4f}")
print("\n  coprime spot checks (expect ~1/49 = 0.020408):")
for (a, b) in [(101, 102), (211, 305), (97, 300), (500, 501)]:
    print(f"    ({a},{b}): {overlap(a,b):.6f}")
print("DONE.")
