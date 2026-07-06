#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S22 -- working the SOLE OPEN PIECE (the density floor) collaboratively:
develop opus's theta-sum (HYP-4446) by RELATION-LENGTH SHELLS and the BEURLING-SELBERG
MAJORANT bound (the residual kps-S31 isolated).

opus's identity: safe(S,beta) = SUM_{a in L(S)} prod_i hhat(a_i), L(S)={a: sum a_i v_i=0},
  hhat(0)=1-2beta, hhat(m)=-sin(2 pi m beta)/(pi m)  (|hhat(m)| ~ 1/m).
Shells by relation length: safe = (1-2beta)^n + [harmonic (1,-2,1)] + [longer].
kps-S30: a family has ALL harmonic relations IFF it is an AP.  opus-S113: Freiman is
NECESSARY-NOT-SUFFICIENT -- the n=7 tiler {1,5,6,11,16,17} has FEWER 3-term relations
than the AP yet still tiles (safe=0), because the WIDER k=6 window absorbs the deficit
(STRUCTURE x WIDTH).

THIS SCRIPT: compute the theta-sum truncated to SPARSE SHORT relations (<=4 nonzero
coords, |a_i|<=2) for the AP, a non-AP non-tiler, and the tilers (AP + n=7 gap member),
to expose HOW MUCH the short shells cancel the main term, and the residual the Selberg
tail bound must control.  BEURLING-SELBERG framework: a band-limited MAJORANT g+ >= g
(support |m|<=N, excess 1/(N+1)) gives safe >= INT prod(1-g+) = a FINITE theta-sum;
positivity for non-AP is a bounded check -- but needs N ~ 2k^2 (excess < window 1/(2k^2)).
"""
from math import sin, pi, gcd
from itertools import combinations, product
from functools import reduce

def hhat(m, beta):
    if m == 0: return 1 - 2*beta
    return -sin(2*pi*m*beta)/(pi*m)

def theta_truncated(S, beta, max_coord=2, max_support=4):
    """sum over relations a with sum a_i v_i = 0, <=max_support nonzero coords,
       |a_i|<=max_coord, of prod hhat(a_i).  Includes the a=0 main term."""
    n = len(S)
    total = hhat(0, beta)**n     # main term
    seen = set()
    for k in range(2, max_support+1):     # k nonzero coords
        for idx in combinations(range(n), k):
            for vals in product(range(-max_coord, max_coord+1), repeat=k):
                if 0 in vals: continue
                if sum(vals[j]*S[idx[j]] for j in range(k)) != 0: continue
                key = tuple(sorted(zip(idx, vals)))
                if key in seen: continue
                seen.add(key)
                term = 1.0
                for j in range(n): term *= hhat(0, beta)
                # replace the k nonzero coords
                for j in range(k):
                    term = term / hhat(0, beta) * hhat(vals[j], beta)
                total += term
    return total, len(seen)

def n_harmonic(S):
    S=sorted(S); return sum(1 for i in range(len(S)-2) if S[i]-2*S[i+1]+S[i+2]==0)

def log(m=""): print(m, flush=True)

# n=13 (12 runners), beta=2/25
b13 = 2/25
fams13 = {
 "AP {1..12} (safe=0, tiles)": list(range(1,13)),
 "single-lift {1..11,24} (2/25 edge)": list(range(1,12))+[24],
 "generic {1..11,23} (loose, safe>0)": list(range(1,12))+[23],
 "far {1,2,3,4,5,6,7,8,9,10,11,50}": [1,2,3,4,5,6,7,8,9,10,11,50],
}
log("=== n=13, beta=2/25: theta-sum truncated to sparse short relations (|a|<=2, <=4 coords) ===")
log(f"  main term (1-2b)^12 = {(1-2*b13)**12:.5f}")
log(f"  {'family':38s} {'#harm':>6} {'trunc-theta':>12} {'#relations':>11}")
for name,S in fams13.items():
    th,nr = theta_truncated(S, b13)
    log(f"  {name:38s} {n_harmonic(S):>6} {th:>12.5f} {nr:>11}")

# n=7 (6 runners), beta=2/13
b7 = 2/13
fams7 = {
 "AP {1..6} (safe=0, tiles)": list(range(1,7)),
 "n=7 GAP MEMBER {1,5,6,11,16,17} (safe=0, tiles!)": [1,5,6,11,16,17],
 "generic n=7 {1,2,3,4,5,20} (loose)": [1,2,3,4,5,20],
}
log(f"\n=== n=7, beta=2/13: main term (1-2b)^6 = {(1-2*b7)**6:.5f} ===")
log(f"  {'family':50s} {'#harm':>6} {'trunc-theta':>12}")
for name,S in fams7.items():
    th,nr = theta_truncated(S, b7)
    log(f"  {name:50s} {n_harmonic(S):>6} {th:>12.5f}")

log("\nREADING: the AP has the most harmonic relations => most short-shell cancellation.")
log("The n=7 tiler cancels to safe=0 with FEWER harmonic relations (opus-S113 structure x width):")
log("its wider window means the LONGER shells finish the cancellation.  The truncated theta")
log("tracks 'how close to tiling'; the Selberg majorant makes the tail rigorous (N~2k^2).")
