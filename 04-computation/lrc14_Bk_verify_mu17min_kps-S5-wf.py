#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_Bk_verify_mu17min_kps-S5-wf.py   (kind-pasteur-2026-06-18-S5, ADVERSARIAL)

THE LOAD-BEARING CLAIM: for k=8..13, CONSECUTIVE E={0,..,k-1} minimizes mu_{1/7}(E)
(over 0 in E, |E|=k, gcd=1). The whole k>=8 global-witness floor rests on this.

We HUNT for any E with mu_{1/7}(E) < mu_{1/7}(consecutive_k):
  (i)  EXHAUSTIVE over all E with spread s in [k-1 .. S_exh] where C(s-1,k-2) is feasible.
  (ii) STRUCTURED adversaries: perforated APs, common-factor / subtorus shapes
       (E = c*A union B), Sidon-ish, geometric, two-scale {0..t} u {M,...}.
  (iii) LARGE random search at many spreads (including spread >> 2k).
A SINGLE hit with mu < consec  =>  the 1/7 spread bound is FALSE.

Also: print the consecutive mu_{1/7}(k) and cross-check vs the quoted values
   k=8..13: 691/735, 247/294, 38/49, 1381/2205, 13823/24255, 477/1078.
"""
import sys, itertools, random
from fractions import Fraction as F
from math import gcd, comb
from functools import reduce
try: sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception: pass
random.seed(20260618)

ONE7 = F(1,7)

def collision_breakpoints(E):
    bp = {F(0), F(1)}; Es = sorted(set(E))
    for i in range(len(Es)):
        for j in range(i+1, len(Es)):
            d = Es[j]-Es[i]
            for m in range(0, d+1): bp.add(F(m, d))
    return sorted(b for b in bp if F(0) <= b <= F(1))

def cell_cyclic_gaps(E, a, b):
    mid = (a+b)/2; Es = sorted(set(E)); pts = []
    for e in Es:
        val = e*mid; n = val - (val % 1)
        pts.append((val - n, e, n))
    pts.sort(key=lambda t: t[0]); m = len(pts); gaps = []
    for i in range(m):
        (_, e_i, n_i) = pts[i]
        if i < m-1:
            (_, e_j, n_j) = pts[i+1]; alpha = e_j - e_i; beta = n_i - n_j
        else:
            (_, e0, n0) = pts[0]; alpha = e0 - e_i; beta = n_i - n0 + 1
        gaps.append((F(alpha), F(beta)))
    return gaps

def union_gt(gaps, a, b, theta):
    ivs = []
    for alpha, beta in gaps:
        if alpha == 0:
            if beta > theta: ivs.append((a, b))
        else:
            xstar = (theta - beta)/alpha
            if alpha > 0: lo, hi = max(a, xstar), b
            else:         lo, hi = a, min(b, xstar)
            if lo < hi: ivs.append((lo, hi))
    if not ivs: return F(0)
    ivs.sort(); tot = F(0); clo, chi = ivs[0]
    for lo, hi in ivs[1:]:
        if lo <= chi: chi = max(chi, hi)
        else: tot += chi-clo; clo, chi = lo, hi
    tot += chi-clo; return tot

def mu_theta(E, theta=ONE7):
    E = sorted(set(E))
    if len(E) == 1: return F(1)
    bps = collision_breakpoints(E); tot = F(0)
    for a, b in zip(bps, bps[1:]):
        if a == b: continue
        tot += union_gt(cell_cyclic_gaps(E, a, b), a, b, theta)
    return tot

def gcd1(E): return reduce(gcd, E) == 1

quoted = {8:F(691,735),9:F(247,294),10:F(38,49),11:F(1381,2205),12:F(13823,24255),13:F(477,1078)}

print("="*88)
print("STEP 1: consecutive mu_1/7(k), k=8..13, vs quoted")
print("="*88)
consec = {}
for k in range(8, 14):
    m = mu_theta(list(range(k)), ONE7); consec[k] = m
    q = quoted.get(k); tag = f"quoted={q} {'OK' if q==m else 'MISMATCH'}" if q else ""
    print(f"  k={k:2d}: mu_1/7(consec) = {m} = {float(m):.6f}   {tag}")

print()
print("="*88)
print("STEP 2: EXHAUSTIVE hunt for mu_1/7(E) < consec, over feasible spreads")
print("="*88)
EXH_CAP = 350000
overall_min_holds = True
for k in range(8, 14):
    base = consec[k]; best = base; bestE = tuple(range(k)); checked = 0
    smax = k - 1
    # grow spread while exhaustive feasible
    s = k - 1
    while comb(s-1, k-2) <= EXH_CAP:
        smax = s; s += 1
    for s in range(k-1, smax+1):
        for interior in itertools.combinations(range(1, s), k-2):
            E = (0,) + interior + (s,)
            if not gcd1(E): continue
            checked += 1
            m = mu_theta(list(E), ONE7)
            if m < best: best = m; bestE = E
    status = "consec IS min" if best == base else f"VIOLATION min={best} at {bestE}"
    if best < base: overall_min_holds = False
    print(f"  k={k:2d}: exhaustive spread<= {smax}  ({checked} sets)  -> {status}")
print(f"  ==> exhaustive (feasible spreads): consecutive-min-holds = {overall_min_holds}")

print()
print("="*88)
print("STEP 3: STRUCTURED + LARGE-SPREAD adversaries (the dangerous regime for any 'min' claim)")
print("="*88)
def hunt_structured(k, base):
    found = []
    cands = set()
    # perforated APs of larger spread: drop j elements from {0..k-1+drop}
    for extra in range(1, 9):
        full = list(range(k+extra))
        for drops in itertools.combinations(range(1, k+extra), extra):
            E = tuple(x for x in full if x not in drops)
            if len(E)==k and E[0]==0: cands.add(E)
    # common-factor / subtorus: E = c*A (then gcd reduces; but mixed c*A u B)
    for c in (2,3,4,5,7):
        for extra in range(0, 4):
            A = list(range(0, (k-extra)))
            E = tuple(sorted(set([c*a for a in A] + list(range(1, extra+1)))))
            if len(E)==k and 0 in E: cands.add(E)
    # two-scale {0..t} u {M..}
    for t in range(2, k):
        rem = k - (t+1)
        if rem < 1: continue
        for M in (k, k+2, 2*k, 3*k, 5*k, 10*k, 50*k):
            tail = tuple(range(M, M+rem))
            E = tuple(range(0, t+1)) + tail
            if len(set(E))==k: cands.add(tuple(sorted(set(E))))
    # geometric-ish / Sidon-ish
    for _ in range(30000):
        sp = random.choice([k-1,k,k+1,k+2,k+3,k+5,k+8,2*k,3*k,5*k,10*k])
        if sp < k-1: continue
        body = random.sample(range(1, sp+1), min(k-1, sp))
        E = tuple(sorted(set([0]+body)))
        if len(E)==k: cands.add(E)
    for E in cands:
        if len(set(E))!=k or 0 not in E: continue
        m = mu_theta(list(E), ONE7)
        if m < base: found.append((m, E))
    return found

struct_holds = True
for k in range(8, 14):
    base = consec[k]
    found = hunt_structured(k, base)
    if found:
        struct_holds = False
        found.sort()
        print(f"  k={k:2d}: *** {len(found)} VIOLATIONS *** lowest: mu={found[0][0]}={float(found[0][0]):.5f} at {found[0][1]}")
    else:
        print(f"  k={k:2d}: no structured/large-spread E beats consec ({float(base):.5f})")
print(f"  ==> structured+large hunt: consecutive-min-holds = {struct_holds}")

print()
print("="*88)
print("VERDICT (1/7 spread bound): consecutive minimizes mu_1/7 for k=8..13?")
print(f"   exhaustive-feasible: {overall_min_holds} ;  structured/large: {struct_holds}")
print(f"   => 1/7-spread-bound consecutive-is-min: {overall_min_holds and struct_holds}")
print("="*88)
