#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Independent adversarial verification of the "gcd-3 stratum carries the danger" REFUTATION.
Own exact code (Fraction). Re-derive dist_p / L_y from scratch (do not import).

Claims under test:
 (A) Every AP/near-AP/d=2-GAP competitor at k=8,9,10 has elements in 0..k-1 < 27,
     so the relation lattice has ZERO TRUE antipodal {a,27-a} relations.
 (B) Residue-preserving shifts (+27,+81) collapse L_y by ~0.2 regardless of stratum.
 (C) CONTROLLED test (k=9, fix removal of central element 4, vary stranger stratum)
     gives FLAT L_y drop in [-0.224,-0.201], gcd-3 interleaving unit -> no separation.
 (D) Caps: cap_8=0.38153, cap_9=0.49426, cap_10=0.6044; consec values.
"""
from fractions import Fraction
from math import gcd
from functools import reduce

C_MOD = 27

def N_at(E, x):
    hit = set()
    for e in E:
        v = e * x
        v = v - (v.numerator // v.denominator)
        s = (v.numerator * 7) // v.denominator
        hit.add(s)
    return sum(1 for j in range(1, 7) if j not in hit)

def dist_p(E):
    E = sorted(set(E))
    bps = {Fraction(0), Fraction(1)}
    for e in E:
        if e == 0:
            continue
        for a in range(0, 7 * e + 1):
            bps.add(Fraction(a, 7 * e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p = [Fraction(0)] * 7
    for i in range(len(bps) - 1):
        lo, hi = bps[i], bps[i + 1]
        if hi == lo:
            continue
        mid = (lo + hi) / 2
        t = N_at(E, mid)
        p[t] += (hi - lo)
    return p

def g_poly(k):
    g = []
    for t in range(7):
        if k == 8:
            g.append(Fraction((t-1)*(t-2)*(t-4)*(t-5), 40))
        elif k in (9, 10):
            g.append(Fraction(-(t-2)*(t-3)*(t-6), 36))
        else:
            g.append(Fraction((t-3)*(t-4), 12))
    return g

def L_y(E, k):
    p = dist_p(E)
    g = g_poly(k)
    return sum(p[t] * g[t] for t in range(7))

def stratum(a):
    """gcd of a with 27=3^3: 1->unit, 3->gcd-3, 9->gcd-9, 0/27->the 0 class."""
    r = a % C_MOD
    if r == 0:
        return 'zero'
    gg = gcd(r, C_MOD)
    return {1: 'unit', 3: 'gcd3', 9: 'gcd9'}[gg]

print("=== CAPS / CONSEC baseline (claim D) ===")
caps = {8: 0.38153, 9: 0.49426, 10: 0.6044}
consec_L = {}
for k in [8, 9, 10]:
    L = L_y(list(range(k)), k)
    consec_L[k] = L
    print(f"  k={k}: L_y(consec)={L} = {float(L):.5f}   cap={caps[k]}  margin={caps[k]-float(L):+.5f}")

print("\n=== CLAIM A: antipodal {a,27-a} relations among AP/GAP competitors ===")
# A true antipodal relation needs two distinct elements a, b in E with a+b ≡ 0 mod 27 (b=27-a).
# Check membership-level: for elements all in 0..k-1<27, the only way is a + (27-a) but 27-a >= 18 > k-1.
fams = {
    8:  [list(range(8)), [0,1,2,3,4,5,6,8], [0,2,4,6,8,10,12,14], [0,1,2,3,5,6,7,8]],
    9:  [list(range(9)), [0,1,2,3,4,5,6,7,9], [0,2,4,6,8,10,12,14,16]],
    10: [list(range(10)), [0,1,2,3,4,5,6,7,8,10], [0,2,4,6,8,10,12,14,16,18]],
}
for k, lst in fams.items():
    for E in lst:
        Em = [e % C_MOD for e in E]
        antip = sum(1 for a in Em for b in Em if a < b and (a + b) % C_MOD == 0)
        print(f"  k={k} E={E}: residues={Em} antipodal-pairs={antip}")

print("\n=== CLAIM B: residue-preserving shift (+27,+81) collapses L_y regardless of stratum ===")
# k=9 consec; push one element by +27 or +81 (preserves its residue/stratum) and see L_y.
base = list(range(9))
for push in [27, 81]:
    print(f"  push +{push} on k=9 consec, replacing element idx i with base[i]+push:")
    for i in range(9):
        E = base[:i] + [base[i] + push] + base[i+1:]
        L = L_y(E, 9)
        st = stratum(base[i])  # residue unchanged
        print(f"    move elt {base[i]} (stratum {st}) -> {base[i]+push}: L_y={float(L):.5f}  dL={float(L-consec_L[9]):+.5f}")

print("\n=== CLAIM C: CONTROLLED test (k=9, remove central element 4, vary stranger stratum) ===")
# Base remove 4 from consec -> {0,1,2,3,5,6,7,8}, then add a stranger s in 28..39, classify stratum.
core = [0,1,2,3,5,6,7,8]
print(f"  core (consec minus 4) = {core}")
results = {'unit': [], 'gcd3': [], 'gcd9': []}
for s in range(28, 40):
    E = sorted(core + [s])
    L = L_y(E, 9)
    dL = float(L - consec_L[9])
    st = stratum(s)
    if st in results:
        results[st].append((s, dL))
    print(f"    stranger {s} (stratum {st}): L_y={float(L):.5f}  dL={dL:+.5f}")
print("\n  Summary by stratum:")
for st, vals in results.items():
    if vals:
        dls = [d for _, d in vals]
        print(f"    {st}: dL range [{min(dls):+.5f}, {max(dls):+.5f}]  values={[f'{d:+.4f}' for _,d in vals]}")

print("\n=== COUNTEREXAMPLE HUNT: any E (k=8,9,10) beating cap? gcd-3 strangers near-AP ===")
# Try gcd-3 strangers placed to mimic AP corrections; check none exceed cap.
best = {8: (None, -1), 9: (None, -1), 10: (None, -1)}
import itertools
for k in [8, 9, 10]:
    # consec with one element bumped to a gcd-3 residue value (small spread, near-AP)
    cset = list(range(k))
    cands = []
    for i in range(1, k):
        for delta in range(1, 8):
            E = sorted(set(cset[:i] + [cset[i] + delta] + cset[i+1:]))
            if len(E) == k and reduce(gcd, E) == 1:
                cands.append(E)
    for E in cands:
        L = float(L_y(E, k))
        if L > best[k][1]:
            best[k] = (E, L)
    bE, bL = best[k]
    print(f"  k={k}: max over near-consec perturbations L_y={bL:.5f} at {bE}  (cap {caps[k]}) {'EXCEEDS CAP!' if bL>caps[k] else 'below cap'} {'EXCEEDS CONSEC!' if bL>float(consec_L[k]) else ''}")
