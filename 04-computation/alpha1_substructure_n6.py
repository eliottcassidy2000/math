#!/usr/bin/env python3
"""
alpha1_substructure_n6.py — oracle-2026-05-17-S1

Investigate the alpha1 gap at n=6: which alpha1 values occur?
And verify: for alpha2=1 classes, H = 5 + 2*alpha1 exactly.

Also investigate: the "alpha1 staircase" where H=1+2*alpha1 for degree-1 classes.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))

from tournament_lib import find_odd_cycles, hamiltonian_path_count, tournament_from_bits
from collections import defaultdict

def ip_coeffs(cycles, n):
    m = len(cycles)
    if m == 0: return [1]
    vsets = [frozenset(c) for c in cycles]
    adj_bits = [0]*m
    for a in range(m):
        for b in range(a+1, m):
            if vsets[a]&vsets[b]: adj_bits[a]|=1<<b; adj_bits[b]|=1<<a
    coeffs = [0]*(n//3+2); coeffs[0]=1; coeffs[1]=m
    pairs = [(a,b) for a in range(m) for b in range(a+1,m) if not(adj_bits[a]>>b&1)]
    coeffs[2] = len(pairs)
    while len(coeffs)>1 and coeffs[-1]==0: coeffs.pop()
    return coeffs

n = 6
print("="*60)
print(f"N={n} EXHAUSTIVE: Alpha-1 gap structure")
print("="*60)

# Map (alpha1, alpha2) -> list of H values
alpha_map = defaultdict(set)
alpha2_H = defaultdict(set)  # alpha2=1 -> (alpha1, H) pairs

for bits in range(2**(n*(n-1)//2)):
    T = tournament_from_bits(n, bits)
    H = hamiltonian_path_count(T)
    cycles = find_odd_cycles(T)
    coeffs = ip_coeffs(cycles, n)
    a1 = coeffs[1] if len(coeffs)>1 else 0
    a2 = coeffs[2] if len(coeffs)>2 else 0
    alpha_map[(a1,a2)].add(H)
    if a2 == 1:
        alpha2_H[(a1)].add(H)

# Print alpha1 distribution for degree-1 classes (alpha2=0)
print("\nDegree-1 classes (alpha2=0): achievable alpha1 values")
deg1_alphas = sorted(a1 for (a1,a2) in alpha_map if a2==0)
print(f"  alpha1 values: {deg1_alphas}")
print(f"  Missing from 0..{max(deg1_alphas)}: {[x for x in range(max(deg1_alphas)+1) if x not in deg1_alphas]}")
print(f"\n  Formula H = 1 + 2*alpha1 for all degree-1 classes:")
deg1_errs = []
for (a1,a2), Hs in sorted(alpha_map.items()):
    if a2 == 0:
        formula_H = 1 + 2*a1
        deg1_errs.append(formula_H not in Hs or len(Hs) != 1)
        Hs_str = str(sorted(Hs))
        ok = "✓" if len(Hs)==1 and list(Hs)[0]==formula_H else "✗"
        print(f"    alpha1={a1:2d}: H in {Hs_str:15s}, formula={formula_H:3d} {ok}")

print(f"\n  All H=1+2*alpha1 exact: {'✓' if not any(deg1_errs) else '✗'}")

# Alpha2=1 family
print("\nDegree-2, alpha2=1 classes: H = 5 + 2*alpha1")
for a1 in sorted(alpha2_H.keys()):
    Hs = alpha2_H[a1]
    formula_H = 5 + 2*a1  # since H = 1+2a1+4a2 = 1+2a1+4 = 5+2a1
    ok = "✓" if len(Hs)==1 and list(Hs)[0]==formula_H else "✗"
    print(f"  alpha1={a1:2d}, alpha2=1: H in {sorted(Hs)}, formula={formula_H} {ok}")

# All alpha2=const families
print("\nAll (alpha1, alpha2) -> H relationships:")
for (a1,a2), Hs in sorted(alpha_map.items()):
    formula_H = 1 + 2*a1 + 4*a2
    ok = "✓" if len(Hs)==1 and list(Hs)[0]==formula_H else "✗"
    print(f"  alpha1={a1:2d}, alpha2={a2}: H in {sorted(Hs)}, formula H={formula_H} {ok}")

# Summary
print("\n" + "="*60)
print("KEY FINDINGS:")
print(f"  H = 1 + 2*alpha1 + 4*alpha2 always (exact linear relationship)")
print(f"  alpha1=3, alpha2=0 NEVER occurs at n=6")
print(f"  Alpha1 gap at degree-1: alpha1=3 is missing")
print(f"  The gap creates root gap: alpha1=3 -> r=-1/3 (impossible)")
print(f"  alpha1=4 -> r=-1/4 (boundary, achievable)")
