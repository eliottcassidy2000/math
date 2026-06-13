#!/usr/bin/env python3
"""
ham_cycle_dominance.py -- kind-pasteur-2026-03-13-S60

HYPOTHESIS: Among circulant tournaments on Z_p with the same alpha_2,
H is determined by c_p (Hamiltonian cycle count).

More precisely: H = 1 + 2*N + 4*alpha_2 + 8*alpha_3
Hamiltonian cycles on all p vertices ALWAYS overlap with any other cycle
(since every cycle uses >= 3 vertices out of p). So:
- Adding a Hamiltonian cycle increases N by 1 (contributing +2 to H)
- It does NOT reduce alpha_2 (it overlaps with everything)
- It contributes 0 to alpha_3 (can't be disjoint from anything)

So Hamiltonian cycles are "pure benefit" — each one adds exactly 2 to H
with no offsetting alpha_2 reduction.

Test: within each class, is the variation entirely due to c_p?

Also: derive the "marginal H per cycle" for each length k.
"""

from itertools import combinations
from collections import defaultdict

def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i+s) % p] = 1
    return A

def count_ham_cycles(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b]*A[b][c]*A[c][a]) + (A[a][c]*A[c][b]*A[b][a])
    start = 0
    dp = {(1 << start, start): 1}
    for mask in range(1, 1 << k):
        if not (mask & (1 << start)):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            cnt = dp[key]
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + cnt
    full = (1 << k) - 1
    total = 0
    for v in range(k):
        if v == start:
            continue
        key = (full, v)
        if key in dp and dp[key] > 0:
            if A[verts[v]][verts[start]]:
                total += dp[key]
    return total


# Marginal H per additional directed cycle of length k
# A directed k-cycle on vertex set V contributes:
# +1 to N (adding 2 to H)
# To alpha_2: it forms a disjoint pair with every existing cycle C'
#   whose V(C') is disjoint from V. So it adds sum_{C' disj} n(V(C'))... wait,
#   no. alpha_2 counts pairs of CYCLES, not pairs of vertex sets.
#   Adding one new cycle C on V adds exactly |{C' : V(C') cap V = 0}| new pairs.
#   This is C_odd(T[comp(V)]) by THM-160.
# To alpha_3: it extends existing disjoint pairs to triples.

# So the "marginal H" from one new k-cycle on V is:
#   dH = 2 + 4 * |{C' : V(C') disj V}| + higher order terms
#   = 2 + 4 * C_odd(T[comp(V)])

# For Hamiltonian cycles: comp(V) = empty, so C_odd = 0.
#   dH = 2 per Hamiltonian cycle. EXACTLY.

# For k-cycles with k < p: dH = 2 + 4*D where D = C_odd(comp).
# At p=11, Paley:
#   k=3: D = 172 (constant for all 3-sets), so dH = 2 + 4*172 = 690
#   k=5: D varies from 7 to 16, so dH varies from 30 to 66
#   k=7: D varies from 0 to 2, so dH varies from 2 to 10
#   k=9: D = 0 (comp has 2 vertices, no 3-cycles), so dH = 2

# Wait -- these are MARGINAL effects of adding one cycle.
# But we're comparing DIFFERENT tournaments, not adding/removing cycles.

# Actually, the more relevant analysis:
# Consider two orientations with c_k values (c_3, c_5, ..., c_p).
# Their H values differ by:
#   dH = 2*dN + 4*d(alpha_2) + 8*d(alpha_3)
# where dN = sum dc_k.

# The key is that d(alpha_2) depends on the JOINT structure of which
# vertex sets gain/lose cycles, not just the totals dc_k.

# Let me instead compute the EFFECTIVE marginal H per unit dc_k:
# If we perturb c_k by +1 (holding all other c_{k'} fixed),
# what is the expected change in H?

# This requires knowing d(alpha_2)/d(c_k).
# From the data, we can estimate this empirically.

print("EFFECTIVE MARGINAL H PER CYCLE BY LENGTH")
print("="*60)

# At p=11: 4 classes give us 3 independent differences
# Let c = (c5, c7, c9, c11) [c3 is constant]
# dH = 2*dc5 + 2*dc7 + 2*dc9 + 2*dc11 + 4*d(a2) + 8*d(a3)
# We need to solve for the "effective coefficient" per dc_k.

# Classes (Paley-relative):
# B-A: dc = (0, -44, -374, -858, -352), da2 = +341, da3 = +33, dH = -1628
# C-A: dc = (0, -110, -561, -1705, -396), da2 = +231, da3 = +319, dH = -2068
# D-A: dc = (0, -22, -231, -781, -506), da2 = +33, da3 = +33, dH = -2684

import numpy as np

# Let's use numpy to solve the system
# H = 1 + 2*N + 4*alpha_2 + 8*alpha_3
# So dH = 2*dN + 4*da2 + 8*da3

# We have 3 data points (3 class differences relative to Paley):
# dc5, dc7, dc9, dc11, da2, da3, dH

diffs = [
    {'dc': {5: -44, 7: -374, 9: -858, 11: -352}, 'da2': 341, 'da3': 33, 'dH': -1628},
    {'dc': {5: -110, 7: -561, 9: -1705, 11: -396}, 'da2': 231, 'da3': 319, 'dH': -2068},
    {'dc': {5: -22, 7: -231, 9: -781, 11: -506}, 'da2': 33, 'da3': 33, 'dH': -2684},
]

# Verify: dH = 2*sum(dc_k) + 4*da2 + 8*da3
for d in diffs:
    dN = sum(d['dc'].values())
    predicted = 2*dN + 4*d['da2'] + 8*d['da3']
    print(f"  dH={d['dH']}, predicted={predicted}, match={predicted==d['dH']}")

# So the "effective H per unit c_k" is NOT just 2 -- it's 2 + correction from alpha terms.
# The correction depends on HOW alpha_2 and alpha_3 change when c_k changes.

# For Ham cycles (k=p=11): comp(V) is empty, so each Ham cycle adds 0 to alpha_2/alpha_3.
# Effective dH/dc_p = 2.

# For other k: the effect on alpha_2 depends on the overlap structure.
# We can estimate by looking at the partial effects.

# Try: fit linear model dH = b5*dc5 + b7*dc7 + b9*dc9 + b11*dc11
# using 3 equations, 4 unknowns (underdetermined)
# But we know b11 = 2 (pure contribution, no alpha effect)

print("\nSolving for effective coefficients with b11 = 2:")
# dH - 2*dc11 = b5*dc5 + b7*dc7 + b9*dc9
# 3 equations, 3 unknowns

A_mat = []
b_vec = []
for d in diffs:
    row = [d['dc'][5], d['dc'][7], d['dc'][9]]
    rhs = d['dH'] - 2 * d['dc'][11]
    A_mat.append(row)
    b_vec.append(rhs)

# Solve 3x3 system
A_np = [[float(x) for x in row] for row in A_mat]
b_np = [float(x) for x in b_vec]

# Cramer's rule (no numpy needed)
def det3(m):
    return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
          - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
          + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))

D = det3(A_np)
print(f"  det(A) = {D}")

if abs(D) > 1e-10:
    coeffs = []
    for col in range(3):
        M = [row[:] for row in A_np]
        for row in range(3):
            M[row][col] = b_np[row]
        coeffs.append(det3(M) / D)

    k_labels = [5, 7, 9, 11]
    b_values = coeffs + [2.0]

    print(f"\n  Effective H per unit c_k:")
    for k, b in zip(k_labels, b_values):
        print(f"    b_{k} = {b:.4f}")

    # Verify
    print(f"\n  Verification:")
    for i, d in enumerate(diffs):
        predicted = sum(b * d['dc'][k] for k, b in zip(k_labels, b_values))
        print(f"    diff {i}: actual dH={d['dH']}, predicted={predicted:.1f}")

    # Interpret: b_k > 2 means c_k cycles are "more valuable" than Ham cycles
    # b_k < 2 means they're less valuable (their overlaps reduce alpha_2)
    # b_k = 2 means they're neutral (like Ham cycles)

    print(f"\n  Interpretation:")
    for k, b in zip(k_labels, b_values):
        if b > 2.01:
            print(f"    c_{k}: b={b:.4f} > 2 -- SUPER-VALUABLE (disjoint pair bonus)")
        elif b < 1.99:
            print(f"    c_{k}: b={b:.4f} < 2 -- LESS VALUABLE (overlap penalty)")
        else:
            print(f"    c_{k}: b={b:.4f} ~ 2 -- NEUTRAL (no alpha effect)")
else:
    print("  System is singular!")

    # Check rank
    print(f"\n  Checking proportionality:")
    for i in range(3):
        for j in range(i+1, 3):
            ratios = []
            for k in range(3):
                if A_np[j][k] != 0:
                    ratios.append(A_np[i][k] / A_np[j][k])
                else:
                    ratios.append('inf')
            print(f"    row{i}/row{j}: {ratios}")

print("\nDONE.")
