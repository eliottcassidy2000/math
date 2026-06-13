#!/usr/bin/env python3
"""real_data_s115.py — Apply the Cayley-Delannoy theory to REAL tournaments"""
from math import sqrt, log, factorial, comb
from fractions import Fraction
from itertools import permutations

def cv_formula(n):
    """CV = sqrt(2/n) approximation."""
    return sqrt(2.0/n)

def E_H(n):
    """Expected Hamiltonian path count for random tournament."""
    return factorial(n) / 2**(n-1)

def H_tournament(T, n):
    """Count Hamiltonian paths in tournament T (adjacency matrix)."""
    count = 0
    for perm in permutations(range(n)):
        ok = True
        for i in range(n-1):
            if T[perm[i]][perm[i+1]] != 1:
                ok = False
                break
        if ok:
            count += 1
    return count

def ranking_significance(H, n):
    """Z-score of H(T) against random tournament null."""
    mu = E_H(n)
    cv = cv_formula(n)
    sigma = mu * cv
    z = (H - mu) / sigma
    return z

# ============================================================
# EXAMPLE 1: A 5-team round-robin (like a sports group stage)
# ============================================================
print("="*60)
print("EXAMPLE 1: FIFA-STYLE 5-TEAM GROUP")
print("="*60)
print()

# Typical group stage: clear hierarchy
# Team 0 (Brazil) > 1 (Germany) > 2 (France) > 3 (Japan) > 4 (Costa Rica)
# But with upsets: Japan beats Germany (like 2022 World Cup!)
T_group = [
    [0, 1, 1, 1, 1],  # Brazil beats everyone
    [0, 0, 1, 0, 1],  # Germany beats France, Costa Rica, loses to Japan
    [0, 0, 0, 1, 1],  # France beats Japan, Costa Rica
    [0, 1, 0, 0, 1],  # Japan beats Germany, Costa Rica
    [0, 0, 0, 0, 0],  # Costa Rica loses all
]
n = 5
H = H_tournament(T_group, n)
z = ranking_significance(H, n)
mu = E_H(n)

print(f"Teams: Brazil > France > Japan > Germany > Costa Rica (with upset)")
print(f"H(T) = {H} Hamiltonian paths")
print(f"E[H] = {mu:.1f} for random 5-team tournament")
print(f"CV = sqrt(2/5) = {cv_formula(5):.3f}")
print(f"Z-score = {z:.2f}")
print(f"Interpretation: {'Significantly ordered' if abs(z) > 2 else 'Within random noise' if abs(z) < 1 else 'Moderately ordered'}")
print()

# Compare: transitive tournament (perfect hierarchy)
T_trans = [[1 if j > i else 0 for j in range(5)] for i in range(5)]
H_trans = H_tournament(T_trans, 5)
z_trans = ranking_significance(H_trans, 5)
print(f"Perfect hierarchy: H={H_trans}, Z={z_trans:.2f}")

# Full cycle (maximum chaos)
T_cycle = [[0]*5 for _ in range(5)]
for i in range(5):
    T_cycle[i][(i+1)%5] = 1
    T_cycle[i][(i+2)%5] = 1
H_cycle = H_tournament(T_cycle, 5)
z_cycle = ranking_significance(H_cycle, 5)
print(f"Cyclic tournament: H={H_cycle}, Z={z_cycle:.2f}")
print()

# ============================================================
# EXAMPLE 2: Chess-style 6-player round-robin
# ============================================================
print("="*60)
print("EXAMPLE 2: CHESS 6-PLAYER CANDIDATES")
print("="*60)
print()

# Candidates tournament style: top 2 clearly better, rest competitive
T_chess = [
    [0, 1, 1, 1, 1, 1],  # Player 0: dominant (5-0)
    [0, 0, 1, 1, 1, 1],  # Player 1: strong (4-1)
    [0, 0, 0, 1, 0, 1],  # Player 2: middle (2-3)
    [0, 0, 0, 0, 1, 1],  # Player 3: middle (2-3)
    [0, 0, 1, 0, 0, 1],  # Player 4: middle (2-3)
    [0, 0, 0, 0, 0, 0],  # Player 5: tail (0-5)
]
n = 6
H = H_tournament(T_chess, n)
z = ranking_significance(H, n)
print(f"Strong hierarchy with competitive middle: H={H}, Z={z:.2f}")
print(f"E[H] = {E_H(6):.1f}, CV = {cv_formula(6):.3f}")
print()

# More competitive version
T_chess2 = [
    [0, 1, 1, 0, 1, 1],  # 4-1
    [0, 0, 1, 1, 0, 1],  # 3-2
    [0, 0, 0, 1, 1, 1],  # 3-2
    [1, 0, 0, 0, 1, 1],  # 3-2
    [0, 1, 0, 0, 0, 1],  # 2-3
    [0, 0, 0, 0, 0, 0],  # 0-5
]
H2 = H_tournament(T_chess2, n)
z2 = ranking_significance(H2, n)
print(f"Competitive with upsets: H={H2}, Z={z2:.2f}")
print()

# ============================================================
# EXAMPLE 3: How many teams before rankings are meaningless?
# ============================================================
print("="*60)
print("EXAMPLE 3: WHEN DO RANKINGS BECOME MEANINGFUL?")
print("="*60)
print()
print("For a ranking to be 'real' (Z > 2), the tournament must have")
print("H(T) > E[H] * (1 + 2*CV) = E[H] * (1 + 2*sqrt(2/n)).")
print()
print("  n    E[H]          CV      Threshold H/E[H]")
print("  " + "-"*50)
for n in [4, 5, 6, 8, 10, 16, 20, 32]:
    mu = E_H(n)
    cv = cv_formula(n)
    threshold = 1 + 2*cv
    print(f"  {n:3d}   {mu:12.1f}   {cv:.4f}   {threshold:.4f}")

print()
print("For n=32 teams: any H within 50% of E[H] is noise.")
print("Only H > 1.50*E[H] indicates real skill differences.")
print()

# ============================================================
# EXAMPLE 4: The 'league competitiveness index'
# ============================================================
print("="*60)
print("EXAMPLE 4: LEAGUE COMPETITIVENESS INDEX")
print("="*60)
print()
print("L(T) = H(T)/E[H]: normalized Hamiltonian path count.")
print("  L >> 1: clear hierarchy (easy to rank)")
print("  L ~ 1: competitive (rankings are noise)")
print("  L << 1: paradoxical (many contradictions)")
print()

# Generate some example tournaments for n=5
import random
random.seed(42)

print("Random 5-team tournaments (1000 samples):")
L_values = []
for _ in range(1000):
    T = [[0]*5 for _ in range(5)]
    for i in range(5):
        for j in range(i+1, 5):
            if random.random() < 0.5:
                T[i][j] = 1
            else:
                T[j][i] = 1
    H = H_tournament(T, 5)
    L = H / E_H(5)
    L_values.append(L)

L_values.sort()
print(f"  Mean L: {sum(L_values)/len(L_values):.3f} (should be ~1)")
print(f"  Std L:  {(sum((l-1)**2 for l in L_values)/len(L_values))**0.5:.3f} (should be ~{cv_formula(5):.3f})")
print(f"  Min L:  {min(L_values):.3f}")
print(f"  Max L:  {max(L_values):.3f}")
print(f"  P(L > 1.5): {sum(1 for l in L_values if l > 1.5)/len(L_values)*100:.1f}%")
print(f"  P(L > 2.0): {sum(1 for l in L_values if l > 2.0)/len(L_values)*100:.1f}%")

# ============================================================
# EXAMPLE 5: Recommendation systems
# ============================================================
print()
print("="*60)
print("EXAMPLE 5: PRODUCT RANKING (A/B TESTING)")
print("="*60)
print()
print("7 products compared pairwise by users.")
print("Each pair: which product is preferred?")
print("The tournament T encodes all pairwise preferences.")
print()
print(f"E[H] = {E_H(7):.1f} for n=7")
print(f"CV = {cv_formula(7):.3f}")
print(f"sigma = {E_H(7)*cv_formula(7):.1f}")
print()
print("If your A/B test gives H = 189 (maximum possible):")
print(f"  Z = {ranking_significance(189, 7):.2f} — extremely significant ranking")
print("If your A/B test gives H = 79 (near random):")
print(f"  Z = {ranking_significance(79, 7):.2f} — no significant ranking")
print("If your A/B test gives H = 120 (moderate):")
print(f"  Z = {ranking_significance(120, 7):.2f} — moderate evidence of ranking")
