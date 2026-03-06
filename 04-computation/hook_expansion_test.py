#!/usr/bin/env python3
"""
Hook Expansion & Symmetric Function Analysis
================================================
Test the Irving-Omar hook expansion [s_{(i,1^{n-i})}]U_D for tournaments.

For a tournament T on [n]:
- [s_{(1^n)}]U_T = ham(T) = H(T) (Hamiltonian paths)
- [s_{(n)}]U_T = ham(T^op) = H(T) (by path reversal)
- For hooks (i,1^{n-i}): counts permutations with D-descent set {i,...,n-1}

Also: the "e-h" duality. For tournaments, U_T in the power sum basis:
  U_T = sum_{sigma in S(T), all cycles odd} 2^{psi(sigma)} * p_{type(sigma)}
where p_{type} is the power sum indexed by cycle type.

This gives: [p_{(n)}]U_T = 2 * c_n(T) (directed n-cycles, counted twice for forward/backward)
            [p_{(1^n)}]U_T = H(T)

kind-pasteur-2026-03-05-S16
"""

import sys
sys.path.insert(0, r'C:\Users\Eliott\Documents\GitHub\math\03-artifacts\code')
from tournament_lib import (tournament_from_bits, all_tournaments,
                             hamiltonian_path_count, find_odd_cycles,
                             opposite_tournament)
from itertools import permutations
from math import factorial

def count_perms_with_descent_set(T, desc_set):
    """Count permutations sigma in S_n whose T-descent set equals desc_set.
    A T-descent at position i means T[sigma[i]][sigma[i+1]] = 0 (i.e., sigma[i+1] -> sigma[i])."""
    n = len(T)
    count = 0
    for perm in permutations(range(n)):
        actual_desc = set()
        for i in range(n-1):
            if T[perm[i+1]][perm[i]]:  # sigma[i+1] -> sigma[i] means descent
                actual_desc.add(i+1)  # 1-indexed position
        if actual_desc == desc_set:
            count += 1
    return count

def hook_coefficients(T):
    """Compute [s_{(i,1^{n-i})}]U_T for i=1,...,n.
    By Prop 26: this is the count of perms with D-descent set {i,...,n-1}."""
    n = len(T)
    coeffs = []
    for i in range(1, n+1):
        desc_set = set(range(i, n))  # {i, i+1, ..., n-1}
        c = count_perms_with_descent_set(T, desc_set)
        coeffs.append(c)
    return coeffs

# ── n=3 ──
print("=" * 60)
print("HOOK EXPANSION ANALYSIS")
print("=" * 60)

# All n=3 tournaments
for bits in range(1 << 3):
    T = tournament_from_bits(3, bits)
    h = hamiltonian_path_count(T)
    hooks = hook_coefficients(T)
    cycles = find_odd_cycles(T)
    nc = len(cycles)
    print(f"bits={bits:03b}: H={h}, hooks={hooks}, [s_(1,1,1)]={hooks[0]}, [s_(3)]={hooks[2]}, cycles={nc}")

print(f"\nVerification: [s_(1^n)] should equal H(T), [s_(n)] should equal H(T^op) = H(T)")

# ── n=4 ──
print(f"\n{'='*60}")
print("n=4: All tournaments")
print("=" * 60)

hook_sums = [0, 0, 0, 0]
for bits in range(1 << 6):
    T = tournament_from_bits(4, bits)
    hooks = hook_coefficients(T)
    for i in range(4):
        hook_sums[i] += hooks[i]

print(f"Sum of hooks over all 64 tournaments: {hook_sums}")
print(f"Average hooks: {[h/64 for h in hook_sums]}")

# Check specific interesting tournaments
print("\nDetailed n=4:")
for bits in range(1 << 6):
    T = tournament_from_bits(4, bits)
    h = hamiltonian_path_count(T)
    if h in [1, 3, 5]:  # representative values
        hooks = hook_coefficients(T)
        print(f"  bits={bits:06b}: H={h}, hooks={hooks}")
        break  # just one example per H

# ── n=5: sample ──
print(f"\n{'='*60}")
print("n=5: 10 representative tournaments")
print("=" * 60)

import random
random.seed(42)

# Transitive
T = tournament_from_bits(5, 0)
hooks = hook_coefficients(T)
print(f"Transitive: H={hamiltonian_path_count(T)}, hooks={hooks}")

# Cyclic (Paley T_5 doesn't exist since 5 ≡ 1 mod 4)
# Use the tournament with max H
max_h = 0
max_bits = 0
for bits in range(1 << 10):
    T = tournament_from_bits(5, bits)
    h = hamiltonian_path_count(T)
    if h > max_h:
        max_h = h
        max_bits = bits

T = tournament_from_bits(5, max_bits)
hooks = hook_coefficients(T)
print(f"Max H: H={max_h}, hooks={hooks}")

# Look at hook symmetry: h_i vs h_{n+1-i}
print(f"\nHook symmetry check (h_i vs h_(n+1-i)):")
for bits in [0, max_bits, 42]:
    T = tournament_from_bits(5, bits)
    h = hamiltonian_path_count(T)
    hooks = hook_coefficients(T)
    rev = list(reversed(hooks))
    sym = all(hooks[i] == rev[i] for i in range(5))
    print(f"  H={h}: hooks={hooks}, reversed={rev}, symmetric={sym}")

# ── Key question: what determines the hook coefficients? ──
print(f"\n{'='*60}")
print("HOOK COEFFICIENT PATTERNS AT n=4")
print("=" * 60)

from collections import defaultdict
hook_by_h = defaultdict(list)
for bits in range(1 << 6):
    T = tournament_from_bits(4, bits)
    h = hamiltonian_path_count(T)
    hooks = hook_coefficients(T)
    hook_by_h[h].append(hooks)

for h_val in sorted(hook_by_h.keys()):
    patterns = set(tuple(h) for h in hook_by_h[h_val])
    print(f"H={h_val}: {len(hook_by_h[h_val])} tournaments, {len(patterns)} distinct patterns")
    for p in sorted(patterns):
        count = sum(1 for h in hook_by_h[h_val] if tuple(h) == p)
        print(f"  {list(p)} x {count}")

print("\nDone.")
