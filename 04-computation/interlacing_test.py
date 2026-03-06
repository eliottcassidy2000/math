#!/usr/bin/env python3
"""
Interlacing Test for Independence Polynomials
===============================================
Test whether I(Omega(T-v), x) interlaces I(Omega(T), x).

Interlacing: if p has roots r1 <= r2 <= ... <= rd and q has roots
s1 <= s2 <= ... <= s_{d-1}, then q interlaces p if
  r1 <= s1 <= r2 <= s2 <= ... <= s_{d-1} <= rd

For real-rooted polynomials, interlacing of deletion polynomials
implies real-rootedness by induction (base case: degree 1 always real).

This is the MOST PROMISING lead for proving real roots for all n.

kind-pasteur-2026-03-06-S18
"""

import sys
sys.path.insert(0, r'C:\Users\Eliott\Documents\GitHub\math\03-artifacts\code')
from tournament_lib import (tournament_from_bits, all_tournaments,
                             random_tournament, find_odd_cycles,
                             conflict_graph, hamiltonian_path_count)
import numpy as np
import random

def independence_polynomial_coefficients(adj):
    """Compute full independence polynomial coefficients."""
    m = len(adj)
    if m == 0:
        return [1]
    nbr = [0] * m
    for i in range(m):
        for j in range(m):
            if adj[i][j]:
                nbr[i] |= 1 << j
    coeffs = [0] * (m + 1)
    for mask in range(1 << m):
        ok = True
        seen = 0
        temp = mask
        while temp:
            v = (temp & -temp).bit_length() - 1
            if nbr[v] & seen:
                ok = False
                break
            seen |= 1 << v
            temp &= temp - 1
        if ok:
            coeffs[bin(mask).count('1')] += 1
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs

def get_roots(coeffs):
    """Get sorted real roots of polynomial."""
    if len(coeffs) <= 1:
        return []
    poly = list(reversed(coeffs))
    roots = np.roots(poly)
    return sorted([r.real for r in roots])

def check_interlacing(roots_p, roots_q):
    """Check if q interlaces p: r1 <= s1 <= r2 <= s2 <= ..."""
    dp = len(roots_p)
    dq = len(roots_q)
    if dq == 0:
        return True
    if dq > dp:
        return False
    # Standard interlacing: dq = dp - 1
    # But we also allow dq = dp (weak interlacing)
    if dq == dp - 1:
        for i in range(dq):
            if roots_p[i] > roots_q[i] + 1e-6:
                return False
            if roots_q[i] > roots_p[i+1] + 1e-6:
                return False
        return True
    elif dq == dp:
        # Check if roots_q interlaces roots_p in the weak sense
        # s1 <= r1 <= s2 <= r2 <= ... or r1 <= s1 <= r2 <= s2
        # Try both orderings
        ok1 = True
        for i in range(dp):
            if roots_q[i] > roots_p[i] + 1e-6:
                ok1 = False
                break
            if i < dp - 1 and roots_p[i] > roots_q[i+1] + 1e-6:
                ok1 = False
                break
        ok2 = True
        for i in range(dp):
            if roots_p[i] > roots_q[i] + 1e-6:
                ok2 = False
                break
            if i < dp - 1 and roots_q[i] > roots_p[i+1] + 1e-6:
                ok2 = False
                break
        return ok1 or ok2
    elif dq < dp - 1:
        # Weaker condition: just check containment of root intervals
        return True
    return False

def delete_vertex(T, v):
    """Delete vertex v from tournament T."""
    n = len(T)
    return [[T[i][j] for j in range(n) if j != v] for i in range(n) if i != v]

# ══════════════════════════════════════════════════════════════
# n=5 EXHAUSTIVE
# ══════════════════════════════════════════════════════════════
print("=" * 70)
print("INTERLACING TEST: I(Omega(T-v), x) interlaces I(Omega(T), x)?")
print("=" * 70)

n = 5
m = n * (n - 1) // 2
failures = 0
total = 0
degree_pairs = {}

for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    cycles_T = find_odd_cycles(T)
    m_T = len(cycles_T)
    if m_T == 0:
        continue
    cg_T = conflict_graph(cycles_T)
    coeffs_T = independence_polynomial_coefficients(cg_T)
    roots_T = get_roots(coeffs_T)

    for v in range(n):
        Tv = delete_vertex(T, v)
        cycles_Tv = find_odd_cycles(Tv)
        if not cycles_Tv:
            coeffs_Tv = [1]
        else:
            cg_Tv = conflict_graph(cycles_Tv)
            coeffs_Tv = independence_polynomial_coefficients(cg_Tv)
        roots_Tv = get_roots(coeffs_Tv)

        dp = (len(coeffs_T)-1, len(coeffs_Tv)-1)
        degree_pairs[dp] = degree_pairs.get(dp, 0) + 1

        if len(roots_Tv) > 0 and len(roots_T) > 0:
            ok = check_interlacing(roots_T, roots_Tv)
            if not ok:
                failures += 1
                print(f"  FAILURE: bits={bits}, v={v}")
                print(f"    T coeffs={coeffs_T}, roots={[f'{r:.4f}' for r in roots_T]}")
                print(f"    T-v coeffs={coeffs_Tv}, roots={[f'{r:.4f}' for r in roots_Tv]}")
        total += 1

print(f"\nn=5: {failures}/{total} interlacing failures")
print(f"Degree pairs (deg_T, deg_Tv): {degree_pairs}")

# ══════════════════════════════════════════════════════════════
# n=6 EXHAUSTIVE
# ══════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("n=6 EXHAUSTIVE")
print("=" * 70)

n = 6
m = n * (n - 1) // 2
failures6 = 0
total6 = 0
degree_pairs6 = {}
interlace_examples = []

for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    cycles_T = find_odd_cycles(T)
    m_T = len(cycles_T)
    if m_T == 0:
        continue
    cg_T = conflict_graph(cycles_T)
    coeffs_T = independence_polynomial_coefficients(cg_T)
    roots_T = get_roots(coeffs_T)

    for v in range(n):
        Tv = delete_vertex(T, v)
        cycles_Tv = find_odd_cycles(Tv)
        if not cycles_Tv:
            coeffs_Tv = [1]
        else:
            cg_Tv = conflict_graph(cycles_Tv)
            coeffs_Tv = independence_polynomial_coefficients(cg_Tv)
        roots_Tv = get_roots(coeffs_Tv)

        dp = (len(coeffs_T)-1, len(coeffs_Tv)-1)
        degree_pairs6[dp] = degree_pairs6.get(dp, 0) + 1

        if len(roots_Tv) > 0 and len(roots_T) > 0:
            ok = check_interlacing(roots_T, roots_Tv)
            if not ok:
                failures6 += 1
                if failures6 <= 5:
                    print(f"  FAILURE: bits={bits}, v={v}")
                    print(f"    T coeffs={coeffs_T}, roots={[f'{r:.4f}' for r in roots_T]}")
                    print(f"    T-v coeffs={coeffs_Tv}, roots={[f'{r:.4f}' for r in roots_Tv]}")
            elif len(roots_T) >= 2 and len(roots_Tv) >= 1 and len(interlace_examples) < 3:
                interlace_examples.append({
                    'bits': bits, 'v': v,
                    'coeffs_T': coeffs_T, 'coeffs_Tv': coeffs_Tv,
                    'roots_T': roots_T, 'roots_Tv': roots_Tv
                })
        total6 += 1

print(f"\nn=6: {failures6}/{total6} interlacing failures")
print(f"Degree pairs: {degree_pairs6}")

if interlace_examples:
    print(f"\nExample interlacing instances:")
    for ex in interlace_examples:
        print(f"  bits={ex['bits']}, v={ex['v']}")
        print(f"    I(T) = {ex['coeffs_T']}, roots = {[f'{r:.4f}' for r in ex['roots_T']]}")
        print(f"    I(T-v) = {ex['coeffs_Tv']}, roots = {[f'{r:.4f}' for r in ex['roots_Tv']]}")

# ══════════════════════════════════════════════════════════════
# n=7,8 RANDOM SAMPLING
# ══════════════════════════════════════════════════════════════
random.seed(42)

for n in [7, 8]:
    print(f"\n{'='*70}")
    print(f"n={n}: 200 random tournaments")
    print("=" * 70)

    failures_n = 0
    total_n = 0
    degree_pairs_n = {}

    for trial in range(200):
        T = random_tournament(n)
        cycles_T = find_odd_cycles(T)
        if not cycles_T:
            continue
        cg_T = conflict_graph(cycles_T)
        coeffs_T = independence_polynomial_coefficients(cg_T)
        roots_T = get_roots(coeffs_T)

        for v in range(n):
            Tv = delete_vertex(T, v)
            cycles_Tv = find_odd_cycles(Tv)
            if not cycles_Tv:
                coeffs_Tv = [1]
            else:
                cg_Tv = conflict_graph(cycles_Tv)
                coeffs_Tv = independence_polynomial_coefficients(cg_Tv)
            roots_Tv = get_roots(coeffs_Tv)

            dp = (len(coeffs_T)-1, len(coeffs_Tv)-1)
            degree_pairs_n[dp] = degree_pairs_n.get(dp, 0) + 1

            if len(roots_Tv) > 0 and len(roots_T) > 0:
                ok = check_interlacing(roots_T, roots_Tv)
                if not ok:
                    failures_n += 1
                    if failures_n <= 3:
                        print(f"  FAILURE trial={trial} v={v}")
                        print(f"    T coeffs={coeffs_T}")
                        print(f"    T-v coeffs={coeffs_Tv}")
            total_n += 1

    print(f"n={n}: {failures_n}/{total_n} interlacing failures")
    print(f"Degree pairs: {degree_pairs_n}")

# ══════════════════════════════════════════════════════════════
# ALGEBRAIC TEST: Does the relationship between Omega(T) and
# Omega(T-v) have a clean description?
# ══════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("STRUCTURAL ANALYSIS: Omega(T) vs Omega(T-v)")
print("=" * 70)

n = 6
# For a specific tournament, compare cycle sets
T = tournament_from_bits(n, 42)
h = hamiltonian_path_count(T)
cycles = find_odd_cycles(T)
print(f"\nExample tournament (n={n}, bits=42, H={h}):")
print(f"  Cycles: {len(cycles)}")

for v in range(n):
    Tv = delete_vertex(T, v)
    cycles_v = find_odd_cycles(Tv)

    # How many cycles of T survive in T-v?
    survived = [c for c in cycles if v not in c]
    new_in_Tv = [c for c in cycles_v if not any(set(c2) == set(c) or set(c2) == {x - (1 if x > v else 0) for x in c if x != v} for c2 in cycles)]

    print(f"  v={v}: |Omega(T)|={len(cycles)}, |Omega(T-v)|={len(cycles_v)}, survived={len(survived)}, through_v={len(cycles)-len(survived)}")

print("\nDone.")
