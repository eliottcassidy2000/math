#!/usr/bin/env python3
"""
Fast Interlacing Test — skip canonical form, just test interlacing.
kind-pasteur-2026-03-06-S18
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import (tournament_from_bits, random_tournament,
                             find_odd_cycles, conflict_graph)
import numpy as np
import random

def indep_poly_coeffs(adj):
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
    if len(coeffs) <= 1:
        return []
    poly = list(reversed(coeffs))
    roots = np.roots(poly)
    return sorted([r.real for r in roots])

def check_interlacing(roots_big, roots_small):
    """Check roots_small interlaces roots_big (d vs d-1)."""
    d = len(roots_big)
    ds = len(roots_small)
    if ds == 0:
        return True
    if ds == d - 1:
        for i in range(ds):
            if roots_big[i] > roots_small[i] + 1e-5:
                return False
            if roots_small[i] > roots_big[i+1] + 1e-5:
                return False
        return True
    elif ds == d:
        # Same degree — check mutual interlacing
        ok1 = all(roots_big[i] <= roots_small[i] + 1e-5 for i in range(d)) and \
              all(roots_small[i] <= roots_big[i+1] + 1e-5 for i in range(d-1))
        ok2 = all(roots_small[i] <= roots_big[i] + 1e-5 for i in range(d)) and \
              all(roots_big[i] <= roots_small[i+1] + 1e-5 for i in range(d-1))
        return ok1 or ok2
    elif ds < d - 1:
        return True  # trivially
    return False

def delete_vertex(T, v):
    n = len(T)
    return [[T[i][j] for j in range(n) if j != v] for i in range(n) if i != v]

def get_ip(T):
    """Get independence polynomial of Omega(T)."""
    cycles = find_odd_cycles(T)
    if not cycles:
        return [1]
    cg = conflict_graph(cycles)
    return indep_poly_coeffs(cg)

# n=5 exhaustive
print("=" * 70)
print("FAST INTERLACING TEST")
print("=" * 70)

for n in [5, 6]:
    m = n * (n - 1) // 2
    failures = 0
    total = 0
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        coeffs_T = get_ip(T)
        roots_T = get_roots(coeffs_T)

        for v in range(n):
            Tv = delete_vertex(T, v)
            coeffs_Tv = get_ip(Tv)
            roots_Tv = get_roots(coeffs_Tv)

            if roots_T and roots_Tv:
                ok = check_interlacing(roots_T, roots_Tv)
                if not ok:
                    failures += 1
                    if failures <= 3:
                        print(f"  FAIL bits={bits} v={v}: T={coeffs_T} Tv={coeffs_Tv}")
                        print(f"    roots_T={[f'{r:.4f}' for r in roots_T]}")
                        print(f"    roots_Tv={[f'{r:.4f}' for r in roots_Tv]}")
            total += 1

        if bits % 1000 == 0 and bits > 0:
            print(f"  n={n}: {bits}/{1<<m} done, {failures} failures so far")

    print(f"\nn={n}: {failures}/{total} interlacing failures")

# n=7,8 random
random.seed(42)
for n in [7, 8]:
    failures = 0
    total = 0
    num_trials = 200 if n <= 7 else 100

    for trial in range(num_trials):
        T = random_tournament(n)
        coeffs_T = get_ip(T)
        roots_T = get_roots(coeffs_T)

        for v in range(n):
            Tv = delete_vertex(T, v)
            coeffs_Tv = get_ip(Tv)
            roots_Tv = get_roots(coeffs_Tv)

            if roots_T and roots_Tv:
                ok = check_interlacing(roots_T, roots_Tv)
                if not ok:
                    failures += 1
                    if failures <= 3:
                        print(f"  FAIL trial={trial} v={v}: T={coeffs_T} Tv={coeffs_Tv}")
            total += 1

    print(f"\nn={n}: {failures}/{total} interlacing failures ({num_trials} random)")

# Also test: does I(Omega_3(T-v), x) interlace I(Omega_3(T), x)?
# (3-cycle subgraph only)
print(f"\n{'='*70}")
print("3-CYCLE SUBGRAPH INTERLACING")
print("=" * 70)

def get_ip_3only(T):
    """Independence polynomial using only 3-cycles."""
    cycles = [c for c in find_odd_cycles(T) if len(c) == 3]
    if not cycles:
        return [1]
    cg = conflict_graph(cycles)
    return indep_poly_coeffs(cg)

for n in [7, 8, 9]:
    failures = 0
    total = 0
    num_trials = 100

    for trial in range(num_trials):
        T = random_tournament(n)
        coeffs_T = get_ip_3only(T)
        roots_T = get_roots(coeffs_T)

        for v in range(n):
            Tv = delete_vertex(T, v)
            coeffs_Tv = get_ip_3only(Tv)
            roots_Tv = get_roots(coeffs_Tv)

            if roots_T and roots_Tv:
                ok = check_interlacing(roots_T, roots_Tv)
                if not ok:
                    failures += 1
                    if failures <= 3:
                        print(f"  FAIL trial={trial} v={v}: T deg {len(coeffs_T)-1} Tv deg {len(coeffs_Tv)-1}")
                        print(f"    roots_T={[f'{r:.4f}' for r in roots_T]}")
                        print(f"    roots_Tv={[f'{r:.4f}' for r in roots_Tv]}")
            total += 1

    print(f"n={n} (3-cycle): {failures}/{total} interlacing failures ({num_trials} random)")

print("\nDone.")
