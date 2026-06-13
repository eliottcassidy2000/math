#!/usr/bin/env python3
"""
Fast test: Does every SC tournament have an involution anti-automorphism?

Avoids canonical form computation. Instead:
1. Enumerate all tournaments via bit-encoding
2. For each, check self-converse by finding ANY anti-automorphism
3. If SC, check if at least one anti-aut is an involution

Uses early termination and avoids storing canonical forms.

Author: opus-2026-03-06-S18
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits
from itertools import permutations

def find_involution_anti_aut(T):
    """Find an involution anti-aut of T, or return None."""
    n = len(T)
    for perm in permutations(range(n)):
        # Check if perm is an anti-aut AND an involution
        # Involution: perm[perm[i]] == i for all i
        is_involution = all(perm[perm[i]] == i for i in range(n))
        if not is_involution:
            continue
        # Anti-aut: T[perm[i]][perm[j]] == 1 - T[i][j] for all i≠j
        ok = True
        for i in range(n):
            for j in range(i+1, n):
                if T[perm[i]][perm[j]] != (1 - T[i][j]):
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return perm
    return None

def find_any_anti_aut(T):
    """Find any anti-aut of T, or return None."""
    n = len(T)
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(i+1, n):
                if T[perm[i]][perm[j]] != (1 - T[i][j]):
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return perm
    return None

for n in [4, 5, 6, 7]:
    m = n * (n - 1) // 2
    total = 1 << m
    print(f"\n{'='*60}")
    print(f"n={n}: Testing {total} tournaments")
    print("=" * 60)

    sc_count = 0
    has_inv = 0
    no_inv = 0

    for bits in range(total):
        T = tournament_from_bits(n, bits)

        # Quick check: find any anti-aut
        aa = find_any_anti_aut(T)
        if aa is None:
            continue  # Not SC
        sc_count += 1

        # Check for involution anti-aut
        inv_aa = find_involution_anti_aut(T)
        if inv_aa is not None:
            has_inv += 1
        else:
            no_inv += 1
            print(f"  NO INVOLUTION ANTI-AUT! bits={bits}")

        if sc_count % 100 == 0 and n >= 7:
            print(f"  ... {sc_count} SC tournaments checked so far ({has_inv} with inv, {no_inv} without)")

    print(f"\nn={n}: {sc_count} SC tournaments (NOT iso classes)")
    print(f"  Has involution anti-aut: {has_inv}")
    print(f"  No involution anti-aut: {no_inv}")
    if no_inv == 0:
        print(f"  CONJECTURE HOLDS at n={n}")
    else:
        print(f"  CONJECTURE FAILS at n={n}!")

print("\nDone.")
