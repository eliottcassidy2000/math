#!/usr/bin/env python3
"""
beta2_doubling_closure_general.py - Generalize doubling-closure to |S| > 2

HYP-217: For |S|=2 circulant digraphs C_n^S with S cap (-S) = empty,
beta_2 = 0 iff S is "doubling-closed" (2s1 = s2 or 2s2 = s1 mod n).

Question: what characterizes beta_2 = 0 for |S| = 3, 4, ...?
For tournaments (|S| = (n-1)/2), we always have beta_2 = 0 (HYP-207).
What about intermediate |S|?

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
from itertools import combinations
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def build_circulant_digraph(n, S):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for s in S:
            A[i][(i + s) % n] = 1
    return A

def compute_betti_2(A, n):
    """Compute just beta_2."""
    allowed = {}
    for p in range(5):
        allowed[p] = enumerate_allowed_paths(A, n, p)
        if not allowed[p]:
            break

    if 2 not in allowed or not allowed[2]:
        return 0

    omega_basis = {}
    for p in range(4):
        if p not in allowed or not allowed[p]:
            omega_basis[p] = np.zeros((0, 0))
            continue
        basis = compute_omega_basis(A, n, p, allowed[p],
                                     allowed[p-1] if p-1 in allowed and p >= 1 else [])
        omega_basis[p] = basis

    bd_omega = {}
    for p in [2, 3]:
        dim_p = omega_basis[p].shape[1] if omega_basis[p].ndim == 2 else 0
        if dim_p == 0:
            continue
        bd = build_full_boundary_matrix(allowed[p], allowed[p-1] if p-1 in allowed else [])
        bd_omega[p] = bd @ omega_basis[p]

    dim_2 = omega_basis[2].shape[1] if omega_basis[2].ndim == 2 else 0
    if dim_2 == 0:
        return 0

    if 2 not in bd_omega:
        ker = dim_2
    else:
        S_v = np.linalg.svd(bd_omega[2], compute_uv=False)
        ker = dim_2 - sum(s > 1e-8 for s in S_v)

    if 3 not in bd_omega:
        im = 0
    else:
        S_v = np.linalg.svd(bd_omega[3], compute_uv=False)
        im = sum(s > 1e-8 for s in S_v)

    return ker - im

def neg_S(S, n):
    return frozenset((n - s) % n for s in S)

def is_doubling_closed(S, n):
    """Check if S is closed under doubling mod n (2*s in S for all s)."""
    for s in S:
        if (2*s) % n not in S:
            return False
    return True

def has_doubling_pair(S, n):
    """Check if S contains any pair (s, 2s mod n)."""
    for s in S:
        if (2*s) % n in S:
            return True
    return False

def doubling_orbit(S, n):
    """Compute the doubling orbit: {2^k * s mod n : s in S, k >= 0}."""
    orbit = set(S)
    changed = True
    while changed:
        changed = False
        new = set()
        for s in orbit:
            d = (2*s) % n
            if d != 0 and d not in orbit:
                new.add(d)
                changed = True
        orbit |= new
    return frozenset(orbit)

print("=" * 70)
print("BETA_2 CHARACTERIZATION FOR GENERAL |S| CIRCULANT DIGRAPHS")
print("=" * 70)

for n in [5, 7, 9, 11]:
    print(f"\n{'='*70}")
    print(f"n = {n}")
    print(f"{'='*70}")

    for s_size in range(2, (n-1)//2 + 1):
        print(f"\n  --- |S| = {s_size} ---")
        tested = set()
        b2_zero = []
        b2_pos = []

        for S_tuple in combinations(range(1, n), s_size):
            S = frozenset(S_tuple)
            nS = neg_S(S, n)
            if S & nS:
                continue
            canon = min(tuple(sorted(S)), tuple(sorted(nS)))
            if canon in tested:
                continue
            tested.add(canon)

            A = build_circulant_digraph(n, S)
            b2 = compute_betti_2(A, n)

            dc = is_doubling_closed(S, n)
            dp = has_doubling_pair(S, n)
            dorb = doubling_orbit(S, n)

            info = {
                'S': tuple(sorted(S)),
                'b2': b2,
                'doubling_closed': dc,
                'has_doubling_pair': dp,
                'doubling_orbit_size': len(dorb),
            }

            if b2 == 0:
                b2_zero.append(info)
            else:
                b2_pos.append(info)

        print(f"  beta_2 = 0: {len(b2_zero)} cases")
        print(f"  beta_2 > 0: {len(b2_pos)} cases")

        # Check doubling-closed correlation
        zero_dc = sum(1 for x in b2_zero if x['doubling_closed'])
        pos_dc = sum(1 for x in b2_pos if x['doubling_closed'])
        zero_dp = sum(1 for x in b2_zero if x['has_doubling_pair'])
        pos_dp = sum(1 for x in b2_pos if x['has_doubling_pair'])

        print(f"  Doubling-closed in b2=0: {zero_dc}/{len(b2_zero)}")
        print(f"  Doubling-closed in b2>0: {pos_dc}/{len(b2_pos)}")
        print(f"  Has-doubling-pair in b2=0: {zero_dp}/{len(b2_zero)}")
        print(f"  Has-doubling-pair in b2>0: {pos_dp}/{len(b2_pos)}")

        # Show details for b2>0 cases
        if b2_pos and len(b2_pos) <= 10:
            for x in b2_pos:
                print(f"    S={x['S']}: b2={x['b2']}, dc={x['doubling_closed']}, dp={x['has_doubling_pair']}, dorb={x['doubling_orbit_size']}")

        # Show details for b2=0 cases
        if b2_zero and len(b2_zero) <= 10:
            for x in b2_zero:
                print(f"    S={x['S']}: b2=0, dc={x['doubling_closed']}, dp={x['has_doubling_pair']}, dorb={x['doubling_orbit_size']}")

    # Tournament case
    if n % 2 == 1:
        print(f"\n  --- |S| = {(n-1)//2} (TOURNAMENT) ---")
        tested = set()
        count = 0
        b2_zero_t = 0
        for S_tuple in combinations(range(1, n), (n-1)//2):
            S = frozenset(S_tuple)
            nS = neg_S(S, n)
            if S & nS:
                continue
            if S | nS != frozenset(range(1, n)):
                continue
            canon = min(tuple(sorted(S)), tuple(sorted(nS)))
            if canon in tested:
                continue
            tested.add(canon)
            count += 1

            A = build_circulant_digraph(n, S)
            b2 = compute_betti_2(A, n)
            if b2 == 0:
                b2_zero_t += 1
            else:
                print(f"    COUNTEREXAMPLE: S={tuple(sorted(S))}, b2={b2}")

            if count >= 50:  # limit for large n
                break

        print(f"  Tournaments tested: {count}, all b2=0: {b2_zero_t == count}")

print("\nDone.")
