# ⚠️ WARNING: This script uses QR mod p for p ≡ 1 (mod 4), which does NOT
# produce a tournament (S ∩ (-S) ≠ ∅ gives bidirectional edges).
# Results for those primes are INVALID. See MISTAKE-011b.
# Valid Paley tournaments require p ≡ 3 (mod 4).

#!/usr/bin/env python3
"""
TEST: Does the complement pairing E_a(S+a)*B_b(R+b) = E_a(R+a)*B_b(S+b)
hold for all position-uniform tournaments?

At odd n, |U| = n-2 is odd, so S and U\S have opposite signs.
If the cross-ratio identity holds, then M[a,b] = 0 for T[a,b]=0.

This would prove NONHAM=0 for all position-uniform tournaments at odd n.

kind-pasteur-2026-03-06-S25c
"""

from itertools import permutations
import time

def circulant_tournament(n, gen_set):
    T = {}
    for i in range(n):
        for j in range(n):
            if i == j: continue
            T[(i,j)] = 1 if (j - i) % n in gen_set else 0
    return T

def E_paths(T, verts, a):
    verts = list(verts)
    if len(verts) == 1:
        return 1 if verts[0] == a else 0
    count = 0
    for p in permutations(verts):
        if p[-1] != a: continue
        valid = True
        for k in range(len(p)-1):
            if T.get((p[k], p[k+1]), 0) != 1:
                valid = False; break
        if valid: count += 1
    return count

def B_paths(T, verts, b):
    verts = list(verts)
    if len(verts) == 1:
        return 1 if verts[0] == b else 0
    count = 0
    for p in permutations(verts):
        if p[0] != b: continue
        valid = True
        for k in range(len(p)-1):
            if T.get((p[k], p[k+1]), 0) != 1:
                valid = False; break
        if valid: count += 1
    return count


# ============================================================
# n=5 Paley: Cross-ratio identity
# ============================================================
print("=" * 70)
print("n=5 Paley: Cross-ratio identity E_a(S+a)*B_b(R+b) = E_a(R+a)*B_b(S+b)")
print("=" * 70)

n = 5
T = circulant_tournament(n, {1, 4})
fails = 0

for a in range(n):
    for b in range(n):
        if a == b: continue
        if T[(a,b)] == 1: continue  # Only check non-edges

        U = [v for v in range(n) if v != a and v != b]

        for mask in range(1 << len(U)):
            S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
            R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
            S_set = sorted(set(S_list) | {a})
            R_set = sorted(set(R) | {b})

            # Also compute complement pairing
            comp_S_set = sorted(set(R) | {a})
            comp_R_set = sorted(set(S_list) | {b})

            lhs = E_paths(T, S_set, a) * B_paths(T, R_set, b)
            rhs = E_paths(T, comp_S_set, a) * B_paths(T, comp_R_set, b)

            if lhs != rhs:
                fails += 1
                if fails <= 5:
                    print(f"  FAIL: ({a},{b}), S={S_list}: LHS={lhs}, RHS={rhs}")

print(f"  n=5 Paley: {fails} cross-ratio failures")


# ============================================================
# n=7 circulant tournaments: Cross-ratio identity
# ============================================================
print("\n" + "=" * 70)
print("n=7 circulant: Cross-ratio identity")
print("=" * 70)

n = 7
half = list(range(1, (n+1)//2))
gen_sets = set()
for mask in range(1 << len(half)):
    gs = set()
    for k, d in enumerate(half):
        if mask & (1 << k):
            gs.add(d)
        else:
            gs.add(n - d)
    gen_sets.add(frozenset(gs))

t0 = time.time()
for gs in sorted(gen_sets):
    T = circulant_tournament(n, gs)
    fails = 0
    total_pairs = 0

    # By circulant symmetry, check only a=0
    a = 0
    for b in range(1, n):
        if T[(a,b)] == 1: continue

        U = [v for v in range(n) if v != a and v != b]

        for mask in range(1 << len(U)):
            S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
            R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
            S_set = sorted(set(S_list) | {a})
            R_set = sorted(set(R) | {b})
            comp_S_set = sorted(set(R) | {a})
            comp_R_set = sorted(set(S_list) | {b})

            lhs = E_paths(T, S_set, a) * B_paths(T, R_set, b)
            rhs = E_paths(T, comp_S_set, a) * B_paths(T, comp_R_set, b)
            total_pairs += 1

            if lhs != rhs:
                fails += 1
                if fails <= 2:
                    print(f"  gen_set={sorted(gs)}: ({a},{b}), S={S_list}")
                    print(f"    LHS: E_a({S_set})={E_paths(T, S_set, a)}, B_b({R_set})={B_paths(T, R_set, b)} => {lhs}")
                    print(f"    RHS: E_a({comp_S_set})={E_paths(T, comp_S_set, a)}, B_b({comp_R_set})={B_paths(T, comp_R_set, b)} => {rhs}")

    elapsed = time.time() - t0
    status = "OK" if fails == 0 else f"FAIL ({fails})"
    print(f"  gen_set={sorted(gs)}: {status} ({total_pairs} pairs, {elapsed:.1f}s)")


# ============================================================
# n=5 ALL position-uniform: Cross-ratio identity
# ============================================================
print("\n" + "=" * 70)
print("n=5 ALL position-uniform: Cross-ratio identity")
print("=" * 70)

import numpy as np

def tournament_from_bits(n, bits):
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    T = {}
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1
    return T

def count_H(T, n):
    count = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            if T.get((perm[k], perm[k+1]), 0) == 0:
                prod = 0; break
        count += prod
    return count

def position_matrix(T, n):
    P = np.zeros((n, n), dtype=int)
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            if T.get((perm[k], perm[k+1]), 0) == 0:
                prod = 0; break
        if prod > 0:
            for k in range(n):
                P[perm[k], k] += 1
    return P

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
total_uniform = 0
cross_ratio_fails = 0
cross_ratio_ok = 0

for bits in range(1 << len(pairs)):
    T = tournament_from_bits(n, bits)
    H = count_H(T, n)
    P = position_matrix(T, n)
    is_uniform = H % n == 0 and all(P[v,k] == H // n for v in range(n) for k in range(n))
    if not is_uniform:
        continue

    total_uniform += 1
    this_fails = 0

    for a in range(n):
        for b in range(n):
            if a == b: continue
            if T[(a,b)] == 1: continue

            U = [v for v in range(n) if v != a and v != b]

            for mask in range(1 << len(U)):
                S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
                R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
                S_set = sorted(set(S_list) | {a})
                R_set = sorted(set(R) | {b})
                comp_S_set = sorted(set(R) | {a})
                comp_R_set = sorted(set(S_list) | {b})

                lhs = E_paths(T, S_set, a) * B_paths(T, R_set, b)
                rhs = E_paths(T, comp_S_set, a) * B_paths(T, comp_R_set, b)

                if lhs != rhs:
                    this_fails += 1

    if this_fails > 0:
        cross_ratio_fails += 1
        if cross_ratio_fails <= 3:
            print(f"  bits={bits}: {this_fails} cross-ratio failures")
    else:
        cross_ratio_ok += 1

print(f"\n  Position-uniform tournaments: {total_uniform}")
print(f"  Cross-ratio holds for all: {cross_ratio_ok}")
print(f"  Cross-ratio fails for some: {cross_ratio_fails}")

if cross_ratio_fails == 0:
    print("""
  THEOREM CANDIDATE: For ALL position-uniform tournaments at odd n,
  the complement pairing identity holds:

    E_a(S+a) * B_b(R+b) = E_a(R+a) * B_b(S+b)

  for all non-edge pairs (a,b) with T[a,b]=0, and all subsets S of U.

  Since |U| = n-2 is odd at odd n, S and U\\S have opposite (-1)^|S| signs.
  The complement pairing then gives M[a,b] = 0 for all non-edges.
  Combined with THM-030 (symmetry), this yields M = (H/n)*I.
""")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
