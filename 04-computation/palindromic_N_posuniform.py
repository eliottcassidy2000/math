#!/usr/bin/env python3
"""
KEY TEST: Does position-uniform imply palindromic N(a,b,j)?

N(a,b,j) = #{Ham paths with {a,b} at consecutive positions {j,j+1}}
(symmetrized: counts both a->b at (j,j+1) AND b->a at (j,j+1))

If palindromic N implies alternating sum = 0 at odd n (THM-051 corollary),
then position-uniform => palindromic N => scalar M would close the proof.

Also test: does self-complementary (T isomorphic to T^op) relate?

kind-pasteur-2026-03-06-S25d
"""

from itertools import permutations
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

def ham_paths_list(T, n):
    paths = []
    for perm in permutations(range(n)):
        ok = True
        for k in range(n-1):
            if T.get((perm[k], perm[k+1]), 0) == 0:
                ok = False; break
        if ok:
            paths.append(perm)
    return paths

def position_matrix(paths, n):
    P = np.zeros((n, n), dtype=int)
    for perm in paths:
        for k in range(n):
            P[perm[k], k] += 1
    return P

def compute_symmetrized_N(paths, n):
    """N[a][b][j] = #{paths with {a,b} at positions {j,j+1}}"""
    N = [[[0]*(n-1) for _ in range(n)] for _ in range(n)]
    for perm in paths:
        for j in range(n-1):
            a, b = perm[j], perm[j+1]
            N[a][b][j] += 1
            N[b][a][j] += 1  # symmetrize
    return N

def is_palindromic_N(N, n):
    for a in range(n):
        for b in range(a+1, n):
            for j in range(n-1):
                if N[a][b][j] != N[a][b][n-2-j]:
                    return False
    return True

def compute_M_from_N(N, n):
    """M[a,b] = sum_j (-1)^j N(a,b,j) for a!=b"""
    M = np.zeros((n, n), dtype=int)
    for a in range(n):
        for b in range(n):
            if a == b:
                continue
            M[a][b] = sum((-1)**j * N[a][b][j] for j in range(n-1))
    return M

# ============================================================
# n=5: Exhaustive test
# ============================================================
print("=" * 70)
print("n=5: Position-uniform => palindromic N(a,b,j)?")
print("=" * 70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
uniform_count = 0
uniform_palindromic = 0
uniform_not_palindromic = 0

for bits in range(1 << len(pairs)):
    T = tournament_from_bits(n, bits)
    paths = ham_paths_list(T, n)
    H = len(paths)
    if H % n != 0:
        continue
    P = position_matrix(paths, n)
    if not all(P[v,k] == H // n for v in range(n) for k in range(n)):
        continue

    uniform_count += 1
    N = compute_symmetrized_N(paths, n)
    if is_palindromic_N(N, n):
        uniform_palindromic += 1
    else:
        uniform_not_palindromic += 1
        if uniform_not_palindromic <= 3:
            print(f"\n  COUNTEREXAMPLE bits={bits}, H={H}")
            for a in range(n):
                for b in range(a+1, n):
                    nab = [N[a][b][j] for j in range(n-1)]
                    if any(nab[j] != nab[n-2-j] for j in range(n-1)):
                        alt = sum((-1)**j * nab[j] for j in range(n-1))
                        print(f"    N({a},{b}) = {nab}, alt_sum = {alt}")

print(f"\n  Position-uniform: {uniform_count}")
print(f"  Palindromic N: {uniform_palindromic}")
print(f"  Non-palindromic N: {uniform_not_palindromic}")

if uniform_not_palindromic == 0:
    print("\n  ==> Position-uniform => Palindromic N at n=5!")
else:
    print("\n  ==> Position-uniform does NOT imply palindromic N at n=5")
    print("  BUT: check if alternating sum still vanishes...")

# ============================================================
# n=7: Circulant tournaments
# ============================================================
print("\n" + "=" * 70)
print("n=7: Circulant tournaments => palindromic N(a,b,j)?")
print("=" * 70)

def circulant_tournament(n, gen_set):
    T = {}
    for i in range(n):
        for j in range(n):
            if i == j: continue
            T[(i,j)] = 1 if (j - i) % n in gen_set else 0
    return T

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

pal_count = 0
non_pal_count = 0

for gs in sorted(gen_sets):
    T = circulant_tournament(n, gs)
    paths = ham_paths_list(T, n)
    H = len(paths)
    N = compute_symmetrized_N(paths, n)

    is_pal = is_palindromic_N(N, n)
    if is_pal:
        pal_count += 1
    else:
        non_pal_count += 1
        print(f"\n  NON-PALINDROMIC: gs={sorted(gs)}, H={H}")
        # Show one example
        for a in range(min(n, 3)):
            for b in range(a+1, min(n, 3)):
                nab = [N[a][b][j] for j in range(n-1)]
                alt = sum((-1)**j * nab[j] for j in range(n-1))
                print(f"    N({a},{b}) = {nab}, alt_sum = {alt}")

print(f"\n  Circulant count: {len(gen_sets)}")
print(f"  Palindromic N: {pal_count}")
print(f"  Non-palindromic N: {non_pal_count}")

if non_pal_count == 0:
    print("\n  ==> ALL circulant n=7 have palindromic N!")
    print("  Combined with position-uniform at n=5:")
    print("  CONJECTURE: Position-uniform => palindromic N at all odd n")
    print("  => M[a,b]=0 for a!=b => M=(H/n)*I")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
