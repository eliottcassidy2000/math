#!/usr/bin/env python3
"""
Test: Are ALL circulant tournaments at n=7 self-complementary?
If not, what other mechanism gives palindromic N?

For a circulant T on Z/nZ with generator S:
- T^op has generator S^op = {n-d : d in S}
- T is self-complementary iff there exists r with r*S = S^op (mod n)

kind-pasteur-2026-03-06-S25d
"""

def circulant_tournament(n, gen_set):
    T = {}
    for i in range(n):
        for j in range(n):
            if i == j: continue
            T[(i,j)] = 1 if (j - i) % n in gen_set else 0
    return T

def is_self_comp_circulant(n, S):
    """Check if circulant with generator S is self-complementary."""
    S_op = frozenset((n - d) % n for d in S)
    # Check if there exists r coprime to n with r*S = S_op (mod n)
    for r in range(1, n):
        if n > 1:
            from math import gcd
            if gcd(r, n) != 1:
                continue
        rS = frozenset((r * d) % n for d in S)
        if rS == S_op:
            return True, r
    return False, None

# ============================================================
# n=7: Check all circulant tournaments
# ============================================================
print("=" * 70)
print("n=7: Self-complementary check for ALL circulant tournaments")
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

sc_count = 0
not_sc_count = 0

for gs in sorted(gen_sets):
    is_sc, r = is_self_comp_circulant(n, gs)
    gs_op = frozenset((n - d) % n for d in gs)
    if is_sc:
        sc_count += 1
        print(f"  gs={sorted(gs)}, gs^op={sorted(gs_op)}: SC via r={r}")
    else:
        not_sc_count += 1
        print(f"  gs={sorted(gs)}, gs^op={sorted(gs_op)}: NOT self-comp")

print(f"\n  Self-comp: {sc_count}, Not SC: {not_sc_count}")

# ============================================================
# n=9: Check all circulant tournaments
# ============================================================
print("\n" + "=" * 70)
print("n=9: Self-complementary check for ALL circulant tournaments")
print("=" * 70)

n = 9
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

sc_count = 0
not_sc_count = 0

for gs in sorted(gen_sets):
    is_sc, r = is_self_comp_circulant(n, gs)
    gs_op = frozenset((n - d) % n for d in gs)
    if is_sc:
        sc_count += 1
    else:
        not_sc_count += 1
        print(f"  gs={sorted(gs)}, gs^op={sorted(gs_op)}: NOT self-comp")

print(f"\n  Total: {len(gen_sets)}")
print(f"  Self-comp: {sc_count}, Not SC: {not_sc_count}")

if not_sc_count > 0:
    print("\n  Some circulant tournaments are NOT self-complementary!")
    print("  But they still have palindromic N. Why?")
    print("  Possible: vertex-transitivity alone suffices (not self-comp).")

# ============================================================
# n=5: Is position-uniform equivalent to vertex-transitive?
# ============================================================
print("\n" + "=" * 70)
print("n=5: Position-uniform vs vertex-transitive")
print("=" * 70)

from itertools import permutations
import numpy as np

def tournament_from_bits(n, bits):
    pairs_list = [(i,j) for i in range(n) for j in range(i+1, n)]
    T = {}
    for idx, (i,j) in enumerate(pairs_list):
        if (bits >> idx) & 1:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1
    return T

def automorphism_group(T, n):
    auts = []
    for perm in permutations(range(n)):
        is_aut = True
        for i in range(n):
            for j in range(n):
                if i == j: continue
                if T[(perm[i], perm[j])] != T[(i,j)]:
                    is_aut = False
                    break
            if not is_aut:
                break
        if is_aut:
            auts.append(perm)
    return auts

def is_vertex_transitive(auts, n):
    for v in range(1, n):
        found = False
        for perm in auts:
            if perm[0] == v:
                found = True
                break
        if not found:
            return False
    return True

n = 5
pairs_list = [(i,j) for i in range(n) for j in range(i+1, n)]
uniform_vt = 0
uniform_not_vt = 0
uniform_count = 0

for bits in range(1 << len(pairs_list)):
    T = tournament_from_bits(n, bits)
    paths = []
    for perm in permutations(range(n)):
        ok = True
        for k in range(n-1):
            if T.get((perm[k], perm[k+1]), 0) == 0:
                ok = False; break
        if ok:
            paths.append(perm)
    H = len(paths)
    if H % n != 0:
        continue
    P = np.zeros((n, n), dtype=int)
    for perm in paths:
        for k in range(n):
            P[perm[k], k] += 1
    if not all(P[v,k] == H // n for v in range(n) for k in range(n)):
        continue

    uniform_count += 1
    auts = automorphism_group(T, n)
    vt = is_vertex_transitive(auts, n)
    if vt:
        uniform_vt += 1
    else:
        uniform_not_vt += 1
        scores = sorted([sum(T[(i,j)] for j in range(n) if j != i) for i in range(n)], reverse=True)
        if uniform_not_vt <= 3:
            print(f"  Not VT but uniform: bits={bits}, H={H}, scores={scores}, |Aut|={len(auts)}")

print(f"\n  Position-uniform: {uniform_count}")
print(f"  Vertex-transitive: {uniform_vt}")
print(f"  Not VT: {uniform_not_vt}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
