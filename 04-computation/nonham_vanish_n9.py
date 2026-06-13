#!/usr/bin/env python3
"""
Test NONHAM=0 for circulant tournaments at n=9.

n=9 has C(4,1)*C(4,2)*... well, let's enumerate properly.
For n=9, half = {1,2,3,4}, complement = {5,6,7,8}.
Each d in {1,2,3,4}: either d or 9-d in gen_set.
So 2^4 = 16 circulant tournaments.

7! = 5040 permutations for sub-path counting.
For n=9 we need 9! = 362880 permutations for H computation, which is slow.
But for NONHAM check on pair (0,b), we need permutations of subsets
of size up to 8, which is up to 8! = 40320. Should be feasible but slow.

We'll use memoization and only check (0,b) pairs by circulant symmetry.

kind-pasteur-2026-03-06-S25c
"""

from itertools import permutations
import numpy as np
import time

def circulant_tournament(n, gen_set):
    T = {}
    for i in range(n):
        for j in range(n):
            if i == j: continue
            T[(i,j)] = 1 if (j - i) % n in gen_set else 0
    return T

def count_H(T, n):
    count = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            if T.get((perm[k], perm[k+1]), 0) == 0:
                prod = 0
                break
        count += prod
    return count

def position_matrix(T, n):
    P = np.zeros((n, n), dtype=int)
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            if T.get((perm[k], perm[k+1]), 0) == 0:
                prod = 0
                break
        if prod > 0:
            for k in range(n):
                P[perm[k], k] += 1
    return P

def E_paths(T, verts, a):
    """Count Ham paths in T[verts] ending at a."""
    verts = list(verts)
    if len(verts) == 1:
        return 1 if verts[0] == a else 0
    count = 0
    for p in permutations(verts):
        if p[-1] != a: continue
        valid = True
        for k in range(len(p)-1):
            if T.get((p[k], p[k+1]), 0) != 1:
                valid = False
                break
        if valid:
            count += 1
    return count

def B_paths(T, verts, b):
    """Count Ham paths in T[verts] starting at b."""
    verts = list(verts)
    if len(verts) == 1:
        return 1 if verts[0] == b else 0
    count = 0
    for p in permutations(verts):
        if p[0] != b: continue
        valid = True
        for k in range(len(p)-1):
            if T.get((p[k], p[k+1]), 0) != 1:
                valid = False
                break
        if valid:
            count += 1
    return count

def check_nonham(T, n, a, b):
    """Check NONHAM(a,b) = M[a,b] when T[a,b]=0."""
    U = [v for v in range(n) if v != a and v != b]
    total = 0

    for mask in range(1 << len(U)):
        S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
        R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
        sign = (-1)**len(S_list)
        S_set = sorted(set(S_list) | {a})
        R_set = sorted(set(R) | {b})

        ea = E_paths(T, S_set, a)
        bb = B_paths(T, R_set, b)

        total += sign * ea * bb

    return total


print("=" * 70)
print("n=9: NONHAM=0 for circulant tournaments?")
print("=" * 70)

n = 9
half = list(range(1, (n+1)//2))  # {1,2,3,4}
gen_sets = []
for mask in range(1 << len(half)):
    gs = set()
    for k, d in enumerate(half):
        if mask & (1 << k):
            gs.add(d)
        else:
            gs.add(n - d)
    gen_sets.append(frozenset(gs))
gen_sets = list(set(gen_sets))

print(f"  {len(gen_sets)} distinct circulant tournaments at n={n}")
t0 = time.time()

fail_count = 0
for idx, gs in enumerate(gen_sets):
    T = circulant_tournament(n, gs)

    # Quick H computation
    H = count_H(T, n)
    elapsed = time.time() - t0
    print(f"  [{idx+1}/{len(gen_sets)}] gen_set={sorted(gs)}, H={H} ({elapsed:.1f}s)", end="", flush=True)

    # Check position uniformity
    P = position_matrix(T, n)
    is_uniform = H % n == 0 and all(P[v,k] == H // n for v in range(n) for k in range(n))
    if not is_uniform:
        print(f" NOT UNIFORM (skipping)")
        continue

    # By circulant symmetry, only check (0, b) for b=1,...,n-1
    # And by M symmetry, M[0,b] = M[b,0], so check both edge directions
    found_fail = False
    for b in range(1, n):
        if T.get((0,b),0) == 1:
            continue  # NONHAM=0 trivially when T[a,b]=1
        nh = check_nonham(T, n, 0, b)
        if nh != 0:
            found_fail = True
            fail_count += 1
            print(f" FAIL: (0,{b}): NONHAM={nh}")
            break

    if not found_fail:
        elapsed = time.time() - t0
        print(f" OK: NONHAM=0 ({elapsed:.1f}s)")

print(f"\n  Checked {len(gen_sets)} circulant tournaments, {fail_count} failures")
if fail_count == 0:
    print(f"  CONFIRMED: NONHAM=0 for ALL circulant n={n} tournaments")

print(f"\n  Total time: {time.time()-t0:.1f}s")
print("=" * 70)
print("DONE")
print("=" * 70)
