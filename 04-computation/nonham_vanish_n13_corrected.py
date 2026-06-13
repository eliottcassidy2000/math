#!/usr/bin/env python3
"""
CORRECTED: NONHAM vanishing test for n=13 circulant tournament.

The previous test (nonham_vanish_n13_paley.py) used QR mod 13 = {1,3,4,9,10,12},
which is NOT a tournament (p=13 = 1 mod 4, MISTAKE-011).

This test uses a valid circulant tournament S = {1,2,3,4,5,6}.

kind-pasteur-2026-03-06-S25d
"""

import time

def circulant_tournament(n, gen_set):
    T = {}
    for i in range(n):
        for j in range(n):
            if i == j: continue
            T[(i,j)] = 1 if (j - i) % n in gen_set else 0
    return T

def count_paths_ending_at(T, verts, target):
    """Count Ham paths on verts ending at target, using DP."""
    verts = list(verts)
    n = len(verts)
    if n == 1:
        return 1 if verts[0] == target else 0
    v_to_idx = {v: i for i, v in enumerate(verts)}
    if target not in v_to_idx:
        return 0
    full = (1 << n) - 1
    dp = {}
    for i, v in enumerate(verts):
        dp[(1 << i, i)] = 1
    for mask in range(1, 1 << n):
        for li in range(n):
            if not ((mask >> li) & 1): continue
            cnt = dp.get((mask, li), 0)
            if cnt == 0: continue
            for ni in range(n):
                if (mask >> ni) & 1: continue
                if T.get((verts[li], verts[ni]), 0) == 0: continue
                nkey = (mask | (1 << ni), ni)
                dp[nkey] = dp.get(nkey, 0) + cnt
    return dp.get((full, v_to_idx[target]), 0)

def count_paths_starting_at(T, verts, source):
    """Count Ham paths on verts starting at source, using DP."""
    verts = list(verts)
    n = len(verts)
    if n == 1:
        return 1 if verts[0] == source else 0
    v_to_idx = {v: i for i, v in enumerate(verts)}
    if source not in v_to_idx:
        return 0
    full = (1 << n) - 1
    dp = {}
    si = v_to_idx[source]
    dp[(1 << si, si)] = 1
    for mask in range(1, 1 << n):
        for li in range(n):
            if not ((mask >> li) & 1): continue
            cnt = dp.get((mask, li), 0)
            if cnt == 0: continue
            for ni in range(n):
                if (mask >> ni) & 1: continue
                if T.get((verts[li], verts[ni]), 0) == 0: continue
                nkey = (mask | (1 << ni), ni)
                dp[nkey] = dp.get(nkey, 0) + cnt
    total = 0
    for li in range(n):
        total += dp.get((full, li), 0)
    return total

n = 13
S = {1, 2, 3, 4, 5, 6}
T = circulant_tournament(n, S)

# Check it's a tournament
for i in range(n):
    for j in range(i+1, n):
        assert T[(i,j)] + T[(j,i)] == 1, f"Not a tournament at ({i},{j})"

print(f"n={n}: Circulant with S={sorted(S)}")
print(f"Verified: IS a tournament")

# By circulant symmetry, check NONHAM(0, b) for non-edges
a = 0
t0 = time.time()
for b in range(1, n):
    if T[(a,b)] == 1:
        continue  # Skip edges (trivially NONHAM=0)

    U = [v for v in range(n) if v != a and v != b]
    M_ab = 0
    for mask in range(1 << len(U)):
        S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
        R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
        sign = (-1) ** len(S_list)
        S_set = sorted(set(S_list) | {a})
        R_set = sorted(set(R) | {b})
        ea = count_paths_ending_at(T, S_set, a)
        bb = count_paths_starting_at(T, R_set, b)
        M_ab += sign * ea * bb

    elapsed = time.time() - t0
    print(f"  NONHAM(0,{b}) = M[0,{b}] = {M_ab} ({'PASS' if M_ab == 0 else 'FAIL'}) [{elapsed:.1f}s]")

    if M_ab != 0:
        print("  *** NONHAM DOES NOT VANISH! ***")
        break

print(f"\nTotal time: {time.time() - t0:.1f}s")
print("DONE")
