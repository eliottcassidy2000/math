#!/usr/bin/env python3
"""
Test NONHAM=0 for Paley tournament T_11.

n=11 is large (11! = ~40M permutations for H), but for NONHAM we only
need to check pair (0,b) with T[0,b]=0, and by circulant symmetry
we only need ONE such pair.

For Paley T_11: gen_set = QR = {1,3,4,5,9} (quadratic residues mod 11).
T[0,b]=1 iff b in QR, T[0,b]=0 iff b in QNR = {2,6,7,8,10}.

So we check NONHAM(0, 2) (since 2 is the smallest QNR).

The subset enumeration is over U = {1,3,4,...,10}\{2} = 9 vertices,
so 2^9 = 512 subsets. Each subset requires path counting on smaller sets.
Max path counting: 10! = 3.6M for the full set. Most subsets are smaller.

This should be feasible but may take minutes.

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

def check_nonham_detailed(T, n, a, b):
    """Check NONHAM(a,b) with detailed output."""
    U = [v for v in range(n) if v != a and v != b]
    total = 0
    nonzero_count = 0

    for mask in range(1 << len(U)):
        S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
        R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
        sign = (-1)**len(S_list)
        S_set = sorted(set(S_list) | {a})
        R_set = sorted(set(R) | {b})

        ea = E_paths(T, S_set, a)
        bb = B_paths(T, R_set, b)

        contrib = sign * ea * bb
        total += contrib

        if ea > 0 and bb > 0:
            nonzero_count += 1

        if mask % 64 == 63:
            elapsed = time.time() - t0
            print(f"    mask {mask+1}/512: total so far = {total} ({elapsed:.1f}s)", flush=True)

    return total, nonzero_count


print("=" * 70)
print("n=11: NONHAM=0 for Paley T_11?")
print("=" * 70)

n = 11
# Quadratic residues mod 11: 1^2=1, 2^2=4, 3^2=9, 4^2=5, 5^2=3
QR = {1, 3, 4, 5, 9}
QNR = {2, 6, 7, 8, 10}
print(f"  QR = {sorted(QR)}")
print(f"  QNR = {sorted(QNR)}")

T = circulant_tournament(n, QR)

# Verify a few edges
print(f"  T[0,1]={T[(0,1)]}, T[0,2]={T[(0,2)]}, T[0,3]={T[(0,3)]}")

a, b = 0, 2  # T[0,2] = 0 (2 is QNR)
print(f"\n  Checking NONHAM({a},{b}) where T[{a},{b}]={T[(a,b)]}")

t0 = time.time()
nh, nz = check_nonham_detailed(T, n, a, b)
elapsed = time.time() - t0

print(f"\n  NONHAM({a},{b}) = {nh}")
print(f"  Nonzero (E*B) pairs: {nz} out of 512 subsets")
print(f"  Time: {elapsed:.1f}s")

if nh == 0:
    print("\n  CONFIRMED: NONHAM=0 for Paley T_11, pair (0,2)")
    print("  By circulant symmetry, NONHAM=0 for ALL pairs (a,b) with T[a,b]=0")
    print("  Therefore M[a,b]=0 for all non-edges, hence M=(H/11)*I")
else:
    print(f"\n  COUNTEREXAMPLE: NONHAM != 0 for Paley T_11!")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
