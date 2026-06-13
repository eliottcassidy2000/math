# ⚠️ WARNING: This script uses QR mod p for p ≡ 1 (mod 4), which does NOT
# produce a tournament (S ∩ (-S) ≠ ∅ gives bidirectional edges).
# Results for those primes are INVALID. See MISTAKE-011b.
# Valid Paley tournaments require p ≡ 3 (mod 4).

#!/usr/bin/env python3
"""
Debug: verify OCF (H = I(Omega, 2)) for specific tournaments.
Check whether the issue is vertex sets vs directed cycles.
"""
import sys
from itertools import combinations, permutations
from collections import Counter

def H_matrix(A, n):
    """Count Hamiltonian paths from adjacency matrix."""
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            c = dp[mask][v]
            if not (mask & (1 << v)) or c == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += c
    return sum(dp[full])


def find_all_directed_odd_cycles(A, n):
    """Find ALL directed odd cycles (as ordered tuples).
    Each cycle (v0, v1, ..., vk) means v0->v1->...->vk->v0.
    Canonical: minimum vertex first (removes rotational duplicates)."""
    cycles = []
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            first = verts[0]  # minimum vertex
            for perm in permutations(verts[1:]):
                path = (first,) + perm
                valid = True
                for i in range(length):
                    if not A[path[i]][path[(i+1) % length]]:
                        valid = False
                        break
                if valid:
                    cycles.append(path)
    return cycles


def find_odd_cycle_vsets(A, n):
    """Find vertex sets that support at least one directed odd cycle."""
    result = []
    for size in range(3, n+1, 2):
        for subset in combinations(range(n), size):
            s = size
            sub = list(subset)
            idx = {v: i for i, v in enumerate(sub)}
            sub_A = [[0]*s for _ in range(s)]
            for a in sub:
                for b in sub:
                    if a != b and A[a][b]:
                        sub_A[idx[a]][idx[b]] = 1
            dp = [[0]*s for _ in range(1 << s)]
            dp[1][0] = 1
            for mask in range(1, 1 << s):
                if not (mask & 1): continue
                for v in range(s):
                    c = dp[mask][v]
                    if not (mask & (1 << v)) or c == 0: continue
                    for u in range(s):
                        if mask & (1 << u): continue
                        if sub_A[v][u]:
                            dp[mask | (1 << u)][u] += c
            full = (1 << s) - 1
            if any(dp[full][v] > 0 and sub_A[v][0] for v in range(1, s)):
                result.append(frozenset(subset))
    return result


def conflict_check(c1, c2):
    """Two cycles conflict iff they share a vertex."""
    return bool(set(c1) & set(c2))


def I_omega_2_from_cycles(cycles):
    """Compute I(Omega, 2) where Omega has one vertex per directed cycle."""
    m = len(cycles)
    if m == 0:
        return 1

    adj_bits = [0] * m
    vsets = [frozenset(c) for c in cycles]
    for a in range(m):
        for b in range(a+1, m):
            if vsets[a] & vsets[b]:
                adj_bits[a] |= 1 << b
                adj_bits[b] |= 1 << a

    total = 0
    for mask in range(1 << m):
        ok = True
        seen = 0
        temp = mask
        while temp:
            v = (temp & -temp).bit_length() - 1
            if adj_bits[v] & seen:
                ok = False
                break
            seen |= 1 << v
            temp &= temp - 1
        if ok:
            total += 2 ** bin(mask).count('1')
    return total


def I_omega_2_from_vsets(vsets):
    """Compute I(Omega, 2) where Omega has one vertex per vertex set."""
    m = len(vsets)
    if m == 0:
        return 1

    adj_bits = [0] * m
    for a in range(m):
        for b in range(a+1, m):
            if vsets[a] & vsets[b]:
                adj_bits[a] |= 1 << b
                adj_bits[b] |= 1 << a

    total = 0
    for mask in range(1 << m):
        ok = True
        seen = 0
        temp = mask
        while temp:
            v = (temp & -temp).bit_length() - 1
            if adj_bits[v] & seen:
                ok = False
                break
            seen |= 1 << v
            temp &= temp - 1
        if ok:
            total += 2 ** bin(mask).count('1')
    return total


# Test on n=3 cyclic: 0->1->2->0
print("=== n=3 cyclic: 0->1->2->0 ===")
A = [[0,1,0],[0,0,1],[1,0,0]]
h = H_matrix(A, 3)
dc = find_all_directed_odd_cycles(A, 3)
vs = find_odd_cycle_vsets(A, 3)
I_dc = I_omega_2_from_cycles(dc)
I_vs = I_omega_2_from_vsets(vs)
print(f"H={h}, directed cycles={len(dc)}, vertex sets={len(vs)}")
print(f"I(Omega,2) from directed cycles = {I_dc}")
print(f"I(Omega,2) from vertex sets = {I_vs}")
print(f"Directed cycles: {dc}")
print()

# Test on n=5 with score (2,2,2,2,2) - regular
print("=== n=5 regular: 0->1->2->3->4->0, 0->2, 1->3, 2->4, 3->0, 4->1 ===")
# Paley T_5: QR mod 5 = {1,4}
A5 = [[0]*5 for _ in range(5)]
qr = {1, 4}
for i in range(5):
    for j in range(5):
        if i != j and (j - i) % 5 in qr:
            A5[i][j] = 1
h5 = H_matrix(A5, 5)
dc5 = find_all_directed_odd_cycles(A5, 5)
vs5 = find_odd_cycle_vsets(A5, 5)
I_dc5 = I_omega_2_from_cycles(dc5)
I_vs5 = I_omega_2_from_vsets(vs5)
print(f"H={h5}, directed cycles={len(dc5)}, vertex sets={len(vs5)}")
print(f"I(Omega,2) from directed cycles = {I_dc5}")
print(f"I(Omega,2) from vertex sets = {I_vs5}")
print(f"3-cycle directed: {sum(1 for c in dc5 if len(c)==3)}")
print(f"5-cycle directed: {sum(1 for c in dc5 if len(c)==5)}")
print(f"3-cycle vsets: {sum(1 for v in vs5 if len(v)==3)}")
print(f"5-cycle vsets: {sum(1 for v in vs5 if len(v)==5)}")
print()

# Now test on several n=6 tournaments
print("=== Systematic test: OCF correctness at n=6 ===")
n = 6
import random
random.seed(42)
errors = 0
for trial in range(200):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    h = H_matrix(A, n)
    dc = find_all_directed_odd_cycles(A, n)
    vs = find_odd_cycle_vsets(A, n)
    I_dc = I_omega_2_from_cycles(dc)
    I_vs = I_omega_2_from_vsets(vs)

    if h != I_dc:
        print(f"  ERROR at trial {trial}: H={h}, I_dc={I_dc}, I_vs={I_vs}, "
              f"#dc={len(dc)}, #vs={len(vs)}")
        errors += 1
    elif h != I_vs:
        print(f"  MISMATCH: H={h}=I_dc, but I_vs={I_vs} differs! "
              f"#dc={len(dc)}, #vs={len(vs)}")

    if h == 21:
        print(f"  H=21 FOUND at trial {trial}!")
        print(f"    #directed_cycles={len(dc)}, #vertex_sets={len(vs)}")
        print(f"    I_dc={I_dc}, I_vs={I_vs}")

if errors == 0:
    print(f"  All {200} tests: H == I(Omega_directed_cycles, 2) -- CORRECT")
else:
    print(f"  {errors} errors found!")

# Check: does the vertex-set version DIFFER from directed-cycle version?
print("\n=== Check: vertex-set I vs directed-cycle I ===")
mismatches = 0
for trial in range(500):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    dc = find_all_directed_odd_cycles(A, n)
    vs = find_odd_cycle_vsets(A, n)
    I_dc = I_omega_2_from_cycles(dc)
    I_vs = I_omega_2_from_vsets(vs)

    if I_dc != I_vs:
        mismatches += 1
        if mismatches <= 5:
            h = H_matrix(A, n)
            print(f"  Mismatch #{mismatches}: H={h}, I_dc={I_dc}, I_vs={I_vs}, "
                  f"#dc={len(dc)}, #vs={len(vs)}")

print(f"\n  Mismatches: {mismatches}/500")
print(f"  NOTE: Omega should use DIRECTED CYCLES as vertices (one node per cycle).")
print(f"  Using vertex sets instead of directed cycles is WRONG if a vertex set")
print(f"  supports multiple cycles.")
sys.stdout.flush()
