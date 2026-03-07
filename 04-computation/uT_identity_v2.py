#!/usr/bin/env python3
"""
Verify identity: u_T(m) = m^n * I(Omega(T), 2/m^2)

Fix: use correct alpha computation (no nc cutoff).
For small n, use exhaustive independent set enumeration.

opus-2026-03-07-S39
"""
from itertools import permutations, combinations
from collections import defaultdict


def make_tournament(n, bits):
    A = [[0]*n for _ in range(n)]
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    for k, (i,j) in enumerate(edges):
        if bits & (1 << k):
            A[j][i] = 1
        else:
            A[i][j] = 1
    return A


def held_karp(n, adj):
    """Count Hamiltonian paths."""
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n):
        dp[1<<v][v] = 1
    for S in range(1, 1<<n):
        for v in range(n):
            if not (S & (1<<v)): continue
            if dp[S][v] == 0: continue
            for u in range(n):
                if S & (1<<u): continue
                if adj[v][u]:
                    dp[S|(1<<u)][u] += dp[S][v]
    full = (1<<n)-1
    return sum(dp[full][v] for v in range(n))


def compute_uT_coeffs(n, A):
    """u_T(m) coefficients via permutation sum."""
    edge_set = set()
    opp_set = set()
    for i in range(n):
        for j in range(n):
            if i != j:
                if A[i][j]:
                    edge_set.add((i, j))
                else:
                    opp_set.add((i, j))

    coeffs = defaultdict(int)
    for sigma in permutations(range(n)):
        visited = [False] * n
        cycles = []
        for s in range(n):
            if visited[s]:
                continue
            cyc = []
            c = s
            while not visited[c]:
                visited[c] = True
                cyc.append(c)
                c = sigma[c]
            cycles.append(tuple(cyc))

        valid = True
        phi = 0
        has_even = False
        for cyc in cycles:
            if len(cyc) == 1:
                continue
            if len(cyc) % 2 == 0:
                has_even = True
                break
            is_T = all((cyc[i], cyc[(i+1)%len(cyc)]) in edge_set for i in range(len(cyc)))
            is_Top = all((cyc[i], cyc[(i+1)%len(cyc)]) in opp_set for i in range(len(cyc)))
            if not is_T and not is_Top:
                valid = False
                break
            if is_T:
                phi += len(cyc) - 1
        if not valid or has_even:
            continue
        sign = (-1) ** phi
        coeffs[len(cycles)] += sign
    return coeffs


def compute_alpha_direct(n, A):
    """Compute alpha_k directly from OCF: count independent sets in Omega(T).

    Returns dict {k: alpha_k}.
    """
    edge_set = set()
    for i in range(n):
        for j in range(n):
            if i != j and A[i][j]:
                edge_set.add((i, j))

    # Find all directed odd cycles
    all_cycles = []
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            seen = set()
            for p in permutations(verts):
                if all((p[i], p[(i+1)%length]) in edge_set for i in range(length)):
                    min_idx = list(p).index(min(p))
                    canon = tuple(list(p)[min_idx:] + list(p)[:min_idx])
                    if canon not in seen:
                        seen.add(canon)
                        all_cycles.append((frozenset(verts), canon))

    nc = len(all_cycles)

    # Count independent sets by size using backtracking
    # Two cycles conflict if they share a vertex
    alpha = defaultdict(int)
    alpha[0] = 1

    if nc == 0:
        return alpha

    # For efficiency at moderate nc, use vertex-set conflicts
    def backtrack(idx, used_verts, size):
        alpha[size] += 1
        for j in range(idx, nc):
            vset = all_cycles[j][0]
            if not (vset & used_verts):  # no shared vertices
                backtrack(j + 1, used_verts | vset, size + 1)

    backtrack(0, frozenset(), 0)
    # We counted the empty set extra at start, so subtract 1 from alpha[0]
    # Actually no: backtrack starts by counting alpha[0]+=1, then explores.
    # Wait, we initialized alpha[0]=1 and backtrack also adds 1 to alpha[0].
    # Fix: remove the initialization.
    alpha = defaultdict(int)

    def backtrack2(idx, used_verts, size):
        alpha[size] += 1
        for j in range(idx, nc):
            vset = all_cycles[j][0]
            if not (vset & used_verts):
                backtrack2(j + 1, used_verts | vset, size + 1)

    backtrack2(0, frozenset(), 0)
    return alpha


# === n=5 exhaustive verification ===
print("=== n=5: u_T(m) = m^n * I(Omega, 2/m^2) ===")
n = 5
mismatches = 0
total = 0

for bits in range(1 << 10):
    A = make_tournament(n, bits)
    coeffs = compute_uT_coeffs(n, A)
    alpha = compute_alpha_direct(n, A)

    # Check identity at m=1
    uT_1 = sum(coeffs.get(k, 0) for k in range(n+1))
    I_2 = sum(alpha.get(k, 0) * 2**k for k in range(n+1))

    if uT_1 != I_2:
        mismatches += 1
        if mismatches <= 3:
            print(f"  u_T(1)={uT_1} vs I(Omega,2)={I_2} at bits={bits}")

    # Check identity at m=2
    uT_2 = sum(coeffs.get(k, 0) * 2**k for k in range(n+1))
    rhs_2 = 2**n * sum(alpha.get(k, 0) * (2.0/4)**k for k in range(n+1))

    if abs(uT_2 - rhs_2) > 0.01:
        if mismatches <= 3:
            # Hmm, let's check the coefficient identity directly
            # u_T(m) = sum_j c_j * m^j
            # RHS = m^n * sum_k alpha_k * (2/m^2)^k = sum_k alpha_k * 2^k * m^{n-2k}
            rhs_coeffs = defaultdict(int)
            for k, ak in alpha.items():
                deg = n - 2*k
                rhs_coeffs[deg] += ak * 2**k

            print(f"  bits={bits}: u_T coeffs = {dict(sorted(coeffs.items()))}")
            print(f"             RHS coeffs = {dict(sorted(rhs_coeffs.items()))}")
            print(f"             alpha = {dict(sorted(alpha.items()))}")

    total += 1

print(f"u_T(1)=I(Omega,2) mismatches: {mismatches}/{total}")

# More detailed comparison at specific bits
print("\n=== Detailed coefficient comparison at n=5 ===")
for bits in [0, 42, 100, 341]:
    A = make_tournament(n, bits)
    coeffs = compute_uT_coeffs(n, A)
    alpha = compute_alpha_direct(n, A)

    rhs_coeffs = {}
    for k, ak in alpha.items():
        deg = n - 2*k
        rhs_coeffs[deg] = rhs_coeffs.get(deg, 0) + ak * 2**k

    H = held_karp(n, A)

    print(f"\nbits={bits}, H={H}:")
    print(f"  u_T coefficients: {dict(sorted(coeffs.items()))}")
    print(f"  alpha: {dict(sorted(alpha.items()))}")
    print(f"  RHS coefficients (m^n*I(2/m^2)): {dict(sorted(rhs_coeffs.items()))}")

    match = all(coeffs.get(k, 0) == rhs_coeffs.get(k, 0) for k in set(list(coeffs.keys()) + list(rhs_coeffs.keys())))
    print(f"  Exact match: {match}")
