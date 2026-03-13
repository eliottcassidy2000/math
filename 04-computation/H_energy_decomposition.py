#!/usr/bin/env python3
"""
H_energy_decomposition.py -- kind-pasteur-2026-03-13-S60

At p=11, E determines H uniquely but NOT monotonically:
  E=65: H=95095 (Paley, max)
  E=69: H=92411 (min!)
  E=73: H=93467
  E=85: H=93027

This script decomposes H via the OCF to understand the non-monotonicity:
  H = I(Omega, 2) = sum alpha_j * 2^j
  H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...

Since alpha_2 = disj3 = (p/4)*E + C (THM-156), and alpha_1 is the total
number of directed odd cycles, the non-monotonicity must come from alpha_1
(or higher alphas) varying NON-linearly with E.

Key question: Is there a formula alpha_1 = f(E, other invariants)?
"""

import numpy as np
from itertools import combinations
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i+s)%p] = 1
    return A


def additive_energy(S, p):
    S_set = set(S)
    energy = 0
    for a in S:
        for b in S:
            for c in S:
                d = (a + b - c) % p
                if d in S_set:
                    energy += 1
    return energy


def compute_all_cycle_counts(A, p):
    """Count directed cycles of each odd length."""
    A_np = np.array(A, dtype=np.float64)
    counts = {}
    for k in range(3, p + 1, 2):
        Ak = np.linalg.matrix_power(A_np, k)
        counts[k] = int(round(np.trace(Ak))) // k
    return counts


def compute_H_and_alpha(A, p):
    """Compute H and alpha decomposition."""
    n = p
    # Enumerate ALL directed odd cycles
    cycles = []
    for k in range(3, n + 1, 2):
        for subset in combinations(range(n), k):
            verts = list(subset)
            nc = count_directed_ham_cycles(A, verts)
            for _ in range(nc):
                cycles.append(frozenset(subset))

    # Build conflict graph
    n_cyc = len(cycles)
    adj = [[False] * n_cyc for _ in range(n_cyc)]
    for i in range(n_cyc):
        for j in range(i + 1, n_cyc):
            if cycles[i] & cycles[j]:
                adj[i][j] = True
                adj[j][i] = True

    # Count independent sets by size
    nbr = [0] * n_cyc
    for i in range(n_cyc):
        for j in range(n_cyc):
            if adj[i][j]:
                nbr[i] |= (1 << j)

    alpha = [0] * (n_cyc + 1)

    def backtrack(v, mask, size):
        alpha[size] += 1
        for w in range(v + 1, n_cyc):
            if not (mask & (1 << w)):
                backtrack(w, mask | nbr[w], size + 1)

    if n_cyc <= 30:
        backtrack(-1, 0, 0)
    else:
        # Fallback: just count alpha_0, alpha_1, alpha_2
        alpha[0] = 1
        alpha[1] = n_cyc
        alpha[2] = 0
        for i in range(n_cyc):
            for j in range(i + 1, n_cyc):
                if not adj[i][j]:
                    alpha[2] += 1

    H = sum(alpha[j] * (2**j) for j in range(len(alpha)))
    return H, alpha, cycles


def count_directed_ham_cycles(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        fwd = A[a][b] * A[b][c] * A[c][a]
        bwd = A[a][c] * A[c][b] * A[b][a]
        return fwd + bwd
    start = 0
    dp = {(1 << start, start): 1}
    for mask in range(1, 1 << k):
        if not (mask & (1 << start)):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            cnt = dp[key]
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + cnt
    full = (1 << k) - 1
    total = 0
    for v in range(k):
        if v == start:
            continue
        key = (full, v)
        if key in dp and dp[key] > 0:
            if A[verts[v]][verts[start]]:
                total += dp[key]
    return total


def full_decomposition(p):
    """Full H = sum alpha_j * 2^j decomposition for all orientations."""
    m = (p - 1) // 2
    N = 1 << m

    print(f"\n{'='*70}")
    print(f"H-ENERGY DECOMPOSITION at p={p}")
    print(f"{'='*70}")

    results = []
    for bits in range(N):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))

        A = build_adj(p, S)
        E = additive_energy(S, p)
        ck = compute_all_cycle_counts(A, p)

        # alpha_1 = total directed odd cycles
        alpha_1 = sum(ck.values())

        H, alpha, cycles = compute_H_and_alpha(A, p)

        results.append({
            'bits': bits, 'S': S, 'E': E, 'H': H,
            'ck': ck, 'alpha': alpha, 'alpha_1': alpha_1
        })

    # Group by E
    by_E = defaultdict(list)
    for r in results:
        by_E[r['E']].append(r)

    print(f"\n  Full decomposition (grouped by E):")
    print(f"  {'E':>4} {'H':>8} {'a0':>4} {'a1':>6} {'a2':>6} {'a3':>5} {'a4':>5}  "
          f"{'c3':>5} {'c5':>5} {'c7':>5} {'c9':>6} {'c11':>6}")

    for E in sorted(by_E.keys()):
        group = by_E[E]
        r = group[0]  # all should have same H
        alpha = r['alpha']
        ck = r['ck']

        a_str = [str(alpha[j]) if j < len(alpha) and alpha[j] > 0 else '0' for j in range(5)]
        ck_str = [str(ck.get(k, 0)) for k in range(3, 12, 2)]

        print(f"  {E:4d} {r['H']:8d} {a_str[0]:>4} {a_str[1]:>6} {a_str[2]:>6} "
              f"{a_str[3]:>5} {a_str[4]:>5}  "
              f"{'  '.join(ck_str)}")

    # Verify H = 1 + 2*a1 + 4*a2 + 8*a3 + 16*a4 + ...
    print(f"\n  Verification of OCF decomposition:")
    for E in sorted(by_E.keys()):
        r = by_E[E][0]
        alpha = r['alpha']
        H_check = sum(alpha[j] * (2**j) for j in range(len(alpha)))
        print(f"    E={E}: H={r['H']}, sum(a_j*2^j) = {H_check}, match: {r['H'] == H_check}")

    # Check: which alpha varies most with E?
    print(f"\n  Alpha variation across energy levels:")
    for j in range(5):
        vals = [by_E[E][0]['alpha'][j] if j < len(by_E[E][0]['alpha']) else 0
                for E in sorted(by_E.keys())]
        print(f"    alpha_{j}: {vals}")

    # The key insight: alpha_1 = sum c_k is NOT a linear function of E
    # because c_k for k >= 7 contributes non-linearly through sum D^{2k}
    print(f"\n  Cycle count decomposition:")
    for E in sorted(by_E.keys()):
        r = by_E[E][0]
        ck = r['ck']
        alpha_1 = sum(ck.values())
        print(f"    E={E}: c3={ck.get(3,0)}, c5={ck.get(5,0)}, c7={ck.get(7,0)}, "
              f"c9={ck.get(9,0)}, c11={ck.get(11,0)}, alpha_1={alpha_1}")

    # Decompose H into contributions
    print(f"\n  H contribution analysis:")
    for E in sorted(by_E.keys()):
        r = by_E[E][0]
        alpha = r['alpha']
        contribs = [alpha[j] * (2**j) for j in range(len(alpha)) if alpha[j] > 0]
        total = sum(contribs)
        print(f"    E={E}: H={total} = " +
              " + ".join(f"{alpha[j]}*{2**j}" for j in range(len(alpha))
                         if alpha[j] > 0))

    # Test: is H a polynomial function of (c3, c5, c7, ...)?
    # Since c3 is constant, the variation in H comes from c5, c7, c9, c11
    print(f"\n  Testing H as function of cycle counts:")
    Es = sorted(by_E.keys())
    if len(Es) >= 2:
        H_arr = np.array([by_E[E][0]['H'] for E in Es], dtype=float)
        for k in range(5, p + 1, 2):
            ck_arr = np.array([by_E[E][0]['ck'].get(k, 0) for E in Es], dtype=float)
            if np.std(ck_arr) > 0:
                r = np.corrcoef(H_arr, ck_arr)[0, 1]
                print(f"    Corr(H, c_{k}) = {r:.6f}")

    return results


# ================================================================
# MAIN
# ================================================================

print("=" * 70)
print("H-ENERGY DECOMPOSITION: Understanding non-monotonicity")
print("=" * 70)

# p=7 first (quick)
full_decomposition(7)

# p=11 (main interest)
full_decomposition(11)

print("\nDONE.")
