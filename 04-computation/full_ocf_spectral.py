"""
full_ocf_spectral.py — Full OCF decomposition with ALL odd cycles,
and its connection to eigenvalue power sums.

KEY INSIGHT (from session S56c + opus S58-S60):
  - H = I(Omega(T), 2) where Omega has vertices = ALL directed odd cycles
  - For circulant tournaments on Z_p, c_3 is CONSTANT (all are regular)
  - The variation in H comes from c_5, c_7, ... (higher odd cycle counts)
  - THM-133 (opus): H = (462 - tr(A^4))/2 at p=7 (exact)
  - THM-134 (opus): Schur-concavity dichotomy between p mod 4 classes

THIS SCRIPT:
  1. Compute FULL Omega (all 3,5,7-cycles) at p=7
  2. Verify H = I(Omega, 2) with the CORRECT Omega
  3. Decompose: which cycle lengths contribute most to H variation?
  4. Find the precise relationship between cycle counts and spectral data

Author: kind-pasteur-2026-03-12-S56c
"""

import sys
import time
import cmath
import math
from itertools import combinations, permutations
from collections import defaultdict

sys.path.insert(0, '04-computation')


def all_circulant_tournaments(n):
    pairs = []
    used = set()
    for a in range(1, n):
        if a not in used:
            b = n - a
            if a == b:
                return []
            pairs.append((a, b))
            used.add(a)
            used.add(b)
    results = []
    for bits in range(2 ** len(pairs)):
        S = []
        for i, (a, b) in enumerate(pairs):
            S.append(a if (bits >> i) & 1 else b)
        results.append(tuple(sorted(S)))
    return results


def ham_count_dp(n, S):
    S_set = set(S)
    adj = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in S_set:
                adj[i][j] = True

    full_mask = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v] == 0 or not (mask & (1 << v)):
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if adj[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    return sum(dp[full_mask][v] for v in range(n))


def find_all_directed_cycles(p, S):
    """Find ALL directed cycles in circulant tournament on Z_p.
    Returns dict: length -> list of vertex frozensets."""
    S_set = set(S)
    adj = [[False] * p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in S_set:
                adj[i][j] = True

    cycles_by_len = defaultdict(list)

    # For each odd k from 3 to p, find all directed k-cycles
    for k in range(3, p + 1, 2):
        # Enumerate k-subsets
        cycle_vsets = set()
        for subset in combinations(range(p), k):
            # Count directed Hamiltonian cycles on this k-vertex subtournament
            for perm in permutations(subset):
                # Check if perm[0]->perm[1]->...->perm[k-1]->perm[0] is valid
                # Only consider cycles starting at min vertex to avoid counting rotations
                if perm[0] != min(subset):
                    continue
                is_cycle = True
                for idx in range(k):
                    if not adj[perm[idx]][perm[(idx + 1) % k]]:
                        is_cycle = False
                        break
                if is_cycle:
                    cycle_vsets.add(frozenset(subset))
                    break  # Found one direction, skip rest for this starting vertex
                    # Actually we should count ALL directed Hamiltonian cycles,
                    # but each vertex set can support multiple (for k>3).
                    # However, for Omega we only care about the VERTEX SET.

        # For Omega, each cycle is identified by its vertex set
        # But a vertex set can support MULTIPLE directed cycles (for k >= 5)!
        # The conflict graph Omega has one vertex per DIRECTED CYCLE, not per vertex set.

        # Let me recount properly: count directed cycles, identified uniquely
        for subset in combinations(range(p), k):
            # Find ALL directed Hamiltonian cycles on this subset
            for perm in permutations(subset):
                if perm[0] != min(subset):
                    continue
                is_cycle = True
                for idx in range(k):
                    if not adj[perm[idx]][perm[(idx + 1) % k]]:
                        is_cycle = False
                        break
                if is_cycle:
                    # This directed cycle visits perm[0]->...->perm[k-1]->perm[0]
                    # Canonical form: start at min, go in the direction
                    cycle_vsets.add(frozenset(subset))
                    # Note: for k=3, each 3-cycle vertex set has exactly 1 directed cycle
                    # For k=5, a 5-vertex subtournament can have 0, 1, 2, or 3 directed 5-cycles

        cycles_by_len[k] = list(cycle_vsets)

    return cycles_by_len


def find_all_directed_cycles_precise(p, S):
    """Find ALL directed cycles, distinguishing multiple cycles on same vertex set.
    Returns dict: length -> list of (vertex_frozenset, cycle_tuple).
    For the conflict graph, two cycles conflict iff they share a vertex."""
    S_set = set(S)
    adj = [[False] * p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in S_set:
                adj[i][j] = True

    all_cycles = []

    for k in range(3, p + 1, 2):
        for subset in combinations(range(p), k):
            min_v = min(subset)
            for perm in permutations(subset):
                if perm[0] != min_v:
                    continue
                is_cycle = True
                for idx in range(k):
                    if not adj[perm[idx]][perm[(idx + 1) % k]]:
                        is_cycle = False
                        break
                if is_cycle:
                    all_cycles.append((frozenset(subset), perm, k))

    return all_cycles


def compute_omega_ip(all_cycles):
    """Compute independence polynomial of the conflict graph Omega.
    Two cycles conflict iff they share a vertex.
    Returns list [alpha_0, alpha_1, alpha_2, ...] and I(Omega, 2)."""
    n_cyc = len(all_cycles)
    if n_cyc > 25:
        # Too many for brute force enumeration
        return None, None

    # Build conflict adjacency
    conflict = [[False] * n_cyc for _ in range(n_cyc)]
    for i in range(n_cyc):
        for j in range(i + 1, n_cyc):
            if all_cycles[i][0] & all_cycles[j][0]:  # share vertex
                conflict[i][j] = True
                conflict[j][i] = True

    # Enumerate all independent sets by bitmask
    alphas = [0] * (n_cyc + 1)
    alphas[0] = 1

    for mask in range(1, 1 << n_cyc):
        verts = [i for i in range(n_cyc) if mask & (1 << i)]
        is_ind = True
        for a in range(len(verts)):
            for b in range(a + 1, len(verts)):
                if conflict[verts[a]][verts[b]]:
                    is_ind = False
                    break
            if not is_ind:
                break
        if is_ind:
            alphas[len(verts)] += 1

    I_at_2 = sum(alphas[k] * (2 ** k) for k in range(n_cyc + 1))
    return alphas, I_at_2


def circulant_eigenvalues(n, S):
    omega = cmath.exp(2j * cmath.pi / n)
    return [sum(omega ** (k * s) for s in S) for k in range(n)]


def main():
    print("=" * 70)
    print("FULL OCF DECOMPOSITION WITH ALL ODD CYCLES")
    print("=" * 70)

    # ================================================================
    # p = 7: Full analysis (manageable: c3=14, c5~20, c7=1 => ~35 cycles)
    # ================================================================
    p = 7
    print(f"\n{'=' * 60}")
    print(f"p = {p}: Full Omega with ALL odd cycles")
    print(f"{'=' * 60}")

    all_S = all_circulant_tournaments(p)
    H_vals_seen = set()

    for S in all_S:
        H = ham_count_dp(p, S)
        if H in H_vals_seen:
            continue
        H_vals_seen.add(H)

        print(f"\n  S = {list(S)}, H = {H}")

        # Find all directed cycles
        all_cycles = find_all_directed_cycles_precise(p, S)

        # Count by length
        by_len = defaultdict(int)
        for vset, perm, k in all_cycles:
            by_len[k] += 1
        print(f"    Cycle counts: {dict(sorted(by_len.items()))}")
        print(f"    Total cycles (|Omega|): {len(all_cycles)}")

        # Compute independence polynomial
        alphas, I_at_2 = compute_omega_ip(all_cycles)
        if alphas is not None:
            # Trim trailing zeros
            while len(alphas) > 1 and alphas[-1] == 0:
                alphas.pop()
            print(f"    alpha_k = {alphas}")
            print(f"    I(Omega, 2) = {I_at_2}")
            print(f"    H = {H}, match: {I_at_2 == H}")
        else:
            print(f"    Too many cycles for brute force IP computation")

        # Spectral data
        eigs = circulant_eigenvalues(p, S)
        y2 = [eigs[k].imag ** 2 for k in range(1, (p - 1) // 2 + 1)]
        print(f"    y^2 = {[round(y, 4) for y in y2]}")

        # Traces
        for k in [3, 4, 5, 6, 7]:
            tr_k = sum(e ** k for e in eigs).real
            print(f"    tr(A^{k}) = {tr_k:.4f}")

    # ================================================================
    # p = 5: Verify
    # ================================================================
    p = 5
    print(f"\n{'=' * 60}")
    print(f"p = {p}: Full Omega")
    print(f"{'=' * 60}")

    all_S = all_circulant_tournaments(p)
    for S in all_S[:1]:  # All have same H
        H = ham_count_dp(p, S)
        print(f"\n  S = {list(S)}, H = {H}")
        all_cycles = find_all_directed_cycles_precise(p, S)
        by_len = defaultdict(int)
        for vset, perm, k in all_cycles:
            by_len[k] += 1
        print(f"    Cycle counts: {dict(sorted(by_len.items()))}")
        print(f"    Total cycles: {len(all_cycles)}")
        alphas, I_at_2 = compute_omega_ip(all_cycles)
        if alphas is not None:
            while len(alphas) > 1 and alphas[-1] == 0:
                alphas.pop()
            print(f"    alpha_k = {alphas}")
            print(f"    I(Omega, 2) = {I_at_2}")
            print(f"    Match: {I_at_2 == H}")

    # ================================================================
    # p = 11: Need to be clever about counting
    # ================================================================
    p = 11
    print(f"\n{'=' * 60}")
    print(f"p = {p}: Cycle counts (Omega too large for brute IP)")
    print(f"{'=' * 60}")

    all_S = all_circulant_tournaments(p)
    H_vals_data = {}

    for S in all_S:
        H = ham_count_dp(p, S)
        if H in H_vals_data:
            continue

        print(f"\n  S = {list(S)}, H = {H}")

        # Find cycle counts by length (but not full enumeration of Omega)
        S_set = set(S)
        adj = [[False] * p for _ in range(p)]
        for i in range(p):
            for j in range(p):
                if i != j and (j - i) % p in S_set:
                    adj[i][j] = True

        for k in [3, 5, 7, 9, 11]:
            if k > p:
                break
            c_k = 0
            for subset in combinations(range(p), k):
                min_v = min(subset)
                for perm in permutations(subset):
                    if perm[0] != min_v:
                        continue
                    is_cycle = True
                    for idx in range(k):
                        if not adj[perm[idx]][perm[(idx + 1) % k]]:
                            is_cycle = False
                            break
                    if is_cycle:
                        c_k += 1
            print(f"    c_{k} = {c_k}")

        # Eigenvalue data
        eigs = circulant_eigenvalues(p, S)
        y2 = [eigs[k].imag ** 2 for k in range(1, (p - 1) // 2 + 1)]
        sum_y4 = sum(y ** 2 for y in y2)
        print(f"    sum(y^4) = {sum_y4:.4f}")

        H_vals_data[H] = S

    # ================================================================
    # KEY ANALYSIS: How do cycle counts vary between circulant tournaments?
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"CYCLE COUNT VARIATION ANALYSIS")
    print(f"{'=' * 60}")
    print("""
    Finding: c_3 is CONSTANT for all circulant tournaments on Z_p
    (because all are regular with degree (p-1)/2).

    Question: Is c_5 also constant? Or does it vary?
    If c_5 varies: it could explain the H variation.
    If c_5 is constant: we need c_7 or the ARRANGEMENT structure.
    """)

    # Check c_5 variation at p=7
    p = 7
    all_S = all_circulant_tournaments(p)
    print(f"\n  p = {p}: c_5 variation across {len(all_S)} circulant tournaments")

    for S in all_S:
        S_set = set(S)
        adj = [[False] * p for _ in range(p)]
        for i in range(p):
            for j in range(p):
                if i != j and (j - i) % p in S_set:
                    adj[i][j] = True

        H = ham_count_dp(p, S)

        # c_5
        c5 = 0
        for subset in combinations(range(p), 5):
            min_v = min(subset)
            for perm in permutations(subset):
                if perm[0] != min_v:
                    continue
                is_cycle = True
                for idx in range(5):
                    if not adj[perm[idx]][perm[(idx + 1) % 5]]:
                        is_cycle = False
                        break
                if is_cycle:
                    c5 += 1

        # c_7 (only 1 possible subset)
        c7 = 0
        subset = tuple(range(p))
        for perm in permutations(subset):
            if perm[0] != 0:
                continue
            is_cycle = True
            for idx in range(p):
                if not adj[perm[idx]][perm[(idx + 1) % p]]:
                    is_cycle = False
                    break
            if is_cycle:
                c7 += 1

        # Also count DISTINCT 5-cycle vertex sets
        c5_vsets = set()
        for subset in combinations(range(p), 5):
            min_v = min(subset)
            for perm in permutations(subset):
                if perm[0] != min_v:
                    continue
                is_cycle = True
                for idx in range(5):
                    if not adj[perm[idx]][perm[(idx + 1) % 5]]:
                        is_cycle = False
                        break
                if is_cycle:
                    c5_vsets.add(frozenset(subset))

        print(f"    S={list(S)}: H={H}, c3=14, c5={c5} ({len(c5_vsets)} vsets), c7={c7}")


if __name__ == '__main__':
    main()
