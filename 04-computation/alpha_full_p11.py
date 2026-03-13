#!/usr/bin/env python3
"""
alpha_full_p11.py -- Full alpha_j decomposition at p=11

Compute alpha_j for ALL j for Interval and Paley at p=11.
The conflict graph Omega has ~848 vertices (Interval) or ~793 (Paley).
This is too large for brute-force enumeration of all independent sets,
but we can use the backtracking algorithm with pruning.

Key question: HOW MUCH do alpha_3, alpha_4, ... contribute to H?
And does Interval's advantage come from alpha_2 (as we thought)
or from higher-order terms?

At p=7: H(Int)=175, H(Pal)=189, alpha_2(Int)=14, alpha_2(Pal)=7
  H(Int) = 1 + 2*36 + 4*14 = 129. So "higher" = 175-129 = 46
  H(Pal) = 1 + 2*36 + 4*7  = 101. So "higher" = 189-101 = 88
  So Paley's HIGHER-ORDER advantage (88-46=42) overwhelms Interval's alpha_2 advantage (14-7=7)*4=28
  Net: Paley wins by 189-175=14 = 88-46-(28) = 14. Checks out.

At p=11: H(Int)=93027, H(Pal)=95095
  alpha_1(Int)=848, alpha_2(Int)=4972
  alpha_1(Pal)=793, alpha_2(Pal)=3861
  H_2(Int) = 1+1696+19888 = 21585
  H_2(Pal) = 1+1586+15444 = 17031
  Remainder(Int) = 93027 - 21585 = 71442
  Remainder(Pal) = 95095 - 17031 = 78064

  So Interval leads by 4554 at level 2, but Paley leads by 6622 at level 3+.
  Net: Paley wins by 95095-93027 = 2068 = 6622-4554.

  PALEY'S HIGHER-ORDER (alpha_3+) ADVANTAGE IS WHAT MAKES IT WIN!

Author: kind-pasteur-2026-03-12-S59c
"""

import time
from itertools import combinations
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def has_ham_cycle(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b] * A[b][c] * A[c][a] + A[a][c] * A[c][b] * A[b][a]) > 0
    dp = set()
    dp.add((1 << 0, 0))
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            if (mask, v) not in dp:
                continue
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    dp.add((mask | (1 << w), w))
    full = (1 << k) - 1
    for v in range(1, k):
        if (full, v) in dp and A[verts[v]][verts[0]]:
            return True
    return False


def enumerate_cycle_vertex_sets(A, p):
    """All odd-cycle vertex sets."""
    by_k = {}
    for k in range(3, p + 1, 2):
        sets_k = []
        for subset in combinations(range(p), k):
            if has_ham_cycle(A, list(subset)):
                sets_k.append(frozenset(subset))
        by_k[k] = sets_k
    return by_k


def independence_polynomial_backtrack(cycles, x=2):
    """Compute I(Omega, x) using backtracking.

    For each independent set, accumulate x^|S|.
    Uses neighbor bitmask for pruning.
    """
    n = len(cycles)

    # Build neighbor bitmasks
    # Use integer array for speed
    nbr = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if cycles[i] & cycles[j]:  # share a vertex
                nbr[i] |= (1 << j)
                nbr[j] |= (1 << i)

    # Count by size
    alpha = defaultdict(int)

    def backtrack(v, excluded, size):
        alpha[size] += 1
        for w in range(v, n):
            if not (excluded & (1 << w)):
                backtrack(w, excluded | nbr[w], size + 1)

    backtrack(-1, 0, 0)

    # Compute I(Omega, x)
    max_j = max(alpha.keys())
    result = sum(alpha[j] * (x ** j) for j in range(max_j + 1))
    return result, dict(sorted(alpha.items()))


def independence_poly_chunked(cycles, x=2):
    """Compute I(Omega, x) for larger cycle sets using integer set tricks.

    Split cycles into groups and use the polynomial multiplication:
    I(G, x) = I(G[A], x) * I(G[B], x) if A and B are disconnected.

    For connected graphs, use the standard backtracking.
    """
    n = len(cycles)

    if n > 25:
        # Too large for bitmask. Use set-based approach.
        nbr = [set() for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                if cycles[i] & cycles[j]:
                    nbr[i].add(j)
                    nbr[j].add(i)

        alpha = defaultdict(int)

        def backtrack(v, excluded, size):
            alpha[size] += 1
            for w in range(v, n):
                if w not in excluded:
                    backtrack(w, excluded | nbr[w], size + 1)

        backtrack(0, set(), 0)

        max_j = max(alpha.keys()) if alpha else 0
        result = sum(alpha[j] * (x ** j) for j in range(max_j + 1))
        return result, dict(sorted(alpha.items()))
    else:
        return independence_polynomial_backtrack(cycles, x)


def compute_alpha_3_by_type(by_k, p):
    """Compute alpha_3 decomposed by cycle-length triple (k1, k2, k3)."""
    all_cycles = []
    for k in sorted(by_k):
        for fs in by_k[k]:
            all_cycles.append((fs, k))

    n = len(all_cycles)
    result = defaultdict(int)
    total = 0

    for i in range(n):
        for j in range(i + 1, n):
            if all_cycles[i][0] & all_cycles[j][0]:
                continue
            used = all_cycles[i][0] | all_cycles[j][0]
            for l in range(j + 1, n):
                if not (used & all_cycles[l][0]) and \
                   not (all_cycles[j][0] & all_cycles[l][0]):  # redundant but clear
                    ks = tuple(sorted([all_cycles[i][1], all_cycles[j][1], all_cycles[l][1]]))
                    result[ks] += 1
                    total += 1

    return total, dict(sorted(result.items()))


def main():
    print("=" * 70)
    print("FULL ALPHA DECOMPOSITION AT p=7 AND p=11")
    print("=" * 70)

    # ====== p=7: EXACT ======
    p = 7
    m = (p - 1) // 2
    S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
    S_int = list(range(1, m + 1))

    print(f"\n{'='*70}")
    print(f"p={p}, m={m}")
    print(f"{'='*70}")

    for name, S in [("Interval", S_int), ("Paley", S_qr)]:
        A = build_adj(p, S)
        by_k = enumerate_cycle_vertex_sets(A, p)

        all_cycles = []
        for k in sorted(by_k):
            for fs in by_k[k]:
                all_cycles.append(fs)

        n = len(all_cycles)
        print(f"\n  {name}: {n} cycle vertex sets")
        for k in sorted(by_k):
            print(f"    c_{k} = {len(by_k[k])}")

        t0 = time.time()
        H, alpha = independence_polynomial_backtrack(all_cycles, x=2)
        t1 = time.time()

        print(f"\n  I(Omega, 2) = H = {H} ({t1-t0:.2f}s)")
        print(f"  Alpha coefficients:")
        for j in sorted(alpha):
            contribution = alpha[j] * (2**j)
            pct = 100 * contribution / H
            print(f"    alpha_{j} = {alpha[j]:>8} (contributes {contribution:>10} = {pct:>6.2f}%)")

    # ====== p=11: COMPUTE alpha_1, alpha_2, alpha_3 by type ======
    p = 11
    m = (p - 1) // 2
    S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
    S_int = list(range(1, m + 1))

    print(f"\n{'='*70}")
    print(f"p={p}, m={m}")
    print(f"{'='*70}")

    for name, S in [("Interval", S_int), ("Paley", S_qr)]:
        A = build_adj(p, S)
        by_k = enumerate_cycle_vertex_sets(A, p)

        all_cycles = []
        for k in sorted(by_k):
            for fs in by_k[k]:
                all_cycles.append(fs)

        n = len(all_cycles)
        print(f"\n  {name}: {n} cycle vertex sets")
        for k in sorted(by_k):
            print(f"    c_{k} = {len(by_k[k])}")

        # alpha_1 = n
        alpha_1 = n

        # alpha_2: count disjoint pairs
        t0 = time.time()
        alpha_2 = 0
        alpha_2_by_type = defaultdict(int)
        for i in range(n):
            for j in range(i + 1, n):
                if not (all_cycles[i] & all_cycles[j]):
                    alpha_2 += 1
        t1 = time.time()

        print(f"\n  alpha_1 = {alpha_1}")
        print(f"  alpha_2 = {alpha_2} ({t1-t0:.1f}s)")

        # alpha_3 by type
        print(f"  Computing alpha_3...")
        t2 = time.time()
        alpha_3, a3_by_type = compute_alpha_3_by_type(by_k, p)
        t3 = time.time()

        print(f"  alpha_3 = {alpha_3} ({t3-t2:.1f}s)")
        print(f"  alpha_3 by cycle type:")
        for ks, cnt in a3_by_type.items():
            print(f"    {ks}: {cnt}")

        # H estimate
        H_known = {'Interval': 93027, 'Paley': 95095}
        H = H_known.get(name, 0)

        H_3 = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3
        remainder = H - H_3

        print(f"\n  H decomposition:")
        print(f"    1                   = 1")
        print(f"    2*alpha_1 = 2*{alpha_1} = {2*alpha_1}")
        print(f"    4*alpha_2 = 4*{alpha_2} = {4*alpha_2}")
        print(f"    8*alpha_3 = 8*{alpha_3} = {8*alpha_3}")
        print(f"    H_3 = {H_3}")
        print(f"    H = {H}")
        print(f"    Remainder (16*alpha_4 + ...) = {remainder}")

        # Percentage breakdown
        print(f"\n  Percentage breakdown:")
        print(f"    Level 0:  {100*1/H:.4f}%")
        print(f"    Level 1:  {100*2*alpha_1/H:.4f}%")
        print(f"    Level 2:  {100*4*alpha_2/H:.4f}%")
        print(f"    Level 3:  {100*8*alpha_3/H:.4f}%")
        print(f"    Level 4+: {100*remainder/H:.4f}%")

    # ====== COMPARISON ======
    print(f"\n{'='*70}")
    print("INTERVAL vs PALEY COMPARISON")
    print("=" * 70)

    # Recompute for comparison
    results = {}
    for name, S in [("Interval", S_int), ("Paley", S_qr)]:
        A = build_adj(11, S)
        by_k = enumerate_cycle_vertex_sets(A, 11)

        all_cycles = []
        for k in sorted(by_k):
            for fs in by_k[k]:
                all_cycles.append(fs)

        n = len(all_cycles)

        alpha_2 = 0
        for i in range(n):
            for j in range(i + 1, n):
                if not (all_cycles[i] & all_cycles[j]):
                    alpha_2 += 1

        alpha_3, _ = compute_alpha_3_by_type(by_k, 11)

        H = {'Interval': 93027, 'Paley': 95095}[name]
        results[name] = {
            'alpha_1': n,
            'alpha_2': alpha_2,
            'alpha_3': alpha_3,
            'H': H,
        }

    I = results['Interval']
    P = results['Paley']

    print(f"\n  {'Quantity':>20} {'Interval':>12} {'Paley':>12} {'Delta (I-P)':>12} {'Effect on H':>12}")
    print(f"  {'-'*68}")

    d1 = I['alpha_1'] - P['alpha_1']
    d2 = I['alpha_2'] - P['alpha_2']
    d3 = I['alpha_3'] - P['alpha_3']

    I_H3 = 1 + 2*I['alpha_1'] + 4*I['alpha_2'] + 8*I['alpha_3']
    P_H3 = 1 + 2*P['alpha_1'] + 4*P['alpha_2'] + 8*P['alpha_3']
    I_rem = I['H'] - I_H3
    P_rem = P['H'] - P_H3
    d_rem = I_rem - P_rem

    print(f"  {'alpha_1':>20} {I['alpha_1']:>12} {P['alpha_1']:>12} {d1:>12} {2*d1:>12}")
    print(f"  {'alpha_2':>20} {I['alpha_2']:>12} {P['alpha_2']:>12} {d2:>12} {4*d2:>12}")
    print(f"  {'alpha_3':>20} {I['alpha_3']:>12} {P['alpha_3']:>12} {d3:>12} {8*d3:>12}")
    print(f"  {'remainder (4+)':>20} {I_rem:>12} {P_rem:>12} {d_rem:>12} {d_rem:>12}")
    print(f"  {'H':>20} {I['H']:>12} {P['H']:>12} {I['H']-P['H']:>12}")

    # Check
    total_effect = 2*d1 + 4*d2 + 8*d3 + d_rem
    print(f"\n  Check: 2*d1 + 4*d2 + 8*d3 + d_rem = {total_effect} (should be {I['H']-P['H']})")

    print(f"\n  INSIGHT: At p=11, Interval advantage at levels 1-3: {2*d1 + 4*d2 + 8*d3}")
    print(f"  Paley advantage at level 4+: {-d_rem}")
    print(f"  Paley wins by {P['H'] - I['H']} because level 4+ dominance.")


if __name__ == '__main__':
    main()
