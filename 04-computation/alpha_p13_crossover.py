#!/usr/bin/env python3
"""
alpha_p13_crossover.py -- Why does Interval win at p=13?

At p=7,11 (3 mod 4): Paley wins via higher directed-cycle multiplicity.
At p=13 (1 mod 4): Interval wins. WHY?

Key insight from p=11: H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 EXACTLY.
At p=13: max independent set size = floor(13/3) = 4? Actually:
  3+3+3+3 = 12 < 13, so 4 disjoint 3-cycles CAN fit.
  3+5+5 = 13 = p, so a (3,5,5) triple uses ALL vertices -- impossible
  unless cycles are exactly disjoint AND cover Z_13.
  Actually 3+5+5=13, so we'd need 3 disjoint cycles covering all 13 vertices.

So at p=13, alpha_4 > 0 is possible! The truncation level increases.

This script computes:
1. Directed cycle counts and multiplicities at p=13 (k=3,5 only, k=7+ too expensive)
2. alpha_1 and alpha_2 for directed cycles at k=3,5
3. Alpha_3 from (3,3,3) triples and (3,3,5) triples
4. Alpha_4 from (3,3,3,3) quadruples (possible since 3*4=12<13)

For p=13, exhaustive Held-Karp is feasible (2^13 * 13 ~ 100k), so we also
compute H directly.

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


def count_ham_paths(A, n):
    """Held-Karp DP for Hamiltonian path count."""
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def count_directed_ham_cycles(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return A[a][b]*A[b][c]*A[c][a] + A[a][c]*A[c][b]*A[b][a]
    start = 0
    dp = {}
    dp[(1 << start, start)] = 1
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
        if key in dp and dp[key] > 0 and A[verts[v]][verts[start]]:
            total += dp[key]
    return total


def main():
    print("=" * 70)
    print("ALPHA DECOMPOSITION AT p=13: WHY INTERVAL WINS")
    print("=" * 70)

    p = 13
    m = (p - 1) // 2

    # No Paley at p=13 (1 mod 4), but let's compare Interval with
    # a few other circulant tournaments
    S_int = list(range(1, m + 1))  # [1,2,3,4,5,6]

    # Also try the "opposite" -- anti-interval
    # S_anti = list(range(m+1, p))  # [7,8,9,10,11,12] = complement
    # But anti-interval gives same H (complement symmetry)

    # Let's find QR set for comparison (not a valid Paley, but still a connection set)
    S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
    print(f"QR at p=13: {S_qr}")

    tournaments = [
        ("Interval", S_int),
        ("QR (not tournament)", S_qr),  # QR has 6 elements, valid connection set
    ]

    # Compute H directly
    print(f"\np={p}, m={m}")

    for name, S in tournaments:
        A = build_adj(p, S)

        # Verify it's a tournament
        is_tournament = True
        for i in range(p):
            for j in range(p):
                if i != j:
                    if A[i][j] + A[j][i] != 1:
                        is_tournament = False
        if not is_tournament:
            print(f"\n  {name}: NOT a tournament (S={S})")
            continue

        t0 = time.time()
        H = count_ham_paths(A, p)
        t1 = time.time()

        print(f"\n{'='*60}")
        print(f"  {name} (S={S}): H = {H} ({t1-t0:.1f}s)")
        print(f"{'='*60}")

        # Directed cycle enumeration (k=3,5,7 only)
        directed_by_k = {}
        vertex_sets_by_k = {}

        for k in [3, 5, 7]:
            vs = []
            total_dc = 0
            mult_dist = defaultdict(int)

            t2 = time.time()
            for subset in combinations(range(p), k):
                d = count_directed_ham_cycles(A, list(subset))
                if d > 0:
                    vs.append((frozenset(subset), d))
                    total_dc += d
                    mult_dist[d] += 1
            t3 = time.time()

            directed_by_k[k] = total_dc
            vertex_sets_by_k[k] = vs

            print(f"\n    k={k}: {len(vs)} vertex sets, {total_dc} directed cycles ({t3-t2:.1f}s)")
            avg_mult = total_dc / len(vs) if vs else 0
            print(f"      avg multiplicity = {avg_mult:.2f}")
            print(f"      multiplicity dist: {dict(sorted(mult_dist.items()))}")

        # Now build directed cycle list for k=3,5
        all_dc_35 = []
        for k in [3, 5]:
            for fs, d in vertex_sets_by_k[k]:
                for _ in range(d):
                    all_dc_35.append(fs)

        n_35 = len(all_dc_35)
        print(f"\n    Total k=3,5 directed cycles: {n_35}")

        # alpha_2 for k=3,5
        t4 = time.time()
        alpha_2_35 = 0
        for i in range(n_35):
            for j in range(i + 1, n_35):
                if not (all_dc_35[i] & all_dc_35[j]):
                    alpha_2_35 += 1
        t5 = time.time()
        print(f"    alpha_2 (k=3,5 only) = {alpha_2_35} ({t5-t4:.1f}s)")

        # alpha_3 for k=3,5
        print(f"    Computing alpha_3 (k=3,5)...")
        t6 = time.time()
        alpha_3_35 = 0
        # Pre-compute non-adjacent
        non_adj = [[] for _ in range(n_35)]
        for i in range(n_35):
            for j in range(i + 1, n_35):
                if not (all_dc_35[i] & all_dc_35[j]):
                    non_adj[i].append(j)

        for i in range(n_35):
            for j in non_adj[i]:
                for l in non_adj[j]:
                    if l <= j:
                        continue
                    if not (all_dc_35[i] & all_dc_35[l]):
                        alpha_3_35 += 1
        t7 = time.time()
        print(f"    alpha_3 (k=3,5 only) = {alpha_3_35} ({t7-t6:.1f}s)")

        # alpha_4 for k=3 only (4 disjoint 3-cycles: 12 vertices out of 13)
        dc_3 = []
        for fs, d in vertex_sets_by_k[3]:
            for _ in range(d):
                dc_3.append(fs)
        n_3 = len(dc_3)

        print(f"    Computing alpha_4 from k=3 only ({n_3} directed 3-cycles)...")
        t8 = time.time()

        # Non-adjacent pairs for 3-cycles
        non_adj_3 = [[] for _ in range(n_3)]
        for i in range(n_3):
            for j in range(i + 1, n_3):
                if not (dc_3[i] & dc_3[j]):
                    non_adj_3[i].append(j)

        # Mutually disjoint quadruples
        alpha_4_333 = 0
        for i in range(n_3):
            for j in non_adj_3[i]:
                used_ij = dc_3[i] | dc_3[j]
                for l_idx, l in enumerate(non_adj_3[j]):
                    if l <= j:
                        continue
                    if dc_3[i] & dc_3[l]:
                        continue
                    used_ijl = used_ij | dc_3[l]
                    for m_idx, m_val in enumerate(non_adj_3[l]):
                        if m_val <= l:
                            continue
                        if dc_3[i] & dc_3[m_val]:
                            continue
                        if dc_3[j] & dc_3[m_val]:
                            continue
                        alpha_4_333 += 1
        t9 = time.time()
        print(f"    alpha_4 (3,3,3,3 only) = {alpha_4_333} ({t9-t8:.1f}s)")

        # H decomposition
        alpha_1 = sum(directed_by_k[k] for k in [3, 5])
        # Add k=7 to alpha_1
        alpha_1_full = alpha_1 + directed_by_k.get(7, 0)
        print(f"\n    alpha_1 (k=3,5) = {alpha_1}")
        print(f"    alpha_1 (k=3,5,7) = {alpha_1_full}")

        H_partial = 1 + 2*alpha_1 + 4*alpha_2_35 + 8*alpha_3_35 + 16*alpha_4_333
        print(f"\n    H estimate (k=3,5 only, through alpha_4):")
        print(f"      1 + 2*{alpha_1} + 4*{alpha_2_35} + 8*{alpha_3_35} + 16*{alpha_4_333}")
        print(f"      = 1 + {2*alpha_1} + {4*alpha_2_35} + {8*alpha_3_35} + {16*alpha_4_333}")
        print(f"      = {H_partial}")
        print(f"    Actual H = {H}")
        print(f"    Missing from k>=7 and cross-terms = {H - H_partial}")

    # ====== ALL CIRCULANT TOURNAMENTS AT p=13 ======
    print(f"\n{'='*70}")
    print("ALL CIRCULANT TOURNAMENTS AT p=13")
    print("=" * 70)

    # Generate all circulant tournaments
    pairs = [(s, p - s) for s in range(1, m + 1)]
    H_values = []

    for bits in range(1 << m):
        S = []
        for i, (a, b) in enumerate(pairs):
            S.append(a if (bits & (1 << i)) else b)
        S = sorted(S)
        A = build_adj(p, S)
        H = count_ham_paths(A, p)
        is_int = (S == S_int)
        H_values.append((H, S, is_int))

    H_values.sort(reverse=True)
    print(f"\n  {len(H_values)} circulant tournaments")
    print(f"\n  Top 10:")
    for i, (H, S, is_int) in enumerate(H_values[:10]):
        label = " <-- INTERVAL" if is_int else ""
        print(f"    {i+1}. H={H:>12}, S={S}{label}")

    print(f"\n  Bottom 5:")
    for i, (H, S, is_int) in enumerate(H_values[-5:]):
        label = " <-- INTERVAL" if is_int else ""
        print(f"    {len(H_values)-4+i}. H={H:>12}, S={S}{label}")

    # Distinct H values
    distinct = sorted(set(h for h, _, _ in H_values), reverse=True)
    print(f"\n  {len(distinct)} distinct H values:")
    for h in distinct:
        count = sum(1 for hh, _, _ in H_values if hh == h)
        has_int = any(is_int for hh, _, is_int in H_values if hh == h)
        label = " [Interval orbit]" if has_int else ""
        print(f"    H={h:>12}: {count} tournaments{label}")

    # For each distinct H, compute directed 3-cycle and 5-cycle counts
    print(f"\n  H vs directed cycle count (k=3,5):")
    print(f"  {'H':>12} {'dc_3':>8} {'dc_5':>8} {'vs_5':>8} {'avg_m5':>8} {'alpha_2(3)':>12}")

    for h in distinct:
        # Pick one representative
        S_rep = [S for hh, S, _ in H_values if hh == h][0]
        A = build_adj(p, S_rep)

        dc_3 = 0
        vs_3 = 0
        for subset in combinations(range(p), 3):
            d = count_directed_ham_cycles(A, list(subset))
            if d > 0:
                vs_3 += 1
                dc_3 += d

        dc_5 = 0
        vs_5 = 0
        for subset in combinations(range(p), 5):
            d = count_directed_ham_cycles(A, list(subset))
            if d > 0:
                vs_5 += 1
                dc_5 += d

        avg_m5 = dc_5 / vs_5 if vs_5 > 0 else 0

        # alpha_2 from 3-cycles only
        c3_list = []
        for subset in combinations(range(p), 3):
            d = count_directed_ham_cycles(A, list(subset))
            for _ in range(d):
                c3_list.append(frozenset(subset))
        a2_3 = 0
        for i in range(len(c3_list)):
            for j in range(i+1, len(c3_list)):
                if not (c3_list[i] & c3_list[j]):
                    a2_3 += 1

        print(f"  {h:>12} {dc_3:>8} {dc_5:>8} {vs_5:>8} {avg_m5:>8.2f} {a2_3:>12}")


if __name__ == '__main__':
    main()
