#!/usr/bin/env python3
"""
disj3_c5_theorem.py -- Prove disj3 = -c5/2 + const(p)

THEOREM (kind-pasteur-2026-03-12-S60):
For any circulant tournament on Z_p, the number of vertex-disjoint 3-cycle
pairs equals -c5_dir/2 + const(p), where c5_dir = # directed 5-cycles.

PROVED NUMERICALLY: p = 7, 11, 13, 17 (zero error, overlap_gauss_bridge.py).

This script:
1. Derives const(p) via combinatorial analysis
2. Exhaustive 5-vertex tournament census to find the underlying identity
3. Tests 4-fold Gauss sum coefficient at p=7,11,13

Author: kind-pasteur-2026-03-12-S60
"""

import math
from itertools import combinations
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def count_c3_vertex_sets(A, p):
    c3_sets = []
    for a, b, c in combinations(range(p), 3):
        if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
            c3_sets.append(frozenset([a, b, c]))
    return list(set(c3_sets))


def count_c5_directed(A, p):
    """Count directed 5-cycles using Held-Karp (FIXED indentation)."""
    total = 0
    for subset in combinations(range(p), 5):
        verts = list(subset)
        n = 5
        start = 0
        dp = {(1 << start, start): 1}
        for mask in range(1, 1 << n):
            if not (mask & (1 << start)):
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                key = (mask, v)
                if key not in dp or dp[key] == 0:
                    continue
                cnt = dp[key]
                for w in range(n):
                    if mask & (1 << w):
                        continue
                    if A[verts[v]][verts[w]]:
                        nkey = (mask | (1 << w), w)
                        dp[nkey] = dp.get(nkey, 0) + cnt
        # Count cycles OUTSIDE the mask loop
        full = (1 << n) - 1
        for v in range(n):
            if v == start:
                continue
            key = (full, v)
            if key in dp and dp[key] > 0:
                if A[verts[v]][verts[start]]:
                    total += dp[key]
    return total


def disjoint_pair_counts(cycle_sets):
    n = len(cycle_sets)
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            if not (cycle_sets[i] & cycle_sets[j]):
                count += 1
    return count


def overlap_counts(cycle_sets):
    n = len(cycle_sets)
    ov = defaultdict(int)
    for i in range(n):
        for j in range(i+1, n):
            o = len(cycle_sets[i] & cycle_sets[j])
            ov[o] += 1
    return ov


def all_orientations_circ(p):
    m = (p - 1) // 2
    orientations = []
    for bits in range(1 << m):
        S = []
        for i in range(m):
            chord = i + 1
            if (bits >> i) & 1:
                S.append(chord)
            else:
                S.append(p - chord)
        orientations.append((bits, S))
    return orientations


def valid_qr_tournament(p):
    """Return a valid QR tournament connection set (only for p≡3 mod 4)."""
    if p % 4 == 3:
        return sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
    else:
        # p≡1 mod 4: use first orientation (all chords forward)
        m = (p - 1) // 2
        return list(range(1, m + 1))


def gauss4(ga, gb, gc, gd, p):
    total = 0
    for t in range(1, p):
        val = 1
        for g in [ga, gb, gc, gd]:
            val *= math.sin(2 * math.pi * g * t / p)
        total += val
    return total


def main():
    print("=" * 70)
    print("DISJ3 = -C5/2 + CONST(P) THEOREM")
    print("=" * 70)

    # ====== PART 1: Verify const(p) using valid circulant tournaments ======
    print("\n" + "=" * 60)
    print("PART 1: Verify const(p) at key primes")
    print("=" * 60)

    constants = {}
    for p in [7, 11, 13, 17, 19, 23]:
        m = (p - 1) // 2
        c3_formula = math.comb(p, 3) - p * math.comb(m, 2)

        # Use first valid orientation to compute const
        S_first = list(range(1, m + 1))  # Interval
        A = build_adj(p, S_first)
        c3_sets = count_c3_vertex_sets(A, p)
        n3 = len(c3_sets)
        c5 = count_c5_directed(A, p)
        disj3 = disjoint_pair_counts(c3_sets)

        const_p = disj3 + c5 / 2
        constants[p] = const_p

        # Cross-check with QR tournament if p≡3 mod 4
        if p % 4 == 3:
            S_qr = valid_qr_tournament(p)
            A_qr = build_adj(p, S_qr)
            c3_qr = count_c3_vertex_sets(A_qr, p)
            c5_qr = count_c5_directed(A_qr, p)
            disj3_qr = disjoint_pair_counts(c3_qr)
            const_check = disj3_qr + c5_qr / 2
            match = abs(const_check - const_p) < 0.01
        else:
            # Use a random other orientation
            S_other = []
            for i in range(m):
                chord = i + 1
                if i % 2 == 0:
                    S_other.append(chord)
                else:
                    S_other.append(p - chord)
            A_other = build_adj(p, S_other)
            c3_oth = count_c3_vertex_sets(A_other, p)
            c5_oth = count_c5_directed(A_other, p)
            disj3_oth = disjoint_pair_counts(c3_oth)
            const_check = disj3_oth + c5_oth / 2
            match = abs(const_check - const_p) < 0.01

        print(f"\n  p={p}, m={m}:")
        print(f"    c3={n3} (formula: {c3_formula})")
        print(f"    Interval: c5={c5}, disj3={disj3}")
        print(f"    const(p) = {const_p}")
        print(f"    Cross-check: {const_check} {'OK' if match else 'MISMATCH'}")

    # ====== PART 2: Pattern in const(p) ======
    print("\n" + "=" * 60)
    print("PART 2: Pattern in const(p)")
    print("=" * 60)

    print(f"\n  {'p':>4} {'m':>4} {'c3':>6} {'C(c3,2)':>10} {'const':>10} {'const/c3':>10} {'const/p':>10}")
    for p in sorted(constants):
        m = (p - 1) // 2
        c3 = math.comb(p, 3) - p * math.comb(m, 2)
        C2 = math.comb(c3, 2)
        C = constants[p]
        print(f"  {p:>4} {m:>4} {c3:>6} {C2:>10} {C:>10.0f} {C/c3:>10.4f} {C/p:>10.4f}")

    # Try to find const(p) formula
    print(f"\n  Trying const(p) = a*C(p,5) + b*C(p,4) + c*C(p,3):")
    # 3-variable system from p=7,11,13
    ps = [7, 11, 13]
    vals = [constants[p] for p in ps]
    C5s = [math.comb(p, 5) for p in ps]
    C4s = [math.comb(p, 4) for p in ps]
    C3s = [math.comb(p, 3) for p in ps]

    # Solve: a*C5 + b*C4 + c*C3 = const
    # 21a + 35b + 35c = 28       (p=7)
    # 462a + 330b + 165c = 792   (p=11)
    # 1287a + 715b + 286c = 2457 (p=13)
    # Use Cramer's rule
    det = (C5s[0] * (C4s[1]*C3s[2] - C4s[2]*C3s[1])
         - C4s[0] * (C5s[1]*C3s[2] - C5s[2]*C3s[1])
         + C3s[0] * (C5s[1]*C4s[2] - C5s[2]*C4s[1]))

    if abs(det) > 0.001:
        a = (vals[0] * (C4s[1]*C3s[2] - C4s[2]*C3s[1])
           - C4s[0] * (vals[1]*C3s[2] - vals[2]*C3s[1])
           + C3s[0] * (vals[1]*C4s[2] - vals[2]*C4s[1])) / det
        b = (C5s[0] * (vals[1]*C3s[2] - vals[2]*C3s[1])
           - vals[0] * (C5s[1]*C3s[2] - C5s[2]*C3s[1])
           + C3s[0] * (C5s[1]*vals[2] - C5s[2]*vals[1])) / det
        c = (C5s[0] * (C4s[1]*vals[2] - C4s[2]*vals[1])
           - C4s[0] * (C5s[1]*vals[2] - C5s[2]*vals[1])
           + vals[0] * (C5s[1]*C4s[2] - C5s[2]*C4s[1])) / det

        print(f"    a = {a:.6f}, b = {b:.6f}, c = {c:.6f}")

        # Verify at all primes
        for p in sorted(constants):
            pred = a * math.comb(p, 5) + b * math.comb(p, 4) + c * math.comb(p, 3)
            err = abs(pred - constants[p])
            print(f"    p={p}: predicted = {pred:.4f}, actual = {constants[p]:.0f}, error = {err:.4f}")

    # Try other basis: const = f(p, m, c3)
    print(f"\n  Try: const = a*c3^2 + b*c3 + c*p:")
    c3_vals = {p: math.comb(p, 3) - p * math.comb((p-1)//2, 2) for p in sorted(constants)}
    ps = sorted(constants)[:3]
    c3s = [c3_vals[p] for p in ps]
    cs = [constants[p] for p in ps]
    pps = list(ps)

    # a*c3^2 + b*c3 + c*p = const
    det2 = (c3s[0]**2 * (c3s[1]*pps[2] - c3s[2]*pps[1])
          - c3s[0] * (c3s[1]**2*pps[2] - c3s[2]**2*pps[1])
          + pps[0] * (c3s[1]**2*c3s[2] - c3s[2]**2*c3s[1]))
    if abs(det2) > 0.001:
        a2 = (cs[0] * (c3s[1]*pps[2] - c3s[2]*pps[1])
            - c3s[0] * (cs[1]*pps[2] - cs[2]*pps[1])
            + pps[0] * (cs[1]*c3s[2] - cs[2]*c3s[1])) / det2
        b2 = (c3s[0]**2 * (cs[1]*pps[2] - cs[2]*pps[1])
            - cs[0] * (c3s[1]**2*pps[2] - c3s[2]**2*pps[1])
            + pps[0] * (c3s[1]**2*cs[2] - c3s[2]**2*cs[1])) / det2
        c2 = (c3s[0]**2 * (c3s[1]*cs[2] - c3s[2]*cs[1])
            - c3s[0] * (c3s[1]**2*cs[2] - c3s[2]**2*cs[1])
            + cs[0] * (c3s[1]**2*c3s[2] - c3s[2]**2*c3s[1])) / det2

        print(f"    a={a2:.6f}, b={b2:.6f}, c={c2:.6f}")
        for p in sorted(constants):
            c3p = c3_vals[p]
            pred = a2 * c3p**2 + b2 * c3p + c2 * p
            err = abs(pred - constants[p])
            print(f"    p={p}: c3={c3p}, pred={pred:.4f}, actual={constants[p]:.0f}, err={err:.4f}")

    # ====== PART 3: Exhaustive 5-vertex tournament census ======
    print("\n" + "=" * 60)
    print("PART 3: Exhaustive 5-vertex tournament census")
    print("=" * 60)

    n = 5
    results_5 = set()

    for bits in range(1 << (n*(n-1)//2)):
        A5 = [[0]*n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if (bits >> idx) & 1:
                    A5[i][j] = 1
                else:
                    A5[j][i] = 1
                idx += 1

        # c3
        c3_5_sets = []
        for a, b, cc in combinations(range(n), 3):
            if (A5[a][b] and A5[b][cc] and A5[cc][a]) or (A5[a][cc] and A5[cc][b] and A5[b][a]):
                c3_5_sets.append(frozenset([a, b, cc]))
        c3_5_sets = list(set(c3_5_sets))
        c3_5 = len(c3_5_sets)

        # c5
        start = 0
        dp = {(1 << start, start): 1}
        for mask in range(1, 1 << n):
            if not (mask & (1 << start)):
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                key = (mask, v)
                if key not in dp or dp[key] == 0:
                    continue
                cnt = dp[key]
                for w in range(n):
                    if mask & (1 << w):
                        continue
                    if A5[v][w]:
                        nkey = (mask | (1 << w), w)
                        dp[nkey] = dp.get(nkey, 0) + cnt
        full = (1 << n) - 1
        c5_5 = 0
        for v in range(n):
            if v == start:
                continue
            key = (full, v)
            if key in dp and dp[key] > 0:
                if A5[v][start]:
                    c5_5 += dp[key]

        # Overlap counts
        ov1_5 = 0
        ov2_5 = 0
        disj_5 = 0
        for i in range(len(c3_5_sets)):
            for j in range(i+1, len(c3_5_sets)):
                overlap = len(c3_5_sets[i] & c3_5_sets[j])
                if overlap == 0:
                    disj_5 += 1
                elif overlap == 1:
                    ov1_5 += 1
                elif overlap == 2:
                    ov2_5 += 1

        # Score sequence
        scores = tuple(sorted([sum(A5[i]) for i in range(n)]))

        results_5.add((c3_5, c5_5, ov1_5, ov2_5, disj_5, scores))

    print(f"\n  Distinct (c3, c5, ov1, ov2, disj, scores) patterns: {len(results_5)}")
    print(f"  {'c3':>3} {'c5':>3} {'ov1':>4} {'ov2':>4} {'disj':>4} {'C(c3,2)':>7} {'c5+2*disj':>10} {'scores':>16}")

    combos = set()
    for c3_5, c5_5, ov1_5, ov2_5, disj_5, scores in sorted(results_5):
        total_pairs = c3_5 * (c3_5 - 1) // 2
        combo = c5_5 + 2 * disj_5
        combos.add(combo)
        print(f"  {c3_5:>3} {c5_5:>3} {ov1_5:>4} {ov2_5:>4} {disj_5:>4} {total_pairs:>7} {combo:>10}  {str(scores):>16}")

    if len(combos) == 1:
        print(f"\n  *** c5 + 2*disj = {list(combos)[0]} IS UNIVERSAL FOR ALL 5-VERTEX TOURNAMENTS! ***")
    else:
        print(f"\n  c5 + 2*disj values: {sorted(combos)}")
        # Check if c5 + 2*ov1 is universal
        combos_ov1 = set()
        for c3_5, c5_5, ov1_5, ov2_5, disj_5, scores in results_5:
            combos_ov1.add(c5_5 + 2 * ov1_5)
        print(f"  c5 + 2*ov1 values: {sorted(combos_ov1)}")

        # What IS universal?
        # Try: a*c5 + b*disj + c*ov1 + d*ov2 = constant?
        # Check c5 + 2*(ov1 + ov2) = c5 + 2*(C(c3,2) - disj)
        combos_alt = set()
        for c3_5, c5_5, ov1_5, ov2_5, disj_5, scores in results_5:
            total_pairs = c3_5 * (c3_5 - 1) // 2
            combos_alt.add(c5_5 + 2 * (total_pairs - disj_5))
        print(f"  c5 + 2*(C(c3,2)-disj) values: {sorted(combos_alt)}")

        # If c5 + 2*disj is NOT constant, is c5 - 2*disj?
        combos_neg = set()
        for c3_5, c5_5, ov1_5, ov2_5, disj_5, scores in results_5:
            combos_neg.add(c5_5 - 2 * disj_5)
        print(f"  c5 - 2*disj values: {sorted(combos_neg)}")

        # Or linear in c3?
        by_c3 = defaultdict(set)
        for c3_5, c5_5, ov1_5, ov2_5, disj_5, scores in results_5:
            by_c3[c3_5].add(c5_5 + 2 * disj_5)
        print(f"\n  c5 + 2*disj grouped by c3:")
        for c3_5 in sorted(by_c3):
            print(f"    c3={c3_5}: {sorted(by_c3[c3_5])}")

    # ====== PART 4: Identity on a single 5-vertex tournament ======
    print("\n" + "=" * 60)
    print("PART 4: Looking for 5-vertex identity")
    print("=" * 60)

    # Let's compute MANY quantities for 5-vertex tournaments
    n = 5
    all_data = []
    for bits in range(1 << (n*(n-1)//2)):
        A5 = [[0]*n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if (bits >> idx) & 1:
                    A5[i][j] = 1
                else:
                    A5[j][i] = 1
                idx += 1

        c3_5_sets = []
        for a, b, cc in combinations(range(n), 3):
            if (A5[a][b] and A5[b][cc] and A5[cc][a]) or (A5[a][cc] and A5[cc][b] and A5[b][a]):
                c3_5_sets.append(frozenset([a, b, cc]))
        c3_5_sets = list(set(c3_5_sets))
        c3_5 = len(c3_5_sets)

        start = 0
        dp = {(1 << start, start): 1}
        for mask in range(1, 1 << n):
            if not (mask & (1 << start)):
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                key = (mask, v)
                if key not in dp or dp[key] == 0:
                    continue
                cnt = dp[key]
                for w in range(n):
                    if mask & (1 << w):
                        continue
                    if A5[v][w]:
                        nkey = (mask | (1 << w), w)
                        dp[nkey] = dp.get(nkey, 0) + cnt
        full = (1 << n) - 1
        c5_5 = 0
        for v in range(n):
            if v == start:
                continue
            key = (full, v)
            if key in dp and dp[key] > 0:
                if A5[v][start]:
                    c5_5 += dp[key]

        disj_5 = 0
        ov1_5 = 0
        ov2_5 = 0
        for i in range(len(c3_5_sets)):
            for j in range(i+1, len(c3_5_sets)):
                overlap = len(c3_5_sets[i] & c3_5_sets[j])
                if overlap == 0:
                    disj_5 += 1
                elif overlap == 1:
                    ov1_5 += 1
                elif overlap == 2:
                    ov2_5 += 1

        scores = tuple(sorted([sum(A5[i]) for i in range(n)]))
        all_data.append((c3_5, c5_5, ov1_5, ov2_5, disj_5, scores))

    # Search for universal identities involving c3, c5, ov1, ov2, disj
    # We know: disj + ov1 + ov2 = C(c3, 2)
    # Search: a*c5 + b*disj + c*ov1 + d*ov2 = constant for ALL tournaments

    print(f"  Searching for a*c5 + b*disj + c*ov1 + d*ov2 = constant...")

    # Group by unique (c3, c5, disj, ov1, ov2) tuples
    unique_data = list(set(all_data))
    print(f"  {len(unique_data)} distinct patterns from {len(all_data)} tournaments")

    # Check: c5 + 2*ov2 = ?
    vals = set()
    for c3_5, c5_5, ov1_5, ov2_5, disj_5, _ in unique_data:
        vals.add(c5_5 + 2 * ov2_5)
    print(f"  c5 + 2*ov2: {sorted(vals)}")

    vals = set()
    for c3_5, c5_5, ov1_5, ov2_5, disj_5, _ in unique_data:
        vals.add(c5_5 + 2 * ov1_5)
    print(f"  c5 + 2*ov1: {sorted(vals)}")

    vals = set()
    for c3_5, c5_5, ov1_5, ov2_5, disj_5, _ in unique_data:
        vals.add(c5_5 - 2 * disj_5)
    print(f"  c5 - 2*disj: {sorted(vals)}")

    vals = set()
    for c3_5, c5_5, ov1_5, ov2_5, disj_5, _ in unique_data:
        vals.add(c5_5 + 4 * ov2_5)
    print(f"  c5 + 4*ov2: {sorted(vals)}")

    # Maybe c5 = f(c3, ov1, ov2) for some f?
    # Since disj = C(c3,2) - ov1 - ov2, we have disj determined by (c3, ov1, ov2)
    # So disj3 + c5/2 = C(c3,2) - ov1 - ov2 + c5/2
    # For this to be constant across orientations of Z_p, we need:
    # c5/2 - ov1 - ov2 = const(p) - C(c3,2) across all orientations

    # On a single 5-vertex tournament:
    # What is c5 - 2*(ov1 + ov2)?
    vals = set()
    for c3_5, c5_5, ov1_5, ov2_5, disj_5, _ in unique_data:
        vals.add(c5_5 - 2 * (ov1_5 + ov2_5))
    print(f"  c5 - 2*(ov1+ov2): {sorted(vals)}")

    # Hmm, maybe the identity works at the level of each 5-subset
    # summed over all C(p,5) subsets, not for individual 5-vertex tournaments

    # Let's verify: sum over C(p,5) subsets of [c5(V) - 2*ov_share1(V)] = 0?
    # Where ov_share1(V) counts pairs of 3-cycles in V sharing 1 vertex
    print(f"\n  Checking sum over 5-subsets:")
    for p in [7, 11]:
        m = (p - 1) // 2
        S_int = list(range(1, m + 1))
        A = build_adj(p, S_int)
        c3_sets = count_c3_vertex_sets(A, p)

        sum_c5 = 0
        sum_ov1 = 0
        sum_ov2 = 0
        sum_disj = 0

        for subset in combinations(range(p), 5):
            V = frozenset(subset)
            verts = list(subset)

            # c5 for this subset
            nn = 5
            start = 0
            dp = {(1 << start, start): 1}
            for mask in range(1, 1 << nn):
                if not (mask & (1 << start)):
                    continue
                for v in range(nn):
                    if not (mask & (1 << v)):
                        continue
                    key = (mask, v)
                    if key not in dp or dp[key] == 0:
                        continue
                    cnt = dp[key]
                    for w in range(nn):
                        if mask & (1 << w):
                            continue
                        if A[verts[v]][verts[w]]:
                            nkey = (mask | (1 << w), w)
                            dp[nkey] = dp.get(nkey, 0) + cnt
            full = (1 << nn) - 1
            c5_local = 0
            for v in range(nn):
                if v == start:
                    continue
                key = (full, v)
                if key in dp and dp[key] > 0:
                    if A[verts[v]][verts[start]]:
                        c5_local += dp[key]

            # 3-cycle pairs within this 5-set
            c3_in = [c for c in c3_sets if c.issubset(V)]
            ov1_local = 0
            ov2_local = 0
            disj_local = 0
            for i in range(len(c3_in)):
                for j in range(i+1, len(c3_in)):
                    overlap = len(c3_in[i] & c3_in[j])
                    if overlap == 0:
                        disj_local += 1
                    elif overlap == 1:
                        ov1_local += 1
                    elif overlap == 2:
                        ov2_local += 1

            sum_c5 += c5_local
            sum_ov1 += ov1_local
            sum_ov2 += ov2_local
            sum_disj += disj_local

        # The global counts
        c5_total = count_c5_directed(A, p)
        disj3_total = disjoint_pair_counts(c3_sets)
        ov_total = overlap_counts(c3_sets)

        print(f"\n  p={p}, Interval:")
        print(f"    sum_c5 (from 5-subsets) = {sum_c5} (global c5 = {c5_total})")
        print(f"    sum_ov1 = {sum_ov1}, sum_ov2 = {sum_ov2}, sum_disj = {sum_disj}")
        print(f"    Global: ov1={ov_total[1]}, ov2={ov_total[2]}, disj={disj3_total}")

        # Each pair of 3-cycles sharing k vertices: the pair sits inside
        # C(p-5+k, 5-5+k) = ... hmm
        # A pair of 3-cycles with overlap=0 uses 6 vertices, sits inside C(p-6, -1)=0 5-subsets (impossible!)
        # A pair with overlap=1 uses 5 vertices, sits inside C(p-5, 0)=1 five-subset
        # A pair with overlap=2 uses 4 vertices, sits inside C(p-4, 1)=p-4 five-subsets

        print(f"    ov1 pairs: each in exactly 1 five-subset. sum_ov1 = ov1_global = {ov_total[1]} {'OK' if sum_ov1 == ov_total[1] else 'FAIL'}")
        print(f"    ov2 pairs: each in {p-4} five-subsets. sum_ov2/(p-4) = {sum_ov2/(p-4):.1f}, ov2_global = {ov_total[2]}")
        print(f"    disj pairs: 0 five-subsets contain a disjoint pair. sum_disj = {sum_disj}")

        # KEY INSIGHT: sum_disj should be 0 because two disjoint 3-cycles use 6 vertices,
        # which can't both fit in a 5-vertex subset!
        if sum_disj == 0:
            print(f"    *** CONFIRMED: no disjoint 3-cycle pairs within any 5-subset ***")

        # So: sum over 5-subsets of c5 = global c5 (each 5-cycle counted once)
        #     sum over 5-subsets of ov1 = global ov1 (each ov1 pair in exactly 1 five-set)
        #     sum over 5-subsets of ov2 = (p-4) * global ov2

        # From the 5-vertex identities:
        # For EACH 5-vertex tournament, c5 and ov1 and ov2 satisfy some identity
        # Let's check: sum c5 = ? as function of sum ov1

        # But c5 and ov1 are summed over ALL 5-subsets.
        # c5_total = sum_{V:|V|=5} c5(V)
        # ov1_total = sum_{V:|V|=5} ov1_in_V(V) = global_ov1 (since ov1 pair in exactly 1 V)
        # ov2_total = sum_{V:|V|=5} ov2_in_V(V) = (p-4) * global_ov2

        # For the identity disj3 = -c5/2 + const:
        # disj3 = C(c3,2) - ov1 - ov2
        # So -c5/2 + const = C(c3,2) - ov1 - ov2
        # => c5 - 2*ov1 - 2*ov2 = 2*C(c3,2) - 2*const (constant!)
        # => c5 - 2*ov1 - 2*ov2 = constant for all orientations

        val = c5_total - 2 * ov_total[1] - 2 * ov_total[2]
        print(f"    c5 - 2*ov1 - 2*ov2 = {val}")

    # ====== PART 5: Verify c5 - 2*ov1 - 2*ov2 = const across all orientations ======
    print("\n" + "=" * 60)
    print("PART 5: Verify c5 - 2*ov1 - 2*ov2 = const for all orientations")
    print("=" * 60)

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        print(f"\n  p={p}, m={m}:")

        vals_combo = set()
        for bits, S in all_orientations_circ(p):
            A = build_adj(p, S)
            c5 = count_c5_directed(A, p)
            c3_sets = count_c3_vertex_sets(A, p)
            ov = overlap_counts(c3_sets)
            combo = c5 - 2 * ov[1] - 2 * ov[2]
            vals_combo.add(combo)

        print(f"    c5 - 2*ov1 - 2*ov2 values: {sorted(vals_combo)}")
        if len(vals_combo) == 1:
            print(f"    *** CONSTANT = {list(vals_combo)[0]} ***")

    print("\nDONE.")


if __name__ == '__main__':
    main()
