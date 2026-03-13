#!/usr/bin/env python3
"""
multiplicity_anatomy.py -- Why does Paley have more directed cycles per vertex set?

The directed-cycle multiplicity (number of directed Ham cycles on a k-vertex subset)
determines alpha_1, which in turn dominates H for p≡3 mod 4.

Key question: what structural property of Paley's sub-tournaments makes them
have more directed Hamiltonian cycles?

Hypothesis: Paley's sub-tournaments are "more regular" (closer to balanced
out-degrees), and regular tournaments maximize the directed Ham cycle count.

Analysis:
1. Score sequences of k-vertex sub-tournaments vs their directed Ham cycle count
2. How the connection set S affects sub-tournament regularity
3. The relationship between directed cycles, regularity, and additive structure

Author: kind-pasteur-2026-03-12-S59c
"""

import time
from itertools import combinations
from collections import defaultdict
from math import comb


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


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


def sub_tournament_score(A, verts):
    """Compute score sequence (sorted out-degrees) of induced sub-tournament."""
    k = len(verts)
    scores = []
    for v in verts:
        out_deg = sum(A[v][w] for w in verts if w != v)
        scores.append(out_deg)
    return tuple(sorted(scores))


def main():
    print("=" * 70)
    print("DIRECTED CYCLE MULTIPLICITY ANATOMY")
    print("=" * 70)

    for p in [7, 11]:
        m = (p - 1) // 2
        S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        S_int = list(range(1, m + 1))

        print(f"\n{'='*70}")
        print(f"p={p}, m={m}")
        print(f"{'='*70}")

        for name, S in [("Interval", S_int), ("Paley", S_qr)]:
            A = build_adj(p, S)

            print(f"\n  {name} (S={S}):")
            print(f"  {'-'*60}")

            for k in range(3, p + 1, 2):
                print(f"\n    k={k}:")

                # For each k-vertex subset: score, #directed Ham cycles
                score_to_mult = defaultdict(list)
                mult_to_sets = defaultdict(list)
                total_directed = 0
                total_sets = 0

                for subset in combinations(range(p), k):
                    d = count_directed_ham_cycles(A, list(subset))
                    if d > 0:
                        total_sets += 1
                        total_directed += d
                        score = sub_tournament_score(A, list(subset))
                        score_to_mult[score].append(d)
                        mult_to_sets[d].append(subset)

                print(f"      {total_sets} vertex sets, {total_directed} directed cycles")

                # Show score → multiplicity mapping
                print(f"      Score sequence analysis:")
                for score in sorted(score_to_mult):
                    mults = score_to_mult[score]
                    n_sets = len(mults)
                    avg_mult = sum(mults) / len(mults)
                    distinct = sorted(set(mults))
                    # Score regularity: variance of scores
                    scores = list(score)
                    mean_s = sum(scores) / len(scores)
                    var_s = sum((s - mean_s)**2 for s in scores) / len(scores)
                    is_regular = (var_s == 0 and k % 2 == 1)

                    print(f"        score={score}: {n_sets} sets, avg mult={avg_mult:.1f}, "
                          f"mults={distinct}, var_s={var_s:.2f}"
                          f"{' REGULAR' if is_regular else ''}")

                # Multiplicity distribution
                print(f"      Multiplicity distribution:")
                for mult in sorted(mult_to_sets):
                    n_sets = len(mult_to_sets[mult])
                    pct = 100 * n_sets / total_sets if total_sets > 0 else 0
                    print(f"        mult={mult}: {n_sets} sets ({pct:.1f}%)")

            # ====== TOTAL DIRECTED CYCLES AND ALPHA_1 ======
            print(f"\n    TOTAL alpha_1 = directed cycles:")
            total_alpha_1 = 0
            for k in range(3, p + 1, 2):
                subtotal = 0
                for subset in combinations(range(p), k):
                    subtotal += count_directed_ham_cycles(A, list(subset))
                total_alpha_1 += subtotal
                print(f"      k={k}: {subtotal} directed cycles")
            print(f"      TOTAL alpha_1 = {total_alpha_1}")

    # ====== KEY INSIGHT: REGULARITY → MULTIPLICITY ======
    print(f"\n{'='*70}")
    print("KEY INSIGHT: SCORE REGULARITY → DIRECTED CYCLE MULTIPLICITY")
    print("=" * 70)

    for p in [7, 11]:
        m = (p - 1) // 2
        S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        S_int = list(range(1, m + 1))

        for name, S in [("Interval", S_int), ("Paley", S_qr)]:
            A = build_adj(p, S)

            # At k=5: count how many sub-tournaments are regular
            k = 5
            reg_count = 0
            nonreg_count = 0
            reg_mult = []
            nonreg_mult = []

            for subset in combinations(range(p), k):
                d = count_directed_ham_cycles(A, list(subset))
                if d == 0:
                    continue
                score = sub_tournament_score(A, list(subset))
                if all(s == score[0] for s in score):
                    reg_count += 1
                    reg_mult.append(d)
                else:
                    nonreg_count += 1
                    nonreg_mult.append(d)

            total = reg_count + nonreg_count
            print(f"\n  p={p}, {name}, k={k}:")
            print(f"    Regular sub-tournaments: {reg_count}/{total} ({100*reg_count/total:.1f}%)")
            print(f"      avg mult = {sum(reg_mult)/len(reg_mult):.2f}" if reg_mult else "")
            print(f"    Non-regular: {nonreg_count}/{total} ({100*nonreg_count/total:.1f}%)")
            print(f"      avg mult = {sum(nonreg_mult)/len(nonreg_mult):.2f}" if nonreg_mult else "")

    # ====== GAP STRUCTURE AND MULTIPLICITY ======
    print(f"\n{'='*70}")
    print("GAP STRUCTURE AND DIRECTED CYCLE COUNT")
    print("=" * 70)

    for p in [7, 11]:
        m = (p - 1) // 2
        S_int = list(range(1, m + 1))
        S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        for name, S in [("Interval", S_int), ("Paley", S_qr)]:
            A = build_adj(p, S)

            k = 5
            print(f"\n  p={p}, {name}, k={k}:")

            # For each 5-vertex set, compute the "gap type" (sorted gaps between vertices)
            gap_to_mult = defaultdict(list)
            for subset in combinations(range(p), k):
                d = count_directed_ham_cycles(A, list(subset))
                if d == 0:
                    continue
                verts = sorted(subset)
                gaps = tuple(sorted([
                    (verts[i+1] - verts[i]) % p if i < k-1 else (verts[0] + p - verts[-1]) % p
                    for i in range(k)
                ]))
                gap_to_mult[gaps].append(d)

            print(f"    Gap types (sorted vertex-to-vertex gaps on Z_p circle):")
            for gaps in sorted(gap_to_mult):
                mults = gap_to_mult[gaps]
                avg = sum(mults) / len(mults)
                print(f"      gaps={gaps}: {len(mults)} sets, avg mult={avg:.1f}, "
                      f"mults={sorted(set(mults))}")

    # ====== PALEY MULTIPLICITY UNIFORMITY ======
    print(f"\n{'='*70}")
    print("PALEY MULTIPLICITY UNIFORMITY")
    print("=" * 70)
    print("At p=7, ALL 21 five-vertex sets in Paley have mult=2.")
    print("At p=11, Paley five-vertex sets have mult ∈ {1,2,3}.")
    print("Interval has more mult=1 (less cyclic sub-tournaments).")
    print()
    print("Key hypothesis: Paley's multiplicative closure (QR is closed under")
    print("products) forces sub-tournaments to be more 'balanced', leading to")
    print("higher directed cycle counts. Interval's additive structure")
    print("(S = {1,...,m}) creates more transitive-like sub-tournaments.")

    # ====== QUANTIFY THE MULTIPLICITY EFFECT ======
    print(f"\n{'='*70}")
    print("QUANTIFYING MULTIPLICITY EFFECT ON H")
    print("=" * 70)

    for p in [7, 11]:
        m = (p - 1) // 2
        S_int = list(range(1, m + 1))
        S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        H_known = {7: {'Interval': 175, 'Paley': 189},
                   11: {'Interval': 93027, 'Paley': 95095}}

        print(f"\n  p={p}:")
        for name, S in [("Interval", S_int), ("Paley", S_qr)]:
            A = build_adj(p, S)

            # For each k: count vertex sets and directed cycles
            vs_count = {}
            dc_count = {}
            for k in range(3, p + 1, 2):
                vs = 0
                dc = 0
                for subset in combinations(range(p), k):
                    d = count_directed_ham_cycles(A, list(subset))
                    if d > 0:
                        vs += 1
                        dc += d
                vs_count[k] = vs
                dc_count[k] = dc

            H = H_known[p][name]
            total_vs = sum(vs_count.values())
            total_dc = sum(dc_count.values())

            print(f"\n    {name}:")
            print(f"      {'k':>4} {'vertex_sets':>12} {'directed':>12} {'avg_mult':>10}")
            for k in sorted(vs_count):
                avg = dc_count[k] / vs_count[k] if vs_count[k] > 0 else 0
                print(f"      {k:>4} {vs_count[k]:>12} {dc_count[k]:>12} {avg:>10.2f}")
            print(f"      {'TOT':>4} {total_vs:>12} {total_dc:>12} {total_dc/total_vs:>10.2f}")
            print(f"      H = {H}, alpha_1 = {total_dc}")

            # What if we normalize: if ALL vertex sets had the SAME multiplicity?
            # Then alpha_1 = total_vs * avg_mult
            # Paley and Interval have the SAME vertex sets at k=3 and k=7,9,11
            # (because c_k(T) depends only on cycle counts, which for DRT are equal)
            # Wait — is that true? Let me check.

        # Check if vertex set counts match
        print(f"\n    Vertex set comparison (same for both?):")
        for k in range(3, p + 1, 2):
            A_int = build_adj(p, S_int)
            A_pal = build_adj(p, S_qr)
            vs_int = sum(1 for subset in combinations(range(p), k)
                        if count_directed_ham_cycles(A_int, list(subset)) > 0)
            vs_pal = sum(1 for subset in combinations(range(p), k)
                        if count_directed_ham_cycles(A_pal, list(subset)) > 0)
            print(f"      k={k}: Int={vs_int}, Pal={vs_pal}, same={vs_int==vs_pal}")


if __name__ == '__main__':
    main()
