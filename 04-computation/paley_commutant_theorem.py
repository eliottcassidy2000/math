#!/usr/bin/env python3
"""
paley_commutant_theorem.py -- kind-pasteur-2026-03-21-S13

THE PALEY COMMUTANT THEOREM (conjecture):

Among all tournaments on n vertices, Paley T_p has the LARGEST commutant
of its signed adjacency matrix, AND the maximum H.

dim(Comm(S_{Paley})) = (p^2 - 2p + 3) / 2

Sequence: p=3: 3, p=7: 19, p=11: 51, p=23: 243

TESTING: Does max comm_dim always coincide with max H?
Is the commutant a sufficient certificate for H-maximality?

Also: the SPECTRAL FLATNESS PRINCIPLE:
Among regular tournaments, tr(S^4) is MINIMIZED by Paley (flattest spectrum).
Does min tr(S^4) among regular tournaments always coincide with max H?

Author: kind-pasteur-2026-03-21-S13
"""

import numpy as np
from itertools import combinations, permutations
from collections import defaultdict


def binary_to_tournament(bits, n):
    A = np.zeros((n, n), dtype=int)
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << pos):
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A


def signed_adjacency(A):
    n = A.shape[0]
    S = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                S[i][j] = 2.0 * A[i][j] - 1.0
    return S


def count_ham_paths(A):
    n = A.shape[0]
    dp = defaultdict(int)
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[(mask, v)] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[(mask | (1 << w), w)] += dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp[(full, v)] for v in range(n))


def commutant_dimension(S):
    """
    Dimension of {X : XS = SX} = sum d_k^2 where d_k are eigenspace dims.
    """
    n = S.shape[0]
    S_kron = np.kron(S, np.eye(n)) - np.kron(np.eye(n), S.T)
    return n * n - np.linalg.matrix_rank(S_kron, tol=0.01)


def spectral_kurtosis(S):
    """
    Spectral kurtosis: tr(S^4) / tr(S^2)^2 * n.
    Measures how spread out the eigenvalue magnitudes are.
    Minimum for DRT (flat spectrum), higher for non-DRT.
    """
    tr2 = np.trace(S @ S)
    tr4 = np.trace(np.linalg.matrix_power(S, 4))
    if abs(tr2) < 0.001:
        return 0
    return float(tr4 / (tr2 * tr2) * S.shape[0])


def test_commutant_maximality():
    """
    Test: Does max commutant dimension coincide with max H
    across ALL tournaments (not just regular)?
    """
    print("=" * 70)
    print("TEST: COMMUTANT DIMENSION vs H MAXIMALITY")
    print("=" * 70)

    for n in [3, 5, 7]:
        m = n * (n - 1) // 2
        total = 1 << m

        if total > 3000000:
            print(f"\nn={n}: too many, sampling")
            continue

        print(f"\nn={n} ({total} tournaments):")

        max_H = 0
        max_comm = 0
        max_H_bits = 0
        max_comm_bits = 0
        H_comm_pairs = []

        count = 0
        for bits in range(total):
            A = binary_to_tournament(bits, n)
            H = count_ham_paths(A)
            S = signed_adjacency(A)
            comm = commutant_dimension(S)

            if H > max_H:
                max_H = H
                max_H_bits = bits
            if comm > max_comm:
                max_comm = comm
                max_comm_bits = bits

            H_comm_pairs.append((H, comm))
            count += 1

        # Check if max H tournament also has max comm
        max_H_comms = [c for h, c in H_comm_pairs if h == max_H]
        max_comm_Hs = [h for h, c in H_comm_pairs if c == max_comm]

        print(f"  Max H = {max_H}, comm at max H: {sorted(set(max_H_comms))}")
        print(f"  Max comm = {max_comm}, H at max comm: {sorted(set(max_comm_Hs))}")
        print(f"  COINCIDE: {'YES' if max_H in max_comm_Hs else 'NO'}")

        # Correlation
        Hs = [h for h, c in H_comm_pairs]
        comms = [c for h, c in H_comm_pairs]
        corr = np.corrcoef(Hs, comms)[0, 1]
        print(f"  Correlation(H, comm_dim) = {corr:.4f}")

        # Distribution of comm_dim
        comm_dist = defaultdict(list)
        for h, c in H_comm_pairs:
            comm_dist[c].append(h)
        print(f"\n  comm_dim distribution:")
        for c in sorted(comm_dist.keys()):
            hs = comm_dist[c]
            print(f"    comm={c}: {len(hs)} tours, "
                  f"H range=[{min(hs)}, {max(hs)}], "
                  f"mean H={np.mean(hs):.1f}")


def test_spectral_flatness():
    """
    Test: Among regular tournaments, does minimum tr(S^4) coincide with max H?
    This would confirm the "spectral flatness principle."
    """
    print("\n" + "=" * 70)
    print("TEST: SPECTRAL FLATNESS PRINCIPLE (REGULAR TOURNAMENTS)")
    print("=" * 70)

    for n in [3, 5, 7]:
        m = n * (n - 1) // 2
        total = 1 << m

        print(f"\nn={n}:")

        regular_data = []
        for bits in range(total):
            A = binary_to_tournament(bits, n)
            scores = A.sum(axis=1)
            if not all(s == (n - 1) // 2 for s in scores):
                continue

            H = count_ham_paths(A)
            S = signed_adjacency(A)
            S4_trace = float(np.trace(np.linalg.matrix_power(S, 4)))
            kurt = spectral_kurtosis(S)
            comm = commutant_dimension(S)

            regular_data.append({
                'bits': bits, 'H': H, 'tr_S4': S4_trace,
                'kurtosis': kurt, 'comm': comm
            })

        if not regular_data:
            print("  No regular tournaments (n even or too small)")
            continue

        # Group by H
        by_H = defaultdict(list)
        for d in regular_data:
            by_H[d['H']].append(d)

        print(f"  {'H':>5s} {'count':>6s} {'tr(S^4)':>10s} {'kurtosis':>10s} {'comm':>6s}")
        for H in sorted(by_H.keys()):
            group = by_H[H]
            rep = group[0]
            print(f"  {H:5d} {len(group):6d} {rep['tr_S4']:10.1f} "
                  f"{rep['kurtosis']:10.4f} {rep['comm']:6d}")

        # Check: min tr(S^4) <=> max H?
        all_data = regular_data
        min_tr4 = min(d['tr_S4'] for d in all_data)
        max_H = max(d['H'] for d in all_data)
        min_tr4_H = [d['H'] for d in all_data if d['tr_S4'] == min_tr4]
        max_H_tr4 = [d['tr_S4'] for d in all_data if d['H'] == max_H]

        print(f"\n  Min tr(S^4) = {min_tr4}: H = {sorted(set(min_tr4_H))}")
        print(f"  Max H = {max_H}: tr(S^4) = {sorted(set(max_H_tr4))}")
        print(f"  SPECTRAL FLATNESS PRINCIPLE: "
              f"{'HOLDS' if max_H == max(min_tr4_H) else 'FAILS'}")


def investigate_comm_formula():
    """
    Verify the commutant dimension formula for Paley tournaments.
    dim(Comm(S_{Paley T_p})) = (p^2 - 2p + 3) / 2
    """
    print("\n" + "=" * 70)
    print("COMMUTANT DIMENSION FORMULA FOR PALEY")
    print("=" * 70)

    for p in [3, 5, 7, 11, 13]:
        if not all(p % i != 0 for i in range(2, int(p**0.5)+1)):
            continue
        if p % 4 != 3:
            print(f"\n  p={p}: p != 3 mod 4, Paley tournament not standard DRT")
            continue

        qr = set()
        for k in range(1, p):
            qr.add((k * k) % p)

        A = np.zeros((p, p), dtype=int)
        for i in range(p):
            for j in range(p):
                if i != j and (j - i) % p in qr:
                    A[i][j] = 1

        S = signed_adjacency(A)
        comm = commutant_dimension(S)
        formula = (p * p - 2 * p + 3) // 2

        H = count_ham_paths(A) if p <= 13 else None

        print(f"\n  p={p}: comm = {comm}, formula = {formula}, "
              f"{'MATCH' if comm == formula else 'MISMATCH'}")
        if H is not None:
            print(f"    H = {H}")

    # For n=5 (p=5, p=1 mod 4): use the unique regular tournament instead
    print(f"\n  n=5 (regular, NOT p=3 mod 4):")
    m = 10
    for bits in range(1 << m):
        A = binary_to_tournament(bits, 5)
        if all(A.sum(axis=1) == 2):
            S = signed_adjacency(A)
            comm = commutant_dimension(S)
            H = count_ham_paths(A)
            print(f"    Regular T_5: comm = {comm}, H = {H}")
            break


def investigate_comm_H_at_n5():
    """
    At n=5, there are only two spectra.
    Check the full comm vs H picture.
    """
    print("\n" + "=" * 70)
    print("FULL COMM vs H PICTURE AT n=5")
    print("=" * 70)

    n = 5
    m = n * (n - 1) // 2
    total = 1 << m

    by_pair = defaultdict(int)  # (H, comm) -> count

    for bits in range(total):
        A = binary_to_tournament(bits, n)
        H = count_ham_paths(A)
        S = signed_adjacency(A)
        comm = commutant_dimension(S)
        by_pair[(H, comm)] += 1

    print(f"\n  {'H':>5s} {'comm':>6s} {'count':>8s}")
    for (H, comm), count in sorted(by_pair.items()):
        print(f"  {H:5d} {comm:6d} {count:8d}")

    print(f"\n  KEY: Higher comm generally => higher H")


if __name__ == "__main__":
    test_commutant_maximality()
    test_spectral_flatness()
    investigate_comm_formula()
    investigate_comm_H_at_n5()
