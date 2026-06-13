#!/usr/bin/env python3
"""
exact_H_comparison.py -- Exact H(T) via Held-Karp DP for Hamiltonian paths

H(T) = number of Hamiltonian paths in tournament T.
Uses bitmask DP: O(2^n * n^2) time.

Computes H for Interval and Paley tournaments at p=7, 11, 13, 17, 19.

Author: kind-pasteur-2026-03-12-S59c
"""

import time


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def count_ham_paths(A, n):
    """Count all Hamiltonian paths in digraph A on n vertices using Held-Karp DP."""
    # dp[mask][v] = number of paths visiting exactly the vertices in mask, ending at v
    dp = [[0] * n for _ in range(1 << n)]

    # Base case: single-vertex paths
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


def count_ham_paths_sparse(A, n):
    """Memory-efficient Held-Karp for larger n using dict."""
    # For n > 20, use dict-based DP
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1

    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp:
                continue
            cnt = dp[key]
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + cnt

    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def main():
    print("=" * 70)
    print("EXACT H(T) COMPARISON: Interval vs Paley")
    print("=" * 70)

    results = {}

    for p in [5, 7, 11, 13]:
        m = (p - 1) // 2
        S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        S_int = list(range(1, m + 1))

        paley_exists = (p % 4 == 3)

        tournaments = []
        if paley_exists:
            tournaments.append(("Paley", S_qr))
        tournaments.append(("Interval", S_int))

        # Also test "reverse interval" (anti-interval)
        S_rint = list(range(m + 1, p))  # complement of Interval
        # Note: S_rint is NOT a valid tournament connection set for p odd
        # since it doesn't have exactly m elements... let me check.
        # S_int has m elements, S_rint has m elements (from m+1 to p-1).
        # But S_rint = {m+1, ..., 2m} = {m+1, ..., p-1}
        # These are the "far" connections.
        # S_rint union -S_int should cover Z_p* ... actually S_rint = p - S_int (mod p)
        # = {p-1, p-2, ..., p-m} = {2m, 2m-1, ..., m+1}. So S_rint = {m+1,...,2m}.
        # And S_int = {1,...,m}. S_int intersect S_rint = empty. S_int union S_rint = Z_p*.
        # So S_rint is the "anti-Interval" = complement tournament.
        # For the complement T^c of a tournament T: H(T^c) = H(T) (by path reversal).
        # Wait, that's not right. The complement reverses ALL arcs, so
        # a Hamiltonian path in T becomes a Hamiltonian path in T^c (same path, reversed).
        # So H(T^c) = H(T). So H(Interval) = H(anti-Interval).

        print(f"\n{'='*70}")
        print(f"p={p}, m={m}")
        print(f"{'='*70}")

        for name, S in tournaments:
            A = build_adj(p, S)

            t0 = time.time()
            if p <= 15:
                H = count_ham_paths(A, p)
            else:
                H = count_ham_paths_sparse(A, p)
            t1 = time.time()

            print(f"  {name} (S={S}): H = {H} ({t1-t0:.2f}s)")
            results[(p, name)] = H

            # Also compute |Aut(T)|
            # For circulant: Aut contains Z_p rotations at minimum
            aut_size = p  # at least p (circulant)
            # Check if i -> -i (mod p) is an automorphism
            # This maps arc (i,j) with j-i in S to arc (-i,-j) with -j-(-i)=i-j=-(j-i)
            # So it's an automorphism iff -S = S mod p, i.e., S is symmetric.
            # For Interval: -S = {p-1, p-2, ..., p-m} = S_rint. So -S != S unless m=0.
            # For Paley: QR is closed under negation iff -1 is a QR, iff p = 1 mod 4.
            # But Paley exists only for p = 3 mod 4 where -1 is NOT a QR.
            # So for Paley p=3 mod 4: -QR = NQR. The map i->-i sends QR arcs to NQR arcs,
            # which reverses the tournament. So it's an ANTI-automorphism, not automorphism.
            # |Aut(Paley)| = p for p prime (Z_p rotations only? Actually larger for some p).

            print(f"    |Aut| >= {aut_size}")
            print(f"    H/|Aut| = {H // aut_size} (remainder {H % aut_size})")

    # ====== SUMMARY ======
    print(f"\n{'='*70}")
    print("SUMMARY: Phase Transition")
    print("=" * 70)

    for p in [5, 7, 11, 13]:
        paley_exists = (p % 4 == 3)
        print(f"\n  p={p}:")
        if paley_exists and (p, "Paley") in results:
            h_pal = results[(p, "Paley")]
            h_int = results[(p, "Interval")]
            winner = "Paley" if h_pal > h_int else "Interval" if h_int > h_pal else "TIE"
            print(f"    Paley:    H = {h_pal}")
            print(f"    Interval: H = {h_int}")
            print(f"    Winner: {winner}")
            print(f"    Ratio Int/Pal: {h_int/h_pal:.6f}")
            print(f"    Delta: {h_int - h_pal}")
        elif (p, "Interval") in results:
            h_int = results[(p, "Interval")]
            print(f"    Interval: H = {h_int}")

    # ====== ALPHA DECOMPOSITION ======
    print(f"\n{'='*70}")
    print("OCF DECOMPOSITION: H = 1 + 2*alpha_1 + 4*alpha_2 + ...")
    print("=" * 70)

    for p in [7, 11]:
        m = (p - 1) // 2
        S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        S_int = list(range(1, m + 1))

        paley_exists = (p % 4 == 3)
        if not paley_exists:
            continue

        print(f"\n  p={p}:")

        from itertools import combinations as combs

        for name, S in [("Paley", S_qr), ("Interval", S_int)]:
            A = build_adj(p, S)

            # Get all odd-cycle vertex sets
            all_cycles = []
            cycle_by_k = {}
            for k in range(3, p + 1, 2):
                sets_k = []
                for subset in combs(range(p), k):
                    verts = list(subset)
                    if k == 3:
                        a, b, c = verts
                        hc = (A[a][b]*A[b][c]*A[c][a] + A[a][c]*A[c][b]*A[b][a]) > 0
                    else:
                        hc = False
                        dp2 = set()
                        dp2.add((1 << 0, 0))
                        for mask in range(1, 1 << k):
                            if not (mask & 1):
                                continue
                            for v in range(k):
                                if not (mask & (1 << v)):
                                    continue
                                if (mask, v) not in dp2:
                                    continue
                                for w in range(k):
                                    if mask & (1 << w):
                                        continue
                                    if A[verts[v]][verts[w]]:
                                        dp2.add((mask | (1 << w), w))
                        full = (1 << k) - 1
                        for v in range(1, k):
                            if (full, v) in dp2 and A[verts[v]][verts[0]]:
                                hc = True
                                break
                    if hc:
                        sets_k.append(frozenset(subset))
                cycle_by_k[k] = sets_k
                all_cycles.extend(sets_k)

            n = len(all_cycles)
            print(f"\n    {name}: {n} total odd-cycle vertex sets")
            for k in sorted(cycle_by_k.keys()):
                print(f"      c_{k} = {len(cycle_by_k[k])}")

            # Compute alpha_1, alpha_2 from Omega
            alpha_1 = n

            # alpha_2: count disjoint pairs
            alpha_2 = 0
            for i in range(n):
                for j in range(i + 1, n):
                    if not (all_cycles[i] & all_cycles[j]):
                        alpha_2 += 1

            # alpha_2 by cycle type
            print(f"      alpha_1 = {alpha_1}")
            print(f"      alpha_2 = {alpha_2}")

            # Break down alpha_2 by cycle type
            for k1 in sorted(cycle_by_k.keys()):
                for k2 in sorted(cycle_by_k.keys()):
                    if k1 > k2:
                        continue
                    if 2*max(k1,k2) > p:
                        continue  # can't be disjoint
                    count = 0
                    for s1 in cycle_by_k[k1]:
                        start = 0 if k1 != k2 else cycle_by_k[k2].index(s1) + 1
                        for idx in range(start if k1 == k2 else 0, len(cycle_by_k[k2])):
                            s2 = cycle_by_k[k2][idx]
                            if k1 == k2 and idx <= cycle_by_k[k1].index(s1):
                                continue
                            if not (s1 & s2):
                                count += 1
                    if count > 0 or k1 + k2 <= p:
                        print(f"        ({k1},{k2})-disjoint: {count}")

            H = results[(p, name)]
            remainder = H - 1 - 2*alpha_1 - 4*alpha_2
            print(f"      H = {H}")
            print(f"      1 + 2*alpha_1 + 4*alpha_2 = {1 + 2*alpha_1 + 4*alpha_2}")
            print(f"      Remainder (8*alpha_3 + ...) = {remainder}")


if __name__ == '__main__':
    main()
