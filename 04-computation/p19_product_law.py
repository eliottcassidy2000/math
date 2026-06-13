#!/usr/bin/env python3
"""
p19_product_law.py -- Optimized product law verification at p=19

Uses symmetry H(sigma) = H(-sigma) to halve computation.
Also tries p=23 if time permits.

Key question: does sign(h_hat[{a,b}]) = chi(a*b) at p=19?
At p=19, chi(3)=-1, so 3-resonant pairs have sign(D^4)*chi(ab) = -1 (wrong sign).
But the full h_hat includes all D^{2n} contributions weighted by OCF.
Does the nonlinear OCF structure restore the product law?

Author: kind-pasteur-2026-03-12-S60
"""

import math
import time


def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1


def count_ham_paths_fast(adj_list, n):
    """Count Hamiltonian paths using bitmask DP with adjacency list.
    O(2^n * n) time, O(2^n * n) space."""
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1

    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            cnt = dp[mask][v]
            if cnt == 0:
                continue
            for w in adj_list[v]:
                if not (mask & (1 << w)):
                    dp[mask | (1 << w)][w] += cnt

    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def classify_resonance(a, b, p):
    """Classify the resonance type of pair (a,b)."""
    resonances = []
    for k in range(1, p):
        q = 2*k - 1
        if q >= p:
            break
        if (q*a - b) % p == 0:
            resonances.append((q, f"{q}a=b"))
        if (q*a + b) % p == 0:
            resonances.append((q, f"{q}a=-b"))
        if (a - q*b) % p == 0 and q != 1:
            resonances.append((q, f"a={q}b"))
        if (a + q*b) % p == 0 and q != 1:
            resonances.append((q, f"a=-{q}b"))
    return resonances


def product_law_test(p):
    """Test sign(h_hat[{a,b}]) = chi(a*b) at prime p."""
    m = (p - 1) // 2
    pairs = [(s, p - s) for s in range(1, m + 1)]

    print(f"\n{'='*70}")
    print(f"PRODUCT LAW TEST at p={p}, m={m}")
    print(f"{'='*70}")
    print(f"  chi(-1) = {legendre(-1, p):+d}  (p mod 4 = {p % 4})")
    print(f"  chi(3)  = {legendre(3, p):+d}")
    print(f"  chi(5)  = {legendre(5, p):+d}")
    print(f"  chi(7)  = {legendre(7, p):+d}")
    print(f"  chi(9)  = {legendre(9, p):+d}")

    # Use symmetry: H(sigma) = H(-sigma), so only compute half
    # For degree-2 Walsh: h_hat[{a,b}] = (2/2^m) * sum_{sigma with sigma_0=+1} H(sigma)*sigma_a*sigma_b
    half = 1 << (m - 1)
    H_plus = {}  # H values for sigma with sigma_0 = +1

    print(f"\n  Computing H for {half} orientations (symmetry reduction from {1<<m})...")
    t0 = time.time()

    for bits_rest in range(half):
        # sigma_0 = +1 always; bits_rest encodes sigma_1, ..., sigma_{m-1}
        sigma = [1]  # sigma_0 = +1
        for i in range(m - 1):
            sigma.append(1 if bits_rest & (1 << i) else -1)

        # Build adjacency list
        S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))
        adj_list = [[] for _ in range(p)]
        for v in range(p):
            for s in S:
                adj_list[v].append((v + s) % p)

        H_plus[bits_rest] = count_ham_paths_fast(adj_list, p)

        if bits_rest > 0 and bits_rest % 32 == 0:
            elapsed = time.time() - t0
            pct = bits_rest / half * 100
            eta = elapsed / bits_rest * (half - bits_rest)
            print(f"    {pct:.0f}% ({bits_rest}/{half}), elapsed={elapsed:.0f}s, ETA={eta:.0f}s")

    t1 = time.time()
    print(f"  Done in {t1-t0:.1f}s")

    # Compute degree-2 Walsh using symmetry
    # h_hat[{a,b}] = (1/2^m) * sum_sigma H(sigma)*sigma_a*sigma_b
    # = (1/2^m) * 2 * sum_{sigma: sigma_0=+1} H(sigma)*sigma_a*sigma_b
    # = (1/2^{m-1}) * sum_{bits_rest} H_plus[bits_rest] * sigma_a * sigma_b

    print(f"\n  Degree-2 Walsh coefficients:")
    match_count = 0
    total_pairs = 0

    for a_idx in range(m):
        for b_idx in range(a_idx + 1, m):
            a, b = a_idx + 1, b_idx + 1
            chi_ab = legendre(a * b, p)

            h_hat = 0
            for bits_rest in range(half):
                # Reconstruct sigma_a and sigma_b
                if a_idx == 0:
                    sa = 1  # sigma_0 always +1
                else:
                    sa = 1 if bits_rest & (1 << (a_idx - 1)) else -1
                if b_idx == 0:
                    sb = 1  # sigma_0 always +1
                else:
                    sb = 1 if bits_rest & (1 << (b_idx - 1)) else -1

                h_hat += H_plus[bits_rest] * sa * sb

            h_hat = h_hat / (1 << (m - 1))  # factor of 2 from symmetry included

            sign_h = 1 if h_hat > 0 else (-1 if h_hat < 0 else 0)
            match = (sign_h == chi_ab)

            res = classify_resonance(a, b, p)
            min_q = min(q for q, t in res) if res else 'inf'
            first_type = res[0][1] if res else 'none'

            total_pairs += 1
            if match:
                match_count += 1

            print(f"    ({a:>2},{b:>2}): h_hat={h_hat:>14.2f}, sign={sign_h:+d}, "
                  f"chi={chi_ab:+d}, {'OK' if match else '**FAIL**'}, "
                  f"q={min_q}, type={first_type}")

    print(f"\n  PRODUCT LAW: sign=chi in {match_count}/{total_pairs} pairs")
    if match_count == total_pairs:
        print(f"  >>> HOLDS at p={p}! <<<")
    else:
        print(f"  >>> FAILS at p={p} ({total_pairs - match_count} mismatches) <<<")

    # Additional analysis: degree-1 Walsh (should be zero for symmetric pairs)
    print(f"\n  Degree-1 Walsh (should reflect chi structure):")
    for a_idx in range(m):
        a = a_idx + 1
        h1 = 0
        for bits_rest in range(half):
            if a_idx == 0:
                sa = 1
            else:
                sa = 1 if bits_rest & (1 << (a_idx - 1)) else -1
            h1 += H_plus[bits_rest] * sa
        h1 = h1 / (1 << (m - 1))
        chi_a = legendre(a, p)
        sign_h1 = 1 if h1 > 0.01 else (-1 if h1 < -0.01 else 0)
        print(f"    a={a}: h1={h1:>12.2f}, sign={sign_h1:+d}, chi(a)={chi_a:+d}")


def overlap_resonance_connection(p):
    """Investigate connection between overlap weights and resonance structure.

    Key question: does the cycle overlap matrix in Omega(T) reflect the
    multiplicative resonance structure of pairs?

    For Paley T_p:
    - 3-cycle vertex sets correspond to zero-sum triples in QR
    - Two 3-cycles conflict iff they share a vertex
    - The co-occurrence matrix should reflect additive structure of QR

    Connection to resonance:
    - Resonance qa=±b means chords a,b have a multiplicative relationship
    - Overlap means cycles share vertices (additive relationship on Z_p)
    - The WALSH transform converts between additive and multiplicative views
    """
    m = (p - 1) // 2
    QR = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

    print(f"\n{'='*70}")
    print(f"OVERLAP-RESONANCE CONNECTION at p={p}")
    print(f"{'='*70}")

    # Build Paley adjacency
    QR_set = set(QR)
    A = [[0]*p for _ in range(p)]
    for v in range(p):
        for s in QR:
            A[v][(v + s) % p] = 1

    # Enumerate 3-cycles
    from itertools import combinations
    c3_sets = []
    for a, b, c in combinations(range(p), 3):
        if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
            c3_sets.append(frozenset([a, b, c]))
    c3_sets = list(set(c3_sets))
    n3 = len(c3_sets)

    print(f"  c3 vertex sets: {n3}")

    # For each pair of chord indices (i,j), count how many pairs of 3-cycles
    # "use" chord i and chord j. A 3-cycle {a,b,c} uses chord s if
    # (b-a)%p = s or (c-b)%p = s or (a-c)%p = s (i.e., one of its gaps is s).

    # Build gap decomposition of each 3-cycle
    chord_usage = []  # for each 3-cycle, which chord values (from 1..m) appear as gaps
    for fs in c3_sets:
        verts = sorted(fs)
        gaps = set()
        for i in range(3):
            for j in range(3):
                if i != j:
                    d = (verts[j] - verts[i]) % p
                    if d <= m:
                        gaps.add(d)
                    else:
                        gaps.add(p - d)  # map to chord index
        chord_usage.append(gaps)

    # Now for each pair of chord indices (a,b), count 3-cycles using BOTH
    # This is the "chord co-participation" matrix
    print(f"\n  Chord co-participation matrix (# 3-cycles using both chords a and b):")
    copart = [[0]*m for _ in range(m)]
    for cu in chord_usage:
        cu_list = sorted(cu)
        for i in range(len(cu_list)):
            for j in range(i+1, len(cu_list)):
                ci, cj = cu_list[i] - 1, cu_list[j] - 1  # 0-indexed
                copart[ci][cj] += 1
                copart[cj][ci] += 1

    for a_idx in range(m):
        for b_idx in range(a_idx + 1, m):
            a, b = a_idx + 1, b_idx + 1
            chi_ab = legendre(a * b, p)
            res = classify_resonance(a, b, p)
            min_q = min(q for q, t in res) if res else 'inf'
            print(f"    ({a},{b}): copart={copart[a_idx][b_idx]:>3d}, "
                  f"chi={chi_ab:+d}, q={min_q}")

    # Co-occurrence of vertices at gap d: how many 3-cycles contain
    # a pair of vertices separated by gap d?
    print(f"\n  Vertex gap co-occurrence vs chord index:")
    gap_co = [0] * p
    for fs in c3_sets:
        verts = sorted(fs)
        for i in range(3):
            for j in range(i+1, 3):
                d = min((verts[j] - verts[i]) % p, (verts[i] - verts[j]) % p)
                gap_co[d] += 1
    for d in range(1, m + 1):
        print(f"    d={d:>2}: {gap_co[d]:>3d} co-occurring pairs, chi(d)={legendre(d,p):+d}")

    # KEY CONNECTION: Is the co-occurrence profile (as function of gap d)
    # related to the Walsh spectrum?
    # Fourier transform of gap_co should relate to eigenvalues of the
    # co-occurrence matrix.

    # For circulant structures, everything decomposes by characters of Z_p
    import cmath
    omega = cmath.exp(2j * cmath.pi / p)
    fourier_co = []
    for t in range(p):
        val = sum(gap_co[d] * omega**(t*d) for d in range(1, p))
        fourier_co.append(val)

    print(f"\n  Fourier of co-occurrence (|c_hat(t)|):")
    for t in range(1, m + 1):
        val = abs(fourier_co[t])
        chi_t = legendre(t, p)
        print(f"    t={t:>2}: |c_hat|={val:>8.2f}, chi(t)={chi_t:+d}")

    # Disjointness rate by chord pair type
    print(f"\n  Disjointness analysis by chord pair resonance:")
    # For each pair of 3-cycles, track which chord pairs they share
    # and whether they're disjoint

    for a_idx in range(m):
        for b_idx in range(a_idx + 1, m):
            a, b = a_idx + 1, b_idx + 1
            # Count 3-cycle pairs where one uses chord a and the other uses chord b
            uses_a = [i for i, cu in enumerate(chord_usage) if a in cu]
            uses_b = [i for i, cu in enumerate(chord_usage) if b in cu]

            disj = 0
            total = 0
            for i in uses_a:
                for j in uses_b:
                    if i == j:
                        continue
                    total += 1
                    if not (c3_sets[i] & c3_sets[j]):
                        disj += 1

            if total > 0:
                chi_ab = legendre(a * b, p)
                res = classify_resonance(a, b, p)
                min_q = min(q for q, t in res) if res else 'inf'
                print(f"    ({a},{b}): disj_rate={disj/total:.4f} ({disj}/{total}), "
                      f"chi={chi_ab:+d}, q={min_q}")


def main():
    # First: overlap-resonance connection at small primes
    for p in [7, 11]:
        overlap_resonance_connection(p)

    # Then: product law at p=19
    product_law_test(19)


if __name__ == '__main__':
    main()
