"""
satake_ndrt_h.py — Compute H and Omega dims for Satake's cyclotomic NDR tournaments.

Satake (arXiv:2502.12090): For prime q ≡ 5 (mod 8) with q = s² + 4,
the cyclotomic tournament CT_q is nearly-doubly-regular (NDR).
Connection set: S = C_4 ∪ g·C_4 where C_4 = {x^4 mod q : x ∈ F*_q},
g = primitive root mod q.

Key question: Do these NDR tournaments MAXIMIZE H(T) among all q-vertex tournaments?
Our PALEY MAXIMIZER result: H(T_p) = max for p ≡ 3 (mod 4).
For q ≡ 5 (mod 8) (i.e. q ≡ 1 mod 4), Paley doesn't apply — these NDRTs are candidates.

Candidates:
  q=5  (s=1):  S = {1,2}              (5-cycle tournament — classic DRT, n≡1 mod 4)
  q=13 (s=3):  S = {1,2,3,5,6,9}
  q=29 (s=5):  ...

Method:
  H(T) via Hamiltonian path enumeration using diff-seqs (circulant approach).
  For vertex-transitive tournament CT_q: H(T) = q × |{Ham paths starting at 0}| / q = |A_{q-1}|
  Wait: actually H(T) = total Ham paths in T = |A_{q-1}| × q (by vertex-transitivity,
  paths starting at each vertex are equinumerous, each contributes |A_{q-1}|).
  Actually |A_{q-1}| for the diff-seq enumeration = # Ham paths starting at vertex 0.
  Total Ham paths H(T) = n × |A_{q-1}|.

  CORRECTED: For circulant tournament, vertex-transitivity gives H(T) = n × (paths from v=0).
  But paths from v=0 = |A_{q-1}| (the diff-seq count starting from 0) since diff-seqs
  enumerate all ways to choose steps of length in S with distinct partial sums mod q.
  So H(CT_q) = q × |A_{q-1}|... wait, need to check: H(T) = total Hamiltonian paths
  = number of (start, sequence) pairs = sum_v (paths from v) = q × (paths from 0) by symmetry.

  From OEIS: H(T_7) = 189, and the Paley T_7 has |A_6| from vertex 0.
  Let's check: 189 / 7 = 27 paths from each vertex. So H = 7 × 27 = 189 ✓.
  And |A_6| in our circulant_homology = 3 (from test: path_count_sequence[-1]=3?).
  Wait, from test_circulant_homology.py: path_count_sequence[-1] > 0 but let me check actual value.

Author: kind-pasteur-2026-03-12-S56
"""
import sys
import time
sys.path.insert(0, '04-computation')


def primitive_root_mod_q(q):
    """Find primitive root mod q (prime)."""
    from sympy import primitive_root
    return primitive_root(q)


def quartic_cosets(q):
    """
    Compute the four cosets of quartic residues mod q.
    Requires q prime with 4 | (q-1).
    Returns list of 4 sorted cosets.
    """
    from sympy import isprime
    assert isprime(q) and (q - 1) % 4 == 0
    g = primitive_root_mod_q(q)
    n_coset = (q - 1) // 4
    C4 = sorted(pow(g, 4 * k, q) for k in range(n_coset))
    cosets = []
    for i in range(4):
        fi = pow(g, i, q)
        cosets.append(sorted((fi * x) % q for x in C4))
    return cosets, g


def satake_connection_set(q):
    """
    Satake's NDR connection set: S = C_4 ∪ g·C_4 (first two cosets).
    Returns sorted list.
    """
    cosets, _ = quartic_cosets(q)
    S = sorted(set(cosets[0] + cosets[1]))
    neg_S = sorted((q - s) % q for s in S)
    assert set(S) & set(neg_S) == set(), f"S ∩ (-S) ≠ ∅ for q={q}!"
    assert set(S) | set(neg_S) == set(range(1, q))
    return S


def ham_paths_from_0(n, S_set):
    """
    Count Hamiltonian paths starting from vertex 0 in circulant tournament C_n^S.
    Uses dynamic programming on (current_vertex, visited_set).

    Returns: number of Hamiltonian paths starting at 0.
    """
    # DP: state = (last_vertex, frozenset of visited)
    # Initialize: at vertex 0, visited = {0}
    dp = {(0, frozenset([0])): 1}
    S_list = sorted(S_set)

    for step in range(n - 1):
        new_dp = {}
        for (v, vis), cnt in dp.items():
            for s in S_list:
                w = (v + s) % n
                if w not in vis:
                    key = (w, vis | frozenset([w]))
                    new_dp[key] = new_dp.get(key, 0) + cnt
        dp = new_dp
        if not dp:
            break

    return sum(dp.values())


def ham_count_circulant(n, S):
    """H(C_n^S) using vertex-transitivity: H = n × paths_from_0."""
    paths = ham_paths_from_0(n, set(S))
    return n * paths, paths


def count_3cycles(n, S_set):
    """Count directed 3-cycles in circulant tournament C_n^S."""
    count = 0
    for a in range(n):
        for b_diff in S_set:
            b = (a + b_diff) % n
            for c_diff in S_set:
                c = (b + c_diff) % n
                if c != a and (a - c) % n in S_set:
                    count += 1
    return count // 3


def all_circulant_H(n):
    """
    Enumerate ALL valid circulant tournament connection sets on n vertices
    and compute H for each. Returns sorted list of (H, S).

    A valid S: S ∩ (-S) = ∅, S ∪ (-S) = Z/nZ \ {0}.
    """
    from itertools import product
    k = (n - 1) // 2
    # Pair structure: (s, n-s) for s=1,...,k-1, and s=k if n=2k+1 (odd n)
    pairs = [(s, n - s) for s in range(1, k + 1) if s < n - s]
    if n % 2 == 1 and k == (n - 1) // 2:
        # No self-paired element for odd n
        pass

    results = []
    for choices in product([0, 1], repeat=len(pairs)):
        S = sorted(pairs[i][c] for i, c in enumerate(choices))
        H, _ = ham_count_circulant(n, S)
        results.append((H, S))

    results.sort(reverse=True)
    return results


def main():
    print("=" * 70)
    print("SATAKE NDR TOURNAMENT INVESTIGATION (S56)")
    print("=" * 70)

    # -----------------------------------------------------------------------
    # Section 1: Satake connection sets
    # -----------------------------------------------------------------------
    print("\n--- Satake connection sets ---")
    for s_val, q in [(1, 5), (3, 13), (5, 29)]:
        S = satake_connection_set(q)
        c3 = count_3cycles(q, set(S))
        print(f"q={q} (s={s_val}): S={S}, c_3={c3}")

    # -----------------------------------------------------------------------
    # Section 2: H for Satake NDRTs at small q
    # -----------------------------------------------------------------------
    print("\n--- H(CT_q) via Hamiltonian path enumeration ---")

    for s_val, q in [(1, 5), (3, 13)]:
        S = satake_connection_set(q)
        t0 = time.time()
        H, paths0 = ham_count_circulant(q, S)
        elapsed = time.time() - t0
        c3 = count_3cycles(q, set(S))
        print(f"q={q}: H(CT_{q}) = {H} (paths from 0: {paths0}, c_3={c3}) [{elapsed:.2f}s]")

    # -----------------------------------------------------------------------
    # Section 3: Compare to ALL circulant tournaments at n=5 and n=13
    # -----------------------------------------------------------------------
    print("\n--- All circulant H at n=5 (exhaustive) ---")
    results5 = all_circulant_H(5)
    for H, S in results5[:5]:
        print(f"  H={H}, S={S}")
    S_satake5 = satake_connection_set(5)
    print(f"Satake CT_5 S={S_satake5} rank: #{next(i+1 for i,(H,S) in enumerate(results5) if S==S_satake5)}")

    print("\n--- All circulant H at n=13 (exhaustive over circulant) ---")
    t0 = time.time()
    results13 = all_circulant_H(13)
    elapsed = time.time() - t0
    print(f"[{elapsed:.1f}s, {len(results13)} distinct circulant tournaments]")
    # Show top 10
    seen_H = {}
    for H, S in results13:
        if H not in seen_H:
            seen_H[H] = S
    print("Top H values and representative S:")
    for H in sorted(seen_H.keys(), reverse=True)[:10]:
        S = seen_H[H]
        c3 = count_3cycles(13, set(S))
        print(f"  H={H}, c_3={c3}, S={S}")

    S_satake13 = satake_connection_set(13)
    rank13 = next(i + 1 for i, (H, S) in enumerate(results13) if S == S_satake13)
    H_satake13 = next(H for H, S in results13 if S == S_satake13)
    H_max13 = results13[0][0]
    print(f"\nSatake CT_13: H={H_satake13}, rank #{rank13}/{len(results13)}")
    print(f"Max H among circulants at n=13: {H_max13}")
    print(f"Satake is {'THE MAXIMIZER' if H_satake13 == H_max13 else f'NOT max (gap = {H_max13 - H_satake13})'}")

    # -----------------------------------------------------------------------
    # Section 4: Compare coset choices at n=13
    # -----------------------------------------------------------------------
    print("\n--- All valid coset-pair tournaments at n=13 ---")
    from itertools import combinations
    cosets13, g13 = quartic_cosets(13)
    print(f"Cosets of C_4 mod 13 (g={g13}):")
    for i, c in enumerate(cosets13):
        print(f"  coset {i}: {c}")

    for (i, j) in combinations(range(4), 2):
        S_test = sorted(set(cosets13[i] + cosets13[j]))
        neg_S = sorted((13 - s) % 13 for s in S_test)
        if set(S_test) & set(neg_S) == set():
            H_test, paths0_test = ham_count_circulant(13, S_test)
            c3_test = count_3cycles(13, set(S_test))
            is_satake = sorted(S_test) == sorted(S_satake13)
            print(f"  cosets {i}+{j}: S={S_test}, H={H_test}, c_3={c3_test} {'<- Satake' if is_satake else ''}")

    # -----------------------------------------------------------------------
    # Section 5: Check isomorphism classes at n=13
    # -----------------------------------------------------------------------
    print("\n--- Unique H values at n=13 ---")
    print(f"Satake H(CT_13) = {H_satake13}")
    print(f"Maximum circulant H at n=13 = {H_max13}")

    # -----------------------------------------------------------------------
    # Section 6: n=9 context (q=9 = 3^2, not prime, but ≡ 1 mod 4)
    # -----------------------------------------------------------------------
    print("\n--- n=9 context (Satake needs prime q, but let's check all circulants) ---")
    t0 = time.time()
    results9 = all_circulant_H(9)
    print(f"All circulant H at n=9 (top 5): {[H for H,_ in results9[:5]]}")
    print(f"Max circulant H at n=9: {results9[0][0]}, S={results9[0][1]}")


if __name__ == '__main__':
    main()
