"""
c5_spectral_formula.py — Derive c_5 = f(eigenvalues) for circulant tournaments.

KEY DISCOVERY (from full_ocf_spectral.py):
  At p=7: c_3 = 14 (constant), but c_5 varies:
    Paley:     c_5 = 42, c_7 = 24   (42 directed 5-cycles on 21 vsets = 2.0/vset)
    Non-Paley: c_5 = 28, c_7 = 17   (28 directed 5-cycles on 21 vsets = 1.33/vset)

  The extra directed cycles in Paley come from REGULAR 5-vertex subtournaments
  (score 2,2,2,2,2) which have 2 directed 5-cycles each, vs near-regular
  subtournaments (score 1,2,2,2,3) which have 1 directed 5-cycle.

  For circulant tournaments, c_5 can be computed from eigenvalues via:
    tr(A^5) = sum lambda_k^5 = 5*c_5 + correction_terms
  where correction terms count non-simple closed 5-walks.

  THIS SCRIPT:
  1. Compute tr(A^5) and c_5 independently, find the correction
  2. Show c_5 is a polynomial in eigenvalue power sums
  3. Derive the FULL spectral formula H = f(c_3, c_5, c_7)

Author: kind-pasteur-2026-03-12-S56c
"""

import sys
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


def circulant_eigenvalues(n, S):
    omega = cmath.exp(2j * cmath.pi / n)
    return [sum(omega ** (k * s) for s in S) for k in range(n)]


def count_directed_k_cycles(p, S, k):
    """Count directed k-cycles in circulant tournament on Z_p."""
    S_set = set(S)
    adj = [[False] * p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in S_set:
                adj[i][j] = True

    count = 0
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
                count += 1
    return count


def main():
    print("=" * 70)
    print("c_5 SPECTRAL FORMULA FOR CIRCULANT TOURNAMENTS")
    print("=" * 70)

    # ================================================================
    # SECTION 1: tr(A^k) vs c_k relationship
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"SECTION 1: TRACE vs CYCLE COUNT DECOMPOSITION")
    print(f"{'=' * 60}")
    print("""
    tr(A^k) counts ALL closed k-walks, not just simple k-cycles.
    For tournaments: no backward edges (a_{ij}*a_{ji} = 0).
    But closed walks CAN revisit vertices through different edges.

    For k=3: tr(A^3) = 3*c_3 (all closed 3-walks are simple cycles)
    For k=4: tr(A^4) = ??? (includes non-simple 4-walks)
    For k=5: tr(A^5) = 5*c_5 + correction_terms
    """)

    for p in [5, 7, 11]:
        print(f"\n--- p = {p} ---")
        all_S = all_circulant_tournaments(p)
        m = (p - 1) // 2

        H_seen = set()
        for S in all_S:
            H = ham_count_dp(p, S)
            if H in H_seen:
                continue
            H_seen.add(H)

            eigs = circulant_eigenvalues(p, S)

            # Traces
            traces = {}
            for k in range(1, p + 1):
                traces[k] = sum(e ** k for e in eigs).real

            # Direct cycle counts
            c3 = count_directed_k_cycles(p, S, 3)
            c5 = count_directed_k_cycles(p, S, 5) if p >= 5 else 0
            c7 = count_directed_k_cycles(p, S, 7) if p >= 7 else 0

            print(f"  S={list(S)}, H={H}")
            print(f"    tr(A^3) = {traces[3]:.0f}, 3*c_3 = {3 * c3}, diff = {traces[3] - 3 * c3:.0f}")
            print(f"    tr(A^4) = {traces[4]:.0f}")
            print(f"    tr(A^5) = {traces[5]:.0f}, 5*c_5 = {5 * c5}, diff = {traces[5] - 5 * c5:.0f}")
            if p >= 7:
                print(f"    tr(A^6) = {traces[6]:.0f}")
                print(f"    tr(A^7) = {traces[7]:.0f}, 7*c_7 = {7 * c7}, diff = {traces[7] - 7 * c7:.0f}")

            # For k=4 in tournament:
            # tr(A^4) = sum_i (A^4)_{ii} = sum_{i,j,k,l} a_{ij}*a_{jk}*a_{kl}*a_{li}
            # Non-simple 4-walks: can visit only 3 vertices (e.g., i->j->k->j is impossible
            # since a_{kj}*a_{jk} can't both be 1; but i->j->k->i->j has 5 steps not 4)
            # Actually for k=4: all 4-walks must visit 3 or 4 distinct vertices
            # 4 distinct: simple directed 4-cycle (but tournaments have no even cycles? NO!)
            # Wait, tournaments CAN have even-length directed cycles!
            # e.g., i->j, j->k, k->l, l->i is a valid 4-cycle

            # So tr(A^4) = 4*c_4 + correction for walks on 3 vertices
            # 3-vertex walks of length 4: i->j->k->i->j (but this has 4 steps: i,j,k,i,j)
            # Wait: a 4-step closed walk starting at i visits i,j,k,l,i (5 positions)
            # If only 3 distinct vertices {i,j,k}: possible walk i->j->k->i->j... no, must return to i
            # i->j->k->i: this is 3 steps, not 4
            # i->j->i->j: impossible since a_{ji}*a_{ij} = 0
            # So 3-vertex closed 4-walks are IMPOSSIBLE in tournaments
            # Therefore tr(A^4) = 4*c_4
            c4 = count_directed_k_cycles(p, S, 4) if p >= 4 else 0
            if c4 > 0 or True:
                print(f"    4*c_4 = {4 * c4}, tr(A^4) = {traces[4]:.0f}, diff = {traces[4] - 4 * c4:.0f}")

            # For k=5:
            # 5-step closed walk on tournament. Can visit 3, 4, or 5 distinct vertices.
            # 3 vertices: i->j->k->i->j->k->i but that's 6 steps...
            # Actually, 5 steps: i->j->k->i->j->i... need a_{ji}=1 and a_{ij}=1, impossible
            # 4 vertices: i->j->k->l->j->i? Need a_{lj}=1 and a_{ji}=1 OK
            # Wait, the walk is positions 0,1,2,3,4,0 with 5 edges
            # If 4 distinct: one vertex visited twice, say vertex j appears at positions 1 and 3
            # Walk: i->j->k->j->l->i. Needs a_{kj}=1 and a_{jl}=1. Possible!
            # These ARE non-simple closed walks.

            # Actually, let me count non-simple 5-walks directly
            if p <= 7:
                S_set = set(S)
                adj = [[False] * p for _ in range(p)]
                for i in range(p):
                    for j in range(p):
                        if i != j and (j - i) % p in S_set:
                            adj[i][j] = True

                # Count all closed 5-walks
                total_walks = 0
                simple_walks = 0
                for i0 in range(p):
                    for i1 in range(p):
                        if not adj[i0][i1]: continue
                        for i2 in range(p):
                            if not adj[i1][i2]: continue
                            for i3 in range(p):
                                if not adj[i2][i3]: continue
                                for i4 in range(p):
                                    if not adj[i3][i4]: continue
                                    if not adj[i4][i0]: continue
                                    total_walks += 1
                                    if len({i0,i1,i2,i3,i4}) == 5:
                                        simple_walks += 1

                non_simple = total_walks - simple_walks
                print(f"    Total 5-walks = {total_walks}, simple = {simple_walks}, "
                      f"non-simple = {non_simple}")
                print(f"    tr(A^5) should = total_walks = {total_walks}: "
                      f"match = {abs(traces[5] - total_walks) < 0.5}")
                print(f"    5*c_5 = {5 * c5}, simple_walks = {simple_walks}: "
                      f"match = {simple_walks == 5 * c5}")

    # ================================================================
    # SECTION 2: Non-simple walk correction
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"SECTION 2: NON-SIMPLE WALK CORRECTION FOR tr(A^5)")
    print(f"{'=' * 60}")

    # At p=7:
    # tr(A^5) = 5*c_5 + non_simple_5_walks
    # The non-simple part depends on c_3 and the graph structure
    # For a tournament: non-simple closed 5-walk visits 4 distinct vertices
    # with one vertex visited twice
    # Pattern: i->j->k->j->l->i (or permutations)
    # This requires: j is the repeated vertex, a_{kj}=1 (k->j) and a_{jl}=1 (j->l)
    # So vertex j has "in-neighbor" k and "out-neighbor" l within the 4-set {i,j,k,l}

    # Count: for each 4-subset and repeated vertex, count valid walks
    p = 7
    print(f"\n  Detailed walk decomposition at p = {p}:")
    for S in all_circulant_tournaments(p):
        H = ham_count_dp(p, S)
        S_set = set(S)
        adj = [[False] * p for _ in range(p)]
        for i in range(p):
            for j in range(p):
                if i != j and (j - i) % p in S_set:
                    adj[i][j] = True

        c5 = count_directed_k_cycles(p, S, 5)
        eigs = circulant_eigenvalues(p, S)
        tr5 = sum(e ** 5 for e in eigs).real

        # Non-simple 5-walks: count walks on 4 distinct vertices
        # with one repeated
        ns_walks = 0
        for i0 in range(p):
            for i1 in range(p):
                if not adj[i0][i1]: continue
                for i2 in range(p):
                    if not adj[i1][i2]: continue
                    for i3 in range(p):
                        if not adj[i2][i3]: continue
                        for i4 in range(p):
                            if not adj[i3][i4]: continue
                            if not adj[i4][i0]: continue
                            n_distinct = len({i0,i1,i2,i3,i4})
                            if n_distinct < 5:
                                ns_walks += 1

        # Also check: walks on 3 distinct vertices
        ns3 = 0
        for i0 in range(p):
            for i1 in range(p):
                if not adj[i0][i1]: continue
                for i2 in range(p):
                    if not adj[i1][i2]: continue
                    for i3 in range(p):
                        if not adj[i2][i3]: continue
                        for i4 in range(p):
                            if not adj[i3][i4]: continue
                            if not adj[i4][i0]: continue
                            if len({i0,i1,i2,i3,i4}) <= 3:
                                ns3 += 1

        print(f"  S={list(S)}: H={H}, c5={c5}, tr5={tr5:.0f}, "
              f"ns_walks={ns_walks}, ns3={ns3}")
        break  # Just one per H value is enough to see the pattern

    # ================================================================
    # SECTION 3: The alpha decomposition — H via OCF coefficients
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"SECTION 3: OCF DECOMPOSITION H = sum alpha_k * 2^k")
    print(f"{'=' * 60}")

    # At p=5: alpha = [1, 7], I(Omega, 2) = 1 + 14 = 15 = H. VERIFIED.
    # The 7 cycles are: 5 three-cycles + 2 five-cycles.
    # All 7 cycles MUST share at least one vertex (since 5 vertices total,
    # 3 per 3-cycle, 5 per 5-cycle). So alpha_2 = 0.
    # I(Omega, 2) = 1 + 2*7 = 15. Beautiful!

    print(f"""
  p=5: 7 cycles in Omega (5 three-cycles + 2 five-cycles)
       All share vertices (7 cycles on 5 vertices)
       alpha = [1, 7, 0] => I(Omega, 2) = 1 + 14 = 15 = H

  p=7: TWO classes of circulant tournaments
       Paley:     c3=14, c5=42, c7=24 => |Omega| = 80
       Non-Paley: c3=14, c5=28, c7=17 => |Omega| = 59

       H(Paley) = 189, H(non-Paley) = 175, delta = 14
       |Omega|(Paley) - |Omega|(non-Paley) = 80 - 59 = 21
       2 * 21 = 42... not exactly 14.
       But alpha_1 enters as 2*alpha_1, and alpha_2 enters as 4*alpha_2...
       The net effect depends on both linear and quadratic terms.
    """)

    # Let me compute: for Paley vs non-Paley at p=7,
    # what's the contribution of 3-cycles, 5-cycles, 7-cycles to H?
    print(f"  CONTRIBUTION BY CYCLE LENGTH at p=7:")

    for S in [tuple(sorted(set(pow(a,2,7) for a in range(1,7)) - {0})),  # QR = Paley
              (4, 5, 6)]:  # Cyclic interval
        S_set = set(S)
        adj = [[False] * 7 for _ in range(7)]
        for i in range(7):
            for j in range(7):
                if i != j and (j - i) % 7 in S_set:
                    adj[i][j] = True
        H = ham_count_dp(7, S)

        # Find all directed cycles by length
        all_cycles = []
        for k in [3, 5, 7]:
            for subset in combinations(range(7), k):
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

        # Compute H restricted to subsets of cycle lengths
        # I(Omega_{3-only}, 2), I(Omega_{3+5}, 2), I(Omega_{all}, 2)
        for lengths_included in [[3], [3, 5], [3, 5, 7]]:
            restricted = [c for c in all_cycles if c[2] in lengths_included]
            n_cyc = len(restricted)

            if n_cyc > 25:
                print(f"  S={list(S)}, H={H}: lengths {lengths_included}: "
                      f"|Omega|={n_cyc} (too large for IP)")
                continue

            # Build conflict and compute IP
            conflict = [[False] * n_cyc for _ in range(n_cyc)]
            for i in range(n_cyc):
                for j in range(i + 1, n_cyc):
                    if restricted[i][0] & restricted[j][0]:
                        conflict[i][j] = True
                        conflict[j][i] = True

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
            trimmed = [a for a in alphas if a > 0] if any(a > 0 for a in alphas) else [0]
            # Find last nonzero
            last_nz = 0
            for i, a in enumerate(alphas):
                if a > 0:
                    last_nz = i
            alpha_str = str(alphas[:last_nz + 1])
            print(f"    S={list(S)}: lengths={lengths_included}, |Omega|={n_cyc}, "
                  f"alpha={alpha_str}, I(Omega,2)={I_at_2}")

    # ================================================================
    # SECTION 4: Paley cycle dominance
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"SECTION 4: PALEY CYCLE DOMINANCE")
    print(f"{'=' * 60}")
    print("""
    At p=7: Paley has MORE cycles at EVERY length:
      c_3: 14 = 14 (same — all regular)
      c_5: 42 > 28 (+50%)
      c_7: 24 > 17 (+41%)

    At p=11: Same pattern:
      c_3:  55 =  55 (same)
      c_5: 594 > 484-572 (max among non-Paley)
      c_7: 3960 > 3399-3729
      c_9: 11055 > 9350-10274
      c_11: 5505 > 4999-5153

    CONJECTURE: Paley tournament MAXIMIZES c_k for ALL odd k >= 5
    among circulant tournaments on Z_p (when p = 3 mod 4).

    This would give: Paley maximizes alpha_1 = sum c_k.
    By OCF: H = 1 + 2*alpha_1 + 4*alpha_2 + ...
    If alpha_1 dominance is sufficient, Paley maximizes H.

    But we need to check: does max alpha_1 => max H?
    (THM-132 says at p=11, alpha_1 ordering does NOT match H ordering!)
    """)

    # Compute total cycles (alpha_1) for all circulant tournaments at p=7
    p = 7
    print(f"  p = {p}: Total cycles per tournament")
    for S in all_circulant_tournaments(p):
        S_set = set(S)
        adj = [[False] * p for _ in range(p)]
        for i in range(p):
            for j in range(p):
                if i != j and (j - i) % p in S_set:
                    adj[i][j] = True
        H = ham_count_dp(p, S)
        c3 = count_directed_k_cycles(p, S, 3)
        c5 = count_directed_k_cycles(p, S, 5)
        c7 = count_directed_k_cycles(p, S, 7)
        alpha_1 = c3 + c5 + c7
        print(f"    S={list(S)}: H={H}, c3={c3}, c5={c5}, c7={c7}, "
              f"alpha_1={alpha_1}, 2*alpha_1+1={2*alpha_1+1}")


if __name__ == '__main__':
    main()
