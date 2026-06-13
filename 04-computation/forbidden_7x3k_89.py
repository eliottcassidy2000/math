#!/usr/bin/env python3
"""
forbidden_7x3k_89.py — opus-2026-03-14-S89
The forbidden value pattern {7, 21, 63} = 7 × 3^k.

Key discovery: At n=7, the permanently forbidden small values are:
  7 = 7 × 1 = 7 × 3^0
  21 = 7 × 3 = 7 × 3^1
  63 = 7 × 9 = 7 × 3^2

Question: Is this pattern fundamental? Does 189 = 7 × 27 = 7 × 3^3 survive?
(Answer: 189 IS achievable at n=7!)

This script analyzes WHY 7×3^k is forbidden for k=0,1,2 but not k=3.

STRUCTURAL ANALYSIS:
  7 = Φ₃(2) = 2²+2+1 = 111₂
  21 = Φ₃(4) = 4²+4+1 = 111₄
  63 = 2⁶-1 = (2³-1)(2³+1) = 7×9

  Alternative:
  7 = (2³-1)/(2-1) = repunit R₃(2)
  21 = (4³-1)/(4-1) = repunit R₃(4)
  63 = (4³-1)/(4-1) × 3 = 3×R₃(4) ??
  No: 63 = (2⁶-1)/(2-1) = R₆(2) = repunit with 6 ones in base 2

  Actually: R₃(2) = 7, R₃(4) = 21
  And R₆(2) = 63 = R₃(2) × R₃(8)/? No, R₃(8) = 73 ≠ 63/7 = 9.

  The factorization: 63 = 7 × 9 = R₃(2) × (2³+1) = (2³-1)(2³+1) = 2⁶-1
"""

from itertools import permutations
from collections import Counter
from fractions import Fraction
import time

def compute_H_dp(adj_bits, n):
    """Compute H using Hamiltonian path DP."""
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            count = dp[mask][v]
            remaining = ((1 << n) - 1) & ~mask
            targets = adj_bits[v] & remaining
            u = targets
            while u:
                bit = u & (-u)
                idx = bit.bit_length() - 1
                dp[mask | bit][idx] += count
                u &= u - 1
    full = (1 << n) - 1
    return sum(dp[full])


def get_spectrum(n):
    """Get full H-spectrum for tournaments on n vertices."""
    m = n * (n - 1) // 2
    total = 1 << m
    edges = []
    for u in range(n):
        for v in range(u+1, n):
            edges.append((u, v))

    h_counts = Counter()
    for idx in range(total):
        adj_bits = [0] * n
        for e, (u, v) in enumerate(edges):
            if idx & (1 << e):
                adj_bits[u] |= (1 << v)
            else:
                adj_bits[v] |= (1 << u)
        h = compute_H_dp(adj_bits, n)
        h_counts[h] += 1
    return h_counts


def main():
    print("="*70)
    print("FORBIDDEN VALUE PATTERN: 7 × 3^k")
    print("opus-2026-03-14-S89")
    print("="*70)

    # Load precomputed spectra for n=3..6, compute for n=7 from saved data
    spectra = {}

    for n in range(3, 8):
        print(f"\n  Computing H-spectrum at n={n}...")
        t0 = time.time()
        if n <= 6:
            spectra[n] = get_spectrum(n)
        else:
            # Use known results from exhaustive computation
            # (already verified in n7_deep_spectrum_89.py)
            spectra[7] = {
                1:5040, 3:8400, 5:20160, 9:31920, 11:15120, 13:15120,
                15:31584, 17:20160, 19:20160, 23:40320, 25:25200,
                27:21840, 29:50400, 31:30240, 33:52080, 35:10080,
                37:60480, 39:10080, 41:35280, 43:30240, 45:25200,
                47:10080, 49:35280, 51:30240, 53:20160, 55:20160,
                57:13440, 59:40320, 61:5040, 65:55440, 67:30240,
                69:40320, 71:35280, 73:40320, 75:28560, 77:50400,
                79:25200, 81:55440, 83:35280, 85:5040, 87:5040,
                89:70560, 91:75600, 93:35280, 95:25200, 97:30240,
                99:12320, 101:50400, 103:50400, 105:20160, 109:70560,
                111:47040, 113:20160, 115:35280, 117:35280, 121:20160,
                123:40320, 125:5040, 127:5040, 129:31920, 131:40320,
                133:30240, 135:29568, 137:30240, 139:30240, 141:15120,
                143:15120, 145:10080, 147:5040, 151:25200, 153:15120,
                155:5040, 157:5040, 159:16800, 171:1680, 175:720, 189:240
            }
        elapsed = time.time() - t0
        h_vals = sorted(spectra[n].keys())
        print(f"    Done in {elapsed:.1f}s, {len(h_vals)} distinct H values")

    # Check 7×3^k pattern across all n
    print(f"\n" + "="*70)
    print(f"PART 1: STATUS OF 7 × 3^k VALUES")
    print("="*70)

    targets = [7 * 3**k for k in range(6)]  # 7, 21, 63, 189, 567, 1701
    for k, val in enumerate(targets):
        print(f"\n  7 × 3^{k} = {val}:")
        for n in range(3, 8):
            count = spectra[n].get(val, 0)
            max_h = max(spectra[n].keys())
            in_range = val <= max_h
            status = "ACHIEVABLE" if count > 0 else ("FORBIDDEN (in range)" if in_range else "out of range")
            print(f"    n={n}: max_H={max_h:5d}, count={count:8d}, {status}")

    # Analyze forbidden structure
    print(f"\n" + "="*70)
    print(f"PART 2: REPUNIT AND CYCLOTOMIC STRUCTURE")
    print("="*70)

    print(f"\n  Forbidden values and their representations:")
    print(f"    7 = R₃(2) = Φ₃(2) = 2²+2+1 = 111₂")
    print(f"    21 = R₃(4) = Φ₃(4) = 4²+4+1 = 111₄ = 3×7")
    print(f"    63 = R₆(2) = 2⁶-1 = 111111₂ = 333₄ = 9×7")
    print(f"    189 = 7×27 = 7×3³ — ACHIEVABLE at n=7!")

    print(f"\n  WHY 189 survives:")
    print(f"    189 in binary: {bin(189)[2:]}")
    print(f"    189 in base 4: ", end="")
    b4 = ""
    h = 189
    while h > 0:
        b4 = str(h % 4) + b4
        h //= 4
    print(f"{b4}₄")

    # 189 = 2*64 + 3*16 + 3*4 + 1 = 128+48+12+1 = 189
    # in base 4: 189 = 2*64 + 3*16 + 3*4 + 1 = 2331₄
    # Not a repdigit!

    print(f"\n  The key: 63 = 333₄ (repdigit in base 4), but 189 = 2331₄ (NOT repdigit)")
    print(f"  The repdigit/repunit connection is the obstruction, not just 7×3^k!")

    # Check all repdigits in base 4
    print(f"\n  Repdigits in base 4:")
    for d in range(1, 4):
        for k in range(1, 6):
            val = d * sum(4**i for i in range(k))
            print(f"    {''.join([str(d)]*k)}₄ = {val}", end="")
            # Check forbidden status
            if val <= 189:
                for nn in range(3, 8):
                    if val in spectra[nn] or val > max(spectra[nn].keys()):
                        continue
                    if val not in spectra[nn]:
                        print(f" — FORBIDDEN at n={nn}", end="")
                        break
            print()

    # Level-2 analysis: what is the "level" structure of 63?
    print(f"\n" + "="*70)
    print(f"PART 3: OCF LEVEL DECOMPOSITION")
    print("="*70)

    # H(T) = n! * Σ c_S * (-1)^|S| / 2^|S|
    # Actually H = 1 + Σ (over distinct Fourier levels)
    # The variance decomposes as Var = Σ E_{2k}
    # E_2 = level-2 energy, E_4 = level-4 energy, etc.

    # At n=7, Var/Mean² = 131/504
    # Level-2: E_2/E_0 = ?
    # We know E_2/E_0 = 4/C(n,2) at n≤4 = 1/3
    # At n=5: E_2/E_0 = 3/10 (from kind-pasteur)
    # At n=6: E_2/E_0 = 4/15 (from kind-pasteur)

    # The pattern for E_2/E_0 = 2/(m-1) where m = C(n,2)
    # No wait. Let me recalculate from the master formula.

    # From level4_counting.out:
    # E_2/E_0 = 4*N_2/(P(n,2))² where N_2 = n-covering 2-edge subsets
    # P(n,2) = n(n-1)

    # N_2(n) = number of 2-edge subsets of K_n that cover all n vertices
    # For n=3: need 2 edges covering 3 vertices → each edge covers 2 verts,
    # 2 edges cover at most 4 verts. For 3 verts: C(3,2)=3 edges, need 2 covering all 3.
    # From 3 edges, 2-subsets: C(3,2)=3. Each pair covers 3 verts iff they share a vertex.
    # All pairs share a vertex (triangle), so N_2(3) = 3.
    # E_2/E_0 = 4*3/(3*2)² = 12/36 = 1/3. ✓

    # N_2(n) = C(m,2) - C(n,2)*C(m-3,0)... too complicated.
    # Actually N_2(n) = C(m,2) - n*C(m-n+1, 2)... no.
    # N_2(n) = total 2-edge subsets minus those missing at least one vertex
    # = C(m,2) - Σ C(m_i, 2) + ... (inclusion-exclusion)
    # where m_i = edges not touching vertex i = C(n-1,2)

    # Let me just compute N_2(n) directly
    print(f"\n  Computing N_2(n) for n=3..7:")
    for n in range(3, 8):
        edges_n = []
        for u in range(n):
            for v in range(u+1, n):
                edges_n.append((u, v))
        mm = len(edges_n)

        count = 0
        for i in range(mm):
            for j in range(i+1, mm):
                # Check if edges i,j cover all n vertices
                verts = set(edges_n[i]) | set(edges_n[j])
                if len(verts) == n:
                    count += 1

        pn2 = n * (n-1)
        e2_e0 = Fraction(4 * count, pn2**2)
        print(f"    n={n}: N_2={count}, P(n,2)={pn2}, E_2/E_0 = {e2_e0} = {float(e2_e0):.10f}")

    # Now compute N_4(n) for n=3..7
    print(f"\n  Computing N_4(n) for n=3..7:")
    for n in range(3, 8):
        edges_n = []
        for u in range(n):
            for v in range(u+1, n):
                edges_n.append((u, v))
        mm = len(edges_n)

        if mm >= 4:
            from itertools import combinations
            count = 0
            for combo in combinations(range(mm), 4):
                verts = set()
                for idx in combo:
                    verts |= set(edges_n[idx])
                if len(verts) == n:
                    count += 1

            pn4 = n * (n-1) * (n-2) * (n-3) if n >= 4 else 1
            e4_e0 = Fraction(4 * count, pn4**2) if n >= 4 else Fraction(0)
            print(f"    n={n}: N_4={count}, P(n,4)={pn4}, E_4/E_0 = {e4_e0} = {float(e4_e0):.10f}")
        else:
            print(f"    n={n}: not enough edges for 4-subsets")

    # N_6 at n=7
    print(f"\n  Computing N_6(n=7) — this may take a moment...")
    n = 7
    edges_n = []
    for u in range(n):
        for v in range(u+1, n):
            edges_n.append((u, v))
    mm = len(edges_n)  # 21

    from itertools import combinations
    t0 = time.time()
    count_6 = 0
    for combo in combinations(range(mm), 6):
        verts = set()
        for idx in combo:
            verts |= set(edges_n[idx])
        if len(verts) == n:
            count_6 += 1
    elapsed = time.time() - t0
    pn6 = n * (n-1) * (n-2) * (n-3) * (n-4) * (n-5)
    e6_e0 = Fraction(4 * count_6, pn6**2)
    print(f"    n=7: N_6={count_6}, P(n,6)={pn6}, E_6/E_0 = {e6_e0} = {float(e6_e0):.10f}")
    print(f"    (computed in {elapsed:.1f}s)")

    # Check: E_2 + E_4 + E_6 + ... should give Var/Mean²
    print(f"\n  Computing master formula check at n=7:")
    # Need E_2, E_4, E_6, ...
    # E_2k/E_0 = 4 * N_{2k}(n) / P(n,2k)^2

    # We have N_2(7), N_4(7), N_6(7)
    # Need N_8, N_10, ... up to N_{2*floor(m/2)} but limited by coverage

    # Actually the maximum useful level is 2k where P(n,2k) is defined, i.e. 2k ≤ n
    # For n=7: 2k ≤ 7, so k ≤ 3, levels 2, 4, 6
    # Wait, P(n,2k) = n!/(n-2k)! requires n ≥ 2k
    # n=7: 2k can be 2, 4, 6 (since 2k ≤ 7)

    # Actually 2k ≤ m = 21 edges, but the covering condition requires 2k ≥ n-1 = 6
    # edges to cover n=7 vertices. So N_{2k}=0 for 2k < 6? No, that's not right.
    # N_{2k} = number of 2k-edge subsets covering all n vertices.
    # For n=7: need at least ceil(7/2) = 4 edges to cover 7 vertices.
    # So N_2(7) = 0 (2 edges can cover at most 4 vertices).
    # N_4(7) = ? (4 edges can cover at most 8 vertices, so could cover 7)

    # Let me recompute N_2 more carefully
    print(f"\n  RECHECK: Does N_2(7)=0?")
    n = 7
    edges_n = []
    for u in range(n):
        for v in range(u+1, n):
            edges_n.append((u, v))
    mm = len(edges_n)
    count_2 = 0
    for i in range(mm):
        for j in range(i+1, mm):
            verts = set(edges_n[i]) | set(edges_n[j])
            if len(verts) == n:
                count_2 += 1
    print(f"    N_2(7) = {count_2}")  # Should be 0 since 2 edges cover at most 4 vertices

    # So the master formula at n=7:
    # Var/Mean² = E_2/E_0 + E_4/E_0 + E_6/E_0 + ...
    # But wait, E_2 at n=7 means level-2 Fourier, which is |S|=2 subsets
    # where S is a set of ARCS (oriented edges), not undirected edges.
    # I might be confusing the formulas. Let me re-derive.

    # Actually from level4_counting.out, the master formula is:
    # Var/Mean² = 4 * sum_k N_{2k}(n) / P(n,2k)^2
    # where the sum is over k=1,2,...

    # But P(n,2k) = n!/(n-2k)! (falling factorial)
    # For n=7: P(7,2)=42, P(7,4)=840, P(7,6)=5040

    # Hmm wait, from level4_counting.out:
    # E_{2k}/E_0 = 4 * N_{2k}(n) / P(n,2k)^2
    # n=5: E_2/E_0 = 4*30/20^2 = 120/400 = 3/10 ✓

    # But N_2(n) = covering 2-edge subsets of K_n
    # At n=7: 2 edges cover at most 4 vertices, so N_2(7) = 0
    # => E_2/E_0 = 0 at n=7?? That can't be right.

    # Let me reread the formula. The "2k-arc" subsets are ordered edges.
    # An arc is an ORDERED pair (i,j) with i→j. There are P(n,2) = n(n-1) arcs.
    # A 2-arc subset = 2 ordered edges. "Covering" means every vertex appears.
    # But 2 arcs have at most 4 endpoints, can't cover 7 vertices.
    # So N_2(7) = 0 is correct for the arc formulation too.

    # This means E_2/E_0 = 0 at n=7, and the FIRST nonzero level is E_4 or higher.
    # Wait, that contradicts the known result. The level-2 coefficient c_2(S)
    # for |S|=1 (single arc) is always nonzero...

    # I think the confusion is: the "level" in the master formula refers to
    # PAIRS of arcs in the Fourier expansion, not just single arcs.
    # E_{2k} is the energy in the 2k-arc correlation.
    # But the variance involves PAIRS of paths, so the correlation is over
    # pairs of arc-subsets.

    # Actually, I think the issue is that N_{2k} counts subsets of 2k ARCS
    # from the complete DIRECTED graph on n vertices (which has n(n-1) arcs),
    # not undirected edges.

    # Let me reconsider. From level4_counting.out:
    # P(n,2) = n(n-1), N_2(n=5)=30, E_2/E_0 = 4*30/20^2 = 3/10
    # P(5,2) = 20, and 30 = ? "2-arc subsets covering all 5 vertices"
    # 2 arcs from K_5 directed: there are P(5,2)*... no.
    # 2 arcs = 2 directed edges. Total = C(20,2) = 190.
    # Those covering all 5 vertices: need 2 arcs whose 4 endpoints include
    # all 5 vertices. But 2 arcs have at most 4 vertices. Can't cover 5.
    # Unless arcs can share an endpoint.
    # 2 arcs: (a→b), (c→d). Vertices involved: {a,b,c,d} ≤ 4.
    # Can't be 5. So N_2(5) should be 0, not 30!

    # Something is wrong with my understanding of the formula.
    # Let me re-read level4_counting.out more carefully.

    # From the output:
    # k=1, n=5: P(n,2)=20, N_2=30, E_2/E_0 = 0.3000000000
    # But 0.3 = 3/10 which IS the level-2 energy.
    # So what is N_2 = 30 counting?

    # Maybe N_2 = "number of UNDIRECTED 2-edge subsets covering all n vertices"?
    # For n=5: 2 edges covering 5 vertices is impossible.

    # OR: N_2 could be the number of 2-ELEMENT subsets of {1,...,m}
    # (where m = C(n,2)) covering all vertices.
    # Same problem: 2 undirected edges cover at most 4 vertices.

    # WAIT: I think "covering" might mean something different here.
    # Maybe it means "all vertices appear as an endpoint" which for
    # 2 edges is impossible if n ≥ 5. But the output says N_2(5) = 30.

    # Let me look at the actual code that computed this.
    print(f"\n  FORMULA INVESTIGATION:")
    print(f"  From level4_counting.out: N_2(5)=30, N_2(6)=60")
    print(f"  But 2 edges can cover at most 4 vertices!")
    print(f"  So 'N_2' must mean something else.")
    print(f"  ")
    print(f"  Possible interpretation: N_2 is not edge-covering,")
    print(f"  but rather the number of 2-element subsets S of [m]")
    print(f"  such that the corresponding flip changes ALL vertices.")
    print(f"  Or: N_2 counts something in the Fourier expansion.")

    # Let me go back to basics. The variance formula is:
    # Var(H) = Σ_{S⊂[m], |S|≥1} ĉ_S²
    # where ĉ_S is the Fourier coefficient of H on subset S.
    # E_{2k}/E_0 = (1/Mean²) Σ_{|S|=k} ĉ_S²

    # At level 1: E_2/E_0 = (1/Mean²) Σ_{|S|=1} ĉ_S²
    # Each single-arc Fourier coefficient ĉ_{e} can be computed.
    # By symmetry, all |ĉ_e| are equal at each n.
    # So E_2/E_0 = m * ĉ_e² / Mean²

    # At n=5: ĉ_e² = ?, m=10
    # E_2/E_0 = 3/10 → ĉ_e² = 3/10 * Mean² / m
    # Mean = 5!/2^4 = 120/16 = 15/2
    # ĉ_e² = 3/10 * (15/2)² / 10 = 3/10 * 225/4 / 10 = 3*225/400 = 675/400 = 27/16
    # ĉ_e = ±√(27/16) = ±3√3/4

    # Hmm, let me reconsider. Maybe N_2 and the formula are about
    # something else entirely. Let me just verify the master formula
    # directly for n=7.

    # Direct computation: for n=7, we have Var/Mean² = 131/504
    # Mean = 315/4
    # Var = 206325/128

    # Fourier decomposition of Var:
    # Var = Σ_{|S|≥1} ĉ_S²
    # Group by |S|: Var = Σ_{k=1}^{21} E_k where E_k = Σ_{|S|=k} ĉ_S²

    # The "level-2" in the master formula seems to be about |S|=1 (single arcs)
    # mapped to "2-arc" because each arc involves 2 vertices.
    # So "N_2" might be m itself, and the formula has a different form.

    # Let me just compute level-by-level for small n to verify.
    print(f"\n  Level-by-level Fourier energy for n=5:")
    n = 5
    m = n*(n-1)//2  # 10
    total = 1 << m
    edges_n = []
    for u in range(n):
        for v in range(u+1, n):
            edges_n.append((u, v))

    # Compute H for all tournaments
    h_vals_n5 = []
    for idx in range(total):
        adj_bits = [0] * n
        for e, (u, v) in enumerate(edges_n):
            if idx & (1 << e):
                adj_bits[u] |= (1 << v)
            else:
                adj_bits[v] |= (1 << u)
        h = compute_H_dp(adj_bits, n)
        h_vals_n5.append(h)

    mean_n5 = sum(h_vals_n5) / total
    print(f"    Mean = {mean_n5}")

    # Compute Fourier coefficients
    # ĉ_S = (1/2^m) Σ_T (-1)^{<S,T>} H(T)
    # where <S,T> = number of edges in S that are oriented the "wrong" way in T
    # Actually in Walsh-Hadamard: ĉ_S = (1/2^m) Σ_x (-1)^{S·x} f(x)
    # where x is the m-bit tournament vector.

    # For each subset S of [m], compute ĉ_S
    # Then group by |S|

    level_energy = [0.0] * (m+1)
    for S in range(1, 1 << m):
        # Compute ĉ_S
        c_S = 0
        for idx in range(total):
            # (-1)^{popcount(S & idx)}
            sign = 1 - 2*(bin(S & idx).count('1') % 2)
            c_S += sign * h_vals_n5[idx]
        c_S /= total
        level = bin(S).count('1')
        level_energy[level] += c_S * c_S

    print(f"    Fourier energy by level:")
    var_sum = 0
    for k in range(1, m+1):
        if level_energy[k] > 0.001:
            ratio = level_energy[k] / mean_n5**2
            print(f"      level {k:2d}: E_{k} = {level_energy[k]:12.4f}, E_{k}/E_0 = {ratio:.10f}")
            var_sum += level_energy[k]

    print(f"    Total Var = {var_sum:.4f}")
    print(f"    Var/Mean² = {var_sum/mean_n5**2:.10f}")
    print(f"    Expected: {float(Fraction(19,60)):.10f}")

    # Now for n=6 (this will be slower: 2^15 = 32768 tournaments, 2^15 subsets)
    # Actually 2^15 subsets × 2^15 tournaments = 2^30 ≈ 10^9, too slow.
    # Skip n=6 and just report n=5 findings.

    print(f"\n" + "="*70)
    print(f"PART 4: CONE AND FORBIDDEN VALUE MECHANISM")
    print("="*70)

    # Key insight: H(dominant_cone(T)) = H(T) (cone H-preservation)
    # So if H=7 is forbidden for n=k, it's forbidden for n=k+1 via cone?
    # No: cones PRESERVE H, but non-cone tournaments could achieve H=7.
    # The question is: can ANY non-cone tournament at n=k+1 have H=7?

    # At n=4: H values are {1, 3} (no 7 since max H=3)
    # At n=5: H values are {1, 3, 5, 9, 11, 13, 15} (no 7)
    # At n=6: H values are {1, 3, 5, 9, ...} — need to check

    print(f"\n  Forbidden value 7 across n:")
    for nn in range(3, 8):
        h_vals_n = sorted(spectra[nn].keys())
        has_7 = 7 in spectra[nn]
        max_h = max(h_vals_n)
        print(f"    n={nn}: max_H={max_h}, H=7 present: {has_7}")

    print(f"\n  Forbidden value 21 across n:")
    for nn in range(3, 8):
        h_vals_n = sorted(spectra[nn].keys())
        has_21 = 21 in spectra[nn]
        max_h = max(h_vals_n)
        print(f"    n={nn}: max_H={max_h}, H=21 present: {has_21}")

    print(f"\n  Forbidden value 63 across n:")
    for nn in range(5, 8):
        h_vals_n = sorted(spectra[nn].keys())
        has_63 = 63 in spectra[nn]
        max_h = max(h_vals_n)
        print(f"    n={nn}: max_H={max_h}, H=63 present: {has_63}")

    # The spectrum sizes
    print(f"\n  H-spectrum sizes: ", end="")
    sizes = []
    for nn in range(3, 8):
        s = len(spectra[nn])
        sizes.append(s)
        print(f"n={nn}→{s} ", end="")
    print()
    print(f"  Sizes: {sizes}")
    print(f"  Ratios: ", end="")
    for i in range(1, len(sizes)):
        print(f"{sizes[i]/sizes[i-1]:.4f} ", end="")
    print()

    # GCD of counts at each n
    print(f"\n  GCD of spectrum counts:")
    from math import gcd
    from functools import reduce
    for nn in range(3, 8):
        counts = list(spectra[nn].values())
        g = reduce(gcd, counts)
        print(f"    n={nn}: GCD = {g}")

    # The 63 = 333₄ insight
    print(f"\n" + "="*70)
    print(f"PART 5: THE REPDIGIT HYPOTHESIS")
    print("="*70)

    # Hypothesis: H is forbidden iff H is a repdigit in base 2 or base 4
    # that falls in [1, max_H] and is odd
    print(f"\n  Repdigits in base 2 (repunits): 1, 3, 7, 15, 31, 63, 127, 255, ...")
    print(f"    = 2^k - 1 for k=1,2,3,4,5,6,7,8,...")
    print(f"    Odd ones: all of them")
    print(f"    Which are forbidden?")

    repunits_2 = [2**k - 1 for k in range(1, 10)]
    for r in repunits_2:
        statuses = []
        for nn in range(3, 8):
            max_h = max(spectra[nn].keys())
            if r > max_h:
                statuses.append(f"n={nn}:OOR")
            elif r in spectra[nn]:
                statuses.append(f"n={nn}:OK")
            else:
                statuses.append(f"n={nn}:FORB")
        print(f"    R_{r} = {r}: {', '.join(statuses)}")

    print(f"\n  Repdigits in base 4:")
    for d in [1, 3]:  # only odd digits
        for k in range(1, 5):
            val = d * sum(4**i for i in range(k))
            statuses = []
            for nn in range(3, 8):
                max_h = max(spectra[nn].keys())
                if val > max_h:
                    statuses.append(f"n={nn}:OOR")
                elif val in spectra[nn]:
                    statuses.append(f"n={nn}:OK")
                else:
                    statuses.append(f"n={nn}:FORB")
            rep = str(d) * k
            print(f"    {rep}₄ = {val}: {', '.join(statuses)}")

    # KEY: at n=7, which repunits are forbidden?
    # 1 = 1₂ → ACHIEVABLE (H=1 always exists)
    # 3 = 11₂ → ACHIEVABLE
    # 7 = 111₂ → FORBIDDEN
    # 15 = 1111₂ → ACHIEVABLE (in n=7 spectrum)
    # 31 = 11111₂ → ACHIEVABLE
    # 63 = 111111₂ → FORBIDDEN!
    # 127 = 1111111₂ → ACHIEVABLE!

    # So 7 and 63 are forbidden repunits, but 15, 31, 127 are not.
    # 7 = R₃(2), 63 = R₆(2). What's special about 3 and 6?
    # 3 and 6 are both divisible by 3!
    # R₃(2) = 7, R₆(2) = 63, R₉(2) = 511
    # Repunits with digit count divisible by 3?

    print(f"\n  Repunits R_k(2) = 2^k - 1 where k ≡ 0 mod 3:")
    for k in range(3, 15, 3):
        val = 2**k - 1
        print(f"    R_{k}(2) = {val} = {' × '.join(str(f) for f in prime_factors(val))}")
        if val <= 189:
            in_spec = val in spectra[7]
            print(f"      At n=7: {'ACHIEVABLE' if in_spec else 'FORBIDDEN'}")

    print(f"\n" + "="*70)
    print(f"DONE — FORBIDDEN 7×3^k ANALYSIS")
    print("="*70)


def prime_factors(n):
    """Return list of prime factors."""
    factors = []
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors.append(d)
            n //= d
        d += 1
    if n > 1:
        factors.append(n)
    return factors


if __name__ == "__main__":
    main()
