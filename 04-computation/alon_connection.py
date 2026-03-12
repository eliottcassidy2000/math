"""
alon_connection.py — Noga Alon's 1990 paper and our findings

CRITICAL DISCOVERY: Alon (Combinatorica, 1990) explicitly constructs the
INTERVAL TOURNAMENT T_n as the best explicit construction for H-maximization:

  T_n on Z_n: (i,j) is edge iff (i-j) mod n < n/2

This is EXACTLY our interval tournament C_p = {1,...,m}.

Alon's results:
  - Upper bound: P(n) ≤ c * n^{3/2} * n!/2^{n-1}  (Theorem 1.1, via Brégman)
  - Lower bound: P(n) ≥ n!/2^{n-1}  (Szele 1943, probabilistic method)
  - Alon's T_n: P(T_n) ≥ n!/(2+o(1))^n  (asymptotically optimal)
  - Szele's conjecture: lim P(n)^{1/n} = n/(2e) — PROVED by Alon

Van der Waerden connection (Remark 1 in Alon):
  For REGULAR tournaments: 2/(n-1) * A_T is doubly stochastic.
  By VdW: Per(A_T) ≥ ((n-1)/2)^n * n!/n^n = (1+o(1)) * (1/e) * n!/2^n

Our contribution: Understanding WHICH tournament achieves maximum at FINITE n.
  - p=7: Paley wins (189 > 175 for interval)
  - p=11: Paley wins (95095 > 93027)
  - p=19: Interval wins (1184212824763 > 1172695746915)
  - OCF mechanism: α_1 vs α_2+ tradeoff through vertex-disjoint cycle packings

The gap: Alon's bound leaves O(n^{3/2}) gap. Can the OCF give tighter bounds?

Author: opus-2026-03-12-S60
"""
import sys
import math
import numpy as np
from collections import defaultdict
from fractions import Fraction
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)


def hamiltonian_paths_dp(A, n):
    dp = defaultdict(lambda: defaultdict(int))
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, full + 1):
        if not dp[mask]:
            continue
        for v in dp[mask]:
            if dp[mask][v] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    return sum(dp[full][v] for v in range(n))


def circulant_adj(n, S):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for s in S:
            A[i][(i+s)%n] = 1
    return A


def paley_set(p):
    return frozenset(pow(x, 2, p) for x in range(1, p))


def interval_set(n):
    m = (n - 1) // 2
    return frozenset(range(1, m + 1))


def main():
    print("ALON'S INTERVAL TOURNAMENT AND OUR OCF ANALYSIS")
    print("=" * 75)

    # Part 1: Verify Alon's bounds
    print("\nPART 1: ALON'S BOUNDS vs ACTUAL VALUES")
    print("-" * 75)

    # OEIS A038375: 1, 1, 3, 5, 15, 45, 189, 661, 3357, 15745, 95095
    known_max = {3: 3, 4: 5, 5: 15, 6: 45, 7: 189, 8: 661, 9: 3357, 10: 15745, 11: 95095}

    print(f"  {'n':>4} {'P(n)':>12} {'n!/2^(n-1)':>12} {'ratio':>8} {'sqrt(pi/2e)*n^1.5':>20} {'UB ratio':>10}")
    for n in range(3, 12):
        P_n = known_max[n]
        szele = math.factorial(n) / 2**(n-1)
        ratio = P_n / szele
        alon_factor = math.sqrt(math.pi / (2 * math.e)) * n**1.5
        ub_ratio = P_n / (alon_factor * szele)
        print(f"  {n:>4} {P_n:>12} {szele:>12.1f} {ratio:>8.4f} {alon_factor:>20.4f} {ub_ratio:>10.6f}")

    # Part 2: Compare Paley, Interval, and Szele bound
    print(f"\n{'='*75}")
    print("PART 2: PALEY vs INTERVAL vs SZELE BOUND")
    print("-" * 75)

    # Known H values
    paley_h = {3: 3, 7: 189, 11: 95095, 19: 1172695746915}
    interval_h = {3: 3, 7: 175, 11: 93027, 19: 1184212824763}

    print(f"  {'p':>4} {'H(Paley)':>20} {'H(Interval)':>20} {'Szele':>20} {'P/Sz':>8} {'I/Sz':>8}")
    for p in [3, 7, 11, 19]:
        szele = math.factorial(p) / 2**(p-1)
        hp = paley_h[p]
        hi = interval_h[p]
        print(f"  {p:>4} {hp:>20} {hi:>20} {szele:>20.1f} {hp/szele:>8.4f} {hi/szele:>8.4f}")

    # Part 3: Brégman bound on 1-factors
    print(f"\n{'='*75}")
    print("PART 3: BRÉGMAN (MINC) BOUND ON 1-FACTORS")
    print("-" * 75)

    # For regular tournament: all row sums = (n-1)/2 = k
    # Brégman: Per(A) ≤ Π (r_i!)^{1/r_i} = (k!)^{n/k}
    # For k = (n-1)/2: (((n-1)/2)!)^{2n/(n-1)}

    for n in [7, 11, 19]:
        k = (n - 1) // 2
        bregman_ub = (math.factorial(k) ** (1/k)) ** n
        szele_avg = math.factorial(n) / 2**(n-1)
        vdw_lb = ((n-1)/2)**n * math.factorial(n) / n**n
        print(f"  n={n}: k={k}")
        print(f"    Brégman UB for 1-factors: {bregman_ub:.2e}")
        print(f"    VdW LB for regular: {vdw_lb:.2e}")
        print(f"    Szele avg HP: {szele_avg:.2e}")
        print(f"    Brégman/Szele ratio: {bregman_ub/szele_avg:.4f}")
        print(f"    Alon UB factor: sqrt(π/2e)*n^1.5 = {math.sqrt(math.pi/(2*math.e))*n**1.5:.4f}")

    # Part 4: OCF perspective on Alon's bound
    print(f"\n{'='*75}")
    print("PART 4: OCF PERSPECTIVE ON THE n^{3/2} GAP")
    print("-" * 75)

    print("""
  Alon's proof chain:
    P(T) ≤ n * C(T)  (any HP becomes a HC by adding one edge)
    C(T) ≤ F(T) = Per(A_T)  (HC ⊂ 1-factors)
    Per(A_T) ≤ Π (r_i!)^{1/r_i}  (Brégman/Minc)
    → P(n) ≤ c * n^{3/2} * n!/2^{n-1}

  OCF perspective:
    H(T) = I(Ω(T), 2) = Σ_{k≥0} α_k * 2^k

  The α_k encode the independence structure of the odd-cycle complex.
  For regular tournaments (both Paley and Interval):
    α_0 = 1
    α_1 = # odd cycles (grows roughly as n!/n for the n-cycle component)
    α_2 = # vertex-disjoint cycle pairs
    ...

  The n^{3/2} factor in Alon's bound comes from:
    - n factor from HP → HC conversion
    - √n factor from Brégman's inequality applied to regular matrix

  Can the OCF give a TIGHTER bound?
  The key constraint: independence number of Ω(T) is at most ⌊n/3⌋
  (since each odd cycle uses ≥ 3 vertices and they must be disjoint).
  So H ≤ 1 + 2^1*α_1 + 2^2*α_2 + ... + 2^{⌊n/3⌋}*α_{⌊n/3⌋}
  where α_k ≤ C(total_cycles, k) but with vertex-disjointness constraint.
    """)

    # Part 5: The EXACT growth rate comparison
    print(f"{'='*75}")
    print("PART 5: GROWTH RATE P(n)^{1/n}")
    print("-" * 75)

    # Szele proved and Alon confirmed: lim P(n)^{1/n} = n/(2e)
    # Both Paley and Interval achieve this rate
    print(f"  Szele-Alon: lim P(n)^{{1/n}} = n/(2e)")
    print()
    print(f"  {'n':>4} {'P(n)':>20} {'P(n)^(1/n)':>12} {'n/(2e)':>8}")
    for n, P_n in sorted(known_max.items()):
        growth = P_n ** (1/n)
        szele_limit = n / (2 * math.e)
        print(f"  {n:>4} {P_n:>20} {growth:>12.6f} {szele_limit:>8.4f}")

    # Part 6: Compare interval with OEIS values
    print(f"\n{'='*75}")
    print("PART 6: IS INTERVAL THE GLOBAL MAX FOR ALL n ≥ 19?")
    print("-" * 75)

    # Compute H(interval) at small n and compare with OEIS
    for n in [3, 4, 5, 6, 7, 8, 9, 10, 11]:
        m = (n - 1) // 2
        S = interval_set(n)
        A = circulant_adj(n, S)
        H_int = hamiltonian_paths_dp(A, n)
        H_max = known_max.get(n, '?')
        match = "= MAX" if H_int == H_max else f"< MAX ({H_max})"
        print(f"  n={n}: H(interval) = {H_int:>8}, P(n) = {str(H_max):>8}  {match}")

    # Part 7: Alon's conjecture about C(T_n) = C(n)
    print(f"\n{'='*75}")
    print("PART 7: ALON'S OPEN CONJECTURE: C(T_n) = C(n)?")
    print("-" * 75)
    print("  Alon conjectures the interval tournament maximizes Hamiltonian CYCLES.")
    print("  'It seems plausible that C(T_n) = C(n), but at the moment we are")
    print("  unable to prove or disprove this statement.' (1990)")
    print()
    print("  Our OCF analysis suggests this IS true for large n:")
    print("  The interval's local structure creates better vertex-disjoint")
    print("  cycle packings, and the OCF's 2^k weighting amplifies this.")
    print("  The crossover at p=19 for paths likely also holds for cycles.")


if __name__ == '__main__':
    main()
    print("\nDONE.")
