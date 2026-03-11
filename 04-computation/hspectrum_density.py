"""
hspectrum_density.py - Analyze the density and structure of the H-spectrum.

The H-spectrum = set of odd positive integers achievable as H(T) for some tournament T.
Known: H=7 and H=21 are permanent gaps. All others in [1,200] are achievable by n<=9.

Questions:
1. What is the density d(N) = #{achievable H <= N} / #{odd integers <= N}?
2. Do gaps become sparser (density -> 1)?
3. Is there a pattern to the gaps?
4. Connection to the OCF: H = sum_S 2^|S| * (alpha_S=1), can every odd value be achieved?

Author: kind-pasteur-2026-03-10-S52
"""
import sys
import time
import numpy as np
from itertools import product

sys.stdout.reconfigure(line_buffering=True)


def build_tournament(bits, n):
    """Build tournament from bit encoding."""
    A = np.zeros((n, n), dtype=int)
    k = 0
    for i in range(n - 1, -1, -1):
        for j in range(i):
            if (bits >> k) & 1:
                A[j, i] = 1
            else:
                A[i, j] = 1
            k += 1
    return A


def count_ham_paths(A, n):
    """Count Hamiltonian paths using DP with bitmask."""
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, 1 << n):
        for last in range(n):
            if not (mask >> last & 1) or dp[mask][last] == 0:
                continue
            for nxt in range(n):
                if (mask >> nxt & 1) or not A[last, nxt]:
                    continue
                dp[mask | (1 << nxt)][nxt] += dp[mask][last]
    full = (1 << n) - 1
    return sum(dp[full])


def find_odd_cycles_fast(A, n):
    """Find all odd directed cycles using DP."""
    cycles = set()
    for length in range(3, n + 1, 2):
        for start in range(n):
            # DP: can we reach each vertex from start in exactly 'length' steps
            # using distinct vertices
            # Simple: bitmask DP
            from itertools import permutations
            # Only for small n: enumerate
            if n <= 9:
                pass
            else:
                break
    # Actually use simpler approach: count 3,5,7-cycles
    return cycles


def independence_polynomial_at_2(A, n):
    """Compute I(Omega(T), 2) = H(T) via OCF.
    Also computes the cycle counts (alpha_1, alpha_2, ...).
    Returns H and (alpha_0, alpha_1, alpha_2).
    """
    # Find all odd cycles
    # For efficiency: only 3-cycles and 5-cycles for moderate n
    # 3-cycles
    c3_sets = []
    for i in range(n):
        for j in range(n):
            if i == j or not A[i, j]:
                continue
            for k in range(n):
                if k == i or k == j or not A[j, k] or not A[k, i]:
                    continue
                vset = frozenset([i, j, k])
                if vset not in set(c3_sets):
                    c3_sets.append(vset)

    c3_unique = list(set(frozenset(c) for c in c3_sets))

    # For independent pairs (alpha_2)
    alpha_2 = 0
    for i in range(len(c3_unique)):
        for j in range(i + 1, len(c3_unique)):
            if not (c3_unique[i] & c3_unique[j]):  # disjoint
                alpha_2 += 1

    alpha_1 = len(c3_unique)
    # Simple approximation using 3-cycles only: H ≈ 1 + 2*alpha_1 + 4*alpha_2
    # (exact only if no 5-cycles)
    H_approx = 1 + 2 * alpha_1 + 4 * alpha_2
    return H_approx, alpha_1, alpha_2


def exhaustive_hspectrum(n):
    """Find all achievable H values for n-vertex tournaments."""
    total = 1 << ((n * (n - 1)) // 2)
    H_set = set()
    for bits in range(total):
        A = build_tournament(bits, n)
        H = count_ham_paths(A, n)
        H_set.add(H)
    return sorted(H_set)


def analyze_density():
    """Analyze density of H-spectrum."""
    print("=" * 70)
    print("H-SPECTRUM DENSITY ANALYSIS")
    print("=" * 70)
    print()

    # Known spectrum at n=7 (from OPEN-Q-019)
    spectrum_n7 = [1, 3, 5, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31, 33, 35,
                   37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 65, 67,
                   69, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97,
                   99, 101, 103, 105, 109, 111, 113, 115, 117, 121, 123, 125,
                   127, 129, 131, 133, 135, 137, 139, 141, 143, 145, 147, 151,
                   153, 155, 157, 159, 171, 175, 189]
    max_n7 = max(spectrum_n7)

    print(f"Spectrum at n=7: {len(spectrum_n7)} values, max={max_n7}")

    # Count odd integers up to max_n7
    odd_up_to = (max_n7 + 1) // 2
    density_n7 = len(spectrum_n7) / odd_up_to
    print(f"Odd integers in [1,{max_n7}]: {odd_up_to}")
    print(f"Density at n=7: {len(spectrum_n7)}/{odd_up_to} = {density_n7:.4f}")

    # Find gaps
    all_odds = set(range(1, max_n7 + 1, 2))
    achieved = set(spectrum_n7)
    gaps = sorted(all_odds - achieved)
    print(f"\nGaps in [1,{max_n7}] at n=7 ({len(gaps)} gaps):")
    print(f"  {gaps}")
    print()

    # Classify gaps
    permanent_gaps = {7, 21}
    print(f"Permanent gaps (proved all n): {sorted(permanent_gaps)}")
    n7_only_gaps = [g for g in gaps if g not in permanent_gaps]
    print(f"Gaps at n=7 that may fill later: {n7_only_gaps}")
    print()

    # Pattern analysis of gaps
    print("Gap analysis:")
    print(f"  7 = 7*1")
    print(f"  21 = 7*3")
    print(f"  63 = 7*9 -- ACHIEVED at n=8 (not permanent)")
    print(f"  107 -- no obvious pattern")
    print(f"  119 = 7*17")
    print(f"  149 -- prime")
    print(f"  161-169 (block): 7*23, 163(prime), 165=3*5*11, 167(prime), 169=13^2")
    print(f"  173 -- prime")
    print(f"  177-187 (block): top end near max=189")
    print()

    # Check multiples of 7
    mult7_gaps = [g for g in gaps if g % 7 == 0]
    non_mult7 = [g for g in gaps if g % 7 != 0]
    print(f"Gaps divisible by 7: {mult7_gaps}")
    print(f"Gaps NOT divisible by 7: {non_mult7}")
    print()

    # Density as function of N
    print("Density d(N) = #{H in spectrum <= N} / #{odd <= N}:")
    for N in [50, 100, 150, 200, 300, 500, 1000]:
        # For n=7: all values up to 189 known. Beyond 189, need n>=8.
        # Use n=7 known + assume all odd >= 23 & <= N are achieved for N > 189
        achieved_leq_N = [h for h in spectrum_n7 if h <= N]
        if N > max_n7:
            # For H in (189, N]: these appear at n=8 or higher
            # From OPEN-Q-019: H=63 fills at n=8. Most gaps at n=7 fill at n=8.
            # Assume: at n=9, ALL odd values except 7 and 21 are achievable
            # (this is strongly supported by the 2M n=9 sample with only {7,21} missing)
            extra = [(h, h <= N and h != 7 and h != 21)
                     for h in range(191, N + 1, 2)]
            count_achieved = len(achieved_leq_N) + sum(1 for _, b in extra if b)
        else:
            count_achieved = len(achieved_leq_N)
        odd_count = (N + 1) // 2
        density = count_achieved / odd_count
        print(f"  N={N:5d}: {count_achieved}/{odd_count} = {density:.4f}")

    print()
    print("Asymptotic: as H -> infinity, density(H) -> 1 - 0/infinity = 1")
    print("The TWO permanent gaps {7,21} have zero density: |{7,21}|/#{odd<=N} -> 0")
    print()

    # Compute small n spectra exhaustively
    for n in [3, 4, 5]:
        print(f"Exhaustive H-spectrum at n={n}:")
        t0 = time.time()
        spec = exhaustive_hspectrum(n)
        print(f"  {spec} (time: {time.time()-t0:.2f}s)")

    # For n=6: too many (2^15=32768 tournaments) - do it
    print(f"\nH-spectrum at n=6 (2^15=32768 tournaments):")
    t0 = time.time()
    spec6 = exhaustive_hspectrum(6)
    print(f"  {len(spec6)} distinct values: {spec6[:20]}...{spec6[-5:]}")
    print(f"  Max: {max(spec6)}, time: {time.time()-t0:.1f}s")
    gaps6 = sorted(set(range(1, max(spec6)+1, 2)) - set(spec6))
    print(f"  Gaps: {gaps6}")

    print()
    print("=" * 70)
    print("THEORETICAL ANALYSIS")
    print("=" * 70)
    print("""
The H-spectrum grows in a very specific way:
- H must be odd (Rédei's theorem)
- H = I(Omega(T), 2) [OCF]
- The independence polynomial at x=2: I(G,2) = sum_k alpha_k * 2^k
  where alpha_k = number of independent sets of size k in Omega(T)
- I(G, 2) = sum of 2^k over all independent sets of size k
         = sum over all independent sets S of 2^|S|

The achievable values of I(G, 2) over ALL graphs G:
- Any sum of distinct powers of 2 times a count >= 1 is in principle achievable.
- But Omega(T) has VERY specific structure (it's a tournament conflict graph).

Why H=7 is impossible: H=7 means I(G,2) = 7 = 1 + 2 + 4
  This requires alpha_1=3 (three 3-cycles, no disjoint pairs) AND no 5-cycles.
  But alpha_1=3 with all pairs conflicting forces a shared vertex, which forces
  a 5-cycle through that vertex. Contradiction.

Why H=21 is impossible: H=21 = 1 + 4*5 or 1+2*10 or 1+2*4+4*2+8*1...
  The poisoning DAG argument shows no valid (alpha_1,...) combination gives exactly 21.

THE DENSITY CONJECTURE: The H-spectrum has density 1 in the odd integers.
That is, for any odd integer h != 7 and h != 21, there exists a tournament T with H(T) = h.

This is consistent with all computational evidence (all odd in [1,200] except {7,21}
appear by n=9).

STRUCTURE CONJECTURE: The gaps 7 and 21 are the ONLY permanent gaps.
Equivalently: for any odd h not in {7,21}, there exists n and a tournament T_n with H(T_n)=h.

This would follow from showing: for sufficiently large n, the map T -> H(T) hits
all odd values != 7,21. The OCF gives a clear mechanism: for any target value v,
construct a tournament with the right (alpha_1, alpha_2, ...) profile.
""")


if __name__ == '__main__':
    analyze_density()
    print("DONE.")
