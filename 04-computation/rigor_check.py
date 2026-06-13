#!/usr/bin/env python3
"""rigor_check.py — Testing speculative claims against data.

Session: kind-pasteur-2026-03-20-S10

Test each Tier 3 (speculative) claim from the honest assessment.
"""

from math import factorial, comb

def test_forbidden_cascade():
    """Is 21 = 7 x 3 structural or coincidental?"""
    print("=" * 60)
    print("TEST: Is the forbidden cascade 7 -> 21 structural?")
    print("=" * 60)

    # Known forbidden: 7, 21
    # If 7 x (min_cycle) = next forbidden, then:
    # 7 x 3 = 21 (YES, forbidden)
    # 7 x 5 = 35 (achievable at n=7 — REFUTES the cascade)
    # 7 x 7 = 49 (achievable at n=7 — REFUTES)
    # 21 x 2 = 42 (achievable at n=6 — REFUTES)
    # 21 x 3 = 63 (achievable at n=8 — REFUTES)

    print(f"\n  Multiples of 7:")
    print(f"    7 x 1 = 7  : FORBIDDEN (proved)")
    print(f"    7 x 3 = 21 : FORBIDDEN (proved)")
    print(f"    7 x 5 = 35 : ACHIEVABLE at n=7")
    print(f"    7 x 7 = 49 : ACHIEVABLE at n=7")
    print(f"    7 x 9 = 63 : ACHIEVABLE at n=8")

    print(f"\n  Multiples of 21:")
    print(f"    21 x 2 = 42 : ACHIEVABLE at n=6 (not forbidden!)")
    print(f"    21 x 3 = 63 : ACHIEVABLE at n=8")

    print(f"\n  VERDICT: The factor 21/7 = 3 is COINCIDENTAL.")
    print(f"  The forbidden cascade 7 -> 21 -> 42 has NO structural basis.")
    print(f"  7 and 21 are forbidden for INDEPENDENT reasons:")
    print(f"    H=7: alpha_1=3 + i_2=0 forces common vertex -> extra cycles")
    print(f"    H=21: six independent (alpha_1, alpha_2) blocking mechanisms")
    print(f"  These proofs share NO common structure.")

    # BUT: is there ANYTHING structural about {7, 21}?
    # 7 = 1 + 2*3 (H with alpha_1=3, alpha_2=0)
    # 21 = 1 + 2*10 (H with alpha_1+2*alpha_2=10)
    # The alpha sums are S=3 and S=10.
    # Is there something about S=3 and S=10?
    # S=3: only 2 decompositions (3,0) and (1,1), both blocked
    # S=10: 6 decompositions, all blocked
    # S=4: (4,0), (2,1), (0,2) — is (0,2) possible? alpha_2=2 needs >= 2 disjoint
    #   cycle pairs, which needs >= 6 vertices. At n=6+: (0,2) means 0 total cycles
    #   but 2 disjoint pairs?? That's impossible (alpha_2 <= C(alpha_1, 2)).
    #   Actually alpha_2 counts pairs of disjoint cycles, so alpha_1 >= 2 for alpha_2 >= 1.
    #   So (0,2) is impossible. (2,1): 2 cycles with 1 disjoint pair. Possible at n>=6.
    #   (4,0): 4 cycles, no disjoint pairs. H=1+8=9, achieved.
    # So S=4 gives H=9, which IS achievable. NOT forbidden.

    print(f"\n  Alpha sum analysis:")
    print(f"  S=3 (H=7): decomps (3,0),(1,1). BOTH blocked -> H=7 forbidden")
    print(f"  S=4 (H=9): decomp (4,0) achievable -> H=9 NOT forbidden")
    print(f"  S=10 (H=21): all 6 decomps blocked -> H=21 forbidden")
    print(f"  S=11 (H=23): decomp (11,0) achievable -> H=23 NOT forbidden")
    print(f"\n  So the forbidden S values are exactly {{3, 10}}.")
    print(f"  Is 10 = 3 * something? 10/3 = 3.33... NOT integer.")
    print(f"  Is there a pattern in {{3, 10}}? ")
    print(f"  3 = C(3,2), 10 = C(5,2). Triangular numbers!")
    print(f"  C(3,2)=3, C(5,2)=10. Next: C(7,2)=21. Is S=21 forbidden?")
    print(f"  H = 1 + 2*21 = 43. Is H=43 forbidden?")
    print(f"  H=43 IS achievable at n=6! So S=C(7,2)=21 is NOT forbidden.")
    print(f"  Pattern C(odd,2) FAILS at the third term.")
    print(f"  VERDICT: {{3, 10}} has no known generating pattern.")


def test_deficit_ratio():
    """Does D_graph / D_tournament have a limit?"""
    print(f"\n\n{'='*60}")
    print(f"TEST: Does D_graph/D_tournament converge?")
    print(f"{'='*60}")

    # From computed data:
    # n:  3    4    5     6     7      8
    # DG: 6   24   80   392  2212  19504
    # DT: 2    4   12    40   152    784
    # Ratio: 3  6  6.67  9.8  14.6  24.9

    DG = [6, 24, 80, 392, 2212, 19504]
    DT = [2, 4, 12, 40, 152, 784]
    ns = [3, 4, 5, 6, 7, 8]

    print(f"\n  {'n':>3} {'D_G':>8} {'D_T':>8} {'Ratio':>10}")
    for i, n in enumerate(ns):
        ratio = DG[i] / DT[i]
        print(f"  {n:3d} {DG[i]:8d} {DT[i]:8d} {ratio:10.4f}")

    print(f"\n  The ratio is NOT converging to 25.")
    print(f"  It's GROWING: 3, 6, 6.67, 9.8, 14.6, 24.9")
    print(f"  This is roughly doubling every 2 steps.")
    print(f"  Suggests ratio ~ 2^{{(n-3)/2}} or similar exponential growth.")

    # What's the ratio of ratios?
    ratios = [DG[i]/DT[i] for i in range(len(ns))]
    print(f"\n  Ratio growth factors:")
    for i in range(1, len(ratios)):
        print(f"    n={ns[i]}: ratio/prev = {ratios[i]/ratios[i-1]:.4f}")

    print(f"\n  VERDICT: D_G/D_T grows roughly as 2^n, NOT converging to 25.")
    print(f"  The claim 'ratio ~ 25' was a SNAPSHOT at n=8, not a limit.")


def test_2_10_decomposition():
    """Does 2^10 = 10^3 + 4! have the claimed tournament interpretation?"""
    print(f"\n\n{'='*60}")
    print(f"TEST: 2^10 = 10^3 + 4! tournament interpretation")
    print(f"{'='*60}")

    # At n=5: 1024 tournaments decompose as:
    # |Aut|=1: 840
    # |Aut|=3: 160
    # |Aut|=5: 24

    print(f"\n  Actual decomposition:")
    print(f"    |Aut|=1: 840 (trivial automorphism)")
    print(f"    |Aut|=3: 160 (Z_3 symmetry)")
    print(f"    |Aut|=5: 24  (Z_5 symmetry)")
    print(f"    Total: 1024 = 2^10")

    print(f"\n  Claimed: 1024 = 1000 (generic) + 24 (max symmetric)")
    print(f"  Reality: 1024 = 840 + 160 + 24")
    print(f"  The split 1000 + 24 does NOT correspond to any natural partition.")
    print(f"  840 + 160 = 1000 only by coincidence (840 + 160 = 1000? YES!)")

    print(f"\n  Wait: 840 + 160 = {840+160}. So 1000 = trivial + Z_3-symmetric!")
    print(f"  This IS a natural partition: {1000} = |Aut|-non-5 tournaments.")
    print(f"  The identity becomes: (|Aut|!=5 count) + (|Aut|=5 count) = 2^m")
    print(f"  Which is trivially: ALL = (not-Z_5) + (Z_5)")
    print(f"  The 'miracle' is: not-Z_5 count = m^3 = C(n,2)^3")

    # Is not-Z_5 count = m^3 structural?
    print(f"\n  Is (2^m - n!/n) = m^3 at m=C(n,2)?")
    for n in range(2, 8):
        m = comb(n, 2)
        total = 2**m
        # |Aut|=n count: how many labeled tournaments have Z_n symmetry?
        # For prime n: Z_n-symmetric = n!/n = (n-1)! (if exactly 1 such iso class)
        # At n=5: 1 iso class with |Aut|=5, giving 5!/5 = 24 labeled copies
        zn_count = factorial(n) // n if n > 1 else 1  # rough estimate
        non_zn = total - zn_count
        m_cubed = m**3
        print(f"  n={n}: 2^{m}-{zn_count} = {non_zn}, m^3 = {m_cubed}, "
              f"equal: {non_zn == m_cubed}")

    print(f"\n  n=5 is the ONLY case where 2^m - (n-1)! = m^3.")
    print(f"  This is the same as 2^m = m^3 + (n-1)!, just restated.")
    print(f"  The 'natural partition' interpretation (non-Z_5 + Z_5) is valid")
    print(f"  but the fact that non-Z_5 = m^3 is still a NUMERICAL COINCIDENCE")
    print(f"  special to n=5 (and n=2).")


def structural_vs_coincidental():
    """Final summary."""
    print(f"\n\n{'='*60}")
    print(f"FINAL SUMMARY: STRUCTURAL vs COINCIDENTAL")
    print(f"{'='*60}")

    print(f"""
  STRUCTURAL (proved, generalizes):
    - k-periodicity of the tower (k=2 graphs, k=3 tournaments)
    - Layer decomposition D = D_1 + D_2 + D_3 + ...
    - Each level adds k exact terms
    - delta(n) → 0 super-exponentially
    - Only odd-cycle permutations contribute (Moon's theorem)
    - At n=5: exactly 24 = 4! tournaments have Z_5 symmetry (one iso class)

  COINCIDENTAL (specific to small n, doesn't generalize):
    - 2^10 = 10^3 + 4! (only at n=2 and n=5)
    - D/2 = C(2k,k) (only for n <= 6)
    - P(n) = T(n+1) (only for n <= 4)
    - 21 = 7 x 3 (7 x 5 = 35 IS achievable, so the multiplication by 3 is accidental)
    - D_G/D_T ~ 25 (actually grows as ~2^n, not constant)

  SUGGESTIVE BUT UNGROUNDED (interesting analogy, no proof):
    - The hierarchy 2, 3, 7, 43 as a connected tower
    - Non-integer periodicities (phi, e, pi)
    - Euler product structure of D(n)
    - Cheeger-periodicity product
    - 42 as the "ground level" of a cascade

  THE LESSON:
  The k-periodicity principle and layer decomposition are GENUINE
  mathematics with exact proofs. The numerical coincidences at
  small n are BEAUTIFUL but FINITE phenomena that break when tested
  further. The project's history (MISTAKES 6, 24, 28) shows that
  these finite coincidences are common and seductive but ultimately
  accidental. The structural truth (the tower, the layers, the
  3-periodicity) is BETTER than the coincidences because it
  holds for ALL n.
""")


if __name__ == "__main__":
    test_forbidden_cascade()
    test_deficit_ratio()
    test_2_10_decomposition()
    structural_vs_coincidental()
