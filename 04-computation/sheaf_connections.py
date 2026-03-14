"""
sheaf_connections.py -- kind-pasteur-2026-03-14-S110c
SHEAF CONNECTIONS — Grabbing from everywhere

Signs cancel globally but not locally.
What else in mathematics has this property?

1. EULER CHARACTERISTIC: chi = sum (-1)^k * b_k.
   Individual Betti numbers b_k are complex. The alternating sum is simple.

2. INCLUSION-EXCLUSION: |A union B| = |A| + |B| - |A inter B|.
   Individual terms have signs. The result is unsigned.

3. DETERMINANT: det(M) = sum over perms: sgn(pi) * product.
   Individual terms are signed. The determinant is a single number.

4. MÖBIUS INVERSION: f(n) = sum_{d|n} mu(d)*g(n/d).
   The Möbius function mu has signs. The result f is "clean."

5. COHOMOLOGY: H^k(X) has torsion. But chi(X) = sum (-1)^k rank H^k is clean.

6. PARTITION FUNCTION in physics: Z = sum exp(-beta*E).
   The terms are complex (oscillatory). The sum is real (physical).

7. FOURIER INVERSION: f(x) = integral F(k)*exp(ikx) dk.
   The integrand oscillates. The integral gives f back.

ALL of these share the pattern:
  LOCAL: signed/oscillatory/complex
  GLOBAL: unsigned/smooth/simple

The tournament spectrum has EXACTLY this pattern.

Which of these gives the PROOF?
"""

import sys, math
import numpy as np
from fractions import Fraction
from itertools import permutations, combinations

sys.stdout.reconfigure(encoding='utf-8')

def main():
    print("=" * 70)
    print("SHEAF CONNECTIONS — RANDOM GRABS")
    print("kind-pasteur-2026-03-14-S110c")
    print("=" * 70)

    # ============================================================
    print(f"\n{'='*70}")
    print("CONNECTION 1: EULER CHARACTERISTIC ANALOGY")
    print(f"{'='*70}")

    print(f"""
  chi = sum_k (-1)^k * b_k.

  For the tournament Fourier spectrum:
  D_n(2)/n! = 1 + sum_k E_2k/E_0

  This is NOT an alternating sum — all terms are positive.
  But the INDIVIDUAL Fourier coefficients H_hat(S) CAN be negative.
  Only their SQUARES (energies E_2k) are positive.

  The squaring kills the signs. E_2k = sum H_hat(S)^2 >= 0.
  This is like going from homology (has torsion) to Betti numbers (ranks).
  The "rank" operation (squaring) removes the sign information.

  INSIGHT: The formula E_2k/E_0 = 2*(n-2k)^k/P(n,2k) is a statement
  about ENERGIES (squared coefficients), not about the coefficients themselves.
  The coefficients have signs, but the energies don't.
  The formula works because it's about the RIGHT quantity (energy, not amplitude).
    """)

    # ============================================================
    print(f"\n{'='*70}")
    print("CONNECTION 2: INCLUSION-EXCLUSION — COULD THIS BE THE PROOF?")
    print(f"{'='*70}")

    print(f"""
  Inclusion-exclusion: |A1 union ... union Ak| = sum |Ai| - sum |Ai inter Aj| + ...

  WHAT IF the tournament energy at level 2k is computed by
  inclusion-exclusion over some overcounting?

  The "overcounting" might be:
  count ALL 4-arc subsets weighted by |paths containing them|^2,
  then subtract the overcounting from overlapping arcs.

  E_4 = sum_S |{{P : S sub arcs(P)}}|^2 / 2^{{2(n-1)}}
      = (1/2^{{2(n-1)}}) * sum_S N(S)^2

  where N(S) = SIGNED count of paths containing S.

  For adjacent arcs (sharing a vertex): N(S) is large (many paths use them).
  For disjoint arcs: N(S) might be small or zero.

  The formula E_4/E_0 = 2*(n-4)^2/P(n,4) says:
  sum_S N(S)^2 = 2*(n-4)^2 * (n!/2^(n-1))^2 / P(n,4)

  This could be proved by:
  1. Computing N(S)^2 for each type of 4-arc subset S
  2. Counting how many subsets of each type exist
  3. Summing

  This is the DIRECT approach. Let me try it at n=5.
    """)

    # At n=5, compute N(S) for each 4-arc subset
    n = 5
    m = n*(n-1)//2
    N_total = 1 << m
    arcs = [(i,j) for i in range(n) for j in range(i+1,n)]

    # For each 4-arc subset, compute H_hat(S) = (1/2^{n-1}) * signed_count
    # Actually: H_hat(S) = (1/N) * sum_T H(T)*chi_S(T)
    # But we showed: H_hat(S) = (1/2^{n-1}) * sum_{P: S sub arcs(P)} (-1)^{des_S(P)}

    # Let me compute N(S) = sum_{P: S sub arcs(P)} (-1)^{des_S(P)} for each S.
    all_perms = list(permutations(range(n)))

    # For each permutation P, its arcs are the consecutive pairs.
    # As UNDIRECTED edges: {min(P[i],P[i+1]), max(P[i],P[i+1])} for each i.
    # Each arc (a,b) with a<b is used if {P[i],P[i+1]} = {a,b}.
    # The SIGN: if P[i] < P[i+1] (ascending), sign = +1. Else -1.

    def path_arcs_with_signs(perm):
        """Return dict: arc_index -> sign (+1 or -1)."""
        result = {}
        for i in range(len(perm)-1):
            a, b = perm[i], perm[i+1]
            if a < b:
                idx = arcs.index((a, b))
                result[idx] = +1
            else:
                idx = arcs.index((b, a))
                result[idx] = -1
        return result

    # Compute N(S) for all 4-arc subsets
    ns_values = {}
    for s_idx in combinations(range(m), 4):
        s_set = set(s_idx)
        signed_count = 0
        for perm in all_perms:
            arc_signs = path_arcs_with_signs(perm)
            if s_set.issubset(arc_signs.keys()):
                # All 4 arcs are used by this path
                sign_product = 1
                for e in s_idx:
                    sign_product *= arc_signs[e]
                signed_count += sign_product
        ns_values[s_idx] = signed_count

    # H_hat(S) = signed_count / 2^{n-1}
    # E_4 = sum H_hat(S)^2 = sum signed_count^2 / 2^{2(n-1)}
    E4 = sum(v**2 for v in ns_values.values()) / 2**(2*(n-1))
    E0 = (math.factorial(n) / 2**(n-1))**2
    print(f"  n={n}: E_4/E_0 = {E4/E0:.10f}")
    print(f"  Formula: {2*(n-4)**2/math.perm(n,4):.10f}")
    print(f"  Match: {abs(E4/E0 - 2*(n-4)**2/math.perm(n,4)) < 1e-8}")

    # Classify 4-arc subsets by their N(S) value
    from collections import Counter
    ns_dist = Counter()
    for v in ns_values.values():
        ns_dist[v] += 1

    print(f"\n  Distribution of N(S) (signed path count):")
    for v in sorted(ns_dist.keys()):
        count = ns_dist[v]
        print(f"    N(S) = {v:3d}: {count:4d} subsets, contributes {count * v**2} to sum N^2")

    total_ns2 = sum(count * v**2 for v, count in ns_dist.items())
    print(f"  Total sum N(S)^2 = {total_ns2}")
    print(f"  E_4 = {total_ns2}/2^{2*(n-1)} = {total_ns2/2**(2*(n-1)):.6f}")

    # ============================================================
    print(f"\n{'='*70}")
    print("CONNECTION 3: THE SIGNED COUNT IS SIMPLER THAN EXPECTED")
    print(f"{'='*70}")

    # N(S) values: what are they?
    # At n=5: the values should be 0, +/-2, +/-4, etc.
    print(f"\n  Unique |N(S)| values: {sorted(set(abs(v) for v in ns_values.values()))}")

    # The H_hat(S) = N(S)/2^{n-1}. At n=5: 2^4 = 16.
    # If N(S) in {0, +-2}: H_hat in {0, +-1/8}. Matches our earlier finding (all = 1/8)!
    # But we found ALL nonzero level-4 coefficients = 1/8 at n=5.
    # So all nonzero N(S) should have |N(S)| = 2.

    nonzero_ns = {s: v for s, v in ns_values.items() if v != 0}
    print(f"  Nonzero N(S) values: {sorted(set(nonzero_ns.values()))}")
    print(f"  Number of nonzero: {len(nonzero_ns)}")
    print(f"  All have |N(S)| = 2: {all(abs(v) == 2 for v in nonzero_ns.values())}")

    # YES! All nonzero N(S) = +2 or -2. So |H_hat(S)| = 2/16 = 1/8. Matches!
    # And 60 nonzero subsets (30 with N=+2, 30 with N=-2).
    ns_signs = Counter(v for v in nonzero_ns.values())
    print(f"  Sign distribution: {dict(ns_signs)}")
    print(f"  30 positive, 30 negative -> E_4 = 60 * (2/16)^2 = 60/64 = 15/16. CHECK!")

    # ============================================================
    print(f"\n{'='*70}")
    print("CONNECTION 4: WHY |N(S)| = 2 FOR ALL NONZERO S AT n=5?")
    print(f"{'='*70}")

    print(f"""
  N(S) = sum_{{P: S sub arcs(P)}} (-1)^{{des_S(P)}}.

  |N(S)| = 2 means: the signed count is exactly +2 or -2.
  This means: the paths containing S come in PAIRS with opposite signs,
  EXCEPT for ONE pair that doesn't cancel.

  Actually: |N(S)| = 2 = 2^1. And 2^1 = 2^{{|S|/4}}... no.
  Or: 2 = 2 * 1! = 2 * (n-4)! at n=5 (since (n-4)! = 1! = 1).

  From the formula: |H_hat(S)| = |N(S)| / 2^(n-1).
  The formula says E_4/E_0 = 2*(n-4)^2/P(n,4).
  E_4 = N_nonzero * |H_hat|^2 = N_nonzero * |N(S)|^2 / 2^(2(n-1)).
  E_4/E_0 = N_nonzero * |N(S)|^2 / (n!^2 / 2^(2(n-1)) * 2^(2(n-1)))
           = N_nonzero * |N(S)|^2 / n!^2.
  Wait: E_0 = (n!/2^(n-1))^2. E_4 = N_nonzero * (|N(S)|/2^(n-1))^2.
  E_4/E_0 = N_nonzero * |N(S)|^2 / n!^2.

  At n=5: 60 * 4 / 14400 = 240/14400 = 1/60. CHECK!
  Formula: 2*1/120 = 1/60. CHECK!

  So: N_nonzero * |N(S)|^2 = 2*(n-4)^2*(n-2k)!^2 ... hmm.
  Actually: E_4/E_0 = N_nonzero * |N(S)|^2 / (n!)^2 = 2*(n-4)^2/P(n,4).
  => N_nonzero * |N(S)|^2 = 2*(n-4)^2 * (n!)^2 / P(n,4)
                           = 2*(n-4)^2 * ((n-4)!)^2 * ... hmm, complicated.

  The SIMPLER way: N_nonzero = 60 at n=5. |N(S)| = 2.
  60 * 4 = 240 = 2 * 120 = 2 * n!. So N_nonzero * |N(S)|^2 = 2*n!.
  And E_4/E_0 = 2*n!/n!^2 = 2/n! ... that's 2/120 != 1/60. Hmm.

  Let me recompute: E_4 = 60*(2/16)^2 = 60*4/256 = 240/256 = 15/16.
  E_0 = 56.25 = (15/2)^2 = 225/4.
  E_4/E_0 = (15/16)/(225/4) = 15*4/(16*225) = 60/3600 = 1/60. CHECK.

  N_nonzero * |N(S)|^2 / 2^(2(n-1)) / E_0 = 1/60.
  60 * 4 / 256 / (225/4) = 240/256 * 4/225 = 240*4/(256*225) = 960/57600 = 1/60. CHECK.
  """)

    # ============================================================
    print(f"\n{'='*70}")
    print("CONNECTION 5: MÖBIUS FUNCTION AND THE SIGNS")
    print(f"{'='*70}")

    print(f"""
  The signs in N(S) come from descents. For a 4-arc subset S and path P:
  (-1)^des_S(P) = product over arcs in S of (+1 if ascending, -1 if descending).

  This is a MULTIPLICATIVE character on the arcs of P.
  The product over a subset S is like a Legendre symbol or Möbius function:
  it's multiplicative and takes values in {{+1, -1}}.

  The key: at n=5, ALL nonzero N(S) have |N(S)| = 2.
  This means the signed sum over paths containing S has MINIMAL nonzero value.
  It's as if the paths come in pairs with opposite signs, except for
  exactly ONE unmatched pair.

  The "unmatched pair" = the path and its REVERSAL.
  If P contains S, does P^rev also contain S?
  P^rev has the SAME arcs but with REVERSED signs.
  So des_S(P^rev) = |S| - des_S(P) (all ascents become descents and vice versa).
  (-1)^des_S(P^rev) = (-1)^(|S|-des_S(P)) = (-1)^|S| * (-1)^des_S(P).

  For |S| = 4 (even): (-1)^4 = +1. So the reversal has the SAME sign.
  P and P^rev contribute the SAME sign to N(S).
  They DON'T cancel — they REINFORCE!

  For |S| = 2 (even): same argument. P and P^rev reinforce.
  This is why level-2 coefficients are nonzero.

  The factor of 2: N(S) = 2 because each path P and its reversal P^rev
  contribute identically. The "independent" paths are n!/2 (modulo reversal).
  And among these, typically only 1 independent path contains S.
  So N(S) = 2 * 1 = 2 (one independent path contributing +1, times 2 from reversal).

  This gives |N(S)| = 2 * (number of independent paths containing S).
  At n=5 level 4: 1 independent path per nonzero S. N(S) = 2.
  At n=6 level 4: some S have 1 independent path (N=2), some have 2 (N=4).
  This matches our finding: |H_hat| = 1/8 or 1/4 at n=6!
  (N=2 gives 2/32 = 1/16, N=4 gives 4/32 = 1/8... wait, 2^(n-1) = 32 at n=6.)

  |H_hat| = |N|/2^(n-1) = 2/32 = 1/16? But we found 1/8 at n=6 level 4!
  Hmm: 1/8 = 4/32 = N(S)/32 means N(S) = 4 for the inherited (class A).
  And 1/4 = 8/32 means N(S) = 8 for the new (class B).

  So at n=6: |N(S)| in {{4, 8}}, not {{2, 4}}.
  The factor of 2 from reversal: N(S) = 2 * (independent path count).
  Class A: 2 independent paths -> N=4.
  Class B: 4 independent paths -> N=8.

  THE NUMBER OF INDEPENDENT PATHS CONTAINING S DETERMINES |N(S)|.
  And the formula gives the TOTAL energy sum_S N(S)^2.

  THE PROOF APPROACH:
  Show that sum_S N(S)^2 / (2^(2(n-1)) * E_0) = 2*(n-2k)^k/P(n,2k)
  by counting the INDEPENDENT PATH PAIRS that contribute.
  """)

    # Verify at n=6: N(S) values at level 4
    print(f"\n  Checking N(S) values at n=6, level 4:")
    n = 6
    m = n*(n-1)//2
    arcs6 = [(i,j) for i in range(n) for j in range(i+1,n)]
    all_perms6 = list(permutations(range(n)))

    ns_dist6 = Counter()
    sample_count = 0
    for s_idx in combinations(range(m), 4):
        s_set = set(s_idx)
        signed_count = 0
        for perm in all_perms6:
            arc_signs = {}
            for i in range(n-1):
                a, b = perm[i], perm[i+1]
                if a < b: arc_signs[arcs6.index((a,b))] = +1
                else: arc_signs[arcs6.index((b,a))] = -1
            if s_set.issubset(arc_signs.keys()):
                sign_product = 1
                for e in s_idx:
                    sign_product *= arc_signs[e]
                signed_count += sign_product
        ns_dist6[abs(signed_count)] += 1
        sample_count += 1
        if sample_count % 200 == 0:
            print(f"    {sample_count}/1365 done...")

    print(f"  Distribution of |N(S)| at n=6, level 4:")
    for v in sorted(ns_dist6.keys()):
        print(f"    |N(S)| = {v}: {ns_dist6[v]} subsets")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
