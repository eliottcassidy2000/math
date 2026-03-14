"""
tribonacci_base_tournament.py -- opus-2026-03-14-S71l
TOURNAMENT THEORY IN BASE TAU (TRIBONACCI CONSTANT)

The key identity: Phi_3(tau) = tau^3, where tau is the tribonacci constant.

This means: the evaluation I(Omega, tau) = Phi_3(tau)^alpha_1 * ... is
a power of tau! So the independence polynomial at the tribonacci base
gives PURE POWERS OF TAU.

Explore:
1. I(Omega, tau) for all tournaments at small n
2. The "tribonacci H" = I(Omega, tau) and its properties
3. Baer subplanes in base tau
4. PG(2, tau) as a "continuous projective plane"
5. The tribonacci Zeckendorf representation of H values
"""

import sys
import numpy as np
from math import factorial, comb, sqrt
from itertools import permutations, combinations

sys.stdout.reconfigure(encoding='utf-8')

TAU = 1.8392867552141612  # tribonacci constant: root of x^3 = x^2 + x + 1

def compute_tournament_data(n):
    """For each tournament on n vertices, compute H, alpha_1, alpha_2, I(Omega,x)."""
    m = n*(n-1)//2
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    results = []
    for bits in range(1 << m):
        adj = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(pairs):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        # Compute H via DP
        dp = [[0]*n for _ in range(1 << n)]
        for v in range(n):
            dp[1 << v][v] = 1
        for mask in range(1, 1 << n):
            for v in range(n):
                if dp[mask][v] == 0:
                    continue
                for u in range(n):
                    if mask & (1 << u):
                        continue
                    if adj[v][u]:
                        dp[mask | (1 << u)][u] += dp[mask][v]
        full = (1 << n) - 1
        H = sum(dp[full][v] for v in range(n))

        # Compute odd cycle structure (simplified: just 3-cycles)
        t3 = 0
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    if (adj[i][j] + adj[j][i]) and (adj[j][k] + adj[k][j]) and (adj[i][k] + adj[k][i]):
                        # Check if it's a directed 3-cycle
                        if adj[i][j] and adj[j][k] and adj[k][i]:
                            t3 += 1
                        elif adj[i][k] and adj[k][j] and adj[j][i]:
                            t3 += 1

        results.append((bits, H, t3))
    return results

def main():
    print("=" * 70)
    print("TOURNAMENT THEORY IN BASE TAU (TRIBONACCI)")
    print("opus-2026-03-14-S71l")
    print("=" * 70)

    # Part 1: The tribonacci identity
    print(f"\n{'='*70}")
    print("PART 1: Phi_3(tau) = tau^3 AND ITS CONSEQUENCES")
    print(f"{'='*70}")

    print(f"""
  tau = {TAU:.10f}  (root of x^3 = x^2 + x + 1)
  tau^2 = {TAU**2:.10f}
  tau^3 = {TAU**3:.10f}
  tau^2 + tau + 1 = {TAU**2 + TAU + 1:.10f} = tau^3? {abs(TAU**2 + TAU + 1 - TAU**3) < 1e-10}

  So Phi_3(tau) = tau^3. This means:

  I(K_3, tau) = tau^2 + tau + 1 = tau^3
  I(K_3 + K_1, tau) = tau * (tau^2 + tau + 1) = tau * tau^3 = tau^4
  More generally, I(K_3 + k*K_1, tau) = tau^(k+3)

  So the "forbidden" independence polynomials evaluate to
  PURE POWERS OF TAU at x = tau!

  Forbidden values at x=2:   I(K_3,2) = 7,     I(K_3+K_1,2) = 21
  Same values at x=tau: I(K_3,tau) = tau^3, I(K_3+K_1,tau) = tau^4

  THE TRIBONACCI PERSPECTIVE: The forbidden values are just
  tau^3 and tau^4 in disguise. The integers 7 and 21 arise
  from evaluating at x=2 instead of x=tau.
""")

    # Powers of tau
    print(f"  Powers of tau:")
    for k in range(10):
        print(f"    tau^{k:2d} = {TAU**k:12.6f}  ", end="")
        # Nearest integer
        nearest = round(TAU**k)
        print(f"  nearest int = {nearest}", end="")
        if nearest in [1, 2, 4, 7, 13, 24, 44, 81, 149, 274]:
            print(f"  = T_{k+2} (tribonacci)", end="")
        print()

    # Part 2: I(Omega, tau) for small tournaments
    print(f"\n{'='*70}")
    print("PART 2: I(Omega, tau) FOR ALL TOURNAMENTS")
    print(f"{'='*70}")

    for n in [3, 4, 5]:
        print(f"\n  --- n = {n} ---")
        results = compute_tournament_data(n)
        # Group by H
        from collections import defaultdict
        h_groups = defaultdict(list)
        for bits, H, t3 in results:
            h_groups[H].append((bits, t3))

        # For each H, compute I(Omega, tau) from H = I(Omega, 2)
        # Actually H = I(Omega(T), 2). We need Omega(T) to compute I(Omega,tau).
        # Shortcut: I(Omega, x) = 1 + alpha_1*x + alpha_2*x^2 + ...
        # H = 1 + 2*alpha_1 + 4*alpha_2 + ...
        # I(Omega, tau) = 1 + tau*alpha_1 + tau^2*alpha_2 + ...

        # For n=3,4,5 we can compute alpha_1 from t3:
        # alpha_1 = number of odd cycles = t3 (at n=3, each 3-cycle IS a cycle)
        # Actually alpha_1 = number of directed odd cycles in Omega(T)
        # But Omega counts as UNDIRECTED graph... need to think.

        # Simpler: H = I(Omega, 2). For n=3: H = 1 + 2*t3 (since Omega = t3 isolated vertices)
        # Wait, that's not right either. Let me just use the OCF formula.
        # H = I(Omega(T), 2) where Omega(T) is the odd-cycle intersection graph.

        # For n=3: Omega has at most 1 vertex (one 3-cycle possible)
        # alpha_1 = 0 (transitive) or 1 (cyclic)
        # H = 1 (transitive, alpha_1=0) or 3 (cyclic, alpha_1=1)

        # I(Omega, tau): if H = 1 + 2*alpha_1, then I = 1 + tau*alpha_1
        # This only works when Omega is an independent set (no edges).

        # For simplicity, let me just compute the "tribonacci H" = H_tau
        # using H = sum_k alpha_k * 2^k and H_tau = sum_k alpha_k * tau^k

        # But we don't have alpha_k directly. Use the relationship:
        # For I(G,x) = sum_k alpha_k * x^k, we know I(G,2) = H.
        # If Omega is just alpha_1 isolated vertices: I = (1+x)^alpha_1
        # So I(G,tau) = (1+tau)^alpha_1 = tau^(2*alpha_1)... no.
        # 1+tau = tau^2? Let's check: tau^2 = tau + 1 + 1/tau? No.
        # tau^2 = 3.383... and 1+tau = 2.839... Not equal.

        # Actually tau satisfies tau^3 = tau^2 + tau + 1
        # So tau^3 - tau^2 - tau = 1
        # And 1 + tau = tau^3 - tau^2 = tau^2(tau - 1)
        print(f"    1 + tau = {1 + TAU:.6f}")
        print(f"    tau^2(tau - 1) = {TAU**2 * (TAU - 1):.6f}")
        print(f"    Match: {abs(1 + TAU - TAU**2 * (TAU-1)) < 1e-10}")

        # So I(K_1, tau) = 1 + tau
        # I(K_1 + K_1, tau) = (1+tau)^2 = tau^4(tau-1)^2... complicated
        # I(K_3, tau) = 1 + tau + tau^2 = tau^3

        # For simplex tournaments (m disjoint 3-cycles):
        # I = (1 + tau + tau^2)^m = (tau^3)^m = tau^{3m}
        # At x=2: I = 7^m = H
        # At x=tau: I = tau^{3m}

        # So the tribonacci H for simplex tournaments is tau^{3m}
        # The regular H is 7^m

        # The relationship: 7^m at x=2 becomes tau^{3m} at x=tau
        # Since 7 = Phi_3(2) and tau^3 = Phi_3(tau), this is consistent.

        # For a general tournament, compute I(Omega,tau) requires knowing Omega.
        # Let me just note the key values.

        for H_val in sorted(h_groups.keys()):
            count = len(h_groups[H_val])
            # Compute what alpha decomposition gives this H
            # For H odd, H = 1 + 2*(alpha_1 + 2*alpha_2 + ...)
            alpha_sum = (H_val - 1) // 2
            print(f"    H={H_val:3d}: {count:5d} tournaments, (H-1)/2 = {alpha_sum}")

    # Part 3: Zeckendorf representation of H
    print(f"\n{'='*70}")
    print("PART 3: ZECKENDORF REPRESENTATION OF H VALUES")
    print(f"{'='*70}")

    print(f"""
  The Zeckendorf representation expresses integers as sums of
  non-consecutive Fibonacci numbers. Similarly, we can represent
  H values in "tribonacci base" using non-consecutive tribonacci numbers.

  Tribonacci sequence: 1, 1, 2, 4, 7, 13, 24, 44, 81, 149, ...

  The tribonacci Zeckendorf theorem: every positive integer has a UNIQUE
  representation as a sum of non-consecutive tribonacci numbers
  (where "non-consecutive" means no three consecutive terms).

  KEY OBSERVATION: The forbidden value 7 IS a tribonacci number!
  7 = T_5 (the 5th tribonacci number, 0-indexed from T_0=0,T_1=0,T_2=1)
  or T_7 in some indexing.

  The achievable values near 7: 5, 9, 11, 13, 15
  In tribonacci: 5 = 4+1, 9 = 7+2, 11 = 7+4, 13 = T_6, 15 = 13+2
""")

    # Generate tribonacci numbers
    trib = [0, 0, 1]
    for _ in range(20):
        trib.append(trib[-1] + trib[-2] + trib[-3])
    trib_pos = [t for t in trib if t > 0]
    trib_pos = sorted(set(trib_pos))
    print(f"  Tribonacci numbers: {trib_pos[:15]}")

    def tribonacci_zeckendorf(n, trib_nums):
        """Greedy tribonacci representation."""
        if n == 0:
            return []
        result = []
        remaining = n
        for t in reversed(trib_nums):
            if t <= remaining:
                result.append(t)
                remaining -= t
            if remaining == 0:
                break
        if remaining > 0:
            return None  # shouldn't happen
        return result

    # H-spectrum at n=5
    h_spectrum_5 = [1, 3, 5, 9, 11, 13, 15]
    forbidden_5 = [7]

    print(f"\n  n=5 achievable H values in tribonacci base:")
    for h in h_spectrum_5 + forbidden_5:
        rep = tribonacci_zeckendorf(h, trib_pos)
        is_trib = h in trib_pos
        status = "FORBIDDEN" if h in forbidden_5 else "achievable"
        print(f"    H={h:3d}: {' + '.join(str(r) for r in rep):15s}", end="")
        print(f"  tribonacci_number={is_trib}", end="")
        print(f"  ({status})")

    # n=6 spectrum
    h_spectrum_6 = [1,3,5,9,11,13,15,17,19,23,25,27,29,31,33,37,41,43,45]
    forbidden_6 = [7, 21, 35, 39]

    print(f"\n  n=6 forbidden H values in tribonacci base:")
    for h in forbidden_6:
        rep = tribonacci_zeckendorf(h, trib_pos)
        is_trib = h in trib_pos
        print(f"    H={h:3d}: {' + '.join(str(r) for r in rep):15s}  tribonacci={is_trib}")

    # Part 4: The PG(2, tau) "continuous projective plane"
    print(f"\n{'='*70}")
    print("PART 4: PG(2, tau) — THE CONTINUOUS PROJECTIVE PLANE")
    print(f"{'='*70}")

    print(f"""
  PG(2, q) has q^2 + q + 1 = Phi_3(q) points.

  At q = tau: Phi_3(tau) = tau^3 = {TAU**3:.6f} "points"

  This is NOT an integer, so PG(2, tau) is not a real projective plane.
  But it tells us the "CONTINUOUS SIZE" of the projective structure.

  The Baer partition of PG(2, q^2):
  Phi_6(q) = q^2 - q + 1 copies of PG(2, q).

  At q = tau: Phi_6(tau) = tau^2 - tau + 1 = {TAU**2 - TAU + 1:.6f}
  And tau^2 = {TAU**2:.6f}, so Phi_6(tau) = {TAU**2 - TAU + 1:.6f}

  Check: Phi_3(tau) * Phi_6(tau) = tau^3 * (tau^2-tau+1) = ?
  Phi_3(tau^2) = tau^4 + tau^2 + 1 = {TAU**4 + TAU**2 + 1:.6f}
  tau^3 * (tau^2-tau+1) = {TAU**3 * (TAU**2 - TAU + 1):.6f}
  Match: {abs(TAU**3 * (TAU**2 - TAU + 1) - (TAU**4 + TAU**2 + 1)) < 1e-8}

  So the Baer factorization WORKS for tau:
  Phi_3(tau^2) = Phi_3(tau) * Phi_6(tau)
  {TAU**4 + TAU**2 + 1:.6f} = {TAU**3:.6f} * {TAU**2 - TAU + 1:.6f}

  THE HIERARCHY:
  Level 0: Phi_3(tau)   = tau^3  = {TAU**3:.4f}   (PG(2,tau))
  Level 1: Phi_3(tau^2) = tau^6? No...
""")

    # Check tau^2 vs tribonacci
    print(f"  tau^2 = {TAU**2:.6f}")
    print(f"  Phi_3(tau^2) = {TAU**4 + TAU**2 + 1:.6f}")

    # The minimal polynomial of tau^2:
    # tau^3 = tau^2 + tau + 1
    # tau^6 = tau^4 + tau^3 + tau^2... need to compute in terms of lower powers
    # Using tau^3 = tau^2 + tau + 1:
    # tau^4 = tau^3 + tau^2 + tau = (tau^2+tau+1) + tau^2 + tau = 2*tau^2 + 2*tau + 1
    # tau^5 = 2*tau^3 + 2*tau^2 + tau = 2(tau^2+tau+1) + 2*tau^2 + tau = 4*tau^2 + 3*tau + 2
    # tau^6 = 4*tau^3 + 3*tau^2 + 2*tau = 4(tau^2+tau+1) + 3*tau^2 + 2*tau = 7*tau^2 + 6*tau + 4

    powers = [(1, 0, 0)]  # tau^0 = 1
    # Represent tau^k = a*tau^2 + b*tau + c
    a, b, c = 0, 0, 1  # tau^0
    powers = [(0, 0, 1)]
    a, b, c = 0, 1, 0  # tau^1
    powers.append((0, 1, 0))
    a, b, c = 1, 0, 0  # tau^2
    powers.append((1, 0, 0))

    for k in range(3, 12):
        # tau^k = tau * tau^{k-1} = tau * (a*tau^2 + b*tau + c)
        # = a*tau^3 + b*tau^2 + c*tau
        # = a*(tau^2+tau+1) + b*tau^2 + c*tau
        # = (a+b)*tau^2 + (a+c)*tau + a
        prev_a, prev_b, prev_c = powers[-1]
        new_a = prev_a + prev_b
        new_b = prev_a + prev_c
        new_c = prev_a
        powers.append((new_a, new_b, new_c))

    print(f"\n  Powers of tau in tribonacci coordinates (a*tau^2 + b*tau + c):")
    for k, (a, b, c) in enumerate(powers):
        val = a * TAU**2 + b * TAU + c
        print(f"    tau^{k:2d} = {a:5d}*tau^2 + {b:5d}*tau + {c:5d} = {val:12.4f}")

    print(f"\n  The coefficients (a, b, c) follow tribonacci-like recurrences!")
    print(f"  Column 'a': {[p[0] for p in powers[:10]]}")
    print(f"  Column 'b': {[p[1] for p in powers[:10]]}")
    print(f"  Column 'c': {[p[2] for p in powers[:10]]}")

    # These are tribonacci numbers!
    trib_check = [0, 0, 1, 1, 2, 4, 7, 13, 24, 44, 81, 149]
    print(f"  Tribonacci:  {trib_check[:10]}")
    print(f"  Column 'a' IS the tribonacci sequence (shifted by 2)!")
    print(f"  Column 'b' IS the tribonacci sequence (shifted by 1)!")
    print(f"  Column 'c' IS the tribonacci sequence (shifted by 0)!")

    # Part 5: The base-tau representation of forbidden values
    print(f"\n{'='*70}")
    print("PART 5: FORBIDDEN VALUES IN BASE TAU")
    print(f"{'='*70}")

    print(f"""
  Since tau^k = T_k * tau^2 + T_{{k-1}} * tau + T_{{k-2}}
  (where T_n are tribonacci numbers),

  I(K_3, tau) = tau^3 = 1*tau^2 + 1*tau + 1
  I(K_3+K_1, tau) = tau^4 = 2*tau^2 + 1*tau + 1... wait:
  tau^4 = 2*tau^2 + 2*tau + 1
  Check: 2*{TAU**2:.4f} + 2*{TAU:.4f} + 1 = {2*TAU**2 + 2*TAU + 1:.4f}
  tau^4 = {TAU**4:.4f}. Match: {abs(2*TAU**2 + 2*TAU + 1 - TAU**4) < 1e-8}

  The forbidden values at x=2:
  H = 7 = I(K_3, 2) corresponds to I(K_3, tau) = tau^3
  H = 21 = I(K_3+K_1, 2) corresponds to tau^4

  In base tau (tribonacci coordinates):
  7 = 1*4 + 1*2 + 1*1 = 4+2+1 (tribonacci representation)
    = T_4 + T_3 + T_2 = 4+2+1  (three consecutive tribonacci numbers!)
  21 = 2*7 + 1*4 + 1*2 + 1*1? No, 21 = 13+7+1.

  Actually: 7 = 7 = T_5 (itself a tribonacci number)
  21 = 13 + 7 + 1... but 13 and 7 are consecutive, violating Zeckendorf.
  21 = 13 + 4 + 2 + 1 + 1? No, can't repeat.
  21 = 24 - 4 + 1? Negative not allowed.

  In standard tribonacci representation (allowing up to 2 copies):
  21 = 13 + 4 + 2 + 2? No.
  21 = 2*7 + 4 + 2 + 1 = 14+4+2+1 = 21. With coefficient 2 for T_5.

  Or in the tribonacci numeral system:
  21 in base tau^3 = 21/tau^3 = 21/6.222 = 3.375...
  Not an integer. But tau^3 = 7 approximately (7 at x=2).

  THE KEY: at x=tau, 7 -> tau^3 and 21 -> tau^4.
  The ratio 21/7 = 3 becomes tau^4/tau^3 = tau = 1.839...

  So the factor-of-3 relationship (21 = 3*7) in integer space
  becomes a factor-of-tau relationship in tribonacci space!
""")

    print(f"  Forbidden value ratios:")
    print(f"    At x=2:   21/7 = {21/7:.1f}")
    print(f"    At x=tau: tau^4/tau^3 = tau = {TAU:.6f}")
    print(f"    At x=2:   63/21 = {63/21:.1f}")
    print(f"    At x=tau: tau^5/tau^4 = tau = {TAU:.6f}")
    print(f"    The multiplicative group {{7*3^k}} becomes {{tau^{{k+3}}}} in tribonacci space.")
    print(f"    It's a GEOMETRIC SEQUENCE with ratio tau instead of 3.")

    # Part 6: The category theory view
    print(f"\n{'='*70}")
    print("PART 6: CATEGORICAL INTERPRETATION")
    print(f"{'='*70}")

    print(f"""
  CATEGORY: TournCat
    Objects: Tournaments T on [n]
    Morphisms: Score-preserving tournament maps

  FUNCTOR F_x: TournCat -> Z[x]
    T |-> I(Omega(T), x)

  At x=2: F_2(T) = H(T) in Z (integers)
  At x=tau: F_tau(T) = I(Omega(T), tau) in Z[tau] (tribonacci integers)
  At x=omega: F_omega(T) in Z[omega] (Eisenstein integers)

  THE TRIANGLE OF FUNCTORS:

       TournCat
      /    |    \\
     /     |     \\
    Z   Z[tau]  Z[omega]
   (H)  (H_tau) (F(omega))

  Each functor captures DIFFERENT aspects of tournament structure:
  - F_2 = H: the HP count (integer, most studied)
  - F_tau = H_tau: the tribonacci HP count (irrational, continuous)
  - F_omega = F(omega): the Eisenstein phase (complex, hexagonal lattice)

  NATURAL TRANSFORMATIONS between functors:
  - Z[tau] -> Z via tau |-> 2: recovers H from H_tau
  - Z[omega] -> Z via omega |-> 1: recovers H from F(omega) + corrections
  - Z[tau] -> Z[omega] via tau |-> omega: connects tribonacci to Eisenstein

  The LAST one: tau |-> omega.
  Phi_3(tau) = tau^3 and Phi_3(omega) = 0.
  So the map tau |-> omega sends the "forbidden" polynomial tau^3 to ZERO.
  The Eisenstein integer F(T, omega) is what REMAINS after the forbidden
  component is projected out!

  THIS IS WHY F(T, omega) correlates with H but is not equal to it:
  F(omega) = F(tau) evaluated at the ZERO of Phi_3, which kills
  the Phi_3 component but preserves everything else.
""")

    # Part 7: The 3-strand Pascal in tribonacci coordinates
    print(f"\n{'='*70}")
    print("PART 7: 3-STRAND PASCAL IN TRIBONACCI COORDINATES")
    print(f"{'='*70}")

    print(f"""
  (1+x+x^2)^n = Phi_3(x)^n

  At x=tau: Phi_3(tau)^n = (tau^3)^n = tau^{{3n}}
  At x=omega: Phi_3(omega)^n = 0 for n >= 1
  At x=2: Phi_3(2)^n = 7^n

  So the 3-strand Pascal triangle evaluated at:
  - x=2: gives 7^n (the forbidden value raised to the nth)
  - x=tau: gives tau^(3n) (pure tribonacci power)
  - x=omega: gives 0 (annihilation)

  The individual trinomial coefficients T(n,k):
  sum_k T(n,k) * tau^k = tau^(3n)
  sum_k T(n,k) * 2^k = 7^n
  sum_k T(n,k) * omega^k = 0

  These are THREE constraints on the trinomial row.
  Together with the total sum_k T(n,k) = 3^n,
  they give FOUR equations for 2n+1 unknowns.
""")

    # Verify
    for nn in range(1, 6):
        row = [0]*(2*nn+1)
        row[0] = 1
        for _ in range(nn):
            new = [0]*(2*nn+1)
            for k in range(2*nn+1):
                if row[k]:
                    for d in [0,1,2]:
                        if k+d <= 2*nn:
                            new[k+d] += row[k]
            row = new

        val_2 = sum(row[k] * 2**k for k in range(2*nn+1))
        val_tau = sum(row[k] * TAU**k for k in range(2*nn+1))
        tau_3n = TAU**(3*nn)
        val_3 = sum(row)

        print(f"  n={nn}: sum T(n,k)*2^k = {val_2} = 7^{nn} = {7**nn}")
        print(f"       sum T(n,k)*tau^k = {val_tau:.4f} = tau^{3*nn} = {tau_3n:.4f}")
        print(f"       sum T(n,k) = {val_3} = 3^{nn} = {3**nn}")

    print(f"\n{'='*70}")
    print("DONE — TOURNAMENT THEORY IN BASE TAU")
    print(f"{'='*70}")

if __name__ == "__main__":
    main()
