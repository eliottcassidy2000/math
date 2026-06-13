#!/usr/bin/env python3
"""
category_recurrence_unify_89.py -- opus-2026-03-14-S89

Unifying category theory, recurrences, and the Grand Theorem.

Key threads:
1. The category of tournaments: objects=tournaments, morphisms=???
2. Recurrences: Fibonacci, Padovan, tribonacci and tournament invariants
3. The functor H: Tournament -> Z (Hamiltonian path count)
4. Natural transformations and the cone

This script explores whether H can be seen as a categorical trace,
and whether the Grand Theorem is a consequence of a deeper categorical identity.
"""

from math import factorial, comb
from fractions import Fraction
from itertools import combinations


def compute_H_dp(adj_bits, n):
    """DP for Hamiltonian path count."""
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            c = dp[mask][v]
            targets = adj_bits[v] & ~mask
            u = targets
            while u:
                bit = u & (-u)
                idx = bit.bit_length() - 1
                dp[mask | bit][idx] += c
                u &= u - 1
    return sum(dp[full])


def tournament_from_bits(bits, n):
    adj = [0] * n
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits & (1 << idx):
                adj[i] |= (1 << j)
            else:
                adj[j] |= (1 << i)
            idx += 1
    return adj


def main():
    print("="*70)
    print("CATEGORICAL RECURRENCE UNIFICATION")
    print("opus-2026-03-14-S89")
    print("="*70)

    # PART 1: H as a categorical trace
    print(f"\n{'='*70}")
    print("PART 1: H AS TRACE / PERMANENT")
    print("="*70)

    # H(T) = number of Hamiltonian paths in T.
    # If A is the adjacency matrix of T (A[i][j] = 1 if i -> j, 0 otherwise):
    # H(T) = sum over permutations sigma in S_n of prod_{k=1}^{n-1} A[sigma(k)][sigma(k+1)]
    #
    # This is NOT the permanent of A (which counts perfect matchings).
    # But it's the PERMANENT of a DIFFERENT matrix.
    # Consider the (n-1) x n matrix B where B[k][j] = A[sigma(k)][j]... no.
    #
    # Actually: H(T) = sum_{sigma in S_n} prod_{i=1}^{n-1} A_{sigma(i), sigma(i+1)}
    # This is the sum over all orderings of vertices, checking if each consecutive
    # pair has the right arc.
    #
    # Equivalently: H = tr(A^{(path)}) where A^{(path)} is the transfer matrix
    # representation on the path graph P_n.
    #
    # The transfer matrix: M = A (the adjacency matrix).
    # A path of length k from v to w: (A^k)[v][w] = number of k-walks.
    # But H counts SIMPLE paths (no repeated vertices), which is harder.
    # H(T) = sum_v (DP result for paths ending at v) = sum_v dp[full][v]
    #
    # In matrix terms: define the inclusion-exclusion permanent:
    # H = sum_S (-1)^{n-|S|} ... no, the DP IS the most direct way.

    # The KEY observation: H(T) is a degree-(n-1) polynomial in the edge variables
    # (each edge is 0 or 1). This is because each Hamiltonian path uses exactly
    # n-1 edges, and H = sum over Ham paths of product of edges.

    # This makes H a MULTILINEAR polynomial of degree n-1 in m = C(n,2) variables.
    # The Walsh expansion IS this polynomial written in the {+1,-1} basis.

    print("""
  H(T) = sum_{sigma in S_n} prod_{i=1}^{n-1} T_{sigma(i), sigma(i+1)}

  This is a multilinear polynomial of degree n-1 in m=C(n,2) variables.
  Each monomial corresponds to a set of n-1 edges forming a Hamiltonian path.

  In the Walsh basis: H = sum_S c_hat_S chi_S
  where chi_S(T) = prod_{e in S} (-1)^{T_e}

  The coefficients c_hat_S are determined by the overlap between
  the Hamiltonian-path polynomial and the Walsh characters.
""")

    # PART 2: The recurrence structure of Var/Mean^2
    print(f"\n{'='*70}")
    print("PART 2: RECURRENCE FOR a(n) = V(n)*n!/2")
    print("="*70)

    # a(n) = sum_{k=1}^K (n-2k)^k * (n-2k)!
    # For odd n: a(n) = (n-2)*(n-2)! + (n-4)^2*(n-4)! + ... + 1
    # For even n: a(n) = (n-2)*(n-2)! + (n-4)^2*(n-4)! + ... + 2^{(n-2)/2}*2!

    # Can we find a recurrence?
    # a(n) - a(n-2) = new term from one more level
    # When n increases by 2: K increases by 1, and all existing terms change.
    # a(n) = (n-2)*(n-2)! + sum_{k=2}^K (n-2k)^k*(n-2k)!
    # a(n-2) = (n-4)*(n-4)! + sum_{k=2}^{K-1} (n-2-2k)^k*(n-2-2k)!
    #         = (n-4)*(n-4)! + sum_{k=2}^{K-1} (n-2(k+1))^k*(n-2(k+1))!
    # Setting j=k+1: a(n-2) = (n-4)*(n-4)! + sum_{j=3}^K (n-2j)^{j-1}*(n-2j)!

    # So a(n) - a(n-2) = (n-2)*(n-2)! - (n-4)*(n-4)!
    #   + sum_{k=2}^K [(n-2k)^k - (n-2k)^{k-1}] * (n-2k)!
    # = (n-2)! * [(n-2) - (n-4)*(n-4)!/(n-2)!] ...
    # This is messy. Let me just compute the recurrence numerically.

    def a_n(n):
        K = (n - 1) // 2
        total = 0
        for k in range(1, K + 1):
            j = n - 2 * k
            total += j**k * factorial(j)
        return total

    vals = {n: a_n(n) for n in range(3, 25)}

    # Try: a(n) = f(n)*a(n-1) + g(n)*a(n-2)
    # For two consecutive n, solve for f, g:
    # a(n) = f*a(n-1) + g*a(n-2)
    # a(n+1) = f'*a(n) + g'*a(n-1)
    # This gives f(n) as a function of n.

    print("  Solving for f(n) assuming a(n) = f(n)*a(n-1) + g(n)*a(n-2):")
    for n in range(5, 20):
        # a(n) = f*a(n-1) + g*a(n-2)
        # a(n-1) = f'*a(n-2) + g'*a(n-3)
        # Two equations, two unknowns (f, g) for fixed n:
        # a(n) = f*a(n-1) + g*a(n-2)
        # a(n+1) = f'*a(n) + g'*a(n-1)
        # But f, g might depend on n.

        # Let's just compute f and g for each n:
        # a(n) = f(n)*a(n-1) + g(n)*a(n-2)
        # a(n-1) = f(n-1)*a(n-2) + g(n-1)*a(n-3)

        # Two equations for f(n), g(n):
        # a(n) = f*vals[n-1] + g*vals[n-2]
        # a(n+1) = ?*vals[n] + ?*vals[n-1]  -- different f, g

        # Just find f(n) given g(n) = some value.
        # Try g=0: f(n) = a(n)/a(n-1)
        r = Fraction(vals[n], vals[n-1])
        # The ratio a(n)/a(n-1):
        decimal = float(r)
        excess = decimal - (n-1)
        print(f"    n={n}: a(n)/a(n-1) = {decimal:.8f}, excess over n-1: {excess:.8f}")

    # The ratio a(n)/a(n-1) is approximately n + c/n for some c.
    # More precisely: a(n)/a(n-1) = n - 1 + correction
    # Let me compute the correction more carefully.

    print(f"\n  Precise excess: a(n)/a(n-1) - (n-1):")
    for n in range(5, 20):
        excess = Fraction(vals[n], vals[n-1]) - (n-1)
        # Is excess = c/(n-1) for some constant c?
        c_estimate = float(excess) * (n-1)
        print(f"    n={n}: excess = {float(excess):.8f}, c_est = excess*(n-1) = {c_estimate:.6f}")

    # c_estimate should converge if the formula is a(n)/a(n-1) = n-1 + c/(n-1)
    # Let's see if c_estimate -> some constant.
    # From the data: c_est grows, so the formula isn't a(n)/a(n-1) = n-1 + c/(n-1).
    # Maybe: a(n)/a(n-1) = n-1 + (n-3)^2/(n-2)! * something?

    # The dominant term is a(n) ~ (n-2)*(n-2)! and a(n-1) ~ (n-3)*(n-3)!
    # Ratio ~ (n-2)*(n-2)!/[(n-3)*(n-3)!] = (n-2)*(n-2)/(n-3) = (n-2)^2/(n-3)
    # For n=10: (8)^2/7 = 9.14, actual = 9.106. Close!
    print(f"\n  Dominant term ratio (n-2)^2/(n-3):")
    for n in range(5, 20):
        predicted = (n-2)**2 / (n-3)
        actual = float(Fraction(vals[n], vals[n-1]))
        print(f"    n={n}: (n-2)^2/(n-3)={predicted:.6f}, actual={actual:.6f}, diff={actual-predicted:.6f}")

    # Very close! The difference is small and decreasing.
    # So: a(n)/a(n-1) ~ (n-2)^2/(n-3) for large n.
    # And (n-2)^2/(n-3) = n - 1 + 1/(n-3).
    # So the correction is 1/(n-3) + higher order.

    # PART 3: The three recurrences and tournament hierarchy
    print(f"\n{'='*70}")
    print("PART 3: THREE RECURRENCE FAMILIES")
    print("="*70)

    print("""
  THREE RECURRENCE FAMILIES in tournament theory:

  1. FIBONACCI (phi: x^2 = x + 1):
     - Counts tilings of [n] by 1-tiles and 2-tiles
     - Counts Hamiltonian paths on path graphs
     - Period pi(4) = 6 governs tournament parity mod 4
     - phi = (1+sqrt(5))/2

  2. PADOVAN (p: x^3 = x + 1):
     - Counts tilings of [n] by 2-tiles and 3-tiles
     - The {2,3} = {edge, triangle} building blocks of tournaments
     - Plastic number p = 1.3247... governs the "cycle complexity growth"
     - Real root of x^3 - x - 1 = 0

  3. TRIBONACCI (tau: x^3 = x^2 + x + 1):
     - Counts tilings of [n] by 1,2,3-tiles
     - tau^3 = Phi_3(tau) connects to cyclotomic structure
     - tau = 1.8393... governs the "full complexity growth"
     - Real root of x^3 - x^2 - x - 1 = 0

  Connection: Phi_3(x) = x^2 + x + 1 evaluates to:
    Phi_3(1) = 3 (the triangle)
    Phi_3(2) = 7 (the Fano plane / first forbidden value)
    Phi_3(phi) = phi^2 + phi + 1 = phi+1+phi+1 = 2phi+2 = 1+sqrt(5)+2 = 3+sqrt(5)
    Phi_3(p) = p^2 + p + 1 (where p^3 = p+1, so p^2 = (p+1)/p = 1 + 1/p)
            = 1 + 1/p + p + 1 = 2 + p + 1/p

  The three constants phi, p, tau are the "metallic means" of tournaments:
  phi governs the parity period,
  p governs the cycle-edge decomposition,
  tau governs the full spectral structure.
""")

    # Compute Phi_3 at phi, p, tau
    from math import sqrt
    phi = (1 + sqrt(5)) / 2
    # p = real root of x^3 - x - 1 = 0
    # Use Newton's method
    x = 1.3
    for _ in range(100):
        x = x - (x**3 - x - 1)/(3*x**2 - 1)
    p = x
    # tau = real root of x^3 - x^2 - x - 1 = 0
    x = 1.8
    for _ in range(100):
        x = x - (x**3 - x**2 - x - 1)/(3*x**2 - 2*x - 1)
    tau = x

    print(f"  phi = {phi:.10f}")
    print(f"  p   = {p:.10f}")
    print(f"  tau = {tau:.10f}")
    print(f"")
    print(f"  Phi_3(phi) = {phi**2 + phi + 1:.10f} = {3 + sqrt(5):.10f} (= 3 + sqrt(5))")
    print(f"  Phi_3(p)   = {p**2 + p + 1:.10f}")
    print(f"  Phi_3(tau) = {tau**2 + tau + 1:.10f} = tau^3 = {tau**3:.10f}")

    # Phi_3(tau) = tau^3! This is the tribonacci defining equation.
    # So tau is the unique real number where Phi_3(x) = x^3.

    # PART 4: The "tournament spectrum" of recurrences
    print(f"\n{'='*70}")
    print("PART 4: RECURRENCE SPECTRUM")
    print("="*70)

    # For each constant (phi, p, tau), what is the associated recurrence
    # sequence starting from H-values?

    # Fibonacci-like: a(n) = a(n-1) + a(n-2), starting from a(0)=1, a(1)=H(T)
    # Padovan-like: a(n) = a(n-2) + a(n-3)
    # Tribonacci-like: a(n) = a(n-1) + a(n-2) + a(n-3)

    # How does H interact with these?
    # If H(T) = k for a tournament T, what is the "Fibonacci transform" of k?

    # Actually, the more interesting question:
    # Does the SEQUENCE of H-means (Mean(H) at n=1,2,3,...) satisfy a recurrence?
    # Mean(H) = n!/2^{n-1}. Ratio: Mean(n+1)/Mean(n) = (n+1)/2.
    # So Mean(n) = Mean(3) * prod_{k=4}^n k/2 = (3/2) * prod k/2
    # This is not Fibonacci but factorial growth.

    print("  Mean(H) sequence and ratios:")
    for n in range(3, 12):
        mean = Fraction(factorial(n), 2**(n-1))
        if n > 3:
            prev_mean = Fraction(factorial(n-1), 2**(n-2))
            ratio = mean / prev_mean
            print(f"    n={n}: Mean = {float(mean):.4f}, ratio = {ratio} = {float(ratio):.4f}")
        else:
            print(f"    n={n}: Mean = {float(mean):.4f}")

    # Mean(n)/Mean(n-1) = n/2 exactly. Simple!

    # PART 5: Deletion-contraction and the DC recurrence
    print(f"\n{'='*70}")
    print("PART 5: DELETION-CONTRACTION TREE")
    print("="*70)

    # H(T) satisfies deletion-contraction on edges:
    # H(T) = H(T\e) + H(T/e)  (for any edge e)
    # where T\e deletes edge e (making it a "non-edge")
    # and T/e contracts edge e (merging the two endpoints).

    # Wait: for TOURNAMENTS, every edge has an orientation. Deletion means
    # removing the arc, which turns two vertices into an "incomparable pair."
    # This takes us outside the tournament world.
    # So DC for tournaments is:
    # For arc i->j: H(T) = H(T with i->j) = ?
    # If we flip i->j to j->i: H(T') might differ.
    # H(T) - H(T') = the contribution of paths using the i->j arc.

    # Actually the standard DC for H on tournaments:
    # Fix an arc e = (u,v). Then:
    # Ham paths NOT using e: these are Ham paths that don't go u->v.
    # Ham paths USING e: these go ...->u->v->... somewhere.
    # H(T) = H_{not using e} + H_{using e}
    # H_{using e} = number of Ham paths of T/e (contracted) times some factor?
    # Not quite standard. Let me just verify for small cases.

    # For n=3 (3-cycle 0->1->2->0):
    # Arc 0->1: paths using it = {0->1->2}. Paths not = {2->0->1}.
    # H = 3 (the 3 cyclic paths). Wait, n=3 cycle has H=3:
    # Paths: 0->1->2, 1->2->0, 2->0->1. All use at least one arc.
    # Paths using arc 0->1: 0->1->2, 2->0->1 (the 0->1 appears). Count = 2.
    # Paths not using arc 0->1: 1->2->0. Count = 1.
    # Total: 2+1=3. Good.

    # For the FLIP: if we flip 0->1 to 1->0:
    # New tournament: 1->0, 0->2->0?? No: 1->0, 2->0, 1->2. Wait.
    # Original: 0->1, 1->2, 2->0.
    # Flip 0->1 to 1->0: now 1->0, 1->2, 2->0. This is transitive (1->0, 1->2, 2->0).
    # H(transitive on 3) = 1. So flipping one arc changes H from 3 to 1.
    # The "edge sensitivity" = 2.

    # For n=4,5,6: compute average edge sensitivity
    print("  Edge sensitivity |H(T) - H(T')| when flipping one arc:")
    for n in range(3, 7):
        m = n * (n - 1) // 2
        total_sensitivity = 0
        count = 0
        for bits in range(1 << m):
            adj = tournament_from_bits(bits, n)
            h = compute_H_dp(adj, n)
            for e in range(m):
                # Flip edge e
                bits_flipped = bits ^ (1 << e)
                adj_f = tournament_from_bits(bits_flipped, n)
                hf = compute_H_dp(adj_f, n)
                total_sensitivity += abs(h - hf)
                count += 1

        avg = total_sensitivity / count
        print(f"    n={n}: avg |dH/de| = {avg:.4f}")

    # PART 6: The categorical viewpoint
    print(f"\n{'='*70}")
    print("PART 6: CATEGORICAL PERSPECTIVE")
    print("="*70)

    print("""
  CATEGORY OF TOURNAMENTS (Tourn):
    Objects: tournaments T on [n] for various n
    Morphisms: tournament homomorphisms f: T -> T'
      (vertex maps preserving arc direction)

  FUNCTORS from Tourn:
    H: Tourn -> Z  (Hamiltonian path count)
      H is NOT a functor (doesn't preserve morphisms).
      But it IS an INVARIANT: H(T) = H(sigma(T)) for all sigma in S_n.

    Cone: Tourn_n -> Tourn_{n+1}  (cone construction)
      This IS a functor: it maps T to Cone(T), and morphisms f to Cone(f).
      THM-205: H . Cone = H (the cone preserves H).
      In categorical language: H is a NATURAL INVARIANT of Cone.

    Complement: Tourn -> Tourn  (reverse all arcs)
      THM-203: H . Complement = H (complement preserves H).
      Complement is an INVOLUTION (Complement^2 = Id).

  NATURAL TRANSFORMATION perspective:
    The cone functor Cone: Tourn_n -> Tourn_{n+1} satisfies:
    H_{n+1} . Cone = H_n (as invariants)
    This means H is "stable" in the sense of algebraic K-theory.
    The stabilization: colim_n H_n gives a "universal" invariant.

  MONOIDAL STRUCTURE:
    Tourn has a product: T x T' = direct product tournament.
    H is NOT multiplicative: H(T x T') != H(T) * H(T').
    But the complement involution IS monoidal: (T x T')^bar = T^bar x T'^bar.

  THE KEY QUESTION:
    Is there a "derived" category of tournaments where H becomes a trace?
    In other words: is there a chain complex C(T) such that
    chi(C(T)) = H(T) (Euler characteristic = Hamiltonian path count)?
""")

    # PART 7: Euler characteristic interpretation
    print(f"\n{'='*70}")
    print("PART 7: H AS EULER CHARACTERISTIC?")
    print("="*70)

    # H(T) is always odd (Redei's theorem applied to each starting vertex).
    # An Euler characteristic is typically: chi = sum (-1)^k beta_k.
    # Can we find beta_k such that sum (-1)^k beta_k = H(T)?

    # From GLMY path homology: beta_k(T) are the Betti numbers.
    # But chi_GLMY(T) = sum (-1)^k beta_k^GLMY is NOT H(T) in general.

    # From the ODD CYCLE COLLECTION: H(T) = sum_S 2^|S| where S ranges
    # over independent sets of odd cycles. This is 2^{beta_1} when
    # all cycles are disjoint (i.e., form a matching in the cycle graph).
    # Not quite Euler characteristic.

    # However: H(T) = I(Omega(T), 2) where I(G, x) is the independence polynomial
    # of the odd-cycle graph Omega(T), evaluated at x=2.
    # The independence polynomial is: I(G, x) = sum_k a_k x^k where a_k = number
    # of independent sets of size k.
    # H(T) = I(Omega, 2) = sum_k a_k 2^k.

    # This is NOT an Euler characteristic but a SPECIALIZATION of the
    # Tutte polynomial or independence polynomial.

    # The connection to Euler characteristic would be:
    # If we could find a CW-complex X(T) such that:
    # X has exactly a_k cells of dimension k (for each k)
    # Then chi(X) = sum (-1)^k a_k, which is I(Omega, -1).
    # But H = I(Omega, 2), not I(Omega, -1).

    # What IS I(Omega, -1)?
    # I(G, -1) = sum_k a_k (-1)^k = number of independent sets with even size
    # minus number with odd size.

    # For n=3 (cycle): Omega has 1 cycle (the 3-cycle). a_0=1, a_1=1.
    # I = 1 + 2 = 3 at x=2. I(-1) = 1 - 1 = 0. Interesting!

    # For n=3 (transitive): Omega empty. a_0=1. I = 1 at x=2. I(-1) = 1.

    # For a general graph: I(G, -1) = 0 if G has an edge (by inclusion-exclusion).
    # Wait no: I(G, -1) = sum over indep sets S of (-1)^|S|.
    # For the empty graph: I = (1+(-1))^n = 0^n = 0 for n>0, 1 for n=0. Hmm.
    # Actually independence polynomial of K_0 (no vertices) = 1.
    # For a graph with at least one vertex: if we include the empty set (size 0),
    # I(G, -1) depends on the structure.

    # Let me just compute for small cases
    for n in range(3, 7):
        m = n * (n - 1) // 2
        i_minus_1_vals = set()
        for bits in range(1 << m):
            adj = tournament_from_bits(bits, n)
            h = compute_H_dp(adj, n)
            # We'd need to compute I(Omega, -1) which requires the cycle graph.
            # Skip for now.

        # Instead, check: is H(T) mod 2 always 1?
        h_mod2 = set()
        for bits in range(1 << m):
            adj = tournament_from_bits(bits, n)
            h = compute_H_dp(adj, n)
            h_mod2.add(h % 2)
        print(f"  n={n}: H mod 2 values: {h_mod2}")

    print(f"\n  H(T) is always ODD (confirmed n=3..6). This is Redei's theorem.")
    print(f"  Write H = 2Q + 1 for some integer Q >= 0.")
    print(f"  Q(T) = (H(T)-1)/2 IS a natural candidate for an Euler-like characteristic.")
    print(f"  Q(transitive) = 0, Q(3-cycle) = 1, Q(max at n=7) = 94.")

    print(f"\n{'='*70}")
    print("DONE -- CATEGORICAL RECURRENCE UNIFICATION")
    print("="*70)


if __name__ == "__main__":
    main()
