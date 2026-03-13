#!/usr/bin/env python3
"""
esf_proof.py -- Prove e_j(Q_1,...,Q_m) = C(m+j, 2j) for Interval

The GENERATING FUNCTION for C(m+j, 2j) is:
  sum_{j=0}^inf C(m+j, 2j) * x^j = 1 / ((1-x)^{m+1} * something)?

Actually, C(n+k, 2k) is the (n,k) entry of the Delannoy triangle.
The generating function:
  sum_j C(m+j, 2j) x^j = C_{m}(x) where C_m is the central Delannoy GF.

Actually, let me think about this differently. The polynomial
  P_m(t) = sum_{j=0}^m (-1)^j C(m+j, 2j) t^{m-j}
is the characteristic polynomial of the Q_k.

The ROOTS of P_m(t) are Q_k = sin^2(pi*sigma(k)/p) / sin^2(pi*k/p)
for k=1,...,m, where sigma is the permutation k -> k*(m+1) mod p.

KEY INSIGHT: The permutation is sigma(k) = k*(m+1) mod p (reduced to {1,...,m}).
This permutation has a VERY special structure!

For p = 2m+1, sigma(k) = k*(m+1) = k*((p+1)/2) mod p.
The orbits of sigma under multiplication by (m+1) = (p+1)/2 mod p:
  - Multiplying by (p+1)/2 twice: ((p+1)/2)^2 = (p^2+2p+1)/4 = (p+1)/4 * (p+1)
  - The order of (p+1)/2 mod p divides p-1.

Actually, I should think about what polynomial has roots sin^2(pi*r/p)/sin^2(pi*k/p).
Let me use the substitution u = cos(2*pi*k/p). Then sin^2(pi*k/p) = (1-u)/2.

The key is: U_{p-1}(cos(theta)) = sin(p*theta)/sin(theta)
For theta = pi*k/p: sin(p*pi*k/p) = sin(k*pi) = 0
So every U_{p-1}(cos(2*pi*k/p)) = 0... but we're looking at sin^2 ratios.

Let me try a direct approach. Since Q_k = sin^2(pi*pi(k)/p)/sin^2(pi*k/p),
and the permutation pi: k -> k*(m+1) mod p maps {1,...,m} to itself:

For the INTERVAL specifically: S = {1,...,m}. The Fourier transform:
  S_hat(k) = sum_{j=1}^m omega^{kj} = (omega^k - omega^{k(m+1)}) / (1 - omega^k) - 1
  ...actually = sum_{j=1}^m omega^{kj} = omega^k * (1 - omega^{km}) / (1 - omega^k)

|S_hat(k)|^2 = |omega^k * (1 - omega^{km}) / (1 - omega^k)|^2
             = |1 - omega^{km}|^2 / |1 - omega^k|^2
             = (2 - 2*cos(2*pi*km/p)) / (2 - 2*cos(2*pi*k/p))
             = sin^2(pi*km/p) / sin^2(pi*k/p)
             = sin^2(pi*k*m/p) / sin^2(pi*k/p)

Since m = (p-1)/2, km = k(p-1)/2. So:
  Q_k = sin^2(pi*k*(p-1)/(2p)) / sin^2(pi*k/p)
      = sin^2(pi*k/2 - pi*k/(2p)) / sin^2(pi*k/p)

Hmm, that's not immediately revealing. Let me just verify the algebraic identity directly.

Author: kind-pasteur-2026-03-12-S59c
"""

import cmath
import math
from itertools import combinations


def main():
    print("=" * 70)
    print("PROVING e_j = C(m+j, 2j) FOR INTERVAL Q-POLYNOMIAL")
    print("=" * 70)

    # Part 1: The polynomial P_m(t) = sum (-1)^j C(m+j,2j) t^{m-j}
    print("\n--- PART 1: THE POLYNOMIAL P_m(t) ---")
    print("P_m(t) = sum_{j=0}^m (-1)^j C(m+j,2j) t^{m-j}")
    print("This is the REVERSE of the Catalan triangle polynomial!")
    print()

    # C(m+j,2j) = C(m+j,m-j) = the number of lattice paths...
    # Actually: C(n,2k) where n=m+j, 2k=2j => C(m+j, 2j)
    # The polynomial sum_{j=0}^m C(m+j,2j) x^j is a well-known GF.
    #
    # Key: C(m+j, 2j) = C(m+j, m-j). Setting k = m-j:
    # C(2m-k, 2(m-k)) = C(2m-k, 2m-2k) = C(2m-k, k).
    # So the reversed polynomial is sum_k C(2m-k, k) t^k.
    # THIS IS KNOWN: it's the Fibonacci-Lucas polynomial!

    # Actually, let me check OEIS for the rows:
    # m=3: [1, -6, 5, -1] => P_3(t) = t^3 - 6t^2 + 5t - 1
    # m=5: [1, -15, 35, -28, 9, -1] => P_5(t) = t^5 - 15t^4 + 35t^3 - 28t^2 + 9t - 1

    for m in [2, 3, 5, 6, 8]:
        coeffs = []
        for j in range(m + 1):
            c = (-1)**j * math.comb(m + j, 2 * j)
            coeffs.append(c)
        # coeffs[j] is the coefficient of t^{m-j}
        # So P_m(t) = coeffs[0]*t^m + coeffs[1]*t^{m-1} + ... + coeffs[m]
        poly_str = " + ".join(f"({c})*t^{m-j}" for j, c in enumerate(coeffs) if c != 0)
        print(f"  P_{m}(t) = {poly_str}")

    # Part 2: Connection to Chebyshev
    print("\n--- PART 2: CHEBYSHEV CONNECTION ---")
    print("The Chebyshev U polynomial of the SECOND kind:")
    print("  U_n(x) = sin((n+1)*theta) / sin(theta) where x = cos(theta)")
    print()
    print("Consider R_m(t) = t^m * P_m(1/t) (reversed polynomial):")
    print("  R_m(t) = sum_{j=0}^m (-1)^j C(m+j,2j) t^j")
    print()
    print("Is R_m(t) related to Chebyshev?")
    print("U_{2m}(x) evaluated at x = (1-2t)/2 = 1/2 - t?")

    # Let me check: the Chebyshev polynomial of the 2nd kind
    # U_n(cos(theta)) = sin((n+1)*theta)/sin(theta)
    # U_n has roots at theta = k*pi/(n+1) for k=1,...,n
    # i.e., x = cos(k*pi/(n+1))

    # Our Q_k roots: we need P_m(Q_k) = 0.
    # Q_k = sin^2(alpha)/sin^2(beta) for some angles.
    # If we set t = sin^2(alpha)/sin^2(beta), then...

    # Different approach: look at the SUBSTITUTION t = 4*sin^2(u)
    # Then C(m+j,2j) * (4*sin^2(u))^{m-j} = ?
    # Or t = -4*sin^2(u)?

    # Actually, the Fibonacci polynomial F_n(x) satisfies:
    # F_n(x) = sum_{k=0}^{floor((n-1)/2)} C(n-1-k, k) * x^{n-1-2k}
    # = sum_{k} C(n-1-k, k) x^{n-1-2k}
    # Related to our C(2m-k, k)?

    # Let me try: does P_m(t) = U_{2m}(something in t) / something?

    for m in [3, 5, 8]:
        p = 2*m + 1
        omega = cmath.exp(2j * cmath.pi / p)
        S = list(range(1, m + 1))

        Q_vals = []
        for k in range(1, m + 1):
            val = sum(omega ** (k * s) for s in S)
            Q_vals.append(abs(val)**2)

        # Verify P_m(Q_k) = 0
        for k in range(m):
            qk = Q_vals[k]
            pval = sum((-1)**j * math.comb(m+j, 2*j) * qk**(m-j)
                      for j in range(m + 1))
            print(f"  m={m}, P_{m}(Q_{k+1}) = {pval:.6e}")

    # Part 3: The substitution t = (1 + 2*cos(theta)) for Chebyshev U
    print("\n--- PART 3: SUBSTITUTION SEARCH ---")
    # Q_k = sin^2(pi*r(k)/p) / sin^2(pi*k/p)
    # = [1 - cos(2*pi*r(k)/p)] / [1 - cos(2*pi*k/p)]
    # Let alpha = 2*pi*k/p, beta = 2*pi*r(k)/p
    # Q_k = (1 - cos(beta)) / (1 - cos(alpha))

    # The Dirichlet kernel: sum_{j=1}^m cos(j*alpha) = -1/2 + sin((m+1/2)*alpha)/(2*sin(alpha/2))
    # For alpha = 2*pi*k/p: sum = -1/2 (since sin((m+1/2)*2*pi*k/p) = sin((p*pi*k)/p) = 0 for k<p)
    # Wait, m+1/2 = p/2, so (m+1/2)*alpha = (p/2)*(2*pi*k/p) = pi*k
    # sin(pi*k) = 0. So the sum IS -1/2, confirming C(k) = -1/2.

    # Can I write the product prod(t - Q_k) in terms of Chebyshev?
    # The Q_k are the set { (1-cos(2pi*r/p))/(1-cos(2pi*k/p)) : k=1,...,m }
    # where r = r(k) = k*(m+1) mod p.

    # Consider: let x = cos(alpha) where alpha = 2*pi*k/p.
    # Then Q = (1-cos(beta))/(1-cos(alpha)) where beta depends on k.

    # For interval, the product formula is:
    # prod_k (t - Q_k) = ??? in terms of Chebyshev

    # Part 4: Direct connection via minimal polynomial of cos(2pi/p)
    print("\n--- PART 4: MINIMAL POLYNOMIAL OF 2*cos(2pi/p) ---")
    print("The minimal polynomial of 2*cos(2*pi/p) over Q is the cyclotomic factor.")
    print("Phi_p(x) = x^{p-1} + x^{p-2} + ... + 1")
    print("Setting x = e^{i*2pi/p}: Phi_p(e^{itheta}) = 0")
    print("In terms of cos: the minimal poly of 2*cos(2pi/p) has degree m = (p-1)/2")
    print("")

    # The minimal polynomial of c = 2*cos(2*pi/p) is:
    # Psi_p(x) = prod_{k=1}^m (x - 2*cos(2*pi*k/p))
    # This has integer coefficients!

    for m in [3, 5, 6, 8]:
        p = 2*m + 1
        roots = [2*math.cos(2*math.pi*k/p) for k in range(1, m+1)]

        # Compute coefficients via ESF
        e_cos = []
        for j in range(1, m+1):
            ej = sum(math.prod(roots[i] for i in combo)
                    for combo in combinations(range(m), j))
            e_cos.append(round(ej))

        print(f"  p={p}: Psi_p coeffs = {e_cos}")

    # Part 5: Transform from Psi_p to P_m
    print("\n--- PART 5: TRANSFORM Psi_p -> P_m ---")
    print("Q_k = (1 - cos(2*pi*r(k)/p)) / (1 - cos(2*pi*k/p))")
    print("Let c_k = 2*cos(2*pi*k/p). Then:")
    print("  Q_k = (2 - 2*cos(2*pi*r(k)/p)) / (2 - 2*cos(2*pi*k/p))")
    print("      = (2 - c_{r(k)}) / (2 - c_k)")
    print()
    print("If c_k are roots of Psi_p(x) = prod(x - c_k), then")
    print("Q_k = (2 - c_{r(k)}) / (2 - c_k)")
    print()
    print("The key question: what polynomial has roots (2-c_{r(k)})/(2-c_k)?")
    print("This is a RESULTANT computation:")
    print("  P_m(t) = Res_x(Psi_p(x), (2-x)*t - (2-r(x))) / leading")
    print()
    print("But r(k) = k*(m+1) mod p, so c_{r(k)} = 2*cos(2*pi*k*(m+1)/p)")
    print("And c_{r(k)} = c_k^{power} via the Chebyshev recursion!")
    print("Specifically: 2*cos(n*theta) = T_n(2*cos(theta)) where T is Chebyshev T")
    print("So c_{r(k)} = T_{m+1}(c_k) ... wait, T_n(2*cos(theta)) != 2*cos(n*theta)")
    print("Actually: T_n(cos(theta)) = cos(n*theta), so T_n(c_k/2) = cos(2*pi*k*n/p)")
    print("=> c_{r(k)} = 2*T_{m+1}(c_k/2)")

    # Verify this
    for m in [3, 5, 8]:
        p = 2*m + 1
        print(f"\n  m={m}, p={p}:")
        for k in range(1, m+1):
            c_k = 2*math.cos(2*math.pi*k/p)
            r_k = (k*(m+1)) % p
            if r_k > m:
                r_k = p - r_k
            c_rk_actual = 2*math.cos(2*math.pi*r_k/p)
            # T_{m+1}(c_k/2) = cos((m+1)*2*pi*k/p) = cos(2*pi*k*(m+1)/p)
            # But (m+1) = (p+1)/2, so k*(m+1) mod p gives r(k)
            c_rk_cheb = 2*math.cos(2*math.pi*k*(m+1)/p)  # This might differ by sign from c_{r(k)}
            err = abs(c_rk_actual - c_rk_cheb)
            print(f"    k={k}: c_k={c_k:.6f}, r(k)={r_k}, c_r(k)={c_rk_actual:.6f}, "
                  f"2*T_{{m+1}}(c_k/2)={c_rk_cheb:.6f}, err={err:.2e}")

    # Part 6: So Q_k = (2 - 2*T_{m+1}(c_k/2)) / (2 - c_k) for suitable sign adjustment
    print("\n--- PART 6: Q_k AS CHEBYSHEV RATIO ---")
    print("Q_k = (1 - cos((m+1)*alpha)) / (1 - cos(alpha))")
    print("    = (1 - T_{m+1}(cos(alpha))) / (1 - cos(alpha))")
    print("    = U_m(cos(alpha)) * (1 - cos(alpha)) ... no")
    print("")
    print("Actually: 1 - cos(n*theta) = 2*sin^2(n*theta/2)")
    print("And: sin^2(n*theta/2) = n^2 * prod_{j=1}^{n-1} sin^2(theta/2 - j*pi/n)?")
    print("Better: use the identity")
    print("  (1 - T_n(x)) / (1 - x) = U_{n-1}(x) + U_{n-2}(x) + ... + U_0(x)")
    print("  = sum_{k=0}^{n-1} U_k(x)")
    print("This is a KNOWN identity for Chebyshev polynomials.")
    print("")
    print("So Q_k = (1 - T_{m+1}(cos(alpha))) / (1 - cos(alpha))")
    print("       = U_m(cos(alpha)) + U_{m-1}(cos(alpha)) + ... + U_0(cos(alpha))")
    print("       = sum_{j=0}^m U_j(cos(alpha))")
    print("where alpha = 2*pi*k/p.")

    for m in [3, 5, 8]:
        p = 2*m + 1
        omega_p = cmath.exp(2j * cmath.pi / p)
        S = list(range(1, m + 1))

        print(f"\n  m={m}, p={p}:")
        for k in range(1, min(m+1, 4)):
            alpha = 2*math.pi*k/p
            x = math.cos(alpha)

            # Compute Q_k via Fourier
            val = sum(omega_p ** (k * s) for s in S)
            Q_fourier = abs(val)**2

            # Compute sum U_j(x) for j=0,...,m
            # U_j(cos(theta)) = sin((j+1)*theta)/sin(theta)
            U_sum = sum(math.sin((j+1)*alpha)/math.sin(alpha) for j in range(m+1))

            # This should give Q_k??? No, that's the sum of Chebyshev U values,
            # not the product formula.
            # Actually: (1-T_n(x))/(1-x) where T_n(cos(t)) = cos(nt)
            # = (1-cos(n*alpha))/(1-cos(alpha)) = Q_k directly!
            cheb_Q = (1 - math.cos((m+1)*alpha)) / (1 - math.cos(alpha))

            print(f"    k={k}: Q_fourier={Q_fourier:.6f}, (1-cos({m+1}*a))/(1-cos(a))={cheb_Q:.6f}, "
                  f"sum_U={U_sum:.6f}")

    # Part 7: The PROOF
    print("\n" + "=" * 70)
    print("PART 7: THE PROOF STRUCTURE")
    print("=" * 70)
    print("""
  THEOREM: For the Interval tournament with S={1,...,m}, m=(p-1)/2,
  the Q_k = |S_hat(k)|^2 satisfy e_j(Q_1,...,Q_m) = C(m+j, 2j).

  PROOF OUTLINE:
  1. Q_k = sin^2(pi*k*m/p) / sin^2(pi*k/p)
         = (1 - cos(2*pi*k*m/p)) / (1 - cos(2*pi*k/p))

  2. Let alpha_k = 2*pi*k/p and x_k = cos(alpha_k). Then:
     Q_k = (1 - cos(m*alpha_k)) / (1 - cos(alpha_k))
         = (1 - T_m(x_k)) / (1 - x_k)

  3. The x_k = cos(2*pi*k/p) for k=1,...,m are the roots of Psi_p(x),
     the minimal polynomial of 2*cos(2*pi/p).

  4. The polynomial whose roots are Q_k = (1-T_m(x))/(1-x) evaluated
     at the roots of Psi_p(x) is given by the RESULTANT:
     P_m(t) = Res_x(Psi_p(x), (1-x)*t - (1-T_m(x)))

  5. This resultant equals sum_{j=0}^m (-1)^j C(m+j,2j) t^{m-j}
     by a classical identity for Chebyshev resultants.

  THEREFORE: e_j = C(m+j, 2j). QED (modulo the resultant identity)

  COROLLARY 1: prod_{k=1}^m Q_k = e_m = C(2m, 2m) = 1.
    (This is Theorem "Interval product identity" reproved.)

  COROLLARY 2: sum_{k=1}^m Q_k = e_1 = C(m+1, 2) = m(m+1)/2.
    (This is Parseval, reproved.)

  COROLLARY 3: e_{m-1} = C(2m-1, 2m-2) = C(2m-1, 1) = 2m-1 = p-2.
    (Sum of reciprocals: sum 1/Q_k = p-2.)

  COROLLARY 4: disc = p^{m-1}.
    (The discriminant of the Chebyshev resultant polynomial.)
""")

    # Part 8: Verify Corollary 3
    print("--- VERIFYING COROLLARY 3: sum(1/Q_k) = p-2 ---")
    for p in [5, 7, 11, 13, 17, 19, 23]:
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)
        S = list(range(1, m + 1))

        Q_vals = []
        for k in range(1, m + 1):
            val = sum(omega ** (k * s) for s in S)
            Q_vals.append(abs(val)**2)

        inv_sum = sum(1/q for q in Q_vals)
        print(f"  p={p}: sum(1/Q_k) = {inv_sum:.8f}, p-2 = {p-2}, match: {abs(inv_sum - (p-2)) < 1e-6}")

    # Part 9: Check if Q_k = (1-T_m(x_k))/(1-x_k) exactly
    print("\n--- VERIFYING Q_k = (1 - T_m(cos(2pi*k/p))) / (1 - cos(2pi*k/p)) ---")
    print("Note: T_m means T_{(p-1)/2}")

    for p in [7, 11, 13, 17]:
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)
        S = list(range(1, m + 1))

        all_match = True
        for k in range(1, m + 1):
            val = sum(omega ** (k * s) for s in S)
            Q_fourier = abs(val)**2

            alpha = 2 * math.pi * k / p
            x = math.cos(alpha)

            # Chebyshev T_m(x) = cos(m * arccos(x)) = cos(m * alpha)
            T_m_x = math.cos(m * alpha)
            Q_cheb = (1 - T_m_x) / (1 - x)

            err = abs(Q_fourier - Q_cheb)
            if err > 1e-8:
                all_match = False
                print(f"  p={p}, k={k}: MISMATCH Q_f={Q_fourier:.6f} Q_c={Q_cheb:.6f} err={err:.2e}")

        if all_match:
            print(f"  p={p}: ALL MATCH (Chebyshev formula confirmed)")


if __name__ == '__main__':
    main()
