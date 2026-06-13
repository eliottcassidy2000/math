#!/usr/bin/env python3
"""
morgan_voyce_connection.py -- The Morgan-Voyce polynomial identity

DISCOVERY: The Interval Q-polynomial is a MORGAN-VOYCE POLYNOMIAL.

The Morgan-Voyce polynomial b(m, x) = sum_{j=0}^m C(m+j, 2j) * x^j
is a well-known polynomial family (OEIS A085478) with:
  - Row sums = Fibonacci(2m+1): b(m, 1) = F_{2m+1}
  - GF: sum_m b(m,x) * z^m = (1-z) / ((1-z)^2 - xz)
  - Riordan array: (1/(1-x), x/(1-x)^2)
  - Connection to Chebyshev: b(m, x) related to U_m via substitution

Our characteristic polynomial:
  P_m(t) = sum_{j=0}^m (-1)^j C(m+j, 2j) t^{m-j} = t^m * b(m, -1/t)

So the REVERSED polynomial with alternating signs is b(m, -t).

This means:
  P_m(t) has roots Q_k that satisfy:
  prod(t - Q_k) = t^m * b(m, -1/t)

CONSEQUENCES:
1. b(m, 1) = F_{2m+1} => sum_j e_j = F_{2m+1} (with signs: P_m(-1) = (-1)^m * b(m,1) = (-1)^m * F_{2m+1})
2. b(m, -1) = ? (relates to P_m(1))
3. The Chebyshev connection gives us the discriminant formula

Author: kind-pasteur-2026-03-12-S59c
"""

import cmath
import math
from itertools import combinations


def fibonacci(n):
    """Compute Fibonacci number F_n."""
    if n <= 0:
        return 0
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def morgan_voyce_b(m, x):
    """Compute Morgan-Voyce polynomial b(m, x) = sum C(m+j,2j) x^j."""
    return sum(math.comb(m + j, 2 * j) * x**j for j in range(m + 1))


def main():
    print("=" * 70)
    print("MORGAN-VOYCE POLYNOMIAL AND THE INTERVAL Q-POLYNOMIAL")
    print("=" * 70)

    # Part 1: Verify Fibonacci connection
    print("\n--- PART 1: FIBONACCI ROW SUMS ---")
    print("b(m, 1) = sum C(m+j, 2j) = F_{2m+1}")

    for m in range(1, 16):
        b_val = morgan_voyce_b(m, 1)
        fib_val = fibonacci(2*m + 1)
        print(f"  m={m:>2}: b(m,1) = {b_val:>10}, F_{2*m+1} = {fib_val:>10}, match: {b_val == fib_val}")

    # Part 2: Other special values
    print("\n--- PART 2: SPECIAL VALUES ---")
    print("b(m, x) at various x values:")

    for m in range(1, 10):
        vals = {}
        for x in [-1, 0, 1, 2, 4]:
            vals[x] = morgan_voyce_b(m, x)
        print(f"  m={m}: b(m,-1)={vals[-1]:>8}, b(m,0)={vals[0]:>3}, b(m,1)={vals[1]:>8}, "
              f"b(m,2)={vals[2]:>10}, b(m,4)={vals[4]:>12}")

    # Part 3: Connection to our polynomial P_m(t)
    print("\n--- PART 3: P_m(t) = t^m * b(m, -1/t) ---")
    print("P_m(t) = sum_{j=0}^m (-1)^j C(m+j,2j) t^{m-j}")
    print("       = t^m * sum_{j=0}^m C(m+j,2j) * (-1/t)^j")
    print("       = t^m * b(m, -1/t)")

    for m in [3, 5, 8]:
        p = 2*m + 1
        omega = cmath.exp(2j * cmath.pi / p)
        S = list(range(1, m + 1))

        Q_vals = []
        for k in range(1, m + 1):
            val = sum(omega ** (k * s) for s in S)
            Q_vals.append(abs(val)**2)

        # Verify P_m(Q_k) = 0 using t^m * b(m, -1/t)
        for k in range(m):
            qk = Q_vals[k]
            p_val = qk**m * morgan_voyce_b(m, -1/qk)
            print(f"  m={m}, k={k+1}: Q_k^m * b(m,-1/Q_k) = {p_val:.6e}")

    # Part 4: The discriminant via Morgan-Voyce
    print("\n--- PART 4: DISCRIMINANT FROM MORGAN-VOYCE ---")
    print("The Morgan-Voyce polynomial b(m,x) is related to Chebyshev U via:")
    print("  b(m, x) = U_m(1 + x/2) / U_m(1)")
    print("Wait, that's not right. Let me check the actual relation.")
    print("")
    print("Actually: b(m, x) is related to Chebyshev U via the substitution")
    print("  x = -4*sin^2(theta/2) = 2*(cos(theta) - 1)")
    print("Then b(m, 2*(cos(theta)-1)) = sin((m+1)*theta) / sin(theta) = U_m(cos(theta))")
    print("Let's verify...")

    for m in [3, 5, 8]:
        print(f"\n  m={m}:")
        for theta in [0.3, 0.7, 1.1, 1.5, 2.0]:
            x = 2 * (math.cos(theta) - 1)  # = -4*sin^2(theta/2)
            bval = morgan_voyce_b(m, x)
            uval = math.sin((m+1)*theta) / math.sin(theta)
            print(f"    theta={theta:.1f}: b(m, 2(cos-1)) = {bval:.6f}, U_m(cos) = {uval:.6f}, "
                  f"ratio = {bval/uval:.6f}" if abs(uval) > 1e-10 else f"    theta={theta:.1f}: degenerate")

    # Part 5: The ACTUAL Chebyshev relation
    print("\n--- PART 5: EXACT CHEBYSHEV RELATION ---")
    print("Testing: b(m, x-2) = U_m(x/2) where U_m is Chebyshev of 2nd kind")

    for m in [3, 5, 8]:
        print(f"\n  m={m}:")
        for x in [0.5, 1.0, 1.5, 2.5, 3.0]:
            bval = morgan_voyce_b(m, x - 2)
            # U_m(cos(theta)) = sin((m+1)*theta)/sin(theta), so U_m(x/2):
            # let cos(theta) = x/2
            if abs(x/2) <= 1:
                theta = math.acos(x/2)
                uval = math.sin((m+1)*theta) / math.sin(theta)
            else:
                # Use the polynomial form
                # U_0 = 1, U_1 = 2x, U_n = 2x*U_{n-1} - U_{n-2}
                u_prev2 = 1.0  # U_0
                u_prev1 = x    # U_1(x/2) * ... hmm
                # Actually: U_m(c) where c = x/2
                c = x/2
                u0, u1 = 1.0, 2*c
                for _ in range(m - 1):
                    u0, u1 = u1, 2*c*u1 - u0
                uval = u1

            print(f"    x={x:.1f}: b(m,x-2) = {bval:.6f}, U_m(x/2) = {uval:.6f}, "
                  f"match: {abs(bval - uval) < 0.01}")

    # Part 6: Actually, the Morgan-Voyce relation is simpler
    print("\n--- PART 6: SIMPLER CHEBYSHEV RELATION ---")
    print("From the OEIS: Morgan-Voyce b(n,x) relates to Chebyshev via:")
    print("  b(n, x^2 - 2) = ... or")
    print("  b(n, t) = (1/sqrt(t+4)) * ((sqrt(t+4) + sqrt(t))/2)^{2n+1} + symm")
    print("")
    print("Actually, let me just use the GF directly.")
    print("GF: sum_m b(m,x) z^m = (1-z)/((1-z)^2 - xz)")
    print("Setting z = 1/(t+1) might give the Q-polynomial generating function.")
    print("")
    print("Key insight: Our polynomial P_m(t) = t^m * b(m, -1/t)")
    print("The ROOTS of P_m are Q_k = sin^2(pi*r(k)/p) / sin^2(pi*k/p)")
    print("These are given by Q_k = -1/(-1/Q_k) where b(m, -1/Q_k) = 0")
    print("So the ROOTS of b(m, -s) = 0 are s_k = 1/Q_k")
    print("And sum(s_k) = e_{m-1}/e_m = (p-2)/1 = p-2")
    print("This checks out: sum 1/Q_k = p-2 = 2m-1")

    # Part 7: Row sums give Fibonacci
    print("\n--- PART 7: FIBONACCI IMPLICATIONS ---")
    print("Since e_j = C(m+j, 2j) and sum_j C(m+j, 2j) = F_{2m+1}:")
    print("  sum e_j = F_{2m+1}")
    print("This means: P_m(-1) = sum (-1)^j e_j * (-1)^{m-j} = (-1)^m * sum e_j = (-1)^m * F_{2m+1}")
    print("I.e., prod(1 + Q_k) = F_{2m+1} (up to sign)")

    for p in [5, 7, 11, 13, 17, 19, 23]:
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)
        S = list(range(1, m + 1))

        Q_vals = []
        for k in range(1, m + 1):
            val = sum(omega ** (k * s) for s in S)
            Q_vals.append(abs(val)**2)

        prod_1pQ = math.prod(1 + q for q in Q_vals)
        fib = fibonacci(2*m + 1)

        print(f"  p={p}, m={m}: prod(1+Q_k) = {prod_1pQ:.4f}, F_{2*m+1} = {fib}, "
              f"match: {abs(prod_1pQ - fib) < 0.01}")

    # Part 8: b(m, -4) and its meaning
    print("\n--- PART 8: b(m, -4) = prod(Q_k - 4) * (-1)^m ---")
    print("P_m(4) = 4^m * b(m, -1/4)")
    print("This equals prod(4 - Q_k)")

    for p in [5, 7, 11, 13, 17, 19]:
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)
        S = list(range(1, m + 1))

        Q_vals = []
        for k in range(1, m + 1):
            val = sum(omega ** (k * s) for s in S)
            Q_vals.append(abs(val)**2)

        prod_4mQ = math.prod(4 - q for q in Q_vals)
        p_at_4 = sum((-1)**j * math.comb(m+j, 2*j) * 4**(m-j) for j in range(m+1))

        # Also: at Paley, Q_k = (p+1)/4, so prod(4 - (p+1)/4) = (3p-1)^m/4^m?
        print(f"  p={p}, m={m}: prod(4-Q_k) = {prod_4mQ:.4f}, P_m(4) = {p_at_4:.4f}, "
              f"b(m,-1/4) = {morgan_voyce_b(m, -0.25):.6f}")

    # Part 9: The big picture — what does Fibonacci mean here?
    print("\n--- PART 9: THE BIG PICTURE ---")
    print("""
  THEOREM (INTERVAL Q-POLYNOMIAL = MORGAN-VOYCE):
  For the Interval circulant tournament T_p with S = {1,...,m}, m=(p-1)/2:

  1. The Fourier magnitudes Q_k = |S_hat(k)|^2 have ESFs:
     e_j(Q_1,...,Q_m) = C(m+j, 2j)  (Morgan-Voyce coefficient)

  2. The characteristic polynomial of Q_1,...,Q_m is:
     P_m(t) = t^m * b(m, -1/t)
     where b(m, x) = sum C(m+j,2j) x^j is the Morgan-Voyce polynomial.

  3. COROLLARIES:
     a) prod Q_k = 1             (constant term of P_m)
     b) sum Q_k = m(m+1)/2       (Parseval)
     c) sum 1/Q_k = p - 2        (sum of reciprocals)
     d) prod(1 + Q_k) = F_{2m+1} (Fibonacci!)
     e) disc = p^{m-1}           (discriminant)

  4. CHEBYSHEV PROOF:
     Q_k = (1 - T_m(cos(2*pi*k/p))) / (1 - cos(2*pi*k/p))
     where T_m is the Chebyshev polynomial of the first kind.
     The cos(2*pi*k/p) are roots of the minimal polynomial Psi_p,
     and the resultant Res(Psi_p(x), (1-x)*t - (1-T_m(x))) = P_m(t).

  5. WHY FIBONACCI?
     prod(1 + Q_k) = prod(1 + sin^2(pi*r(k)/p)/sin^2(pi*k/p))
                   = prod (sin^2(pi*k/p) + sin^2(pi*r(k)/p)) / sin^2(pi*k/p)
     This equals F_{2m+1} by the Morgan-Voyce identity b(m,1) = F_{2m+1}.
     The Fibonacci number counts NONDECREASING DYCK PATHS of semilength m+1.
""")

    # Part 10: Table of all identities
    print("--- COMPLETE TABLE ---")
    print(f"{'p':>4} {'m':>3} {'sum Q':>8} {'prod Q':>7} {'sum 1/Q':>8} {'prod(1+Q)':>10} {'F(2m+1)':>10} {'disc':>14}")
    print("-" * 70)

    for p in [5, 7, 11, 13, 17, 19, 23, 29, 31]:
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)
        S = list(range(1, m + 1))

        Q_vals = []
        for k in range(1, m + 1):
            val = sum(omega ** (k * s) for s in S)
            Q_vals.append(abs(val)**2)

        sum_Q = sum(Q_vals)
        prod_Q = math.prod(Q_vals)
        sum_inv_Q = sum(1/q for q in Q_vals)
        prod_1pQ = math.prod(1 + q for q in Q_vals)
        fib = fibonacci(2*m + 1)

        disc = 1.0
        for i in range(m):
            for j in range(i + 1, m):
                disc *= (Q_vals[i] - Q_vals[j])**2

        print(f"{p:>4} {m:>3} {sum_Q:>8.1f} {prod_Q:>7.1f} {sum_inv_Q:>8.1f} {prod_1pQ:>10.1f} {fib:>10} {disc:>14.0f}")


if __name__ == '__main__':
    main()
