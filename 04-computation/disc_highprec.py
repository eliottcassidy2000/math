#!/usr/bin/env python3
"""
disc_highprec.py -- High-precision discriminant verification

Test disc = p^{m-1} for Interval Q-polynomial using mpmath.

Also investigate the ESF pattern: e_j for Interval.
The ESFs [3,1], [6,5,1], [15,35,28,9,1], [21,70,84,45,11,1], [36,210,462,495,286,91,15,1]
look like they might be C(2m,2j)/something.

Let me check: C(2m, 2j) for m=2: C(4,2)=6, C(4,4)=1. Close but not exact.
Actually: e_1 = m(m+1)/2 = C(m+1,2).
  m=2: C(3,2)=3 YES
  m=3: C(4,2)=6 YES
  m=5: C(6,2)=15 YES
  m=6: C(7,2)=21 YES
  m=8: C(9,2)=36 YES

And e_{m-1} = 2m-1:
  m=2: e_1=3, 2*2-1=3 YES (but e_1 already, so last non-trivial)
  m=3: e_2=5, 2*3-1=5? NO, 2*3-1=5 YES
  m=5: e_4=9, 2*5-1=9 YES
  m=6: e_5=11, 2*6-1=11 YES
  m=8: e_7=15, 2*8-1=15 YES

So e_{m-1} = p - 2 = 2m - 1.

Actually... let me try C(2m, 2j):
  m=5: C(10,2)=45, C(10,4)=210, C(10,6)=252, C(10,8)=45, C(10,10)=1
  e_j = [15, 35, 28, 9, 1]
  Nope.

What about the pattern? Let me tabulate more carefully.

Author: kind-pasteur-2026-03-12-S59c
"""

import cmath
import math
from itertools import combinations
from decimal import Decimal, getcontext

getcontext().prec = 50


def main():
    print("=" * 70)
    print("HIGH-PRECISION DISCRIMINANT AND ESF ANALYSIS")
    print("=" * 70)

    # Part 1: ESF pattern for Interval
    print("\n--- PART 1: ESF PATTERN ---")

    for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]:
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)
        S = list(range(1, m + 1))

        Q_vals = []
        for k in range(1, m + 1):
            val = sum(omega ** (k * s) for s in S)
            Q_vals.append(abs(val)**2)

        e_vals = []
        for j in range(1, m + 1):
            ej = sum(math.prod(Q_vals[i] for i in combo)
                    for combo in combinations(range(m), j))
            e_vals.append(round(ej))

        print(f"\n  p={p}, m={m}: e_j = {e_vals}")
        print(f"    e_1 = {e_vals[0]} = C({m+1},2) = {math.comb(m+1,2)} {'YES' if e_vals[0]==math.comb(m+1,2) else 'NO'}")
        print(f"    e_m = {e_vals[-1]} (should be 1)")
        if m >= 2:
            print(f"    e_{{m-1}} = {e_vals[-2]} = 2m-1 = {2*m-1} {'YES' if e_vals[-2]==2*m-1 else 'NO'}")

    # Part 2: OEIS lookup of ESF sequences
    print("\n--- PART 2: KNOWN SEQUENCES? ---")
    print("The e_j for different p:")
    print("  p=5:  [3, 1]")
    print("  p=7:  [6, 5, 1]")
    print("  p=11: [15, 35, 28, 9, 1]")
    print("  p=13: [21, 70, 84, 45, 11, 1]")
    print("  p=17: [36, 210, 462, 495, 286, 91, 15, 1]")
    print("  p=19: [45, 330, 924, 1287, 1001, 455, 120, 17, 1]")

    # These look like C(2m, 2j) / p ??
    print("\n  Checking if e_j = C(2m, 2j+2) / p ... (trying various formulas)")

    for p in [5, 7, 11, 13, 17, 19, 23, 29]:
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)
        S = list(range(1, m + 1))

        Q_vals = []
        for k in range(1, m + 1):
            val = sum(omega ** (k * s) for s in S)
            Q_vals.append(abs(val)**2)

        e_vals = []
        for j in range(1, m + 1):
            ej = sum(math.prod(Q_vals[i] for i in combo)
                    for combo in combinations(range(m), j))
            e_vals.append(round(ej))

        print(f"\n  p={p}, m={m}:")
        for j in range(m):
            # Try C(p, 2*(j+1)) / p
            binom_2j2 = math.comb(p, 2*(j+1)) // p if p > 0 else 0
            # Try C(p-1, 2*(j+1)) / something
            # Try C(p-1, 2*j+1)
            c1 = math.comb(p-1, 2*j+1) if 2*j+1 <= p-1 else 0
            c2 = math.comb(p-1, 2*(j+1)) if 2*(j+1) <= p-1 else 0
            # Try C(2m, 2j+2)/(2j+2)
            c3 = math.comb(2*m, 2*j+2) // (2*j+2) if 2*j+2 > 0 else 0

            print(f"    e_{j+1} = {e_vals[j]}, C(p,{2*(j+1)})/p = {binom_2j2}, "
                  f"C(p-1,{2*j+1}) = {c1}, C(2m,{2*j+2})/({2*j+2}) = {c3}")

    # Part 3: The generating function approach
    print("\n--- PART 3: GENERATING FUNCTION ---")
    print("prod_{k=1}^m (1 + t*Q_k) = sum_{j=0}^m e_j * t^j")
    print("If we can identify this generating function, we identify the polynomial.")

    for p in [5, 7, 11, 13, 17, 19, 23]:
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)
        S = list(range(1, m + 1))

        Q_vals = []
        for k in range(1, m + 1):
            val = sum(omega ** (k * s) for s in S)
            Q_vals.append(abs(val)**2)

        e_vals = [1.0]
        for j in range(1, m + 1):
            ej = sum(math.prod(Q_vals[i] for i in combo)
                    for combo in combinations(range(m), j))
            e_vals.append(round(ej))

        # Check: is e_j = C(p-1, 2j) / C(j, j)?
        # Actually, let me just compare e_j with C(p-1, 2j)/j! * something

        # More productively: compare with Chebyshev
        # Q_k = sin^2(pi*sigma(k)/p) / sin^2(pi*k/p)
        # The Dirichlet kernel connection: sum_{k=1}^m Q_k = m(m+1)/2 (Parseval)

        # What is sum Q_k^n (power sums)?
        p_sums = []
        for n in range(1, min(m+1, 8)):
            ps = sum(q**n for q in Q_vals)
            p_sums.append(round(ps))

        print(f"\n  p={p}, m={m}:")
        print(f"    Power sums p_n = {p_sums}")

        # Check if p_n = m(m+1)/2 * f(n)?
        # p_1 = m(m+1)/2 always
        if len(p_sums) >= 2:
            print(f"    p_2/p_1 = {p_sums[1]/p_sums[0]:.6f}")
        if len(p_sums) >= 3:
            print(f"    p_3/p_1 = {p_sums[2]/p_sums[0]:.6f}")

    # Part 4: Verify disc = p^{m-1} at larger p with arbitrary precision
    print("\n--- PART 4: DISC VERIFICATION VIA SIN PRODUCT ---")
    print("disc = prod_{i<j} (Q_i - Q_j)^2")
    print("     = prod_{i<j} (sin^2(pi*sigma(i)/p)/sin^2(pi*i/p) - sin^2(pi*sigma(j)/p)/sin^2(pi*j/p))^2")

    # Actually, let me try a different approach. Factor out:
    # Q_i - Q_j = (sin^2(alpha_i)*sin^2(beta_j) - sin^2(alpha_j)*sin^2(beta_i)) / (sin^2(beta_i)*sin^2(beta_j))
    # where alpha_k = pi*sigma(k)/p, beta_k = pi*k/p

    # The denominator of disc is (prod sin^2(pi*k/p))^{2(m-1)}
    # We know prod_{k=1}^m sin(pi*k/p) = sqrt(p)/2^m
    # So (prod sin^2(pi*k/p))^{2(m-1)} = (p/4^m)^{2(m-1)} = p^{2(m-1)}/4^{2m(m-1)}

    # For disc = p^{m-1}, we need:
    # numerator = p^{m-1} * p^{2(m-1)} / 4^{2m(m-1)} = p^{3(m-1)} / 4^{2m(m-1)}
    # That doesn't simplify nicely. Let me just verify numerically.

    for p in [23, 29, 31]:
        m = (p - 1) // 2

        # Use higher precision via mpmath-style computation
        # Actually, let's compute disc via the resultant / Sylvester formula
        # disc(f) = (-1)^{m(m-1)/2} * (1/a_m) * Res(f, f')

        # Or just use the product of differences with extended precision
        import decimal
        decimal.getcontext().prec = 100

        # Compute Q_k with high precision using mpmath if available
        try:
            import mpmath
            mpmath.mp.dps = 50

            Q_vals = []
            for k in range(1, m + 1):
                val = mpmath.mpf(0)
                for s in range(1, m + 1):
                    val += mpmath.exp(2j * mpmath.pi * k * s / p)
                Q_vals.append(abs(val)**2)

            disc = mpmath.mpf(1)
            for i in range(m):
                for j in range(i + 1, m):
                    disc *= (Q_vals[i] - Q_vals[j])**2

            disc_int = int(mpmath.nint(disc))

            n = disc_int
            p_power = 0
            while n % p == 0:
                n //= p
                p_power += 1

            print(f"\n  p={p}, m={m}: disc = p^{p_power} * {n}")
            print(f"    m-1 = {m-1}")
            print(f"    disc = p^(m-1)? {'YES' if p_power == m-1 and n == 1 else 'NO'}")

        except ImportError:
            # Fallback to standard precision with careful computation
            omega = cmath.exp(2j * cmath.pi / p)

            Q_vals = []
            for k in range(1, m + 1):
                val = sum(cmath.exp(2j * cmath.pi * k * s / p) for s in range(1, m + 1))
                Q_vals.append(abs(val)**2)

            # Instead of computing the full discriminant (overflow-prone),
            # compute log(disc) and check if it equals (m-1)*log(p)
            log_disc = 0.0
            for i in range(m):
                for j in range(i + 1, m):
                    log_disc += 2 * math.log(abs(Q_vals[i] - Q_vals[j]))

            expected_log = (m - 1) * math.log(p)
            print(f"\n  p={p}, m={m}:")
            print(f"    log(disc) = {log_disc:.10f}")
            print(f"    (m-1)*log(p) = {expected_log:.10f}")
            print(f"    ratio = {log_disc/expected_log:.10f}")
            print(f"    disc = p^(m-1)? {'YES' if abs(log_disc - expected_log) < 0.01 else 'NO'}")


if __name__ == '__main__':
    main()
