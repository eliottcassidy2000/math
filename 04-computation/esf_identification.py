#!/usr/bin/env python3
"""
esf_identification.py -- Identify the Interval ESF sequence

For the Interval tournament T_p with connection set {1,...,m}, m=(p-1)/2,
the Q_k = |S_hat(k)|^2 have elementary symmetric functions e_j.

Known:
  e_m = 1 (product = 1)
  e_1 = m(m+1)/2
  e_{m-1} = p - 2 = 2m - 1

Full sequences:
  m=2:  [3, 1]
  m=3:  [6, 5, 1]
  m=5:  [15, 35, 28, 9, 1]
  m=6:  [21, 70, 84, 45, 11, 1]
  m=8:  [36, 210, 462, 495, 286, 91, 15, 1]
  m=9:  [45, 330, 924, 1287, 1001, 455, 120, 17, 1]
  m=11: [66, 715, 3003, 6435, 8008, 6188, 3060, 969, 190, 21, 1]
  m=14: [105, 1820, 12376, 43758, 92378, 125970, 116280, 74613, 33649, 10626, 2300, 325, 27, 1]

These are EXACTLY C(p-1, 2j) = C(2m, 2j)! Let me verify:
  m=3: C(6,2)=15, C(6,4)=15, C(6,6)=1
  Actual: [6, 5, 1] — NO, not C(2m, 2j).

Wait — power sums had ratio p_2/p_1 = (2p-3)/3 at p=5 (7/3=2.33).
Actually no. Let me be more careful about what these numbers are.

Maybe C(2m+1, 2j+1) / (2j+1) (ballot / Catalan-like)?
  m=3: C(7,3)/3 = 35/3 NO

Or C(p, 2j) / p?
  m=3: C(7,2)/7 = 3, C(7,4)/7 = 5, C(7,6)/7 = 1 => [3, 5, 1] YES!

Let me check all:
  p=5: C(5,2)/5=2, C(5,4)/5=1 => [2, 1] vs actual [3, 1] NO
  p=7: C(7,2)/7=3, C(7,4)/7=5, C(7,6)/7=1 => [3, 5, 1] vs actual [6, 5, 1] NO

Hmm. The MIDDLE entries sometimes match but not the edges.

Author: kind-pasteur-2026-03-12-S59c
"""

import cmath
import math
from itertools import combinations


def compute_esf(p):
    """Compute ESFs for Interval Q-polynomial at prime p."""
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
    return e_vals, Q_vals


def main():
    print("=" * 70)
    print("ESF IDENTIFICATION FOR INTERVAL Q-POLYNOMIAL")
    print("=" * 70)

    # Collect all ESF sequences
    all_esf = {}
    for p in [5, 7, 11, 13, 17, 19, 23, 29, 31]:
        m = (p - 1) // 2
        e_vals, Q_vals = compute_esf(p)
        all_esf[p] = e_vals
        print(f"  p={p}, m={m}: {e_vals}")

    # Part 1: Look at the reversed sequence (since e_m=1, e_{m-1}=p-2)
    print("\n--- PART 1: REVERSED ESFs ---")
    for p, e_vals in all_esf.items():
        rev = list(reversed(e_vals))
        print(f"  p={p}: reversed = {rev}")

    # Part 2: Central binomial / Catalan comparison
    print("\n--- PART 2: FORMULA SEARCH ---")
    # The key observation from earlier: at p=7, C(7,4)/7 = 5 = e_2.
    # And C(7,6)/7 = 1 = e_3. But C(7,2)/7 = 3 != 6 = e_1.
    # However, e_1 = m(m+1)/2 = 6, and C(7,2)/7 = 3.
    # So e_1 = 2 * C(p,2)/p? 2*3 = 6 YES!
    # Check: p=11, e_1=15, C(11,2)/11 = 5, 3*5=15 YES!
    # p=13, e_1=21, C(13,2)/13 = 6, 3.5*6=21? NO, 21/6 = 3.5

    # Maybe e_j = C(p, 2j) * A(j) / p for some sequence A?
    print("Trying: e_j * p = C(p, 2j) * ratio(j)")
    for p, e_vals in all_esf.items():
        m = (p - 1) // 2
        print(f"\n  p={p}:")
        for j in range(1, m + 1):
            c = math.comb(p, 2*j)
            if c > 0:
                ratio = e_vals[j-1] * p / c
                print(f"    e_{j} * p / C(p,{2*j}) = {e_vals[j-1]} * {p} / {c} = {ratio:.6f}")

    # Part 3: Try e_j = C(2m, 2j-1) / (2j-1) or similar
    print("\n--- PART 3: FORMULA CANDIDATES ---")
    for p, e_vals in all_esf.items():
        m = (p - 1) // 2
        print(f"\n  p={p}, m={m}:")

        for j in range(1, m + 1):
            # Various candidates
            candidates = {}
            candidates[f'C(2m,2j-1)/(2j-1)'] = math.comb(2*m, 2*j-1) / (2*j-1) if 2*j-1 <= 2*m else 0
            candidates[f'C(2m+1,2j)/(2j)'] = math.comb(2*m+1, 2*j) / (2*j) if 2*j <= 2*m+1 else 0
            candidates[f'C(2m+1,2j+1)/(2j+1)'] = math.comb(2*m+1, 2*j+1) / (2*j+1) if 2*j+1 <= 2*m+1 else 0
            candidates[f'C(p,2j)/p'] = math.comb(p, 2*j) / p if 2*j <= p else 0
            candidates[f'C(p-2,2j-2)'] = math.comb(p-2, 2*j-2) if 2*j-2 <= p-2 else 0

            matches = []
            for name, val in candidates.items():
                if abs(val - e_vals[j-1]) < 0.001:
                    matches.append(name)

            if matches:
                print(f"    e_{j} = {e_vals[j-1]} = {', '.join(matches)}")
            else:
                # Show all values for debugging
                vals_str = ", ".join(f"{name}={val:.1f}" for name, val in candidates.items())
                print(f"    e_{j} = {e_vals[j-1]} ~ {vals_str}")

    # Part 4: Check C(2m, 2j-1)/(2j-1) systematically
    print("\n--- PART 4: TESTING C(2m, 2j-1)/(2j-1) ---")
    for p, e_vals in all_esf.items():
        m = (p - 1) // 2
        match = True
        for j in range(1, m + 1):
            expected = math.comb(2*m, 2*j-1) // (2*j-1) if (2*j-1) <= 2*m else 0
            actual_float = math.comb(2*m, 2*j-1) / (2*j-1)
            if abs(actual_float - e_vals[j-1]) > 0.001:
                match = False
                break
        print(f"  p={p}: C(2m,2j-1)/(2j-1) = {'YES' if match else 'NO'}")

    # Part 5: Central Delannoy / Narayana
    print("\n--- PART 5: NARAYANA / MOTZKIN CHECK ---")
    # Narayana N(n,k) = C(n,k) * C(n,k-1) / n
    # These grow like our sequence...

    for p, e_vals in all_esf.items():
        m = (p - 1) // 2
        print(f"\n  p={p}, m={m}:")
        for j in range(1, m + 1):
            # Narayana N(m+1, j)
            nar = math.comb(m+1, j) * math.comb(m+1, j-1) // (m+1) if j <= m+1 else 0
            # Catalan convolution
            catconv = math.comb(2*j, j) * math.comb(2*(m-j), m-j) // ((m-j+1)*(j+1)) if j <= m else 0
            # Self-conjugate partitions related
            # Actually, let's try: number of lattice paths
            print(f"    e_{j} = {e_vals[j-1]}, Narayana({m+1},{j}) = {nar}")

    # Part 6: Power sums and their structure
    print("\n--- PART 6: POWER SUM FORMULAS ---")
    print("p_n = sum_{k=1}^m Q_k^n")
    print("p_1 = m(m+1)/2 = e_1")
    print("p_2 = e_1^2 - 2*e_2 (Newton)")
    print("")
    print("Testing: p_2 / m = (2p^2 - 6p + 1) / (3p) ?")

    for p, e_vals in all_esf.items():
        m = (p - 1) // 2
        p2 = e_vals[0]**2 - 2*e_vals[1]
        expected_p2_formula = m * (2*m + 1) * (2*m - 1) / 3  # = m(4m^2-1)/3
        print(f"  p={p}: p_2 = {p2}, m*(4m^2-1)/3 = {expected_p2_formula:.1f}, match: {abs(p2 - expected_p2_formula) < 0.01}")

    # Part 7: Q_k as sin^2 ratios — exact algebraic formula
    print("\n--- PART 7: EXACT ALGEBRAIC FORM ---")
    print("Q_k = sin^2(pi*pi_k/p) / sin^2(pi*k/p)")
    print("where pi is the permutation: pi_k = k*(m+1) mod p (reduced to {1,...,m})")

    for p in [5, 7, 11, 13, 17]:
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)

        # Find the permutation
        perm = []
        for k in range(1, m + 1):
            r = (k * (m + 1)) % p
            if r > m:
                r = p - r
            perm.append(r)

        print(f"\n  p={p}, m={m}: permutation = {perm}")

        # The Q_k via sin^2 ratio
        Q_sin = []
        for k in range(1, m + 1):
            r = perm[k-1]
            q = math.sin(math.pi * r / p)**2 / math.sin(math.pi * k / p)**2
            Q_sin.append(q)

        # Verify against Fourier
        S = list(range(1, m + 1))
        Q_fourier = []
        for k in range(1, m + 1):
            val = sum(omega ** (k * s) for s in S)
            Q_fourier.append(abs(val)**2)

        print(f"    Q_sin    = {[f'{q:.6f}' for q in Q_sin]}")
        print(f"    Q_fourier = {[f'{q:.6f}' for q in Q_fourier]}")
        max_err = max(abs(Q_sin[i] - Q_fourier[i]) for i in range(m))
        print(f"    max_err = {max_err:.2e}")

        # The characteristic polynomial has integer coefficients:
        # prod(t - Q_k) = sum_{j=0}^m (-1)^j e_j t^{m-j}
        # Roots are sin^2(pi*r/p)/sin^2(pi*k/p) for k=1,...,m
        # This is related to the Chebyshev polynomial of the second kind!

        # Actually, U_{p-1}(cos(theta)) = sin(p*theta)/sin(theta)
        # So sin(p*theta)/sin(theta) = U_{p-1}(cos(theta))
        # And Q_k = [sin(pi*r/p)/sin(pi*k/p)]^2

        # But the permutation pi complicates things. Let me think about this differently.

    # Part 8: The generating function approach
    print("\n--- PART 8: GENERATING FUNCTION IDENTITY ---")
    print("prod_{k=1}^m (1 + t * Q_k) = 1 + sum e_j t^j")
    print("")
    print("Take log: sum_{k=1}^m log(1 + t*Q_k) = sum_{n>=1} (-1)^{n+1}/n * t^n * p_n")
    print("where p_n = sum Q_k^n is the n-th power sum.")
    print("")
    print("If p_n = m*(4m^2-1)*(16m^4-24m^2+9)*.../... then the gen function might factor.")

    # Actually, let me try: is the char poly related to Chebyshev polynomials?
    # The roots Q_k = sin^2(pi*pi_k/p)/sin^2(pi*k/p) where pi is a KNOWN permutation.
    # Let x_k = cos(2*pi*k/p) so sin^2(pi*k/p) = (1-x_k)/2.
    # Then Q_k = (1-x_{pi_k})/(1-x_k).

    print("\n  Chebyshev connection: Q_k = (1 - cos(2*pi*pi_k/p)) / (1 - cos(2*pi*k/p))")

    for p in [7, 11, 13]:
        m = (p - 1) // 2

        perm = []
        for k in range(1, m + 1):
            r = (k * (m + 1)) % p
            if r > m:
                r = p - r
            perm.append(r)

        # Compute (1-cos(2*pi*pi_k/p))/(1-cos(2*pi*k/p))
        Q_cheb = []
        for k in range(1, m + 1):
            r = perm[k-1]
            num = 1 - math.cos(2*math.pi*r/p)
            den = 1 - math.cos(2*math.pi*k/p)
            Q_cheb.append(num/den)

        e_vals, Q_vals = compute_esf(p)
        max_err = max(abs(Q_cheb[i] - Q_vals[i]) for i in range(m))
        print(f"  p={p}: max_err = {max_err:.2e}")

    # Part 9: Spread polynomial (Wildberger)
    print("\n--- PART 9: SPREAD POLYNOMIAL CONNECTION ---")
    print("The 'spread' s = sin^2(theta) satisfies the 'spread polynomial'")
    print("S_n(s) = sin^2(n*theta)/sin^2(theta) where s = sin^2(theta)")
    print("S_n is a degree n-1 polynomial in s with integer coefficients!")
    print("")
    print("Q_k = S_r(s_k) where s_k = sin^2(pi*k/p), r = pi(k)")
    print("Actually Q_k = sin^2(pi*r/p)/sin^2(pi*k/p) = S_r(s_k)... NO,")
    print("S_n(s) = sin^2(n*arcsin(sqrt(s)))/s, which is different.")
    print("")
    print("Better: let theta = pi*k/p. Then Q_k = sin^2(r*theta/k)/(sin^2(theta))...")
    print("This doesn't factor nicely since r/k is different for each k.")

    # Part 10: Direct formula search via OEIS-style
    print("\n--- PART 10: SEARCH FOR SEQUENCE IN TABLES ---")
    # Let me concatenate and search
    # p=11: [15, 35, 28, 9, 1]
    # OEIS: search "15, 35, 28, 9, 1" — this is C(6,2), ?, ?, ?, ?
    #
    # Actually, I notice: for p=13, the sequence is [21, 70, 84, 45, 11, 1]
    # C(7,2)=21, C(8,2)=28? No. C(7,2)=21, C(7,3)=35 not 70.
    # 70 = C(8,4)? Yes! C(8,4)=70. And 84 = C(9,3)? C(9,3)=84 YES!
    # 45 = C(10,2)? Yes! 11 = C(5,1) or just 11.
    #
    # Hmm: 21 = C(7,2), 70 = C(8,4), 84 = C(9,3), 45 = C(10,2), 11 = ?, 1 = C(?,0)
    # That's: C(2j+5, j+1) for j=0,1,2,3? C(5,1)=5 no.
    # Or: C(2m+1-2j, m-j) for j=0,...,m-1?
    # m=6: C(13-0,6)=C(13,6)=1716 NO.

    # Wait — for m=8 (p=17): [36, 210, 462, 495, 286, 91, 15, 1]
    # 36 = C(9,2), 210 = C(10,4), 462 = C(11,5)? C(11,5)=462 YES!
    # 495 = C(12,4)? C(12,4)=495 YES! 286 = C(13,3)? C(13,3)=286 YES!
    # 91 = C(14,2)? C(14,2)=91 YES! 15 = C(15,1)? C(15,1)=15 YES!
    # 1 = C(16,0) = 1.

    # Pattern: e_j = C(m+j, 2j) / something?
    # e_1 = C(m+1, 2) = C(9,2) = 36 YES
    # e_2 = C(m+2, 4)? C(10,4) = 210 YES
    # e_3 = C(m+3, 6)? C(11,6) = 462 YES!!
    # e_4 = C(m+4, 8)? C(12,8) = 495 YES!!!
    # e_5 = C(m+5, 10)? C(13,10) = 286 YES!!!!
    # e_6 = C(m+6, 12)? C(14,12) = 91 YES!!!!!
    # e_7 = C(m+7, 14)? C(15,14) = 15 YES!!!!!!
    # e_8 = C(m+8, 16)? C(16,16) = 1 YES!!!!!!!

    print("\n  TESTING: e_j = C(m+j, 2j)")
    for p, e_vals in all_esf.items():
        m = (p - 1) // 2
        match = True
        for j in range(1, m + 1):
            expected = math.comb(m + j, 2 * j)
            if expected != e_vals[j-1]:
                match = False
                break
        print(f"  p={p}, m={m}: e_j = C(m+j, 2j)? {'YES!!!' if match else 'NO'}")
        if match:
            print(f"    Verified: {[math.comb(m+j, 2*j) for j in range(1, m+1)]}")
            print(f"    Actual:   {e_vals}")


if __name__ == '__main__':
    main()
