#!/usr/bin/env python3
"""
discriminant_analysis.py -- The Q-polynomial discriminant identity

DISCOVERY: The discriminant of the characteristic polynomial of
{Q_1,...,Q_m} for the Interval tournament equals p^{m(m-1)/2} EXACTLY.

Q_k = |S_hat(k)|^2 where S={1,...,m} is the Interval connection set.
The Q_k are roots of a degree-m polynomial with integer coefficients.
disc = prod_{i<j} (Q_i - Q_j)^2.

For Interval: disc = p^{C(m,2)} = p^{m(m-1)/2}.

p=5: m=2, disc = 5^1 = 5
p=7: m=3, disc = 7^3 = 343? No, we got 49 = 7^2 = 7^{3*2/2} = 7^3?
Wait: C(3,2)=3 but disc=49=7^2. Let me recheck.

Actually: disc(p=7) = 49 = 7^2, and C(3,2) = 3, so disc ≠ p^{C(m,2)}.
disc(p=11) = 14641 = 11^4, and C(5,2) = 10, so disc ≠ p^{C(m,2)} either.
disc(p=13) = 371293 = 13^5, and C(6,2) = 15.

The pattern is disc = p^{m-1}:
  p=7: 7^2 = 49, m-1=2. YES
  p=11: 11^4 = 14641, m-1=4. YES
  p=13: 13^5 = 371293, m-1=5. YES

Let me verify at p=5 and extend to larger primes.

Author: kind-pasteur-2026-03-12-S59c
"""

import cmath
import math
from itertools import combinations


def main():
    print("=" * 70)
    print("Q-POLYNOMIAL DISCRIMINANT ANALYSIS")
    print("=" * 70)

    for p in [5, 7, 11, 13, 17, 19, 23, 29, 31]:
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)

        # Interval: S = {1, ..., m}
        S = list(range(1, m + 1))

        Q_vals = []
        for k in range(1, m + 1):
            val = sum(omega ** (k * s) for s in S)
            Q_vals.append(abs(val)**2)

        # Discriminant = prod_{i<j} (Q_i - Q_j)^2
        disc = 1.0
        for i in range(m):
            for j in range(i + 1, m):
                disc *= (Q_vals[i] - Q_vals[j])**2

        disc_int = round(disc)

        # Factor out powers of p
        if disc_int > 0:
            n = disc_int
            p_power = 0
            while n % p == 0:
                n //= p
                p_power += 1
            remainder = n
        else:
            p_power = 0
            remainder = disc_int

        print(f"\n  p={p}, m={m}:")
        print(f"    disc = {disc_int} = p^{p_power} * {remainder}")
        print(f"    m-1 = {m-1}")
        print(f"    disc = p^(m-1)? {'YES' if p_power == m-1 and remainder == 1 else 'NO'}")

        # Also compute ESFs
        e_vals = []
        for j in range(1, m + 1):
            ej = sum(math.prod(Q_vals[i] for i in combo)
                    for combo in combinations(range(m), j))
            e_vals.append(round(ej))

        print(f"    e_j = {e_vals}")
        print(f"    prod Q_k = {e_vals[-1]}")

        # Check if e_j = C(2m, 2j) / C(m, j)? Or some binomial?
        print(f"    Binomial check: ", end="")
        binom_match = True
        for j in range(m):
            expected = math.comb(2*(m-j), 2) if j == 0 else None
            if j == 0:
                expected_c = math.comb(2*m, 2*1) // math.comb(m, 1) if m > 0 else 0
            # Try C(m+1, j+1)
            c_val = math.comb(m+1, j+1)
            print(f"e_{j+1}={e_vals[j]} vs C({m+1},{j+1})={c_val}", end="; ")
        print()

        # Actually, let me just check if e_j = C(m+1, j+1) for interval
        all_binom = True
        for j in range(m):
            if e_vals[j] != math.comb(m+1, j+1):
                all_binom = False
                break
        print(f"    e_j = C(m+1, j+1)? {'YES' if all_binom else 'NO'}")

    # Part 2: Analytical formula for discriminant of Interval Q-polynomial
    print("\n" + "=" * 70)
    print("PART 2: ANALYTICAL FORMULA FOR INTERVAL DISCRIMINANT")
    print("=" * 70)

    print("\nThe Q-polynomial for Interval has roots Q_k = sin^2(k*(m+1)*pi/p) / sin^2(k*pi/p)")
    print("where k=1,...,m and m=(p-1)/2, m+1=(p+1)/2.")
    print("\nSince {k*(m+1) mod p} is a permutation of {1,...,m} (up to sign),")
    print("Q_k = sin^2(pi*sigma(k)/p) / sin^2(pi*k/p) for some permutation sigma.")
    print("\ndisc = prod_{i<j} (Q_i - Q_j)^2")

    # Part 3: Non-Interval discriminants
    print("\n" + "=" * 70)
    print("PART 3: DISCRIMINANT FOR OTHER ORIENTATIONS")
    print("=" * 70)

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)

        pairs = [(s, p - s) for s in range(1, m + 1)]
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        disc_vals = {}
        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))

            Q_vals = []
            for k in range(1, m + 1):
                val = sum(omega ** (k * s) for s in S)
                Q_vals.append(abs(val)**2)

            disc = 1.0
            for i in range(m):
                for j in range(i + 1, m):
                    disc *= (Q_vals[i] - Q_vals[j])**2

            disc_int = round(disc)
            if disc_int not in disc_vals:
                disc_vals[disc_int] = 0
            disc_vals[disc_int] += 1

        print(f"\n  p={p}, m={m}:")
        for d, count in sorted(disc_vals.items()):
            if d > 0:
                n = d
                p_power = 0
                while n % p == 0:
                    n //= p
                    p_power += 1
                print(f"    disc = {d} = p^{p_power} * {n}, count = {count}")
            else:
                print(f"    disc = {d}, count = {count}")

    # Part 4: e_j = C(m+1, j+1) identity proof
    print("\n" + "=" * 70)
    print("PART 4: e_j = C(m+1, j+1) IDENTITY FOR INTERVAL")
    print("=" * 70)

    print("\nIf the char poly of Q_1,...,Q_m is prod_{k=1}^m (t - Q_k),")
    print("and e_j = C(m+1, j+1), then the char poly is:")
    print("  t^m - C(m+1,2)*t^{m-1} + C(m+1,3)*t^{m-2} - ... + (-1)^m * C(m+1,m+1)")
    print("  = t^m - C(m+1,2)*t^{m-1} + ... + (-1)^m * 1")
    print("\nThis is the Chebyshev-like polynomial arising from sin^2 ratios!")

    for p in [5, 7, 11, 13, 17]:
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

        binom_vals = [math.comb(m+1, j+1) for j in range(m)]

        match = (e_vals == binom_vals)
        print(f"\n  p={p}, m={m}:")
        print(f"    e_j    = {e_vals}")
        print(f"    C(m+1) = {binom_vals}")
        print(f"    Match: {'YES' if match else 'NO'}")

        if match:
            # The char poly is sum_{j=0}^m (-1)^j C(m+1,j+1) t^{m-j}
            # = (1/t) * sum_{j=0}^m (-1)^j C(m+1,j+1) t^{m+1-j}
            # Note: sum_{j=0}^m (-1)^j C(m+1,j+1) = sum_{k=1}^{m+1} (-1)^{k-1} C(m+1,k) = 1
            print(f"    Char poly: sum_j (-1)^j * C({m+1},j+1) * t^({m}-j)")
            print(f"    = [sum_k (-1)^{{k-1}} C({m+1},k) t^{{m+1-k}}] / t")
            print(f"    = [t^{{m+1}} - (1-t)^{{m+1}}] / t")
            print(f"    = [t^{m+1} - (1-t)^{m+1}] / t")

            # Verify: roots of t^{m+1} - (1-t)^{m+1} = 0
            # => (t/(1-t))^{m+1} = 1
            # => t/(1-t) = omega_k where omega_k = exp(2pi*i*k/(m+1))
            # => t = omega_k / (1 + omega_k) = 1/(1 + omega_k^{-1})
            print(f"\n    Roots of t^{{{m+1}}} = (1-t)^{{{m+1}}}:")
            print(f"    => (t/(1-t))^{{{m+1}}} = 1")
            print(f"    => t/(1-t) = exp(2*pi*i*k/{m+1}), k=1,...,{m}")
            print(f"    => t_k = 1/(1 + exp(-2*pi*i*k/{m+1}))")
            print(f"         = 1/2 + i*cot(pi*k/{m+1})/2...wait")

            # Actually: t/(1-t) = e^{i*theta} where theta = 2*pi*k/(m+1)
            # t = e^{i*theta} / (1 + e^{i*theta}) = 1/(1 + e^{-i*theta})
            # = e^{i*theta/2} / (e^{-i*theta/2} + e^{i*theta/2}) = e^{i*theta/2} / (2*cos(theta/2))
            # |t|^2 = 1 / (4*cos^2(theta/2)) = 1/(2*(1+cos(theta))) = 1/(2+2*cos(2*pi*k/(m+1)))

            # Hmm, but our Q_k should be REAL. Let's compute:
            roots = []
            for k in range(1, m + 1):
                theta = 2 * math.pi * k / (m + 1)
                # t/(1-t) = e^{i*theta}
                # t = e^{i*theta}/(1 + e^{i*theta}) = 1/(1 + e^{-i*theta})
                z = cmath.exp(1j * theta)
                t = z / (1 + z)
                roots.append(t)

            # These are complex. But the polynomial has real coefficients.
            # The roots come in conjugate pairs. The Q_k are the actual roots
            # which are REAL.
            # Let me re-derive. poly(t) = [t^{m+1} - (1-t)^{m+1}] / t
            # For t real, this equals 0 when t^{m+1} = (1-t)^{m+1}
            # => |t/(1-t)| = 1 => t and 1-t have same absolute value => t = 1/2
            # But that gives only 1 real root! Unless m+1 is even.
            # For m+1 even (p=3 mod 4): t=1/2 is a root (double?)
            # For m+1 odd (p=1 mod 4): t^{m+1}=(1-t)^{m+1} => t=1-t => t=1/2 is the only real solution

            # Wait, our Q_k ARE real but they're NOT roots of this polynomial??
            # Let me recheck.
            poly_vals = []
            for qk in Q_vals:
                pv = qk**(m+1) - (1-qk)**(m+1)
                poly_vals.append(pv / qk if abs(qk) > 1e-10 else float('inf'))
            print(f"\n    Poly values at Q_k (should be ~0):")
            for k, pv in enumerate(poly_vals):
                print(f"      Q_{k+1} = {Q_vals[k]:.8f}, poly/t = {pv:.6e}")


if __name__ == '__main__':
    main()
