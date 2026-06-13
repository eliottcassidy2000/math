#!/usr/bin/env python3
"""
resonance_fourier_bridge.py -- Connect resonance levels to Fourier structure

The resonance level q of a chord pair (a,b) is the minimum odd integer
such that qa = +/-b or a = +/-qb mod p.

CONJECTURE: The resonance level determines which Gauss sum phase
contributes to the Walsh coefficient, and this is why the two magnitude
classes exist.

Key relationships to test:
1. q determines the "distance" between gaps in the multiplicative group
2. For Paley, the Fourier coefficients split by chi(t) into two phases
3. The Walsh coefficient magnitude depends on how the phases combine

Author: kind-pasteur-2026-03-12-S60
"""

import cmath
import math
from itertools import combinations
from collections import defaultdict


def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def resonance_level(a, b, p):
    for q in range(1, p, 2):
        if (q * a - b) % p == 0 or (q * a + b) % p == 0:
            return q
        if q > 1 and ((a - q * b) % p == 0 or (a + q * b) % p == 0):
            return q
    return p


def multiplicative_distance(a, b, p):
    """Find the minimum generator power that maps a to +/-b in (Z/pZ)*."""
    # Find a primitive root
    for g in range(2, p):
        if pow(g, (p-1)//2, p) != 1:  # g is a primitive root candidate
            # Check all powers
            log_a = None
            log_b = None
            val = 1
            for k in range(p-1):
                if val == a % p:
                    log_a = k
                if val == b % p:
                    log_b = k
                val = (val * g) % p
            if log_a is not None and log_b is not None:
                diff = (log_b - log_a) % (p - 1)
                neg_diff = (log_b - log_a + (p-1)//2) % (p - 1)
                return min(diff, p - 1 - diff, neg_diff, p - 1 - neg_diff)
    return None


def main():
    for p in [7, 11, 13, 17, 19]:
        m = (p - 1) // 2
        print("=" * 70)
        print(f"RESONANCE-FOURIER BRIDGE at p={p}")
        print("=" * 70)

        # Gauss sum
        g = sum(legendre(t, p) * cmath.exp(2j * cmath.pi * t / p)
                for t in range(1, p))
        print(f"  Gauss sum: g = {g.real:.4f} + {g.imag:.4f}i")
        print(f"  |g|^2 = {abs(g)**2:.4f} (should be {p})")
        print(f"  g^2 = {(g**2).real:.4f} + {(g**2).imag:.4f}i "
              f"(should be {'-%d' % p if p%4==3 else '+%d' % p})")

        # Paley S_hat(t) = (-1 + chi(t)*g_1)/2 where g_1 = sum chi(s)*omega^s
        # For Paley, S = QR. hat(1_QR)(t) = sum_{s in QR} omega^{st}
        # = (1/2)(sum_{s=1}^{p-1} omega^{st} + sum_{s=1}^{p-1} chi(s)*omega^{st})
        # = (1/2)(-1 + chi(t)*g)  [using Gauss sum property]
        omega = cmath.exp(2j * cmath.pi / p)
        qr = [s for s in range(1, p) if legendre(s, p) == 1]

        S_hat = []
        for t in range(p):
            val = sum(omega**(t*s) for s in qr)
            S_hat.append(val)

        print(f"\n  S_hat(t) for Paley (QR={qr[:4]}...):")
        for t in range(min(p, 8)):
            chi_t = legendre(t, p)
            predicted = (-1 + chi_t * g) / 2 if t > 0 else complex(m, 0)
            match = abs(S_hat[t] - predicted) < 0.001
            print(f"    t={t}: S_hat={S_hat[t].real:>8.4f}+{S_hat[t].imag:>8.4f}i, "
                  f"chi(t)={chi_t:+d}, "
                  f"|S_hat|^2={abs(S_hat[t])**2:>8.4f}, "
                  f"pred_match={'Y' if match else 'N'}")

        # For general circulant tournament with S, S_hat(t) = sum_{s in S} omega^{st}
        # The key: when we flip chord j (gap g_j <-> p-g_j), the Fourier transform changes:
        # S_hat(t) -> S_hat(t) + (omega^{-g_j*t} - omega^{g_j*t}) = S_hat(t) - 2i*sin(2pi*g_j*t/p)
        # Wait, we replace g_j with p-g_j in S:
        # Delta = omega^{(p-g_j)*t} - omega^{g_j*t} = omega^{-g_j*t} - omega^{g_j*t} = -2i*sin(2pi*g_j*t/p)

        # The Walsh coefficient h_hat[{a,b}] measures the INTERACTION between flipping chord a and chord b.
        # In Fourier space, this is related to the product of sin terms.

        # For a chord pair with gaps (g_a, g_b):
        # The degree-2 Walsh of H involves sum over sigma of sigma_a*sigma_b*H(sigma)
        # In Fourier space, H depends on S_hat via the cycle counts and independence polynomial.

        # The resonance level q connects gaps g_a and g_b: q*g_a = +/-g_b mod p
        # This means that in Fourier space, the t-values where sin(2pi*g_a*t/p) and
        # sin(2pi*g_b*t/p) have a specific relationship determined by q.

        print(f"\n  Resonance structure:")
        pairs = [(a, b) for a in range(m) for b in range(a+1, m)]
        for a, b in pairs:
            ga, gb = a + 1, b + 1
            q = resonance_level(ga, gb, p)
            chi_q = legendre(q, p)
            chi_ab = legendre(ga * gb, p)

            # The key: what is chi(q)?
            # q*ga = +/-gb mod p, so chi(q)*chi(ga) = chi(+/-gb) = chi(gb) (since chi(-1)=chi(-1))
            # For p=3 mod 4: chi(-1) = -1, so:
            # If q*ga = gb: chi(q) = chi(gb)/chi(ga) = chi(gb*ga^{-1}) = chi(gb/ga)
            # If q*ga = -gb: chi(q) = chi(-gb)/chi(ga) = -chi(gb/ga)

            # Compute chi(gb/ga) = chi(gb * ga^{-1})
            ga_inv = pow(ga, p - 2, p)
            ratio = (gb * ga_inv) % p
            chi_ratio = legendre(ratio, p)

            # Check which condition holds
            if (q * ga - gb) % p == 0:
                rel = f"{q}*{ga}={gb}"
                chi_q_pred = chi_ratio
            elif (q * ga + gb) % p == 0:
                rel = f"{q}*{ga}=-{gb}"
                chi_q_pred = -chi_ratio if p % 4 == 3 else chi_ratio
            elif (ga - q * gb) % p == 0:
                rel = f"{ga}={q}*{gb}"
                chi_q_pred = legendre(pow(gb, p-2, p) * ga % p, p)
            elif (ga + q * gb) % p == 0:
                rel = f"{ga}=-{q}*{gb}"
                if p % 4 == 3:
                    chi_q_pred = -legendre(pow(gb, p-2, p) * ga % p, p)
                else:
                    chi_q_pred = legendre(pow(gb, p-2, p) * ga % p, p)
            else:
                rel = "?"
                chi_q_pred = 0

            print(f"    ({a},{b}) gaps=({ga},{gb}): q={q}, chi(q)={chi_q:+d}, "
                  f"chi(ab)={chi_ab:+d}, chi(b/a)={chi_ratio:+d}, "
                  f"rel={rel}")

        # At p=11: what separates q=3 from q=5?
        if p == 11:
            print(f"\n  q=3 vs q=5 at p=11:")
            for q_target in [3, 5]:
                print(f"\n    q={q_target} pairs:")
                for a, b in pairs:
                    ga, gb = a + 1, b + 1
                    q = resonance_level(ga, gb, p)
                    if q == q_target:
                        chi_ab = legendre(ga * gb, p)
                        # Sum of sin products
                        sin_prod_sum = 0
                        for t in range(1, p):
                            sa = math.sin(2 * math.pi * ga * t / p)
                            sb = math.sin(2 * math.pi * gb * t / p)
                            sin_prod_sum += sa * sb
                        print(f"      ({a},{b}) gaps=({ga},{gb}), chi(ab)={chi_ab:+d}, "
                              f"sum(sin*sin)={sin_prod_sum:.4f}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
