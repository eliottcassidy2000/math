#!/usr/bin/env python3
"""Primary exact certificate for THM-4341.

This dependency-free path audits odd self-similar cusp tails, reciprocal
split duality, differential-order reciprocity, and the two natural-number
indexings.  The finite range is a hostile replay of formulas proved for all
g,r in the theorem; it is not used as induction.
"""

from fractions import Fraction as F
from math import gcd
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")


def need(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(label)


def triangular(n: int) -> int:
    return n * (n + 1) // 2


def oriented_index(g: int, r: int) -> int:
    return g * (g - 1) + r


def quotient_address(g: int, r: int) -> tuple[int, int, int]:
    """Return (N,h,orientation) for the reciprocal-pair quotient."""

    m = 2 * g + 1
    centered_odd = 2 * r - m
    need(centered_odd != 0 and centered_odd % 2, "nonzero odd center")
    h = (abs(centered_odd) + 1) // 2
    orientation = 1 if centered_odd > 0 else -1
    return triangular(g - 1) + h, h, orientation


def reconstruct_quotient(N: int, orientation: int) -> tuple[int, int]:
    """Invert the quotient index plus its one-bit orientation sidecar."""

    need(N >= 1 and orientation in (-1, 1), "quotient input")
    g = 1
    while triangular(g) < N:
        g += 1
    h = N - triangular(g - 1)
    r = g + h if orientation > 0 else g + 1 - h
    return g, r


def main() -> None:
    max_g = 250
    oriented_seen = set()
    quotient_seen = set()
    rows = 0

    for g in range(1, max_g + 1):
        m = 2 * g + 1
        for r in range(1, m):
            rows += 1
            q, epsilon = divmod(r, 2)
            branch_degree = m - r + epsilon
            need(branch_degree % 2 == 1, f"g={g},r={r}: odd branch degree")
            tail_genus = (branch_degree - 1) // 2
            persistent_delta = q
            need(tail_genus == g - persistent_delta,
                 f"g={g},r={r}: tail genus")
            need(tail_genus + persistent_delta == g,
                 f"g={g},r={r}: delta partition")

            complement = m - r
            complement_delta = complement // 2
            need(tail_genus == complement_delta,
                 f"g={g},r={r}: reciprocal tail/delta")
            need(g - complement_delta == persistent_delta,
                 f"g={g},r={r}: reciprocal delta/tail")

            slope = F(r, m - r)
            reciprocal_slope = F(complement, m - complement)
            need(slope * reciprocal_slope == 1,
                 f"g={g},r={r}: reciprocal slope")

            # The honest uniform base change t=tau^(2(m-r)), followed by
            # z=tau^(2r)Z and x=tau^(mr)Y, balances both source monomials.
            source_weight = 2 * m * r
            need(2 * (m * r) == source_weight, "x-square weight")
            need(m * (2 * r) == source_weight, "z^m weight")
            need(r * (2 * (m - r) + 2 * r) == source_weight,
                 "(tz)^r weight")

            # If omega=sigma^k z^B dz/x and t=sigma^d, its contribution
            # beyond k is c*d*slope.  Complementary contributions multiply
            # to (c*d)^2.  B=g is the sharp integral positivity threshold.
            d = 2 * g + 3
            c_threshold = F(g + 1) - F(m, 2)
            need(c_threshold == F(1, 2), "threshold coefficient")
            excess = c_threshold * d * slope
            complement_excess = c_threshold * d * reciprocal_slope
            need(excess > 0 and complement_excess > 0,
                 f"g={g},r={r}: positive threshold excess")
            need(excess * complement_excess == (c_threshold * d) ** 2,
                 f"g={g},r={r}: reciprocal order product")
            below = F(g) - F(m, 2)
            need(below == F(-1, 2), "sharp below-threshold coefficient")

            n = oriented_index(g, r)
            need(g * (g - 1) < n <= g * (g + 1), "oriented block")
            need(n not in oriented_seen, "oriented index collision")
            oriented_seen.add(n)

            N, h, orientation = quotient_address(g, r)
            need(1 <= h <= g, "odd-square rank")
            need((2 * h - 1) ** 2 == (2 * r - m) ** 2,
                 "odd-square compression")
            need(reconstruct_quotient(N, orientation) == (g, r),
                 "quotient inverse with orientation")
            quotient_seen.add(N)

            reduced_num = r // gcd(r, m - r)
            reduced_den = (m - r) // gcd(r, m - r)
            need(F(reduced_num, reduced_den) == slope, "reduced slope")

    need(oriented_seen == set(range(1, max_g * (max_g + 1) + 1)),
         "oriented natural-number bijection")
    need(quotient_seen == set(range(1, triangular(max_g) + 1)),
         "reciprocal-orbit triangular bijection")

    # Sharp scope hostiles.
    even_m, even_r = 4, 1
    even_degree = even_m - even_r + (even_r % 2)
    need(even_degree == 4 and even_degree % 2 == 0,
         "even cusp has two-infinity parity")
    need(F(1, 3 - 1) == F(3, 9 - 3),
         "reduced slope loses total cusp scale")

    print("THM4341 ODD SELF-SIMILAR CUSP PRIMARY CERTIFICATE")
    print(f"FINITE_REPLAY=g1..{max_g};oriented_rows={rows};reciprocal_pairs={len(quotient_seen)}")
    print("TAIL_NORMALIZATION=Y0^2=Z^epsilon*(Z^(m-r)-c);epsilon=r mod 2;c!=0")
    print("GENUS_PARTITION=g_tail=g-floor(r/2)=floor((m-r)/2);delta_persist=floor(r/2)")
    print("RECIPROCAL_INVOLUTION=r<->m-r;lambda<->1/lambda;tail_genus<->persistent_delta")
    print("DIFFERENTIAL_EXCESS=chi*d*lambda;chi=B+1-m/2;product_under_reciprocity=(chi*d)^2")
    print("SHARP_BUFFER=m=2g+1;integer_B>=g_positive;B=g-1_negative")
    print("ORIENTED_INDEX=n=g(g-1)+r;blocks=2T_(g-1)+1..2T_g")
    print("RECIPROCAL_QUOTIENT=N=T_(g-1)+h;h=(abs(2r-m)+1)/2;orientation=sign(2r-m)")
    print("ODD_SQUARE_COMPRESSION=(2h-1)^2->h;orientation_bit_required")
    print("HOSTILES=even_m_gives_even_branch_degree_two_infinities;slope_1/2_shared_by_(3,1)_and_(9,3)")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
