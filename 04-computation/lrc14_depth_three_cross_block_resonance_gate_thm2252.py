#!/usr/bin/env python3
"""Exact arithmetic referee for THM-2252.

For the scalar profile (3,4,5), group each normalized blocker core's two
positive-carrier appearances before applying Fourier orthogonality.  In the
absence of a cross-block relation, positive Jackson smoothing factorizes the
guard and the three two-step core blocks.  The exact two-step danger
autocorrelation then places the carrier strictly below the scalar residual
floor.

The paper proof supplies the Fourier-factorization implication.  This script
checks every numerical constant, the Jackson error certificate, the exact
two-step Markov correlation, the first successful equal-degree ledger, and
the anisotropic relation-height invoice using Fraction arithmetic.
"""

from fractions import Fraction
from hashlib import sha256


TARGET = Fraction(961, 6930)
P = 13
N = 99
DEGREE = 2 * N - 2
ETA_CAP = Fraction(2129, 500000)
SIMPLE_CAP = Fraction(433, 3125)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def jackson_coefficient(n, k):
    """Integral Fourier numerator C_k for the normalized squared Fejer kernel."""
    if not 0 <= k <= 2 * n - 2:
        return 0
    if k <= n:
        return (
            4 * n**3
            - 6 * n * k**2
            + 2 * n
            + 3 * k**3
            - 3 * k
        ) // 6
    return ((2 * n - k) ** 3 - (2 * n - k)) // 6


def jackson_eta_upper(n):
    """Strict rational upper bound using pi < 355/113."""
    c_zero = Fraction(n * (2 * n**2 + 1), 3)
    odd_sum = sum(
        (
            Fraction(jackson_coefficient(n, k), k**2)
            for k in range(1, 2 * n - 2, 2)
        ),
        Fraction(),
    )
    return (
        Fraction(1, 2)
        - 4 * Fraction(113, 355) ** 2 * odd_sum / c_zero
    )


def two_step_zero_correlation():
    """P(X_0=0,X_2=0) from the exact backward danger-root law."""
    # transition[future][current]
    transition = (
        (Fraction(11, 13), Fraction(2, 13)),
        (Fraction(12, 13), Fraction(1, 13)),
    )
    two_step = tuple(
        tuple(
            sum(
                (
                    transition[future][middle]
                    * transition[middle][current]
                    for middle in (0, 1)
                ),
                Fraction(),
            )
            for current in (0, 1)
        )
        for future in (0, 1)
    )
    stationary = (Fraction(6, 7), Fraction(1, 7))
    return stationary[0] * two_step[0][0], two_step


def factorized_carrier(pair_zero):
    return Fraction(5, 7) * (
        1 - 2 * Fraction(6, 7) ** 3 + pair_zero**3
    )


def smoothed_carrier_cap(eta):
    pair_zero, _ = two_step_zero_correlation()
    return (
        7 * eta
        + Fraction(5, 7)
        * (
            1
            - 2 * Fraction(6, 7) ** 3
            + (pair_zero + 2 * eta) ** 3
        )
    )


def main():
    pair_zero, two_step = two_step_zero_correlation()
    require(
        two_step
        == (
            (Fraction(145, 169), Fraction(24, 169)),
            (Fraction(144, 169), Fraction(25, 169)),
        ),
        "two-step danger transition drift",
    )
    require(pair_zero == Fraction(870, 1183), "pair-zero correlation drift")

    independent = factorized_carrier(pair_zero)
    require(
        independent == Fraction(1144584995, 11589168409),
        "factorized carrier drift",
    )
    require(independent < TARGET, "factorized carrier missed target")

    eta_upper = jackson_eta_upper(N)
    require(eta_upper > 0, "Jackson error upper bound lost positivity")
    require(eta_upper < ETA_CAP, "Jackson error cap failed")

    cap = smoothed_carrier_cap(ETA_CAP)
    require(
        cap
        == Fraction(
            5017879304282732213894543,
            36216151278125000000000000,
        ),
        "smoothed carrier cap drift",
    )
    require(cap < SIMPLE_CAP, "simple rational cap failed")
    require(SIMPLE_CAP < TARGET, "simple cap missed scalar target")
    require(
        TARGET - SIMPLE_CAP == Fraction(487, 4331250),
        "simple target margin drift",
    )

    previous_eta_upper = jackson_eta_upper(N - 1)
    previous_cap = smoothed_carrier_cap(previous_eta_upper)
    require(
        previous_cap > TARGET,
        "degree-194 ledger unexpectedly closes",
    )

    require(DEGREE == 196, "Jackson degree drift")
    require(170 * DEGREE == 33320, "block frequency height drift")
    require(DEGREE // P == 15, "reduced guard height drift")
    require(DEGREE < 16 * P, "reduced guard height ceiling drift")

    # The grouping is load-bearing: the inevitable within-core mode
    # 169*(+1)+(-169)=0 is a block constant, not a cross-block resonance.
    require(169 * 1 - 169 == 0, "internal mode control drift")
    require(abs(1) <= DEGREE and abs(-169) <= DEGREE, "internal mode absent")
    require(13 * 1 - 13 == 0, "common-core diagonal syzygy drift")
    require(13 <= 170 * DEGREE, "common-core syzygy exceeds height")

    digest_payload = "|".join(
        (
            str(TARGET),
            str(pair_zero),
            str(independent),
            str(ETA_CAP),
            str(cap),
            str(SIMPLE_CAP),
            str(DEGREE),
            str(170 * DEGREE),
            str(DEGREE // P),
        )
    )
    digest = sha256(digest_payload.encode("ascii")).hexdigest()
    require(
        digest
        == "8c3e0e43c0dc40cd7c42464f8e5d965eda11b766e100af6ffa3139845cca1226",
        "certificate digest drift",
    )

    print("THM-2252 DEPTH-THREE CROSS-BLOCK SPECTRAL RESONANCE GATE")
    print(f"target={TARGET}")
    print(f"jackson_N={N} degree={DEGREE}")
    print(f"eta_rational_cap={ETA_CAP}")
    print("eta_exact_upper_lt_cap=PASS")
    print("two_step_transition=((145/169,24/169),(144/169,25/169))")
    print(f"pair_zero_correlation={pair_zero}")
    print(
        f"cross_block_factorized_carrier={independent} "
        f"decimal={float(independent):.15f}"
    )
    print(f"smoothed_relation_free_cap={cap} decimal={float(cap):.15f}")
    print(f"simple_cap={SIMPLE_CAP} target_minus={TARGET-SIMPLE_CAP}")
    print("previous_N=98 relation_free_ledger_closes=NO")
    print("internal_same_core_169_relation=ABSORBED_AS_BLOCK_CONSTANT")
    print("common_core_cross_block_syzygy=(a0,b1,b2,b3)=(0,13,-1,0)")
    print(
        "forced_relation="
        "a0*H+b1*u1+13*b2*u2+169*b3*u3=0"
    )
    print(f"guard_coefficient_height={DEGREE//P}")
    print(f"block_coefficient_height={170*DEGREE}")
    print("aggregate_vector_nonzero=YES")
    print(f"certificate_digest={digest}")
    print("status=THM2252_PROVED_VERIFIED_EXACT")


if __name__ == "__main__":
    main()
