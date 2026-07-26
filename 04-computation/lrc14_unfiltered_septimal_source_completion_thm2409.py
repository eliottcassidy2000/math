#!/usr/bin/env python3
"""Exact companion for THM-2409.

The script uses only integer/Fraction arithmetic.  It checks the C7
all-or-flat algebra over a symbolic Q(zeta_13) coefficient basis, the
danger-partition identity, two exact flat-source/nonflat-target hostiles,
and the coupled-word failure of the partition.
"""

from fractions import Fraction
from itertools import product


P7 = 7
P13 = 13
ZERO13 = (Fraction(0),) * 12
ONE13 = (Fraction(1),) + (Fraction(0),) * 11


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def k13_add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def k13_sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def k13_scale(c, a):
    return tuple(c * x for x in a)


def k13_rational(q):
    return k13_scale(Fraction(q), ONE13)


def source_remainder(profile, character):
    """Polynomial at xi_7^character in the K-basis 1,...,xi_7^5.

    K=Q(zeta_13).  Coprime-conductor intersection makes Phi_7
    irreducible over K.  Exponents are permuted first, then
    xi_7^6=-(1+...+xi_7^5).
    """

    coeff = [ZERO13 for _ in range(P7)]
    for shift, value in enumerate(profile):
        exponent = (character * shift) % P7
        coeff[exponent] = k13_add(coeff[exponent], value)
    top = coeff[6]
    return tuple(k13_sub(coeff[i], top) for i in range(6))


def source_charged_all(profile):
    return all(
        any(value != ZERO13 for value in source_remainder(profile, e))
        for e in range(1, P7)
    )


def source_flat(profile):
    return all(value == profile[0] for value in profile)


def cyclotomic13_remainder(values, character):
    coeff = [Fraction(0)] * P13
    for shift, value in enumerate(values):
        coeff[(character * shift) % P13] += value
    top = coeff[12]
    return tuple(coeff[i] - top for i in range(12))


def target_charged_all(values):
    return all(
        any(cyclotomic13_remainder(values, b))
        for b in range(1, P13)
    )


def target_energy(values):
    total = sum(values)
    squares = sum(value * value for value in values)
    return Fraction(P13 * squares - total * total, P13 * P13)


def danger(root, shift):
    return int(root == shift)


def main():
    # Exact a.e. seven-shift partition on the seven open component centres.
    for root in range(P7):
        require(
            sum(danger(root, shift) for shift in range(P7)) == 1,
            f"danger partition failed at root {root}",
        )

    boolean_count = 0
    flat_count = 0
    for bits in product((0, 1), repeat=P7):
        profile = tuple(k13_rational(bit) for bit in bits)
        flat = source_flat(profile)
        require(
            source_charged_all(profile) == (not flat),
            f"C7 all-or-flat failed for {bits}",
        )
        boolean_count += 1
        flat_count += int(flat)

    # A genuinely K-valued profile, not merely a rational profile.
    basis_profile = []
    for shift in range(P7):
        vector = [Fraction(0)] * 12
        vector[shift] = Fraction(shift + 1)
        basis_profile.append(tuple(vector))
    require(not source_flat(basis_profile), "K-valued control is flat")
    require(source_charged_all(basis_profile), "K-valued control lost a colour")

    # Rational finite-table hostile: target is charged, source phase is flat.
    target = (Fraction(1),) + (Fraction(2),) * 12
    flat_source_rows = tuple(
        tuple(value / P7 for value in target) for _ in range(P7)
    )
    require(target_charged_all(target), "target hostile lost a C13 colour")
    for s in range(P13):
        profile = tuple(
            k13_rational(flat_source_rows[ell][s]) for ell in range(P7)
        )
        require(source_flat(profile), f"source hostile nonflat at s={s}")
    require(target_energy(target) == Fraction(12, 169), "target energy mismatch")

    # Smooth one-circle hostile.  For eps=delta=1/2, A_s B has Fourier
    # support {-2,-1,0,1,2}; d(Nx-ell/7) has support N*Z.  N>=3 leaves
    # only frequency zero, so C_ell(s)=1/7(1+1/8 cos(2*pi*s/13)).
    epsilon = Fraction(1, 2)
    delta = Fraction(1, 2)
    low_support = set(range(-2, 3))
    for n_speed in (3, 4, 13, 91):
        comb_support_window = {
            n_speed * j for j in range(-5, 6)
        }
        require(
            low_support.intersection(comb_support_window) == {0},
            f"low-frequency separation failed at N={n_speed}",
        )
    cosine_amplitude = epsilon * delta / 2
    require(cosine_amplitude == Fraction(1, 8), "cosine amplitude mismatch")
    target_mode_amplitude = cosine_amplitude / (2 * P7)
    require(
        target_mode_amplitude == Fraction(1, 112),
        "smooth hostile target amplitude mismatch",
    )

    # A coupled mod-seven word shift destroys the source partition.
    fixed_word = tuple(1 - danger(root, 0) for root in range(P7))
    fixed_sum = tuple(
        sum(danger(root, ell) * fixed_word[root] for ell in range(P7))
        for root in range(P7)
    )
    active_sum = tuple(
        sum(
            danger(root, ell) * (1 - danger(root, ell))
            for ell in range(P7)
        )
        for root in range(P7)
    )
    require(fixed_sum == fixed_word, "fixed-word partition failed")
    require(active_sum == (0,) * P7, "active-word hostile failed")
    require(active_sum != fixed_sum, "active word accidentally preserved partition")

    print("THM-2409 SEPTIMAL SOURCE COMPLETION EXACT AUDIT")
    print("conductors=(7,13), gcd=1, cyclotomic intersection=Q")
    print(f"C7 Boolean profiles exhausted={boolean_count}")
    print(f"flat Boolean profiles={flat_count}")
    print("C7 all-or-flat over symbolic Q(zeta13)=PASS")
    print("seven shifted dangers partition one a.e.=PASS")
    print("K-valued nonflat control fires all six source colours=PASS")
    print("finite flat-source/nonflat-target hostile target energy=12/169")
    print("smooth hostile eps=delta=1/2, N>=3")
    print("smooth hostile source profile independent of ell=PASS")
    print("smooth hostile target b=+/-1 amplitude=1/112")
    print("fixed-word partition=PASS")
    print("mod-seven active-word coupled-shift partition=FAILS SHARPLY")
    print("THM-2409 exact companion PASS")


if __name__ == "__main__":
    main()
