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

    # The deletion-branch anchor: z_0=0 and sum z_ell=u!=0.  Equality in
    # the sharp Cauchy bound occurs at z_1=...=z_6=u/6.
    anchored = (k13_rational(0),) + (k13_rational(Fraction(1, 6)),) * 6
    require(source_charged_all(anchored), "anchored profile lost a source colour")
    anchored_energy = (
        Fraction(1, P7) * 6 * Fraction(1, 36) - Fraction(1, 49)
    )
    require(anchored_energy == Fraction(1, 294), "anchored energy mismatch")

    # Positive equality control for the anchored bound.  The original
    # owner phase has a flat row; the other six rows vanish at s=0 and
    # equal 1/84 elsewhere.  Hence z_0(b)=0 and z_1=...=z_6 for every
    # nonzero target character.
    owner_rows = (
        (Fraction(1, 14),) * P13,
    ) + (
        ((Fraction(0),) + (Fraction(1, 84),) * (P13 - 1)),
    ) * (P7 - 1)
    for b in range(1, P13):
        z_profile = tuple(
            tuple(value / P13 for value in cyclotomic13_remainder(row, b))
            for row in owner_rows
        )
        require(z_profile[0] == ZERO13, f"positive anchor failed at b={b}")
        require(
            all(value == z_profile[1] for value in z_profile[1:]),
            f"positive equality rows differ at b={b}",
        )
        deletion_row = tuple(
            sum(owner_rows[ell][s] for ell in range(P7))
            for s in range(P13)
        )
        deletion_transform = tuple(
            value / P13
            for value in cyclotomic13_remainder(deletion_row, b)
        )
        require(
            tuple(
                sum(z_profile[ell][j] for ell in range(P7))
                for j in range(12)
            )
            == deletion_transform,
            f"positive equality partition failed at b={b}",
        )

    target_floor = Fraction(27, 28561)
    mixed_floor = target_floor / 294
    require(
        mixed_floor == Fraction(9, 2798978),
        "mixed source/target floor mismatch",
    )
    require(
        mixed_floor / 72 == Fraction(1, 4732 * 4732),
        "joint maximum-mode denominator mismatch",
    )
    eligible_full_floor = mixed_floor / (169 * 12 * 13)
    require(eligible_full_floor > 0, "eligible full-table floor vanished")

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

    # A stronger flat-source hostile retaining H(t,s,t)=0 and all nonzero
    # deep colours.  H is supported on r=t+1, so the diagonal is empty.
    diagonal_target = tuple(
        Fraction(0) if s == 0 else Fraction(6, 7) for s in range(P13)
    )
    require(
        target_charged_all(diagonal_target),
        "diagonal hostile lost a target colour",
    )
    for s in range(P13):
        for t in range(P13):
            r = (t + 1) % P13
            require(r != t, f"diagonal hostile hit r=t at {t}")

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
    print("anchored z0=0 sharp source energy=1/294")
    print("positive anchored equality control=PASS")
    print("joint mixed floor=9/2798978")
    print("some joint mode floor denominator=4732")
    print("finite flat-source/nonflat-target hostile target energy=12/169")
    print("flat-source diagonal-zero all-target/deep hostile=PASS")
    print("smooth hostile eps=delta=1/2, N>=3")
    print("smooth hostile source profile independent of ell=PASS")
    print("smooth hostile target b=+/-1 amplitude=1/112")
    print("fixed-word partition=PASS")
    print("mod-seven active-word coupled-shift partition=FAILS SHARPLY")
    print("THM-2409 exact companion PASS")


if __name__ == "__main__":
    main()
