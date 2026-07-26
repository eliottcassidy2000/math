#!/usr/bin/env python3
"""Exact finite gates for THM-2442.

The analytic BV estimate and the LRC identifications are proved in the
theorem text.  This companion independently exhausts the finite
reflection/visibility and cyclotomic gates.
"""

from fractions import Fraction
from itertools import product


P = 7
NONZERO = set(range(1, P))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def reflection_closed(mask: set[int]) -> bool:
    return all((-x) % P in mask for x in mask)


def phi7_remainder(coefficients: tuple[tuple[int, ...], ...]) -> tuple:
    """Remainder at zeta_7, coordinatewise over an arbitrary Q-basis.

    Since zeta_7^6=-(1+...+zeta_7^5), the remainder coefficients are
    a_i-a_6 for i=0,...,5.
    """

    last = coefficients[6]
    return tuple(
        tuple(coefficients[i][j] - last[j] for j in range(len(last)))
        for i in range(6)
    )


def main() -> None:
    # Reflection-paired nonflat phases versus one invisible word phase.
    reflected_masks = []
    visibility_checks = 0
    for bits in range(1, 1 << 6):
        mask = {i + 1 for i in range(6) if bits & (1 << i)}
        if not reflection_closed(mask):
            continue
        reflected_masks.append(mask)
        for missing in range(P):
            visible = set(range(P)) - {missing}
            require(bool(mask & visible), "reflection visibility failure")
            visibility_checks += 1

    require(len(reflected_masks) == 7, "reflection-mask count drift")
    require(visibility_checks == 49, "reflection-check count drift")

    # Exact mass gate: 0<=q_l<=tau and sum q_l=6 tau.
    mass_vectors = 0
    for tau in range(1, 6):
        for q in product(range(tau + 1), repeat=P):
            if sum(q) != 6 * tau:
                continue
            mass_vectors += 1
            require(
                sum(value == 0 for value in q) <= 1,
                "bounded mass vector has two zero phases",
            )
    require(mass_vectors == 791, "bounded mass-vector count drift")

    # Sharp hostile without reflection: the sole charged phase can equal
    # the sole missing word phase.
    unique_nonflat = {1}
    visible = set(range(P)) - {1}
    require(not (unique_nonflat & visible), "hostile unexpectedly visible")
    require(
        not reflection_closed(unique_nonflat),
        "singleton hostile unexpectedly reflection closed",
    )

    # Phi_7 kernel over Q, checked exhaustively.  The same coordinatewise
    # calculation applies over K=Q(zeta_13) in its Q-basis.
    scalar_kernel = 0
    for coeffs in product((-1, 0, 1), repeat=P):
        lifted = tuple((c,) for c in coeffs)
        if not any(any(row) for row in phi7_remainder(lifted)):
            scalar_kernel += 1
            require(len(set(coeffs)) == 1, "nonconstant scalar Phi7 kernel")
    require(scalar_kernel == 3, "scalar Phi7 kernel count drift")

    # A two-coordinate K-control checks that no hidden coordinate mixing
    # enters the reduction.
    vector_kernel = 0
    alphabet = ((0, 0), (1, 0), (0, 1), (1, 1))
    for coeffs in product(alphabet, repeat=P):
        if not any(any(row) for row in phi7_remainder(coeffs)):
            vector_kernel += 1
            require(len(set(coeffs)) == 1, "nonconstant vector Phi7 kernel")
    require(vector_kernel == 4, "vector Phi7 kernel count drift")

    # Positive rational control.  Phases +/-1 are nonflat; whichever
    # single word phase is missing, at least one survives.  Phase zero is
    # flat.  The thirteen-vectors are represented exactly.
    flat = (Fraction(1),) * 13
    bump = (Fraction(0),) + (Fraction(1),) * 12
    source_rows = [flat for _ in range(P)]
    source_rows[1] = bump
    source_rows[6] = tuple(bump[(-s) % 13] for s in range(13))
    controls = 0
    for missing in range(P):
        q = [Fraction(1) for _ in range(P)]
        q[missing] = Fraction(0)
        weighted = [
            tuple(q[ell] * x for x in source_rows[ell])
            for ell in range(P)
        ]
        require(
            any(weighted[ell] != weighted[0] for ell in (1, 6)),
            "positive control lost both reflected phases",
        )
        controls += 1
    require(controls == 7, "positive-control count drift")

    # Clock parity is only a permutation ell -> -ell.
    plus = tuple(range(P))
    minus = tuple((-ell) % P for ell in range(P))
    require(set(plus) == set(minus), "clock reflection is not a permutation")
    require(13 % 7 == -1 % 7, "13 is not -1 modulo seven")

    print("THM-2442 exact finite companion")
    print(f"reflection_closed_nonempty_phase_sets={len(reflected_masks)}")
    print(f"reflection_visibility_checks={visibility_checks}")
    print(f"bounded_mass_vectors_tau_1_to_5={mass_vectors}")
    print("unique_nonflat_without_reflection_hostile=PASS")
    print(f"phi7_scalar_kernel={scalar_kernel}")
    print(f"phi7_two_coordinate_kernel={vector_kernel}")
    print(f"positive_missing_phase_controls={controls}")
    print("clock_parity_permutation=PASS")
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
