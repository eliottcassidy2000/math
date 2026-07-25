#!/usr/bin/env python3
"""Exact first-depth audit for THM-2198's scalar 5+3 transition.

At valuation profile (1,1,2), normalize the guard modulo 13^3.  The deepest
positive-valuation mask is empty on the primitive guard-danger annulus.  Fix
the two remaining depth-one masks, remove their union, and ask whether five
unit masks can cover the residual.

The program exhausts:

* all 78 depth-one sign classes modulo 13^2;
* all 78*79/2 = 3,081 unordered pairs, with repetition allowed;
* all 1,014 unit sign classes modulo 13^3.

For a fixed shallow pair, repeated unit sign classes do not enlarge a union.
After duplicates are removed, any five unit masks have union size at most the
sum of the five largest individual intersections with the residual.  The
audit computes exactly this deliberately generous conditional upper bound.

All predicates use integer cross-products.  Runtime checks are ordinary
exceptions and remain active under ``python -O``.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256


P = 13
P2 = P * P
N = P**3


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def norm_mod(value: int, modulus: int) -> int:
    residue = value % modulus
    return min(residue, modulus - residue)


def sign_classes(modulus: int) -> tuple[int, ...]:
    return tuple(
        value
        for value in range(1, (modulus + 1) // 2)
        if value % P != 0
    )


ANNULUS = tuple(
    z
    for z in range(N)
    if z % P != 0 and 7 * norm_mod(z, N) > N
)
POSITION = {z: index for index, z in enumerate(ANNULUS)}
FULL = (1 << len(ANNULUS)) - 1


def mask_from_reduced_coefficient(
    coefficient: int,
    modulus: int,
) -> int:
    """Danger mask on ANNULUS for a coefficient reduced to ``modulus``."""
    mask = 0
    for z in ANNULUS:
        if 14 * norm_mod(coefficient * z, modulus) < modulus:
            mask |= 1 << POSITION[z]
    return mask


UNIT_CLASSES = sign_classes(N)
SHALLOW_CLASSES = sign_classes(P2)

UNIT_MASKS = tuple(
    (coefficient, mask_from_reduced_coefficient(coefficient, N))
    for coefficient in UNIT_CLASSES
)
SHALLOW_MASKS = tuple(
    (coefficient, mask_from_reduced_coefficient(coefficient, P2))
    for coefficient in SHALLOW_CLASSES
)


def audit_reductions() -> int:
    """Check reduced predicates against the actual 13a and 169a masks."""
    checks = 0
    for coefficient, reduced_mask in SHALLOW_MASKS:
        direct_mask = mask_from_reduced_coefficient(P * coefficient, N)
        require(
            direct_mask == reduced_mask,
            f"depth-one reduction failed at {coefficient}",
        )
        checks += 1

    # A depth-two coefficient is 169a.  At a primitive z its value modulo
    # 13^3 is a nonzero thirteenth root, whose norm is at least 1/13.
    for coefficient in range(1, P):
        direct_mask = mask_from_reduced_coefficient(P2 * coefficient, N)
        require(
            direct_mask == 0,
            f"depth-two mask is not empty at {coefficient}",
        )
        checks += 1
    return checks


def main() -> None:
    require(N == 2197, "modulus changed")
    require(N % 7 != 0 and N % 14 != 0, "torsion endpoints appeared")
    require(len(ANNULUS) == 1450, "annulus size changed")
    require(len(UNIT_CLASSES) == 1014, "unit sign-class count changed")
    require(len(SHALLOW_CLASSES) == 78, "shallow sign-class count changed")

    reduction_checks = audit_reductions()
    require(reduction_checks == 90, "reduction check count changed")

    pair_count = 0
    margin_histogram: Counter[int] = Counter()
    minimum_margin: int | None = None
    maximum_margin: int | None = None
    worst_rows: list[
        tuple[
            int,
            int,
            int,
            int,
            tuple[tuple[int, int], ...],
        ]
    ] = []
    transcript = sha256()

    for left_index, (left, left_mask) in enumerate(SHALLOW_MASKS):
        for right, right_mask in SHALLOW_MASKS[left_index:]:
            pair_count += 1
            residual = FULL & ~(left_mask | right_mask)
            residual_size = residual.bit_count()

            intersections = sorted(
                (
                    ((unit_mask & residual).bit_count(), unit)
                    for unit, unit_mask in UNIT_MASKS
                ),
                key=lambda row: (-row[0], row[1]),
            )
            top_five = tuple(intersections[:5])
            conditional_cap = sum(size for size, _ in top_five)
            margin = residual_size - conditional_cap

            # The sum is an upper bound on the union: duplicates among the
            # five actual coefficients contribute no new mask, and distinct
            # masks satisfy |union| <= sum of their residual intersections.
            require(
                margin > 0,
                f"conditional five-mask capacity failed at {(left, right)}",
            )

            margin_histogram[margin] += 1
            maximum_margin = (
                margin
                if maximum_margin is None
                else max(maximum_margin, margin)
            )
            row = (
                left,
                right,
                residual_size,
                conditional_cap,
                top_five,
            )
            if minimum_margin is None or margin < minimum_margin:
                minimum_margin = margin
                worst_rows = [row]
            elif margin == minimum_margin:
                worst_rows.append(row)

            transcript.update(
                (
                    f"{left},{right},{residual_size},{conditional_cap},"
                    f"{margin}:"
                ).encode("ascii")
            )
            for size, unit in intersections:
                transcript.update(f"{unit}={size},".encode("ascii"))
            transcript.update(b"\n")

    require(pair_count == 3081, "shallow-pair count changed")
    require(minimum_margin == 86, "minimum margin changed")
    require(maximum_margin == 218, "maximum margin changed")
    expected_worst = [
        (
            14,
            46,
            1046,
            960,
            (
                (204, 183),
                (200, 799),
                (190, 599),
                (186, 1000),
                (180, 1007),
            ),
        )
    ]
    require(worst_rows == expected_worst, "worst row changed")

    digest = transcript.hexdigest()
    require(
        digest
        == "d31e5c874d8b5893ff33fc35095c18dcdf865c7f8296285c1b3441a8d8d679d9",
        "full conditional-capacity table changed",
    )

    print("THM-2198 scalar five-plus-three first-depth audit")
    print(
        f"N={N}; annulus={len(ANNULUS)}; "
        f"unit_classes={len(UNIT_CLASSES)}; "
        f"shallow_classes={len(SHALLOW_CLASSES)}"
    )
    print(
        f"shallow_pairs_with_repetition={pair_count}; "
        f"reduction_checks={reduction_checks}"
    )
    print(
        f"conditional_margin_range={minimum_margin}..{maximum_margin}; "
        f"worst_rows={worst_rows}"
    )
    print(
        f"margin_histogram_size={len(margin_histogram)}; "
        f"margin_histogram_total={sum(margin_histogram.values())}"
    )
    print(f"table_sha256={digest}")
    print("strict_torsion_witnesses_thicken_to_open_intervals=yes")
    print("PASS")


if __name__ == "__main__":
    main()
