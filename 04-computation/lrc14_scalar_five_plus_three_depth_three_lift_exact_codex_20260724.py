#!/usr/bin/env python3
"""Exact finite-lift audit for the scalar valuation profile (1,1,3).

The primitive unit group modulo 13^m, quotiented by sign, is cyclic.  At
m=4 the two valuation-one blockers are pullbacks from the quotient at m=3.
Over each quotient phase there are thirteen rooted sheets.  This script
stores, for every unit mask and quotient phase, the *integer* number (0, 1,
or 2) of its endpoints that survive the guard.  Conditional capacities for
all 514,605 shallow pairs then follow by exact inclusion--exclusion.

There is no floating point and no FFT in the certificate.  NumPy is used
only for fixed-width nonnegative integer arrays, integer matrix products,
integer sums, and selection of the five largest integer entries.  Every
possible accumulator is bounded explicitly below its dtype limit.

The m=3 computation is rerun as a positive control and must reproduce
THM-2198's margin 86 before the new m=4 audit is accepted.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from hashlib import sha256
from struct import pack

import numpy as np


P = 13
GENERATOR = 2


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def norm_mod(value: int, modulus: int) -> int:
    residue = value % modulus
    return min(residue, modulus - residue)


def primitive_root_audit(modulus: int) -> None:
    """Prove that 2 generates the unit group modulo the given 13-power."""
    phi = 12 * (modulus // P)
    require(
        all(pow(GENERATOR, phi // prime, modulus) != 1 for prime in (2, 3, 13)),
        f"{GENERATOR} is not primitive modulo {modulus}",
    )
    require(pow(GENERATOR, phi, modulus) == 1, "bad Euler order")


def sign_label(exponent: int, modulus: int) -> int:
    value = pow(GENERATOR, exponent, modulus)
    return min(value, modulus - value)


@dataclass(frozen=True)
class AuditResult:
    exponent: int
    modulus: int
    quotient: int
    group_size: int
    quotient_group_size: int
    universe_half: int
    fibre_histogram: tuple[tuple[int, int], ...]
    shallow_mask_size: int
    pair_count: int
    minimum_half_margin: int
    worst_pair_exponents: tuple[int, int]
    worst_pair_labels: tuple[int, int]
    worst_residual_half: int
    worst_top_half: tuple[int, ...]
    worst_unit_exponents: tuple[int, ...]
    worst_unit_labels: tuple[int, ...]
    root_matrix_sha256: str
    loss_matrix_sha256: str
    pair_summary_sha256: str
    margin_histogram_size: int


def audit_depth(exponent: int) -> AuditResult:
    """Run the exact quotient-fibre capacity audit at N=13^exponent."""
    require(exponent >= 3, "the lift audit starts at depth three")
    modulus = P**exponent
    quotient = modulus // P
    group_size = 6 * P ** (exponent - 1)
    quotient_group_size = group_size // P
    primitive_root_audit(modulus)
    primitive_root_audit(quotient)

    # Exponents 0,...,group_size-1 represent the unit classes modulo sign.
    values = np.fromiter(
        (pow(GENERATOR, e, modulus) for e in range(group_size)),
        dtype=np.int64,
        count=group_size,
    )
    norms = np.minimum(values, modulus - values)
    guard_safe = 7 * norms > modulus
    terminal_danger = 14 * norms < modulus

    quotient_values = np.fromiter(
        (
            pow(GENERATOR, e, quotient)
            for e in range(quotient_group_size)
        ),
        dtype=np.int64,
        count=quotient_group_size,
    )
    quotient_norms = np.minimum(
        quotient_values,
        quotient - quotient_values,
    )
    quotient_danger = 14 * quotient_norms < quotient

    # Reduction modulo 13^(m-1) is exponent reduction modulo |G_(m-1)|.
    # The 13 entries r+j|G_(m-1)| are precisely the rooted lift fibre.
    fibre_heights = guard_safe.reshape(
        P,
        quotient_group_size,
    ).sum(axis=0, dtype=np.uint8)
    require(
        int(fibre_heights.min()) == 9
        and int(fibre_heights.max()) == 10,
        "guard fibre is not a 9/10-sheet block",
    )
    fibre_histogram_counter = Counter(map(int, fibre_heights))
    fibre_histogram = tuple(sorted(fibre_histogram_counter.items()))

    # W[q,r] is the number of q-mask sheets inside the guard-safe part of
    # the fibre over r.  The root-sheet law predicts W in {0,1,2}.
    root_matrix = np.empty(
        (group_size, quotient_group_size),
        dtype=np.uint8,
    )
    exponent_axis = np.arange(group_size)
    batch_size = 128
    for lower in range(0, group_size, batch_size):
        unit_exponents = np.arange(
            lower,
            min(lower + batch_size, group_size),
        )
        batch = (
            guard_safe[None, :]
            & terminal_danger[
                (exponent_axis[None, :] + unit_exponents[:, None])
                % group_size
            ]
        )
        root_matrix[lower : lower + len(unit_exponents)] = batch.reshape(
            len(unit_exponents),
            P,
            quotient_group_size,
        ).sum(axis=1, dtype=np.uint8)
    require(
        int(root_matrix.min()) == 0
        and int(root_matrix.max()) == 2,
        "unit root-sheet entries escaped {0,1,2}",
    )
    root_matrix_transpose = np.ascontiguousarray(root_matrix.T)
    full_capacities = root_matrix.sum(axis=1, dtype=np.uint16)
    root_matrix_digest = sha256(root_matrix.tobytes()).hexdigest()

    # A[a,r] says that the valuation-one blocker with quotient unit part
    # 2^a is active on the entire fibre over r.
    quotient_axis = np.arange(quotient_group_size)
    shallow = np.empty(
        (quotient_group_size, quotient_group_size),
        dtype=np.uint8,
    )
    for a in range(quotient_group_size):
        shallow[a] = quotient_danger[
            (quotient_axis + a) % quotient_group_size
        ]
    shallow_size = int(shallow[0].sum())
    require(
        np.all(shallow.sum(axis=1, dtype=np.uint16) == shallow_size),
        "shallow translates changed size",
    )

    # The largest possible loss entry is 2*shallow_size, which is 288 at
    # m=4.  Thus the uint16 matrix product is exact without overflow.
    require(2 * quotient_group_size < 2**16, "uint16 capacity bound failed")
    loss_matrix = (
        shallow.astype(np.uint16)
        @ root_matrix_transpose.astype(np.uint16)
    )
    loss_matrix_digest = sha256(loss_matrix.tobytes()).hexdigest()

    shallow_height = (
        shallow.astype(np.uint16)
        * fibre_heights[None, :].astype(np.uint16)
    ).sum(axis=1, dtype=np.uint16)
    universe_half = int(fibre_heights.sum())
    full_capacities_i32 = full_capacities.astype(np.int32)
    loss_i32 = loss_matrix.astype(np.int32)

    pair_digest = sha256()
    margin_histogram: Counter[int] = Counter()
    minimum_margin: int | None = None
    worst_pairs: list[
        tuple[
            int,
            int,
            int,
            tuple[int, ...],
        ]
    ] = []
    pair_count = 0

    for left in range(quotient_group_size):
        left_mask = shallow[left].astype(bool)
        for right in range(left, quotient_group_size):
            pair_count += 1
            intersection_indices = np.flatnonzero(
                left_mask & shallow[right].astype(bool)
            )

            # Exact inclusion--exclusion.  The last term restores the
            # fibres removed by both shallow masks.
            intersection_capacity = root_matrix_transpose[
                intersection_indices
            ].sum(axis=0, dtype=np.uint16)
            capacities = (
                full_capacities_i32
                - loss_i32[left]
                - loss_i32[right]
                + intersection_capacity.astype(np.int32)
            )
            require(
                int(capacities.min()) >= 0,
                "a conditional capacity became negative",
            )

            top_five = np.partition(capacities, -5)[-5:]
            top_five_sorted = tuple(
                sorted(map(int, top_five), reverse=True)
            )
            top_five_sum = sum(top_five_sorted)
            residual = (
                universe_half
                - int(shallow_height[left])
                - int(shallow_height[right])
                + int(fibre_heights[intersection_indices].sum())
            )
            margin = residual - top_five_sum
            require(margin > 0, f"capacity failed at {(left, right)}")

            margin_histogram[margin] += 1
            pair_digest.update(
                pack(
                    "<9I",
                    left,
                    right,
                    residual,
                    *top_five_sorted,
                    margin,
                )
            )

            record = (
                left,
                right,
                residual,
                top_five_sorted,
            )
            if minimum_margin is None or margin < minimum_margin:
                minimum_margin = margin
                worst_pairs = [record]
            elif margin == minimum_margin:
                worst_pairs.append(record)

    require(
        pair_count
        == quotient_group_size * (quotient_group_size + 1) // 2,
        "bad shallow-pair count",
    )
    require(minimum_margin is not None, "no shallow pairs were audited")
    require(len(worst_pairs) == 1, "minimum row is not unique")

    worst_left, worst_right, worst_residual, worst_top = worst_pairs[0]
    worst_intersection = np.flatnonzero(
        shallow[worst_left].astype(bool)
        & shallow[worst_right].astype(bool)
    )
    worst_capacities = (
        full_capacities_i32
        - loss_i32[worst_left]
        - loss_i32[worst_right]
        + root_matrix_transpose[worst_intersection]
        .sum(axis=0, dtype=np.uint16)
        .astype(np.int32)
    )
    canonical_order = np.lexsort(
        (np.arange(group_size), -worst_capacities)
    )
    worst_unit_exponents = tuple(
        map(int, canonical_order[:5])
    )
    canonical_top = tuple(
        map(int, worst_capacities[canonical_order[:5]])
    )
    require(canonical_top == worst_top, "canonical top five changed")

    # Direct full-residue hostile control for the unique minimum row.
    worst_shallow_labels = (
        sign_label(worst_left, quotient),
        sign_label(worst_right, quotient),
    )
    worst_unit_labels = tuple(
        sign_label(e, modulus) for e in worst_unit_exponents
    )
    full_universe = tuple(
        z
        for z in range(modulus)
        if z % P != 0 and 7 * norm_mod(z, modulus) > modulus
    )
    full_residual = tuple(
        z
        for z in full_universe
        if all(
            14 * norm_mod(P * a * z, modulus) >= modulus
            for a in worst_shallow_labels
        )
    )
    direct_capacities = tuple(
        sum(
            14 * norm_mod(q * z, modulus) < modulus
            for z in full_residual
        )
        for q in worst_unit_labels
    )
    require(
        len(full_universe) == 2 * universe_half,
        "sign quotient lost universe points",
    )
    require(
        len(full_residual) == 2 * worst_residual,
        "direct worst residual disagrees with the lift",
    )
    require(
        direct_capacities == tuple(2 * value for value in worst_top),
        "direct worst capacities disagree with the lift",
    )

    # A depth-(m-1) blocker is safe on every primitive depth-m numerator.
    deepest_sizes = tuple(
        sum(
            14
            * norm_mod(P ** (exponent - 1) * unit * int(z), modulus)
            < modulus
            for z in values
        )
        for unit in range(1, P)
    )
    require(
        deepest_sizes == (0,) * (P - 1),
        "the unique deepest blocker is not invisible",
    )

    return AuditResult(
        exponent=exponent,
        modulus=modulus,
        quotient=quotient,
        group_size=group_size,
        quotient_group_size=quotient_group_size,
        universe_half=universe_half,
        fibre_histogram=fibre_histogram,
        shallow_mask_size=shallow_size,
        pair_count=pair_count,
        minimum_half_margin=minimum_margin,
        worst_pair_exponents=(worst_left, worst_right),
        worst_pair_labels=worst_shallow_labels,
        worst_residual_half=worst_residual,
        worst_top_half=worst_top,
        worst_unit_exponents=worst_unit_exponents,
        worst_unit_labels=worst_unit_labels,
        root_matrix_sha256=root_matrix_digest,
        loss_matrix_sha256=loss_matrix_digest,
        pair_summary_sha256=pair_digest.hexdigest(),
        margin_histogram_size=len(margin_histogram),
    )


def main() -> None:
    control = audit_depth(3)
    require(control.minimum_half_margin == 43, "m=3 margin changed")
    require(control.worst_pair_labels == (14, 46), "m=3 pair changed")
    require(control.worst_residual_half == 523, "m=3 residual changed")
    require(
        control.worst_top_half == (102, 100, 95, 93, 90),
        "m=3 top five changed",
    )
    require(
        control.worst_unit_labels == (183, 799, 599, 1000, 1007),
        "m=3 unit labels changed",
    )

    result = audit_depth(4)
    require(result.pair_count == 514605, "m=4 pair count changed")
    require(result.minimum_half_margin == 643, "m=4 margin changed")
    require(
        result.worst_pair_exponents == (30, 599),
        "m=4 exponent pair changed",
    )
    require(
        result.worst_pair_labels == (183, 799),
        "m=4 shallow labels changed",
    )
    require(result.worst_residual_half == 6973, "m=4 residual changed")
    require(
        result.worst_top_half == (1350, 1326, 1253, 1217, 1184),
        "m=4 top five changed",
    )
    require(
        result.worst_unit_labels == (2380, 5193, 10386, 1190, 7139),
        "m=4 unit labels changed",
    )

    print("Scalar (1,1,3) exact cyclic-lift audit")
    print(
        "control_m3="
        f"pairs:{control.pair_count};"
        f"worst_pair:{control.worst_pair_labels};"
        f"residual:{2 * control.worst_residual_half};"
        f"top5:{tuple(2 * x for x in control.worst_top_half)};"
        f"margin:{2 * control.minimum_half_margin}"
    )
    print(
        f"m4=N:{result.modulus};Q:{result.quotient};"
        f"unit_sign_classes:{result.group_size};"
        f"shallow_sign_classes:{result.quotient_group_size}"
    )
    print(
        f"fibre_histogram_half={result.fibre_histogram};"
        f"universe_full={2 * result.universe_half};"
        f"shallow_mask_base_size={result.shallow_mask_size}"
    )
    print(
        f"shallow_pairs_with_repetition={result.pair_count};"
        f"margin_histogram_size={result.margin_histogram_size}"
    )
    print(
        f"worst_pair_exponents={result.worst_pair_exponents};"
        f"worst_pair_labels={result.worst_pair_labels}"
    )
    print(
        f"worst_residual_full={2 * result.worst_residual_half};"
        f"worst_top5_full={tuple(2 * x for x in result.worst_top_half)};"
        f"worst_unit_labels={result.worst_unit_labels};"
        f"minimum_margin_full={2 * result.minimum_half_margin}"
    )
    print(f"root_matrix_sha256={result.root_matrix_sha256}")
    print(f"loss_matrix_sha256={result.loss_matrix_sha256}")
    print(f"pair_summary_sha256={result.pair_summary_sha256}")
    print("arithmetic=float_free_uint8_uint16_int32")
    print("deepest_depth_three_mask_empty_on_primitive_layer=yes")
    print("PASS")


if __name__ == "__main__":
    main()
