#!/usr/bin/env python3
"""Exact hinge/Gram certificate for the scalar depth-(1,1,3) branch.

The script reconstructs the singleton residual-capacity table directly
from the primitive 13^4 torsion geometry.  It then proves one simultaneous
upper bound for all unordered pairs of depth-one blockers.

For x >= 0 and SCALE=S, ceil(sqrt(S^2 x))/S is an upper bound for sqrt(x).
Thus the integer Gram matrix below rigorously upper-bounds the Hellinger
tail kernel.  No floating-point arithmetic is used.
"""

import hashlib
import math

import numpy as np


P = 13
N = P**4
Q = P**3
MASK_COUNT = 5
THETA = 2612
SCALE = 100_000


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(array):
    return hashlib.sha256(array.tobytes()).hexdigest()


def sign_classes(modulus):
    representatives = np.arange(
        1,
        (modulus + 1) // 2,
        dtype=np.int64,
    )
    return representatives[representatives % P != 0]


def norm_numerators(values, modulus):
    residues = values % modulus
    return np.minimum(residues, modulus - residues)


def rows_to_bits(rows):
    packed = np.packbits(rows, axis=1, bitorder="little")
    return tuple(
        int.from_bytes(row.tobytes(), "little")
        for row in packed
    )


def ceiling_square_root(n):
    root = math.isqrt(n)
    return root if root * root == n else root + 1


def singleton_capacity_table(phases, roots, guard, shallow, units):
    """Return C_u(q) for every shallow blocker u and terminal label q."""
    full = (1 << len(phases)) - 1
    danger_rows = (
        14
        * norm_numerators(
            shallow[:, None] * phases[None, :],
            Q,
        )
        < Q
    )
    danger_bits = rows_to_bits(danger_rows)
    safe_bits = tuple(full ^ row for row in danger_bits)

    # On one thirteen-root fibre a unit terminal mask has zero, one, or two
    # guard-safe dangerous sheets.  Store its >=1 and >=2 superlevel sets
    # as exact Python integers.
    height_one = []
    height_two = []
    for start in range(0, len(units), 24):
        labels = units[start : start + 24]
        residues = labels[:, None, None] * roots[None, :, :] % N
        active = (
            14 * np.minimum(residues, N - residues) < N
        )
        heights = np.sum(
            active & guard[None, :, :],
            axis=2,
            dtype=np.uint8,
        )
        require(
            int(heights.min()) >= 0 and int(heights.max()) <= 2,
            "root-sheet height escaped {0,1,2}",
        )
        height_one.extend(rows_to_bits(heights >= 1))
        height_two.extend(rows_to_bits(heights >= 2))

    capacities = np.empty(
        (len(shallow), len(units)),
        dtype=np.uint16,
    )
    for index, residual in enumerate(safe_bits):
        capacities[index] = np.fromiter(
            (
                (residual & first).bit_count()
                + (residual & second).bit_count()
                for first, second in zip(height_one, height_two)
            ),
            dtype=np.uint16,
            count=len(units),
        )

    safe_rows = ~danger_rows
    return capacities, safe_rows


def run():
    phases = np.arange(1, Q, dtype=np.int64)
    phases = phases[phases % P != 0]
    shallow = sign_classes(Q)
    units = sign_classes(N)
    require(
        (len(phases), len(shallow), len(units))
        == (2028, 1014, 13182),
        "sign-class universe changed",
    )

    sheets = np.arange(P, dtype=np.int64)
    roots = phases[:, None] + Q * sheets[None, :]
    guard = 7 * norm_numerators(roots, N) > N
    guard_counts = guard.sum(axis=1, dtype=np.int64)
    require(
        int(guard_counts.min()) == 9
        and int(guard_counts.max()) == 10
        and int(guard_counts.sum()) == 18830,
        "primitive guard-safe universe changed",
    )

    capacities, safe_rows = singleton_capacity_table(
        phases,
        roots,
        guard,
        shallow,
        units,
    )
    capacity_digest = digest(capacities)
    require(
        capacity_digest
        == "f93144ef35b6b733ad08390a4d133d99c785fa92efde821fc909efd67bf0dc5a",
        "singleton-capacity table digest changed",
    )

    # W_uv is the exact number of guard-safe root sheets whose base phase
    # avoids both shallow blockers.
    safe_i64 = safe_rows.astype(np.int64)
    residual_weights = (
        safe_i64 * guard_counts[None, :]
    ) @ safe_i64.T

    # The Ky--Fan hinge tail of min(C_u,C_v) is bounded by
    # sum_q sqrt((C_u(q)-theta)_+ (C_v(q)-theta)_+).
    # Round every square-root feature upward before taking the Gram product.
    excess = np.maximum(
        capacities.astype(np.int64) - THETA,
        0,
    )
    maximum_excess = int(excess.max())
    root_upper = np.fromiter(
        (
            ceiling_square_root(n * SCALE * SCALE)
            for n in range(maximum_excess + 1)
        ),
        dtype=np.int64,
        count=maximum_excess + 1,
    )
    features = root_upper[excess]
    require(
        int(features.max()) ** 2 * len(units) < 2**63,
        "int64 Gram accumulator bound failed",
    )
    gram_upper_numerators = features @ features.T

    scale_squared = SCALE * SCALE
    margin_numerators = (
        (residual_weights - MASK_COUNT * THETA) * scale_squared
        - gram_upper_numerators
    )
    upper_triangle = np.triu_indices(len(shallow))
    upper_margins = margin_numerators[upper_triangle]
    minimum = int(upper_margins.min())
    minimizer_positions = np.flatnonzero(upper_margins == minimum)
    require(minimum > 0, "Hellinger-Gram certificate failed")

    records = []
    for position in minimizer_positions:
        left = int(upper_triangle[0][position])
        right = int(upper_triangle[1][position])
        records.append(
            (
                int(shallow[left]),
                int(shallow[right]),
                int(residual_weights[left, right]),
                int(gram_upper_numerators[left, right]),
                int(margin_numerators[left, right]),
            )
        )

    require(
        minimum == 786_095_718_332,
        f"minimum margin numerator changed: {minimum}",
    )
    require(
        records
        == [(5, 1098, 13580, 4_413_904_281_668, 786_095_718_332)],
        f"minimizer record changed: {records}",
    )

    print("depth-(1,1,3) fixed-threshold Hellinger-Gram certificate")
    print("arithmetic=exact integers; floating_point=none")
    print(
        f"phases={len(phases)} shallow={len(shallow)} "
        f"terminals={len(units)}"
    )
    print(
        f"pairs={len(upper_margins)} p={MASK_COUNT} "
        f"theta={THETA} scale={SCALE}"
    )
    print(f"capacity_table_digest={capacity_digest}")
    print(f"maximum_excess={maximum_excess}")
    print(
        "criterion_numerator="
        "(W-p*theta)*scale^2-(V_int V_int^T)"
    )
    print(f"minimum_margin_numerator={minimum}")
    print(f"minimum_margin={minimum}/{scale_squared}")
    print(f"minimizers={records}")
    print(f"margin_matrix_digest={digest(margin_numerators)}")


if __name__ == "__main__":
    run()
