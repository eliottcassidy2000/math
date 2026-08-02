#!/usr/bin/env python3
"""Exact referee for THM-3083's exceptional F2/F3 tomography clutch."""

from __future__ import annotations

import argparse
from collections import deque
from fractions import Fraction as F
from itertools import combinations, permutations
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "05-knowledge/results/exceptional_binary_point_ternary_direction_s4_tomography_clutch_thm3083.out"
)


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def compose_perm(left, right):
    """Return left after right."""
    return tuple(left[index] for index in right)


def matmul(left, right, p):
    return tuple(
        tuple(
            sum(left[i][k] * right[k][j] for k in range(2)) % p
            for j in range(2)
        )
        for i in range(2)
    )


def normalize_projective(matrix, p):
    flat = [entry % p for row in matrix for entry in row]
    pivot = next(entry for entry in flat if entry)
    inverse = pow(pivot, -1, p)
    flat = [(inverse * entry) % p for entry in flat]
    return ((flat[0], flat[1]), (flat[2], flat[3]))


def matvec(matrix, vector, p):
    return tuple(
        sum(matrix[i][j] * vector[j] for j in range(2)) % p
        for i in range(2)
    )


def inv2(matrix, p):
    a, b = matrix[0]
    c, d = matrix[1]
    det = (a * d - b * c) % p
    invdet = pow(det, -1, p)
    return (
        ((d * invdet) % p, (-b * invdet) % p),
        ((-c * invdet) % p, (a * invdet) % p),
    )


def determinant(matrix):
    size = len(matrix)
    if size == 0:
        return 1
    if size == 1:
        return matrix[0][0]
    total = 0
    for column in range(size):
        minor = tuple(
            tuple(row[index] for index in range(size) if index != column)
            for row in matrix[1:]
        )
        total += (-1) ** column * matrix[0][column] * determinant(minor)
    return total


def smith_from_minors(matrix):
    rows = len(matrix)
    columns = len(matrix[0])
    determinantal = [1]
    for size in range(1, min(rows, columns) + 1):
        divisor = 0
        for row_set in combinations(range(rows), size):
            for column_set in combinations(range(columns), size):
                minor = tuple(
                    tuple(matrix[row][column] for column in column_set)
                    for row in row_set
                )
                divisor = gcd(divisor, abs(determinant(minor)))
        determinantal.append(divisor)
    require(all(determinantal[index] for index in range(1, len(determinantal))), determinantal)
    smith = tuple(
        determinantal[index] // determinantal[index - 1]
        for index in range(1, len(determinantal))
    )
    return tuple(determinantal), smith


B = ((0, 0), (1, 0), (0, 1), (1, 1))
B_INDEX = {value: index for index, value in enumerate(B)}
D = ("I", "0", "1", "2")
D_VECTOR = {"I": (0, 1), "0": (1, 0), "1": (1, 1), "2": (1, 2)}
PHI = {"1": (0, 0), "I": (1, 0), "0": (0, 1), "2": (1, 1)}


def line_label(vector):
    x, y = vector
    if x % 3 == 0:
        return "I"
    return str((y * pow(x, -1, 3)) % 3)


def projective_perm(matrix):
    return tuple(
        D.index(line_label(matvec(matrix, D_VECTOR[label], 3))) for label in D
    )


def affine_perm(linear, shift):
    rows = []
    for point in B:
        linear_image = matvec(linear, point, 2)
        image = tuple(
            (value + delta) % 2 for value, delta in zip(linear_image, shift)
        )
        rows.append(B_INDEX[image])
    return tuple(rows)


I2 = ((1, 0), (0, 1))
S3 = ((0, 2), (1, 0))
T3 = ((1, 0), (1, 1))
R3 = matmul(S3, T3, 3)
J3 = ((1, 0), (0, 2))

SB = affine_perm(I2, (1, 1))
RB = affine_perm(((0, 1), (1, 1)), (0, 0))
JB = affine_perm(((0, 1), (1, 0)), (1, 1))
GENERATORS = ((SB, S3, "S"), (RB, R3, "R"), (JB, J3, "J"))


def generate_exceptional_group():
    identity = (tuple(range(4)), I2)
    queue = deque([identity])
    seen = {identity}
    while queue:
        bperm, matrix = queue.popleft()
        for generator_b, generator_m, _name in GENERATORS:
            candidate = (
                compose_perm(generator_b, bperm),
                normalize_projective(matmul(generator_m, matrix, 3), 3),
            )
            if candidate not in seen:
                seen.add(candidate)
                queue.append(candidate)
    return tuple(sorted(seen, key=repr))


GROUP = generate_exceptional_group()
require(len(GROUP) == 24, ("group order", len(GROUP)))
require(
    len({projective_perm(matrix) for _bperm, matrix in GROUP}) == 24,
    "PGL action not faithful",
)
for bperm, matrix in GROUP:
    dperm = projective_perm(matrix)
    for index, label in enumerate(D):
        require(
            B[bperm[B_INDEX[PHI[label]]]] == PHI[D[dperm[index]]],
            (label, bperm, matrix),
        )

equivariant_bijections = 0
for image_indices in permutations(range(4)):
    candidate = {D[index]: B[image_indices[index]] for index in range(4)}
    if all(
        B[bperm[B_INDEX[candidate[label]]]] == candidate[D[projective_perm(matrix)[index]]]
        for bperm, matrix in GROUP
        for index, label in enumerate(D)
    ):
        equivariant_bijections += 1
require(equivariant_bijections == 1, equivariant_bijections)


V3 = tuple((x, y) for x in range(3) for y in range(3))
V3_INDEX = {value: index for index, value in enumerate(V3)}
LINES = {
    label: frozenset(
        tuple((scalar * coordinate) % 3 for coordinate in D_VECTOR[label])
        for scalar in range(3)
    )
    for label in D
}
require(all(len(line) == 3 for line in LINES.values()), LINES)


def dot(first, second):
    return sum(a * b for a, b in zip(first, second))


ZERO = tuple(F(0) for _point in V3)
ONE = tuple(F(1) for _point in V3)
DELTA_ZERO = tuple(F(point == (0, 0)) for point in V3)
H = {
    label: tuple(F(2 if point in LINES[label] else -1) for point in V3)
    for label in D
}
require(
    [[dot(H[left], H[right]) for right in D] for left in D]
    == [[F(18 if i == j else 0) for j in range(4)] for i in range(4)],
    "H Gram",
)
require(all(dot(ONE, H[label]) == 0 for label in D), "H mean")
require(
    all(
        tuple(value * value for value in H[label]) != H[label]
        for label in D
    ),
    "channel idempotence hostile",
)

# The projector-normalized channels are not idempotents, but puncturing the
# common origin produces an exact multiplicative S4 clutch.
PUNCTURED = {
    label: tuple(
        F(point in LINES[label] and point != (0, 0)) for point in V3
    )
    for label in D
}
for left in D:
    for right in D:
        product = tuple(
            a * b for a, b in zip(PUNCTURED[left], PUNCTURED[right])
        )
        require(
            product == (PUNCTURED[left] if left == right else ZERO),
            (left, right, "punctured idempotents"),
        )
require(
    tuple(sum(PUNCTURED[label][index] for label in D) for index in range(9))
    == tuple(one - delta for one, delta in zip(ONE, DELTA_ZERO)),
    "punctured ideal unit",
)


def add3(left, right):
    return tuple((a + b) % 3 for a, b in zip(left, right))


def p_line(function, label):
    return tuple(
        sum(
            function[V3_INDEX[add3(point, shift)]] for shift in LINES[label]
        )
        / 3
        for point in V3
    )


def p_all(function):
    mean = sum(function) / 9
    return tuple(mean for _point in V3)


def e_line(function, label):
    return tuple(
        value - mean for value, mean in zip(p_line(function, label), p_all(function))
    )


for left in D:
    for right in D:
        require(
            e_line(H[right], left) == (H[right] if left == right else ZERO),
            (left, right, "projector"),
        )
for basis in (ONE, *(H[label] for label in D)):
    require(
        tuple(
            sum(e_line(basis, label)[index] for label in D)
            for index in range(9)
        )
        == tuple(value - mean for value, mean in zip(basis, p_all(basis))),
        (basis, "even reconstruction"),
    )


def act_v3(function, matrix):
    inverse = inv2(matrix, 3)
    return tuple(
        function[V3_INDEX[matvec(inverse, point, 3)]] for point in V3
    )


for _bperm, matrix in GROUP:
    dperm = projective_perm(matrix)
    for index, label in enumerate(D):
        require(
            act_v3(H[label], matrix) == H[D[dperm[index]]],
            (matrix, label, "H action"),
        )


def phi_map(binary_function):
    require(sum(binary_function) == 0, ("not centered", binary_function))
    return tuple(
        sum(
            F(binary_function[B_INDEX[PHI[label]]]) * H[label][point_index]
            for label in D
        )
        for point_index in range(9)
    )


def psi_map(binary_function):
    return tuple(
        sum(
            F(binary_function[B_INDEX[PHI[label]]])
            * PUNCTURED[label][point_index]
            for label in D
        )
        for point_index in range(9)
    )


def act_b(binary_function, permutation):
    inverse = [0] * 4
    for source, target in enumerate(permutation):
        inverse[target] = source
    return tuple(binary_function[inverse[index]] for index in range(4))


centered_basis = tuple(
    tuple(1 if index == basis_index else -1 if index == 3 else 0 for index in range(4))
    for basis_index in range(3)
)
for binary in centered_basis:
    image = phi_map(binary)
    require(
        image == tuple(3 * value for value in psi_map(binary)),
        (binary, "Phi=3Psi on augmentation"),
    )
    require(dot(image, image) == 18 * dot(binary, binary), (binary, "norm"))
    for bperm, matrix in GROUP:
        require(
            phi_map(act_b(binary, bperm)) == act_v3(image, matrix),
            (binary, bperm, matrix, "equivariance"),
        )

# Psi is an equivariant algebra isomorphism from the four-point function
# algebra onto the punctured ideal.  Its ambient integral matrix is a
# permutation matrix, so its A3 augmentation image is saturated.
binary_deltas = tuple(
    tuple(1 if index == source else 0 for index in range(4))
    for source in range(4)
)
for binary in binary_deltas:
    image = psi_map(binary)
    for bperm, matrix in GROUP:
        require(
            psi_map(act_b(binary, bperm)) == act_v3(image, matrix),
            (binary, bperm, matrix, "Psi equivariance"),
        )
psi_integral_matrix = tuple(
    tuple(
        int(psi_map(binary_deltas[column])[V3_INDEX[D_VECTOR[label]]])
        for column in range(4)
    )
    for label in D
)
require(abs(determinant(psi_integral_matrix)) == 1, psi_integral_matrix)

# A unital equivariant map to the full five-orbit even algebra would have to
# send the fixed ternary origin to an S4-fixed binary point.  None exists.
fixed_binary_points = tuple(
    point_index
    for point_index in range(4)
    if all(bperm[point_index] == point_index for bperm, _matrix in GROUP)
)
require(fixed_binary_points == (), fixed_binary_points)


MATCHING_ROWS = []
for character in ((1, 0), (0, 1), (1, 1)):
    walsh = tuple(
        1
        if (character[0] * point[0] + character[1] * point[1]) % 2 == 0
        else -1
        for point in B
    )
    signs = tuple(walsh[B_INDEX[PHI[label]]] for label in D)
    require(signs.count(1) == 2 and signs.count(-1) == 2, (character, signs))
    require(
        phi_map(walsh)
        == tuple(
            sum(F(signs[index]) * H[D[index]][point] for index in range(4))
            for point in range(9)
        ),
        (character, "matching"),
    )
    MATCHING_ROWS.append(
        (
            character,
            tuple(D[index] for index, sign in enumerate(signs) if sign == 1),
            signs,
        )
    )


RADIAL = tuple(sum(H[label][point] for label in D) for point in range(9))
require(
    tuple(int(value) for value in RADIAL)
    == tuple(8 if point == (0, 0) else -1 for point in V3),
    "radial",
)
require(all(dot(phi_map(binary), RADIAL) == 0 for binary in centered_basis), "radial complement")


channel_to_even = tuple(
    tuple(3 * (row == column) - 1 for column in range(4))
    for row in range(4)
)
channel_divisors, channel_smith = smith_from_minors(channel_to_even)
require(channel_divisors == (1, 1, 3, 9, 27), channel_divisors)
require(channel_smith == (1, 3, 3, 3), channel_smith)

composite = tuple(
    tuple(
        (-1 if column == 0 else 3 * ((row == column - 1) - (row == 3)))
        for column in range(4)
    )
    for row in range(4)
)
composite_divisors, composite_smith = smith_from_minors(composite)
require(composite_divisors == (1, 1, 3, 9, 108), composite_divisors)
require(composite_smith == (1, 3, 3, 12), composite_smith)


def matching_squares(roots):
    return tuple(
        sorted(
            {
                (roots[0] + roots[1]) ** 2,
                (roots[0] + roots[2]) ** 2,
                (roots[0] + roots[3]) ** 2,
            }
        )
    )


roots_plus = (0, 1, 2, -3)
roots_minus = tuple(-value for value in roots_plus)
require(matching_squares(roots_plus) == (1, 4, 9), matching_squares(roots_plus))
require(matching_squares(roots_minus) == (1, 4, 9), matching_squares(roots_minus))
require(
    sum(roots_plus[i] * roots_plus[j] * roots_plus[k] for i, j, k in combinations(range(4), 3))
    == -6,
    "plus cubic elementary sum",
)
require(
    sum(roots_minus[i] * roots_minus[j] * roots_minus[k] for i, j, k in combinations(range(4), 3))
    == 6,
    "minus cubic elementary sum",
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    lines = [
        "Exceptional binary-point / ternary-direction S4 tomography clutch",
        f"group_order={len(GROUP)};equivariant_bijections_after_generator_fix={equivariant_bijections}",
        "generators=S:translation_(1,1)<->z|-1/z;R:linear_order3<->z|-1/(z+1);J:affine_reflection<->z|-z",
        "ternary_even_split=constant_rank1+four_centered_direction_channels_rank1",
        "channel_vectors=h_L=3*1_L-1;gram=18*I4;sum=radial(8_at_0,-1_elsewhere)",
        "binary_centered_to_ternary_standard=Phi(f)=sum_d f(phi(d))*h_d;gram_scale=18;S4_equivariant",
        f"channel_lattice_smith={channel_smith};index={channel_divisors[-1]};quotient=F3^3",
        f"radial_plus_standard_smith={composite_smith};index={composite_divisors[-1]};quotient=Z4+F3^3",
        "punctured_line_algebra=Psi(delta_phi(L))=1_(L\\{0});orthogonal_idempotents;sum=1-delta0;S4_equivariant",
        "augmentation_relation=Phi=3*Psi;saturated_integral_Psi_A3;projector_normalization_has_F3^3_defect",
    ]
    for character, positive, signs in MATCHING_ROWS:
        lines.append(
            f"binary_character={character};positive_directions={positive};signs={signs}"
        )
    lines += [
        "quartic_phase_hostile=f_plus:T4-7T2+6T;f_minus:T4-7T2-6T;matching_squares:(1,4,9);linear_signs:opposite",
        "normalized_multiplicative_hostile=delta_b^2=delta_b but the displayed extension delta_phi(L)|h_L has h_L^2!=h_L",
        "unital_boundary=no S4-equivariant unital algebra map k[B]->W_even because ternary origin is fixed and B has no S4-fixed point",
        "first_failed_implication=no individual-Walsh-to-single-line or unital full-even algebra clutch;punctured-line ideal clutch survives",
        "scope=representation-and-integral-lattice clutch only;no physical quartic,Cardano,Weil,Keller,owner,current,or LRC consequence",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
