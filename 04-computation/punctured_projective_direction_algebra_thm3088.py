#!/usr/bin/env python3
"""Exact controls for THM-3088's punctured direction algebra."""

from __future__ import annotations

import argparse
from itertools import product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "05-knowledge/results/punctured_projective_direction_algebra_thm3088.out"
)


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


class Field:
    """Prime fields and the two quadratic controls F4/F9."""

    def __init__(self, p, modulus=None):
        self.p = p
        self.modulus = modulus
        if modulus is None:
            self.elements = tuple(range(p))
            self.zero = 0
            self.one = 1
        else:
            self.elements = tuple(product(range(p), repeat=2))
            self.zero = (0, 0)
            self.one = (1, 0)

    @property
    def q(self):
        return len(self.elements)

    def add(self, left, right):
        if self.modulus is None:
            return (left + right) % self.p
        return tuple((a + b) % self.p for a, b in zip(left, right))

    def neg(self, value):
        if self.modulus is None:
            return (-value) % self.p
        return tuple((-entry) % self.p for entry in value)

    def sub(self, left, right):
        return self.add(left, self.neg(right))

    def mul(self, left, right):
        if self.modulus is None:
            return (left * right) % self.p
        # modulus X^2+a1 X+a0, so X^2=-a1 X-a0
        a0, a1 = self.modulus
        c0 = left[0] * right[0]
        c1 = left[0] * right[1] + left[1] * right[0]
        c2 = left[1] * right[1]
        return (
            (c0 - a0 * c2) % self.p,
            (c1 - a1 * c2) % self.p,
        )

    def inv(self, value):
        require(value != self.zero, (self.q, "inverse zero"))
        for candidate in self.elements:
            if self.mul(value, candidate) == self.one:
                return candidate
        raise RuntimeError((self.q, value, "no inverse"))


FIELDS = (
    Field(2),
    Field(3),
    Field(2, (1, 1)),  # X^2+X+1
    Field(5),
    Field(7),
    Field(3, (1, 0)),  # X^2+1
)


def vector_neg(field, vector):
    return tuple(field.neg(value) for value in vector)


def vector_scale(field, scalar, vector):
    return tuple(field.mul(scalar, value) for value in vector)


def field_lines(field):
    lines = []
    labels = []
    for slope in field.elements:
        lines.append(
            frozenset(
                (scalar, field.mul(slope, scalar))
                for scalar in field.elements
            )
        )
        labels.append(("slope", slope))
    lines.append(frozenset((field.zero, scalar) for scalar in field.elements))
    labels.append(("infinity",))
    return tuple(labels), tuple(lines)


def parity_orbits(field, vectors):
    remaining = set(vectors)
    rows = []
    while remaining:
        value = min(remaining, key=repr)
        orbit = frozenset((value, vector_neg(field, value)))
        rows.append(orbit)
        remaining -= orbit
    return tuple(sorted(rows, key=repr))


def determinant_2(field, matrix):
    return field.sub(
        field.mul(matrix[0][0], matrix[1][1]),
        field.mul(matrix[0][1], matrix[1][0]),
    )


def matrix_vector(field, matrix, vector):
    return tuple(
        field.add(
            field.mul(matrix[row][0], vector[0]),
            field.mul(matrix[row][1], vector[1]),
        )
        for row in range(2)
    )


def audit_field(field):
    q = field.q
    zero_vector = (field.zero, field.zero)
    vectors = tuple(product(field.elements, repeat=2))
    labels, lines = field_lines(field)
    require(len(lines) == q + 1 and len(set(lines)) == q + 1, (q, "lines"))
    require(all(len(line) == q and zero_vector in line for line in lines), q)

    punctured = tuple(frozenset(line - {zero_vector}) for line in lines)
    incidences = {}
    for index, line in enumerate(punctured):
        for point in line:
            require(point not in incidences, (q, point, "two directions"))
            incidences[point] = index
    require(set(incidences) == set(vectors) - {zero_vector}, (q, "partition"))

    # Each punctured line is one F_q^*-orbit.  Thus the direction algebra is
    # exactly the scalar-invariant algebra on V\{0}, not merely a subalgebra
    # with the right dimension.
    nonzero_scalars = tuple(value for value in field.elements if value != field.zero)
    for line in punctured:
        point = next(iter(line))
        scalar_orbit = frozenset(
            vector_scale(field, scalar, point) for scalar in nonzero_scalars
        )
        require(scalar_orbit == line, (q, point, "scalar orbit"))

    # Orthogonal Boolean idempotents and punctured unit.
    for i, left in enumerate(punctured):
        for j, right in enumerate(punctured):
            require((left & right) == (left if i == j else frozenset()), (q, i, j))

    orbits = parity_orbits(field, vectors)
    zero_orbit = next(orbit for orbit in orbits if zero_vector in orbit)
    nonzero_orbits = tuple(orbit for orbit in orbits if orbit != zero_orbit)
    orbit_blocks = []
    for line in punctured:
        block = tuple(index for index, orbit in enumerate(nonzero_orbits) if orbit <= line)
        orbit_blocks.append(block)
    require(all(orbit_blocks), (q, "empty orbit block"))
    require(
        sorted(index for block in orbit_blocks for index in block)
        == list(range(len(nonzero_orbits))),
        (q, "orbit partition"),
    )
    block_sizes = {len(block) for block in orbit_blocks}
    expected_block = q - 1 if field.p == 2 else (q - 1) // 2
    require(block_sizes == {expected_block}, (q, block_sizes, expected_block))

    direction_rank = q + 1
    parity_ideal_rank = len(nonzero_orbits)
    internal_contrast_rank = parity_ideal_rank - direction_rank
    expected_internal = (
        (q + 1) * (q - 2)
        if field.p == 2
        else (q + 1) * (q - 3) // 2
    )
    require(internal_contrast_rank == expected_internal, (q, internal_contrast_rank))
    require(
        (direction_rank == parity_ideal_rank) == (q in (2, 3)),
        (q, "exceptional saturation"),
    )

    # Choosing one parity orbit in each block gives an identity maximal minor,
    # so the integral direction lattice is saturated.
    pivots = tuple(block[0] for block in orbit_blocks)
    pivot_matrix = tuple(
        tuple(int(pivots[row] in orbit_blocks[column]) for column in range(q + 1))
        for row in range(q + 1)
    )
    require(
        pivot_matrix
        == tuple(
            tuple(int(row == column) for column in range(q + 1))
            for row in range(q + 1)
        ),
        (q, "integral saturation"),
    )

    # h_L=q 1_L-1 differs by q times the punctured idempotent on the
    # augmentation lattice.
    point_index = {point: index for index, point in enumerate(vectors)}
    h_vectors = tuple(
        tuple(q * int(point in line) - 1 for point in vectors) for line in lines
    )
    e_vectors = tuple(
        tuple(int(point in line) for point in vectors) for line in punctured
    )
    for index in range(q):
        require(
            tuple(
                h_vectors[index][point] - h_vectors[-1][point]
                for point in range(q * q)
            )
            == tuple(
                q * (e_vectors[index][point] - e_vectors[-1][point])
                for point in range(q * q)
            ),
            (q, index, "projector scaling"),
        )
    require(len(point_index) == q * q, q)

    # Exhaust the GL2 action for the three smallest fields.
    gl_count = 0
    if q <= 4:
        line_lookup = {line: index for index, line in enumerate(lines)}
        for entries in product(field.elements, repeat=4):
            matrix = ((entries[0], entries[1]), (entries[2], entries[3]))
            if determinant_2(field, matrix) == field.zero:
                continue
            gl_count += 1
            images = []
            for line in lines:
                image = frozenset(matrix_vector(field, matrix, point) for point in line)
                require(image in line_lookup, (q, matrix, image))
                images.append(line_lookup[image])
            require(len(set(images)) == q + 1, (q, matrix, "line action"))
    expected_gl = (q * q - 1) * (q * q - q) if q <= 4 else 0
    require(gl_count == expected_gl, (q, gl_count, expected_gl))

    return (
        q,
        field.p,
        q + 1,
        expected_block,
        parity_ideal_rank,
        internal_contrast_rank,
        q in (2, 3),
        gl_count,
    )


ROWS = tuple(audit_field(field) for field in FIELDS)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    lines = [
        "Punctured projective-direction algebra and exceptional parity saturation",
        "identity=e_L=1_(L\\{0});e_L*e_M=delta_LM*e_L;sum_L e_L=1-delta0",
        "algebra=Fun(P1(Fq))=Fun(Fq^2\\{0})^(Fq*) is GL2-equivariant and integrally saturated",
    ]
    for q, p, directions, block, parity_rank, contrast, saturated, gl_count in ROWS:
        lines.append(
            f"q={q};char={p};directions={directions};parity_orbits_per_direction={block};"
            f"parity_ideal_rank={parity_rank};internal_contrast_rank={contrast};"
            f"full_parity_saturation={int(saturated)};GL2_exhausted={gl_count}"
        )
    lines += [
        "exceptional_fields=q2(singleton punctured lines),q3(antipodal-pair punctured lines)",
        "first_failures=q4 leaves residual C3 within each line;q5 leaves residual F5*/{+-1}=C2",
        "projector_relation=on augmentation,h_L=q*1_L-1 gives Phi_q=q*Psi_q;Smith=q repeated q",
        "larger_field_boundary=internal line contrasts survive with rank (q+1)(q-2) in char2 and (q+1)(q-3)/2 in odd characteristic",
        "scope=pointwise finite-field algebra only;no convolution,Farey-tree,quartic-owner,Keller,or LRC intertwiner",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
