#!/usr/bin/env python3
"""Exact standard-library certificate for THM-4298.

The universe is k[p,y] with wt(p)=2, wt(y)=3 under the source-normal seam
substitution p=t(1+x^2*t), y=x*t*p.  This script audits the integral
binomial transform from every weight face to the diagonal 2*n-d=M, its closed
inverse, exact row interval, face isolation, and the M=12 wall dictionary.
"""

from __future__ import annotations

from fractions import Fraction as F
from math import comb


Monomial = tuple[int, int]
XYTerm = tuple[int, int]
Polynomial = dict[Monomial, F]


def require(condition: bool, label: str) -> None:
    if not condition:
        raise AssertionError(label)


def weight(i: int, j: int) -> int:
    return 2 * i + 3 * j


def face(M: int) -> list[Monomial]:
    result: list[Monomial] = []
    for j in range(M // 3 + 1):
        remainder = M - 3 * j
        if remainder >= 0 and remainder % 2 == 0:
            result.append((remainder // 2, j))
    return result


def onset(i: int, j: int) -> XYTerm:
    return j, i + 2 * j


def transform(M: int) -> list[list[int]]:
    monomials = face(M)
    result: list[list[int]] = []
    for s, _ in enumerate(monomials):
        row: list[int] = []
        for r, (i_r, j_r) in enumerate(monomials):
            q = s - r
            row.append(comb(i_r + j_r, q) if 0 <= q <= i_r + j_r else 0)
        result.append(row)
    return result


def closed_inverse(M: int) -> list[list[int]]:
    monomials = face(M)
    result: list[list[int]] = []
    for s, _ in enumerate(monomials):
        row: list[int] = []
        for r, (i_r, j_r) in enumerate(monomials):
            q = s - r
            if 0 <= q <= i_r + j_r:
                row.append((-1) ** q * comb(i_r + j_r, q))
            else:
                row.append(0)
        result.append(row)
    return result


def invert_integer_matrix(matrix: list[list[int]]) -> list[list[int]]:
    n = len(matrix)
    augmented = [
        [F(value) for value in row] + [F(int(i == j)) for j in range(n)]
        for i, row in enumerate(matrix)
    ]
    for column in range(n):
        pivot = next((row for row in range(column, n) if augmented[row][column]), None)
        require(pivot is not None, f"pivot column {column}")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [value / scale for value in augmented[column]]
        for row in range(n):
            if row == column:
                continue
            scale = augmented[row][column]
            if scale:
                augmented[row] = [
                    left - scale * right
                    for left, right in zip(augmented[row], augmented[column])
                ]
    inverse_q = [row[n:] for row in augmented]
    require(
        all(value.denominator == 1 for row in inverse_q for value in row),
        "integral matrix inverse",
    )
    return [[int(value) for value in row] for row in inverse_q]


def matmul(left: list[list[int]], right: list[list[int]]) -> list[list[int]]:
    if not left:
        return []
    return [
        [
            sum(left[i][k] * right[k][j] for k in range(len(right)))
            for j in range(len(right[0]))
        ]
        for i in range(len(left))
    ]


def identity(n: int) -> list[list[int]]:
    return [[int(i == j) for j in range(n)] for i in range(n)]


def seam_expand(poly: Polynomial) -> dict[XYTerm, F]:
    result: dict[XYTerm, F] = {}
    for (i, j), coefficient in poly.items():
        for q in range(i + j + 1):
            term = (j + 2 * q, i + 2 * j + q)
            result[term] = result.get(term, F(0)) + coefficient * comb(i + j, q)
    return {term: coefficient for term, coefficient in result.items() if coefficient}


def diagonal_channels(expanded: dict[XYTerm, F], M: int) -> list[F]:
    return [expanded.get(onset(i, j), F(0)) for i, j in face(M)]


def apply_matrix(matrix: list[list[int]], vector: list[F]) -> list[F]:
    return [sum(F(value) * entry for value, entry in zip(row, vector)) for row in matrix]


def audit_termwise() -> int:
    checked = 0
    for i in range(41):
        for j in range(41):
            M = weight(i, j)
            for q in range(i + j + 1):
                d = j + 2 * q
                n = i + 2 * j + q
                require(2 * n - d == M, "termwise diagonal invariant")
                checked += 1
    return checked


def audit_transforms() -> tuple[int, int]:
    weights = 0
    matrix_entries = 0
    for M in range(201):
        monomials = face(M)
        if not monomials:
            require(M == 1, "unique Frobenius gap")
            continue
        rows = [onset(i, j)[1] for i, j in monomials]
        require(rows == list(range(rows[0], rows[0] + len(rows))), "consecutive rows")
        require(rows[0] == (M + 1) // 2, "first row")
        require(rows[-1] == (2 * M) // 3, "last row")
        forward = transform(M)
        require(all(forward[r][r] == 1 for r in range(len(forward))), "unit diagonal")
        require(
            all(forward[r][s] == 0 for r in range(len(forward)) for s in range(r + 1, len(forward))),
            "lower triangular",
        )
        inverse = invert_integer_matrix(forward)
        require(inverse == closed_inverse(M), "closed inverse")
        require(matmul(forward, inverse) == identity(len(forward)), "right inverse")
        require(matmul(inverse, forward) == identity(len(forward)), "left inverse")
        weights += 1
        matrix_entries += len(forward) ** 2
    return weights, matrix_entries


def audit_mixed_isolation() -> int:
    poly: Polynomial = {}
    for i in range(51):
        for j in range(35):
            if weight(i, j) <= 100:
                sign = -1 if (i + j) % 2 else 1
                poly[(i, j)] = F(sign * (1 + 3 * i + 5 * j))
    expanded = seam_expand(poly)
    checked = 0
    for M in range(101):
        monomials = face(M)
        if not monomials:
            continue
        coefficients = [poly.get(monomial, F(0)) for monomial in monomials]
        require(
            diagonal_channels(expanded, M) == apply_matrix(transform(M), coefficients),
            f"mixed face isolation M={M}",
        )
        checked += len(monomials)
    return checked


def audit_m12() -> tuple[list[list[int]], list[list[int]], int]:
    require(face(12) == [(6, 0), (3, 2), (0, 4)], "M12 face")
    forward = transform(12)
    inverse = closed_inverse(12)
    require(forward == [[1, 0, 0], [6, 1, 0], [15, 5, 1]], "M12 transform")
    require(inverse == [[1, 0, 0], [-6, 1, 0], [15, -5, 1]], "M12 inverse")
    checked = 0
    for a in range(-8, 9):
        for b in range(-8, 9):
            for c in range(-8, 9):
                U = a
                W = b - 6 * a
                Z = c - 5 * b + 15 * a
                require(U + W + Z == c - 4 * b + 10 * a, "Lambda dictionary")
                require(
                    W * W - 4 * U * Z == b * b + 8 * a * b - 24 * a * a - 4 * a * c,
                    "D dictionary",
                )
                checked += 1
    p6 = {(6, 0): F(1)}
    p3y2 = {(3, 2): F(1)}
    require(diagonal_channels(seam_expand(p6), 12) == [1, 6, 15], "p6 hostile")
    require(diagonal_channels(seam_expand(p3y2), 12) == [0, 1, 5], "p3y2 hostile")
    return forward, inverse, checked


def main() -> None:
    term_count = audit_termwise()
    weight_count, matrix_entries = audit_transforms()
    channel_count = audit_mixed_isolation()
    forward, inverse, wall_count = audit_m12()

    print("THM-4298 WEIGHTED-FACE VISIBILITY PRIMARY AUDIT: PASS")
    print(f"termwise_invariant_checks={term_count}")
    print(f"unimodular_weights={weight_count}; matrix_entries={matrix_entries}; range=0..200_except_1")
    print(f"mixed_face_channels={channel_count}; mixed_weight_range=0..100")
    print(f"M12_transform={forward}")
    print(f"M12_inverse={inverse}")
    print("M12_channels: h0=U; h1=6U+W; h2=15U+5W+Z")
    print("M12_walls: U=h0; Z=h2-5h1+15h0; Lambda=h2-4h1+10h0")
    print("M12_discriminant: D=h1^2+8h0h1-24h0^2-4h0h2")
    print(f"M12_wall_triples={wall_count}")
    print("hostile: p^6_packet=(1,6,15); p^3y^2_packet=(0,1,5)")
    print("scope: lossless face observer; no Darboux lift, seam entry, or JC2 conclusion")


if __name__ == "__main__":
    main()
