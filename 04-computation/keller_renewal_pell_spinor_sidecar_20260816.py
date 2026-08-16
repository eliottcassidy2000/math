#!/usr/bin/env python3
"""Exact Pell/spinor sidecar for the THM-3522 renewal matrix.

The abstract packet recurrence is

    (e',m') = (7e-2m, 3e-2m).

With X=6e-9m and Y=m it becomes multiplication by the quadratic integer
alpha=(5+sqrt(57))/2.  The script checks the resulting norm conic, Cassini
determinant, harmonic K4 face, primitive Pythagorean current, and induced
Lorentz similitude.  Only indices 0..6 are currently realized by proved
polynomial packets in the fixed Keller tower; later rows are arithmetic
predictions conditional on future polynomiality gates.
"""

from __future__ import annotations

from hashlib import sha256
import json
from math import gcd


EXPECTED_SEMANTIC_SHA256 = (
    "c5758a441d60c45f254edbfea3d2ec06c34886fdcd0fc9bc27685dc0aaa4a5af"
)

M = ((7, -2), (3, -2))
PELL_TWICE = ((5, 57), (1, 5))
CHANGE = ((6, -9), (0, 1))
TRIPLE_MAP = (
    (20, -8, 20),
    (17, -20, 25),
    (25, -20, 33),
)
LORENTZ = ((-1, 0, 0), (0, -1, 0), (0, 0, 1))
FIXED_NAMES = ("L", "H", "J", "G", "R5", "R6", "R7")
EXPECTED_FIXED = (
    (1, 0),
    (7, 3),
    (43, 15),
    (271, 99),
    (1699, 615),
    (10663, 3867),
    (66907, 24255),
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def matmul(left, right):
    return tuple(
        tuple(
            sum(left[row][inner] * right[inner][column]
                for inner in range(len(right)))
            for column in range(len(right[0]))
        )
        for row in range(len(left))
    )


def transpose(matrix):
    return tuple(zip(*matrix))


def scale(matrix, scalar: int):
    return tuple(tuple(scalar * value for value in row) for row in matrix)


def matvec(matrix, vector):
    return tuple(sum(row[index] * vector[index] for index in range(len(vector)))
                 for row in matrix)


def determinant_2(matrix) -> int:
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def determinant_3(matrix) -> int:
    a, b, c = matrix[0]
    d, e, f = matrix[1]
    g, h, i = matrix[2]
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def next_packet(pair):
    return matvec(M, pair)


def pell_coordinates(pair):
    e, m = pair
    return 6 * e - 9 * m, m


def euclid_triple(pair):
    e, m = pair
    return (e * e - m * m) // 2, e * m, (e * e + m * m) // 2


def harmonic_face(pair):
    e, m = pair
    x = e - m
    return m * x, e * x, e * m


def iterate(count: int):
    pairs = [(1, 0)]
    for _ in range(count - 1):
        pairs.append(next_packet(pairs[-1]))
    return tuple(pairs)


def main() -> None:
    # The half-integral Pell action is recorded without floating arithmetic:
    # 2*C*M = [[5,57],[1,5]]*C.
    require(scale(matmul(CHANGE, M), 2) == matmul(PELL_TWICE, CHANGE),
            "Pell semiconjugacy")
    require(determinant_2(M) == -8, "renewal determinant")
    require(determinant_3(TRIPLE_MAP) == -512, "spinor determinant")
    require(
        matmul(transpose(TRIPLE_MAP), matmul(LORENTZ, TRIPLE_MAP))
        == scale(LORENTZ, 64),
        "Lorentz multiplier",
    )

    pairs = iterate(65)
    require(pairs[:len(EXPECTED_FIXED)] == EXPECTED_FIXED,
            ("fixed packet prefix", pairs[:len(EXPECTED_FIXED)]))
    rows = []
    for index, pair in enumerate(pairs):
        e, m = pair
        x_pell, y_pell = pell_coordinates(pair)
        norm = x_pell * x_pell - 57 * y_pell * y_pell
        quadratic = 3 * e * e - 9 * e * m + 2 * m * m
        require(norm == 36 * (-8) ** index,
                ("Pell norm", index, pair, norm))
        require(quadratic == 3 * (-8) ** index,
                ("packet conic", index, pair, quadratic))
        if index >= 2:
            require(e == 5 * pairs[index - 1][0] + 8 * pairs[index - 2][0],
                    ("e recurrence", index))
            require(m == 5 * pairs[index - 1][1] + 8 * pairs[index - 2][1],
                    ("m recurrence", index))
        if index + 1 < len(pairs):
            e_next, m_next = pairs[index + 1]
            require(e * m_next - m * e_next == 3 * (-8) ** index,
                    ("Cassini determinant", index))

        if index == 0:
            rows.append((index, pair, (x_pell, y_pell), norm, None, None))
            continue

        require((e % 6, m % 6) == (1, 3),
                ("mod-six packet", index, pair))
        require(gcd(e, m) == 1 and e > m > 0,
                ("primitive packet", index, pair))
        triple = euclid_triple(pair)
        a, b, c = triple
        require(a * a + b * b == c * c and gcd(gcd(a, b), c) == 1,
                ("primitive Pythagorean current", index, triple))
        require(a - 9 * b + 5 * c == 3 * (-8) ** index,
                ("moving Pythagorean plane", index, triple))
        if index + 1 < len(pairs):
            require(matvec(TRIPLE_MAP, triple) == euclid_triple(pairs[index + 1]),
                    ("spinor map", index))

        u, v, z = harmonic_face(pair)
        require(v * z == u * (v + z),
                ("harmonic K4 face", index, (u, v, z)))
        require((gcd(u, v), gcd(u, z), gcd(v, z)) == (e - m, m, e),
                ("harmonic gcd decoder", index, pair))
        rows.append((index, pair, (x_pell, y_pell), norm, triple, (u, v, z)))

    # Hostiles: neither dropping the -2m term nor treating the triple map as
    # a Berggren isometry preserves the exact structures.
    wrong_m = ((7, -1), (3, -2))
    wrong_pair = matvec(wrong_m, pairs[1])
    wrong_x, wrong_y = pell_coordinates(wrong_pair)
    hostile = (
        wrong_pair,
        wrong_x * wrong_x - 57 * wrong_y * wrong_y,
        36 * (-8) ** 2,
        matmul(transpose(TRIPLE_MAP), matmul(LORENTZ, TRIPLE_MAP)) != LORENTZ,
        abs(determinant_3(TRIPLE_MAP)) != 1,
    )
    require(hostile[1] != hostile[2] and hostile[3:] == (True, True),
            ("Pell/Berggren hostiles", hostile))

    fixed_records = tuple(
        (FIXED_NAMES[index],) + rows[index]
        for index in range(len(FIXED_NAMES))
    )
    record = (
        M,
        PELL_TWICE,
        CHANGE,
        TRIPLE_MAP,
        LORENTZ,
        fixed_records,
        rows[-1],
        hostile,
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("renewal Pell semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== THM-3522 renewal Pell/spinor sidecar ==")
    print("renewal=(e,m)->(7e-2m,3e-2m);scalar_recurrence=u[n+2]=5u[n+1]+8u[n]")
    print("pell=(X,Y)=(6e-9m,m);2*(X',Y')=((5,57),(1,5))*(X,Y)")
    print("quadratic_integer=alpha=(5+sqrt(57))/2;Norm(alpha)=-8;(X+Ysqrt57)=6*alpha^n")
    print("invariants=X^2-57Y^2=36*(-8)^n;3e^2-9em+2m^2=3*(-8)^n")
    print("cassini=e[n]m[n+1]-m[n]e[n+1]=3*(-8)^n")
    print("pythagorean_plane=A-9B+5C=3*(-8)^n;triple_map_det=-512;Lorentz_multiplier=64")
    print(f"proved_fixed_prefix_(name,n,e,m,triple)={tuple((name, index, pair[0], pair[1], triple) for name, index, pair, _pell, _norm, triple, _harmonic in fixed_records)}")
    print(f"next_conditional_rows={tuple((index, pair) for index, pair in enumerate(pairs[7:10], start=7))}")
    print(f"hostiles_(wrong_pair,wrong_norm,required_norm,not_isometry,not_unimodular)={hostile}: PASS")
    print(f"semantic_sha256={semantic}")
    print("status=PROVED-ALGEBRAIC arithmetic sidecar;fixed Keller realization proved only through R7")
    print("scope=no next polynomiality,no image,no irreducibility,no all-level Keller theorem,no general JC")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
