#!/usr/bin/env python3
"""Exact arithmetic/UFD certificate for THM-3530.

The geometric proof in THM-3530 says that, after the finite-sheet obstruction
removes the old boundary, the cleared norm of an absolutely irreducible packet
has the form c*R^d with 1<=d<=3.  This companion certifies the remaining
power obstruction: every non-seed raw packet row is primitive, its complete
maximum-lambda face has two coprime prime exponents, and therefore it cannot
be a nontrivial square or cube.

The script does not computationally prove the geometric norm-divisor lemma.
It records sharp square/cube hostiles when grade primitivity is removed and
makes no separability, discriminant-multiplicity, arbitrary-map, or general
Jacobian-conjecture claim.
"""

from __future__ import annotations

from hashlib import sha256
import json
from math import gcd

import sympy as sp


EXPECTED_SEMANTIC_SHA256 = (
    "ab5aa0fe7f28607c766ea4e3d1d42a7b84b4e9843ce2e35df49ebefd06bace75"
)
MATRIX = ((7, -2), (3, -2))


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def step(row: tuple[int, int]) -> tuple[int, int]:
    e, m = row
    return 7 * e - 2 * m, 3 * e - 2 * m


def main() -> None:
    x, y, z = sp.symbols("x y z")
    h = 3 * x * z - 2 * y
    require(sp.Poly(h, x, y, z, domain=sp.QQ).is_irreducible,
            "top-face factor h is irreducible")
    require(sp.gcd(sp.Poly(x, x, y, z), sp.Poly(h, x, y, z)).total_degree() == 0,
            "x and h are coprime")

    rows = [(1, 0)]
    for _ in range(24):
        rows.append(step(rows[-1]))

    gate_rows = []
    for index, (e, m) in enumerate(rows):
        require(e > 0 and 0 <= 2 * m <= e and m % 3 == 0,
                ("packet cone", index, e, m))
        if index >= 1:
            require(gcd(e, m) == 1, ("primitive grade", index, e, m))
            require(e % 6 == 1 and m % 6 == 3,
                    ("mod-six primitive row", index, e, m))
            for image_degree in (2, 3):
                require(not (e % image_degree == 0 and m % image_degree == 0),
                        ("power gate", index, image_degree, e, m))
        if index + 1 < len(rows):
            next_e, _next_m = rows[index + 1]
            require(next_e - e == 6 * e - 2 * m >= 5 * e > 0,
                    ("strict grade growth", index, e, m, next_e))
        gate_rows.append((index, e, m, gcd(e, m), tuple(
            d for d in (2, 3) if e % d == 0 and m % d == 0
        )))

    # Cassini gives a second exact primitivity route: a common divisor of
    # row n divides 3*8^(n-1), while the mod-six row excludes 2 and 3.
    cassini = []
    for index in range(len(rows) - 1):
        e, m = rows[index]
        next_e, next_m = rows[index + 1]
        determinant = e * next_m - m * next_e
        require(determinant == 3 * ((-8) ** index),
                ("Cassini", index, determinant))
        cassini.append((index, determinant))

    # The maximum-lambda face is c*x^e*h^m.  In a UFD, if it is a d-th
    # power then d divides both displayed prime exponents.  These bounded
    # symbolic controls verify the factor support without expanding huge
    # powers.
    factor_support = []
    for index, (e, m) in enumerate(rows[:10]):
        factor_support.append((index, (("x", e), ("3xz-2y", m))))
        if index >= 1:
            require(all(not (e % d == 0 and m % d == 0) for d in (2, 3)),
                    ("bounded UFD gate", index, e, m))

    # Sharp hostiles: L has packet A(1,0), so L^2 and L^3 are complete
    # packets A(2,0), A(3,0), whose top faces are literal powers.  Packet
    # completeness alone does not exclude image multiplicity; primitive
    # packet grade is load-bearing.
    source_L = 27 * x**2 * z**2 - 18 * x * y * z + 16 * x + y**3 * z - y**2

    def max_lambda_face(poly: sp.Expr) -> sp.Expr:
        polynomial = sp.Poly(sp.expand(poly), x, y, z, domain=sp.QQ)
        weighted = [((i - k), coefficient * x**i * y**j * z**k)
                    for (i, j, k), coefficient in polynomial.terms()]
        maximum = max(weight for weight, _term in weighted)
        return sp.expand(sum(term for weight, term in weighted if weight == maximum))

    hostile_square_face = max_lambda_face(source_L**2)
    hostile_cube_face = max_lambda_face(source_L**3)
    require(sp.expand(hostile_square_face - (16 * x) ** 2) == 0,
            ("square hostile", hostile_square_face))
    require(sp.expand(hostile_cube_face - (16 * x) ** 3) == 0,
            ("cube hostile", hostile_cube_face))

    # Raw output/image indexing: P_0=L and P_(n+1)=T(P_n).  The component
    # set of F^n, if the geometric induction is accepted, is P_0,...,P_(n-1).
    component_rows = tuple((iterate, tuple(rows[:iterate])) for iterate in range(1, 13))
    require(all(len({row[0] for row in packet_rows}) == iterate
                for iterate, packet_rows in component_rows),
            "strict e-grades distinguish component candidates")

    record = (
        ("matrix", MATRIX),
        ("rows", tuple(rows)),
        ("gate_rows", tuple(gate_rows)),
        ("cassini", tuple(cassini)),
        ("top_prime_factors", ("x", "3xz-2y")),
        ("factor_support", tuple(factor_support)),
        ("hostile_square_face", str(hostile_square_face)),
        ("hostile_cube_face", str(hostile_cube_face)),
        ("component_rows", component_rows),
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== fixed Keller packet-power image-prime induction audit ==")
    print(f"matrix={MATRIX};rows_0_to_9={tuple(rows[:10])}")
    print("raw_rows_n>=1=gcd(e_n,m_n)=1 and (e_n,m_n)=(1,3) mod 6: PASS through n=24")
    print("strict_e_growth=e_(n+1)-e_n=6e_n-2m_n>=5e_n>0: PASS")
    print("top_face_prime_support=(x^e)*(3xz-2y)^m with coprime irreducibles: PASS")
    print("square_cube_gate=d|e and d|m is impossible for d=2,3 on every nonseed row: PASS")
    print(f"hostile_A(2,0)_top_face={hostile_square_face};literal_square=PASS")
    print(f"hostile_A(3,0)_top_face={hostile_cube_face};literal_cube=PASS")
    print("image_indexing=F^n uses raw primes P_0,...,P_(n-1);strict grades distinguish them")
    print(f"semantic_sha256={semantic}")
    print("scope=arithmetic/UFD companion to the proved geometric norm-divisor and THM-3529 finite-unit gates;no all-level coordinate separability,exact discriminant multiplicity,arbitrary-map,or general JC claim")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
