#!/usr/bin/env python3
"""Finite exact integral-point census in the THM-3646 rank-three subgroup.

This scans the coefficient cube |a|,|b|,|c| <= 20 in <P,Q,R>.  It does not
assume or assert that P,Q,R form a Mordell--Weil basis.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path


B = 1_225_041
LIMIT = 20
P = (Fraction(232), Fraction(3_703))
Q = (Fraction(4_960), Fraction(349_321))
R = (
    Fraction(8_279_053_120, 216_766_729),
    Fraction(-3_611_785_597_108_493, 3_191_456_551_067),
)
Point = tuple[Fraction, Fraction] | None
CHECKS = 0


def require(condition: bool, payload: object) -> None:
    global CHECKS
    if condition is not True:
        raise RuntimeError(payload)
    CHECKS += 1


def digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def on_curve(point: Point) -> bool:
    return point is None or point[1] ** 2 == point[0] ** 3 + B


def negate(point: Point) -> Point:
    return None if point is None else (point[0], -point[1])


def add(left: Point, right: Point) -> Point:
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2 and y1 == -y2:
        return None
    slope = (3 * x1 * x1 / (2 * y1) if x1 == x2
             else (y2 - y1) / (x2 - x1))
    x3 = slope * slope - x1 - x2
    return x3, -y1 + slope * (x1 - x3)


def multiples(point: Point) -> dict[int, Point]:
    positive: list[Point] = [None]
    for _ in range(LIMIT):
        positive.append(add(positive[-1], point))
    require(all(on_curve(value) for value in positive), "multiple table")
    return ({i: positive[i] for i in range(LIMIT + 1)}
            | {-i: negate(positive[i]) for i in range(1, LIMIT + 1)})


def main() -> None:
    require(B == 107**3 - 2, "curve constant")
    require(all(on_curve(point) for point in (P, Q, R)), "basis points")
    mp, mq, mr = (multiples(point) for point in (P, Q, R))
    pq = {(a, b): add(mp[a], mq[b])
          for a in range(-LIMIT, LIMIT + 1)
          for b in range(-LIMIT, LIMIT + 1)}

    integral: dict[tuple[int, int], tuple[int, int, int]] = {}
    evaluated = 0
    for a in range(-LIMIT, LIMIT + 1):
        for b in range(-LIMIT, LIMIT + 1):
            for c in range(-LIMIT, LIMIT + 1):
                if a == b == c == 0:
                    continue
                evaluated += 1
                point = add(pq[a, b], mr[c])
                if point is None:
                    raise RuntimeError(("unexpected torsion relation", a, b, c))
                x, y = point
                if x.denominator != 1 or y.denominator != 1:
                    continue
                coeff = (a, b, c)
                if y < 0:
                    y = -y
                    coeff = tuple(-v for v in coeff)
                key = (int(x), int(y))
                require(on_curve((x, y)), ("integral output", key, coeff))
                old = integral.get(key)
                if old is None or (sum(abs(v) for v in coeff), coeff) < (
                        sum(abs(v) for v in old), old):
                    integral[key] = coeff

    rows = tuple(sorted((x, y, (y - 1) // 2, integral[x, y])
                        for x, y in integral))
    expected = (
        (232, 3_703, 1_851, (1, 0, 0)),
        (4_960, 349_321, 174_660, (0, 1, 0)),
    )
    require(evaluated == (2 * LIMIT + 1) ** 3 - 1 == 68_920,
            "coefficient universe")
    require(rows == expected, ("integral census", rows))
    require(all(x > 0 and y > 0 and y % 2 == 1 for x, y, _r, _c in rows),
            "Berggren positivity")

    source = Path(__file__).resolve()
    source_bytes = source.read_bytes()
    require(b"\r\n" not in source_bytes, "source raw LF")
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source_bytes.decode("utf-8")))),
            "Python assert node present")

    semantic = (B, LIMIT, evaluated, rows)
    print("== fixed-107 rank-three subgroup integral lattice box ==")
    print(f"curve=y^2=x^3+{B};basis=P,Q,R_from_THM-3646")
    print(f"box=[-{LIMIT},{LIMIT}]^3\\{{0}};coefficients={evaluated}")
    print("integral_absolute_points=" + json.dumps(rows, separators=(",", ":")))
    print("positive_odd_ordinates=2;new_fixed107_two_cube_depths=0")
    print(f"semantic_sha256={digest(semantic)}")
    print(f"CHECKS={CHECKS}")
    print("status=FINITE-EXACT SUBGROUP CENSUS")
    print("scope=coefficient box only;not a Mordell-Weil basis or integral-point classification")


if __name__ == "__main__":
    main()
