#!/usr/bin/env python3
"""Exact hostile probe for prescribed-coordinate orders at the newest Keller prime.

The script has two independent lanes.

1.  Over Q it uses the rational point q0=(2/27,1,1) on V(L), factors the
    other two points in F^{-1}(F(q0)) as one quadratic algebra, and computes
    the three last-step coordinate blocks by exact resultants.  For each of
    y, z, and u=1/x, the degree-nine special-fibre polynomial has exactly one
    repeated linear factor: the double shadow on the unique L=0 block.  Thus
    no cross-block residue collision is forced at the level-two newest prime.

2.  Over finite fields it searches for complete split inverse trees rooted
    at a point q0 on V(L).  It then multiplies all last-step y-, z-, and
    reciprocal-x cubics and checks that gcd(P,P') has degree exactly one.
    These are finite exact witnesses against an identically vanishing
    cross-block resultant at the corresponding characteristic-zero divisor.

The geometric all-level reduction and its scope belong in the companion
reflection/theorem candidate.  Finite-field witnesses are not promoted here
to an all-level noncollision theorem.
"""

from __future__ import annotations

from hashlib import sha256
import json

import sympy as sp


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def l_value(a, b, c):
    return 27 * a**2 * c**2 - 18 * a * b * c + 16 * a + b**3 * c - b**2


def f_map(point):
    x, y, z = point
    unit = 1 + x * y
    return (
        unit**3 * z + y**2 * unit * (4 + 3 * x * y),
        y + 3 * x * unit**2 * z + 3 * x * y**2 * (4 + 3 * x * y),
        2 * x - 3 * x**2 * y - x**3 * z,
    )


def inverse_section(target, x):
    """The rational inverse section from THM-3519, with x left algebraic."""
    a, b, c = target
    denominator = (12 * a - b**2) * x**2 + b * x + 2
    y = b - 3 * a * x * ((9 * a * c - b) * x + 2) / denominator
    z = (2 * x - c - 3 * x**2 * y) / x**3
    return x, sp.cancel(y), sp.cancel(z)


def y_block(target, variable):
    """Constant-leading cubic for the last-step source y-coordinate."""
    a, b, c = target
    r = b - variable
    return sp.expand(27 * a**2 * c - 18 * a * r + 3 * b * r**2 - 2 * r**3)


def z_block(target, variable):
    """Constant-leading cubic Q_z of THM-2546."""
    a, b, c = target
    L = l_value(a, b, c)
    q2 = 324 * a**2 * c**2 - 216 * a * b * c + 408 * a - 15 * b**3 * c + 6 * b**2
    Sz = 27 * a**2 * c**2 - 18 * a * b * c + 52 * a + b**3 * c + 14 * b**2
    Tz = (
        729 * a**4 * c**4
        - 972 * a**3 * b * c**3
        + 2322 * a**3 * c**2
        + 54 * a**2 * b**3 * c**3
        + 270 * a**2 * b**2 * c**2
        - 3735 * a**2 * b * c
        - 338 * a**2
        - 36 * a * b**4 * c**2
        + 122 * a * b**3 * c
        + 1372 * a * b**2
        + b**6 * c**2
        - 2 * b**5 * c
        - 80 * b**4
    )
    return sp.expand(8 * variable**3 + q2 * variable**2 + 6 * L * Sz * variable + L * Tz)


def reciprocal_x_block(target, variable):
    """Integral reciprocal cubic L+(4-3bc)u^2-2cu^3."""
    a, b, c = target
    return sp.expand(l_value(a, b, c) + (4 - 3 * b * c) * variable**2 - 2 * c * variable**3)


def primitive_qq_poly(expression, variable):
    polynomial = sp.Poly(sp.cancel(expression), variable, domain=sp.QQ)
    return polynomial.monic()


def rational_level_two_witness():
    X, Y, Z, U = sp.symbols("X Y Z U")
    q0 = (sp.Rational(2, 27), sp.Integer(1), sp.Integer(1))
    require(l_value(*q0) == 0, "q0 must lie on L")
    eta = tuple(sp.factor(value) for value in f_map(q0))
    outer = sp.factor(l_value(*eta) * X**3 + (4 - 3 * eta[1] * eta[2]) * X - 2 * eta[2])
    quotient, remainder = sp.div(sp.Poly(outer, X), sp.Poly(X - q0[0], X))
    require(remainder.is_zero and quotient.degree() == 2, (outer, quotient, remainder))
    h = quotient.monic().as_expr()

    qx = inverse_section(eta, X)
    require(all(sp.cancel(qx[i].subs(X, q0[0]) - q0[i]) == 0 for i in range(3)), qx)

    rows = {}
    for name, variable, constructor, expected_gcd in (
        ("y", Y, y_block, Y - sp.Rational(1, 3)),
        ("z", Z, z_block, Z),
        ("u", U, reciprocal_x_block, U),
    ):
        block0 = primitive_qq_poly(constructor(q0, variable), variable)
        generic_numerator = sp.together(constructor(qx, variable)).as_numer_denom()[0]
        other_product_raw = sp.resultant(h, generic_numerator, X)
        other_product = primitive_qq_poly(other_product_raw, variable)
        total = primitive_qq_poly(block0.as_expr() * other_product.as_expr(), variable)
        gcd = sp.gcd(total, total.diff()).monic()
        require(total.degree() == 9, (name, "degree", total.degree()))
        require(gcd == sp.Poly(expected_gcd, variable, domain=sp.QQ).monic(), (name, gcd))
        require(sp.gcd(block0, other_product).degree() == 0, (name, "cross-block collision"))
        require(sp.gcd(other_product, other_product.diff()).degree() == 0, (name, "other blocks collide"))
        rows[name] = {
            "block0_factorization": str(sp.factor(block0.as_expr())),
            "other_degree": other_product.degree(),
            "total_degree": total.degree(),
            "derivative_gcd": str(gcd.as_expr()),
            "coefficient_digest": digest([str(c) for c in total.all_coeffs()]),
        }

    return {
        "q0": tuple(str(v) for v in q0),
        "eta": tuple(str(v) for v in eta),
        "outer_factorization": str(outer),
        "quadratic_factor": str(sp.factor(h)),
        "rows": rows,
    }


def f_mod(point, prime):
    x, y, z = point
    unit = (1 + x * y) % prime
    return (
        (unit**3 * z + y**2 * unit * (4 + 3 * x * y)) % prime,
        (y + 3 * x * unit**2 * z + 3 * x * y**2 * (4 + 3 * x * y)) % prime,
        (2 * x - 3 * x**2 * y - x**3 * z) % prime,
    )


def l_mod(point, prime):
    a, b, c = point
    return (27 * a**2 * c**2 - 18 * a * b * c + 16 * a + b**3 * c - b**2) % prime


def poly_mod(expression, variable, prime):
    return sp.Poly(expression, variable, modulus=prime).monic()


def coordinate_collision_row(blocks, prime):
    Y, Z, U = sp.symbols("Y Z U")
    result = {}
    for name, variable, constructor in (
        ("y", Y, y_block),
        ("z", Z, z_block),
        ("u", U, reciprocal_x_block),
    ):
        total = sp.Poly(1, variable, modulus=prime)
        block_polynomials = []
        for block in blocks:
            polynomial = poly_mod(constructor(block, variable), variable, prime)
            require(polynomial.degree() == 3, (prime, name, block, polynomial))
            block_polynomials.append(polynomial)
            total *= polynomial
        total = total.monic()
        gcd = sp.gcd(total, total.diff()).monic()
        pairwise = True
        for i in range(len(block_polynomials)):
            for j in range(i):
                if sp.gcd(block_polynomials[i], block_polynomials[j]).degree() != 0:
                    pairwise = False
                    break
            if not pairwise:
                break
        result[name] = {
            "degree": total.degree(),
            "derivative_gcd_degree": gcd.degree(),
            "derivative_gcd": str(gcd.as_expr()),
            "pairwise_coprime": pairwise,
            "coefficient_digest": digest([int(c) % prime for c in total.all_coeffs()]),
        }
    return result


def split_tree(reverse, target, depth):
    level = [target]
    for _ in range(depth):
        next_level = []
        for point in level:
            fibre = reverse.get(point, ())
            if len(fibre) != 3:
                return None
            next_level.extend(fibre)
        if len(set(next_level)) != len(next_level):
            return None
        level = next_level
    return tuple(sorted(level))


def finite_field_witnesses(prime, depths):
    reverse = {}
    l_points = []
    for x in range(prime):
        for y in range(prime):
            for z in range(prime):
                point = (x, y, z)
                reverse.setdefault(f_mod(point, prime), []).append(point)
                if l_mod(point, prime) == 0:
                    l_points.append(point)

    pending = set(depths)
    witnesses = {}
    for q0 in l_points:
        if not pending:
            break
        # Generic one-step ramified-block gates for y and reciprocal x.
        if q0[0] == 0 or q0[2] == 0 or (4 - 3 * q0[1] * q0[2]) % prime == 0:
            continue
        for depth in sorted(tuple(pending)):
            eta = q0
            for _ in range(depth):
                eta = f_mod(eta, prime)
            blocks = split_tree(reverse, eta, depth)
            if blocks is None or q0 not in blocks:
                continue
            l_blocks = [point for point in blocks if l_mod(point, prime) == 0]
            if l_blocks != [q0]:
                continue
            if any(point[2] == 0 for point in blocks):
                continue
            collision = coordinate_collision_row(blocks, prime)
            if not all(
                row["degree"] == 3 ** (depth + 1)
                and row["derivative_gcd_degree"] == 1
                and row["pairwise_coprime"]
                for row in collision.values()
            ):
                continue
            witnesses[depth] = {
                "level": depth + 1,
                "q0": q0,
                "eta": eta,
                "preceding_blocks": len(blocks),
                "unique_L_block": True,
                "collision": collision,
            }
            pending.remove(depth)
    return {
        "prime": prime,
        "L_point_count": len(l_points),
        "witnesses": witnesses,
        "missing_depths": sorted(pending),
    }


def main() -> None:
    rational = rational_level_two_witness()
    finite_rows = []
    missing = {1, 2, 3}
    for prime in (31, 37, 41, 43, 47, 53, 59, 61):
        if not missing:
            break
        row = finite_field_witnesses(prime, sorted(missing))
        finite_rows.append(row)
        missing -= set(row["witnesses"])

    record = {
        "verdict": "exact cross-block noncollision witnesses; all-level criterion remains separate",
        "rational_level_two": rational,
        "finite_field_rows": finite_rows,
        "unwitnessed_depths": sorted(missing),
    }
    record["semantic_sha256"] = digest(record)
    print(json.dumps(record, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
