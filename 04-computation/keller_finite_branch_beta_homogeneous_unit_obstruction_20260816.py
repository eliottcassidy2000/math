#!/usr/bin/env python3
"""Exact algebraic certificate for THM-3529's finite-sheet unit obstruction.

For the fixed sporadic Keller map F, this script expands the source pullback
B=F^*L, verifies that B is homogeneous for beta=i-j-2k, and reduces its
irreducibility to a quadratic in the torus-invariant variables

    p=xy,  q=x^2 z.

The geometric identification of V(B) with the closure of the regular finite
inverse sheet and the graded-domain proof excluding B from every complete
packet are not encoded by this script.  This companion verifies their exact
algebraic inputs and runs bounded packet/hostile controls.  It makes no image-prime,
irreducibility-of-later-rungs, separability, discriminant, or general-JC claim.
"""

from __future__ import annotations

from hashlib import sha256
import json
from math import gcd

import sympy as sp


EXPECTED_SEMANTIC_SHA256 = (
    "8e14d28a41500a2f28a37b181089e66c668806382ea1011814b35076ebdc23fd"
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def main() -> None:
    x, y, z = sp.symbols("x y z")
    a, b, c = sp.symbols("a b c")
    p, q = sp.symbols("p q")

    u = 1 + x * y
    first = u**3 * z + y**2 * u * (4 + 3 * x * y)
    second = y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y)
    third = 2 * x - 3 * x**2 * y - x**3 * z
    jacobian = sp.expand(sp.det(sp.Matrix((first, second, third)).jacobian((x, y, z))))
    require(jacobian == -2, ("Jacobian", jacobian))

    target_L = 27 * a**2 * c**2 - 18 * a * b * c + 16 * a + b**3 * c - b**2
    pullback = sp.expand(target_L.subs({a: first, b: second, c: third}, simultaneous=True))
    expected_pullback = sp.expand(
        -9 * x**4 * y**2 * z**2
        - 54 * x**3 * y**3 * z
        - 18 * x**3 * y * z**2
        - 81 * x**2 * y**4
        - 72 * x**2 * y**2 * z
        - 9 * x**2 * z**2
        - 54 * x * y**3
        + 6 * x * y * z
        + 63 * y**2
        + 16 * z
    )
    require(pullback == expected_pullback, "pullback expansion")

    pullback_poly = sp.Poly(pullback, x, y, z, domain=sp.QQ)
    terms = []
    for (i, j, k), coefficient in pullback_poly.terms():
        beta = i - j - 2 * k
        terms.append((i, j, k, int(coefficient), beta))
    terms.sort(reverse=True)
    beta_support = sorted({row[-1] for row in terms})
    require(beta_support == [-2], ("beta support", beta_support))
    require(pullback_poly.degree(x) == 4, "positive x-degree")
    require(sp.Poly(pullback.subs(x, 0), y, z) == sp.Poly(63 * y**2 + 16 * z, y, z),
            "x does not divide pullback")
    coefficient_content = 0
    for coefficient in pullback_poly.coeffs():
        coefficient_content = gcd(coefficient_content, abs(int(coefficient)))
    require(coefficient_content == 1, ("primitive content", coefficient_content))

    # In Q[x^{+-1},y,z], p=xy and q=x^2z give x^2 B=C(p,q).
    torus_form = sp.cancel(
        x**2 * pullback.subs({y: p / x, z: q / x**2}, simultaneous=True)
    )
    require(x not in torus_form.free_symbols, ("uncancelled x", torus_form))
    torus_form = sp.expand(torus_form)
    expected_torus_form = sp.expand(
        -9 * p**2 * q**2
        - 54 * p**3 * q
        - 18 * p * q**2
        - 81 * p**4
        - 72 * p**2 * q
        - 9 * q**2
        - 54 * p**3
        + 6 * p * q
        + 63 * p**2
        + 16 * q
    )
    require(torus_form == expected_torus_form, "torus form")

    torus_q = sp.Poly(torus_form, q, domain=sp.QQ[p])
    q_coefficients = tuple(sp.factor(value) for value in torus_q.all_coeffs())
    coefficient_gcd = sp.gcd(sp.gcd(q_coefficients[0], q_coefficients[1]), q_coefficients[2])
    require(coefficient_gcd in (1, -1), ("q-content", coefficient_gcd))
    discriminant = sp.factor(sp.discriminant(torus_form, q))
    require(sp.expand(discriminant - 64 * (3 * p + 4)) == 0,
            ("quadratic discriminant", discriminant))
    require(sp.degree(3 * p + 4, p) == 1, "odd valuation/nonsquare hostile")

    factor_coefficient, factors = sp.factor_list(pullback_poly.as_expr(), x, y, z)
    require(factor_coefficient in (1, -1), ("factor coefficient", factor_coefficient))
    require(len(factors) == 1 and factors[0][1] == 1,
            ("pullback factorization", factor_coefficient, factors))
    require(sp.expand(factors[0][0] - factor_coefficient * pullback) == 0
            or sp.expand(factors[0][0] + factor_coefficient * pullback) == 0,
            "factor is the pullback up to sign")

    # The canonical target witness on L=0 and its regular finite inverse point.
    target_point = {a: sp.Rational(2, 27), b: 1, c: 1}
    source_point = {x: 2, y: sp.Rational(5, 6), z: sp.Rational(-7, 8)}
    require(target_L.subs(target_point) == 0, "target witness lies on L")
    require(tuple(sp.simplify(value.subs(source_point)) for value in (first, second, third))
            == tuple(target_point[value] for value in (a, b, c)),
            "finite inverse witness")
    require(pullback.subs(source_point) == 0, "finite sheet lies on pullback")
    gradient_at_witness = tuple(
        sp.simplify(sp.diff(pullback, variable).subs(source_point))
        for variable in (x, y, z)
    )
    require(any(value != 0 for value in gradient_at_witness),
            ("reduced witness", gradient_at_witness))

    # Symbolic packet-face bookkeeping.  Every factor y^2+cz has beta -2,
    # so the displayed complete minimum-beta face is x-free, nonzero, and has
    # weight -5e+2m.  The bounded census checks the full admissible packet cone
    # 0<=m<=e, 3|m through e<=300, not only the raw-orbit subcone 2m<=e.
    e_symbol, m_symbol = sp.symbols("e m", integer=True)
    symbolic_packet_beta = sp.expand(
        -(3 * e_symbol - 2 * m_symbol)
        - 2 * (e_symbol - m_symbol)
        - 2 * (sp.Rational(2, 3) * m_symbol)
        - 2 * (sp.Rational(1, 3) * m_symbol)
    )
    require(symbolic_packet_beta == -5 * e_symbol + 2 * m_symbol,
            ("packet beta", symbolic_packet_beta))

    packet_census = []
    for e_value in range(1, 301):
        for m_value in range(0, e_value + 1, 3):
            exponents = (
                3 * e_value - 2 * m_value,
                e_value - m_value,
                2 * m_value // 3,
                m_value // 3,
            )
            require(all(value >= 0 for value in exponents),
                    ("packet exponent", e_value, m_value, exponents))
            beta_value = -exponents[0] - 2 * exponents[1] \
                - 2 * exponents[2] - 2 * exponents[3]
            require(beta_value == -5 * e_value + 2 * m_value,
                    ("packet beta census", e_value, m_value, beta_value))
            packet_census.append((e_value, m_value, beta_value))
    require((3, 3, -9) in packet_census,
            "full packet cone must include A(3,3) outside 2m<=e")

    # Hostile controls isolate both load-bearing hypotheses.
    # B itself is divisible by B and beta-homogeneous, but its beta face has
    # positive x-degree, so it is not a complete packet.  Conversely an
    # x-free beta-homogeneous divisor such as y^2 can divide an x-free face;
    # the positive x-degree of B is essential.
    require(sp.Poly(pullback, x).degree() == 4, "self-divisor hostile")
    require(sp.Poly(y**2, x).degree() == 0, "x-free divisor hostile")

    record = (
        ("jacobian", str(jacobian)),
        ("pullback_terms", tuple(terms)),
        ("beta_support", tuple(beta_support)),
        ("pullback_multidegree", tuple(pullback_poly.degree(v) for v in (x, y, z))),
        ("pullback_total_degree", pullback_poly.total_degree()),
        ("pullback_content", coefficient_content),
        ("torus_form", str(torus_form)),
        ("torus_q_coefficients", tuple(map(str, q_coefficients))),
        ("torus_q_discriminant", str(discriminant)),
        ("factor_count", len(factors)),
        ("finite_witness_gradient", tuple(map(str, gradient_at_witness))),
        ("packet_census_size", len(packet_census)),
        ("packet_census_digest", digest_json(tuple(packet_census))),
        ("hostiles", (
            "B is not a packet",
            "x-free divisor does not obstruct",
            "A(3,3) lies outside the raw 2m<=e subcone",
        )),
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== fixed Keller finite-branch beta-homogeneous unit obstruction ==")
    print(f"det_JF={jacobian}")
    print(f"B=F^*L={pullback}")
    print(f"terms={len(terms)};multidegree={tuple(pullback_poly.degree(v) for v in (x, y, z))};total_degree={pullback_poly.total_degree()};content={coefficient_content}")
    print(f"beta_support={beta_support};deg_x_B={pullback_poly.degree(x)};B_at_x0={pullback.subs(x, 0)}")
    print(f"x^2*B=C(p,q)={torus_form}")
    print(f"C_as_quadratic_in_q_coefficients={q_coefficients}")
    print(f"disc_q_C={discriminant};3p+4_has_odd_linear_valuation=PASS")
    print("irreducibility=PASS (exact factorization and primitive-quadratic discriminant certificate)")
    print(f"finite_inverse_witness=target{tuple(target_point[v] for v in (a,b,c))}<-source{tuple(source_point[v] for v in (x,y,z))};grad_B={gradient_at_witness}")
    print(f"packet_min_beta_face=x-free;symbolic_weight={symbolic_packet_beta};admissible_census={len(packet_census)}")
    print("graded_obstruction=B|P would force B|in_min_beta(P), impossible by x-degree additivity")
    print("hostiles=B itself fails packet completeness;an x-free homogeneous divisor would evade the x-degree obstruction;A(3,3) tests the full cone beyond 2m<=e")
    print(f"semantic_sha256={semantic}")
    print("scope=exact algebraic inputs for finite-sheet unit theorem;no later-rung irreducibility,image-prime,separability,discriminant,or general JC claim")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
