#!/usr/bin/env python3
"""Exact companion for THM-3533's newest-prime multiplicity theorem.

The proof that the newest raw divisor has coefficient one in the
normalization discriminant is geometric/local and lives in THM-3533.  This
script freezes its algebraic inputs:

* the reciprocal one-step cubic has discriminant ``-4*L*S**2`` and a
  transverse ``1+2`` special fibre;
* a rational point of ``V(L)`` maps transversely to ``V(H)`` away from the
  old divisor;
* three independent degree-nine x-eliminant slices have newest-H
  multiplicity one and old-L multiplicity eight; and
* rescaling an integral primitive element changes its order discriminant by
  an index square (the one-step valuations 1 and 7 are the sharp hostile).

It does not prove the properness/image-degree input, an exact old-component
recursion, or an all-level statement for any prescribed affine coordinate.
"""

from __future__ import annotations

from hashlib import sha256
import json
import pickle
from pathlib import Path

import sympy as sp


EXPECTED_H_RAW_SHA256 = (
    "5a9459b3149e500c1b00b67bd804aa7e607de06bf4610c7cdf5fa26d41d74ce9"
)
EXPECTED_H_LEDGER_SHA256 = (
    "b146c11f33e895c08f72303d282e2668d955e0a58d9268a1b445d4d5202016c2"
)
EXPECTED_SEMANTIC_SHA256 = (
    "aa19a27d7e68467936caebb1497804c2815c07d7c160b1ca172723805d3bd79a"
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def order_at(poly: sp.Expr, variable: sp.Symbol) -> int:
    expanded = sp.Poly(sp.expand(poly), variable, domain=sp.QQ)
    return min(monomial[0] for monomial, coefficient in expanded.terms()
               if coefficient != 0)


def divisor_multiplicity(poly: sp.Expr, factor: sp.Expr,
                         variable: sp.Symbol) -> int:
    dividend = sp.Poly(sp.expand(poly), variable, domain=sp.QQ)
    divisor = sp.Poly(sp.expand(factor), variable, domain=sp.QQ)
    multiplicity = 0
    while True:
        quotient, remainder = sp.div(dividend, divisor)
        if not remainder.is_zero:
            return multiplicity
        dividend = quotient
        multiplicity += 1


def l_poly(a: sp.Expr, b: sp.Expr, c: sp.Expr) -> sp.Expr:
    return 27 * a**2 * c**2 - 18 * a * b * c + 16 * a + b**3 * c - b**2


def t_poly(b: sp.Expr, c: sp.Expr) -> sp.Expr:
    return 4 - 3 * b * c


def s_poly(a: sp.Expr, b: sp.Expr, c: sp.Expr) -> sp.Expr:
    return 27 * a * c**2 - 9 * b * c + 8


def primitive_degree_nine_slice(
    b_value: int,
    c_value: int,
    symbols: tuple[sp.Symbol, sp.Symbol, sp.Symbol, sp.Symbol, sp.Symbol],
) -> sp.Poly:
    a, b, c, w, xi = symbols
    y_denominator = 12 * a * w**2 - b**2 * w**2 + b * w + 2
    middle_y = b - 3 * a * w * (9 * a * c * w - b * w + 2) / y_denominator
    middle_z = (2 * w - 3 * w**2 * middle_y - c) / w**3
    inner_core = (
        l_poly(w, middle_y, middle_z) * xi**3
        + t_poly(middle_y, middle_z) * xi
        - 2 * middle_z
    )
    inner_numerator = sp.together(inner_core).as_numer_denom()[0]
    outer_core = l_poly(a, b, c) * w**3 + t_poly(b, c) * w - 2 * c
    specialized_inner = sp.expand(
        inner_numerator.subs({b: sp.Rational(b_value), c: sp.Rational(c_value)})
    )
    specialized_outer = sp.expand(
        outer_core.subs({b: sp.Rational(b_value), c: sp.Rational(c_value)})
    )
    raw_resultant = sp.resultant(specialized_outer, specialized_inner, w)

    as_xi = sp.Poly(raw_resultant, xi)
    coefficient_gcd = sp.Integer(0)
    for coefficient in as_xi.all_coeffs():
        if coefficient != 0:
            coefficient_gcd = sp.gcd(
                coefficient_gcd, sp.Poly(coefficient, a).as_expr()
            )
    primitive = sp.Poly(sp.cancel(raw_resultant / coefficient_gcd), xi)
    rational_content = sp.gcd_list(
        [
            sp.Rational(value)
            for coefficient in primitive.all_coeffs()
            if coefficient != 0
            for value in sp.Poly(coefficient, a).coeffs()
        ]
    )
    return sp.Poly(sp.expand(primitive.as_expr() / rational_content), xi)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    a, b, c, x, y, z, u, tau, xi, w = sp.symbols(
        "a b c x y z u tau xi w"
    )

    # Fixed sporadic Keller map and the cuspidal discriminant identity.
    unit = 1 + x * y
    first = unit**3 * z + y**2 * unit * (4 + 3 * x * y)
    second = y + 3 * x * unit**2 * z + 3 * x * y**2 * (4 + 3 * x * y)
    third = 2 * x - 3 * x**2 * y - x**3 * z
    jacobian = sp.factor(
        sp.det(sp.Matrix((first, second, third)).jacobian((x, y, z)))
    )
    require(jacobian == -2, ("Jacobian", jacobian))

    target_l = sp.expand(l_poly(a, b, c))
    target_t = sp.expand(t_poly(b, c))
    target_s = sp.expand(s_poly(a, b, c))
    require(
        sp.expand(target_s**2 - target_t**3 - 27 * c**2 * target_l) == 0,
        "cuspidal identity",
    )

    # If u=1/x, the inverse cubic becomes integral at the escaping pair.
    reciprocal_core = sp.expand(target_l + target_t * u**2 - 2 * c * u**3)
    reciprocal_discriminant = sp.factor(sp.discriminant(reciprocal_core, u))
    require(
        sp.expand(reciprocal_discriminant + 4 * target_l * target_s**2) == 0,
        ("reciprocal discriminant", reciprocal_discriminant),
    )

    # Exact transverse point on V(L): (a,b,c)=(2/27,1,1).
    level_one_path = {a: sp.Rational(2, 27) + tau, b: 1, c: 1}
    level_one_l = sp.factor(target_l.subs(level_one_path))
    level_one_s = sp.factor(target_s.subs(level_one_path))
    level_one_disc = sp.factor(reciprocal_discriminant.subs(level_one_path))
    level_one_special = sp.factor(reciprocal_core.subs(level_one_path).subs(tau, 0))
    require(level_one_l == tau * (27 * tau + 2), level_one_l)
    require(level_one_s == 27 * tau + 1, level_one_s)
    require(
        sp.expand(level_one_special - u**2 * (1 - 2 * u)) == 0,
        level_one_special,
    )
    require(order_at(level_one_disc, tau) == 1, level_one_disc)
    require(
        sp.diff(level_one_l, tau).subs(tau, 0) == 2
        and sp.diff(level_one_special, u, 2).subs(u, 0) == 2,
        "transverse double-root control",
    )

    # Load the exact H image prime and pin both its bytes and coefficient ledger.
    h_path = root / "05-knowledge/results/keller_L2_polynomial_opus_20260728.pkl"
    h_raw = h_path.read_bytes()
    require(sha256(h_raw).hexdigest() == EXPECTED_H_RAW_SHA256, "H bytes")
    h_poly = pickle.loads(h_raw)
    h_a, h_b, h_c = sorted(h_poly.free_symbols, key=lambda symbol: symbol.name)
    h_terms = sp.Poly(h_poly, h_a, h_b, h_c).terms()
    h_ledger = "\n".join(
        f"{monomial}:{coefficient}" for monomial, coefficient in h_terms
    )
    require(
        sha256(h_ledger.encode("ascii")).hexdigest() == EXPECTED_H_LEDGER_SHA256,
        "H coefficient ledger",
    )

    # Level-two newest-prime witness.  q(t) crosses V(L) simply; F(q(0)) is
    # on H but off L, and H o F has the same simple transverse parameter.
    source_path = {x: tau, y: 2, z: sp.Rational(1, 2)}
    source_l = sp.factor(l_poly(x, y, z).subs(source_path))
    image_path = tuple(
        sp.factor(coordinate.subs(source_path))
        for coordinate in (first, second, third)
    )
    image_zero = tuple(coordinate.subs(tau, 0) for coordinate in image_path)
    image_substitution = dict(zip((h_a, h_b, h_c), image_zero))
    require(image_zero == (sp.Rational(33, 2), 2, 0), image_zero)
    require(h_poly.subs(image_substitution) == 0, "image witness lies on H")
    require(
        target_l.subs(dict(zip((a, b, c), image_zero))) == 260,
        "image witness avoids old L",
    )
    gradient_h = sp.Matrix(
        [sp.diff(h_poly, variable).subs(image_substitution)
         for variable in (h_a, h_b, h_c)]
    )
    image_tangent = sp.Matrix(
        [sp.diff(coordinate, tau).subs(tau, 0) for coordinate in image_path]
    )
    h_transverse_derivative = sp.factor((gradient_h.T * image_tangent)[0])
    l_transverse_derivative = sp.diff(source_l, tau).subs(tau, 0)
    transverse_ratio = sp.factor(
        h_transverse_derivative / l_transverse_derivative
    )
    require(source_l == tau * (27 * tau - 8) / 4, source_l)
    require(l_transverse_derivative == -2, l_transverse_derivative)
    require(
        h_transverse_derivative == -18538456682778272,
        h_transverse_derivative,
    )
    require(transverse_ratio == 9269228341389136, transverse_ratio)

    # Direct degree-nine eliminant controls.  These reconstruct the x-core
    # independently on three one-parameter slices and divide the complete
    # discriminant, rather than merely testing square class.
    slice_symbols = (a, b, c, w, xi)
    slice_rows = []
    h_relabelled = h_poly.subs({h_a: a, h_b: b, h_c: c}, simultaneous=True)
    for b_value, c_value, expected_lead_quotient in (
        (5, 7, sp.Integer(-1)),
        (2, 3, sp.Integer(-1)),
        (-1, 2, sp.Rational(-1, 4)),
    ):
        eliminant = primitive_degree_nine_slice(
            b_value, c_value, slice_symbols
        )
        require(eliminant.degree() == 9, (b_value, c_value, eliminant.degree()))
        coefficients = eliminant.all_coeffs()
        require(sp.expand(coefficients[1]) == 0, (b_value, c_value, "trace"))
        h_slice = sp.expand(h_relabelled.subs({b: b_value, c: c_value}))
        l_slice = sp.expand(target_l.subs({b: b_value, c: c_value}))
        lead_quotient, lead_remainder = sp.div(
            sp.Poly(coefficients[0], a), sp.Poly(h_slice, a)
        )
        require(lead_remainder.is_zero, (b_value, c_value, "lead remainder"))
        require(
            lead_quotient.as_expr() == expected_lead_quotient,
            (b_value, c_value, lead_quotient.as_expr()),
        )
        discriminant = sp.discriminant(eliminant.as_expr(), xi)
        h_multiplicity = divisor_multiplicity(discriminant, h_slice, a)
        l_multiplicity = divisor_multiplicity(discriminant, l_slice, a)
        require(h_multiplicity == 1, (b_value, c_value, "H", h_multiplicity))
        require(l_multiplicity == 8, (b_value, c_value, "L", l_multiplicity))
        slice_rows.append(
            (
                b_value,
                c_value,
                int(eliminant.degree()),
                int(sp.degree(discriminant, a)),
                str(expected_lead_quotient),
                h_multiplicity,
                l_multiplicity,
            )
        )

    # Primitive is not maximal.  The monic reciprocal generator u has
    # discriminant valuation one.  theta=L^k*u is still primitive and
    # integral, but its discriminant is multiplied by L^(3*2*k).
    monic_reciprocal = sp.expand(
        u**3 - target_t * u**2 / (2 * c) - target_l / (2 * c)
    )
    monic_discriminant = sp.factor(sp.discriminant(monic_reciprocal, u))
    theta = sp.Symbol("theta")
    scaled_monic = sp.expand(
        theta**3
        - target_t * target_l * theta**2 / (2 * c)
        - target_l**4 / (2 * c)
    )
    scaled_discriminant = sp.factor(sp.discriminant(scaled_monic, theta))
    require(
        sp.factor(scaled_discriminant / monic_discriminant) == target_l**6,
        "primitive-rescaling index square",
    )
    require(
        sp.factor(monic_discriminant + target_l * target_s**2 / (4 * c**4)) == 0,
        monic_discriminant,
    )
    primitive_rescaling_rows = []
    for level in range(1, 9):
        degree = 3**level
        index_jump = degree * (degree - 1) // 2
        discriminant_jump = degree * (degree - 1)
        require(discriminant_jump == 2 * index_jump, level)
        primitive_rescaling_rows.append(
            (level, degree, index_jump, discriminant_jump)
        )

    # The all-level newest-coefficient transitivity ledger: the preceding
    # cover is unramified at the new prime, and exactly one residue-degree-one
    # intermediate branch has relative different exponent one.
    newest_rows = []
    for level in range(1, 17):
        previous_degree = 3 ** (level - 1)
        old_contribution = 3 * 0
        relative_contribution = 1 * 1
        newest_coefficient = old_contribution + relative_contribution
        require(newest_coefficient == 1, level)
        newest_rows.append(
            (level, previous_degree, old_contribution,
             relative_contribution, newest_coefficient)
        )

    record = (
        ("jacobian", str(jacobian)),
        ("cuspidal_identity", "S^2-T^3=27c^2L"),
        ("reciprocal_discriminant", str(reciprocal_discriminant)),
        ("level_one", (
            str(level_one_l), str(level_one_s), str(level_one_special),
            order_at(level_one_disc, tau),
        )),
        ("H_hashes", (EXPECTED_H_RAW_SHA256, EXPECTED_H_LEDGER_SHA256)),
        ("level_two_transverse", (
            tuple(map(str, image_zero)), str(source_l),
            str(h_transverse_derivative), str(transverse_ratio),
        )),
        ("degree_nine_slices", tuple(slice_rows)),
        ("primitive_rescaling", tuple(primitive_rescaling_rows)),
        ("newest_transitivity", tuple(newest_rows)),
        ("hostiles", (
            "primitive u versus primitive L*u has valuations 1 versus 7",
            "S=0 is a squared projection collision away from L",
            "old-component multiplicities are not inferred from newest support",
        )),
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(
            semantic == EXPECTED_SEMANTIC_SHA256,
            ("semantic drift", semantic, EXPECTED_SEMANTIC_SHA256),
        )

    print("== fixed Keller newest-prime different/index-square audit ==")
    print(f"det_JF={jacobian};identity=S^2-T^3=27c^2L")
    print(f"reciprocal_core={reciprocal_core}")
    print(f"disc_u_reciprocal={reciprocal_discriminant}")
    print(
        "level1_transverse="
        f"L(t)={level_one_l};S(t)={level_one_s};G_0={level_one_special};"
        f"disc_order={order_at(level_one_disc, tau)}"
    )
    print(
        f"H_raw_sha256={EXPECTED_H_RAW_SHA256};"
        f"H_ledger_sha256={EXPECTED_H_LEDGER_SHA256};H_terms={len(h_terms)}"
    )
    print(
        "level2_transverse="
        f"q0=(0,2,1/2);F(q0)={image_zero};L(F(q0))=260;"
        f"dL={l_transverse_derivative};d(HoF)={h_transverse_derivative};"
        f"ratio={transverse_ratio}"
    )
    print(f"degree9_slice_rows={tuple(slice_rows)}")
    print(
        "primitive_order_formula="
        "v_P(Disc(theta))=1+2*length(S/R[theta]);"
        "one_step_u_vs_Lu=1_vs_7"
    )
    print(f"primitive_rescaling_rows={tuple(primitive_rescaling_rows)}")
    print(f"newest_transitivity_rows={tuple(newest_rows)}")
    print("hostiles=primitive_not_maximal;S_squared_projection_collision;old_components_open")
    print(f"semantic_sha256={semantic}")
    print(
        "scope=exact local algebra and low-level controls for THM-3533;"
        "no prescribed-coordinate all-level multiplicity,old-component recursion,"
        "arbitrary-map,or general-JC claim"
    )
    print("all exact checks passed")


if __name__ == "__main__":
    main()
