#!/usr/bin/env python3
"""Exact companion for THM-2387's H4 elliptic three-isogeny atlas.

The script checks the algebraic identities used by the proof:

* the Cardano norm and reconstruction constants;
* Tate normal form, its order-three flex, the Velu quotient, and both
  j-invariants;
* the binary-quartic invariant convention and sector polynomial; and
* the four-line Weil-pairing tournament together with its switching gauge.

The divisor and connected-cover arguments are geometric and are proved in the
theorem.  Every executable validity check uses ``require`` so optimized mode
cannot remove it.
"""

from __future__ import annotations

from itertools import product

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def weierstrass_invariants(
    a1: sp.Expr,
    a2: sp.Expr,
    a3: sp.Expr,
    a4: sp.Expr,
    a6: sp.Expr,
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Return c4,c6,Delta for a general Weierstrass equation."""

    b2 = sp.expand(a1**2 + 4 * a2)
    b4 = sp.expand(2 * a4 + a1 * a3)
    b6 = sp.expand(a3**2 + 4 * a6)
    b8 = sp.expand(
        a1**2 * a6
        + 4 * a2 * a6
        - a1 * a3 * a4
        + a2 * a3**2
        - a4**2
    )
    c4 = sp.expand(b2**2 - 24 * b4)
    c6 = sp.expand(-b2**3 + 36 * b2 * b4 - 216 * b6)
    discriminant = sp.expand(
        -b2**2 * b8
        - 8 * b4**3
        - 27 * b6**2
        + 9 * b2 * b4 * b6
    )
    require(
        sp.expand(c4**3 - c6**2 - 1728 * discriminant) == 0,
        "Weierstrass invariant identity changed",
    )
    return c4, c6, discriminant


def reduce_on_tate_curve(
    expression: sp.Expr,
    xvar: sp.Symbol,
    yvar: sp.Symbol,
    avar: sp.Symbol,
    bvar: sp.Symbol,
) -> sp.Expr:
    """Reduce a rational expression modulo y^2+a*x*y+b*y-x^3."""

    numerator = sp.together(expression).as_numer_denom()[0]
    relation = yvar**2 + avar * xvar * yvar + bvar * yvar - xvar**3
    remainder = sp.rem(
        sp.Poly(numerator, yvar),
        sp.Poly(relation, yvar),
    )
    return sp.expand(remainder.as_expr())


def symplectic_pair(first: tuple[int, int], second: tuple[int, int]) -> int:
    return (first[0] * second[1] - first[1] * second[0]) % 3


def oriented_edge(
    first_name: str,
    first: tuple[int, int],
    second_name: str,
    second: tuple[int, int],
) -> tuple[str, str]:
    pairing = symplectic_pair(first, second)
    require(pairing in {1, 2}, "distinct projective lines became orthogonal")
    return (
        (first_name, second_name)
        if pairing == 1
        else (second_name, first_name)
    )


def main() -> None:
    # Cardano norm and reconstruction.
    pcap, qcap, scap, hcap, tcap = sp.symbols("P Q S H t")
    acap, abarcap, zcap, wcap = sp.symbols("A Abar z w")
    norm_left = sp.expand(
        (7 * qcap + scap * tcap) * (7 * qcap - scap * tcap)
        + 4 * pcap**3
    )
    norm_relations = sp.expand(
        -scap**2 * (tcap**2 - hcap)
        - (hcap * scap**2 - 4 * pcap**3 - 49 * qcap**2)
    )
    require(
        sp.expand(norm_left - norm_relations) == 0,
        "Cardano norm identity changed",
    )

    cardano_sum = zcap + wcap
    cardano_identity = sp.expand(
        cardano_sum**3
        + 12 * pcap * cardano_sum
        + 4 * (acap + abarcap)
        - (
            (zcap**3 + 4 * acap)
            + (wcap**3 + 4 * abarcap)
            + 3 * cardano_sum * (zcap * wcap + 4 * pcap)
        )
    )
    require(cardano_identity == 0, "Cardano cubic identity changed")
    cube_compatibility = sp.expand(
        (-4 * pcap) ** 3
        + 4 * abarcap * (-4 * acap)
        + 16 * (acap * abarcap + 4 * pcap**3)
    )
    require(
        cube_compatibility == 0,
        "conjugate Cardano cube compatibility changed",
    )

    reconstruction_scale = sp.Rational(2, 1701)
    depressed_p_scale = sp.Rational(16, 964467)
    depressed_q_scale = sp.Rational(64, 703096443)
    require(
        depressed_p_scale * reconstruction_scale
        == 12 * reconstruction_scale**3
        and depressed_q_scale == 56 * reconstruction_scale**3,
        "degree-eighteen Cardano reconstruction constants changed",
    )

    # Tate normal form and the exact 3-isogeny quotient.
    avar, bvar, xvar, yvar, lam = sp.symbols("a b x y lambda")
    source_c4, source_c6, source_delta = weierstrass_invariants(
        avar,
        0,
        bvar,
        0,
        0,
    )
    require(
        sp.expand(source_c4 - avar * (avar**3 - 24 * bvar)) == 0
        and sp.expand(
            source_c6
            - (-avar**6 + 36 * avar**3 * bvar - 216 * bvar**2)
        )
        == 0
        and sp.expand(source_delta - bvar**3 * (avar**3 - 27 * bvar))
        == 0,
        "Tate normal-form invariants changed",
    )
    # Y=0 is the tangent at (0,0), and it meets the cubic with order three.
    tangent_intersection = sp.expand(
        (yvar**2 + avar * xvar * yvar + bvar * yvar - xvar**3).subs(
            yvar,
            0,
        )
    )
    require(
        tangent_intersection == -xvar**3,
        "the marked Tate point stopped being a flex",
    )

    quotient_a4 = -5 * avar * bvar
    quotient_a6 = -bvar * (avar**3 + 7 * bvar)
    quotient_c4, quotient_c6, quotient_delta = weierstrass_invariants(
        avar,
        0,
        bvar,
        quotient_a4,
        quotient_a6,
    )
    require(
        sp.expand(quotient_c4 - avar * (avar**3 + 216 * bvar)) == 0
        and sp.expand(
            quotient_delta - bvar * (avar**3 - 27 * bvar) ** 3
        )
        == 0,
        "Velu quotient invariants changed",
    )

    lambda_expression = avar**3 / bvar
    source_j = sp.cancel(source_c4**3 / source_delta)
    quotient_j = sp.cancel(quotient_c4**3 / quotient_delta)
    require(
        sp.cancel(
            source_j
            - lambda_expression
            * (lambda_expression - 24) ** 3
            / (lambda_expression - 27)
        )
        == 0
        and sp.cancel(
            quotient_j
            - lambda_expression
            * (lambda_expression + 216) ** 3
            / (lambda_expression - 27) ** 3
        )
        == 0,
        "Tate source or quotient j-map changed",
    )

    # Derive the x-coordinate of Velu's quotient from the two translations.
    slope_plus = yvar / xvar
    slope_minus = (yvar + bvar) / xvar
    translated_plus_raw = slope_plus**2 + avar * slope_plus - xvar
    translated_minus_raw = slope_minus**2 + avar * slope_minus - xvar
    translated_plus = -bvar * yvar / xvar**2
    translated_minus = (
        bvar * (yvar + bvar + avar * xvar) / xvar**2
    )
    require(
        reduce_on_tate_curve(
            translated_plus_raw - translated_plus,
            xvar,
            yvar,
            avar,
            bvar,
        )
        == 0
        and reduce_on_tate_curve(
            translated_minus_raw - translated_minus,
            xvar,
            yvar,
            avar,
            bvar,
        )
        == 0,
        "translation-by-three-torsion formulas changed",
    )
    lattes_map = sp.cancel(
        xvar + translated_plus + translated_minus
    )
    expected_lattes = (xvar**3 + avar * bvar * xvar + bvar**2) / xvar**2
    require(
        sp.cancel(lattes_map - expected_lattes) == 0,
        "Velu/Lattes x-map changed",
    )

    # Binary quartic invariants and the incoming-isogeny sector equation.
    h0, h1, h2, h3, h4 = sp.symbols("h0 h1 h2 h3 h4")
    binary_quartic = (
        h0 * xvar**4
        + h1 * xvar**3
        + h2 * xvar**2
        + h3 * xvar
        + h4
    )
    invariant_i = sp.expand(12 * h0 * h4 - 3 * h1 * h3 + h2**2)
    invariant_j = sp.expand(
        72 * h0 * h2 * h4
        + 9 * h1 * h2 * h3
        - 27 * h0 * h3**2
        - 27 * h1**2 * h4
        - 2 * h2**3
    )
    quartic_delta_invariant = sp.expand(4 * invariant_i**3 - invariant_j**2)
    quartic_discriminant = sp.discriminant(binary_quartic, xvar)
    require(
        sp.expand(
            27 * quartic_discriminant - quartic_delta_invariant
        )
        == 0,
        "binary-quartic discriminant convention changed",
    )
    quartic_j_numerator = 6912 * invariant_i**3
    sector_polynomial = sp.expand(
        quartic_j_numerator * (lam - 27) ** 3
        - quartic_delta_invariant * lam * (lam + 216) ** 3
    )
    sector_poly = sp.Poly(sector_polynomial, lam)
    require(
        sector_poly.degree() == 4
        and sp.expand(sector_poly.LC() + quartic_delta_invariant) == 0,
        "smooth binary-quartic sector polynomial stopped having degree four",
    )
    require(
        sp.expand(
            sector_polynomial.subs(lam, 27)
            + quartic_delta_invariant * 27 * 243**3
        )
        == 0,
        "lambda=27 sector guard changed",
    )

    # Positive and hostile controls for the invariant normalization.
    smooth_control = {
        h0: 1,
        h1: 0,
        h2: 0,
        h3: 0,
        h4: -1,
    }
    singular_control = {
        h0: 1,
        h1: 0,
        h2: -2,
        h3: 0,
        h4: 1,
    }
    smooth_i = invariant_i.subs(smooth_control)
    smooth_j = invariant_j.subs(smooth_control)
    smooth_delta = quartic_delta_invariant.subs(smooth_control)
    require(
        smooth_i == -12
        and smooth_j == 0
        and smooth_delta == -6912
        and sp.cancel(6912 * smooth_i**3 / smooth_delta) == 1728,
        "smooth binary-quartic control changed",
    )
    require(
        quartic_delta_invariant.subs(singular_control) == 0
        and sp.Poly(
            sector_polynomial.subs(singular_control),
            lam,
        ).degree()
        < 4,
        "singular binary-quartic hostile control stopped tripping the guard",
    )

    # E[3] has eight oriented nonzero points and four projective lines.
    nonzero_vectors = [
        vector
        for vector in product(range(3), repeat=2)
        if vector != (0, 0)
    ]
    projective_lines = {
        frozenset(
            {
                vector,
                ((-vector[0]) % 3, (-vector[1]) % 3),
            }
        )
        for vector in nonzero_vectors
    }
    require(
        len(nonzero_vectors) == 8 and len(projective_lines) == 4,
        "F3 projective sector count changed",
    )

    vertices = [
        ("infinity", (0, 1)),
        ("0", (1, 0)),
        ("1", (1, 1)),
        ("-1", (1, 2)),
    ]
    edges: list[tuple[str, str]] = []
    for first_index in range(len(vertices)):
        for second_index in range(first_index + 1, len(vertices)):
            edges.append(
                oriented_edge(
                    vertices[first_index][0],
                    vertices[first_index][1],
                    vertices[second_index][0],
                    vertices[second_index][1],
                )
            )
    expected_edges = [
        ("0", "infinity"),
        ("1", "infinity"),
        ("-1", "infinity"),
        ("0", "1"),
        ("-1", "0"),
        ("1", "-1"),
    ]
    require(edges == expected_edges, "Weil-pairing tournament changed")

    edge_set = set(edges)
    for switched_name, switched_vector in vertices:
        switched_vertices = [
            (
                name,
                (
                    ((-vector[0]) % 3, (-vector[1]) % 3)
                    if name == switched_name
                    else vector
                ),
            )
            for name, vector in vertices
        ]
        switched_edges: set[tuple[str, str]] = set()
        for first_index in range(len(switched_vertices)):
            for second_index in range(first_index + 1, len(switched_vertices)):
                switched_edges.add(
                    oriented_edge(
                        switched_vertices[first_index][0],
                        switched_vertices[first_index][1],
                        switched_vertices[second_index][0],
                        switched_vertices[second_index][1],
                    )
                )
        for first_name, second_name in edge_set:
            expected_edge = (
                (second_name, first_name)
                if switched_name in {first_name, second_name}
                else (first_name, second_name)
            )
            require(
                expected_edge in switched_edges,
                "generator negation stopped being a vertex switch",
            )
    inverse_root_edges: set[tuple[str, str]] = set()
    for first_index in range(len(vertices)):
        for second_index in range(first_index + 1, len(vertices)):
            first_name, first_vector = vertices[first_index]
            second_name, second_vector = vertices[second_index]
            pairing = symplectic_pair(first_vector, second_vector)
            inverse_root_edges.add(
                (first_name, second_name)
                if pairing == 2
                else (second_name, first_name)
            )
    reversed_edges = {(second, first) for first, second in edge_set}
    require(
        inverse_root_edges == reversed_edges and len(reversed_edges) == 6,
        "primitive-root inversion stopped reversing the tournament",
    )
    scores = {
        name: sum(first == name for first, _ in edges)
        for name, _ in vertices
    }
    require(
        sorted(scores.values()) == [0, 2, 2, 2],
        "Weil tournament score profile changed",
    )

    # Stable transcript.
    print("THM-2387 degree-18 H4 elliptic three-isogeny atlas exact companion")
    print("Cardano norm: (7Q+St)(7Q-St)=-4P^3")
    print("Cardano reconstruction: v=(2/1701)(z-4P/z)")
    print("divisor ledger: four selected P-lifts, two double infinities, degree 0")
    print("Tate discriminant: b^3(a^3-27b)")
    print("Velu quotient discriminant: b(a^3-27b)^3")
    print("Lattes x-map: (x^3+abx+b^2)/x^2")
    print("j_X=lambda(lambda-24)^3/(lambda-27)")
    print("j_E=lambda(lambda+216)^3/(lambda-27)^3")
    print("binary-quartic convention: 27 Disc(H)=4I^3-J^2")
    print("sector polynomial: degree 4; lambda=27 excluded when H is smooth")
    print("sector counts: 8 oriented points; 4 projective lines")
    print(
        "Weil tournament edges: "
        "0>infinity,1>infinity,-1>infinity,0>1,-1>0,1>-1"
    )
    print("Weil tournament score multiset: 0,2,2,2")
    print("gauge: generator negation switches one star; zeta inversion reverses all")
    print("controls: y^4-1 smooth with j=1728; (y^2-1)^2 rejected as singular")
    print("conclusion: exact isogeny atlas, not an H4 obstruction")


if __name__ == "__main__":
    main()
