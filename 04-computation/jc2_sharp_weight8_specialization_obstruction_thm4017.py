#!/usr/bin/env python3
"""Exact checks for THM-4017's sharp weight-eight specialization obstruction."""

from fractions import Fraction
from itertools import combinations
from math import gcd


def pick_interior(vertices):
    twice_area = abs(
        sum(
            vertices[i][0] * vertices[(i + 1) % len(vertices)][1]
            - vertices[(i + 1) % len(vertices)][0] * vertices[i][1]
            for i in range(len(vertices))
        )
    )
    boundary = sum(
        gcd(
            abs(vertices[(i + 1) % len(vertices)][0] - vertices[i][0]),
            abs(vertices[(i + 1) % len(vertices)][1] - vertices[i][1]),
        )
        for i in range(len(vertices))
    )
    interior = Fraction(twice_area - boundary + 2, 2)
    assert interior.denominator == 1
    return twice_area, boundary, int(interior)


def lower_facets(points):
    """Return projected lower faces of lifted points (x,y,height,label)."""
    faces = {}
    for triple in combinations(range(len(points)), 3):
        p1, p2, p3 = (points[i] for i in triple)
        x1, y1, z1, _ = p1
        x2, y2, z2, _ = p2
        x3, y3, z3, _ = p3
        det = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
        if det == 0:
            continue
        a = Fraction(
            (z2 - z1) * (y3 - y1) - (z3 - z1) * (y2 - y1), det
        )
        b = Fraction(
            (x2 - x1) * (z3 - z1) - (x3 - x1) * (z2 - z1), det
        )
        c = Fraction(z1) - a * x1 - b * y1
        gaps = [Fraction(z) - a * x - b * y - c for x, y, z, _ in points]
        if min(gaps) >= 0:
            on = tuple(i for i, gap in enumerate(gaps) if gap == 0)
            faces[on] = (a, b, c)
    return faces


def main():
    # THM-4016 at a=1, so gamma=-1/2 and A5=1. Keep normalized
    # coefficients of R/gamma separate from the raw coefficients of R.
    gamma = Fraction(-1, 2)
    c40_tilde = Fraction(-512, 9)
    c02_tilde = Fraction(-8128, 135)
    c40_raw = gamma * c40_tilde
    epsilon = Fraction(-1376, 135)
    kappa = Fraction(4064, 135)
    assert kappa == gamma * c02_tilde
    assert c40_raw == Fraction(256, 9)
    assert c40_tilde != 0 and c02_tilde != 0 and epsilon * kappa != 0
    assert epsilon + kappa == Fraction(2688, 135)
    assert -epsilon / (epsilon + kappa) == Fraction(43, 84)
    assert kappa / (epsilon + kappa) == Fraction(127, 84)
    print("PASS sharp coefficients and six-point ratios")

    # Under q=rho^-6, p=rho^-2 P, y=rho^-3 SP, a residual monomial
    # p^a y^b enters rho^6 H with exponent 6-(2a+3b).
    def rho_exponent(a, b):
        return 6 - (2 * a + 3 * b)

    assert rho_exponent(3, 0) == 0
    assert rho_exponent(0, 2) == 0
    assert rho_exponent(4, 0) == -2
    print("PASS weight-six scaling: p^3,y^2 integral; forced p^4 has rho exponent -2")

    weight6_polygon = [(0, 1), (2, 0), (4, 2), (2, 3), (0, 4)]
    weight8_polygon = [(0, 1), (2, 0), (4, 2), (2, 4), (0, 5)]
    assert pick_interior(weight6_polygon) == (18, 8, 6)
    assert pick_interior(weight8_polygon) == (24, 10, 8)
    print("PASS Newton polygons: weight-six interior=6; truncated p^4+y^2 interior=8")

    # For Q=q^-1, F_Q=(s^2-p)(1-QH)+gamma*Q*s^2. The lifted
    # coefficient support for H=lambda*p+alpha*p^2+epsilon*p^3
    # +kappa*s^2*p^2+c4*p^4 is listed below. The shared (2,3)
    # coefficient is nonzero for the sharp values, so it remains in support.
    points = [
        (2, 0, 0, "base_s2"),
        (0, 1, 0, "base_p"),
        (2, 1, 1, "lambda_s2p"),
        (0, 2, 1, "lambda_p2"),
        (2, 2, 1, "alpha_s2p2"),
        (0, 3, 1, "alpha_p3"),
        (2, 3, 1, "epsilon/kappa_shared"),
        (0, 4, 1, "epsilon_p4"),
        (4, 2, 1, "kappa_s4p2"),
        (2, 4, 1, "c4_s2p4"),
        (0, 5, 1, "c4_p5"),
    ]
    facets = lower_facets(points)
    labelled = {
        tuple(points[i][3] for i in on): plane for on, plane in facets.items()
    }
    expected = {
        ("base_s2", "base_p", "c4_s2p4", "c4_p5"): (
            Fraction(1, 8),
            Fraction(1, 4),
            Fraction(-1, 4),
        ),
        ("base_s2", "kappa_s4p2", "c4_s2p4"): (
            Fraction(1, 4),
            Fraction(1, 4),
            Fraction(-1, 2),
        ),
    }
    assert labelled == expected
    print("PASS exact lower Newton subdivision has two facets")
    for labels, plane in labelled.items():
        print("  FACET", labels, "PLANE", plane)

    # The second face polynomial is
    # s^2(1-kappa*s^2*p^2-c4*p^4); after T=sp its nontrivial component is
    # kappa*T^2=1-c4*p^4. Its quartic invariant J vanishes and I is nonzero,
    # which is the standard exact j=1728 criterion in characteristic zero.
    # For f(p)=a*p^4+b*p^3+c*p^2+d*p+e:
    # I=12ae-3bd+c^2, J=72ace+9bcd-27ad^2-27b^2e-2c^3.
    a4, b3, c2, d1, e0 = -c40_raw, 0, 0, 0, 1
    quartic_I = 12 * a4 * e0 - 3 * b3 * d1 + c2 * c2
    quartic_J = (
        72 * a4 * c2 * e0
        + 9 * b3 * c2 * d1
        - 27 * a4 * d1 * d1
        - 27 * b3 * b3 * e0
        - 2 * c2 * c2 * c2
    )
    assert quartic_I != 0 and quartic_J == 0
    print("PASS truncated honest side-face elliptic component has j=1728 (I!=0,J=0)")

    # Dual-graph ledger for the truncated two-face model. Face A has the
    # rational component s^2=p and four rational vertical components
    # 1-c4*p^4=0. Each vertical component meets the first one twice.
    # The shared edge has four simple roots, hence four gluing edges to the
    # one elliptic component from face B.
    face_a_vertices = 1 + 4
    face_a_edges = 4 * 2
    face_a_rank = face_a_edges - face_a_vertices + 1
    total_vertices = face_a_vertices + 1
    total_edges = face_a_edges + 4
    total_toric_rank = total_edges - total_vertices + 1
    assert face_a_rank == 4
    assert total_toric_rank == 7
    assert total_toric_rank + 1 == 8
    print("PASS truncated dual graph: face-A toric rank=4; total toric rank=7 plus elliptic genus=1")

    # A later weight-eight p*y^2 term adds the lifted exponent (4,3,1),
    # which lies below the old side-face plane z=(x+y-2)/4 and destroys it.
    a, b, c = expected[("base_s2", "kappa_s4p2", "c4_s2p4")]
    new_point_gap = Fraction(1) - (a * 4 + b * 3 + c)
    assert new_point_gap == Fraction(-1, 4)
    print("PASS hostile: a p*y^2 term destroys the j=1728 side facet; gap=-1/4")

    print("ALL THM-4017 PRIMARY CHECKS PASSED")


if __name__ == "__main__":
    main()
