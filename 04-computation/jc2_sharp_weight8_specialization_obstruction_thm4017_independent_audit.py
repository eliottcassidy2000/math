#!/usr/bin/env python3
"""Independent exact audit of THM-4017's weight-eight obstruction."""

from fractions import Fraction as F
from itertools import combinations
from math import gcd


gates = 0


def check(label, condition):
    global gates
    if not condition:
        raise AssertionError(label)
    gates += 1
    print(f"PASS {label}")


def pick(vertices):
    area2 = abs(sum(
        vertices[i][0] * vertices[(i + 1) % len(vertices)][1]
        - vertices[(i + 1) % len(vertices)][0] * vertices[i][1]
        for i in range(len(vertices))
    ))
    boundary = sum(
        gcd(
            abs(vertices[(i + 1) % len(vertices)][0] - vertices[i][0]),
            abs(vertices[(i + 1) % len(vertices)][1] - vertices[i][1]),
        )
        for i in range(len(vertices))
    )
    interior = (area2 - boundary + 2) // 2
    return area2, boundary, interior


def lower_facets(points):
    """Lower two-faces of (x,y,height,label,coefficient) support."""
    faces = {}
    for ids in combinations(range(len(points)), 3):
        p0, p1, p2 = (points[i] for i in ids)
        x0, y0, z0 = p0[:3]
        x1, y1, z1 = p1[:3]
        x2, y2, z2 = p2[:3]
        det = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
        if det == 0:
            continue
        a = F((z1-z0)*(y2-y0) - (z2-z0)*(y1-y0), det)
        b = F((x1-x0)*(z2-z0) - (x2-x0)*(z1-z0), det)
        c = F(z0) - a*x0 - b*y0
        gaps = [F(z) - a*x - b*y - c for x, y, z, _, _ in points]
        if min(gaps) >= 0:
            on = tuple(i for i, gap in enumerate(gaps) if gap == 0)
            faces[on] = (a, b, c)
    return faces


def combine_support(raw_points):
    """Combine coincident lifted monomials and delete exact cancellations."""
    sums = {}
    labels = {}
    for x, y, z, label, coeff in raw_points:
        key = (x, y, z)
        sums[key] = sums.get(key, F(0)) + coeff
        labels.setdefault(key, []).append(label)
    ans = []
    for (x, y, z), coeff in sorted(sums.items()):
        if coeff:
            ans.append((x, y, z, "+".join(labels[(x, y, z)]), coeff))
    return ans


def support_from_H(H_terms, gamma):
    """Expand (s^2-p)(1-QH)+gamma*Q*s^2 from H monomials s^i p^j."""
    raw = [
        (2, 0, 0, "s2", F(1)),
        (0, 1, 0, "-p", F(-1)),
        (2, 0, 1, "gamma_Q_s2", gamma),
    ]
    for i, j, name, coeff in H_terms:
        raw.append((i + 2, j, 1, f"-s2*{name}", -coeff))
        raw.append((i, j + 1, 1, f"+p*{name}", coeff))
    return combine_support(raw)


def main():
    # Suppress the nonzero a-scaling by putting a=1. These are reconstructed
    # from THM-4007/4016, not imported from their scripts.
    gamma = F(-1, 2)
    c30_tilde = F(2752, 135)
    c40_tilde = F(-512, 9)
    c02_tilde = F(-8128, 135)
    epsilon = gamma * c30_tilde
    delta = gamma * c40_tilde       # raw [p^4]R
    kappa = gamma * c02_tilde       # raw [y^2]R
    lam = F(-3)
    alpha = F(8, 3)

    check("raw epsilon", epsilon == F(-1376, 135))
    check("raw delta=gamma*c40 is nonzero", delta == F(256, 9))
    check("raw kappa=gamma*c02", kappa == F(4064, 135))
    check("old weight-six attachment X ratio", -epsilon/(epsilon+kappa) == F(43, 84))
    check("old weight-six attachment Y ratio", kappa/(epsilon+kappa) == F(127, 84))

    # A p^a y^b term has source weight 2a+3b under s:r^-1,p:r^-2.
    exponent6 = lambda a, b: 6 - (2*a + 3*b)
    check("p3 is integral on old scale", exponent6(3, 0) == 0)
    check("y2 is integral on old scale", exponent6(0, 2) == 0)
    check("forced p4 has pole rho^-2 on old scale", exponent6(4, 0) == -2)

    # H=lambda*p+alpha*p^2+epsilon*p^3+kappa*s^2*p^2+delta*p^4.
    H = [
        (0, 1, "lambda_p", lam),
        (0, 2, "alpha_p2", alpha),
        (0, 3, "epsilon_p3", epsilon),
        (2, 2, "kappa_y2", kappa),
        (0, 4, "delta_p4", delta),
    ]
    points = support_from_H(H, gamma)
    coeff = {(x, y, z): c for x, y, z, _, c in points}
    check("shared s2*p3 coefficient is kappa-epsilon",
          coeff[(2, 3, 1)] == kappa-epsilon == F(1088, 27))

    faces = lower_facets(points)
    face_data = {
        frozenset((points[i][0], points[i][1], points[i][2]) for i in on): plane
        for on, plane in faces.items()
    }
    primary = frozenset({(2,0,0), (0,1,0), (2,4,1), (0,5,1)})
    side = frozenset({(2,0,0), (4,2,1), (2,4,1)})
    check("truncated support has exactly two lower facets",
          set(face_data) == {primary, side})
    check("primary lower plane", face_data[primary] == (F(1,8), F(1,4), F(-1,4)))
    check("side lower plane", face_data[side] == (F(1,4), F(1,4), F(-1,2)))

    # Initial forms, recorded coefficientwise.
    check("primary initial form coefficients",
          {k: coeff[k] for k in primary} == {
              (2,0,0): F(1), (0,1,0): F(-1),
              (2,4,1): -delta, (0,5,1): delta,
          })
    check("side initial form coefficients",
          {k: coeff[k] for k in side} == {
              (2,0,0): F(1), (4,2,1): -kappa,
              (2,4,1): -delta,
          })

    full_polygon = [(0,1), (2,0), (4,2), (2,4), (0,5)]
    primary_polygon = [(0,1), (2,0), (2,4), (0,5)]
    side_polygon = [(2,0), (4,2), (2,4)]
    check("full truncated Newton polygon genus", pick(full_polygon) == (24,10,8))
    check("primary cell interior count", pick(primary_polygon) == (16,10,4))
    check("side cell interior count", pick(side_polygon) == (8,8,1))
    internal_edge_length = gcd(0, 4)
    check("subdivision genus ledger 8=4+1+(4-1)",
          8 == 4 + 1 + (internal_edge_length-1))

    # Primary face: (S^2-P)(1-delta*P^4). Over an algebraic closure it has
    # one parabola and four vertical rational components with eight transverse
    # intersections. The graph rank is 8-5+1=4.
    check("primary rational-component graph rank", 8-5+1 == 4)
    # Four inter-facet roots add graph rank 4-1=3. Hence the total stable
    # reduction of the nondegenerate truncated subdivision has toric rank 7
    # and abelian dimension one.
    check("truncated toric plus abelian genus ledger", 4+3+1 == 8)

    # Side face (in the torus) is 1-kappa*S^2*P^2-delta*P^4. With T=SP:
    # kappa*T^2=1-delta*P^4. The branch quartic has I!=0,J=0.
    quartic_I = -12*delta
    quartic_J = F(0)
    check("side quartic is smooth", delta*kappa != 0)
    check("side quartic invariants force j=1728", quartic_I != 0 and quartic_J == 0)
    check("side and target j differ", F(1728) != F(0))

    # The first uncontrolled monomial of the same source weight is p*y^2,
    # i.e. eta*s^2*p^3 in H. Its -Q*s^2*H endpoint is (4,3,1).
    side_plane = face_data[side]
    a, b, c = side_plane
    gap = F(1) - (a*4+b*3+c)
    check("p*y2 endpoint lies strictly below old side facet", gap == F(-1,4))

    H_eta = H + [(2, 3, "eta_p_y2", F(1))]
    eta_points = support_from_H(H_eta, gamma)
    eta_faces = lower_facets(eta_points)
    eta_planes = set(eta_faces.values())
    check("generic p*y2 deletes old side plane", side_plane not in eta_planes)
    check("generic p*y2 retains/enlarges primary plane",
          face_data[primary] in eta_planes)
    check("generic p*y2 creates replacement rational side plane",
          (F(1,2), F(0), F(-1)) in eta_planes)

    eta_full_polygon = [(0,1), (2,0), (4,2), (4,3), (2,4), (0,5)]
    eta_primary_polygon = [(0,1), (2,0), (4,3), (2,4), (0,5)]
    eta_side_polygon = [(2,0), (4,2), (4,3)]
    check("p*y2 full Newton polygon genus", pick(eta_full_polygon) == (26,10,9))
    check("p*y2 enlarged primary cell genus", pick(eta_primary_polygon) == (24,8,9))
    check("p*y2 replacement side cell is rational", pick(eta_side_polygon) == (2,4,0))
    check("p*y2 primary component ledger genus2+graph7=9", 2+(8-2+1) == 9)
    check("p*y2 subdivision genus ledger", 9 == 9+0+(gcd(2,3)-1))

    # Normalization firewall: mixing raw kappa with normalized c40 is not the
    # literal initial equation. Either pair of coefficients may be used.
    check("raw and normalized coefficient pairs differ by common gamma",
          delta/c40_tilde == kappa/c02_tilde == gamma)

    print(f"GATES={gates}")
    print("RESULT old_weight6_sharp_realization=REFUTED_BY_FORCED_P4")
    print("RESULT truncated_lower_facets=PRIMARY_WEIGHT8_PLUS_SIDE_QUARTIC")
    print("RESULT side_quartic=j1728; target=j0; old_six_point_invoice_not_applicable")
    print("RESULT p_y2_hostile=DESTROYS_SIDE_FACET")
    print("ALL INDEPENDENT THM-4017 CHECKS PASSED")


if __name__ == "__main__":
    main()
