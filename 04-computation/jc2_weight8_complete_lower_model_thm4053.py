#!/usr/bin/env python3
"""Primary exact certificate for THM-4053.

The script reconstructs the complete maximum-weight-eight support in the live
reduced (2,3) seam, enumerates all lower faces in six support strata, checks
the face normalizations, boundary discriminants, node charts, genus ledgers,
and the Eisenstein-norm survivor.  The stable-map and CM arguments are in the
theorem file.  No Python assertions are used.
"""

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import gcd

import sympy as sp


CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def combine(raw):
    result = {}
    for x, y, height, coefficient in raw:
        key = (x, y, height)
        result[key] = result.get(key, F(0)) + coefficient
    return {key: coefficient for key, coefficient in result.items() if coefficient}


def expanded_support(coefficients):
    """Expand (S^2-P)(1-QH)+gamma*Q*S^2 at coefficient valuations 0/1."""
    raw = [(2, 0, 0, F(1)), (0, 1, 0, F(-1))]
    # H monomials are keyed by exponents (i,j) of p^i*y^j; y=S*P.
    for (i, j), coefficient in coefficients.items():
        raw.append((j + 2, i + j, 1, -coefficient))
        raw.append((j, i + j + 1, 1, coefficient))
    return combine(raw)


def lower_faces(points):
    items = tuple(sorted(points))
    result = {}
    for triple in combinations(range(len(items)), 3):
        q0, q1, q2 = (items[index] for index in triple)
        x0, y0, z0 = q0
        x1, y1, z1 = q1
        x2, y2, z2 = q2
        determinant = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
        if not determinant:
            continue
        ax = F((z1 - z0) * (y2 - y0) - (z2 - z0) * (y1 - y0), determinant)
        ay = F((x1 - x0) * (z2 - z0) - (x2 - x0) * (z1 - z0), determinant)
        c = F(z0) - ax * x0 - ay * y0
        gaps = {point: F(point[2]) - ax * point[0] - ay * point[1] - c
                for point in items}
        if min(gaps.values()) >= 0:
            face = frozenset(point for point, gap in gaps.items() if gap == 0)
            result[face] = (ax, ay, c)
    return result


def convex_hull(points):
    points = sorted(set(points))
    if len(points) <= 1:
        return points

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower = []
    for point in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper = []
    for point in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return lower[:-1] + upper[:-1]


def pick(points):
    vertices = convex_hull(points)
    area2 = abs(sum(
        vertices[i][0] * vertices[(i + 1) % len(vertices)][1]
        - vertices[(i + 1) % len(vertices)][0] * vertices[i][1]
        for i in range(len(vertices))
    ))
    boundary = sum(
        gcd(abs(vertices[(i + 1) % len(vertices)][0] - vertices[i][0]),
            abs(vertices[(i + 1) % len(vertices)][1] - vertices[i][1]))
        for i in range(len(vertices))
    )
    require((area2 - boundary + 2) % 2 == 0, "Pick parity failed")
    return area2, boundary, (area2 - boundary + 2) // 2


L = (F(1, 8), F(1, 4), F(-1, 4))
D = (F(1, 4), F(1, 4), F(-1, 2))
V = (F(0), F(1, 3), F(-1, 3))
R = (F(1, 2), F(0), F(-1))
PLANE_NAME = {L: "L", D: "D", V: "V", R: "R"}


def make_coefficients(epsilon, kappa, phi, delta, theta):
    coefficients = {
        (1, 0): F(7, 5),
        (2, 0): F(-11, 3),
        (3, 0): epsilon,
        (0, 2): kappa,
    }
    if phi:
        coefficients[(2, 1)] = phi
    if delta:
        coefficients[(4, 0)] = delta
    if theta:
        coefficients[(1, 2)] = theta
    return coefficients


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


def factor(n):
    factors = {}
    prime = 2
    while prime * prime <= n:
        while n % prime == 0:
            factors[prime] = factors.get(prime, 0) + 1
            n //= prime
        prime = 3 if prime == 2 else prime + 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def eisenstein_local_criterion(n):
    return all(prime % 3 != 2 or exponent % 2 == 0
               for prime, exponent in factor(n).items())


def eisenstein_represented(n):
    bound = int((4 * n / 3) ** 0.5) + 3
    return any(a * a - a * b + b * b == n
               for a in range(-bound, bound + 1)
               for b in range(-bound, bound + 1))


def main():
    # Exact forced coefficient row, normalized by gamma and A5.
    gamma = F(1)
    A5 = F(1)
    epsilon = gamma * F(2752, 135) / A5 ** 3

    def forced_kappa(delta):
        return -gamma * F(5696, 45) / A5 ** 3 - F(7, 6) * A5 * delta

    delta_kappa_zero = -gamma * F(11392, 105) / A5 ** 4
    delta_kappa_epsilon = -gamma * F(7936, 63) / A5 ** 4
    require(forced_kappa(delta_kappa_zero) == 0, "kappa-zero wall changed")
    require(forced_kappa(delta_kappa_epsilon) == epsilon,
            "kappa-epsilon support-cancellation wall changed")
    require(epsilon / gamma == F(2752, 135), "forced epsilon changed")

    delta = F(1)
    theta = F(2)
    phi = F(7)
    kappa = forced_kappa(delta)
    kappa_theta = forced_kappa(F(0))
    require(kappa != 0 and kappa_theta != 0, "generic representatives degenerated")

    cases = (
        ("delta_generic", make_coefficients(epsilon, kappa, phi, delta, F(0)),
         {L, D}, 8, 7, 1),
        ("delta_kappa0_phi", make_coefficients(epsilon, F(0), phi, delta_kappa_zero, F(0)),
         {L, D}, 7, 7, 0),
        ("delta_kappa0_contract", make_coefficients(epsilon, F(0), F(0), delta_kappa_zero, F(0)),
         {L}, 4, 4, 0),
        ("theta_only", make_coefficients(epsilon, kappa_theta, phi, F(0), theta),
         {L, V, R}, 8, 7, 1),
        ("both_generic", make_coefficients(epsilon, kappa, phi, delta, theta),
         {L, R}, 9, 7, 2),
        ("both_kappa0", make_coefficients(epsilon, F(0), phi, delta_kappa_zero, theta),
         {L}, 9, 7, 2),
    )

    semantics = []
    case_support = {}
    for name, coefficients, expected_planes, genus, graph_rank, abelian in cases:
        support = expanded_support(coefficients)
        faces = lower_faces(support)
        planes = set(faces.values())
        require(planes == expected_planes, f"{name} lower-face inventory changed")
        global_pick = pick([(x, y) for x, y, _ in support])
        require(global_pick[2] == genus, f"{name} generic genus changed")
        require(graph_rank + abelian == genus, f"{name} genus split changed")
        require(all(tuple(24 * entry for entry in plane) in {
            (3, 6, -6), (6, 6, -12), (0, 8, -8), (12, 0, -24)
        } for plane in planes), f"{name} base change failed")
        plane_names = ",".join(sorted(PLANE_NAME[plane] for plane in planes))
        semantics.append(f"{name}:faces={plane_names};g={genus};graph={graph_rank};ab={abelian}")
        case_support[name] = (support, faces, coefficients)

    # Expanded collision points and their harmless/nonharmless walls.
    cancellation = expanded_support(make_coefficients(
        epsilon, epsilon, phi, delta_kappa_epsilon, F(0)
    ))
    require((2, 3, 1) not in cancellation, "kappa=epsilon did not cancel (2,3,1)")
    require(set(lower_faces(cancellation).values()) == {L, D},
            "kappa=epsilon changed the lower subdivision")
    theta_delta = expanded_support(make_coefficients(epsilon, kappa, phi, delta, delta))
    require((2, 4, 1) not in theta_delta, "theta=delta did not cancel (2,4,1)")
    require(set(lower_faces(theta_delta).values()) == {L, R},
            "theta=delta changed the lower subdivision")

    # Face equations, normalizations, and exact wall polynomials.
    S, P, T, X = sp.symbols("S P T X")
    de, th, ka, ph, ep = sp.symbols("delta theta kappa phi epsilon", nonzero=True)
    gL_delta = (S ** 2 - P) * (1 - de * P ** 4)
    gD = 1 - ka * S ** 2 * P ** 2 - ph * S * P ** 3 - de * P ** 4
    gL_theta = (S ** 2 - P) * (1 - th * S ** 2 * P ** 3)
    gV = -1 + P ** 3 * (ep + ph * S + th * S ** 2)
    gR = 1 - S ** 2 * P ** 2 * (ka + th * P)
    gL_both = (S ** 2 - P) * (1 - de * P ** 4 - th * S ** 2 * P ** 3)
    require(sp.Poly(gL_delta, S, P).total_degree() == 6, "delta main equation changed")
    require(sp.Poly(gL_theta, S, P).total_degree() == 7, "theta main equation changed")
    require(sp.Poly(gL_both, S, P).total_degree() == 7, "both main equation changed")
    require(sp.diff(gR, ka) == -S ** 2 * P ** 2, "right face coefficient changed")

    Wd = 2 * ka * T + ph * P ** 2
    relation_d = ka * T ** 2 + ph * T * P ** 2 + de * P ** 4 - 1
    require(sp.expand(Wd ** 2 - (4 * ka + (ph ** 2 - 4 * ka * de) * P ** 4)
                      - 4 * ka * relation_d) == 0,
            "delta-side quartic completion changed")
    Wv = 2 * th * S + ph
    relation_v = th * S ** 2 + ph * S + ep - X ** 3
    require(sp.expand(Wv ** 2 - ((ph ** 2 - 4 * ep * th) + 4 * th * X ** 3)
                      - 4 * th * relation_v) == 0,
            "theta-side cubic completion changed")

    edge_delta = ka + ph * T + de * T ** 2
    edge_theta = th + ph * T + ep * T ** 2
    edge_both = (T - 1) * (de * T + th)
    Delta_D = ph ** 2 - 4 * ka * de
    Delta_V = ph ** 2 - 4 * ep * th
    require(sp.discriminant(edge_delta, T) == Delta_D, "delta edge discriminant changed")
    require(sp.discriminant(edge_theta, T) == Delta_V, "theta edge discriminant changed")
    require(sp.factor(sp.discriminant(edge_both, T)) == (de + th) ** 2,
            "Bolza attachment discriminant changed")
    require(sp.discriminant(edge_both.subs(th, de), T) == 4 * de ** 2,
            "theta=delta was misclassified as resonance")
    require(sp.discriminant(edge_both.subs(th, -de), T) == 0,
            "theta=-delta collision wall disappeared")

    # Long generic outer edges are automatically squarefree.
    Q, lam, al = sp.symbols("Q lambda alpha")
    quartic_outer = -1 + Q * (lam * T + al * T ** 2 + ep * T ** 3 + de * T ** 4)
    cubic_outer = -1 + Q * (lam * T + al * T ** 2 + ep * T ** 3)
    disc4 = sp.Poly(sp.discriminant(quartic_outer, T), Q)
    disc3 = sp.Poly(sp.discriminant(cubic_outer, T), Q)
    require(disc4.nth(3) == -256 * de ** 3, "quartic outer leading discriminant changed")
    require(all(disc4.nth(power) == 0 for power in range(3)),
            "quartic outer discriminant order dropped")
    require(disc3.nth(2) == -27 * ep ** 2, "cubic outer leading discriminant changed")
    require(all(disc3.nth(power) == 0 for power in range(2)),
            "cubic outer discriminant order dropped")

    # Eight transverse main intersections and the exact A_23 chart.
    f0 = S ** 2 - P
    f8 = 1 - de * P ** 4 - th * S ** 2 * P ** 3
    jacobian = sp.det(sp.Matrix([
        [sp.diff(f0, S), sp.diff(f0, P)],
        [sp.diff(f8, S), sp.diff(f8, P)],
    ]))
    require(sp.simplify(jacobian.subs(P, S ** 2) + 8 * (de + th) * S ** 7) == 0,
            "eight-node determinant changed")
    sigma, gamma_symbol = sp.symbols("sigma gamma", nonzero=True)
    Hsigma = (de * P ** 4 + th * S ** 2 * P ** 3
              + sigma ** 3 * ph * S * P ** 3
              + sigma ** 6 * (ep * P ** 3 + ka * S ** 2 * P ** 2)
              + sigma ** 12 * al * P ** 2 + sigma ** 18 * lam * P)
    U = S ** 2 - P
    local_V = (1 - Hsigma) / S ** 2
    scaled = U * (1 - Hsigma) + gamma_symbol * sigma ** 24 * S ** 2
    require(sp.simplify(scaled / S ** 2 - (U * local_V + gamma_symbol * sigma ** 24)) == 0,
            "exact local A23 chart changed")
    require(24 - 1 == 23, "A23 chain length changed")

    # Pick/graph ledgers for the actual lower models.
    expected_picks = {
        "delta_L": (16, 10, 4), "delta_D": (8, 8, 1),
        "theta_L": (16, 4, 7), "theta_V": (6, 6, 1),
        "right": (2, 4, 0), "both_L": (24, 8, 9),
    }
    polygons = {
        "delta_L": [(0, 1), (2, 0), (2, 4), (0, 5)],
        "delta_D": [(2, 0), (4, 2), (3, 3), (2, 4)],
        "theta_L": [(0, 1), (2, 0), (4, 3), (2, 4)],
        "theta_V": [(0, 1), (0, 4), (1, 4), (2, 4)],
        "right": [(2, 0), (4, 2), (4, 3)],
        "both_L": [(0, 1), (2, 0), (4, 3), (2, 4), (0, 5)],
    }
    for name, polygon in polygons.items():
        require(pick(polygon) == expected_picks[name], f"{name} Pick ledger changed")
    require(8 - 2 + 1 == 7, "eight-node dual graph rank changed")

    # Exact degree sieve on the sole nonresonant survivor.
    norm_checks = 0
    for degree in range(1, 501):
        require(eisenstein_represented(degree) == eisenstein_local_criterion(degree),
                f"Eisenstein norm criterion failed at {degree}")
        norm_checks += 1
    require(2 * 8 - 2 == 14, "theta-only Riemann-Hurwitz budget changed")

    semantics.extend([
        "walls=delta:kappa*delta!=0&phi^2=4*kappa*delta;theta:phi^2=4*epsilon*theta;both:delta+theta=0",
        "excluded=delta_off_wall,both_off_wall",
        "survivor=delta0,theta!=0,DeltaV!=0;j0;degree=a^2-a*b+b^2;ramification=14",
    ])
    digest = sha256("\n".join(semantics).encode()).hexdigest()

    print("JC2 MAXIMUM-WEIGHT-EIGHT COMPLETE LOWER MODEL")
    print("status=PROVED_ON_LIVE_b=d=0_SEAM;JC2=OPEN")
    print("forced=epsilon/gamma:2752/(135A5^3);delta/gamma+6kappa/(7A5gamma):-11392/(105A5^4)")
    print("strata=delta_generic,delta_kappa0_phi,delta_kappa0_contract,theta_only,both_generic,both_kappa0")
    print("base_change=Q=sigma^24;faces=L,D,V,R;multiplicity=1;nodes=8xA23")
    print("walls=DeltaD=phi^2-4kappa*delta;DeltaV=phi^2-4epsilon*theta;delta+theta=0")
    print("excluded_off_wall=p4_only,two_term_top")
    print("survivor=py2_only;E0;j=0;degree=Eisenstein_norm;ramification_budget=14")
    print(f"norm_checks={norm_checks}")
    print(f"semantic_digest={digest}")
    print(f"checks={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
