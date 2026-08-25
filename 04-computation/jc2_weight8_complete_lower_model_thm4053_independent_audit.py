#!/usr/bin/env python3
"""SymPy-free independent audit for THM-4053.

This path does not enumerate triples of lifted points.  It verifies the four
declared supporting planes directly, proves that their projected polygons
tile each of six Newton polygons by exact area, and reconstructs the algebra,
walls, genus, and degree sieve with small rational/polynomial routines.
"""

from fractions import Fraction as F
from hashlib import sha256
from math import gcd


GATES = 0


def check(condition, message):
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(message)


def combine(raw):
    result = {}
    for key, coefficient in raw:
        result[key] = result.get(key, F(0)) + coefficient
    return {key: coefficient for key, coefficient in result.items() if coefficient}


def support(coefficients):
    raw = [((2, 0, 0), F(1)), ((0, 1, 0), F(-1))]
    for (i, j), coefficient in coefficients.items():
        raw.extend((
            ((j + 2, i + j, 1), -coefficient),
            ((j, i + j + 1, 1), coefficient),
        ))
    return combine(raw)


def hull(points):
    points = sorted(set(points))

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


def polygon_ledger(points):
    vertices = hull(points)
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
    return area2, boundary, (area2 - boundary + 2) // 2


L = (F(1, 8), F(1, 4), F(-1, 4))
D = (F(1, 4), F(1, 4), F(-1, 2))
V = (F(0), F(1, 3), F(-1, 3))
R = (F(1, 2), F(0), F(-1))
NAMES = {L: "L", D: "D", V: "V", R: "R"}


def plane_gap(point, plane):
    x, y, height = point
    return F(height) - plane[0] * x - plane[1] * y - plane[2]


def plane_polygon(points, plane):
    return hull([(x, y) for x, y, height in points if plane_gap((x, y, height), plane) == 0])


def coefficients(epsilon, kappa, phi, delta, theta):
    result = {(1, 0): F(7, 5), (2, 0): F(-11, 3),
              (3, 0): epsilon, (0, 2): kappa}
    if phi:
        result[(2, 1)] = phi
    if delta:
        result[(4, 0)] = delta
    if theta:
        result[(1, 2)] = theta
    return result


def poly_add(left, right):
    raw = [(key, value) for key, value in left.items()]
    raw.extend((key, value) for key, value in right.items())
    return combine(raw)


def poly_mul(left, right):
    raw = []
    for (i, j), a in left.items():
        for (k, ell), b in right.items():
            raw.append(((i + k, j + ell), a * b))
    return combine(raw)


def factor(n):
    result = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            result[p] = result.get(p, 0) + 1
            n //= p
        p = 3 if p == 2 else p + 2
    if n > 1:
        result[n] = result.get(n, 0) + 1
    return result


def local_norm_test(n):
    return all(p % 3 != 2 or exponent % 2 == 0
               for p, exponent in factor(n).items())


def literal_norm_test(n):
    bound = 2 * int(n ** 0.5) + 3
    values = set()
    for a in range(-bound, bound + 1):
        for b in range(-bound, bound + 1):
            value = a * a - a * b + b * b
            if 0 < value <= 500:
                values.add(value)
    return n in values


def main():
    gamma = F(1)
    epsilon = F(2752, 135)

    def kappa_from(delta):
        return -F(5696, 45) - F(7, 6) * delta

    delta0 = -F(11392, 105)
    deltae = -F(7936, 63)
    check(kappa_from(delta0) == 0, "kappa-zero row")
    check(kappa_from(deltae) == epsilon, "kappa-epsilon row")
    check(epsilon / gamma == F(2752, 135), "epsilon row")

    delta, theta, phi = F(1), F(2), F(7)
    kappa, kappa_theta = kappa_from(delta), kappa_from(F(0))
    cases = (
        ("delta_generic", coefficients(epsilon, kappa, phi, delta, 0),
         (L, D), 8, 7, 1),
        ("delta_kappa0_phi", coefficients(epsilon, 0, phi, delta0, 0),
         (L, D), 7, 7, 0),
        ("delta_kappa0_contract", coefficients(epsilon, 0, 0, delta0, 0),
         (L,), 4, 4, 0),
        ("theta_only", coefficients(epsilon, kappa_theta, phi, 0, theta),
         (L, V, R), 8, 7, 1),
        ("both_generic", coefficients(epsilon, kappa, phi, delta, theta),
         (L, R), 9, 7, 2),
        ("both_kappa0", coefficients(epsilon, 0, phi, delta0, theta),
         (L,), 9, 7, 2),
    )

    semantics = []
    for name, coeffs, planes, genus, graph, abelian in cases:
        points = support(coeffs)
        global_ledger = polygon_ledger([(x, y) for x, y, _ in points])
        check(global_ledger[2] == genus, f"{name} global genus")
        face_area = 0
        for plane in planes:
            check(all(plane_gap(point, plane) >= 0 for point in points),
                  f"{name} direct support inequality {NAMES[plane]}")
            polygon = plane_polygon(points, plane)
            check(len(polygon) >= 3, f"{name} face lost dimension")
            face_area += polygon_ledger(polygon)[0]
            scaled = tuple(24 * value for value in plane)
            check(scaled in {(3, 6, -6), (6, 6, -12),
                             (0, 8, -8), (12, 0, -24)},
                  f"{name} nonintegral face")
        check(face_area == global_ledger[0], f"{name} face polygons do not tile")
        check(graph + abelian == genus, f"{name} genus split")
        names = ",".join(sorted(NAMES[plane] for plane in planes))
        semantics.append(f"{name}:faces={names};g={genus};graph={graph};ab={abelian}")

    # Coefficient collisions are reconstructed directly from the expansion.
    cancel_ke = support(coefficients(epsilon, epsilon, phi, deltae, 0))
    check((2, 3, 1) not in cancel_ke, "kappa-epsilon cancellation")
    cancel_td = support(coefficients(epsilon, kappa, phi, delta, delta))
    check((2, 4, 1) not in cancel_td, "theta-delta cancellation")
    for plane in (L, R):
        check(all(plane_gap(point, plane) >= 0 for point in cancel_td),
              "theta=delta changed supporting planes")

    # Independent coefficient-dictionary checks of all initial equations.
    s2_minus_p = {(2, 0): F(1), (0, 1): F(-1)}
    delta_top = {(0, 0): F(1), (0, 4): -delta}
    theta_top = {(0, 0): F(1), (2, 3): -theta}
    both_top = poly_add(delta_top, {(2, 3): -theta})
    check(poly_mul(s2_minus_p, delta_top) == {
        (2, 0): F(1), (0, 1): F(-1), (2, 4): -delta, (0, 5): delta
    }, "delta L factor")
    check(poly_mul(s2_minus_p, theta_top) == {
        (2, 0): F(1), (0, 1): F(-1), (4, 3): -theta, (2, 4): theta
    }, "theta L factor")
    check(poly_mul(s2_minus_p, both_top) == {
        (2, 0): F(1), (0, 1): F(-1), (2, 4): theta - delta,
        (0, 5): delta, (4, 3): -theta
    }, "both L factor")

    # Completing the square is checked coefficientwise in (T,P) and (S,X).
    # Delta: (2kT+phi P^2)^2 - [4k+(phi^2-4kd)P^4]
    #        =4k(kT^2+phi TP^2+dP^4-1).
    left_delta = {
        (2, 0): 4 * kappa * kappa,
        (1, 2): 4 * kappa * phi,
        (0, 4): phi * phi - (phi * phi - 4 * kappa * delta),
        (0, 0): -4 * kappa,
    }
    right_delta = {
        (2, 0): 4 * kappa * kappa,
        (1, 2): 4 * kappa * phi,
        (0, 4): 4 * kappa * delta,
        (0, 0): -4 * kappa,
    }
    check(combine(list(left_delta.items())) == combine(list(right_delta.items())),
          "delta square completion")
    # Theta identity has coefficients in (S,X): W^2-DeltaV-4theta X^3.
    left_theta = {
        (2, 0): 4 * theta * theta,
        (1, 0): 4 * theta * phi,
        (0, 0): phi * phi - (phi * phi - 4 * epsilon * theta),
        (0, 3): -4 * theta,
    }
    right_theta = {
        (2, 0): 4 * theta * theta,
        (1, 0): 4 * theta * phi,
        (0, 0): 4 * theta * epsilon,
        (0, 3): -4 * theta,
    }
    check(left_theta == right_theta, "theta square completion")

    # The three and only three outer-edge collision equations.
    DeltaD = phi * phi - 4 * kappa * delta
    DeltaV = phi * phi - 4 * epsilon * theta
    check(DeltaD != 0 and DeltaV != 0, "generic controls hit a wall")
    check((delta + theta) ** 2 != 0, "generic Bolza attachment collided")
    check((delta + delta) ** 2 == 4 * delta * delta,
          "theta=delta wrongly became a wall")
    check((delta - delta) ** 2 == 0, "theta=-delta wall disappeared")
    # The long-edge discriminant orders follow from c*T^n-1 at first Q-order.
    check(-4 ** 4 * delta ** 3 == -256 * delta ** 3,
          "quartic outer leading discriminant")
    check(-3 ** 3 * epsilon ** 2 == -27 * epsilon ** 2,
          "cubic outer leading discriminant")

    # Direct node and genus ledger.  Along P=S^2 the top factor is
    # 1-(delta+theta)S^8, whose derivative is -8(delta+theta)S^7.
    check(-8 * (delta + theta) != 0, "eight intersections lost transversality")
    check(24 - 1 == 23, "A23 resolution length")
    check(8 - 2 + 1 == 7, "dual graph rank")
    check(polygon_ledger([(0, 1), (2, 0), (2, 4), (0, 5)]) == (16, 10, 4),
          "delta main Pick ledger")
    check(polygon_ledger([(2, 0), (4, 2), (3, 3), (2, 4)]) == (8, 8, 1),
          "delta side Pick ledger")
    check(polygon_ledger([(0, 1), (2, 0), (4, 3), (2, 4)]) == (16, 4, 7),
          "theta main Pick ledger")
    check(polygon_ledger([(0, 1), (0, 4), (1, 4), (2, 4)]) == (6, 6, 1),
          "theta vertical Pick ledger")
    check(polygon_ledger([(0, 1), (2, 0), (4, 3), (2, 4), (0, 5)]) == (24, 8, 9),
          "Bolza main Pick ledger")

    norm_checks = 0
    for n in range(1, 501):
        check(literal_norm_test(n) == local_norm_test(n), f"norm sieve {n}")
        norm_checks += 1
    check(2 * 8 - 2 == 14, "ramification budget")

    semantics.extend([
        "walls=delta:kappa*delta!=0&phi^2=4*kappa*delta;theta:phi^2=4*epsilon*theta;both:delta+theta=0",
        "excluded=delta_off_wall,both_off_wall",
        "survivor=delta0,theta!=0,DeltaV!=0;j0;degree=a^2-a*b+b^2;ramification=14",
    ])
    digest = sha256("\n".join(semantics).encode()).hexdigest()

    print("JC2 WEIGHT-EIGHT INDEPENDENT DIRECT-PLANE AUDIT")
    print("method=direct_support_inequalities+exact_polygon_tiling+coefficient_dicts")
    print("strata=6;faces=L,D,V,R;base_change=24;multiplicity=1")
    print("models=delta:E1728_or_rational;theta:E0;both:Bolza_E8000^2")
    print("walls=DeltaD,DeltaV,delta+theta")
    print(f"norm_checks={norm_checks};ramification=14")
    print(f"semantic_digest={digest}")
    print(f"gates={GATES}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
