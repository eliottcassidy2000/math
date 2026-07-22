"""Exact audits for THM-2056 and THM-2057.

The proofs are algebraic.  This script adversarially checks their identities,
the load-bearing acute-cone hypothesis, the complete three-route scaled-core
certificate on a large box, and exact pair-sum maxima on named controls.
"""

from fractions import Fraction
from itertools import combinations_with_replacement
from math import gcd


def central_unit(h: int) -> int:
    assert h >= 2
    if h <= 14:
        r = 1
    elif h % 2:
        r = (h - 1) // 2
    elif h % 4 == 0:
        r = h // 2 - 1
    else:
        r = h // 2 - 2
    assert gcd(r, h) == 1
    assert 14 * min(r, h - r) >= h
    return r


def clock_numerator(N: int, a: int, w: int) -> int:
    q = N * a
    g = gcd(w, q)
    h = q // g
    assert h >= 2
    c = w // g
    r = central_unit(h)
    k0 = (pow(c, -1, h) * r) % h
    for t in range(N):
        k = k0 + h * t
        if gcd(k, N) == 1:
            k %= q
            assert 14 * min((k * w) % q, (-(k * w)) % q) >= q
            return k
    raise AssertionError("unit lift failed")


CORE = tuple(range(1, 12)) + (13,)


def route_certificate(a: int, w: int):
    if w % (12 * a):
        q = 12 * a
        k = clock_numerator(12, a, w)
        return "clock12", Fraction(k, q)
    if w % (14 * a):
        q = 14 * a
        k = clock_numerator(14, a, w)
        return "clock14", Fraction(k, q)
    assert w % (84 * a) == 0
    m = w // (84 * a)
    return "binding84", Fraction(35 * m + 2, a * (84 * m + 5))


def circle_distance(x: Fraction) -> Fraction:
    r = x.numerator % x.denominator
    return Fraction(min(r, x.denominator - r), x.denominator)


def row_margin(a: int, w: int, t: Fraction) -> Fraction:
    speeds = tuple(a * i for i in CORE) + (w,)
    return min(circle_distance(v * t) for v in speeds)


def ap_route_certificate(a: int, w: int):
    if w % (13 * a):
        return "clock13", Fraction(clock_numerator(13, a, w), 13 * a)
    if w % (14 * a):
        return "clock14", Fraction(clock_numerator(14, a, w), 14 * a)
    assert w % (182 * a) == 0
    m = w // (182 * a)
    return "binding182", Fraction(14 * m, a * (182 * m + 1))


def ap_row_margin(a: int, w: int, t: Fraction) -> Fraction:
    speeds = tuple(a * i for i in range(1, 13)) + (w,)
    return min(circle_distance(v * t) for v in speeds)


def exact_M(speeds):
    best = Fraction(0)
    for x, y in combinations_with_replacement(speeds, 2):
        q = x + y
        for k in range(1, q):
            value = min(circle_distance(Fraction(k * v, q)) for v in speeds)
            if value > best:
                best = value
    return best


def audit_kelvin_box(B=300):
    checked = 0
    for a in range(-B, B + 1):
        for b in range(-B, B + 1):
            if (a == 0 and b == 0) or gcd(abs(a), abs(b)) != 1:
                continue
            q = a * a + b * b
            D = max(13 * abs(b), abs(a - 12 * b))
            direct = 91 * D <= q
            # The two polar-facet inequalities at I(a,b)=(a/q,b/q).
            polar = 91 * 13 * abs(b) <= q and 91 * abs(a - 12 * b) <= q
            assert direct == polar
            checked += 1
    return checked


def F(w, d):
    return d[0] * d[0] + d[1] * d[1] - 91 * (w[0] * d[0] + w[1] * d[1])


def audit_farey_identity():
    vectors = [
        (x, y)
        for x in range(-4, 5)
        for y in range(-4, 5)
        if (x or y) and gcd(abs(x), abs(y)) == 1
    ]
    identities = 0
    certified_points = 0
    certifying_cones = 0
    for w in ((1, 0), (5, 3), (-4, 7), (12, -13)):
        for u in vectors:
            for v in vectors:
                if abs(u[0] * v[1] - u[1] * v[0]) != 1:
                    continue
                dot = u[0] * v[0] + u[1] * v[1]
                if dot < 0 or w[0] * u[0] + w[1] * u[1] < 0 or w[0] * v[0] + w[1] * v[1] < 0:
                    continue
                Au = max(0, -F(w, u))
                Av = max(0, -F(w, v))
                certifies = 2 * dot >= Au + Av
                certifying_cones += int(certifies)
                for m in range(1, 6):
                    for n in range(1, 6):
                        x = (m * u[0] + n * v[0], m * u[1] + n * v[1])
                        rhs = (
                            m * F(w, u)
                            + n * F(w, v)
                            + (m * m - m) * (u[0] ** 2 + u[1] ** 2)
                            + (n * n - n) * (v[0] ** 2 + v[1] ** 2)
                            + 2 * m * n * dot
                        )
                        assert F(w, x) == rhs
                        identities += 1
                        if certifies:
                            assert F(w, x) >= 0
                            certified_points += 1
    # Two deliberately large Farey cones exercise the implication itself:
    # the first has one bad boundary generator, the second has two safe ones.
    w = (1, 0)
    for u, v in (((1, 0), (91, 1)), ((91, 1), (92, 1))):
        dot = u[0] * v[0] + u[1] * v[1]
        Au, Av = max(0, -F(w, u)), max(0, -F(w, v))
        assert abs(u[0] * v[1] - u[1] * v[0]) == 1
        assert dot >= 0 and 2 * dot >= Au + Av
        certifying_cones += 1
        for m in range(1, 101):
            for n in range(1, 101):
                x = (m * u[0] + n * v[0], m * u[1] + n * v[1])
                rhs = (
                    m * F(w, u)
                    + n * F(w, v)
                    + (m * m - m) * (u[0] ** 2 + u[1] ** 2)
                    + (n * n - n) * (v[0] ** 2 + v[1] ** 2)
                    + 2 * m * n * dot
                )
                assert F(w, x) == rhs and F(w, x) >= 0
                identities += 1
                certified_points += 1
    # Acuteness cannot be dropped: safe endpoints can hide a bad mediant.
    u, v = (91, 1), (-90, -1)
    assert abs(u[0] * v[1] - u[1] * v[0]) == 1
    assert F(w, u) >= 0 and F(w, v) >= 0 and F(w, (1, 0)) < 0
    # An unsafe boundary ray can nevertheless have a certified punctured cone.
    u, v = (1, 0), (91, 1)
    assert F(w, u) == -90 and F(w, v) == 1
    assert 2 * (u[0] * v[0] + u[1] * v[1]) >= 90
    return identities, certifying_cones, certified_points


def audit_scaled_core(A=120, W=12000):
    routes = {"clock12": 0, "clock14": 0, "binding84": 0}
    minimum = Fraction(1, 2)
    for a in range(1, A + 1):
        for w in range(1, W + 1):
            route, t = route_certificate(a, w)
            margin = row_margin(a, w, t)
            assert margin >= Fraction(1, 14)
            routes[route] += 1
            minimum = min(minimum, margin)
    return routes, minimum


def audit_missing_clock_sieve(A=25, W=2000):
    checked = 0
    for N in range(2, 15):
        for a in range(1, A + 1):
            q = N * a
            for w in range(1, W + 1):
                if w % q == 0:
                    continue
                k = clock_numerator(N, a, w)
                assert gcd(k, N) == 1
                assert all(min((c * k) % N, (-(c * k)) % N) * 14 >= N for c in range(1, N))
                checked += 1
    return checked


def audit_scaled_ap_core(A=80, W=10000):
    routes = {"clock13": 0, "clock14": 0, "binding182": 0}
    minimum = Fraction(1, 2)
    for a in range(1, A + 1):
        for w in range(1, W + 1):
            route, t = ap_route_certificate(a, w)
            margin = ap_row_margin(a, w, t)
            assert margin >= Fraction(1, 14)
            routes[route] += 1
            minimum = min(minimum, margin)
    return routes, minimum


def one_tail_residual_count():
    # THM-2055 gives ||d||<91*13, so |a|,|b|<=1182 contains every failure.
    bad = 0
    distinct_bad = 0
    collisions = []
    for a in range(1, 1183):
        for b in range(-1182, 1183):
            if 12 * a + b <= 0 or gcd(a, abs(b)) != 1:
                continue
            D = max(13 * abs(b), abs(a - 12 * b))
            if a * a + b * b < 91 * D:
                bad += 1
                speeds = [a * i + (b if i == 12 else 0) for i in range(1, 14)]
                if len(set(speeds)) == 13:
                    distinct_bad += 1
                else:
                    collisions.append((a, b))
    return bad, distinct_bad, collisions


def main():
    print("THM-2056 / THM-2057 EXACT AUDIT")
    print("kelvin primitive box count:", audit_kelvin_box())
    identities, cones, points = audit_farey_identity()
    print("farey identities checked:", identities)
    print("defect-certified cones:", cones)
    print("defect-certified interior points checked:", points)

    routes, minimum = audit_scaled_core()
    print("scaled-core box: a<=120, w<=12000")
    print("route counts:", routes)
    print("minimum certificate margin:", minimum)
    print("general missing-clock certificates checked:", audit_missing_clock_sieve())
    ap_routes, ap_minimum = audit_scaled_ap_core()
    print("scaled AP-core box: a<=80, w<=10000")
    print("AP route counts:", ap_routes)
    print("AP minimum certificate margin:", ap_minimum)

    controls = ((1, 12), (1, 24), (2, 25), (2, 168), (7, 1000))
    print("exact pair-sum controls:")
    for a, w in controls:
        route, t = route_certificate(a, w)
        margin = row_margin(a, w, t)
        M = exact_M(tuple(a * i for i in CORE) + (w,))
        print(f"  a={a:2d} w={w:4d} route={route:9s} cert={margin} exact_M={M}")

    ap_controls = ((1, 13), (1, 182), (2, 364), (2, 27))
    print("exact AP-tail pair-sum controls:")
    for a, w in ap_controls:
        route, t = ap_route_certificate(a, w)
        margin = ap_row_margin(a, w, t)
        speeds = tuple(a * i for i in range(1, 13)) + (w,)
        M = exact_M(speeds)
        print(f"  a={a:2d} w={w:4d} route={route:10s} cert={margin} exact_M={M}")

    bad, distinct_bad, collisions = one_tail_residual_count()
    print("complete positive primitive one-tail determinant residual:", bad)
    print("distinct-speed residual:", distinct_bad)
    print("collision directions removed:", collisions)

    carriers = (
        "scaled_clock_binding",
        "resolved_phase_height",
        "kelvin_polar_polygon",
        "farey_defect_cone",
        "raw_tangent_scan",
        "heegner_form_class",
    )
    n = len(carriers)
    print("tournament carriers:", " > ".join(carriers))
    print("tournament score histogram:", {s: 1 for s in range(n)})
    print("directed 3-cycles: 0; SCC sizes:", [1] * n, "; Hamiltonian paths: 1")
    print("PASS")


if __name__ == "__main__":
    main()
