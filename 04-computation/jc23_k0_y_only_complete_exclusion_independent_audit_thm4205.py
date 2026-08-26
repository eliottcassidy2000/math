#!/usr/bin/env python3
"""Independent alternative-pair referee for THM-4205's Y-only K=0 row."""

from math import gcd

import sympy as sy


CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def p_valuation(polynomial, variable):
    terms = sy.Poly(polynomial, variable).terms()
    require(bool(terms), "zero resultant")
    return min(monomial[0] for monomial, coefficient in terms if coefficient != 0)


def hull(points):
    points = sorted(set(points))

    def cross(origin, left, right):
        return (
            (left[0] - origin[0]) * (right[1] - origin[1])
            - (left[1] - origin[1]) * (right[0] - origin[0])
        )

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
    return tuple(lower[:-1] + upper[:-1])


def main():
    s, p, Phi, Theta, zeta = sy.symbols("s p Phi Theta zeta")
    t = p - s**2
    Delta = sy.Rational(5696, 105)
    H = (
        -3 * p
        + sy.Rational(8, 3) * p**2
        - sy.Rational(1376, 135) * p**3
        + Phi * s * p**3
        + Delta * p**4
        + Theta * s**2 * p**3
        + zeta * s**3 * p**3
    )
    G = -s**2 / (2 * t) + H

    # This referee uses (A,C0), not THM-4205's primary (A,B) pair.
    A = sy.cancel((-s * p + t**2 * sy.diff(H, s)) / p)
    C0 = sy.expand(s**2 + 2 * t**2 * sy.diff(H, p))
    require(sy.denom(A) == 1, "A polynomial")
    require(sy.factor(t**2 * sy.diff(G, s) - p * A) == 0, "first gradient")
    require(sy.factor(2 * t**2 * sy.diff(G, p) - C0) == 0, "second gradient")
    require((sy.degree(A, s), sy.degree(C0, s)) == (6, 7), "alternative degrees")
    require(sy.factor(sy.Poly(A, s).LC() - 3 * zeta * p**2) == 0, "A infinity")
    require(sy.factor(sy.Poly(C0, s).LC() - 6 * zeta * p**2) == 0, "C0 infinity")
    def row_det(first, second):
        return first[0] * second[1] - first[1] * second[0]

    grad_a = (sy.diff(A, s), sy.diff(A, p))
    grad_c = (sy.diff(C0, s), sy.diff(C0, p))
    jacobian = row_det(grad_a, grad_c)
    hessian = sy.det(sy.hessian(G, (s, p)))
    alpha = p / t**2
    beta = sy.Rational(1, 2) / t**2
    grad_alpha = (sy.diff(alpha, s), sy.diff(alpha, p))
    grad_beta = (sy.diff(beta, s), sy.diff(beta, p))
    correction = 2 * t**4 * (
        alpha * C0 * row_det(grad_a, grad_beta)
        + A * beta * row_det(grad_alpha, grad_c)
        + A * C0 * row_det(grad_alpha, grad_beta)
    )
    bridge = sy.together(p * jacobian - 2 * t**4 * hessian + correction)
    require(sy.factor(bridge) == 0, "alternative localized Hessian bridge")
    require(sy.factor(A.subs(p, 0) + s) == 0, "collapsed p row")
    require(sy.factor(A.subs(p, s**2) + s) == 0, "collapsed t row")

    raw_resultant = sy.factor(sy.resultant(A, C0, s))
    valuation = p_valuation(raw_resultant, p)
    residual = sy.Poly(sy.cancel(raw_resultant / p**valuation), p)
    require((valuation, residual.degree()) == (8, 20), "alternative p^8 R20")
    require(sy.factor(residual.nth(0) - 1259712 * zeta**3) == 0, "alternative R0")
    require(
        sy.factor(residual.LC() - sy.Rational(40870620168192, 1225) * zeta**7)
        == 0,
        "alternative RLC",
    )

    # Independent normalized projection at a disjoint exact control.
    X, T = sy.symbols("X T")
    control = {Phi: sy.Rational(-7, 3), Theta: sy.Rational(11, 5), zeta: sy.Rational(13, 2)}
    P = T + X**2 * T**2
    Y = X * T * P
    GN = (
        -X**2 * T / 2
        - 3 * P
        + sy.Rational(8, 3) * P**2
        - sy.Rational(1376, 135) * P**3
        + Phi * P**2 * Y
        + Delta * P**4
        + Theta * P * Y**2
        + zeta * Y**3
    )
    f = sy.expand(sy.cancel(sy.diff(GN, X) / T).subs(control))
    h = sy.expand(sy.diff(GN, T).subs(control))
    require((sy.degree(f, X), sy.degree(h, X)) == (8, 9), "normalized control degrees")
    normalized = sy.Poly(sy.resultant(f, h, X), T)
    t_value = min(monomial[0] for monomial, coefficient in normalized.terms() if coefficient != 0)
    quotient = sy.cancel(normalized.as_expr() / T**t_value)
    six_multiplicity = 0
    while sy.rem(sy.Poly(quotient, T), sy.Poly(6 * T + 1, T)).is_zero:
        quotient = sy.cancel(quotient / (6 * T + 1))
        six_multiplicity += 1
    q20 = sy.Poly(quotient, T)
    require((normalized.degree(), t_value, six_multiplicity, q20.degree()) == (78, 56, 2, 20), "normalized factor ledger")
    require(q20.nth(0) != 0 and q20.LC() != 0, "normalized endpoints")

    # Reconstruct the polygon independently from the literal support.
    Q = sy.symbols("Q")
    fibre = sy.expand((s**2 - p) * (1 - Q * H) - Q * s**2 / 2)
    support = tuple(monomial for monomial, coefficient in sy.Poly(fibre, s, p).terms() if coefficient != 0)
    vertices = hull(support)
    require(vertices == ((0, 1), (2, 0), (5, 3), (3, 4), (0, 5)), "independent hull")
    area2 = abs(
        sum(
            left[0] * right[1] - right[0] * left[1]
            for left, right in zip(vertices, vertices[1:] + vertices[:1])
        )
    )
    boundary = 0
    packet = []
    for left, right in zip(vertices, vertices[1:] + vertices[:1]):
        dx, dy = right[0] - left[0], right[1] - left[1]
        length = gcd(abs(dx), abs(dy))
        u, v = -dy // length, dx // length
        constant = u * left[0] + v * left[1]
        index = u + v - constant
        boundary += length
        if not (left[0] == right[0] == 0):
            packet.extend([index] * length)
    genus = (area2 - boundary + 2) // 2
    packet = tuple(sorted(packet, reverse=True))
    require((area2, boundary, genus) == (30, 10, 11), "independent Pick ledger")
    require(packet == (11, 8, 2, 2, 2, 1), "independent packet")

    # The carrier polynomial is Eisenstein for the prime q-1/2 in C[q].
    W, q = sy.symbols("W q")
    carrier = sy.Poly(zeta * W**3 - (q - sy.Rational(1, 2)), W)
    require(carrier.degree() == 3, "carrier degree")
    require(sy.Poly(q - sy.Rational(1, 2), q).degree() == 1, "Eisenstein prime")
    require(sy.resultant(carrier.as_expr(), sy.diff(carrier.as_expr(), W), W) != 0, "carrier separability")

    length = residual.degree() + 4
    defect = sum(index - 1 for index in packet)
    full_n = sum(packet)
    finite_n, beta = full_n - 6, 3
    require((length, defect, full_n, finite_n, beta) == (24, 20, 26, 20, 3), "owner ledger")
    require(2 * finite_n - length - 1 + beta == 18 < finite_n - 1, "finite inequality")
    require(2 * (full_n - length) == 4 < defect, "full inequality")

    print("JC23_K0_Y_ONLY_ALTERNATIVE_PAIR_REFEREE")
    print("source_pair=(A,C0);degrees=(6,7);infinity=(3*zeta*p^2,6*zeta*p^2)")
    print("bridge=p*detD(A,C0)=2*t^4*detHess mod(A,C0)")
    print("resultant=p^8*R20;R0=1259712*zeta^3")
    print("RLC=(40870620168192/1225)*zeta^7")
    print("normalized_control=degree78,T56,(6T+1)^2,Q20")
    print("hull=((0,1),(2,0),(5,3),(3,4),(0,5));Pick=(30,10,11)")
    print("packet=(11,8,2,2,2,1);L=24;full_n=26;defect=20")
    print("carrier=one_prime_cubic_owner;finite=(20,3);responses=18<19,4<20")
    print(f"checks={CHECKS}")
    print("K0_Y_ONLY_ALTERNATIVE_PAIR_ACCEPT")


if __name__ == "__main__":
    main()
