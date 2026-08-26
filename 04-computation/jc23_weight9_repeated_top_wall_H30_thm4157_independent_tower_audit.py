#!/usr/bin/env python3
"""Clean-room exact audit of the H30 critical degeneration.

This file deliberately imports no producer module and reads no producer
artifact.  It rebuilds row D from the displayed source, derives H30 from the
primary projection discriminant, and works in

    K = Q[u]/h15(u),       eta^2 = u,

after the symmetry-adapted substitution X = eta*r.
"""

from hashlib import sha256
from math import isqrt
import os
import sys

import sympy as sp
from sympy.polys.domains import QQ


CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def progress(stage):
    if os.environ.get("H30_AUDIT_PROGRESS"):
        print("progress=" + stage, file=sys.stderr, flush=True)


def monic_digest(poly):
    p = poly.monic()
    payload = "degree=" + str(p.degree()) + ";" + ",".join(
        f"{sp.Rational(c).p}/{sp.Rational(c).q}" for c in p.all_coeffs()
    )
    return sha256(payload.encode()).hexdigest()


def term_digest(poly):
    payload = ";".join(
        ",".join(str(i) for i in monomial) + ":"
        + str(sp.Rational(value).p) + "/" + str(sp.Rational(value).q)
        for monomial, value in sorted(poly.terms())
    )
    return sha256(payload.encode()).hexdigest()


def eta_odd_to_u(expression, eta, u):
    """Return expression/eta with eta^2 replaced by u."""
    p = sp.Poly(sp.expand(expression), eta)
    require(all(e[0] % 2 == 1 for e, _ in p.terms()),
            "primary residual is not eta-odd")
    return sp.expand(sum(c * u ** ((e[0] - 1) // 2)
                         for e, c in p.terms()))


def tower_reduce(expression, X, T, eta, r, u, parity):
    """Substitute X=eta*r, divide eta^parity, replace eta^2 by u."""
    p = sp.Poly(sp.expand(expression), X, T, eta, domain=QQ)
    output = 0
    observed = set()
    for (xd, td, ed), coefficient in p.terms():
        observed.add((xd + ed) % 2)
        exponent = xd + ed - parity
        require(exponent >= 0 and exponent % 2 == 0,
                "invalid tower exponent")
        output += coefficient * r ** xd * T ** td * u ** (exponent // 2)
    require(observed == {parity}, "combined (X,eta)-parity changed")
    return sp.expand(output)


def field_poly_in_T(expression, T, u, field, u0):
    p = sp.Poly(expression, T, u, domain=QQ)
    coefficients = {}
    up = [field.one]
    for _ in range(p.degree(u)):
        up.append(up[-1] * u0)
    for (td, ud), coefficient in p.terms():
        coefficients[td] = coefficients.get(td, field.zero) \
            + field.convert(coefficient) * up[ud]
    return sp.Poly.from_dict({(d,): c for d, c in coefficients.items()},
                             gens=T, domain=field)


def specialize_rT(expression, r, T, u, t0, field, u0):
    """Specialize u=u0,T=t0, retaining a polynomial in r."""
    p = sp.Poly(expression, r, T, u, domain=QQ)
    tp = [field.one]
    up = [field.one]
    for _ in range(p.degree(T)):
        tp.append(tp[-1] * t0)
    for _ in range(p.degree(u)):
        up.append(up[-1] * u0)
    coefficients = {}
    for (rd, td, ud), coefficient in p.terms():
        coefficients[rd] = coefficients.get(rd, field.zero) \
            + field.convert(coefficient) * tp[td] * up[ud]
    return sp.Poly.from_dict({(d,): c for d, c in coefficients.items()},
                             gens=r, domain=field)


def evaluate_rT(expression, r, T, u, r0, t0, field, u0):
    p = sp.Poly(expression, r, T, u, domain=QQ)
    rp = [field.one]
    tp = [field.one]
    up = [field.one]
    for _ in range(p.degree(r)):
        rp.append(rp[-1] * r0)
    for _ in range(p.degree(T)):
        tp.append(tp[-1] * t0)
    for _ in range(p.degree(u)):
        up.append(up[-1] * u0)
    answer = field.zero
    for (rd, td, ud), coefficient in p.terms():
        answer += (field.convert(coefficient) * rp[rd] * tp[td] * up[ud])
    return answer


def field_payload(element):
    return ",".join(str(c) for c in element.to_list())


def main():
    X, T, eta, u, r = sp.symbols("X T eta u r")

    # Exact row-D source, reconstructed directly rather than imported.
    P = T + X ** 2 * T ** 2
    Y = X * T * P
    epsilon = -sp.Rational(1376, 135)
    K = sp.Rational(1376, 135)
    Delta = sp.Rational(2048, 45)
    Theta = -Delta
    G = sp.expand(
        -X ** 2 * T / 2 - 3 * P + sp.Rational(8, 3) * P ** 2
        + epsilon * P ** 3 + K * Y ** 2 + Delta * P ** 4
        + Theta * P * Y ** 2 + eta * P ** 3 * Y - eta * Y ** 3
    )
    # Independently expose the cancellations special to row D.
    G_contracted = sp.expand(
        -X ** 2 * T / 2 - 3 * P + sp.Rational(8, 3) * P ** 2
        + epsilon * T * P ** 2 + Delta * T * P ** 3
        + eta * X * T ** 2 * P ** 3
    )
    require(G == G_contracted, "row-D contraction failed")

    f = sp.cancel(sp.diff(G, X) / T)
    h = sp.diff(G, T)
    require(sp.rem(sp.diff(G, X), T, T) == 0, "G_X is not divisible by T")
    require(sp.degree(f, X) == 6 and sp.degree(h, X) == 7,
            "source derivative X-degrees changed")
    require(sp.factor(sp.Poly(f, X).LC() - 7 * T ** 7 * eta) == 0,
            "f leading-X coefficient changed")
    require(sp.factor(sp.Poly(h, X).LC() - 8 * T ** 7 * eta) == 0,
            "h leading-X coefficient changed")

    jac = sp.expand(sp.diff(f, X) * sp.diff(h, T)
                    - sp.diff(f, T) * sp.diff(h, X))
    hessian = sp.expand(sp.diff(G, X, 2) * sp.diff(G, T, 2)
                        - sp.diff(G, X, T) ** 2)
    require(sp.cancel(T * jac - hessian - f * sp.diff(G, X, T)) == 0,
            "Hessian/Jacobian identity failed")

    # Symmetry reduction used later at the actual point.
    F = tower_reduce(f, X, T, eta, r, u, 1)
    H = tower_reduce(h, X, T, eta, r, u, 0)
    Jac = tower_reduce(jac, X, T, eta, r, u, 0)
    Hess = tower_reduce(hessian, X, T, eta, r, u, 0)
    require(sp.cancel(Jac - (sp.diff(F, r) * sp.diff(H, T)
                             - sp.diff(F, T) * sp.diff(H, r))) == 0,
            "reduced Jacobian chain rule failed")
    require(sp.cancel(T * Jac - Hess - u * F * tower_reduce(
        sp.diff(G, X, T), X, T, eta, r, u, 1)) == 0,
        "reduced Hessian/Jacobian identity failed")
    progress("source")

    # Rebuild the primary projection independently.
    resultant_X = sp.resultant(f, h, X)
    quotient = sp.cancel(resultant_X / (T ** 35 * (6 * T + 1) ** 2))
    require(sp.denom(quotient) == 1, "residual resultant is not polynomial")
    Q_eta = sp.Poly(sp.expand(quotient), T, eta, domain=QQ)
    require(Q_eta.degree(T) == 12, "residual resultant degree changed")
    Q_u_expr = eta_odd_to_u(Q_eta.as_expr(), eta, u)
    Q_u = sp.Poly(Q_u_expr, T, u, domain=QQ)
    require(Q_u.degree(T) == 12, "u-reduced residual degree changed")
    progress("resultant")

    # Derive h15, rather than reading the producer's factor ledger.
    discriminant = sp.Poly(sp.discriminant(Q_u_expr, T), u, domain=QQ)
    _, factors = sp.factor_list(discriminant)
    factor_degrees = sorted((factor.degree(), multiplicity)
                            for factor, multiplicity in factors)
    degree15 = [(factor, multiplicity) for factor, multiplicity in factors
                if factor.degree() == 15]
    require(len(degree15) == 1 and degree15[0][1] == 1,
            "unique simple degree-15 discriminant factor missing")
    h15 = degree15[0][0].primitive()[1]
    if h15.LC() < 0:
        h15 = -h15
    require(h15.is_irreducible, "h15 is reducible over Q")
    progress("h15")

    # Prove eta^2-u0 irreducible over K with a norm obstruction.  If u0 were
    # a square in K, Norm_K/Q(u0) would be a rational square.
    hm = h15.monic()
    norm_u = sp.Rational((-1) ** hm.degree()) * sp.Rational(hm.nth(0))
    norm_num = abs(int(norm_u.p))
    norm_den = int(norm_u.q)
    norm_is_square = (norm_u >= 0 and isqrt(norm_num) ** 2 == norm_num
                      and isqrt(norm_den) ** 2 == norm_den)
    require(not norm_is_square, "Norm(u0) is a rational square")
    require(sp.Poly(hm.as_expr().subs(u, eta ** 2), eta,
                    domain=QQ).is_irreducible,
            "h15(eta^2) is reducible over Q")
    progress("tower_irreducibility")
    if os.environ.get("H30_AUDIT_STAGE1_ONLY"):
        print("stage1_h15_monic_sha256=" + monic_digest(h15))
        print("stage1_Q_u_sha256=" + term_digest(Q_u))
        print("stage1_factor_degrees=" + str(factor_degrees))
        print("stage1_tower_irreducible=yes")
        print("stage1_checks=" + str(CHECKS))
        return

    field = QQ.alg_field_from_poly(hm, alias="u0")
    u0 = field.unit
    QK = field_poly_in_T(Q_u_expr, T, u, field, u0)
    require(QK.degree() == 12 and QK.LC() != field.zero,
            "Q loses degree at H30")
    q_subresultants = sp.subresultants(QK, QK.diff(), T)
    q_subresultant_degrees = [item.degree() for item in q_subresultants]
    q_common = sp.gcd(QK, QK.diff()).monic()
    require(q_common.degree() == 1,
            "specialized Q does not have exactly one double root")
    require(q_subresultants[-1].degree() == 1,
            "specialized Q subresultant does not terminate linearly")
    t0 = field.convert(-q_common.nth(0))
    require(t0 != field.zero and 6 * t0 + 1 != field.zero,
            "repeated T is universal")
    require(QK.eval(field.zero) != field.zero
            and QK.eval(field.convert(-sp.Rational(1, 6))) != field.zero,
            "H30 meets a residual endpoint wall")
    require(sp.gcd(h15, sp.Poly(u * (22143375 * u + 15510536192), u,
                                  domain=QQ)).degree() == 0,
            "H30 meets eta=0 or J=0")
    progress("specialized_Q")

    # At the repeated T-root, compute the common finite r-root directly.
    Ft = specialize_rT(F, r, T, u, t0, field, u0)
    Ht = specialize_rT(H, r, T, u, t0, field, u0)
    require(Ft.degree() == 6 and Ht.degree() == 7,
            "an X-leading coefficient vanishes at t0")
    require(Ft.LC() != field.zero and Ht.LC() != field.zero,
            "finite-X infinity sidecar is not a unit")
    point_subresultants = sp.subresultants(Ft, Ht, r)
    point_degrees = [item.degree() for item in point_subresultants]
    common_r = sp.gcd(Ft, Ht).monic()
    require(common_r.degree() == 1,
            "repeated T fibre does not have one finite point")
    require(point_subresultants[-1].degree() == 1,
            "point subresultant does not terminate linearly")
    r0 = field.convert(-common_r.nth(0))
    require(Ft.eval(r0) == field.zero and Ht.eval(r0) == field.zero,
            "computed point is not critical")
    progress("point")

    # Direct independent degeneration checks.
    jac0 = evaluate_rT(Jac, r, T, u, r0, t0, field, u0)
    hess0 = evaluate_rT(Hess, r, T, u, r0, t0, field, u0)
    require(jac0 == field.zero, "Jacobian of (f,h) is nonzero")
    require(hess0 == field.zero, "Hessian of G is nonzero")
    progress("hessian")

    h_digest = monic_digest(h15)
    q_digest = term_digest(Q_u)
    t_digest = sha256(field_payload(t0).encode()).hexdigest()
    r_digest = sha256(field_payload(r0).encode()).hexdigest()
    semantic = "\n".join((
        "h15=" + ",".join(str(c) for c in h15.all_coeffs()),
        "Q_u=" + q_digest,
        "q_subresultants=" + ",".join(map(str, q_subresultant_degrees)),
        "point_subresultants=" + ",".join(map(str, point_degrees)),
        "t0=" + field_payload(t0),
        "r0=" + field_payload(r0),
        "jacobian=hessian=0",
    )) + "\n"

    print("CLEAN-ROOM H30 TOWER/NONREDUCEDNESS AUDIT")
    print("source=direct-row-D;producer_imports=none;producer_reads=none")
    print("primary_factor_degrees=" + str(factor_degrees))
    print("h15_degree=15;h15_irreducible=yes")
    print("Norm(u0)_square=no;eta^2-u0_irreducible=yes")
    print("H30=h15(eta^2);H30_irreducible=yes")
    print("Q_degree=12;gcd(Q,Q')=linear;double_root=unique")
    print("Q_subresultant_degrees=" + str(q_subresultant_degrees))
    print("repeated_T_avoids=0,-1/6;Q_infinity_unit=yes")
    print("point_subresultant_degrees=" + str(point_degrees))
    print("unique_finite_point=(eta*r0,t0);X_infinity_unit=yes")
    print("intersection_multiplicity=2;critical_scheme_nonreduced=yes")
    print("Jacobian(f,h)=0;Hessian(G)=0;identity=verified")
    print("gcd(H30,eta*J)=1")
    print("h15_monic_sha256=" + h_digest)
    print("Q_u_sha256=" + q_digest)
    print("t0_sha256=" + t_digest)
    print("r0_sha256=" + r_digest)
    print("semantic_sha256=" + sha256(semantic.encode()).hexdigest())
    print("checks=" + str(CHECKS))
    print("verdict=PASS")


if __name__ == "__main__":
    main()
