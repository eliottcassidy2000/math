#!/usr/bin/env python3
"""Exact H30 nonreduced-critical-point certificate via subresultants.

The repeated T-resultant root has exactly one finite X-point.  Hence its
intersection multiplicity is at least two, so the Jacobian of (f,h) vanishes.
On the critical scheme T*det D(f,h)=det Hess(G), and T is a unit.
"""

import ast
from hashlib import sha256
import importlib.util
import os
from pathlib import Path
import sys

import sympy as sp
from sympy.polys.domains import QQ


ROOT = Path(__file__).resolve().parents[1]
PRIMARY = ROOT / "04-computation/jc23_weight9_repeated_top_wall_thm4157.py"
FACTOR_OUTPUT = (
    ROOT / "05-knowledge/results/"
    "jc23_weight9_repeated_top_wall_factors_thm4157.out"
)
FACTOR_OUTPUT_SHA256 = (
    "91b91809c99ebeecf26af27fdb7168311fa5124a90f71a8f942c1495133cf483"
)

spec = importlib.util.spec_from_file_location("repeated_wall_primary", PRIMARY)
primary = importlib.util.module_from_spec(spec)
spec.loader.exec_module(primary)

CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def progress(stage):
    if os.environ.get("H30_PROGRESS"):
        print(f"progress={stage}", file=sys.stderr, flush=True)


def eta_parity_to_u(expression, eta, u, parity):
    polynomial = sp.Poly(expression, eta)
    require(all(exponent[0] % 2 == parity
                for exponent, _ in polynomial.terms()),
            "eta parity changed")
    return sp.expand(sum(
        coefficient * u ** ((exponent[0] - parity) // 2)
        for exponent, coefficient in polynomial.terms()
    ))


def polynomial_digest(polynomial):
    payload = ";".join(
        f"{','.join(str(index) for index in monomial)}:{sp.Rational(value).p}/"
        f"{sp.Rational(value).q}"
        for monomial, value in sorted(polynomial.terms())
    )
    return sha256(payload.encode()).hexdigest()


def main():
    X, T, eta, u = sp.symbols("X T eta u")
    P = T + X ** 2 * T ** 2
    Y = X * T * P
    G = primary.source_polynomial(
        X, T, P, Y,
        sp.Rational(2048, 45), -sp.Rational(2048, 45),
        sp.Rational(0), eta,
    )
    f_expr = sp.cancel(sp.diff(G, X) / T)
    h_expr = sp.diff(G, T)
    f = sp.Poly(f_expr, X, T, eta, domain=QQ)
    h = sp.Poly(h_expr, X, T, eta, domain=QQ)
    hessian_expr = (
        sp.diff(G, X, 2) * sp.diff(G, T, 2)
        - sp.diff(G, X, T) ** 2
    )
    require(sp.factor(sp.Poly(f_expr, X).LC() - 7 * T ** 7 * eta) == 0,
            "f leading-X sidecar changed")
    require(sp.factor(sp.Poly(h_expr, X).LC() - 8 * T ** 7 * eta) == 0,
            "h leading-X sidecar changed")
    jacobian_expr = (
        sp.diff(f_expr, X) * sp.diff(h_expr, T)
        - sp.diff(f_expr, T) * sp.diff(h_expr, X)
    )
    require(sp.cancel(
        T * jacobian_expr - hessian_expr - f_expr * sp.diff(G, X, T)
    ) == 0, "critical Hessian/Jacobian identity changed")

    x_subresultants = sp.subresultants(f_expr, h_expr, X)
    require([sp.degree(item, X) for item in x_subresultants]
            == [7, 6, 5, 4, 3, 2, 1, 0],
            "symbolic X-subresultant ladder changed")
    linear_raw = sp.Poly(x_subresultants[-2], X)
    linear_content = sp.gcd(linear_raw.nth(1), linear_raw.nth(0))
    require(sp.factor(
        linear_content - T ** 35 * eta * (6 * T + 1) / 6
    ) == 0, "linear X-subresultant content changed")
    linear_primitive = sp.Poly(
        sp.cancel(linear_raw.as_expr() / linear_content), X
    )
    A_u = eta_parity_to_u(linear_primitive.nth(1), eta, u, 0)
    B_u = eta_parity_to_u(linear_primitive.nth(0), eta, u, 1)
    A_poly = sp.Poly(A_u, T, u, domain=QQ)
    B_poly = sp.Poly(B_u, T, u, domain=QQ)
    require(A_poly.degree(T) == 12 and B_poly.degree(T) == 11,
            "primitive linear subresultant degrees changed")
    progress("source_and_linear_subresultant")

    resultant = sp.resultant(f_expr, h_expr, X)
    Q12 = sp.Poly(sp.cancel(
        resultant / (T ** 35 * (6 * T + 1) ** 2)
    ), T)
    require(Q12.degree() == 12, "primary residual degree changed")
    Q_u = eta_parity_to_u(Q12.as_expr(), eta, u, 1)
    Q_u_poly = sp.Poly(Q_u, T, u, domain=QQ)
    require(Q_u_poly.degree(T) == 12,
            "symmetry-reduced primary residual degree changed")
    progress("primary_resultant")

    factor_bytes = Path(FACTOR_OUTPUT).read_bytes()
    require(sha256(factor_bytes).hexdigest() == FACTOR_OUTPUT_SHA256,
            "deepest-factor output hash changed")
    factor_lines = factor_bytes.decode().splitlines()
    h30_index = factor_lines.index("factor=H30;degree_eta=30;degree_u=15")
    prefix = "  primitive_integer_coefficients="
    require(factor_lines[h30_index + 1].startswith(prefix),
            "H30 coefficient ledger is missing")
    h15_coefficients = ast.literal_eval(
        factor_lines[h30_index + 1][len(prefix):]
    )
    require(len(h15_coefficients) == 16,
            "h15 coefficient count changed")
    h15 = sp.Poly.from_list(h15_coefficients, gens=u, domain=QQ).monic()
    progress("h15_ledger")
    QA_resultant = sp.Poly(sp.resultant(Q_u, A_u, T), u, domain=QQ)
    require(sp.gcd(h15, QA_resultant).degree() == 0,
            "a residual Q-root loses the linear X coefficient on H30")
    progress("linear_coefficient_unit")

    field = QQ.alg_field_from_poly(h15, alias="u0")
    u0 = field.unit

    def u_element(expression):
        polynomial = sp.Poly(expression, u, domain=QQ)
        answer = field.zero
        for coefficient in polynomial.all_coeffs():
            answer = answer * u0 + field.convert(coefficient)
        return answer

    def univariate_over_field(expression):
        polynomial = sp.Poly(expression, T)
        return sp.Poly.from_list(
            [u_element(coefficient)
             for coefficient in polynomial.all_coeffs()],
            gens=T,
            domain=field,
        )

    Q_field = univariate_over_field(Q_u)
    require(Q_field.degree() == 12,
            "h15-specialized residual degree changed")
    common_T = sp.gcd(Q_field, Q_field.diff()).monic()
    progress("common_T")
    require(common_T.degree() == 1,
            "H30 does not give one repeated T value")
    Q_at_zero = sp.together(Q_u.subs(T, 0)).as_numer_denom()[0]
    Q_at_universal = sp.together(
        Q_u.subs(T, -sp.Rational(1, 6))
    ).as_numer_denom()[0]
    require(sp.gcd(h15, sp.Poly(Q_at_zero, u)).degree() == 0
            and sp.gcd(h15, sp.Poly(Q_at_universal, u)).degree() == 0,
            "repeated T value can meet a universal divisor")

    # The primitive penultimate subresultant is A(T,u)X+eta*B(T,u).
    # Coprimality of A with the common T-factor says it remains a genuine
    # nonzero linear polynomial at t0, hence the fibre has exactly one finite
    # common X-root.  Unit X-leading coefficients exclude an infinity point.
    # Since common_T divides Q and Q', the resultant has multiplicity >=2 at
    # t0.  With one finite point over t0, that point has local intersection
    # multiplicity >=2.  Thus det D(f,h)=0 there.  The checked symbolic
    # identity gives det Hess(G)=T*det D(f,h) on f=0, and T is a unit.
    J_u = 22143375 * u + 15510536192
    require(sp.gcd(h15, sp.Poly(u * J_u, u)).degree() == 0,
            "H30 meets eta=0 or J=0")

    q_digest = polynomial_digest(Q_u_poly)
    a_digest = polynomial_digest(A_poly)
    b_digest = polynomial_digest(B_poly)
    semantic = "\n".join((
        "h15=" + ",".join(str(c) for c in h15.all_coeffs()),
        "Q_u=" + q_digest,
        "A_u=" + a_digest,
        "B_u=" + b_digest,
        "gcdT=linear;gcd(h15,ResT(Q,A))=1;nonreduced=yes;hessian=0",
    )) + "\n"
    print("DEEPEST REPEATED-WALL H30 NONREDUCEDNESS CERTIFICATE")
    print("field=Q[u]/h15(u);tower_relation=eta^2-u")
    print("h15_degree=15;H30(eta)=h15(eta^2);irreducible=yes")
    print("gcd(Q12/eta,d(Q12/eta))_degree=1")
    print("repeated_T_avoids=T*(6T+1)")
    print("primitive_linear_X_subresultant=A(T,u)X+eta*B(T,u)")
    print("gcd(h15,Res_T(Q_u,A_u))=1;unique_finite_X_point=yes")
    print("resultant_multiplicity_at_point>=2;critical_scheme_nonreduced=yes")
    print("T*det_D(f,h)=det_Hess_G_on_f=0;det_Hess_G=0")
    print("gcd(H30,eta*J)=1")
    print(f"Q_u_sha256={q_digest}")
    print(f"A_u_sha256={a_digest}")
    print(f"B_u_sha256={b_digest}")
    print(f"semantic_sha256={sha256(semantic.encode()).hexdigest()}")
    print(f"checks={CHECKS}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
