#!/usr/bin/env python3
"""Extract exact row-D projection-discriminant factors in u=eta^2."""

from hashlib import sha256
import importlib.util
from pathlib import Path

import sympy as sp


PRIMARY = Path(__file__).with_name("jc23_weight9_repeated_top_wall_thm4157.py")
spec = importlib.util.spec_from_file_location("repeated_wall_primary", PRIMARY)
primary = importlib.util.module_from_spec(spec)
spec.loader.exec_module(primary)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def monic_digest(poly):
    monic = poly.monic()
    payload = "degree=" + str(monic.degree()) + ";" + ",".join(
        f"{sp.Rational(coefficient).p}/{sp.Rational(coefficient).q}"
        for coefficient in monic.all_coeffs()
    )
    return sha256(payload.encode()).hexdigest()


def even_factor_to_u(factor, eta, u):
    eta_poly = sp.Poly(factor, eta)
    require(all(exponent[0] % 2 == 0 for exponent in eta_poly.as_dict()),
            "factor is not even in eta")
    u_expression = sum(
        coefficient * u ** (exponent[0] // 2)
        for exponent, coefficient in eta_poly.terms()
    )
    primitive = sp.Poly(u_expression, u).primitive()[1]
    if primitive.LC() < 0:
        primitive = -primitive
    return primitive, primitive.monic()


def main():
    X, T, s, p, eta, u = sp.symbols("X T s p eta u")
    P = T + X ** 2 * T ** 2
    Y = X * T * P
    Delta = sp.Rational(2048, 45)
    Theta = -Delta
    K = sp.Rational(1376, 135)

    G = primary.source_polynomial(
        X, T, P, Y, Delta, Theta, sp.Rational(0), eta
    )
    resultant_X = sp.resultant(
        sp.cancel(sp.diff(G, X) / T), sp.diff(G, T), X
    )
    Q12 = sp.Poly(sp.cancel(
        resultant_X / (T ** 35 * (6 * T + 1) ** 2)
    ), T)
    primary_factors = {
        sp.degree(factor, eta): (factor, multiplicity)
        for factor, multiplicity in sp.factor_list(
            sp.discriminant(Q12.as_expr(), T), eta
        )[1]
        if sp.degree(factor, eta) > 1
    }
    require({degree: multiplicity for degree, (_, multiplicity)
             in primary_factors.items()} == {30: 1, 36: 2},
            "primary discriminant factor ledger changed")

    t = p - s ** 2
    H = (
        -3 * p + sp.Rational(8, 3) * p ** 2
        - sp.Rational(1376, 135) * p ** 3
        + K * s ** 2 * p ** 2 + Delta * p ** 4
        + Theta * s ** 2 * p ** 3 + eta * s * p ** 3 * t
    )
    Gsp = -s ** 2 / (2 * t) + H
    A = sp.cancel(t ** 2 * sp.diff(Gsp, s) / p)
    C = sp.cancel(2 * t ** 2 * sp.diff(Gsp, p))
    resultant_p = sp.resultant(A, C, s)
    R12 = sp.Poly(sp.cancel(resultant_p / p ** 8), p)
    independent_factors = {
        sp.degree(factor, eta): (factor, multiplicity)
        for factor, multiplicity in sp.factor_list(
            sp.discriminant(R12.as_expr(), p), eta
        )[1]
        if sp.degree(factor, eta) > 1
    }
    require({degree: multiplicity for degree, (_, multiplicity)
             in independent_factors.items()} == {24: 2, 30: 1},
            "independent discriminant factor ledger changed")

    common_primary = sp.Poly(primary_factors[30][0], eta).monic()
    common_independent = sp.Poly(independent_factors[30][0], eta).monic()
    require(common_primary == common_independent,
            "degree-thirty factors are not associates")

    factors = (
        ("H30", primary_factors[30][0]),
        ("H36", primary_factors[36][0]),
        ("H24", independent_factors[24][0]),
    )
    extracted = []
    for name, factor in factors:
        primitive, monic = even_factor_to_u(factor, eta, u)
        extracted.append((name, primitive, monic, monic_digest(monic)))

    q_at_universal = sp.together(
        Q12.eval(-sp.Rational(1, 6))
    ).as_numer_denom()[0]
    J = 22143375 * eta ** 2 + 15510536192
    require(sp.gcd(independent_factors[24][0], primary_factors[36][0]) == 1,
            "H24 and H36 acquired a common root")
    require(sp.gcd(independent_factors[24][0], q_at_universal) == 1,
            "H24 and Q12(-1/6) acquired a common root")
    require(sp.gcd(J, primary_factors[30][0]) == 1,
            "J and H30 acquired a common root")

    semantic = "\n".join(
        f"{name}:{','.join(str(c) for c in primitive.all_coeffs())}"
        for name, primitive, _, _ in extracted
    ) + "\n"
    print("DEEPEST REPEATED-WALL PROJECTION FACTOR EXTRACTOR")
    print("substitution=u=eta^2")
    for name, primitive, monic, digest in extracted:
        print(f"factor={name};degree_eta={2 * primitive.degree()};"
              f"degree_u={primitive.degree()}")
        print(f"  primitive_integer_coefficients={primitive.all_coeffs()}")
        print(f"  monic_rational_coefficients={monic.all_coeffs()}")
        print(f"  monic_sha256={digest}")
    print("gcd(H24,H36)=1")
    print("gcd(H24,numer(Q12(-1/6)))=1")
    print("gcd(J,H30)=1")
    print(f"semantic_sha256={sha256(semantic.encode()).hexdigest()}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
