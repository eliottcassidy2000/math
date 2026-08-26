#!/usr/bin/env python3
"""Exact independent certificate for the deepest repeated-wall J=0 slice.

The calculation is performed over K=Q[eta]/(J), where
J=22143375*eta^2+15510536192.  It reconstructs the row-D source rather than
importing a saved resultant, reduces the primary resultant coefficientwise,
and factors the resulting degree-eleven polynomial over K.
"""

from hashlib import sha256

import sympy as sp


CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def rational_token(value):
    value = sp.Rational(value)
    return f"{value.p}/{value.q}"


def main():
    X, T, eta = sp.symbols("X T eta")
    P = T + X ** 2 * T ** 2
    Y = X * T * P
    Delta = sp.Rational(2048, 45)
    Theta = -Delta
    K = sp.Rational(1376, 135)

    # Exact row-D source on zeta=-eta.
    G = (
        -X ** 2 * T / 2 - 3 * P + sp.Rational(8, 3) * P ** 2
        - sp.Rational(1376, 135) * P ** 3 + K * Y ** 2
        + Delta * P ** 4 + Theta * P * Y ** 2
        + eta * P ** 3 * Y - eta * Y ** 3
    )
    J = 22143375 * eta ** 2 + 15510536192
    Jpoly = sp.Poly(J, eta, domain=sp.QQ)
    require(Jpoly.degree() == 2 and sp.factor(J) == J,
            "J is not the expected irreducible quadratic")

    def reduce_mod_J(expression):
        polynomial = sp.Poly(sp.cancel(expression), eta, domain=sp.QQ)
        return polynomial.rem(Jpoly).as_expr()

    def basis_pair(expression):
        polynomial = sp.Poly(reduce_mod_J(expression), eta, domain=sp.QQ)
        require(polynomial.degree() <= 1,
                "coefficient did not reduce to the basis 1,eta")
        return sp.Rational(polynomial.nth(0)), sp.Rational(polynomial.nth(1))

    GX = sp.diff(G, X)
    f_expression = sp.cancel(GX / T)
    f_denominator = f_expression.as_numer_denom()[1]
    require(not (f_denominator.free_symbols & {X, T, eta}),
            "G_X/T ceased to be polynomial")
    f = sp.Poly(f_expression, X).as_expr()
    h = sp.Poly(sp.diff(G, T), X).as_expr()
    fpoly = sp.Poly(f, X)
    hpoly = sp.Poly(h, X)
    require(fpoly.degree() == 6 and hpoly.degree() == 7,
            "finite-T X degrees changed")
    require(sp.factor(fpoly.LC() - 7 * T ** 7 * eta) == 0,
            "leading coefficient of f changed")
    require(sp.factor(hpoly.LC() - 8 * T ** 7 * eta) == 0,
            "leading coefficient of h changed")

    resultant = sp.resultant(f, h, X)
    quotient = sp.cancel(resultant / (T ** 35 * (6 * T + 1) ** 2))
    require(sp.denom(quotient) == 1,
            "primary resultant quotient ceased to be polynomial")
    Q12 = sp.Poly(quotient, T)
    require(Q12.degree() == 12, "generic primary quotient degree changed")

    reduced_coefficients = [
        reduce_mod_J(Q12.nth(exponent))
        for exponent in range(Q12.degree() + 1)
    ]
    QJ_expression = sp.expand(sum(
        coefficient * T ** exponent
        for exponent, coefficient in enumerate(reduced_coefficients)
    ))
    QJ = sp.Poly(QJ_expression, T)
    require(QJ.degree() == 11 and reduce_mod_J(Q12.LC()) == 0,
            "J specialization did not drop Q12 exactly to degree eleven")

    expected_leading = (
        sp.Rational(1887494744535166760726994878464,
                    12258226409765625) * eta
    )
    expected_constant = -sp.Rational(
        7518022905104433152, 2767921875
    ) * eta
    expected_at_universal = -sp.Rational(
        2983872572026179618984706441216,
        402131117372361328125
    ) * eta
    require(reduce_mod_J(QJ.LC() - expected_leading) == 0,
            "degree-eleven leading coefficient changed")
    require(reduce_mod_J(QJ.nth(0) - expected_constant) == 0,
            "degree-eleven constant coefficient changed")
    require(reduce_mod_J(QJ.eval(-sp.Rational(1, 6))
                         - expected_at_universal) == 0,
            "Q11(-1/6) changed")

    # The quotient field is Q(sqrt(-30)); eta is the displayed nonzero scalar
    # multiple of its generator.  This is a second reduction path and lets
    # SymPy invoke exact algebraic-field factorization.
    theta = sp.sqrt(-30)
    eta_value = sp.Rational(88064, 18225) * theta
    require(sp.simplify(J.subs(eta, eta_value)) == 0,
            "chosen algebraic generator does not satisfy J")
    Q_extension = sp.Poly(
        QJ_expression.subs(eta, eta_value), T, extension=theta
    )
    require(Q_extension.degree() == 11,
            "algebraic-field copy has the wrong degree")
    factor_unit, factor_list = Q_extension.factor_list()
    factor_ledger = tuple(
        (factor.degree(), multiplicity)
        for factor, multiplicity in factor_list
    )
    require(factor_ledger == ((11, 1),),
            "Q11 is reducible over Q(eta)")
    require(sp.gcd(Q_extension, Q_extension.diff()).degree() == 0,
            "Q11 is not squarefree over Q(eta)")

    # The universal T=-1/6 pair belongs to the resultant; T=0 was removed
    # when G_X was divided by T and is restored in the actual critical ideal.
    gcd_minus_sixth = sp.Poly(
        f.subs(T, -sp.Rational(1, 6)), X
    ).gcd(sp.Poly(h.subs(T, -sp.Rational(1, 6)), X)).monic()
    require(sp.expand(gcd_minus_sixth.as_expr() - (X ** 2 - 6)) == 0,
            f"universal T=-1/6 critical pair changed: {gcd_minus_sixth}")
    require(sp.factor(f.subs(T, 0) + X) == 0,
            "f(X,0) changed")
    require(sp.factor(h.subs(T, 0) + (X ** 2 + 6) / 2) == 0,
            "h(X,0) changed")

    hessian = sp.det(sp.hessian(G, (X, T)))
    hessian_minus_sixth = sp.Poly(
        sp.cancel(hessian.subs(T, -sp.Rational(1, 6)) + 6), X
    ).rem(sp.Poly(X ** 2 - 6, X))
    hessian_zero = sp.Poly(
        sp.cancel(hessian.subs(T, 0) - 6), X
    ).rem(sp.Poly(X ** 2 + 6, X))
    require(hessian_minus_sixth.is_zero,
            "Hessian at T=-1/6 is no longer -6")
    require(hessian_zero.is_zero,
            "Hessian at T=0 is no longer 6")

    # Because eta and T are nonzero at every residual root, the simultaneous
    # leading coefficients 7*T^7*eta and 8*T^7*eta exclude X=infinity.
    # Squarefreeness then makes the residual intersection length eleven.
    residual_length = QJ.degree()
    universal_length = 2 + 2
    affine_length = residual_length + universal_length
    require(affine_length == 15, "affine critical length changed")

    packet = (7, 4, 2, 2, 2, 1)
    packet_defect = sum(index - 1 for index in packet)
    finite_degree = 12
    carrier_index = 3
    full_degree = 18
    both_handle_capacity = (
        2 * finite_degree - affine_length - 2 + carrier_index
    )
    merger_threshold = finite_degree - 1
    origin_finite_index = (7 - 1) + (4 - 1) + (1 - 1)
    carrier_product_capacity = carrier_index
    overlap_bound = full_degree - affine_length
    commutator_bound = 2 * overlap_bound
    require(packet_defect == 12, "packet defect changed")
    require(both_handle_capacity == 10 < merger_threshold == 11,
            "finite two-handle merger inequality changed")
    require(origin_finite_index == 9 > carrier_product_capacity == 3,
            "finite identity-handle inequality changed")
    require(overlap_bound == 3 and commutator_bound == 6 < packet_defect,
            "full commutator-overlap inequality changed")

    basis_coefficients = [
        basis_pair(QJ.nth(exponent))
        for exponent in range(QJ.degree(), -1, -1)
    ]
    coefficient_payload = "degree=11;" + ";".join(
        f"{rational_token(a)},{rational_token(b)}"
        for a, b in basis_coefficients
    )
    coefficient_digest = sha256(coefficient_payload.encode()).hexdigest()
    semantic_payload = "\n".join((
        "J=22143375*eta^2+15510536192",
        "Q_degree=11",
        f"Q_digest={coefficient_digest}",
        f"factor_ledger={factor_ledger}",
        "affine_length=15",
        "finite_capacity=10<11",
        "finite_identity=9>3",
        "full_commutator=6<12",
    )) + "\n"

    print("DEEPEST REPEATED-WALL J-SLICE CERTIFICATE")
    print("field=Q[eta]/(22143375*eta^2+15510536192)=Q(sqrt(-30))")
    print("eta_in_sqrt_minus_30_basis=(88064/18225)*sqrt(-30)")
    print("source=row_D,zeta=-eta,Delta=2048/45,Theta=-2048/45,Phi=0")
    print("map=critical_ideal_(G_X,G_T)->primary_T_resultant")
    print("lost_data=X_coordinates;sidecar=finite_T_X_leading_coefficients")
    print("deg_X(f),deg_X(h)=6,7")
    print("LC_X(f)=7*T^7*eta;LC_X(h)=8*T^7*eta")
    print("primary_generic=T^35*(6T+1)^2*Q12")
    print("specialized_degree=11")
    print(f"Q11_LC={expected_leading}")
    print(f"Q11_constant={expected_constant}")
    print(f"Q11_at_-1/6={expected_at_universal}")
    print(f"Q11_basis_coefficients_desc_1_eta={basis_coefficients}")
    print(f"Q11_basis_coefficient_sha256={coefficient_digest}")
    print(f"factor_ledger_over_Q_eta={factor_ledger}")
    print("squarefree_over_Q_eta=True")
    print("universal_T=-1/6:X^2=6,count=2,Hessian=-6")
    print("restored_T=0:X^2=-6,count=2,Hessian=6")
    print("finite_T_infinity_loss=0")
    print("affine_critical_length=11+2+2=15")
    print("packet=(7,4,2,2,2,1);defect=12")
    print("finite_response=(n,beta)=(12,3)")
    print("finite_both_handles=2*n-L-2+beta=10<n-1=11")
    print("finite_identity_handle=origin_index_9>carrier_product_index_3")
    print("full_response=n=18;overlap<=3;commutator_index<=6<12")
    print(f"semantic_sha256={sha256(semantic_payload.encode()).hexdigest()}")
    print(f"checks={CHECKS}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
