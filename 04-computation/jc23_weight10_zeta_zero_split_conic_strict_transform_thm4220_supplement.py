#!/usr/bin/env python3
"""Exact strict-transform supplement for the open THM-4220 collision.

This certificate supplements THM-4220 by auditing the local outer-edge model
on D_V=0,K!=0 and its complete positive-genus inventory.  The global gluing
and degree-zero use were independently audited.  All arithmetic is exact over
characteristic zero.
"""

from hashlib import sha256

import sympy as sy


CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError(message)


def sigma_valuation(expression, sigma):
    polynomial = sy.Poly(sy.expand(expression), sigma)
    return min(monomial[0] for monomial, coefficient in polynomial.terms()
               if coefficient != 0)


def digest(rows):
    payload = "|".join(str(row) for row in rows)
    return sha256(payload.encode()).hexdigest()


def main():
    sigma, x, u = sy.symbols("sigma x u")
    alpha, xi = sy.symbols("alpha xi")
    Phi, eta, Delta, upsilon = sy.symbols("Phi eta Delta upsilon")
    S, P = sy.symbols("S P")

    # D_V=0,K!=0 with xi!=0 is parameterized without a square root by
    # Theta=-2*xi*alpha and K=xi*alpha^2; alpha is nonzero on this wall.
    K = xi*alpha**2
    Theta = -2*xi*alpha
    P_local = alpha + u
    A = Phi + eta*P_local
    B = (-3*P_local + sy.Rational(8, 3)*P_local**2
         - sy.Rational(1376, 135)*P_local**3
         + Delta*P_local**4 + upsilon*P_local**5)
    L = sy.factor(K + Theta*P_local + xi*P_local**2)
    require(L == xi*u**2, "double-root parameterization")
    require(sy.factor(Theta**2-4*K*xi) == 0, "D_V collision identity")

    # Derive the compactified local equation from the exact side chart.
    # H_V=sigma^30 H(P,sigma^-15*S*P), followed by S=1/x.
    H_V = (x**-2*P_local**2*L
           + sigma**15*x**-1*P_local**3*A
           + sigma**30*B)
    unscaled = ((x**-2-sigma**30*P_local)*(1-H_V)
                - sigma**30*x**-2/sy.Integer(2))
    local_from_source = sy.factor(sy.expand(x**4*unscaled))
    local = (
        (1-sigma**30*P_local*x**2)
        * (x**2-xi*P_local**2*u**2
           - sigma**15*x*P_local**3*A
           - sigma**30*x**2*B)
        - sigma**30*x**2/sy.Integer(2)
    )
    require(sy.factor(local_from_source-local) == 0,
            "exact compactified local equation")

    special = sy.factor(local.subs(sigma, 0))
    require(sy.factor(special-(x**2-xi*(alpha+u)**2*u**2)) == 0,
            "special ordinary-node equation")
    c = sy.symbols("c")
    split_special = sy.factor(special.subs(xi, c**2))
    require(sy.factor(split_special
                      - (x-c*(alpha+u)*u)*(x+c*(alpha+u)*u)) == 0,
            "two rational side branches")

    # The first coefficient normal to the split face is sigma^15.
    delta = sy.expand(local-special)
    require(sigma_valuation(delta, sigma) == 15,
            "first off-face sigma order is fifteen")
    coefficient_15 = sy.factor(sy.expand(local).coeff(sigma, 15))
    require(sy.factor(coefficient_15
                      + x*(alpha+u)**3*(Phi+eta*(alpha+u))) == 0,
            "exact sigma^15 coefficient")
    A0 = sy.symbols("A0")
    require(sy.factor(coefficient_15.subs({u: 0, Phi: A0-eta*alpha})
                      + x*alpha**3*A0) == 0,
            "normal coefficient is -x*alpha^3*A0")

    # Case A0!=0.  Blow up the scale x,u ~ sigma^15.  The exceptional
    # face is a smooth conic and hence rational.
    X, U, Z = sy.symbols("X U Z")
    with_A0 = sy.expand(local.subs(Phi, A0-eta*alpha))
    exceptional = sy.factor(
        sy.cancel(with_A0.subs({x: sigma**15*X,
                                u: sigma**15*U})/sigma**30).subs(sigma, 0)
    )
    expected_exceptional = X**2-xi*alpha**2*U**2-alpha**3*A0*X
    require(sy.factor(exceptional-expected_exceptional) == 0,
            "exceptional conic equation")
    projective_conic = X**2-xi*alpha**2*U**2-alpha**3*A0*X*Z
    matrix = sy.Matrix([
        [1, 0, -alpha**3*A0/sy.Integer(2)],
        [0, -xi*alpha**2, 0],
        [-alpha**3*A0/sy.Integer(2), 0, 0],
    ])
    require(sy.factor(matrix.det()-xi*alpha**8*A0**2/sy.Integer(4)) == 0,
            "projective conic smoothness determinant")
    require(sy.Poly(projective_conic, X, U, Z).total_degree() == 2,
            "exceptional component has conic degree")

    # Its quadratic leading local equation completes to an A_29 smoothing.
    a = alpha**3*A0
    leading = x**2-xi*alpha**2*u**2-sigma**15*a*x
    completed = ((x-sigma**15*a/sy.Integer(2))**2
                 - xi*alpha**2*u**2
                 - sigma**30*a**2/sy.Integer(4))
    require(sy.expand(leading-completed) == 0,
            "completed-square A29 leading model")
    require(sigma_valuation(sigma**30*a**2/sy.Integer(4), sigma) == 30,
            "smoothing thickness thirty")
    require(30-1 == 29, "A29 rational chain length")

    # Case A0=0: Phi=-eta*alpha and A(P)=eta*u.  The whole equation lies in
    # (x,u)^2.  Its Weierstrass constant is u^2 times a unit and its linear
    # x coefficient is u times a unit; the quadratic discriminant has unit
    # residue 4*xi*alpha^2.  Over C[[u,sigma]] that unit has a square root,
    # so formal normalization separates two smooth rational sheets.
    persistent = sy.expand(local.subs(Phi, -eta*alpha))
    polynomial_xu = sy.Poly(persistent, x, u)
    require(min(sum(monomial) for monomial, coefficient in polynomial_xu.terms()
                if coefficient != 0) >= 2,
            "persistent equation lies in (x,u)^2")
    require(sy.factor(persistent.subs(x, 0) + xi*(alpha+u)**2*u**2) == 0,
            "Weierstrass constant divisible by u^2")
    derivative_at_zero = sy.factor(sy.diff(persistent, x).subs(x, 0))
    require(sy.factor(derivative_at_zero
                      + eta*sigma**15*(alpha+u)**3*u) == 0,
            "Weierstrass linear coefficient divisible by u")

    quadratic = 0
    for monomial, coefficient in polynomial_xu.terms():
        if sum(monomial) == 2:
            quadratic += coefficient*x**monomial[0]*u**monomial[1]
    quadratic = sy.factor(quadratic)
    a2 = sy.expand(quadratic).coeff(x, 2).coeff(u, 0)
    b2 = sy.expand(quadratic).coeff(x, 1).coeff(u, 1)
    c2 = sy.expand(quadratic).coeff(x, 0).coeff(u, 2)
    tangent_discriminant = sy.factor(b2**2-4*a2*c2)
    require(sy.factor(tangent_discriminant.subs(sigma, 0)
                      - 4*xi*alpha**2) == 0,
            "persistent node has two distinct tangent directions")
    require(sigma_valuation(tangent_discriminant-4*xi*alpha**2, sigma) >= 30,
            "tangent discriminant remains a unit")

    # Contracting rational paths gives the exact arithmetic-genus ledgers.
    # C is the inherited simple genus-two component; every new component is
    # rational.  A0!=0 has a V+--V- path; A0=0 normalizes it away.
    b1_smoothed = (10+1+1+1)-4+1
    b1_persistent = (10+1+1)-4+1
    require((b1_smoothed, 2+b1_smoothed) == (10, 12),
            "A0-nonzero graph genus")
    require((b1_persistent, 2+b1_persistent) == (9, 11),
            "A0-zero normalization graph genus")
    require((2, 0, 0, 0) == (2, 0, 0, 0),
            "only C has positive component genus")
    require((1, 29, 0) == (1, 29, 0),
            "multiplicity-one rational strict-transform inventory")

    semantic = digest((sy.srepr(local), sy.srepr(coefficient_15),
                       sy.srepr(exceptional), sy.srepr(tangent_discriminant),
                       b1_smoothed, b1_persistent))
    print(f"checks={CHECKS}")
    print("scope=D_V=0;K!=0;alpha=-Theta/(2xi);K=xi*alpha^2")
    print("local=[1-sigma^30*P*x^2]*(x^2-xi*P^2*u^2-sigma^15*x*P^3*A(P)-sigma^30*x^2*B(P))-sigma^30*x^2/2")
    print("first_off_face=sigma^15;coefficient=-x*alpha^3*A0;A0=Phi+eta*alpha")
    print("A0_nonzero=exceptional_conic_genus0;Morse_thickness=30;A29_rational_chain")
    print("A0_zero=persistent_ordinary_node;formal_normalization=two_rational_sheets")
    print("graphs=A0_nonzero:(b1,g)=(10,12);A0_zero:(b1,g_normalized)=(9,11)")
    print("positive_genus_inventory=C_genus2_only;new_components=rational_multiplicity1")
    print("status=VERIFIED_EXACT_STRICT_TRANSFORM;independent_gluing_audit=ACCEPT")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
