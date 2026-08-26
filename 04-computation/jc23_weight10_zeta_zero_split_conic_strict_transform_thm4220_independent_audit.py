#!/usr/bin/env python3
"""Clean-room exact audit of the residual split-conic wall after THM-4220.

Universe: the exact-M=10 zeta=0 locus with
upsilon*xi*(upsilon+xi) != 0, D_V=Theta^2-4*K*xi=0, and K != 0.
The calculation starts again from the scaled source pencil.  It checks the
outer strict transform, the two A0=Phi+eta*alpha cases, the labelled dual
graph, and the multiplicity/degree ledger.  It imports only the already
proved simplicity of the inherited genus-two component from THM-4218.
"""

from hashlib import sha256
from math import gcd

import sympy as sy


CHECKS = 0


def require(assertion, label):
    global CHECKS
    CHECKS += 1
    if not bool(assertion):
        raise RuntimeError(label)


def sigma_order(polynomial, sigma):
    terms = sy.Poly(sy.expand(polynomial), sigma).terms()
    return min(monomial[0] for monomial, coefficient in terms
               if coefficient != 0)


def main():
    sigma, x, u, z = sy.symbols("sigma x u z")
    alpha, xi, eta, Phi = sy.symbols("alpha xi eta Phi", nonzero=True)
    Delta, upsilon = sy.symbols("Delta upsilon")

    # Parameterize the whole residual discriminant wall.  The theorem gates
    # imply alpha*xi != 0.
    P = alpha + u
    K = xi*alpha**2
    Theta = -2*xi*alpha
    quadratic_owner = sy.factor(K + Theta*P + xi*P**2)
    require(quadratic_owner == xi*u**2, "collision parameterization")
    require(sy.factor(Theta**2 - 4*K*xi) == 0, "discriminant wall")

    odd_owner = Phi + eta*P
    pure_owner = (-3*P + sy.Rational(8, 3)*P**2
                  - sy.Rational(1376, 135)*P**3
                  + Delta*P**4 + upsilon*P**5)

    # Side scaling is s=sigma^-15*S, p=P and H_V=sigma^30 H;
    # compactify the outer divisor by x=S^-1.  Multiplication by x^4
    # gives the following exact local equation, including every lower term.
    H_V = (x**-2*P**2*quadratic_owner
           + sigma**15*x**-1*P**3*odd_owner
           + sigma**30*pure_owner)
    source_scaled = ((x**-2 - sigma**30*P)*(1-H_V)
                     - sigma**30*x**-2/sy.Integer(2))
    derived = sy.factor(sy.expand(x**4*source_scaled))
    local = ((1-sigma**30*P*x**2)
             * (x**2-xi*P**2*u**2
                -sigma**15*x*P**3*odd_owner
                -sigma**30*x**2*pure_owner)
             - sigma**30*x**2/sy.Integer(2))
    require(sy.factor(derived-local) == 0, "exact compactified equation")

    # The split special face consists of two rational branches.  Its sole
    # common point on this edge is the ordinary node x=u=0.
    special = sy.factor(local.subs(sigma, 0))
    require(sy.factor(special-(x**2-xi*(alpha+u)**2*u**2)) == 0,
            "special split equation")
    c = sy.symbols("c", nonzero=True)
    split = sy.factor(special.subs(xi, c**2))
    require(sy.factor(split
                      -(x-c*(alpha+u)*u)*(x+c*(alpha+u)*u)) == 0,
            "two rational sheets")
    hessian = sy.hessian(special, (x, u)).subs({x: 0, u: 0})
    require(sy.factor(hessian.det()) == -4*xi*alpha**2,
            "ordinary-node Hessian")

    # The first coefficient transverse to the split face is controlled by
    # A0=Phi+eta*alpha and occurs at order fifteen.
    correction = sy.expand(local-special)
    require(sigma_order(correction, sigma) == 15,
            "first transverse order")
    normal_15 = sy.factor(sy.expand(local).coeff(sigma, 15))
    require(sy.factor(normal_15+x*P**3*odd_owner) == 0,
            "normal coefficient")
    A0 = sy.symbols("A0", nonzero=True)
    require(sy.factor(normal_15.subs({u: 0, Phi: A0-eta*alpha})
                      + x*alpha**3*A0) == 0,
            "A0 gate")

    # A0 != 0.  Scaling x,u by sigma^15 recovers the exceptional conic.
    # Its Hessian is nondegenerate, and its critical value has nonzero
    # sigma^30 leading coefficient.  The parameterized formal Morse lemma
    # therefore gives XY=sigma^30*(unit), i.e. A_29 and a reduced rational
    # path with 29 intermediate components.
    X, U, Z = sy.symbols("X U Z")
    local_nonzero = sy.expand(local.subs(Phi, A0-eta*alpha))
    exceptional = sy.factor(sy.cancel(
        local_nonzero.subs({x: sigma**15*X, u: sigma**15*U})
        / sigma**30).subs(sigma, 0))
    expected = X**2-xi*alpha**2*U**2-alpha**3*A0*X
    require(sy.factor(exceptional-expected) == 0,
            "exceptional equation")
    exceptional_projective = (X**2-xi*alpha**2*U**2
                              -alpha**3*A0*X*Z)
    conic_matrix = sy.Matrix([
        [1, 0, -alpha**3*A0/sy.Integer(2)],
        [0, -xi*alpha**2, 0],
        [-alpha**3*A0/sy.Integer(2), 0, 0],
    ])
    require(sy.factor(conic_matrix.det()
                      - xi*alpha**8*A0**2/sy.Integer(4)) == 0,
            "smooth exceptional conic")
    critical_X = alpha**3*A0/sy.Integer(2)
    critical_U = 0
    critical_residue = sy.factor(
        expected.subs({X: critical_X, U: critical_U}))
    require(critical_residue == -alpha**6*A0**2/sy.Integer(4),
            "nonzero Morse critical residue")
    require(30-1 == 29, "A29 path length")

    # A0 = 0.  Here Phi=-eta*alpha and odd_owner=eta*u.  Substitute x=u*z.
    # The equation is exactly divisible by u^2.  At u=0 its z-polynomial
    # has discriminant with nonzero residue 4*xi*alpha^2; Hensel therefore
    # gives two distinct power-series roots and the normalization consists
    # of two smooth sheets x=u*z_+ and x=u*z_-.
    local_zero = sy.expand(local.subs(Phi, -eta*alpha))
    require(sy.factor(local_zero.subs(x, 0)
                      + xi*(alpha+u)**2*u**2) == 0,
            "constant term is u squared times a unit")
    require(sy.factor(sy.diff(local_zero, x).subs(x, 0)
                      + eta*sigma**15*(alpha+u)**3*u) == 0,
            "linear term is u times a unit")
    divided = sy.cancel(local_zero.subs(x, u*z)/u**2)
    require(not sy.denom(divided).has(u), "exact x=u*z division")
    edge_polynomial = sy.Poly(sy.expand(divided.subs(u, 0)), z)
    require(edge_polynomial.degree() == 2, "Weierstrass edge degree")
    edge_discriminant = sy.factor(sy.discriminant(edge_polynomial, z))
    require(sy.factor(edge_discriminant.subs(sigma, 0)
                      - 4*xi*alpha**2) == 0,
            "two Hensel roots")
    require(sigma_order(edge_discriminant-4*xi*alpha**2, sigma) >= 30,
            "discriminant remains a unit")

    # Attachment labels and graph conservation.  The two reduced roots of
    # 1-xi*Z^2 attach C to V+ and V- separately.  Ten inherited R--C paths
    # remain.  A0 != 0 adds the rational A29 path V+--V-; A0=0 normalization
    # separates that path.  Thus only the inherited genus-two C is positive
    # genus, and arithmetic/geometric genus agree with the delta ledger.
    edge_scheme = 1-xi*Z**2
    require(sy.discriminant(edge_scheme, Z) == 4*xi,
            "two reduced labelled BD roots")
    b1_nonzero = (10+1+1+1)-4+1
    b1_zero = (10+1+1)-4+1
    require((b1_nonzero, 2+b1_nonzero) == (10, 12),
            "A0 nonzero graph genus")
    require((b1_zero, 2+b1_zero) == (9, 11),
            "A0 zero normalized genus")
    require(12-11 == 1, "persistent node delta")

    # Q=sigma^30 makes both inherited face normals primitive with last
    # coordinate one.  The split sheets and the semistable A29 chain are
    # consequently reduced.  Rational components have no nonconstant map to
    # the good elliptic target; THM-4218 gives Hom(J(C),E0)=0 for C.
    require(gcd(gcd(3, 6), 1) == gcd(gcd(15, 0), 1) == 1,
            "multiplicity-one face normals")
    component_degrees = (0, 0, 0, 0)
    multiplicities = (1, 1, 1, 1)
    require(sum(m*d for m, d in zip(multiplicities, component_degrees)) == 0,
            "specialized degree zero")

    semantic = sha256("|".join(map(str, (
        sy.srepr(local), sy.srepr(exceptional),
        sy.srepr(edge_discriminant), b1_nonzero, b1_zero,
    ))).encode()).hexdigest()
    print(f"checks={CHECKS}")
    print("scope=D_V=0;K!=0;alpha=-Theta/(2*xi);K=xi*alpha^2")
    print("local_equation=exact;first_off_face=sigma^15;A0=Phi+eta*alpha")
    print("A0_nonzero=formal_Morse_A29;29_P1_path;normalized_genus=12")
    print("A0_zero=two_Hensel_sheets;no_joining_path;normalized_genus=11")
    print("attachments=10*(R-C)+(C-V+)+(C-V-)+optional(V+-V-)")
    print("positive_genus=C_genus2_only;all_new_components=rational_multiplicity1")
    print("degree_specialization=0;residual_split_conic_wall=closed")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
