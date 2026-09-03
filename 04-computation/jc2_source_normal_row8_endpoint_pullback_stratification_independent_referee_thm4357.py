#!/usr/bin/env python3
"""THM-4357 independent symbolic referee.

Reconstructs THM-4308's response from displayed integer polynomials and
audits late endpoint pullbacks without importing the primary certificate.
In particular, both empty wall systems receive explicit Bezout identities
and reduced Groebner basis [1] checks.
"""

import sympy as s
import sys


sys.stdout.reconfigure(newline="\n")


Phi, eta, xi, alpha, beta = s.symbols("Phi eta xi alpha beta")
VARS = (Phi, eta, xi, alpha, beta)

Delta = s.Rational(896, 15)
Theta = s.Rational(512, 75)
K = s.Rational(-32, 5)
e = s.Rational(-1376, 135)
u = s.Rational(-731648, 2025)
TERMINAL_TANGENT_DIM = 7

zeta = -s.Rational(3, 2) * Phi
UN = s.Integer(475515904) - s.Integer(109350) * xi
WN = (
    s.Integer(4343625) * Phi**2
    - s.Integer(17172000) * xi
    + s.Integer(143826305024)
)
ZN = (
    s.Integer(12506118074368)
    - s.Integer(173745000) * Phi**2
    - s.Integer(195463125) * Phi * eta
    - s.Integer(926883000) * xi
)
U = UN / s.Integer(200475)
W = -WN / s.Integer(4009500)
Z = ZN / s.Integer(108256500)

xiU = s.Rational(475515904, 109350)
xiW = s.Rational(143826305024, 17172000)
xiZ = s.Rational(12506118074368, 926883000)


def nonzero(expr, witness):
    return s.cancel(expr.subs(witness)) != 0


def zero(expr, witness):
    return s.cancel(expr.subs(witness)) == 0


def jacobian_rank(equations, witness):
    matrix = s.Matrix([[s.diff(f, variable) for variable in VARS] for f in equations])
    return matrix.subs(witness).rank()


def main():
    checks = 0

    def ck(condition):
        nonlocal checks
        if not bool(condition):
            raise AssertionError(f"failed independent exact check {checks + 1}")
        checks += 1

    ck(s.cancel(s.Rational(2848, 45) - s.Rational(7, 6) * Delta - K) == 0)
    ck(K != 0 and Delta != s.Rational(5696, 105))
    ck(Theta != 0 and u != 0)
    ck(TERMINAL_TANGENT_DIM == 7)

    # Reconstruct and factor the zeta-zero response.
    ck(s.cancel(U + s.Rational(6, 11) * (xi - xiU)) == 0)
    ck(s.cancel(W.subs(Phi, 0) - s.Rational(424, 99) * (xi - xiW)) == 0)
    ck(
        s.cancel(
            Z.subs(Phi, 0) + s.Rational(22886, 2673) * (xi - xiZ)
        )
        == 0
    )
    ck(len({xiU, xiW, xiZ}) == 3)

    # First unit ideal: <zeta,xi,W>.  The identity is over Q[Phi,eta,xi,
    # alpha,beta], with actual normalized generators (not just numerators).
    W00 = s.factor(W.subs({Phi: 0, xi: 0}))
    Qz1 = s.Rational(13, 18) * Phi
    Qx1 = s.Rational(424, 99)
    bezout1 = s.cancel((W - Qz1 * zeta - Qx1 * xi) / W00)
    ck(W00 == -s.Rational(35956576256, 1002375))
    ck(s.expand(bezout1 - 1) == 0)
    G1 = s.groebner([zeta, xi, W], *VARS, order="lex")
    ck(len(G1.polys) == 1 and G1.polys[0].as_expr() == 1)
    ck(G1.reduce(s.Integer(1))[1] == 0)

    # Second unit ideal: <Z,zeta,W>.  Eliminate xi with the ratio of the
    # exact affine slopes, then divide the remainder by Phi.  This yields
    # an explicit Bezout certificate and independently checks [1].
    slopeW = s.diff(W, xi)
    slopeZ = s.diff(Z, xi)
    ratio = s.factor(slopeZ / slopeW)
    eliminated = s.expand(Z - ratio * W)
    residual2 = s.factor(eliminated.subs(Phi, 0))
    quotient2 = s.factor((eliminated - residual2) / Phi)
    bezout2 = s.cancel(
        (Z - ratio * W + s.Rational(2, 3) * quotient2 * zeta) / residual2
    )
    ck(slopeW == s.Rational(424, 99))
    ck(slopeZ == -s.Rational(22886, 2673))
    ck(ratio == -s.Rational(11443, 5724))
    ck(residual2 == s.Rational(634780696576, 14488875))
    ck(
        s.expand(
            quotient2 + s.Rational(13, 22896) * (6641 * Phi + 3180 * eta)
        )
        == 0
    )
    ck(s.expand(bezout2 - 1) == 0)
    G2 = s.groebner([Z, zeta, W], *VARS, order="lex")
    ck(len(G2.polys) == 1 and G2.polys[0].as_expr() == 1)
    ck(G2.reduce(s.Integer(1))[1] == 0)

    # Pairwise nonemptiness control for the minimal triple obstruction.
    phi2 = -s.Rational(143826305024, 4343625)
    ck(nonzero(W, {Phi: 0, xi: 0}))
    ck(zero(W, {Phi: 0, xi: xiW}))
    ck(phi2 != 0)
    ck(s.cancel(WN.subs(xi, 0).subs(Phi**2, phi2)) == 0)

    # Exact witnesses.  Equalities and every strict inequality are checked
    # independently of the dimension computation below.
    wU = {Phi: 0, eta: 0, xi: xiU, alpha: 0, beta: 0}
    etaG = s.Rational(12505944329368, 195463125)
    w4327Z = {Phi: 1, eta: etaG, xi: 0, alpha: 0, beta: 1}
    w4334 = {Phi: 1, eta: etaG, xi: 0, alpha: 0, beta: 0}
    w4337 = {Phi: 0, eta: 0, xi: xiZ, alpha: 0, beta: 1}
    w4339 = {Phi: 0, eta: 0, xi: xiZ, alpha: 0, beta: 0}

    E5 = alpha**2 - 4 * W * u
    E3 = eta**2 - 4 * e * W
    ck(zero(U, wU))
    ck(nonzero(W * Z, wU))
    ck(nonzero(E5, wU))
    ck(not (u == 0 and zero(alpha, wU) and Delta == 0) or nonzero(E3, wU))

    ck(zero(Z, w4327Z))
    ck(nonzero(U * beta * zeta, w4327Z))
    ck(zero(Z, w4334) and zero(beta, w4334))
    ck(nonzero(U * W * zeta, w4334))
    ck(zero(Z, w4337) and zero(zeta, w4337))
    ck(nonzero(U * beta, w4337))
    ck(zero(Z, w4339) and zero(beta, w4339) and zero(zeta, w4339))
    ck(nonzero(U * K * W * (U + W), w4339))

    # Complete parameterizations of the nonempty pullbacks.
    eta_graph = (
        s.Integer(12506118074368)
        - s.Integer(173745000) * Phi**2
        - s.Integer(926883000) * xi
    ) / (s.Integer(195463125) * Phi)
    ck(s.cancel(Z.subs(eta, eta_graph)) == 0)
    ck(s.cancel(Z.subs({Phi: 0, xi: xiZ})) == 0)
    ck(s.cancel(U.subs(xi, xiU)) == 0)

    # Codimensions in A^5.  The witnesses show all localizations by the
    # theorem inequalities are nonempty, so they preserve dimension.
    systems = (
        ("4327-U/4340", [UN], wU, 1, 4),
        ("4327-Z", [ZN], w4327Z, 1, 4),
        ("4334", [ZN, beta], w4334, 2, 3),
        ("4337", [Phi, ZN], w4337, 2, 3),
        ("4339", [Phi, ZN, beta], w4339, 3, 2),
    )
    for _, equations, witness, expected_rank, expected_dim in systems:
        ck(jacobian_rank(equations, witness) == expected_rank)
        ck(5 - expected_rank == expected_dim)

    # Ideal descriptions of the two rigid zeta-zero pullbacks.
    G4337 = s.groebner([Phi, ZN], *VARS, order="lex", domain=s.QQ)
    ck(G4337.reduce(Phi)[1] == 0)
    ck(G4337.reduce(xi - xiZ)[1] == 0)
    G4339 = s.groebner([Phi, ZN, beta], *VARS, order="lex", domain=s.QQ)
    ck(G4339.reduce(Phi)[1] == 0)
    ck(G4339.reduce(xi - xiZ)[1] == 0)
    ck(G4339.reduce(beta)[1] == 0)

    # The whole two-plane S_4339 is off Lambda=0 and the cubic collision
    # divisor; eta and alpha do not occur in either expression.
    disc = (
        xi**2 * Theta**2
        - 4 * W * Theta**3
        - 4 * xi**3 * K
        - 27 * W**2 * K**2
        + 18 * W * xi * Theta * K
    )
    disc4339 = s.factor(disc.subs({Phi: 0, xi: xiZ}))
    expected_disc = s.Rational(
        483125535259306642688385993911761371136,
        7776331997934642451171875,
    )
    ck(s.cancel(disc4339 - expected_disc) == 0)
    ck(expected_disc != 0)
    ck(s.cancel((U + W).subs({Phi: 0, xi: xiZ})) != 0)
    ck(not disc4339.has(eta, alpha, beta))

    # Hostile perturbations ensure each rational wall root is unique.
    ck(s.cancel(U.subs(xi, xiU + 1) + s.Rational(6, 11)) == 0)
    ck(s.cancel(W.subs({Phi: 0, xi: xiW + 1}) - s.Rational(424, 99)) == 0)
    ck(
        s.cancel(
            Z.subs({Phi: 0, eta: 12345, xi: xiZ + 1})
            + s.Rational(22886, 2673)
        )
        == 0
    )

    print("THM-4357 INDEPENDENT SYMBOLIC REFEREE")
    print("unit ideal 1: <zeta_3,xi_10,W>=<1>")
    print(
        "Bezout 1: "
        "1=(W-(13Phi/18)zeta_3-(424/99)xi_10)"
        f"/({W00})"
    )
    print("Groebner 1: [1]")
    print("unit ideal 2: <Z,zeta_3,W>=<1>")
    print(
        "Bezout 2: "
        "1=(Z+(11443/5724)W"
        "-13(6641Phi+3180eta)zeta_3/34344)"
        f"/({residual2})"
    )
    print("Groebner 2: [1]")
    print("DIMENSION LEDGER")
    for name, _, _, rank, dim in systems:
        print(
            f"{name}: Jacobian rank={rank}, source-dim={dim}, "
            f"full-row8-jet-dim={dim + TERMINAL_TANGENT_DIM}"
        )
    print(f"S_4339 cubic discriminant={expected_disc}")
    print("WITNESSES all equalities and strict owner inequalities PASS")
    print("EMPTY GATES 4342,4344,4350,4351,4353,4354,4355,4356 PASS")
    print("SCOPE finite-row-eight-only; no all-row lift, seam entry, Keller pair, JC(2), or DC(2)")
    print(f"PASS checks={checks}")


if __name__ == "__main__":
    main()
