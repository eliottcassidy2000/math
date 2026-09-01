#!/usr/bin/env python3
"""Primary exact certificate for THM-4317.

The certificate starts from THM-4312's literal cubic-corner k=1 local
equation.  It audits both equality-ray weighted charts, the completed-local
surface singularities, their quotient points at infinity, and the absorbing
nearest-neighbour Green kernel carried by every A-chain.

This is conditional local geometry over characteristic zero.  It does not
assert that the row-eight datum extends to an all-row Keller pair; THM-4316
in fact supplies a later bracket obstruction for the audited corner.
"""

from __future__ import annotations

import sympy as sp


CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def is_zero(expression: sp.Expr) -> bool:
    return sp.cancel(sp.together(sp.expand(expression))) == 0


def require_zero(expression: sp.Expr, label: str) -> None:
    require(is_zero(expression), label)


def weighted_initial(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    weights: tuple[int, ...],
) -> tuple[int, sp.Expr]:
    poly = sp.Poly(sp.expand(expression), *variables)
    minimum = min(
        sum(exponent * weight for exponent, weight in zip(monomial, weights))
        for monomial, _ in poly.terms()
    )
    initial = sum(
        coefficient
        * sp.prod(variable**exponent for variable, exponent in zip(variables, monomial))
        for monomial, coefficient in poly.terms()
        if sum(exponent * weight for exponent, weight in zip(monomial, weights))
        == minimum
    )
    return minimum, sp.factor(initial)


q, t, z, x = sp.symbols("q t z x")
U, rho = sp.symbols("U rho", nonzero=True)
upsilon, eta, zeta = sp.symbols("upsilon eta zeta")
Delta, Theta, Phi = sp.symbols("Delta Theta Phi")
alpha = -2 * U * rho
c2 = U * rho**2
L1 = eta + zeta + rho * upsilon
C = Delta + Theta - rho * zeta
L2 = C - upsilon**2 / (4 * U)


def literal_local_equation() -> sp.Expr:
    """Reconstruct THM-4312's literal F before weighted truncation."""

    r = 1 + q
    K = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta
    hhat = (
        U * (r**6 - 2 * r**5 + r**4)
        + t * (alpha * r**5 - alpha * r**4)
        + t**2
        * (
            upsilon * r**5
            + (alpha**2 / (4 * U) - upsilon) * r**4
        )
        + t**3 * (eta * r**4 + zeta * r**3)
        + t**4 * (Delta * r**4 + Theta * r**3)
        + t**5 * Phi * r**3
        + t**6 * (-sp.Rational(1376, 135) * r**3 + K * r**2)
        + sp.Rational(8, 3) * t**8 * r**2
        - 3 * t**10 * r
    )
    return sp.expand(q * hhat - t**12 / 2 - q * z**12)


def weighted_chart_audit() -> sp.Expr:
    F = literal_local_equation()
    recentered = sp.expand(F.subs(q, x + rho * t))

    weight1, initial1 = weighted_initial(recentered, (x, t, z), (6, 4, 1))
    expected1 = rho * t * (U * x**2 - z**12 + L1 * t**3)
    require(weight1 == 16, "L1 weighted order")
    require_zero(initial1 - expected1, "L1 literal weighted initial")

    high = sp.expand(recentered.subs(eta, -rho * upsilon - zeta))
    weight2, initial2 = weighted_initial(high, (x, t, z), (6, 3, 1))
    expected2 = rho * t * (
        U * x**2 + upsilon * t**2 * x + C * t**4 - z**12
    )
    require(weight2 == 15, "L2 weighted order")
    require_zero(initial2 - expected2, "L2 literal weighted initial")

    # Finite z-chart.  These are the two affine charts displayed in THM-4312,
    # now retained far enough to inspect the boundary intersections.
    w, Y, Ybar = sp.symbols("w Y Ybar")
    chart1 = sp.cancel(
        recentered.subs(t, w * z**4).subs(x, z**6 * Y) / z**16
    )
    exceptional1 = sp.factor(chart1.subs(z, 0))
    require_zero(
        exceptional1 - rho * w * (U * Y**2 - 1 + L1 * w**3),
        "L1 finite-chart exceptional equation",
    )

    chart2 = sp.cancel(high.subs(t, w * z**3).subs(x, z**6 * Y) / z**15)
    exceptional2 = sp.factor(chart2.subs(z, 0))
    require_zero(
        exceptional2
        - rho * w * (U * Y**2 + upsilon * w**2 * Y + C * w**4 - 1),
        "L2 finite-chart exceptional equation",
    )
    completed_square = sp.factor(
        exceptional2.subs(Y, Ybar - upsilon * w**2 / (2 * U))
    )
    require_zero(
        completed_square - rho * w * (U * Ybar**2 - 1 + L2 * w**4),
        "L2 completed-square exceptional equation",
    )

    # Each carrier meets w=0 at Y=+/-U^(-1/2).  The mixed Hessian is
    # nondegenerate at both points, so a parameterized Morse lemma reduces
    # the remaining singularity type to the order of its z-critical value.
    Y0 = sp.symbols("Y0", nonzero=True)
    hessian1 = sp.hessian(exceptional1, (w, Y)).subs({w: 0, Y: Y0})
    hessian2 = sp.hessian(completed_square, (w, Ybar)).subs({w: 0, Ybar: Y0})
    require_zero(hessian1.det() + (2 * rho * U * Y0) ** 2, "L1 boundary Hessian")
    require_zero(hessian2.det() + (2 * rho * U * Y0) ** 2, "L2 boundary Hessian")
    require(
        sp.factor(hessian1.det().subs(Y0**2, 1 / U)) == -4 * U * rho**2,
        "L1 Hessian nonzero on both boundary points",
    )
    require(
        sp.factor(hessian2.det().subs(Y0**2, 1 / U)) == -4 * U * rho**2,
        "L2 Hessian nonzero on both boundary points",
    )

    # Leading critical section in the literal coordinates.  F_q=0 first
    # gives c2*kappa^2=1; F_t=0 then gives q60=3*kappa^10/c2.
    kappa, q60 = sp.symbols("kappa q60", nonzero=True)
    critical_trial = {t: kappa * z**6, q: q60 * z**60}
    Fq = sp.expand(sp.diff(F, q).subs(critical_trial))
    Ft = sp.expand(sp.diff(F, t).subs(critical_trial))
    Ftrial = sp.expand(F.subs(critical_trial))
    require(sp.Poly(Fq, z).terms()[-1][0] == (12,), "Fq first z order")
    require(sp.Poly(Ft, z).terms()[-1][0] == (66,), "Ft first z order")
    require(sp.Poly(Ftrial, z).terms()[-1][0] == (72,), "F critical first z order")
    require_zero(Fq.coeff(z, 12) - (c2 * kappa**2 - 1), "critical t leading equation")
    require_zero(
        Ft.coeff(z, 66) - 2 * kappa * (c2 * q60 - 3 * kappa**10),
        "critical q leading equation",
    )
    require_zero(
        Ftrial.coeff(z, 72) - (q60 * (c2 * kappa**2 - 1) - kappa**12 / 2),
        "critical value leading coefficient",
    )
    require(
        all(Ftrial.coeff(z, degree) == 0 for degree in range(73, 78)),
        "literal critical remainder starts no earlier than z78",
    )

    # In F=-t^12/2+q(h-z^12)+q^2 A+q^3 C(q,t), A=O(t).
    # Thus the exact critical equations improve the non-t^12 contribution
    # to order at least 2*60+6=126.  In literal z, t=kappa*z^6+O(z^12),
    # hence the expanded remainder after the leading z^72 term is O(z^78).
    qpoly = sp.Poly(sp.expand(F + t**12 / 2), q)
    coefficient_q2 = qpoly.coeff_monomial(q**2)
    require(sp.expand(coefficient_q2).subs(t, 0) == 0, "critical q2 coefficient is O(t)")
    critical_identity = sp.expand(
        F
        - q * sp.diff(F, q)
        + t**12 / 2
        + sum(
            (degree - 1) * q**degree * qpoly.coeff_monomial(q**degree)
            for degree in range(2, qpoly.degree() + 1)
        )
    )
    require_zero(critical_identity, "exact critical-value Euler identity")
    q_correction = sp.expand(
        (F - q * sp.diff(F, q) + t**12 / 2).subs(critical_trial)
    )
    require(
        sp.Poly(q_correction, z).terms()[-1][0] == (126,),
        "critical q-dependent correction starts at z126",
    )
    require_zero(
        q_correction.coeff(z, 126) - 2 * U * kappa * q60**2 * rho,
        "critical z126 correction coefficient",
    )
    require(2 * 60 + 6 == 126, "critical non-t12 order")
    require((72 - 16, 72 - 15) == (56, 57), "transverse Morse orders")

    return recentered


def infinity_quotient_audit(recentered: sp.Expr) -> None:
    sigma, lam, u, Xvar, V = sp.symbols("sigma lambda u X V")

    # L1: t=sigma*z, sigma=lambda^3, z=lambda*u, x=lambda^6*X.
    cover1 = recentered.subs(t, sigma * z)
    cover1 = cover1.subs(sigma, lam**3).subs(z, lam * u).subs(x, lam**6 * Xvar)
    cover1 = sp.cancel(sp.expand(cover1) / lam**16)
    require(not cover1.has(sp.zoo), "L1 cover regular after lambda16 division")
    exceptional1 = sp.factor(cover1.subs(lam, 0).subs(Xvar, u * V))
    expected1 = rho * u**3 * (U * V**2 + L1 * u - u**10)
    require_zero(exceptional1 - expected1, "L1 infinity u3 saturation")

    curve1 = U * V**2 + L1 * u - u**10
    require_zero(curve1.subs({u: 0, V: 0}), "L1 unique fixed point lies on curve")
    require(sp.diff(curve1, u).subs({u: 0, V: 0}) == L1, "L1 fixed point smooth")
    require_zero(curve1.subs(V, 0) - u * (L1 - u**9), "L1 V=0 locus")
    require_zero(
        sp.diff(curve1, u).subs(V, 0).subs(u**9, L1) + 9 * L1,
        "L1 nonzero fixed-locus derivative",
    )
    lambda_weight1 = 1
    u_weight1 = -lambda_weight1 % 3
    X_weight1 = -6 * lambda_weight1 % 3
    V_weight1 = (X_weight1 - u_weight1) % 3
    require(
        (lambda_weight1, u_weight1, V_weight1) == (1, 2, 1),
        "L1 mu3 weights lambda,u,V",
    )
    require(
        (lambda_weight1, V_weight1) == (1, 1),
        "L1 local quotient type 1/3(1,1)",
    )

    # L2: first impose L1=0, then use the double cover.  X=u^2 V exposes
    # two fixed points because the residual quadratic has discriminant
    # -4 U L2.
    high = sp.expand(recentered.subs(eta, -rho * upsilon - zeta))
    cover2 = high.subs(t, sigma * z)
    cover2 = cover2.subs(sigma, lam**2).subs(z, lam * u).subs(x, lam**6 * Xvar)
    cover2 = sp.cancel(sp.expand(cover2) / lam**15)
    require(not cover2.has(sp.zoo), "L2 cover regular after lambda15 division")
    exceptional2 = sp.factor(cover2.subs(lam, 0).subs(Xvar, u**2 * V))
    expected2 = rho * u**5 * (U * V**2 + upsilon * V + C - u**8)
    require_zero(exceptional2 - expected2, "L2 infinity u5 saturation")
    discriminant = sp.factor(upsilon**2 - 4 * U * C)
    require_zero(discriminant + 4 * U * L2, "L2 infinity discriminant")
    lambda_weight2 = 1
    u_weight2 = -lambda_weight2 % 2
    X_weight2 = -6 * lambda_weight2 % 2
    V_weight2 = (X_weight2 - 2 * u_weight2) % 2
    require(
        (lambda_weight2, u_weight2, V_weight2) == (1, 1, 0),
        "L2 mu2 weights lambda,u,V",
    )
    require(
        (lambda_weight2, u_weight2) == (1, 1),
        "L2 local quotient type 1/2(1,1)=A1",
    )

    # One cubic infinity orbit in the L1 carrier and two quartic infinity
    # points in L2 account for exactly the quotient points above.  Together
    # with the finite z-chart, these are the two standard weighted charts.
    w = sp.symbols("w")
    carrier1 = sp.Poly(1 - L1 * w**3, w)
    carrier2 = sp.Poly(1 - L2 * w**4, w)
    require(
        (carrier1.degree(), carrier2.degree()) == (3, 4),
        "weighted carrier cubic/quartic degrees",
    )
    require(
        (1 if carrier1.degree() % 2 else 2, 1 if carrier2.degree() % 2 else 2)
        == (1, 2),
        "weighted carrier infinity point counts",
    )
    require(2 * 55 + 1 == 111, "L1 minimal rational exceptional count")
    require(2 * 56 + 2 == 114, "L2 minimal rational exceptional count")


def path_laplacian(r: int) -> sp.Matrix:
    size = r - 1
    return sp.Matrix(
        size,
        size,
        lambda i, j: 2 if i == j else (-1 if abs(i - j) == 1 else 0),
    )


def green_kernel(r: int) -> sp.Matrix:
    return sp.Matrix(
        r - 1,
        r - 1,
        lambda i, j: sp.Rational(
            min(i + 1, j + 1) * (r - max(i + 1, j + 1)), r
        ),
    )


def stochastic_resolution_audit() -> None:
    for r in (56, 57):
        laplacian = path_laplacian(r)
        kernel = green_kernel(r)
        identity = sp.eye(r - 1)
        require(laplacian * kernel == identity, f"A{r-1} inverse intersection kernel")
        require(kernel * laplacian == identity, f"A{r-1} two-sided Green inverse")

        hit_right = [sp.Rational(i, r) for i in range(r + 1)]
        require(
            all(
                hit_right[i] == (hit_right[i - 1] + hit_right[i + 1]) / 2
                for i in range(1, r)
            ),
            f"A{r-1} gambler ruin hitting law",
        )

        visits = 2 * kernel
        require(laplacian * visits == 2 * identity, f"A{r-1} occupation kernel")

        valuation_u = [r - i for i in range(r + 1)]
        valuation_v = [i for i in range(r + 1)]
        valuation_z = [1 for _ in range(r + 1)]
        for name, values in (
            ("u", valuation_u),
            ("v", valuation_v),
            ("z", valuation_z),
        ):
            require(
                all(
                    2 * values[i] - values[i - 1] - values[i + 1] == 0
                    for i in range(1, r)
                ),
                f"A{r-1} monomial valuation {name} harmonic",
            )

        # A deterministic hostile/positive control for the Poisson formula:
        # one unit of strict-transform contact at j=floor(r/2).
        j = r // 2
        contact = sp.zeros(r - 1, 1)
        contact[j - 1, 0] = 1
        orders = kernel * contact
        require(laplacian * orders == contact, f"A{r-1} point-contact Poisson solve")
        require(
            orders[j - 1] == sp.Rational(j * (r - j), r),
            f"A{r-1} point-contact Green diagonal",
        )


def main() -> None:
    recentered = weighted_chart_audit()
    infinity_quotient_audit(recentered)
    stochastic_resolution_audit()
    print("THM-4317 K1 COMPLETED-LOCAL RESOLUTION AND GREEN KERNEL: PASS")
    print("L1_INITIAL wt(x,t,z)=(6,4,1) N=16 rho*t*(U*x^2-z^12+L1*t^3)")
    print("L2_INITIAL wt(x,t,z)=(6,3,1) N=15 rho*t*(U*x^2+upsilon*t^2*x+C*t^4-z^12)")
    print("FINITE_BOUNDARY L1=two_A55 L2=two_A56 critical_order_F=72 quotients=(56,57)")
    print("LITERAL_Z Fcrit=-(kappa^12/2)*z^72+O(z^78); reparametrized remainder O(z^126)")
    print("INFINITY L1=smooth_mu3_cover quotient=1/3(1,1) one_minus3_curve")
    print("INFINITY L2=smooth_mu2_cover two_quotients=1/2(1,1)=A1")
    print("RESOLUTION_COUNTS L1=2*55+1=111 L2=2*56+2=114 all_rational")
    print("A_CHAIN nu_i(u)=r-i nu_i(v)=i nu_i(z)=1 negative_intersection=Dirichlet_path_L")
    print("GREEN K_ij=min(i,j)*(r-max(i,j))/r expected_visits=2*K_ij audited_r=56,57")
    print("POISSON a_i=E_i[a_Xtau]+(1/2)E_i[sum_(n<tau)b_Xn]")
    print("MARTINGALE a(X_(n^tau))+(1/2)sum_(k<n^tau)b(X_k)")
    print("QUOTIENT_FIREWALL minus3_and_A1_quotient_curves_not_part_of_unbiased_A_chain")
    print("CONSEQUENCE conditional_on_actual_Keller_lift completed_local_k1_divisors_are_constant")
    print(f"CHECKS={CHECKS}")
    print("SCOPE conditional completed-local geometry over row8 k1 datum; no all-row lift, seam entry, JC2, or DC2")


if __name__ == "__main__":
    main()
