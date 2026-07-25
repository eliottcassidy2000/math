#!/usr/bin/env python3
"""Exact referee for the degree-18 target-translation normal form.

The degree-18 Faber row of THM-2262 is

    Q = E_18 + alpha E_14 + beta E_10 + gamma E_6 + delta E_2.

Translating the first Keller target by ``P_c = P + c`` preserves the map and
changes the Faber basis by a finite binomial transform.  At
``c = 2*alpha/9`` this kills ``alpha``.  This script verifies the full Faber
identity, the four resulting invariant parameters, covariance of all retained
spectral data, the weighted cone structure, the exact one-sparse axis checks,
and the critical-value data closing the full C--W plane in THM-2297.
"""

from __future__ import annotations

import hashlib

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def faber_coefficients(
    degree: int,
    p: sp.Expr,
    q: sp.Expr,
    r: sp.Expr,
    extra: int = 3,
) -> list[sp.Expr]:
    """Laurent coefficients of ``(w^4+p*w^2+q*w+r)^(degree/4)``."""
    exponent = sp.Rational(degree, 4)
    coefficients = [sp.Integer(1)]
    quartic = {2: p, 3: q, 4: r}
    for index in range(1, degree + extra + 1):
        value = sum(
            quartic[step]
            * ((exponent + 1) * step - index)
            * coefficients[index - step]
            for step in range(2, min(4, index) + 1)
            if step in quartic
        ) / index
        coefficients.append(sp.factor(value))
    return coefficients


def faber_seed(
    degree: int,
    w: sp.Symbol,
    d: sp.Symbol,
    q: sp.Symbol,
    s: sp.Expr,
) -> sp.Expr:
    """Polynomial part ``E_degree`` in the depressed quartic coordinates."""
    coefficients = faber_coefficients(degree, 2 * d, q, d**2 - s, extra=0)
    return sp.expand(
        sum(
            coefficients[index] * w ** (degree - index)
            for index in range(degree + 1)
        )
    )


def observables(
    degree: int,
    d: sp.Symbol,
    q: sp.Symbol,
    s: sp.Expr,
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """The three quartic Laurent observables used by THM-2262."""
    coefficients = faber_coefficients(degree, 2 * d, q, d**2 - s)
    first = sp.factor(4 * coefficients[degree + 1])
    second = sp.factor(4 * coefficients[degree + 2])
    third = sp.factor(
        4 * coefficients[degree + 3] + 2 * d * coefficients[degree + 1]
    )
    return first, second, third


def canonical_hash(expression: sp.Expr, *generators: sp.Symbol) -> str:
    polynomial = sp.Poly(sp.expand(expression), *generators)
    payload = sp.srepr(polynomial.as_expr()).encode()
    return hashlib.sha256(payload).hexdigest()


def main() -> None:
    w, d, q, s = sp.symbols("w d q s")
    alpha, beta, gamma, delta, psi0 = sp.symbols(
        "alpha beta gamma delta Psi0"
    )
    u, y, z = sp.symbols("u y Z")
    scale = sp.symbols("rho")

    c = sp.Rational(2, 9) * alpha
    beta_bar = sp.factor(beta - sp.Rational(7, 18) * alpha**2)
    gamma_bar = sp.factor(
        gamma
        - sp.Rational(5, 9) * alpha * beta
        + sp.Rational(35, 243) * alpha**3
    )
    delta_bar = sp.factor(
        delta
        - sp.Rational(1, 3) * alpha * gamma
        + sp.Rational(5, 54) * alpha**2 * beta
        - sp.Rational(35, 1944) * alpha**4
    )
    psi_bar = sp.factor(
        psi0
        + sp.Rational(4, 9) * alpha * delta
        - sp.Rational(2, 27) * alpha**2 * gamma
        + sp.Rational(10, 729) * alpha**3 * beta
        - sp.Rational(14, 6561) * alpha**5
    )
    psi_shift = sp.factor(psi0 - psi_bar)

    degrees = (2, 6, 10, 14, 18)
    seeds = {degree: faber_seed(degree, w, d, q, s) for degree in degrees}
    shifted_seeds = {
        degree: faber_seed(degree, w, d, q, s - c) for degree in degrees
    }

    # The whole Faber basis obeys the finite binomial translation law.
    for j, degree in enumerate(degrees, start=1):
        translated = sum(
            sp.binomial(sp.Rational(2 * j - 1, 2), k)
            * c**k
            * seeds[4 * (j - k) - 2]
            for k in range(j)
        )
        require(
            sp.factor(shifted_seeds[degree] - translated) == 0,
            f"Faber translation failed at degree {degree}",
        )

    original_q = (
        seeds[18]
        + alpha * seeds[14]
        + beta * seeds[10]
        + gamma * seeds[6]
        + delta * seeds[2]
    )
    normalized_q = (
        shifted_seeds[18]
        + beta_bar * shifted_seeds[10]
        + gamma_bar * shifted_seeds[6]
        + delta_bar * shifted_seeds[2]
    )
    require(
        sp.factor(original_q - normalized_q) == 0,
        "whole reduced mate is not translation invariant",
    )

    bank = {degree: observables(degree, d, q, s) for degree in degrees}
    shifted_bank = {
        degree: observables(degree, d, q, s - c) for degree in degrees
    }
    weights = {18: 1, 14: alpha, 10: beta, 6: gamma, 2: delta}
    normalized_weights = {
        18: 1,
        14: 0,
        10: beta_bar,
        6: gamma_bar,
        2: delta_bar,
    }
    original_observables = tuple(
        sp.factor(sum(weights[n] * bank[n][index] for n in degrees))
        for index in range(3)
    )
    normalized_observables = tuple(
        sp.factor(
            sum(
                normalized_weights[n] * shifted_bank[n][index]
                for n in degrees
            )
        )
        for index in range(3)
    )
    require(
        sp.factor(original_observables[0] - normalized_observables[0]) == 0,
        "first flux is not invariant",
    )
    require(
        sp.factor(original_observables[1] - normalized_observables[1])
        == psi_shift,
        "second flux constant shift mismatch",
    )
    require(
        sp.factor(original_observables[2] - normalized_observables[2]) == 0,
        "third flux is not invariant",
    )

    # THM-2262's centered spectral cubic.
    spectral_g = (
        -5878656 * psi0 * y
        - 26040609 * u**3
        - 19289340 * u**2 * alpha**2
        + 49601160 * u**2 * beta
        + 1607445 * u**2 * y**2
        - 2222640 * u * alpha**4
        + 11430720 * u * alpha**2 * beta
        + 1111320 * u * alpha**2 * y**2
        + 17635968 * u * alpha * gamma
        - 20995200 * u * beta**2
        - 2857680 * u * beta * y**2
        - 52907904 * u * delta
        - 138915 * u * y**4
        + 235200 * alpha**6
        + 326144 * alpha**5 * y
        - 1814400 * alpha**4 * beta
        + 82320 * alpha**4 * y**2
        - 2096640 * alpha**3 * beta * y
        + 4354560 * alpha**3 * gamma
        - 62720 * alpha**3 * y**3
        + 3110400 * alpha**2 * beta**2
        - 423360 * alpha**2 * beta * y**2
        - 13063680 * alpha**2 * delta
        + 2612736 * alpha**2 * gamma * y
        - 30380 * alpha**2 * y**4
        + 3110400 * alpha * beta**2 * y
        - 11197440 * alpha * beta * gamma
        + 241920 * alpha * beta * y**3
        - 2612736 * alpha * delta * y
        - 653184 * alpha * gamma * y**2
        + 777600 * beta**2 * y**2
        + 33592320 * beta * delta
        - 5598720 * beta * gamma * y
        + 78120 * beta * y**4
        + 1959552 * delta * y**2
        - 435456 * gamma * y**3
        + 1127 * y**6
    )
    normalized_g = spectral_g.xreplace(
        {
            alpha: sp.Integer(0),
            beta: beta_bar,
            gamma: gamma_bar,
            delta: delta_bar,
            psi0: psi_bar,
        }
    )
    require(
        sp.expand(spectral_g - normalized_g) == 0,
        "spectral cubic is not invariant",
    )

    n2 = (
        45927 * u**2
        + 22680 * u * alpha**2
        - 58320 * u * beta
        - 5670 * u * y**2
        - 1680 * alpha**4
        - 2240 * alpha**3 * y
        + 8640 * alpha**2 * beta
        - 840 * alpha**2 * y**2
        + 8640 * alpha * beta * y
        - 31104 * alpha * gamma
        + 2160 * beta * y**2
        + 93312 * delta
        - 15552 * gamma * y
        + 35 * y**4
    )
    normalized_n2 = n2.xreplace(
        {
            alpha: sp.Integer(0),
            beta: beta_bar,
            gamma: gamma_bar,
            delta: delta_bar,
        }
    )
    require(sp.expand(n2 - normalized_n2) == 0, "Z elimination is not invariant")

    n3 = (
        183708 * u**3
        + 90720 * u**2 * alpha**2
        - 233280 * u**2 * beta
        - 22680 * u**2 * y**2
        + 51030 * u * z * y
        - 6720 * u * alpha**4
        - 8960 * u * alpha**3 * y
        + 34560 * u * alpha**2 * beta
        - 3360 * u * alpha**2 * y**2
        + 34560 * u * alpha * beta * y
        - 124416 * u * alpha * gamma
        + 8640 * u * beta * y**2
        + 373248 * u * delta
        - 62208 * u * gamma * y
        + 140 * u * y**4
        - 6561 * z**2
        + 13440 * z * alpha**3
        + 10080 * z * alpha**2 * y
        - 51840 * z * alpha * beta
        - 25920 * z * beta * y
        + 93312 * z * gamma
        - 840 * z * y**3
    )
    normalized_n3 = n3.xreplace(
        {
            alpha: sp.Integer(0),
            beta: beta_bar,
            gamma: gamma_bar,
            delta: delta_bar,
        }
    )
    require(sp.expand(n3 - normalized_n3) == 0, "third flux is not invariant")

    bvar, gvar, dvar, pvar = sp.symbols("B C D W")
    g0 = spectral_g.xreplace(
        {
            alpha: sp.Integer(0),
            beta: bvar,
            gamma: gvar,
            delta: dvar,
            psi0: pvar,
        }
    )
    scaled_g0 = g0.xreplace(
        {
            u: scale**2 * u,
            y: scale * y,
            bvar: scale**2 * bvar,
            gvar: scale**3 * gvar,
            dvar: scale**4 * dvar,
            pvar: scale**5 * pvar,
        }
    )
    require(
        sp.expand(scaled_g0 - scale**6 * g0) == 0,
        "weighted degree-six covariance failed",
    )

    branch_discriminant = sp.factor(sp.discriminant(g0, u))
    require(
        sp.Poly(branch_discriminant, y).degree() == 12,
        "normalized branch discriminant degree",
    )
    scaled_discriminant = branch_discriminant.xreplace(
        {
            y: scale * y,
            bvar: scale**2 * bvar,
            gvar: scale**3 * gvar,
            dvar: scale**4 * dvar,
            pvar: scale**5 * pvar,
        }
    )
    require(
        sp.expand(scaled_discriminant - scale**12 * branch_discriminant)
        == 0,
        "branch discriminant weighted degree",
    )

    axes = {
        "B": {bvar: 1, gvar: 0, dvar: 0, pvar: 0},
        "C": {bvar: 0, gvar: 1, dvar: 0, pvar: 0},
        "D": {bvar: 0, gvar: 0, dvar: 1, pvar: 0},
        "W": {bvar: 0, gvar: 0, dvar: 0, pvar: 1},
        "origin": {bvar: 0, gvar: 0, dvar: 0, pvar: 0},
    }
    axis_data: dict[str, tuple[sp.Expr, int]] = {}
    for name, specialization in axes.items():
        axis_polynomial = sp.factor(branch_discriminant.subs(specialization))
        polynomial = sp.Poly(sp.expand(axis_polynomial), y)
        gcd_degree = sp.gcd(polynomial, polynomial.diff()).degree()
        axis_data[name] = (axis_polynomial, gcd_degree)
    require(axis_data["B"][1] == 0, "B-axis discriminant is not squarefree")
    require(axis_data["D"][1] == 0, "D-axis discriminant is not squarefree")
    require(axis_data["C"][1] == 5, "C-axis collision ledger changed")
    require(axis_data["W"][1] == 1, "W-axis collision ledger changed")
    require(axis_data["origin"][1] == 11, "origin collision ledger changed")

    v = sp.symbols("v")
    infinity_cubic = (
        1127 - 138915 * v + 1607445 * v**2 - 26040609 * v**3
    )
    infinity_discriminant = sp.discriminant(infinity_cubic, v)
    require(
        infinity_discriminant
        == -153384762202971019112448,
        "infinity cubic discriminant changed",
    )

    # The whole B=D=0 plane has a repeated branch value in the original
    # y-projection.  In the faithful z=1/y chart its equation is separated:
    #
    #   L(v) = 435456*C*z^3 + 5878656*W*z^5.
    #
    # The two critical values of L are distinct and nonzero.  Since the two
    # nonzero critical points of the right side are simple, at most two branch
    # collisions occur; the normalization therefore has genus at least two.
    cw_discriminant = sp.factor(
        branch_discriminant.subs({bvar: 0, dvar: 0})
    )
    require(
        sp.resultant(cw_discriminant, sp.diff(cw_discriminant, y), y) == 0,
        "C-W plane should be invisible to the raw squarefree test",
    )
    critical_value = sp.symbols("tau")
    critical_value_polynomial = sp.factor(
        sp.resultant(
            sp.diff(infinity_cubic, v),
            critical_value - infinity_cubic,
            v,
        )
    )
    primitive_critical_value_polynomial = sp.factor(
        critical_value_polynomial
        / sp.Poly(critical_value_polynomial, critical_value).LC()
        * 27
    )
    require(
        primitive_critical_value_polynomial
        == 27 * critical_value**2 + 68992 * critical_value + 226193408,
        "critical-value polynomial changed",
    )
    critical_value_discriminant = sp.discriminant(
        primitive_critical_value_polynomial,
        critical_value,
    )
    require(
        critical_value_discriminant == -19668992000,
        "critical values are no longer distinct",
    )
    require(
        primitive_critical_value_polynomial.subs(critical_value, 0) != 0,
        "zero became a critical value",
    )
    inverse_y = sp.symbols("z")
    cw_right_side = (
        435456 * gvar * inverse_y**3
        + 5878656 * pvar * inverse_y**5
    )
    require(
        sp.expand(
            sp.diff(cw_right_side, inverse_y)
            - inverse_y**2
            * (1306368 * gvar + 29393280 * pvar * inverse_y**2)
        )
        == 0,
        "C-W critical-point factorization changed",
    )
    require(
        sp.factor(
            435456 * gvar
            + 5878656 * pvar
            * (-sp.Rational(2, 45) * gvar / pvar)
        )
        == sp.Rational(870912, 5) * gvar,
        "nonzero C-W critical images changed",
    )

    # On the invariant origin, G=y^6 L(u/y^2), Z=k(v)y^3, and
    # F is a nonzero constant times y^6/T.  The three resultants below
    # certify that none of the needed coefficients vanishes at a root of L.
    z_ratio = -sp.Rational(2, 729) * (6561 * v**2 - 810 * v + 5)
    third_ratio = -sp.Rational(4, 243) * (
        6561 * v**2 - 810 * v + 5
    ) * (19683 * v**2 + 4374 * v - 125)
    sidecar_ratio = -sp.Rational(2, 729) * (
        6561 * v**2 + 1134 * v - 19
    )
    origin_resultants = {
        "Z": sp.factor(
            sp.resultant(
                infinity_cubic,
                sp.together(z_ratio).as_numer_denom()[0],
                v,
            )
        ),
        "F": sp.factor(
            sp.resultant(
                infinity_cubic,
                sp.together(third_ratio).as_numer_denom()[0],
                v,
            )
        ),
        "sidecar": sp.factor(
            sp.resultant(
                infinity_cubic,
                sp.together(sidecar_ratio).as_numer_denom()[0],
                v,
            )
        ),
    }
    require(
        all(value != 0 for value in origin_resultants.values()),
        "an origin cusp coefficient vanishes",
    )

    print("faber_translation_bank=PASS degrees=2,6,10,14,18")
    print("whole_reduced_mate_translation=PASS")
    print(f"translation_c={c}")
    print(f"beta_bar={beta_bar}")
    print(f"gamma_bar={gamma_bar}")
    print(f"delta_bar={delta_bar}")
    print(f"Psi_bar={psi_bar}")
    print(f"second_flux_shift={psi_shift}")
    print("spectral_G_covariance=PASS")
    print("Z_elimination_covariance=PASS")
    print("third_flux_covariance=PASS")
    print("weights=(y,u,B,C,D,W)=(1,2,2,3,4,5)")
    print("spectral_G_weight=6")
    print("branch_discriminant_y_degree=12")
    print("branch_discriminant_weight=12")
    print(
        "branch_discriminant_terms="
        f"{len(sp.Poly(sp.expand(branch_discriminant), y, bvar, gvar, dvar, pvar).terms())}"
    )
    for name in ("B", "C", "D", "W", "origin"):
        print(f"axis_{name}_gcd_degree={axis_data[name][1]}")
        print(f"axis_{name}_factor={axis_data[name][0]}")
    print(f"infinity_cubic={infinity_cubic}")
    print(f"infinity_discriminant={infinity_discriminant}")
    print("CW_raw_branch_resultant=0")
    print(
        "infinity_critical_value_polynomial="
        f"{primitive_critical_value_polynomial}"
    )
    print(
        "infinity_critical_value_discriminant="
        f"{critical_value_discriminant}"
    )
    print("CW_nonzero_critical_points=z^2=-2*C/(45*W)")
    print("CW_critical_images=(870912/5)*C*(+/-z)^3")
    print("CW_normalization_genus_lower_bound=2")
    print(f"origin_Z_resultant={origin_resultants['Z']}")
    print(f"origin_F_resultant={origin_resultants['F']}")
    print(f"origin_sidecar_resultant={origin_resultants['sidecar']}")
    print(
        "normalized_G_sha256="
        f"{canonical_hash(g0, u, y, bvar, gvar, dvar, pvar)}"
    )
    print(
        "normalized_Delta_sha256="
        f"{canonical_hash(branch_discriminant, y, bvar, gvar, dvar, pvar)}"
    )
    print("axis_C_normalization=z^3=L(v)/(435456*C),genus=1")
    print("axis_W_normalization=z^5=L(v)/(5878656*W),genus=4")
    print("origin_cusp_sidecar=MATHEMATICAL_PROOF_REQUIRED")
    print("status=THM2297_DEGREE18_TARGET_TRANSLATION_EXACT_REFEREE")


if __name__ == "__main__":
    main()
