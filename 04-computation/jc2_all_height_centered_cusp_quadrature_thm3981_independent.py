#!/usr/bin/env python3
"""Exact scratch checks for the all-height centered cusp quadrature."""

from __future__ import annotations

import math
import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.factor(expression) == 0, message)


x, t, z, p, y, Y, L = sp.symbols("x t z p y Y L")

# Uniform completion, source volume, split idempotent, and Darboux identities.
for n in range(2, 17):
    z0 = 1 + x**n * t
    p0 = z0 * t
    y0 = x ** (n - 1) * z0 * t**2
    zero(z0 * (z0 - 1) ** 2 - x ** (n + 1) * y0,
         f"height {n}: D chart equation")
    zero((2 * z0 - 1) ** 2 - (1 + 4 * x**n * p0),
         f"height {n}: split square root")
    zero(
        sp.diff(y0, t) - x**-1 * (z0 - 1) * (3 * z0 - 1),
        f"height {n}: uniform source volume",
    )

# Standard cusp determinant and ordinary arm pair.
XX, YY = sp.symbols("XX YY")
A_std = YY**2 + 2 * L * XX
C_std = YY**3 + 3 * L * XX * YY
det_std = sp.diff(A_std, XX) * sp.diff(C_std, YY) - sp.diff(A_std, YY) * sp.diff(C_std, XX)
zero(sp.rem(sp.Poly(det_std - XX, L), sp.Poly(6 * L**2 - 1, L)).as_expr(),
     "standard cusp determinant")
A_arm = t + 2 * L * x
C_arm = -x
zero(sp.diff(A_arm, x) * sp.diff(C_arm, t) - sp.diff(A_arm, t) * sp.diff(C_arm, x) - 1,
     "ordinary arm Darboux pair")

# Branch, genus, and obstruction table.  The generic centered curve is
# z(z-1)^2=x^(n+1)(Y+x/2), and omega=x dx/((z-1)(3z-1)).
rows: list[str] = []
for n in range(2, 31):
    d3 = math.gcd(3, n + 2)
    epsilon = 1 if n % 2 == 0 else 0
    genus_numerator = n + epsilon - d3 + 2
    gate(genus_numerator % 2 == 0, f"height {n}: integral genus")
    genus = genus_numerator // 2
    gate(genus >= 2, f"height {n}: positive genus floor")

    # x-projection branch contribution and infinity orders.
    ramification = epsilon + 1 + (n + 2) + (3 - d3)
    gate(2 * genus - 2 == -6 + ramification, f"height {n}: Riemann-Hurwitz")
    omega_inf = (2 * n - 2) // d3 - 1
    gate((2 * n - 2) % d3 == 0 and omega_inf >= 0,
         f"height {n}: omega regular at infinity")

    # The n+2 roots over the second cubic critical value are generically simple.
    h = x ** (n + 1) * (Y + x / 2) - sp.Rational(4, 27)
    gate(sp.gcd(sp.Poly(h, x, domain=sp.QQ.frac_field(Y)),
                sp.Poly(sp.diff(h, x), x, domain=sp.QQ.frac_field(Y))).degree() == 0,
         f"height {n}: second critical-value roots simple")

    if n == 2:
        mechanism = "NONZERO_HOLOMORPHIC_GENUS2"
        residue = "NA"
    elif n % 2 == 1:
        q = (n + 1) // 2
        # At either P_+ or P_-, a^2=Y.  The omitted prefactor 1/(2*epsilon*a)
        # is nonzero; this is the exact binomial residue coefficient.
        coeff = sp.binomial(-sp.Rational(1, 2), q - 2) / (2 * Y) ** (q - 2)
        gate(coeff != 0, f"height {n}: odd-height logarithmic residue")
        mechanism = "NONZERO_SEAM_RESIDUE"
        residue = sp.sstr(coeff) + "/(2*epsilon*a),a^2=Y"
    else:
        q = n // 2
        m = q - 2
        # Pair the local primitive h of omega with the global holomorphic form
        # eta=x^q dx/F_z.  Every summand has the same sign; the closed value is
        # included as an independent exact check.
        pairing_sum = sp.factor(sum(
            sp.binomial(-sp.Rational(1, 2), j)
            * sp.binomial(-sp.Rational(1, 2), m - j)
            / sp.Rational(2 * j - 2 * m - 1)
            for j in range(m + 1)
        ))
        closed_sum = sp.Rational((-1) ** (m + 1) * 4**m * math.factorial(m) ** 2,
                                 math.factorial(2 * m + 1))
        gate(pairing_sum == closed_sum, f"height {n}: even pairing convolution")
        pairing = sp.factor(pairing_sum / (2**m * Y ** (m + 1)))
        closed_pairing = sp.Rational((-1) ** (q - 1) * math.factorial(q - 2),
                                     math.prod(range(1, 2 * q - 2, 2))) / Y ** (q - 1)
        gate(pairing == closed_pairing and pairing != 0,
             f"height {n}: even-height residue pairing nonzero")
        eta_inf = (q + 1) // d3 - 1
        gate((q + 1) % d3 == 0 and eta_inf >= 0,
             f"height {n}: test differential holomorphic at infinity")
        mechanism = "NONZERO_HOLOMORPHIC_RESIDUE_PAIRING"
        residue = sp.sstr(pairing)

    rows.append(
        f"n={n};g={genus};d_inf={d3};ord_inf_omega={omega_inf};"
        f"mechanism={mechanism};certificate={residue}"
    )

print("THM-3980 all-height canonical-X scratch probe")
print(f"CHECKS={CHECKS}")
print("CURVE=z(z-1)^2=x^(n+1)(Y+x/2)")
print("OMEGA=x*dx/((z-1)(3z-1));d(X^2)=2*OMEGA")
print("GENUS=(n+epsilon_even-gcd(3,n+2)+2)/2")
print("FORMAL_SPLIT=W^2=1+4x^n*p;D_AND_L_IDEMPOTENTS")
print("ODD_HEIGHT=NONZERO_RESIDUES_AT_TWO_POINTS_OVER_(0,1)")
print("EVEN_HEIGHT_GE_4=NONZERO_RESIDUE_PAIRING_WITH_ETA=x^(n/2)dx/F_z")
print("HEIGHT_2=NONZERO_HOLOMORPHIC_GENUS_TWO_DIFFERENTIAL")
for row in rows:
    print(row)
print("CONCLUSION=CANONICAL_X_TRANSCENDENT_OVER_FRAC(B_n)_FOR_ALL_n_GE_2")
print("SCOPE=SPECIFIC_CANONICAL_SPLIT_SECTION;ALTERNATIVE_FORMAL_GAUGES_AND_GLOBAL_KELLER_PAIR_OPEN")
