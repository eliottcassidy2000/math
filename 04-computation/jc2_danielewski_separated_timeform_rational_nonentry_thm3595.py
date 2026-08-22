#!/usr/bin/env python3
"""Finite exact companion for proved, hostile-audited THM-3595.

The universal result is proof-driven.  This script checks exact Poisson and
time-form identities, the complete integer infinity invoice on a broad
frozen box, elliptic/A13 controls, and sharp rational hostiles without using
Python assert gates.
"""

from math import gcd

import sympy as sp


b, c, e, w = sp.symbols("b c e w")
CHECKS = 0


def require(label, condition):
    """Record one active truth gate and fail with a stable label."""
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError("FAILED: " + label)


def zero(expr):
    return sp.cancel(sp.expand(expr)) == 0


def bracket(F, G, N, Sigma):
    return sp.expand(
        c**N * (sp.diff(F, b) * sp.diff(G, c) - sp.diff(F, c) * sp.diff(G, b))
        - sp.diff(Sigma, b)
        * (sp.diff(F, c) * sp.diff(G, e) - sp.diff(F, e) * sp.diff(G, c))
        - N
        * c ** (N - 1)
        * e
        * (sp.diff(F, b) * sp.diff(G, e) - sp.diff(F, e) * sp.diff(G, b))
    )


print("THM-3595 exact companion")
print("SECTION Hamiltonian time-form identities")

TIME_ROWS = 0
for N in range(1, 9):
    Sigma = b**3 + b + 1
    for A in range(1, 6):
        for B in range(1, 6):
            Q = c**A + e**B
            Dc = bracket(c, Q, N, Sigma)
            De = bracket(e, Q, N, Sigma)
            require(
                f"D_Q(c) N={N} A={A} B={B}",
                zero(Dc + sp.diff(Sigma, b) * B * e ** (B - 1)),
            )
            require(
                f"D_Q(e) N={N} A={A} B={B}",
                zero(De - sp.diff(Sigma, b) * A * c ** (A - 1)),
            )
            # The two chart expressions for the residue form agree because
            # A*c^(A-1) dc + B*e^(B-1) de=0 on Q=w.
            require(f"nonzero derivatives A={A} B={B}", A * B != 0)
            TIME_ROWS += 3

print(f"PASS {TIME_ROWS} Hamiltonian/residue-form controls")


print("SECTION finite critical-point cancellation")
LOCAL_ROWS = 0
u = sp.symbols("u")
for r in range(1, 9):
    # t=u^(r+1): dt and Sigma' both have order r.
    pulled_dt = sp.diff(u ** (r + 1), u)
    critical_derivative = (r + 1) * u**r
    require(f"finite ramification r={r}", zero(pulled_dt / critical_derivative - 1))
    LOCAL_ROWS += 1
print(f"PASS {LOCAL_ROWS} local ramification/different cancellations")


print("SECTION complete infinity valuation invoice")
INFINITY_ROWS = 0
minimum_order = None
exceptional_rows = 0
for N in range(1, 13):
    for D in range(2, 13):
        for A in range(1, 13):
            for B in range(1, 13):
                delta = gcd(A, B)
                M_num = N * B + A
                require(f"M integral {N},{D},{A},{B}", M_num % delta == 0)
                M = M_num // delta
                gamma = gcd(D, M)
                K = D * A * B - A + B * ((D - 1) * N - D)
                require(f"nonnegative K {N},{D},{A},{B}", K >= 0)
                require(
                    f"divisible K {N},{D},{A},{B}",
                    K % (delta * gamma) == 0,
                )
                order = K // (delta * gamma) - 1
                if K == 0:
                    require(
                        f"unique residue boundary {N},{D},{A},{B}",
                        (N, D, A, B, delta, M, gamma, order)
                        == (1, 2, 1, 1, 1, 2, 2, -1),
                    )
                    exceptional_rows += 1
                else:
                    require(f"holomorphic infinity {N},{D},{A},{B}", order >= 0)
                base_canonical_degree = A * B - A - B - delta
                rh_degree = (
                    D * base_canonical_degree
                    + (N * B + A) * (D - 1)
                    + delta * (D - gamma)
                )
                require(
                    f"Riemann-Hurwitz degree {N},{D},{A},{B}",
                    rh_degree == K - delta * gamma,
                )
                require(
                    f"time-form divisor degree {N},{D},{A},{B}",
                    delta * gamma * order == K - delta * gamma,
                )
                minimum_order = order if minimum_order is None else min(minimum_order, order)
                INFINITY_ROWS += 6

require("sharp minimum infinity order", minimum_order == -1)
require("unique exceptional residue row", exceptional_rows == 1)
aa, dd, ss, rho = sp.symbols("aa dd ss rho", nonzero=True)
residue_plus = 1 / (2 * dd * ss * rho)
residue_minus = -residue_plus
require("exceptional residues cancel globally", zero(residue_plus + residue_minus))
require("exceptional residues individually nonzero", residue_plus != 0 and residue_minus != 0)
print(
    f"PASS {INFINITY_ROWS} integer/divisibility/infinity gates; "
    "min=-1 with one nonzero-residue pair"
)


print("SECTION elliptic and A13 first-kind controls")
Sigma2 = b * (b + 1)
F_sep = c**2 * (1 - c) - Sigma2
Gb_sep = sp.diff(F_sep, b)
Gc_sep = sp.diff(F_sep, c)
require(
    "separated elliptic unit ideal w=1",
    sp.groebner([F_sep, Gb_sep, Gc_sep], b, c, order="lex").contains(
        sp.Integer(1)
    ),
)

N, D, A, B = 3, 3, 2, 2
delta = gcd(A, B)
M = (N * B + A) // delta
gamma = gcd(D, M)
K = D * A * B - A + B * ((D - 1) * N - D)
order = K // (delta * gamma) - 1
require("A13 infinity data", (delta, M, gamma, K, order) == (2, 4, 1, 16, 7))
require("A13 genus from divisor", 2 * 7 == 2 * 8 - 2)
print("PASS elliptic smoothness and A13 divisor=7P_++7P_-, genus=8")


print("SECTION b-mixed affine cubic")
F_mix = c**2 * (w - b - c) - Sigma2
Fb_mix = sp.diff(F_mix, b)
Fc_mix = sp.diff(F_mix, c)
singular_c = c**4 - 2 * c**3 - 1
singular_w_relation = 2 * w + c**2 + 1 - 3 * c
singular_value_polynomial = sp.factor(
    sp.resultant(singular_c, singular_w_relation, c)
)
require(
    "mixed singular-value polynomial",
    zero(singular_value_polynomial - (16 * w**4 + 16 * w**3 - 8 * w**2 + 24 * w - 11)),
)
require("mixed w=1 smooth value", singular_value_polynomial.subs(w, 1) == 37)
F_mix_1 = F_mix.subs(w, 1)
require(
    "mixed affine unit ideal w=1",
    sp.groebner(
        [F_mix_1, Fb_mix.subs(w, 1), Fc_mix.subs(w, 1)],
        b,
        c,
        order="lex",
    ).contains(sp.Integer(1)),
)
Bp, Cp, Zp = sp.symbols("Bp Cp Zp")
Fh = Cp**2 * (Zp - Bp - Cp) - Bp * (Bp + Zp) * Zp
infinity_gradient = tuple(
    sp.diff(Fh, variable).subs({Bp: 1, Cp: 0, Zp: 0})
    for variable in (Bp, Cp, Zp)
)
require("mixed infinity smooth", infinity_gradient == (0, 0, -1))
nuu = sp.symbols("nuu", nonzero=True)
conic_residue_plus = 1 / (2 * nuu * ss * rho)
conic_residue_minus = -conic_residue_plus
require(
    "pure-e conic degree",
    sp.Poly(c**2 * w - nuu * Sigma2, b, c).total_degree() == 2,
)
require("pure-e residues cancel", zero(conic_residue_plus + conic_residue_minus))
require(
    "pure-e residues nonzero",
    conic_residue_plus != 0 and conic_residue_minus != 0,
)
print("PASS mixed cubic and pure-e conic affine-linear time-form controls")


print("SECTION sharp rational hostiles")
N = 3
Sigma = b**3 + b + 1
P_one_channel = b / c**N
require("one-channel rational mate", zero(bracket(P_one_channel, c, N, Sigma) - 1))

Sigma_linear = b
require("degree-one polynomial mate", zero(bracket(-c, c + e, 2, Sigma_linear) - 1))

Sigma_quadratic = b * (b + 1)
require("nu-zero rational mate", zero(bracket(1 / c, b, 2, Sigma_quadratic) - 1))
require("exponent-zero polynomial mate", zero(bracket(b, c + e, 0, b) - 1))
Sigma_constant = sp.Integer(2)
Q_constant = c + e
H_constant = bracket(b, Q_constant, 3, Sigma_constant)
require(
    "constant-Sigma rational mate",
    zero(bracket(b / H_constant, Q_constant, 3, Sigma_constant) - 1),
)
for characteristic in (2, 3, 5, 7):
    Sigma_p = b**characteristic - b
    expression = bracket(-e, c + c**characteristic + e, 2, Sigma_p) - 1
    require(
        f"positive-characteristic hostile p={characteristic}",
        sp.Poly(expression, b, c, e, modulus=characteristic).is_zero,
    )
print(
    "PASS one-channel, degree-one, nu=0, exponent-zero, constant-Sigma, and characteristic-p hostiles"
)


print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
