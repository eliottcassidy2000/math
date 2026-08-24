#!/usr/bin/env python3
"""Exact primary certificate for THM-4012.

The theorem is a conditional stable-specialization observer, not a claim that
the weighted leading face alone determines the stable model.  This script
checks the algebraic part of that observer: the common weighted degeneration,
the target good model, every singleton face through weight 18, the exact
weight-six elliptic survivor, its exact boundary normalization and two-prime
attachment no-go, and the self-contained weight-eight Bolza split and
CM-field mismatch.

Universe: rational function rings of characteristic zero, with the explicitly
named nonzero coefficients inverted.  Algebraic roots used for changes of
variables are represented by their exact defining relations.  Stable-model
admissibility and conservation of degree are geometric hypotheses/proofs in
THM-4012 and are not inferred from these symbolic identities.
"""

from __future__ import annotations

import sys
from math import gcd, isqrt

import sympy as sp


sys.stdout.reconfigure(newline="\n")
CHECKS = 0


def simp(expr):
    return sp.factor(sp.cancel(sp.expand(expr)))


def require(condition, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError(label)
    print(f"PASS  {label}")


def equal(left, right, label: str) -> None:
    require(simp(left - right) == 0, label)


def rho_exponents(expr, rho):
    """Return the rho exponent of each expanded nonzero term."""
    return [term.as_powers_dict().get(rho, 0) for term in sp.Add.make_args(sp.expand(expr))]


def ec_add(left, right, prime: int):
    """Add points on y^2=x^3+1 over F_prime; None is infinity."""
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2 and (y1 + y2) % prime == 0:
        return None
    if left == right:
        slope = (3 * x1 * x1) * pow(2 * y1, -1, prime) % prime
    else:
        slope = (y2 - y1) * pow(x2 - x1, -1, prime) % prime
    x3 = (slope * slope - x1 - x2) % prime
    y3 = (slope * (x1 - x3) - y1) % prime
    return x3, y3


def ec_mul(multiplier: int, point, prime: int):
    result = None
    addend = point
    n = multiplier
    while n:
        if n & 1:
            result = ec_add(result, addend, prime)
        addend = ec_add(addend, addend, prime)
        n >>= 1
    return result


def valuation(integer: int, prime: int) -> int:
    value = abs(integer)
    answer = 0
    while value and value % prime == 0:
        answer += 1
        value //= prime
    return answer


print("STATUS=THM-4012;PROVED_CONDITIONAL_OBSERVER;EXACT_MAX6_NO_GO;VERIFIED-EXACT_PRIMARY")
print("SCOPE=general_face_stability_required;exact_max6_face_stability_proved;JC(2)_OPEN")

# -------------------------------------------------------------------------
# I. Exact weighted degeneration and simultaneous target good reduction.
# -------------------------------------------------------------------------

rho, S, P, Z = sp.symbols("rho S P Z")
p, y, s, q = sp.symbols("p y s q")
gamma = sp.symbols("gamma", nonzero=True)

for M in range(2, 13):
    pairs = [
        (i, j)
        for i in range(M // 2 + 1)
        for j in range(M // 3 + 1)
        if 0 < 2 * i + 3 * j <= M
    ]
    top_pairs = [(i, j) for i, j in pairs if 2 * i + 3 * j == M]
    require(bool(top_pairs), f"weight M={M} has a nonempty monomial face")
    coeff = {(i, j): sp.Symbol(f"c_{M}_{i}_{j}", nonzero=True) for i, j in pairs}
    H = sum(coeff[i, j] * p**i * y**j for i, j in pairs)
    H_top = sum(coeff[i, j] * P**i * (S * P) ** j for i, j in top_pairs)
    H_rho = sp.expand(
        rho ** (6 * M)
        * H.subs({p: rho**-12 * P, y: rho**-18 * S * P})
    )
    require(
        min(rho_exponents(H_rho, rho)) >= 0,
        f"weighted H_rho is integral M={M}",
    )
    equal(H_rho.subs(rho, 0), H_top, f"weighted top face survives M={M}")

    source = s**2 * (q + gamma - H) - p * (q - H)
    source_scaled = sp.expand(
        rho ** (6 * M + 12)
        * source.subs(
            {
                q: rho ** (-6 * M),
                s: rho**-6 * S,
                p: rho**-12 * P,
                y: rho**-18 * S * P,
            }
        )
    )
    expected = (S**2 - P) * (1 - H_rho) + gamma * rho ** (6 * M) * S**2
    equal(source_scaled, expected, f"source weighted scaling M={M}")
    equal(
        source_scaled.subs(rho, 0),
        (S**2 - P) * (1 - H_top),
        f"source central factorization M={M}",
    )

    H_bar = sum(
        coeff[i, j]
        * rho ** (6 * (M - 2 * i - 3 * j))
        * S**j
        * P ** (i + j)
        * Z ** (M - 2 * i - 3 * j)
        for i, j in pairs
    )
    equal(H_bar.subs(Z, 1), H_rho, f"weighted homogenization chart M={M}")
    weighted_degrees = []
    for term in sp.Add.make_args(sp.expand(H_bar)):
        powers = term.as_powers_dict()
        weighted_degrees.append(powers.get(S, 0) + 2 * powers.get(P, 0) + powers.get(Z, 0))
    require(set(weighted_degrees) == {M}, f"all homogenized terms have weighted degree M={M}")
    special_closure = sp.expand((S**2 - P) * (Z**M - H_bar.subs(rho, 0)))
    expected_closure = (S**2 - P) * (Z**M - H_top)
    equal(special_closure, expected_closure, f"proper weighted special divisor M={M}")

A, C, X, Y, a = sp.symbols("A C X Y a", nonzero=True)
target = C**2 - A**3 + sp.Rational(3, 4) * a**2 * A - q + a**3 / 4
for M in range(2, 13):
    target_scaled = sp.expand(
        rho ** (6 * M)
        * target.subs(
            {
                q: rho ** (-6 * M),
                A: rho ** (-2 * M) * X,
                C: rho ** (-3 * M) * Y,
            }
        )
    )
    expected = (
        Y**2
        - X**3
        - 1
        + sp.Rational(3, 4) * a**2 * rho ** (4 * M) * X
        + a**3 * rho ** (6 * M) / 4
    )
    equal(target_scaled, expected, f"target simultaneous good scaling M={M}")
    equal(target_scaled.subs(rho, 0), Y**2 - X**3 - 1, f"target j=0 special fibre M={M}")

print("RESULT central_support=(S^2-P)*(Z^M-H_M(P,SP))")
print("RESULT target_special=Y^2-X^3-1 (smooth j=0)")
print("CONTROL stable_specialization_is_a_named_geometric_hypothesis_not_a_symbolic_conclusion")

# -------------------------------------------------------------------------
# II. All singleton top faces are rational at face level.
# -------------------------------------------------------------------------

c = sp.symbols("c", nonzero=True)
singleton_count = 0
for M in range(2, 19):
    for i in range(M // 2 + 1):
        for j in range(M // 3 + 1):
            if 2 * i + 3 * j != M:
                continue
            singleton_count += 1
            exponent_S = j
            exponent_P = i + j
            components = gcd(exponent_S, exponent_P)
            f0 = S**2 - P
            f1 = 1 - c * S**exponent_S * P**exponent_P
            determinant = sp.diff(f0, S) * sp.diff(f1, P) - sp.diff(f0, P) * sp.diff(f1, S)
            intersection = sp.Poly(c * S**M - 1, S, domain=sp.QQ.frac_field(c))
            require(
                components >= 1
                and gcd(M, gcd(exponent_S, exponent_P)) == components
                and simp(determinant.subs(P, S**2) + M * c * S ** (M - 1)) == 0
                and intersection.degree() == M
                and sp.gcd(intersection, intersection.diff()).degree() == 0,
                f"singleton rational-component/attachment ledger ({i},{j})",
            )
            print(
                f"SINGLETON weight={M} monomial=p^{i}y^{j} "
                f"rational_components={components} attachments={M}"
            )

require(singleton_count == 36, "complete singleton census through weight 18")
print("RESULT singleton_face_normalizations=all_P1;observer_Hom=0_under_stable_admissibility")

# -------------------------------------------------------------------------
# III. Weight six is the exact j=0 survivor boundary.
# -------------------------------------------------------------------------

epsilon, kappa = sp.symbols("epsilon kappa", nonzero=True)
T = sp.symbols("T")
H6 = epsilon * P**3 + kappa * S**2 * P**2
equal(
    (S**2 - P) * (1 - H6),
    S**2 * (1 - H6) - P * (1 - H6),
    "weight-six central factorization",
)
equal(
    (1 - H6).subs(S**2, T**2 / P**2),
    1 - epsilon * P**3 - kappa * T**2,
    "weight-six face becomes elliptic under T=SP",
)

u6, v6 = sp.symbols("u6 v6", nonzero=True)
face6 = kappa * T**2 + epsilon * P**3 - 1
equal(
    face6.subs({epsilon: -u6**3, kappa: v6**2, P: X / u6, T: Y / v6}),
    Y**2 - X**3 - 1,
    "weight-six elliptic face is the target j=0 curve",
)

f0 = S**2 - P
f6 = 1 - H6
det6 = sp.diff(f0, S) * sp.diff(f6, P) - sp.diff(f0, P) * sp.diff(f6, S)
equal(
    det6.subs(P, S**2),
    -6 * (epsilon + kappa) * S**5,
    "weight-six six-attachment transversality",
)

A5 = a**5
gamma_live = -a**3 / 2
epsilon_live = -sp.Rational(1376, 135) / a**12
epsilon_tilde_live = sp.Rational(2752, 135) / A5**3
e_tilde_p4_zero = -sp.Rational(5696, 45) / A5**3
kappa_live = simp(gamma_live * e_tilde_p4_zero)
equal(gamma_live * epsilon_tilde_live, epsilon_live, "live normalized/raw p^3 conversion")
equal(kappa_live, sp.Rational(2848, 45) / a**12, "live p4-zero branch raw y^2 coefficient")
equal(
    epsilon_live + kappa_live,
    sp.Rational(7168, 135) / a**12,
    "live p4-zero branch is nonresonant",
)
equal(
    epsilon_tilde_live + e_tilde_p4_zero,
    -sp.Rational(14336, 135) / A5**3,
    "live normalized p4-zero branch is nonresonant",
)
print("RESULT weight6=elliptic_j0_face;Hom_to_target_nonzero;observer_survives_not_suffices")
print("RESULT live_max6_Rtilde_face=(2752*p^3/135-5696*y^2/45)/A5^3;sum=-14336/(135*A5^3)")
print("RESULT live_max6_raw_face=gamma*Rtilde_face;raw_sum=7168/(135*a^12)")

# -------------------------------------------------------------------------
# IIIb. Exact max-weight-six boundary normalization and face stability.
# -------------------------------------------------------------------------

rbd, zbd, wbd = sp.symbols("rbd zbd wbd")
lambda2, coeff_y3, coeff_p4, coeff_py5 = sp.symbols(
    "lambda2 coeff_y3 coeff_p4 coeff_py5"
)
H6_boundary_chart = (
    epsilon * rbd**3
    + kappa * rbd**2
    + coeff_py5 * rho**6 * rbd**2 * zbd
    + coeff_p4 * rho**12 * rbd**2 * zbd**2
    + coeff_y3 * rho**18 * rbd * zbd**3
    + lambda2 * rho**24 * rbd * zbd**4
)
F6_boundary = sp.expand(
    (1 - rbd) * (zbd**6 - H6_boundary_chart)
    + gamma * rho**36 * zbd**6
)
equal(
    F6_boundary.subs(rho, 0),
    (1 - rbd) * (zbd**6 - epsilon * rbd**3 - kappa * rbd**2),
    "max-weight-six S=1 boundary central equation",
)
F6_homogeneous_central = (S**2 - P) * (Z**6 - epsilon * P**3 - kappa * S**2 * P**2)
equal(
    F6_homogeneous_central.subs({S: 0, P: 1, Z: 0}),
    epsilon,
    "max-weight-six closure value at [0:1:0] is epsilon",
)
require(epsilon != 0, "max-weight-six closure avoids the [0:1:0] orbifold point")
equal(
    sp.diff(F6_boundary.subs(rho, 0), rbd).subs({rbd: 1, zbd: 0}),
    epsilon + kappa,
    "r=1 boundary point is smooth off resonance",
)
r_infinity = -kappa / epsilon
equal(
    sp.diff(F6_boundary.subs(rho, 0), rbd).subs({rbd: r_infinity, zbd: 0}),
    -kappa**2 * (epsilon + kappa) / epsilon**2,
    "r=-kappa/epsilon elliptic infinity is smooth off resonance",
)

Phi6 = simp(F6_boundary.subs(rbd, zbd**3 * wbd) / zbd**6)
Phi6_expected = sp.expand(
    (1 - zbd**3 * wbd)
    * (
        1
        - kappa * wbd**2
        - coeff_y3 * rho**18 * wbd
        - zbd * (coeff_py5 * rho**6 * wbd**2 + lambda2 * rho**24 * wbd)
        - coeff_p4 * rho**12 * zbd**2 * wbd**2
        - epsilon * zbd**3 * wbd**3
    )
    + gamma * rho**36
)
equal(Phi6, Phi6_expected, "persistent boundary A5 strict-transform equation")
equal(
    Phi6.subs(zbd, 0),
    1 + gamma * rho**36 - coeff_y3 * rho**18 * wbd - kappa * wbd**2,
    "persistent boundary quadratic at z=0",
)
hensel_discriminant = sp.discriminant(sp.Poly(Phi6.subs(zbd, 0), wbd), wbd)
equal(
    hensel_discriminant,
    coeff_y3**2 * rho**36 + 4 * kappa * (1 + gamma * rho**36),
    "persistent boundary Hensel discriminant",
)
equal(hensel_discriminant.subs(rho, 0), 4 * kappa, "boundary Hensel discriminant is a unit")
print("RESULT exact_nonresonant_max6_boundary=two_smooth_Hensel_branches;no_positive_genus_tail")
print("RESULT exact_nonresonant_max6_face_stability=unique_elliptic_plus_six_UV_rho36_chains")

# -------------------------------------------------------------------------
# IV. Weight eight is the Bolza genus-two face.
# -------------------------------------------------------------------------

c8, d8 = sp.symbols("c8 d8", nonzero=True)
H8 = c8 * P**4 + d8 * S**2 * P**3
T8, Y8 = sp.symbols("T8 Y8")
equal(
    (1 - H8).subs(S**2, T8**2 / P**2),
    1 - c8 * P**4 - d8 * P * T8**2,
    "weight-eight first birational coordinate T=SP",
)
equal(
    P * (1 - c8 * P**4 - d8 * P * T8**2).subs(T8, Y8 / P),
    P * (1 - c8 * P**4) - d8 * Y8**2,
    "weight-eight hyperelliptic equation dY^2=P(1-cP^4)",
)

u8, v8, xb, yb = sp.symbols("u8 v8 xb yb", nonzero=True)
weight8_hyperelliptic = d8 * Y8**2 - P * (1 - c8 * P**4)
equal(
    u8
    * weight8_hyperelliptic.subs(
        {c8: -u8**4, d8: v8**2 / u8, P: xb / u8, Y8: yb / v8}
    ),
    yb**2 - xb * (xb**4 + 1),
    "weight-eight face scales to the Bolza curve",
)

f8 = 1 - H8
det8 = sp.diff(f0, S) * sp.diff(f8, P) - sp.diff(f0, P) * sp.diff(f8, S)
equal(
    det8.subs(P, S**2),
    -8 * (c8 + d8) * S**7,
    "weight-eight eight-attachment transversality",
)

bolza = yb**2 - xb * (xb**4 + 1)
require(sp.gcd(sp.Poly(xb * (xb**4 + 1), xb), sp.Poly(sp.diff(xb * (xb**4 + 1), xb), xb)).degree() == 0,
        "Bolza quintic is squarefree")
require(sp.Poly(xb * (xb**4 + 1), xb).degree() == 5, "Bolza smooth projective genus is two")

# -------------------------------------------------------------------------
# V. Self-contained Bolza quotient split Jac(C) ~ E^2.
# -------------------------------------------------------------------------

inv1_curve = bolza.subs({xb: 1 / xb, yb: yb / xb**3}, simultaneous=True)
equal(xb**6 * inv1_curve, bolza, "Bolza involution i1 preserves the curve")

uq = xb + 1 / xb
vq = yb * (xb + 1) / xb**2
equal(
    vq**2 - (uq + 2) * (uq**2 - 2),
    (xb + 1) ** 2 * bolza / xb**4,
    "i1 quotient equation v^2=(u+2)(u^2-2)",
)

a2, a4, a6 = sp.Integer(2), sp.Integer(-2), sp.Integer(-4)
b2 = 4 * a2
b4 = 2 * a4
b6 = 4 * a6
b8 = 4 * a2 * a6 - a4**2
c4 = b2**2 - 24 * b4
delta = -b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
jE = sp.cancel(c4**3 / delta)
require((b2, b4, b6, b8) == (8, -4, -16, -36), "quotient Weierstrass b-invariants")
require(c4 == 160 and delta == 512 and jE == 8000, "quotient elliptic invariants c4=160 Delta=512 j=8000")

inv2_curve = bolza.subs({xb: -1 / xb, yb: sp.I * yb / xb**3}, simultaneous=True)
equal(xb**6 * inv2_curve, -bolza, "Bolza involution i2 preserves the curve")

scale = -sp.I
equal(scale**2, -1, "conjugating Bolza scaling sends x-involution to -1/x")
equal(scale**3, sp.I, "conjugating Bolza scaling sends y-involution to i*y/x^3")

M_i1 = sp.Matrix([[0, -1], [-1, 0]])
M_i2 = sp.Matrix([[0, sp.I], [-sp.I, 0]])
fixed1 = sp.Matrix([1, -1])
fixed2 = sp.Matrix([1, -sp.I])
require(M_i1 * fixed1 == fixed1, "i1 fixes differential omega0-omega1")
require(M_i2 * fixed2 == fixed2, "i2 fixes differential omega0-i*omega1")
require(sp.det(sp.Matrix.hstack(fixed1, fixed2)) == 1 - sp.I, "quotient differentials span H0(C,Omega)")
print("RESULT Bolza_Jacobian_isogenous_to=E8000^2 via two independent involution quotients")

# -------------------------------------------------------------------------
# VI. Self-contained CM mismatch E8000 versus target E0.
# -------------------------------------------------------------------------

xe, ye = sp.symbols("xe ye", nonzero=True)
fE = xe**3 - 4 * xe**2 + 2 * xe
equal(
    (Y**2 - (X + 2) * (X**2 - 2)).subs({X: xe - 2, Y: ye}),
    ye**2 - fE,
    "shift quotient E to y^2=x^3-4x^2+2x",
)

Xphi = xe - 4 + 2 / xe
Yphi = ye * (2 - xe**2) / xe**2
Eprime_relation = Yphi**2 - (Xphi**3 + 8 * Xphi**2 + 8 * Xphi)
equal(
    Eprime_relation.subs(ye**2, fE),
    0,
    "explicit degree-two quotient E to Eprime",
)

Xp, Yp = sp.symbols("Xp Yp")
Eprime = Yp**2 - (Xp**3 + 8 * Xp**2 + 8 * Xp)
scaled_back = -sp.Rational(1, 8) * Eprime
target_after_scale = (-sp.Rational(1, 8) * Yp**2) - (
    (-Xp / 2) ** 3 - 4 * (-Xp / 2) ** 2 + 2 * (-Xp / 2)
)
equal(target_after_scale, scaled_back, "Eprime isomorphism x=-X/2, y=tau*Y, tau^2=-1/8")

alpha_x = -Xphi / 2
alpha_x_twice = simp(-(alpha_x - 4 + 2 / alpha_x) / 2)
slope_yfree = (3 * xe**2 - 8 * xe + 2) * ye / (2 * fE)
double_x = simp((3 * xe**2 - 8 * xe + 2) ** 2 / (4 * fE) + 4 - 2 * xe)
double_y = simp(-ye + slope_yfree * (xe - double_x))
alpha_y_twice = simp(
    -sp.Rational(1, 8)
    * ye
    * (2 - xe**2)
    / xe**2
    * (2 - alpha_x**2)
    / alpha_x**2
)
equal(alpha_x_twice, double_x, "degree-two self-map alpha squared has x([2])")
equal(alpha_y_twice, -double_y, "degree-two self-map alpha^2=[-2]")

zeta = sp.symbols("zeta")
require(sp.rem(zeta**3 - 1, zeta**2 + zeta + 1, domain=sp.QQ) == 0,
        "target order-three automorphism has CM polynomial zeta^2+zeta+1")
require(
    isqrt(2) ** 2 != 2 and isqrt(3) ** 2 != 3,
    "2/3 is not a rational square, so CM fields -2 and -3 differ",
)
print("RESULT End0(E8000)=Q(sqrt(-2)) from alpha^2=[-2]")
print("RESULT End0(E0)=Q(sqrt(-3)) from the order-three automorphism")
print("RESULT Hom(E8000,E0)=0 and therefore Hom(Jac(Bolza),E0)=0")

# -------------------------------------------------------------------------
# VII. The unique live weight-six attachment is provably nontorsion.
# -------------------------------------------------------------------------

attachment_x_cube = simp(-epsilon_tilde_live / (epsilon_tilde_live + e_tilde_p4_zero))
attachment_y_square = simp(e_tilde_p4_zero / (epsilon_tilde_live + e_tilde_p4_zero))
equal(attachment_x_cube, sp.Rational(43, 224), "live attachment X^3=43/224")
equal(attachment_y_square, sp.Rational(267, 224), "live attachment Y^2=267/224")
equal(attachment_y_square, attachment_x_cube + 1, "live attachment lies on Y^2=X^3+1")

reduction_controls = [(11, (2, 3), 6), (17, (7, 2), 9)]
for prime, point, order in reduction_controls:
    xred, yred = point
    require((xred**3 - int(attachment_x_cube.p) * pow(int(attachment_x_cube.q), -1, prime)) % prime == 0,
            f"attachment cubic root reduces to x={xred} mod {prime}")
    require((yred**2 - int(attachment_y_square.p) * pow(int(attachment_y_square.q), -1, prime)) % prime == 0,
            f"attachment square root reduces to y={yred} mod {prime}")
    require((yred**2 - xred**3 - 1) % prime == 0, f"reduced point lies on E mod {prime}")
    require((-432) % prime != 0, f"E has good reduction at residue characteristic {prime}")
    require(ec_mul(order, point, prime) is None, f"reduced attachment killed by {order} mod {prime}")
    for proper_prime_divisor in sorted(set(sp.factorint(order))):
        require(
            ec_mul(order // proper_prime_divisor, point, prime) is not None,
            f"order {order} mod {prime} survives divisor {proper_prime_divisor}",
        )
    print(
        f"REDUCTION p={prime} point={point} exact_order={order} "
        f"multiples={[ec_mul(n, point, prime) for n in range(1, order + 1)]}"
    )

require(
    valuation(6, 2) == 1
    and valuation(11, 2) == 0
    and valuation(9, 2) == 0
    and valuation(17, 2) == 0,
    "no N can satisfy N/6 an 11-power and N/9 a 17-power",
)
print("RESULT live_attachment_is_nontorsion_by_good_reduction_orders_6_and_9")
print("RESULT exact_max6_six_attachment_invoice_excludes_the_forced_survivor_unconditionally")

# -------------------------------------------------------------------------
# VIII. Live degree-at-most-eight truncation ledger.
# -------------------------------------------------------------------------

r40_tilde, e_tilde = sp.symbols("r40_tilde e_tilde")
row_coupling = r40_tilde + sp.Rational(6, 7) * e_tilde / A5 + sp.Rational(11392, 105) / A5**4
equal(
    row_coupling.subs({r40_tilde: 0, e_tilde: e_tilde_p4_zero}),
    0,
    "THM-4007 p4-zero branch satisfies the exact affine lock",
)
require(
    simp(row_coupling.subs({r40_tilde: 0, e_tilde: 0})) != 0,
    "THM-4007 forbids simultaneous p4=y2=0",
)
print("RESULT b=d=0_exact_max6_candidate_is_nontorsion_excluded")
print("RESULT face_stable_b=d=0_conclusion=total_wtdeg(H)>=9")
print("CONTROL higher_weight_terms_change_the_top_face_and_are_not_erased")

print(f"CHECKS={CHECKS}")
print("THEOREM_ID=THM-4012")
print("NOTE=Hom-factor obstruction is necessary only; face admissibility is explicit")
print("ALL THM-4012 EXACT CHECKS PASSED")
