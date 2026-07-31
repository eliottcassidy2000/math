#!/usr/bin/env python3
"""Exact structural scout for the split degree-22 even-Faber lambda subchart.

THM-2692 closes the ``Phi_Q=0`` fibre in the inherited genuine nonsplit
degree-22 branch.  On the split quadratic deck the anti-invariant constant
line survives.  Restrict here to the even-Faber bank used by THM-2411 and
write ``Phi_Q=lambda`` on one chosen sheet.  This companion reconstructs that
bank and verifies

    q*N1 = -7496192*lambda,              N2 = 0,

together with the invariantized square/fourth-power forms.  After the
THM-2636 weighted scaling, the nonzero-lambda branch becomes

    f2(t,v,zeta)=0,       zeta*f1(t,v,zeta)^4=eta*t^23.

The fixed fibre is ``G3(v)*L5(v)^4``.  Exact transversality checks give five
old branches with ``ord(t)=4`` and three new branches with ``ord(t)=1``.
Moreover the putative fourth-root extension is globally trivial on this
subchart: the chosen-sheet coordinate is the rational function

    q = -7496192*lambda*t^5/(rho^5*f1),       q^4=Z.

Thus the nonsplit Kummer genus gate does not survive the nonzero-lambda
deformation even though the prime-23 fixed fibre remains ramified over ``t``.

The full split deck also permits odd Faber seeds; they are deliberately set
to zero here.  This is a VERIFIED-EXACT structural scout, not a closure of the
full split deck, degree twenty two, or JC(2).
"""

from __future__ import annotations

import hashlib
import importlib.util
from pathlib import Path

import sympy as s


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
BASE_PATH = HERE / "jc2_degree22_first_flux_pole_divisor_thm2411.py"
BASE_SHA256 = "7d14a16aab791db2da9dc2749117db6bdfe539096fe8d7670a9c852d44d956e3"
require(
    hashlib.sha256(BASE_PATH.read_bytes()).hexdigest() == BASE_SHA256,
    "audited THM-2411 dependency changed",
)

spec = importlib.util.spec_from_file_location("thm2411_split_lambda_base", BASE_PATH)
require(spec is not None and spec.loader is not None, "cannot load THM-2411")
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)


def normalized_primitive(expression: s.Expr, variable: s.Symbol) -> s.Expr:
    """Primitive part with positive leading coefficient."""
    polynomial = s.Poly(s.expand(expression), variable, domain=s.QQ)
    _, primitive = s.primitive(polynomial.as_expr(), variable)
    primitive_poly = s.Poly(primitive, variable, domain=s.QQ)
    if primitive_poly.LC() < 0:
        primitive = -primitive
    return s.expand(primitive)


def digest(expression: s.Expr, *variables: s.Symbol) -> str:
    polynomial = s.Poly(s.expand(expression), *variables, domain=s.QQ)
    return hashlib.sha256(s.srepr(polynomial.as_expr()).encode()).hexdigest()


# Depressed-quartic/Faber variables and fixed degree-22 coefficients.
x, d0, q, s0, T = s.symbols("x d q s T")
B, C, D, E, W = s.symbols("B C D E W")
y, u, Z = s.symbols("y u Z")
lam = s.symbols("lambda")
degrees = (2, 6, 10, 14, 18, 22)

bank = {
    degree: base.faber_observables(degree, 2 * d0, q, d0**2 - s0)
    for degree in degrees
}
weights = {22: 1, 14: B, 10: C, 6: D, 2: E}

# The first flux is anti-invariant and therefore divisible by q.  The second
# flux is invariant; W is the translated constant Psi_0.
phi_total = s.factor(sum(weights[n] * bank[n][0] for n in weights))
require(s.rem(s.Poly(phi_total, q), s.Poly(q, q)) == 0, "Phi_Q lost q parity")
phi_over_q = base.replace_even_q(s.factor(phi_total / q), q, T)
phi_base = base.replace_t_square(
    s.factor(phi_over_q.subs({d0: u / T, s0: y / 11})), T, Z
)

psi_total = s.factor(sum(weights[n] * bank[n][1] for n in weights) - W)
psi_even = base.replace_even_q(psi_total, q, T)
psi_base = base.replace_t_square(
    s.factor(psi_even.subs({d0: u / T, s0: y / 11})), T, Z
)

A = 616 * B - 1089 * u + 63 * y**2
K = (
    -745360 * B * u * y
    + 6160 * B * y**3
    + 2342560 * C * u
    - 58080 * C * y**2
    + 511104 * D * y
    - 3748096 * E
    + 922383 * u**2 * y
    - 25410 * u * y**3
    + 63 * y**5
)
N1 = 1331 * A * Z + 4 * K
N2 = (
    15944049 * Z**2
    + 65591680 * B * Z * y
    - 206145280 * C * Z
    - 162339408 * Z * u * y
    + 2236080 * Z * y**3
    + 1443016960 * B * u**2
    - 71554560 * B * u * y**2
    + 98560 * B * y**4
    + 449771520 * C * u * y
    - 1239040 * C * y**3
    - 1978994688 * D * u
    + 16355328 * D * y**2
    - 239878144 * E * y
    - 1319329792 * W
    - 1190488992 * u**3
    + 147581280 * u**2 * y**2
    - 1219680 * u * y**4
    + 672 * y**6
)
DEN = s.Integer(7496192)
require(s.factor(phi_base + N1 / DEN) == 0, "first-flux N1 normalization changed")
require(
    s.factor(psi_base - N2 / s.Integer(1319329792)) == 0,
    "second-flux N2 normalization changed",
)

# Exact ideal-membership certificates for forgetting first q and then T.
X = s.symbols("X")
square_certificate = s.expand(
    (q * X + DEN * lam) * (q * X - DEN * lam)
    - X**2 * (q**2 - T)
    - (T * X**2 - DEN**2 * lam**2)
)
fourth_certificate = s.expand(
    (T * X**2 - DEN**2 * lam**2) * (T * X**2 + DEN**2 * lam**2)
    - X**4 * (T**2 - Z)
    - (Z * X**4 - DEN**4 * lam**4)
)
require(square_certificate == 0, "q-sign invariantization certificate failed")
require(fourth_certificate == 0, "T-sign invariantization certificate failed")

# Weighted scaling from fixed coefficients to the THM-2636 signed scale.
tau, rho, v, zeta = s.symbols("t rho v zeta")
c, dd, e, w = s.symbols("c d e w")
scale = {
    B: rho**2,
    C: c * rho**3,
    D: dd * rho**4,
    E: e * rho**5,
    W: w * rho**6,
    y: rho / tau,
    u: rho**2 * v / tau**2,
    Z: rho**3 * zeta / tau**3,
}
f1 = s.factor(s.cancel(N1.subs(scale) * tau**5 / rho**5))
f2 = s.factor(s.cancel(N2.subs(scale) * tau**6 / rho**6))
require(
    not any(symbol in f1.free_symbols or symbol in f2.free_symbols for symbol in (rho,)),
    "weighted scaling retained rho",
)

# The fourth-power invariant equation becomes zeta*f1^4=eta*t^23 with
# eta=DEN^4*lambda^4/rho^23.  Verify the exponent directly before rho is
# absorbed into the nonzero constant eta.
scaled_left = s.factor(s.cancel((Z * N1**4).subs(scale)))
require(
    s.factor(scaled_left - rho**23 * zeta * f1**4 / tau**23) == 0,
    "prime-23 weighted relation changed",
)

# The chosen split sheet is not a further Kummer cover of the invariantized
# curve.  Indeed q=-DEN*lambda/N1, and in the scaled coordinates this is the
# displayed rational function.  The following identity proves q^4=Z modulo
# the invariantized curve equation, without taking a root or choosing a
# numerical branch.
eta = DEN**4 * lam**4 / rho**23
split_curve_equation = zeta * f1**4 - eta * tau**23
q_reconstructed = -DEN * lam * tau**5 / (rho**5 * f1)
physical_Z = rho**3 * zeta / tau**3
require(
    s.factor(
        q_reconstructed**4
        - physical_Z
        + rho**3 * split_curve_equation / (f1**4 * tau**3)
    )
    == 0,
    "global rational q reconstruction failed",
)

# Fixed section and its exact 3+5 branch decomposition.
f10 = s.factor(f1.subs(tau, 0))
f20 = s.factor(f2.subs(tau, 0))
G3 = normalized_primitive(f20.subs(zeta, 0), v)
G3_EXPECTED = s.expand((121 * v - 1) * (14641 * v**2 - 1694 * v + 1))
require(G3 == G3_EXPECTED, "new split fixed cubic changed")
require(s.gcd(s.Poly(G3, v), s.Poly(s.diff(G3, v), v)).degree() == 0,
        "split fixed cubic is not squarefree")

L5 = (
    155624547606 * v**5
    + 3215383215 * v**4
    - 1700698560 * v**3
    + 58124770 * v**2
    - 855470 * v
    + 2583
)
require(s.gcd(s.Poly(L5, v), s.Poly(s.diff(L5, v), v)).degree() == 0,
        "old fixed quintic is not squarefree")
require(s.gcd(s.Poly(G3, v), s.Poly(L5, v)).degree() == 0,
        "old and new fixed branches meet")

old_resultant = normalized_primitive(s.resultant(f10, f20, zeta), v)
require(old_resultant == L5, "old fixed resultant is not L5")
split_resultant = s.resultant(f20, zeta * f10**4, zeta)
split_product = s.resultant(f20, zeta, zeta) * s.resultant(f20, f10, zeta) ** 4
require(s.expand(split_resultant - split_product) == 0,
        "fixed resultant multiplicativity failed")
split_primitive = normalized_primitive(split_resultant, v)
require(split_primitive == s.expand(G3 * L5**4),
        "fixed split fibre is not G3*L5^4")
require(s.Poly(split_primitive, v).degree() == 23, "fixed fibre degree is not 23")

# Exact local controls.  On L5, f1 is a transverse coordinate and zeta is a
# unit.  On G3, f1 and d(f2|zeta=0)/dv are units.  These are precisely the
# hypotheses for the local models t^23~f1^4 and zeta~t^23.
f1_poly_zeta = s.Poly(f10, zeta)
require(f1_poly_zeta.degree() == 1, "fixed first flux is not linear in zeta")
f1_A = f1_poly_zeta.coeff_monomial(zeta)
f1_B = f1_poly_zeta.coeff_monomial(1)
zeta_solution = -f1_B / f1_A
jacobian = s.det(
    s.Matrix(
        [
            [s.diff(f10, v), s.diff(f10, zeta)],
            [s.diff(f20, v), s.diff(f20, zeta)],
        ]
    )
)
jacobian_numerator = s.together(jacobian.subs(zeta, zeta_solution)).as_numer_denom()[0]
require(s.gcd(s.Poly(L5, v), s.Poly(f1_A, v)).degree() == 0,
        "old fixed zeta reconstruction denominator vanished")
require(s.gcd(s.Poly(L5, v), s.Poly(f1_B, v)).degree() == 0,
        "zeta vanished on an old fixed point")
require(s.gcd(s.Poly(L5, v), s.Poly(jacobian_numerator, v)).degree() == 0,
        "old fixed intersections are not transverse")

new_f1 = s.Poly(f10.subs(zeta, 0), v)
new_f2_derivative = s.Poly(s.diff(f20.subs(zeta, 0), v), v)
require(s.gcd(s.Poly(G3, v), new_f1).degree() == 0,
        "f1 vanished on a new fixed point")
require(s.gcd(s.Poly(G3, v), new_f2_derivative).degree() == 0,
        "new fixed cubic is not transverse in v")

# Integral C2 augmentation control behind the cheapest order-raising test.
Delta = s.Matrix([[-1, 1], [1, -1]])
for exponent in range(1, 9):
    require(Delta**exponent == (-2) ** (exponent - 1) * Delta,
            "quadratic augmentation law changed")

print("degree-22 split even-Faber first-flux lambda / prime-23 scout")
print(f"base_thm2411_sha256={BASE_SHA256}")
print("faber_degrees=2,6,10,14,18,22")
print("first_flux_exact=Phi_Q=q*(-N1/7496192)")
print("second_flux_exact=Psi_Q-Psi0=N2/1319329792")
print("chosen_split_sheet=q*N1=-7496192*lambda")
print("q_sign_quotient=T*N1^2=7496192^2*lambda^2")
print("T_sign_quotient=Z*N1^4=7496192^4*lambda^4")
print("lambda_zero_with_q_unit=N1=0:recovers_THM2692_flux_fibre")
print("scaled_split_system=f2=0;zeta*f1^4=eta*t^23")
print("eta=7496192^4*lambda^4/rho^23:nonzero_on_split_lambda_branch")
print("chosen_sheet_q=-7496192*lambda*t^5/(rho^5*f1):rational_on_curve")
print("global_fourth_root=q^4=rho^3*zeta/t^3:Kummer_extension_trivial")
print(f"scaled_f1_digest={digest(f1, tau, v, zeta, c, dd, e, w)}")
print(f"scaled_f2_digest={digest(f2, tau, v, zeta, c, dd, e, w)}")
print("fixed_new_G3=(121v-1)(14641v^2-1694v+1)")
print("fixed_new_roots=1/121,(7+4sqrt(3))/121,(7-4sqrt(3))/121")
print("fixed_old_L5_degree=5:squarefree=True")
print("fixed_G3_degree=3:squarefree=True:coprime_to_L5=True")
print("fixed_split_eliminant=G3*L5^4:degree=23")
print("old_local_controls=zeta_unit=True:intersection_transverse=True")
print("new_local_controls=f1_unit=True:v_transverse=True")
print("normalization_t_orders=4,4,4,4,4,1,1,1")
print("t_zero_divisor_degree=5*4+3*1=23")
print("physical_Z_orders_at_t0=old:-12,new:20:both_divisible_by_4")
print("q_orders_at_t0=old:-3,new:5")
print("augmentation_law=Delta^j=(-2)^(j-1)Delta:j=1..8_checked")
print("cheapest_order_raising_test=degree6_odd_seed_boundary_Phi_Psi_on_I^j/I^(j+1)")
print("loss_q_to_T=forgets_split_sheet_sign")
print("loss_T_to_Z=forgets_square_root_T_and_enlarges_to_fourth-power_shadow")
print("loss_even_bank=forgets_full_split_odd_Faber_seed_coefficients")
print("scope=VERIFIED_EXACT_EVEN_FABER_SUBCHART_NOT_FULL_SPLIT_NOT_JC2")
print("ALL CHECKS PASSED")
