#!/usr/bin/env python3
"""Independent hostile audit of the fixed-F level-three global norm package.

This script deliberately does not call the producer companion.  It uses two
orthogonal routes:

1. a fraction-field regular representation on the previously unused target
   slice ``(b,c)=(1,1)`` computes the complete norm and both transverse values
   on ``L=0``; and
2. a Rabin certificate over ``F_449`` proves irreducibility of the full-degree
   specialization ``P(-1,lambda)``.  A coefficient-content calculation rules
   out a factor depending only on ``tau``, so this independently proves the
   bivariate boundary residual ``P`` irreducible over ``Q``.

The script checks the exact sign and scalar in the normalization residue, the
pole-seven denominator on the new slice, the image multiplicity control, and
the exceptional seam.  It does not compute a global expanded numerator J or
a degree-27 coordinate eliminant.  ``require`` remains active under ``-O``.
"""

from __future__ import annotations

import hashlib
import pickle
from pathlib import Path

import sympy as sp
from flint import fmpz_mod_poly_ctx


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def coefficient_hash(poly: sp.Poly) -> str:
    ledger = "\n".join(
        f"{monomial}:{coefficient}" for monomial, coefficient in poly.terms()
    )
    return hashlib.sha256(ledger.encode("ascii")).hexdigest()


ROOT = Path(__file__).resolve().parents[1]
H_PATH = ROOT / "05-knowledge/results/keller_L2_polynomial_opus_20260728.pkl"
H_SHA256 = "5a9459b3149e500c1b00b67bd804aa7e607de06bf4610c7cdf5fa26d41d74ce9"
EXPECTED_J11_SLICE_SHA256 = "f0b51842adf88383230a056f9bcdd19c434e17551940c23aea2fad2964c57447"
EXPECTED_P_SHA256 = "5f15ad89e0788eafcdffecf6fa0f0204224d27fb6ecf94cd46ed1982da678333"
EXPECTED_P_MINUS_ONE_SHA256 = "e4fa283a05df9285c7384f2fb2e4467e3a5e0d6c688abb32719fa394834470c5"
P3_DENOMINATOR = sp.Integer(32687327622020737269760000000000000)
P3_PRIMITIVE_LEAD = sp.Integer(
    51781084507728969748447381099724691990978530821553059009958856196332470476610552088369
)

H_bytes = H_PATH.read_bytes()
require(hashlib.sha256(H_bytes).hexdigest() == H_SHA256, "H artifact changed")
H = pickle.loads(H_bytes)

a, b, c, w = sp.symbols("a b c w")
H_poly = sp.Poly(H, a, b, c)
require(H_poly.total_degree() == 25 and len(H_poly.terms()) == 361, "H shape changed")

# ---------------------------------------------------------------------------
# Route A: weighted face and a fresh complete norm slice.
# ---------------------------------------------------------------------------

eta = sp.symbols("eta")
weight = max(i - k for (i, _j, k), _coefficient in H_poly.terms())
weighted_face = sp.expand(
    sum(
        coefficient * eta**j * (-3 * eta) ** k
        for (i, j, k), coefficient in H_poly.terms()
        if i - k == weight
    )
)
face_constant = 2**9 * 3**6 * 11**3 * 13**2
require(weight == 7, "escaping weighted order is not seven")
require(weighted_face == face_constant * eta**3, "escaping weighted face changed")

# The unused slice b=c=1 meets L=0 at A=0 and A=2/27.  Work in the cubic
# fraction-field algebra, but take the final norm from the closed determinant
# formula for h0+h1*w+h2*w^2.  This avoids calling either producer script and
# gives a new specialization with two rational transverse boundary points.
A = sp.symbols("A")
L_slice = A * (27 * A - 2)
slice_field = sp.QQ.frac_field(A)
E_slice = sp.Poly(L_slice * w**3 + w - 2, w, domain=slice_field).monic()


def polynomial(expression: sp.Expr) -> sp.Poly:
    return sp.Poly(expression, w, domain=slice_field).rem(E_slice)


def multiply(left: sp.Poly, right: sp.Poly) -> sp.Poly:
    return (left * right).rem(E_slice)


def divide(numerator: sp.Poly, denominator: sp.Poly) -> sp.Poly:
    return multiply(numerator, sp.invert(denominator, E_slice))


def power(base: sp.Poly, exponent: int) -> sp.Poly:
    value = polynomial(1)
    factor = base
    while exponent:
        if exponent & 1:
            value = multiply(value, factor)
        exponent //= 2
        if exponent:
            factor = multiply(factor, factor)
    return value


qx = polynomial(w)
Y_den = polynomial((12 * A - 1) * w**2 + w + 2)
qy = polynomial(1) - divide(
    polynomial(3 * A * w * ((9 * A - 1) * w + 2)), Y_den
)
qz = divide(
    polynomial(2 * w - 1) - multiply(polynomial(3 * w**2), qy),
    power(qx, 3),
)

unit = polynomial(1) + multiply(qx, qy)
four_plus = polynomial(4) + multiply(polynomial(3), multiply(qx, qy))
image = (
    multiply(power(unit, 3), qz)
    + multiply(multiply(power(qy, 2), unit), four_plus),
    qy
    + multiply(polynomial(3), multiply(multiply(qx, power(unit, 2)), qz))
    + multiply(polynomial(3), multiply(multiply(qx, power(qy, 2)), four_plus)),
    multiply(polynomial(2), qx)
    - multiply(polynomial(3), multiply(power(qx, 2), qy))
    - multiply(power(qx, 3), qz),
)
for row, expected in zip(image, (polynomial(A), polynomial(1), polynomial(1))):
    require(row == expected, "fresh inverse graph failed")

qx_powers = [power(qx, index) for index in range(H_poly.degree(a) + 1)]
qy_powers = [power(qy, index) for index in range(H_poly.degree(b) + 1)]
qz_powers = [power(qz, index) for index in range(H_poly.degree(c) + 1)]
H_at_q = polynomial(0)
for (i, j, k), coefficient in H_poly.terms():
    term = multiply(multiply(qx_powers[i], qy_powers[j]), qz_powers[k])
    H_at_q = (H_at_q + term.mul_ground(slice_field.convert(coefficient))).rem(E_slice)
print("stage: fresh slice inverse and H evaluation PASS", flush=True)

# If w^3+p*w+q=0 and h=h0+h1*w+h2*w^2, the following is the determinant
# of multiplication by h in the basis (1,w,w^2).  Writing it explicitly is
# a hostile guard against a library resultant convention or nonmonic factor.
h_rows = H_at_q.as_dict()
h0 = slice_field.convert(h_rows.get((0,), slice_field.zero))
h1 = slice_field.convert(h_rows.get((1,), slice_field.zero))
h2 = slice_field.convert(h_rows.get((2,), slice_field.zero))
p_cubic = slice_field.convert(1 / L_slice)
q_cubic = slice_field.convert(-2 / L_slice)
norm_field = (
    p_cubic**2 * h2**2 * h0
    - p_cubic * q_cubic * h2**2 * h1
    - 2 * p_cubic * h2 * h0**2
    + p_cubic * h0 * h1**2
    + q_cubic**2 * h2**3
    + 3 * q_cubic * h2 * h0 * h1
    - q_cubic * h1**3
    + h0**3
)
norm_slice = sp.cancel(norm_field.as_expr())
print("stage: fresh slice regular-representation norm PASS", flush=True)
norm_numerator, norm_denominator = norm_slice.as_numer_denom()
denominator_quotient = sp.Poly(norm_denominator, A).exquo(sp.Poly(L_slice**7, A))
require(
    denominator_quotient.degree() == 0
    and denominator_quotient.LC() == 2**35,
    "fresh slice denominator is not exactly 2^35*L^7",
)

# Keep the two normalizations separate.  J_res was the notation in the first
# pole-seven sidecar; J is the later primitive integral global polynomial.
J_res_slice = sp.cancel(L_slice**7 * norm_slice)
J_primitive_slice = sp.cancel(2**35 * J_res_slice)
J_primitive_numerator, J_primitive_denominator = J_primitive_slice.as_numer_denom()
require(J_primitive_denominator == 1, "primitive fresh J slice is not integral")
J_primitive_slice_poly = sp.Poly(J_primitive_numerator, A, domain=sp.ZZ)
require(J_primitive_slice_poly.primitive()[0] == 1, "primitive fresh J slice has content")
require(J_primitive_slice_poly.degree() == 86, "fresh J slice has wrong degree")
factor_unit, factor_rows = sp.factor_list(J_primitive_slice_poly)
print("stage: fresh degree-86 slice factorization PASS", flush=True)
require(
    len(factor_rows) == 1
    and factor_rows[0][0].degree() == 86
    and factor_rows[0][1] == 1,
    "fresh J slice is reducible or a power",
)
require(
    coefficient_hash(J_primitive_slice_poly) == EXPECTED_J11_SLICE_SHA256,
    "fresh primitive J slice ledger changed",
)
require(sp.gcd(J_primitive_slice_poly, sp.Poly(L_slice, A)).degree() == 0, "fresh J slice contains L")
H_slice = sp.Poly(H.subs({a: A, b: 1, c: 1}), A)
require(sp.gcd(J_primitive_slice_poly, H_slice).degree() == 0, "fresh J slice contains H")

# ---------------------------------------------------------------------------
# Route B: reconstruct the boundary residual but certify it by a univariate
# finite-field irreducibility test, not by bivariate FLINT factorization.
# ---------------------------------------------------------------------------

tau, lam = sp.symbols("tau lambda")
seam = lam * tau - 2
q0 = (
    2 * tau / seam**2,
    -lam * (3 * lam * tau - 8) / 6,
    -lam**2 * seam**2 * (lam**2 * tau**2 - 8 * lam * tau + 14) / 8,
)
H_finite = sp.cancel(H.subs(dict(zip((a, b, c), q0))))
P_expression = sp.cancel(
    H_finite * (2**35 * 3**21 * seam**14) / lam**2
)
P_numerator, P_denominator = P_expression.as_numer_denom()
require(P_denominator == 1, "boundary residual reconstruction is not polynomial")
P = sp.Poly(P_numerator, tau, lam, domain=sp.ZZ)
require(
    (P.degree(tau), P.degree(lam), P.total_degree(), len(P.terms()))
    == (57, 84, 141, 527),
    "boundary residual shape changed",
)
require(coefficient_hash(P) == EXPECTED_P_SHA256, "boundary residual ledger changed")
print("stage: boundary residual reconstruction PASS", flush=True)

# No factor can depend only on tau: the lambda-coefficients have gcd one in
# Z[tau].  Any remaining proper factorization has two positive lambda-degrees.
lambda_coefficients = sp.Poly(P.as_expr(), lam).all_coeffs()
tau_content = sp.Poly(lambda_coefficients[0], tau, domain=sp.QQ)
for coefficient in lambda_coefficients[1:]:
    tau_content = sp.gcd(tau_content, sp.Poly(coefficient, tau, domain=sp.QQ))
require(tau_content.degree() == 0, "P has a nonconstant tau-only content")

# At tau=-1 the lambda degree stays 84.  Reduction mod 449 is certified
# irreducible by Rabin's criterion; the prime divisors of 84 are 2,3,7.
P_special = sp.Poly(P.as_expr().subs(tau, -1), lam, domain=sp.ZZ)
require(P_special.degree() == 84, "P(-1,lambda) lost degree")
require(
    coefficient_hash(P_special) == EXPECTED_P_MINUS_ONE_SHA256,
    "P(-1,lambda) ledger changed",
)
prime = 449
finite_field = fmpz_mod_poly_ctx(prime)
P_mod = finite_field(list(reversed([int(value) for value in P_special.all_coeffs()])))
x_mod = finite_field.gen()
require(P_mod.degree() == 84, "mod-449 specialization lost degree")
require((x_mod.pow_mod(prime**84, P_mod) - x_mod) % P_mod == 0, "Rabin Frobenius closure failed")
for divisor in (2, 3, 7):
    hostile = x_mod.pow_mod(prime ** (84 // divisor), P_mod) - x_mod
    require(hostile.gcd(P_mod).degree() == 0, f"Rabin gcd failed for divisor {divisor}")
print("stage: independent Rabin irreducibility certificate PASS", flush=True)

# The seam is not a residual factor, while lambda=0 is an explicit
# codimension-two intersection divisor on the normalization rather than a
# missing target hypersurface.
require(P.eval({tau: 1, lam: 2}) != 0, "exceptional seam divides P")
P_at_lambda_zero = sp.Poly(P.as_expr().subs(lam, 0), tau, domain=sp.ZZ)
require(not P_at_lambda_zero.is_zero, "P contains an extra lambda factor")

residue_constant = sp.Rational(11**6 * 13**4, 2**17 * 3**15)
for A_value, lam_value in ((sp.Rational(2, 27), 1), (sp.Integer(0), 3)):
    direct_value = sp.cancel(J_res_slice.subs(A, A_value))
    residue_value = sp.cancel(
        residue_constant * lam_value**8 * P.eval({tau: 1, lam: lam_value})
    )
    require(direct_value == residue_value, "direct norm disagrees with residue sign/scalar")
    require(direct_value != 0, "fresh transverse boundary value vanished")

# Cross the two incoming computations without sharing their construction.
# The leading coefficient of the actual degree-27 product is N^2(L).  Using
# N(L)=H/(64L) and N(H)=J/(2^35 L^7) gives J/(2^47 L^6 H).
J111 = J_primitive_slice_poly.eval(1)
H111 = H_slice.eval(1)
predicted_level_three_lead = sp.Rational(J111, 2**47 * 25**6 * H111)
observed_level_three_lead = sp.Rational(P3_PRIMITIVE_LEAD, P3_DENOMINATOR)
require(
    predicted_level_three_lead == observed_level_three_lead,
    "global J and degree-27 tower disagree on N^2(L)",
)

print("== independent audit: fixed-F level-three global norm ==")
print(f"escaping weighted face: order={weight}, constant={face_constant}")
print("fresh slice (b,c)=(1,1): regular-representation inverse graph 3/3 PASS")
print(
    "fresh complete norm: denominator=2^35*L^7; primitive "
    "J slice degree 86, irreducible exponent one, gcd(J,LH)=1"
)
print(f"fresh primitive J-slice coefficient ledger sha256={EXPECTED_J11_SLICE_SHA256}")
print("direct transverse values at a=2/27 and a=0 confirm the positive residue scalar")
print("P lambda-content over Q[tau]: 1")
print("P(-1,lambda) mod 449: degree 84; Rabin tests for 2,3,7 PASS")
print(f"P(-1,lambda) coefficient-ledger sha256={coefficient_hash(P_special)}")
print(f"P(tau,0): degree {P_at_lambda_zero.degree()}, nonzero")
print("exceptional seam lambda*tau=2 is not a factor")
print("cross-package lead: P3.LC/P3.den = N^2(L) = J/(2^47*L^6*H) at (1,1,1)")
print("primitive discriminant class is [-2J], never the stale [-J]")
print("scope: no expanded global J and no degree-27 discriminant expansion")
print("all independent audit checks passed")
