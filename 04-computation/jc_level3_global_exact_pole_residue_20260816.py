#!/usr/bin/env python3
"""Exact global L-pole and boundary-residue probe for the fixed F norm tower.

The degree-three function-field norm is the norm attached to the fixed
sporadic Keller map F.  This companion does *not* expand its global level-three
numerator.  Instead it certifies the local calculation at the generic point of
the old Jelonek divisor L:

    v_L(Norm(H)) = -7.

Together with finiteness of F off V(L) (THM-2473), this is the decisive global
denominator statement: L^7 Norm(H) is a polynomial and is not divisible by L.
The script also computes the exact pullback of its L-boundary residue to the
THM-2576 normalization, factors the 527-term residual with python-flint, and
pins the inherited irreducible one-variable slice used to type the generic
image multiplicity as one.

No degree-27 eliminant, Jacobian-counterexample classification, planar-JC, or
higher-family claim is made.  ``require`` remains active under ``python -O``.
"""

from __future__ import annotations

import hashlib
import pickle
from pathlib import Path

import sympy as sp
from flint import ctx as flint_ctx
from flint import fmpz_mpoly_ctx


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
H_PATH = ROOT / "05-knowledge/results/keller_L2_polynomial_opus_20260728.pkl"
SLICE_SCRIPT = ROOT / "04-computation/keller_level_three_norm_slice_probe_20260815.py"
SLICE_OUTPUT = ROOT / "05-knowledge/results/keller_level_three_norm_slice_probe_20260815.out"

H_ARTIFACT_SHA256 = "5a9459b3149e500c1b00b67bd804aa7e607de06bf4610c7cdf5fa26d41d74ce9"
H_LEDGER_SHA256 = "b146c11f33e895c08f72303d282e2668d955e0a58d9268a1b445d4d5202016c2"
SLICE_SCRIPT_LF_SHA256 = "f2b8725341caea3bc2235dea9b69e7d33c5870a89db96a8604ef4521b2154659"
SLICE_OUTPUT_LF_SHA256 = "1086ccb03a69f1c92eb945a6ad0118976ae7f6ef2763d4e10e3973256a53d503"
SLICE_12_K_SHA256 = "de18f1f38b29b92cbbe46c913a1446fe017a4566615d60e37a96d82307ae84a7"
BOUNDARY_P_LEDGER_SHA256 = "5f15ad89e0788eafcdffecf6fa0f0204224d27fb6ecf94cd46ed1982da678333"
EXPECTED_SEMANTIC_SHA256 = "5c12ecc88284a89667bcbbbbcaa505d221c04a4ccf70a1f9dfdc0456439f570c"


def lf_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def coefficient_hash(poly: sp.Poly) -> str:
    ledger = "\n".join(
        f"{monomial}:{coefficient}" for monomial, coefficient in poly.terms()
    )
    return hashlib.sha256(ledger.encode("ascii")).hexdigest()


H_bytes = H_PATH.read_bytes()
require(hashlib.sha256(H_bytes).hexdigest() == H_ARTIFACT_SHA256, "H artifact changed")
H = pickle.loads(H_bytes)

a, b, c, w = sp.symbols("a b c w")
H_poly = sp.Poly(H, a, b, c)
require(coefficient_hash(H_poly) == H_LEDGER_SHA256, "H coefficient ledger changed")
require(H_poly.total_degree() == 25 and len(H_poly.terms()) == 361, "H shape changed")

L = 27 * a**2 * c**2 - 18 * a * b * c + 16 * a + b**3 * c - b**2
T = 4 - 3 * b * c
A = 12 * a - b**2
E = L * w**3 + T * w - 2 * c
require(sp.factor(L) == L, "old Jelonek divisor L is no longer irreducible")


def F(point: tuple[sp.Expr, sp.Expr, sp.Expr]) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    x, y, z = point
    unit = 1 + x * y
    return (
        unit**3 * z + y**2 * unit * (4 + 3 * x * y),
        y + 3 * x * unit**2 * z + 3 * x * y**2 * (4 + 3 * x * y),
        2 * x - 3 * x**2 * y - x**3 * z,
    )


# The global inverse graph, typed in the cubic algebra Q(a,b,c)[w]/(E).
y_denominator = A * w**2 + b * w + 2
qy = sp.cancel(b - 3 * a * w * ((9 * a * c - b) * w + 2) / y_denominator)
qz = sp.cancel(((2 * w - c) - 3 * w**2 * qy) / w**3)
coefficient_field = sp.QQ.frac_field(a, b, c)
E_monic = sp.Poly(E, w, domain=coefficient_field).monic()
for row, target in zip(F((w, qy, qz)), (a, b, c)):
    numerator = sp.cancel(row - target).as_numer_denom()[0]
    remainder = sp.Poly(numerator, w, domain=coefficient_field).rem(E_monic)
    require(remainder.is_zero, "rational inverse row is not typed modulo E")


# On L=0 with A*T nonzero, two roots escape.  Their inverse coordinates obey
# qy -> eta and qz ~ -3*eta/w.  The largest exponent i-k in a^i b^j c^k
# is seven, so this exact four-term face controls each escaping H-value.
eta = sp.symbols("eta")
escape_order = max(i - k for (i, _j, k), _coefficient in H_poly.terms())
escape_face = sp.expand(
    sum(
        coefficient * eta**j * (-3 * eta) ** k
        for (i, j, k), coefficient in H_poly.terms()
        if i - k == escape_order
    )
)
escape_constant = sp.Integer(2) ** 9 * 3**6 * 11**3 * 13**2
require(escape_order == 7, "escaping H-order changed")
require(escape_face == escape_constant * eta**3, "escaping H-face changed")

eta_infinity = sp.factor(b - 3 * a * (9 * a * c - b) / A)
require(
    sp.cancel(eta_infinity - (15 * a * b - b**3 - 27 * a**2 * c) / A) == 0,
    "escaping qy limit changed",
)


# Pull the generic boundary calculation to THM-2576's set-bijective
# normalization nu(tau,lambda) of V(L).
tau, lam = sp.symbols("tau lambda")
nu_a = lam**2 * (3 - tau * lam) / 27
nu_b = lam * (4 - tau * lam) / 3
nu_c = tau
nu = {a: nu_a, b: nu_b, c: nu_c}
seam = lam * tau - 2

require(sp.factor(L.subs(nu)) == 0, "nu left V(L)")
require(sp.factor(A.subs(nu)) == -lam**2 * seam**2 / 9, "A pullback changed")
require(sp.factor(T.subs(nu)) == seam**2, "T pullback changed")
require(sp.factor(eta_infinity.subs(nu)) == lam / 3, "escaping qy pullback changed")

# The unique finite root modulo L and its surviving inverse point.
w0 = 2 * nu_c / seam**2
qy0 = sp.factor(qy.subs(nu).subs(w, w0))
qz0 = sp.factor(qz.subs(nu).subs(w, w0))
expected_qy0 = -lam * (3 * lam * tau - 8) / 6
expected_qz0 = -lam**2 * seam**2 * (lam**2 * tau**2 - 8 * lam * tau + 14) / 8
require(qy0 == expected_qy0, "finite-sheet qy changed")
require(qz0 == expected_qz0, "finite-sheet qz changed")

# This is the only moderately large operation.  SymPy performs the exact
# substitution; FLINT then factors the primitive bivariate numerator.
H_finite = sp.cancel(H.subs({a: w0, b: qy0, c: qz0}))
finite_numerator, finite_denominator = H_finite.as_numer_denom()
expected_finite_denominator = 2**35 * 3**21 * seam**14
require(
    sp.expand(finite_denominator - expected_finite_denominator) == 0,
    "finite-sheet H denominator changed",
)

finite_numerator_poly = sp.Poly(finite_numerator, tau, lam)
P_poly = finite_numerator_poly.exquo(sp.Poly(lam**2, tau, lam))
content, P_primitive = P_poly.primitive()
require(content == 1 and P_primitive == P_poly, "boundary P lost primitivity")
require(
    (P_poly.degree(tau), P_poly.degree(lam), P_poly.total_degree(), len(P_poly.terms()))
    == (57, 84, 141, 527),
    "boundary P shape changed",
)
require(coefficient_hash(P_poly) == BOUNDARY_P_LEDGER_SHA256, "boundary P ledger changed")

flint_ctx.threads = 4
P_context = fmpz_mpoly_ctx.get(["tau", "lambda"], ordering="lex")
P_flint = P_context.from_dict(
    {monomial: int(coefficient) for monomial, coefficient in P_poly.terms()}
)
factor_content, factors = P_flint.factor()
require(
    factor_content == 1
    and len(factors) == 1
    and factors[0][1] == 1
    and factors[0][0] == P_flint,
    "FLINT found a nontrivial factor of boundary P",
)

hostile_finite_value = sp.Rational(
    3393794313700733412883215882425216567,
    359414999291950792704,
)
require(H_finite.subs({tau: 1, lam: 1}) == hostile_finite_value, "finite-sheet hostile changed")
require(P_poly.eval({tau: 1, lam: 2}) != 0, "exceptional seam unexpectedly killed P")

# Since L*w1*w2 -> T for the two escaping roots, the exact residue is
# C^2*T^7*eta^6*H(q0).  The seam powers cancel on the normalization.
residue_constant = sp.cancel(escape_constant**2 / (2**35 * 3**27))
expected_residue_constant = sp.Rational(11**6 * 13**4, 2**17 * 3**15)
require(residue_constant == expected_residue_constant, "boundary residue constant changed")


# A small exact rank/distinctness control for the newest image surface.
source_boundary_point = (sp.Integer(0), sp.Integer(-1), sp.Integer(-1))
H_point = tuple(sp.factor(value) for value in F(source_boundary_point))
J_point = tuple(sp.factor(value) for value in F(H_point))
require(sp.expand(L.subs(dict(zip((a, b, c), source_boundary_point)))) == 0, "source point left L")
require(H_point == (3, -1, 0), "small H point changed")
require(H.subs(dict(zip((a, b, c), H_point))) == 0, "small point left H")
gradient_H = tuple(
    sp.diff(H, variable).subs(dict(zip((a, b, c), H_point))) for variable in (a, b, c)
)
require(
    gradient_H == (-1152495599044, -10372460391396, 43218584964150),
    "H smooth-point gradient changed",
)
require(J_point == (10, -46, 33), "small newest-image point changed")
require(L.subs(dict(zip((a, b, c), J_point))) == -504, "new point met old L")
require(
    H.subs(dict(zip((a, b, c), J_point))) == -1402696598666966597632,
    "new point met old H",
)
jacobian = sp.Matrix(F((a, b, c))).jacobian((a, b, c))
require(sp.factor(jacobian.det()) == -2, "Keller determinant changed")


# Pin (rather than silently restate) the inherited exact slice.  Its
# irreducible exponent-one specialization rules out J0^e with e>1 after the
# global geometry has already shown that only one prime image support occurs.
require(lf_sha256(SLICE_SCRIPT) == SLICE_SCRIPT_LF_SHA256, "slice script changed")
require(lf_sha256(SLICE_OUTPUT) == SLICE_OUTPUT_LF_SHA256, "slice output changed")
slice_text = SLICE_OUTPUT.read_text(encoding="utf-8")
for fragment in (
    "slice (b,c)=(1,2)",
    "Norm(H)=K/(2^21 L^7); deg(K)=86; K irreducible; gcd(K,LH)=1",
    f"K coefficient-ledger sha256={SLICE_12_K_SHA256}",
    "all exact checks passed",
):
    require(fragment in slice_text, "inherited slice certificate lost a pinned statement")


semantic_payload = "\n".join(
    (
        f"H={H_ARTIFACT_SHA256}",
        f"H_ledger={H_LEDGER_SHA256}",
        f"escape_order={escape_order}",
        f"escape_constant={escape_constant}",
        f"P_shape=57,84,141,527",
        f"P_ledger={BOUNDARY_P_LEDGER_SHA256}",
        f"residue_constant={residue_constant}",
        f"image_point={J_point}",
        "image_point_L=-504",
        "image_point_H=-1402696598666966597632",
        f"slice_script={SLICE_SCRIPT_LF_SHA256}",
        f"slice_output={SLICE_OUTPUT_LF_SHA256}",
        f"slice_12_K={SLICE_12_K_SHA256}",
        "global_v_L_NH=-7",
        "global_shape=Norm(H)=u*J/L^7 with J prime and image multiplicity one",
    )
)
semantic_sha256 = hashlib.sha256(semantic_payload.encode("ascii")).hexdigest()
if EXPECTED_SEMANTIC_SHA256 != "UNFROZEN":
    require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256, "semantic payload changed")


print("== fixed-F level-three global norm: exact L-pole and residue ==")
print("global inverse graph: 3/3 rows vanish modulo the typed cubic E")
print(f"escaping H-face: order {escape_order}, coefficient ({escape_constant})*eta^3")
print("two escaping roots + one finite root force v_L(Norm(H))=-7")
print("THM-2473 finiteness off L then gives L^7*Norm(H) in Q[a,b,c], not divisible by L")
print(
    "boundary finite-sheet identity: H(q0)=lambda^2*P/"
    "(2^35*3^21*(lambda*tau-2)^14)"
)
print("P: bidegree (57,84), total degree 141, 527 terms, FLINT-irreducible")
print(f"P coefficient-ledger sha256={BOUNDARY_P_LEDGER_SHA256}")
print(
    "normalization residue: nu^*(L^7*Norm(H)|_L)="
    "(11^6*13^4/(2^17*3^15))*lambda^8*P"
)
print("exceptional T=A=0 seam is excluded from the root split but survives after cancellation")
print("rank/distinctness hostile: H-point (3,-1,0) -> J-point (10,-46,33)")
print("  at J-point: L=-504 and H=-1402696598666966597632")
print("pinned (b,c)=(1,2) slice: irreducible degree-86 exponent-one numerator")
print("therefore the unique global image prime occurs with multiplicity one and differs from L,H")
print(f"semantic sha256={semantic_sha256}")
print("scope: no expanded J, no degree-27 separability, no family/classification/JC(2) claim")
print("all exact checks passed")
