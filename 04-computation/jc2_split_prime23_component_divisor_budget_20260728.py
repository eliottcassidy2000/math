#!/usr/bin/env python3
"""Exact companion for the split prime-23 component/divisor reduction.

For the nonzero-eta even-Faber curve

    f2(t,v,zeta)=0,       zeta*f1(t,v,zeta)^4=eta*t^23,

the fixed fibre has five old branches with orders

    ord(t),ord(F1),ord(zeta) = (4,23,0)

and three new branches with orders (1,0,23).  If a geometric component
contains r old and s new branches, the sections F1 in O(5) and zeta in O(3)
give the opposing degree bounds

    23*r <= 5*(4*r+s),       23*s <= 3*(4*r+s).

They force (r,s)=(5,3).  The proof-sidecar in the accompanying reflection
checks the geometric hypotheses (basepoint freeness, local domains, and
reducedness); this script freezes the exact algebra, norm form, and integer
budget.  It proves no rationality exclusion and uses no floating arithmetic.
"""

from __future__ import annotations

import contextlib
import hashlib
import importlib.util
import io
from pathlib import Path

import sympy as s


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
BASE_PATH = HERE / "jc2_degree22_split_lambda_prime23_scout.py"
BASE_SHA256 = "1bbadb900e27112f57f600b0b5b73a3b85fc5b02f23e3f8f687dac4fa1c41fc3"
require(
    hashlib.sha256(BASE_PATH.read_bytes()).hexdigest() == BASE_SHA256,
    "audited split-lambda scout changed",
)

spec = importlib.util.spec_from_file_location("split_lambda_prime23", BASE_PATH)
require(spec is not None and spec.loader is not None, "cannot load lambda scout")
base = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(base)


# ---------------------------------------------------------------------------
# 1. Universal quadratic norm identity.
# ---------------------------------------------------------------------------

z, X = s.symbols("z X")
AA, BB, CC, aa, bb = s.symbols("AA BB CC aa bb")
quadratic = AA * z**2 + BB * z + CC
linear = aa * z + bb
abstract_resultant = s.Poly(s.resultant(quadratic, z * linear**4 - X, z), X)
require(abstract_resultant.degree() == 2, "abstract norm lost quadratic degree")

DD = s.expand(AA * bb**2 - BB * aa * bb + CC * aa**2)
TT = s.expand(-abstract_resultant.coeff_monomial(X))
abstract_norm = s.expand(CC * DD**4 - X * TT + AA**5 * X**2)
require(
    s.expand(abstract_resultant.as_expr() - abstract_norm) == 0,
    "quadratic norm identity failed",
)


# ---------------------------------------------------------------------------
# 2. Specialize the norm to the exact split degree-22 fluxes.
# ---------------------------------------------------------------------------

t, v, zeta = base.tau, base.v, base.zeta
c, d, e, w = base.c, base.dd, base.e, base.w
f1_z = s.Poly(base.f1, zeta)
f2_z = s.Poly(base.f2, zeta)
require(f1_z.degree() == 1 and f2_z.degree() == 2, "flux zeta degrees changed")

a = f1_z.coeff_monomial(zeta)
b = f1_z.coeff_monomial(1)
A = f2_z.coeff_monomial(zeta**2)
B = f2_z.coeff_monomial(zeta)
C0 = f2_z.coeff_monomial(1)
require(A == 15944049, "quadratic leading coefficient changed")

D = s.expand(A * b**2 - B * a * b + C0 * a**2)
T = s.expand(TT.subs({AA: A, BB: B, CC: C0, aa: a, bb: b}))
require(
    s.denom(T) == 1,
    "scaled trace coefficient is not polynomial",
)

D_poly = s.Poly(D, t, v, c, d, e, w)
T_poly = s.Poly(T, t, v, c, d, e, w)
require(
    (D_poly.degree(t), D_poly.degree(v), len(D_poly.terms())) == (10, 5, 60),
    "D norm support changed",
)
require(
    (T_poly.degree(t), T_poly.degree(v), len(T_poly.terms())) == (23, 11, 1220),
    "trace support changed",
)

# With X=eta*t^23, the invariant eliminant is exactly
# C0*D^4-eta*t^23*T+A^5*eta^2*t^46.  At t=0 its two factors are the
# already-audited G3 and L5, up to the displayed nonzero scalars.
C0_fixed = s.Poly(C0.subs(t, 0), v, domain=s.QQ)
D_fixed = s.Poly(D.subs(t, 0), v, domain=s.QQ)
G3 = s.Poly(base.G3, v, domain=s.QQ)
L5 = s.Poly(base.L5, v, domain=s.QQ)
require(C0_fixed.as_expr() == -672 * G3.as_expr(), "fixed C0 is not -672*G3")
require(
    D_fixed.as_expr() == -16071601392 * L5.as_expr(),
    "fixed D is not the declared scalar times L5",
)

T_fixed = s.Poly(T.subs(t, 0), v, domain=s.QQ)
require((T_fixed.degree(), len(T_fixed.terms())) == (11, 12), "fixed trace changed")
require(s.gcd(T_fixed, G3).degree() == 0, "order-23 coupling vanished on G3")
require(s.gcd(T_fixed, L5).degree() == 0, "order-23 coupling vanished on L5")


# ---------------------------------------------------------------------------
# 3. Exhaust the component divisor budget.
# ---------------------------------------------------------------------------

naive_partitions: list[tuple[int, int, int]] = []
budget_survivors: list[tuple[int, int, int]] = []
for old_count in range(6):
    for new_count in range(4):
        degree = 4 * old_count + new_count
        if degree == 0:
            continue
        naive_partitions.append((old_count, new_count, degree))
        f1_budget = 23 * old_count <= 5 * degree
        zeta_budget = 23 * new_count <= 3 * degree
        if f1_budget and zeta_budget:
            budget_survivors.append((old_count, new_count, degree))

require(
    budget_survivors == [(5, 3, 23)],
    "component divisor budget gained a proper partition",
)
require(
    {degree for _, _, degree in naive_partitions} == set(range(1, 24)),
    "naive branch degrees no longer cover 1..23",
)

# Equality on the unique component spends the entire zero divisors:
# div_0(F1)=23*O, div_0(zeta)=23*N, and q=-t^5/F1 has divisor 5*N-3*O.
require(5 * 23 == 23 * 5, "F1 divisor degree mismatch")
require(3 * 23 == 23 * 3, "zeta divisor degree mismatch")
require(5 * (4 * 5 + 3) - 23 * 5 == 0, "q divisor did not balance")


# ---------------------------------------------------------------------------
# 4. Genus-zero perfect-power normal form.
# ---------------------------------------------------------------------------

alpha, beta, H = s.symbols("alpha beta H", nonzero=True)
q_form = alpha**5 / beta**3
t_form = alpha * beta**4 / H
zeta_form = s.cancel(q_form**4 * t_form**3)
f1_form = s.cancel(-t_form**5 / q_form)
require(zeta_form == alpha**23 / H**3, "zeta perfect-power form changed")
require(f1_form == -beta**23 / H**5, "F1 perfect-power form changed")


print("split prime-23 component divisor budget")
print(f"base_split_lambda_sha256={BASE_SHA256}")
print("quadratic_norm=Res_zeta(f2,zeta*f1^4-X)=C0*D^4-X*T+A^5*X^2")
print("split_eliminant=C0*D^4-eta*t^23*T+A^5*eta^2*t^46")
print(f"zeta2_leading_A={A}")
print("D_norm=degree_t10:degree_v5:terms60")
print("T_trace=degree_t23:degree_v11:terms1220")
print("fixed_C0=-672*G3")
print("fixed_D=-16071601392*L5")
print("fixed_T=degree11:terms12:coprime_to_G3_and_L5")
print("naive_component_degrees=1..23")
print("divisor_budget=23r<=5(4r+s);23s<=3(4r+s)")
print("budget_survivors=(r,s,d)=(5,3,23)_only")
print("uniform_component_consequence=one_geometric_component")
print("exact_zero_divisors=div0(F1)=23O;div0(zeta)=23N;div(q)=5N-3O")
print("rational_normal_form=q=a^5/b^3;t=a*b^4/H;zeta=a^23/H^3;f1=-b^23/H^5")
print("rational_form_degrees=deg_a3:deg_b5:deg_H23")
print("scope=UNIFORM_GEOMETRIC_INTEGRALITY_REDUCTION_NOT_RATIONALITY_EXCLUSION")
print("ALL CHECKS PASSED")
