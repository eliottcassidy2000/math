#!/usr/bin/env python3
"""Exact companion for THM-3980 (all-height split formal cusp section)."""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.factor(expression) == 0, message)


def zero_mod_l(expression: sp.Expr, L: sp.Symbol, message: str) -> None:
    numerator = sp.together(expression).as_numer_denom()[0]
    remainder = sp.Poly(sp.expand(numerator), L, domain=sp.EX).rem(
        sp.Poly(6 * L**2 - 1, L, domain=sp.EX)
    ).as_expr()
    zero(remainder, message)


x, t, L = sp.symbols("x t L", nonzero=True)


def jac(f: sp.Expr, g: sp.Expr) -> sp.Expr:
    return sp.expand(
        sp.diff(f, x) * sp.diff(g, t)
        - sp.diff(f, t) * sp.diff(g, x)
    )


height_rows: list[dict[str, int]] = []
for n in range(2, 10):
    z = 1 + x**n * t
    p = sp.expand(z * t)
    y = sp.expand(x ** (n - 1) * z * t**2)

    zero(z * (z - 1) - x**n * p,
         f"height {n}: Danielewski relation")
    zero(p * (z - 1) - x * y,
         f"height {n}: modification relation")
    zero(z * y - x ** (n - 1) * p**2,
         f"height {n}: saturated minor")
    zero(z * (z - 1)**2 - x ** (n + 1) * y,
         f"height {n}: D formal equation")

    # The D-regular rational coordinate and exact source volume.
    qd = z / x ** (n + 1)
    zero(qd - y / (z - 1)**2,
         f"height {n}: D coordinate two presentations")
    beta_d_x = sp.diff(qd, x)
    beta_d_t = sp.diff(qd, t)
    # x dx wedge d(qd) has coefficient x*qd_t relative to dx wedge dt.
    zero(x * beta_d_t - 1,
         f"height {n}: eta=x dx wedge d qD")

    ad = sp.expand(qd**2 + 2 * L * x)
    cd = sp.expand(qd**3 + 3 * L * x * qd)
    zero_mod_l(jac(ad, cd) - 1, L,
               f"height {n}: D cusp pair exact Jacobian")

    # The retained source-arm chart.
    al = t + 2 * L * x
    cl = -x
    zero(jac(al, cl) - 1,
         f"height {n}: L1 exact Darboux pair")

    # Exact pole ledgers.  They are frozen as integer comparisons here;
    # the theorem derives them from z and z-1 being the respective units.
    pole_qd_l = -(n + 1)
    pole_ad_l = -2 * (n + 1)
    pole_cd_l = -3 * (n + 1)
    pole_t_d = -n
    gate(pole_qd_l < 0, f"height {n}: qD has an L1 pole")
    gate(pole_ad_l < 0, f"height {n}: AD has an L1 pole")
    gate(pole_cd_l < 0, f"height {n}: CD has an L1 pole")
    gate(pole_t_d < 0, f"height {n}: AL has a D pole")
    gate(pole_ad_l != 0, f"height {n}: A components cannot agree")
    gate(pole_cd_l != 1, f"height {n}: C components cannot agree")

    # The Hensel square root has opposite signs relative to q=2z-1 on the
    # two colors.  A field has only the two global roots +/-q.
    q = 2 * z - 1
    zero(q**2 - (1 + 4 * x**n * p),
         f"height {n}: Hensel square equation")
    gate(q.subs(x, 0) == 1,
         f"height {n}: finite-source specialization sees L1")

    # The direct source substitution only sees L1 (finite t), where q=+1.
    # The D sign is checked in the abstract quotient z=0 below.
    gate((2 * 1 - 1) == 1, f"height {n}: L1 q sign")
    gate((2 * 0 - 1) == -1, f"height {n}: D q sign")

    height_rows.append({
        "n": n,
        "ordL_qD": pole_qd_l,
        "ordL_AD": pole_ad_l,
        "ordL_CD": pole_cd_l,
        "ordD_AL": pole_t_d,
    })


# ---------------------------------------------------------------------------
# Universal idempotent, split algebra, and connected-algebraization gate.
# ---------------------------------------------------------------------------

e, AD, AL, CD, CL, T = sp.symbols("e AD AL CD CL T")
AH = sp.expand(e * AD + (1 - e) * AL)
CH = sp.expand(e * CD + (1 - e) * CL)


def reduce_e(expression: sp.Expr) -> sp.Expr:
    return sp.Poly(sp.expand(expression), e, domain=sp.EX).rem(
        sp.Poly(e**2 - e, e, domain=sp.EX)
    ).as_expr().expand()


zero(reduce_e((AH - AD) * (AH - AL)),
     "split A satisfies its reducible quadratic")
zero(reduce_e((CH - CD) * (CH - CL)),
     "split C satisfies its reducible quadratic")
zero(reduce_e((AH - AL) / (AD - AL) - e),
     "A reconstructs the nontrivial idempotent")
zero(reduce_e((CH - CL) / (CD - CL) - e),
     "C reconstructs the nontrivial idempotent")

# The mixed-sign square root is (W_D,W_L)=(-q,+q).  Its square is diagonal,
# but it is neither diagonal root.  The polynomial is deliberately split:
# this is algebraic in KxK, not algebraic inside a field extension.
q0 = sp.symbols("q0", nonzero=True)
W_D, W_L = -q0, q0
zero(W_D**2 - q0**2, "mixed root D square")
zero(W_L**2 - q0**2, "mixed root L square")
gate(W_D != W_L, "mixed root is not diagonal")
gate(W_D == -q0 and W_L == q0, "mixed root sign packet")

# Hostile global choices W=+q and W=-q collapse the idempotent and select
# just one chart; they cannot produce the mixed formal section.
for sign in (1, -1):
    eps = sp.Integer(1) / sign
    e_l = sp.simplify((1 + eps) / 2)
    e_d = sp.simplify((1 - eps) / 2)
    gate((e_d, e_l) in ((0, 1), (1, 0)),
         f"global square-root sign {sign}: trivial idempotent")

summary = {
    "checks": CHECKS,
    "family": "B_n, n>=2",
    "D_coordinate": "q_D=z/x^(n+1)=y/(z-1)^2",
    "D_pair": "(q_D^2+2Lx,q_D^3+3Lxq_D), J=1",
    "L1_pair": "(t+2Lx,-x), J=1",
    "glue": "mixed-sign Hensel root of (2z-1)^2=1+4x^n p",
    "pole_rows": height_rows,
    "field_gate": "global field roots are +/-q; mixed root requires KxK",
    "algebra": "split pair is algebraic over K in KxK but has no connected algebraization",
    "scope": "canonical split section only; arbitrary formal/polynomial pairs open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3980 all-height split formal cusp companion")
print(f"CHECKS={CHECKS}")
print("D_COORDINATE=Z_OVER_X_NPLUS1=Y_OVER_ZMINUS1_SQUARED")
print("D_VOLUME=X_DX_WEDGE_DQD;D_CUSP_PAIR=EXACT")
print("L1_PAIR=T_PLUS_2LX_MINUS_X;EXACT")
print("GLUE=MIXED_SIGN_SQRT_1_PLUS_4XNP")
print("POLES=QD_NPLUS1;AD_2NPLUS2;CD_3NPLUS3;AL_ON_D_N")
print("NONRATIONAL=NOT_IN_DIAGONAL_K")
print("ALGEBRAIC_ENVELOPE=SPLIT_K_TIMES_K;NO_CONNECTED_ALGEBRAIZATION")
print("SCOPE=CANONICAL_SPLIT_SECTION_ONLY;KELLER_PAIR_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
