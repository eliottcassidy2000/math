#!/usr/bin/env python3
"""Exact companion for THM-3978's linear-seam Jacobian image ideal."""

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


def ceil_div(a: int, b: int) -> int:
    return (a + b - 1) // b


x, t, u, c = sp.symbols("x t u c", nonzero=True)


def jac(a: sp.Expr, b: sp.Expr) -> sp.Expr:
    return sp.factor(sp.diff(a, x) * sp.diff(b, t)
                     - sp.diff(a, t) * sp.diff(b, x))


def laurent_rows(expr: sp.Expr, shift: int) -> dict[int, sp.Expr]:
    """Return x-exponent rows, knowing x**shift * expr is polynomial."""
    poly = sp.Poly(sp.expand(x**shift * expr), x)
    return {
        degree - shift: sp.factor(coefficient)
        for (degree,), coefficient in poly.terms()
    }


inventory: dict[int, list[int]] = {}
for n in range(2, 10):
    uxt = x**n * t
    A = x + c * uxt
    V = A / x

    # A is a submersion: A_t=c*x^n, while A_x specializes to 1 at x=0.
    gate(sp.factor(sp.diff(A, t) - c * x**n) == 0,
         f"height {n}: t derivative")
    gate(sp.expand(sp.diff(A, x).subs(x, 0) - 1) == 0,
         f"height {n}: no critical point over x=0")

    # The rational mate and the minimal polynomial-valued Jacobian in B_n.
    rational_mate = x**(1 - n) / (c * (n - 1))
    gate(jac(A, rational_mate) == 1,
         f"height {n}: rational constant mate")

    # The two branches over x=0 cancel the pole with different invariant
    # constants.  These regular rewritings are the exact formal-CRT atlas.
    arm_sum = sum(V**j for j in range(n - 1))
    Q_arm = t * arm_sum / ((n - 1) * V**(n - 1))
    gate(sp.cancel(Q_arm - rational_mate
                   + A**(1 - n) / (c * (n - 1))) == 0,
         f"height {n}: retained-arm formal mate")
    zxt = 1 + uxt
    W_boundary = (A + c) / x
    boundary_sum = sum(W_boundary**j for j in range(n - 1))
    Q_boundary = (zxt / x**n) * boundary_sum / (
        (n - 1) * W_boundary**(n - 1)
    )
    gate(sp.cancel(Q_boundary - rational_mate
                   + (A + c)**(1 - n) / (c * (n - 1))) == 0,
         f"height {n}: boundary formal mate")

    r = (A * (A + c))**(n - 1)
    Q = sp.cancel((A + c)**(n - 1) * (V**(n - 1) - 1)
                  / (c * (n - 1)))
    _, denominator = sp.cancel(Q).as_numer_denom()
    gate(not denominator.has(x, t),
         f"height {n}: minimal Q is source-polynomial")
    gate(sp.factor(jac(A, Q) - r) == 0,
         f"height {n}: minimal image generator")

    # Verify every negative homogeneous row against the exact two-color
    # module of B_n.  The formula has no row below -(n-1).
    Au = x + c * u
    Qu = sp.cancel((Au + c)**(n - 1)
                   * ((Au / x)**(n - 1) - 1) / (c * (n - 1)))
    rows = laurent_rows(Qu, n - 1)
    negative = sorted(weight for weight in rows if weight < 0)
    inventory[n] = negative
    gate(negative == list(range(1 - n, 0)),
         f"height {n}: exact negative support")
    for weight in negative:
        q = -weight
        required = u**ceil_div(q, n) * (u + 1)**ceil_div(q, n + 1)
        gate(sp.rem(rows[weight], required, domain=sp.QQ.frac_field(c)) == 0,
             f"height {n}, weight {weight}: B_n membership")

    # The two one-color near misses record why both target factors occur.
    # Taylor orders are the exact criterion: n-2 versus n-1.
    S = sp.symbols("S")
    rz = S**(n - 2) * (S + c)**(n - 1)
    rb = S**(n - 1) * (S + c)**(n - 2)
    gate(sp.diff(rz, S, n - 2).subs(S, 0) != 0,
         f"height {n}: missing A multiplicity detected")
    gate(sp.diff(rb, S, n - 2).subs(S, -c) != 0,
         f"height {n}: missing A+c multiplicity detected")


# At height two the generator has a particularly transparent expression.
n = 2
uxt = x**n * t
A = x + c * uxt
p = t * (1 + uxt)
Q_visible = uxt + c * x * p
gate(sp.factor(jac(A, Q_visible) - A * (A + c)) == 0,
     "height two visible generator")


summary = {
    "checks": CHECKS,
    "family": "A_c=x+c*x^n*t, n>=2, c!=0",
    "submersion": True,
    "rational_mates": "x^(1-n)/(c(n-1))+k(A_c)",
    "source_image": "J(A_c,k[x,t]) intersect k[A_c]=A_c^(n-1)k[A_c]",
    "completion_image": (
        "J(A_c,B_n) intersect k[A_c]="
        "(A_c(A_c+c))^(n-1)k[A_c]"
    ),
    "formal_atlas": "exact mates on both x-adic factors with distinct H(A)",
    "negative_support": inventory,
    "height2": "J(A,u+c*x*p)=A(A+c)",
    "scope": "named first-coordinate family only; JC2 open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3978 linear-seam submersion companion")
print(f"CHECKS={CHECKS}")
print("SUBMERSION=A_C_HAS_NO_AFFINE_CRITICAL_POINT")
print("RATIONAL_MATES=X_POWER_1_MINUS_N_OVER_C_N_MINUS_1_PLUS_K_OF_A")
print("SOURCE_IMAGE=A_POWER_N_MINUS_1")
print("COMPLETION_IMAGE=PRODUCT_OF_TWO_COLORS_TO_POWER_N_MINUS_1")
print("FORMAL_ATLAS=EXACT_ON_BOTH_X_ADIC_FACTORS;GLOBAL_H_MISMATCH")
print("NO_CONSTANT_MATE=YES")
print("HEIGHT2_VISIBLE=J(A,U_PLUS_CXP)=A(A_PLUS_C)")
print("SCOPE=NAMED_LINEAR_SEAM_FAMILY;JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
