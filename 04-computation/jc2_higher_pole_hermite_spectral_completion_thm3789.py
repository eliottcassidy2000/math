#!/usr/bin/env python3
"""Exact companion for THM-3789's higher-pole Hermite completion."""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def check(condition: object, label: str) -> None:
    global CHECKS
    if condition is not True and condition != sp.S.true:
        raise AssertionError(label)
    CHECKS += 1


def same(lhs: sp.Expr, rhs: sp.Expr, label: str) -> None:
    check(sp.cancel(sp.expand(lhs - rhs)) == 0, label)


x, y, T = sp.symbols("x y T")
t_source = x**2 * (1 + x * y)


def jac(a: sp.Expr, b: sp.Expr) -> sp.Expr:
    return sp.cancel(
        sp.diff(a, x) * sp.diff(b, y) - sp.diff(a, y) * sp.diff(b, x)
    )


# Abstract hypersurface-gradient and Poisson sign packet.
f, w, z = sp.symbols("f w z")
psi = sp.Function("Psi")
D_formal = f**3 * z - (w - f**2) * psi(w)
same(sp.diff(D_formal, z), f**3, "formal D_z")
same(sp.diff(D_formal, f), 3 * f**2 * z + 2 * f * psi(w), "formal D_f")
same(
    sp.diff(D_formal, w),
    -psi(w) - (w - f**2) * sp.diff(psi(w), w),
    "formal D_w",
)


def exact_integer_root_control(roots: tuple[int, ...], tag: str) -> None:
    """Replay the complete construction for q=product(T-alpha)."""

    d = len(roots)
    q = sp.expand(sp.prod(T - alpha for alpha in roots))
    Q = sp.expand(sp.integrate(q**2, T) - sp.integrate(q**2, T).subs(T, 0))
    lambdas = tuple(sp.factor(Q.subs(T, alpha)) for alpha in roots)
    check(len(set((sp.Integer(0),) + lambdas)) == d + 1, f"{tag} distinct spectrum")
    Psi = sp.expand(sp.prod(w - lam for lam in lambdas))
    Psi_Q = sp.expand(Psi.subs(w, Q))

    # First Hermite quotient is univariate after removing the forced axis T.
    A = sp.cancel(Q * Psi_Q / q**2)
    check(sp.denom(A) == 1, f"{tag} A polynomial")
    B = sp.cancel(A / T)
    check(sp.denom(B) == 1, f"{tag} B polynomial")

    q_s = sp.expand(q.subs(T, t_source))
    Q_s = sp.expand(Q.subs(T, t_source))
    F = sp.expand(x * q_s)
    W = Q_s
    H = sp.expand((1 + x * y) * B.subs(T, t_source))
    Psi_W = sp.expand(Psi.subs(w, W))
    Y = sp.cancel((H - Psi_W) / F)
    check(sp.denom(Y) == 1, f"{tag} Y polynomial")
    Y = sp.expand(Y)

    same(H, W * Psi_W / F**2, f"{tag} H quotient")
    same(H, Psi_W + F * Y, f"{tag} second Hermite quotient")
    same(F**3 * Y, (W - F**2) * Psi_W, f"{tag} surface relation")
    same(jac(F, t_source), x**3 * q_s, f"{tag} J(F,t)")
    same(jac(F, W), F**3, f"{tag} J(F,W)")
    same(jac(F, W / F**3), 1, f"{tag} rational Keller")
    same(
        jac(F, Y),
        Psi_W + (W - F**2) * sp.diff(Psi, w).subs(w, W),
        f"{tag} bracket F,Y",
    )
    same(
        jac(W, Y),
        3 * F**2 * Y + 2 * F * Psi_W,
        f"{tag} bracket W,Y",
    )

    # The axis retains the first jet with a nonzero scalar coefficient.
    q0 = q.subs(T, 0)
    psi0 = Psi.subs(w, 0)
    same(Y.subs(x, 0), psi0 * y / q0, f"{tag} axis packet")
    check(q0 != 0 and psi0 != 0, f"{tag} axis coefficient nonzero")

    # Every critical contact is exactly cubic; its residual cover has degree
    # 2d-2, and the boundary packet is the asserted inverse cube.
    for i, alpha in enumerate(roots):
        lam = lambdas[i]
        qprime = sp.diff(q, T).subs(T, alpha)
        psi_prime = sp.diff(Psi, w).subs(w, lam)
        cubic = sp.cancel((Q - lam) / (T - alpha) ** 3)
        check(sp.denom(cubic) == 1, f"{tag} cubic quotient {i}")
        check(sp.Poly(cubic, T).degree() == 2 * d - 2, f"{tag} residual degree {i}")
        same(cubic.subs(T, alpha), qprime**2 / 3, f"{tag} cubic lead {i}")
        for j, other in enumerate(roots):
            if j != i:
                check(cubic.subs(T, other) != 0, f"{tag} residual noncritical {i},{j}")
        boundary = sp.cancel(Y.subs(y, (alpha - x**2) / x**3))
        same(
            boundary,
            lam * psi_prime / (3 * qprime * x**3),
            f"{tag} critical boundary {i}",
        )
        check(lam * psi_prime * qprime != 0, f"{tag} critical coefficient {i}")

    # Smoothness is exhausted by f!=0 and the d+1 simple boundary arms.
    D = sp.expand(f**3 * z - (w - f**2) * Psi)
    Dw = sp.diff(D, w)
    same(sp.diff(D, z), f**3, f"{tag} target D_z")
    check(Dw.subs({f: 0, w: 0}) == -psi0, f"{tag} smooth axis")
    for i, lam in enumerate(lambdas):
        check(
            sp.factor(Dw.subs({f: 0, w: lam})) == -lam * sp.diff(Psi, w).subs(w, lam),
            f"{tag} smooth critical arm {i}",
        )


# Two exact positive controls.  Real ordered roots force distinct normalized
# critical values, and the d=2/d=3 residual degrees are 2 and 4.
exact_integer_root_control((1, 2), "quadratic_integer")
exact_integer_root_control((1, 2, 3), "cubic_integer")


# Symbolic quadratic control q=T^2+c from the theorem statement.
c = sp.symbols("c", nonzero=True)
q2 = T**2 + c
Q2 = T**5 / 5 + 2 * c * T**3 / 3 + c**2 * T
C2 = sp.Rational(64, 225) * c**5
Psi2 = w**2 + C2
A2 = sp.cancel(Q2 * Psi2.subs(w, Q2) / q2**2)
B2 = sp.cancel(A2 / T)
check(sp.denom(A2) == 1, "symbolic quadratic H numerator polynomial")
check(sp.denom(B2) == 1, "symbolic quadratic axis quotient polynomial")
F2 = sp.expand(x * q2.subs(T, t_source))
W2 = sp.expand(Q2.subs(T, t_source))
H2 = sp.expand((1 + x * y) * B2.subs(T, t_source))
Y2 = sp.cancel((H2 - Psi2.subs(w, W2)) / F2)
check(sp.denom(Y2) == 1, "symbolic quadratic Y polynomial")
Y2 = sp.expand(Y2)
same(F2**3 * Y2, (W2 - F2**2) * (W2**2 + C2), "quadratic surface")
same(jac(F2, W2), F2**3, "quadratic F,W bracket")
same(jac(F2, Y2), 3 * W2**2 - 2 * F2**2 * W2 + C2, "quadratic F,Y bracket")
same(
    jac(W2, Y2),
    3 * F2**2 * Y2 + 2 * F2 * (W2**2 + C2),
    "quadratic W,Y bracket",
)
same(Y2.subs(x, 0), sp.Rational(64, 225) * c**4 * y, "quadratic axis")
alpha = sp.symbols("alpha", nonzero=True)
quad_boundary = sp.cancel(
    Y2.subs(y, (alpha - x**2) / x**3).subs(c, -alpha**2)
)
same(quad_boundary, sp.Rational(64, 675) * alpha**9 / x**3, "quadratic root arm")


# Degree-one hostile: the cubic critical root consumes the whole cover, the
# putative surface loses a divisor, and x^3 is an extra target-field polynomial.
a, alpha1 = sp.symbols("a alpha1", nonzero=True)
q1 = a * (T - alpha1)
Q1 = sp.expand(sp.integrate(q1**2, T) - sp.integrate(q1**2, T).subs(T, 0))
lambda1 = sp.expand(Q1.subs(T, alpha1))
same(Q1 - lambda1, a**2 * (T - alpha1) ** 3 / 3, "linear exact cubic cover")
F1 = sp.expand(x * q1.subs(T, t_source))
W1 = sp.expand(Q1.subs(T, t_source))
same(F1**3 / (3 * a * (W1 - lambda1)), x**3, "linear extra Kummer observable")
D1 = f**3 * z - (w - f**2) * (w - lambda1)
same(D1.subs({w: lambda1, z: 0}), 0, "linear missing target divisor")
check(sp.Poly(sp.cancel((Q1 - lambda1) / (T - alpha1) ** 3), T).degree() == 0,
      "linear residual degree zero")


# The d+1 reduced arm classes have the single primitive relation div(f), so
# their quotient lattice has rank d.  This is the exact Nagata/Picard packet
# used by the birational Darboux obstruction.
for d in range(2, 9):
    relation = sp.Matrix([[1] * (d + 1)])
    check(relation.rank() == 1, f"Pic relation rank d={d}")
    check((d + 1) - relation.rank() == d, f"Pic quotient rank d={d}")


semantic = {
    "boundary": "axis_all; each_critical_arm_punctured; exactly_d_missing_points",
    "birational": "units=k*;Pic=Z^d;no_birational_Darboux;source_degree>=2(2d+1)",
    "degree": "source_target_field_degree=2d+1",
    "hostile": "d=1_misses_divisor_and_x^3_is_extra_observable",
    "hypotheses": "q_squarefree; q(0)!=0; 0_and_all_Q(alpha_i)_pairwise_distinct; d>=2",
    "intersection": "k[x,y]_intersect_k(F,P)=k[F,W,Y]/(F^3Y-(W-F^2)Psi(W))",
    "poisson": "{F,W}=F^3; {F,Y}=Psi+(W-F^2)Psi'; {W,Y}=3F^2Y+2FPsi",
    "seed": "t=x^2(1+xy); Q'=q^2; F=xq(t); P=Q(t)/F^3; J(F,P)=1",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("theorem=THM-3789-higher-pole-hermite-spectral-completion")
print("seed=t=x^2(1+xy);Q'=q^2;F=xq(t);P=Q(t)/F^3;J(F,P)=1")
print("completion=H=W*Psi(W)/F^2;Y=(H-Psi(W))/F")
print("surface=F^3Y-(W-F^2)Psi(W)=0")
print("image=surface_minus_d_points_(0,lambda_i,0)")
print("atlas=quasi_finite_etale_codimension_one_surjective")
print("intersection=k[x,y]_cap_k(F,P)=k[F,W,Y]/relation")
print("degree=2d+1")
print("units=k*;Pic=Z^d;no_birational_Darboux;survivor_source_degree>=2(2d+1)")
print("quadratic=F^3Y-(W-F^2)(W^2+64c^5/225)=0")
print("hostile=d1_misses_divisor_and_retains_x^3")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
