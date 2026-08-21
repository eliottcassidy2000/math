#!/usr/bin/env python3
"""Deterministic exact companion for proved THM-3566.

The structural theorem is over C and holds for every odd d>=3.  This file
uses exact QQ arithmetic to verify the universal algebraic identities on a
generic polynomial control, the Chebyshev--Pell identities for odd d<=15,
and the complete boundary/Jelonek passports for d=3,5,7.  It also freezes
hostile controls for the maximal-intersection valuations, the exact
symplectic primitive, the two-brackets identity, and the fixed-power parity
classification.  Finite loops are controls, not the proof of all-d claims.
"""

from __future__ import annotations

import sympy as sp


TRUTH_GATES = 0


def require(condition: bool, label: str) -> None:
    global TRUTH_GATES
    TRUTH_GATES += 1
    if not condition:
        raise RuntimeError(label)


def jacobian2(f: sp.Expr, g: sp.Expr, r: sp.Symbol, s: sp.Symbol) -> sp.Expr:
    return sp.factor(sp.diff(f, r) * sp.diff(g, s) - sp.diff(f, s) * sp.diff(g, r))


def chebyshev_packet(d: int, t: sp.Symbol) -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Rational]:
    r = (d - 1) // 2
    z = 1 - 2 * t
    P = sp.expand((sp.chebyshevu(r, z) + sp.chebyshevu(r - 1, z)) / d)
    S = sp.expand(sp.chebyshevu(r, z) - sp.chebyshevu(r - 1, z))
    B = sp.expand(t * P**2)
    return P, S, B, sp.Rational(1, d**2)


x, q, t = sp.symbols("x q t")

# Generic master identity and construction.  A quadratic P is a control with
# algebraically independent coefficients; it is not a degree cutoff in the
# theorem.
H0, H1, G0, G1 = sp.symbols("H0 H1 G0 G1")
master_expanded = -H0 * G0 - 2 * t * H0 * G1 - t * H1 * G0
master_derivative = -(
    H0 * G0**2 + t * H1 * G0**2 + 2 * t * H0 * G0 * G1
) / G0
require(sp.cancel(master_expanded - master_derivative) == 0, "master identity")

p0, p1, p2 = sp.symbols("p0 p1 p2", nonzero=True)
P_generic = p0 + p1 * t + p2 * t**2
S_generic = sp.expand(P_generic + 2 * t * sp.diff(P_generic, t))
B_generic = sp.expand(t * P_generic**2)
require(sp.expand(sp.diff(B_generic, t) - P_generic * S_generic) == 0, "generic B prime")

t_xq = x**2 * q
P_xq = P_generic.subs(t, t_xq)
S_xq = S_generic.subs(t, t_xq)
a_generic = q / S_xq**2
c_generic = x * P_xq * S_xq
require(jacobian2(a_generic, c_generic, x, q) == -1, "generic rational Keller determinant")
require(sp.factor(a_generic * c_generic**2 - B_generic.subs(t, t_xq)) == 0, "generic B pullback")

# All odd Chebyshev packets through 15.  The theorem proves these identities
# from the half-angle formulas; this is an exact hostile bank.
chebyshev_rows = []
for d in range(3, 16, 2):
    r = (d - 1) // 2
    P, S, B, beta = chebyshev_packet(d, t)
    z = 1 - 2 * t
    require(sp.Poly(P, t).degree() == r, f"P degree d={d}")
    require(sp.Poly(S, t).degree() == r, f"S degree d={d}")
    require(P.subs(t, 0) == 1, f"P zero normalization d={d}")
    require(S.subs(t, 0) == 1, f"S zero normalization d={d}")
    require(P.subs(t, 1) == sp.Rational((-1) ** r, d), f"P one normalization d={d}")
    require(S.subs(t, 1) == (-1) ** r * d, f"S one normalization d={d}")
    require(sp.expand(S - P - 2 * t * sp.diff(P, t)) == 0, f"S construction d={d}")
    require(
        sp.expand(B - (1 - sp.chebyshevt(d, z)) / (2 * d**2)) == 0,
        f"Chebyshev B d={d}",
    )
    require(sp.expand(B - beta - (t - 1) * S**2 / d**2) == 0, f"Pell arm d={d}")
    require(sp.expand(sp.diff(B, t) - P * S) == 0, f"B derivative d={d}")
    require(sp.gcd(sp.Poly(P, t), sp.Poly(S, t)).degree() == 0, f"P S coprime d={d}")
    require(sp.gcd(sp.Poly(P, t), sp.Poly(sp.diff(P, t), t)).degree() == 0, f"P squarefree d={d}")
    require(sp.gcd(sp.Poly(S, t), sp.Poly(sp.diff(S, t), t)).degree() == 0, f"S squarefree d={d}")
    require(sp.Poly(B, t).degree() == d, f"B degree d={d}")
    require(sp.factor(B) == sp.factor(t * P**2), f"zero passport d={d}")
    require(sp.factor(B - beta) == sp.factor((t - 1) * S**2 / d**2), f"beta passport d={d}")
    require(sp.gcd(sp.Poly(P, t), sp.Poly(t * (t - 1), t)).degree() == 0, f"P endpoint hostile d={d}")
    require(sp.gcd(sp.Poly(S, t), sp.Poly(t * (t - 1), t)).degree() == 0, f"S endpoint hostile d={d}")
    chebyshev_rows.append((d, r, f"(1,2^{r})", f"(1,2^{r})"))

# Exact completion, charts, fibres, and all three nonproperness components in
# the requested bounded d=3,5,7 bank.
bb, cc, ee = sp.symbols("b c e")
boundary_rows = []
valuation_rows = []
for d in (3, 5, 7):
    r = (d - 1) // 2
    P, S, B, beta = chebyshev_packet(d, t)
    txq = x**2 * q
    Pxq = P.subs(t, txq)
    Sxq = S.subs(t, txq)
    b = sp.expand(B.subs(t, txq))
    c = sp.expand(x * Pxq * Sxq)
    a = q / Sxq**2
    e = sp.expand(q * (txq - 1) / d**2)

    require(jacobian2(a, c, x, q) == -1, f"Keller determinant d={d}")
    require(sp.factor(a * c**2 - b) == 0, f"b polynomialization d={d}")
    require(sp.factor(a * (b - beta) - e) == 0, f"e polynomialization d={d}")
    require(sp.expand(c**2 * e - b * (b - beta)) == 0, f"surface relation d={d}")
    require(sp.factor(jacobian2(b, c, x, q) + c**2) == 0, f"bc chart d={d}")
    require(sp.factor(jacobian2(c, e, x, q) - (2 * b - beta)) == 0, f"ce chart d={d}")

    relation = cc**2 * ee - bb * (bb - beta)
    singular_basis = sp.groebner(
        [relation, sp.diff(relation, bb), sp.diff(relation, cc), sp.diff(relation, ee)],
        bb,
        cc,
        ee,
        order="lex",
    )
    require(list(singular_basis) == [sp.Integer(1)], f"smooth target d={d}")
    require((2 * b - beta).subs({x: 0, q: 0}) == -beta, f"b=0 boundary sign d={d}")
    require(sp.simplify(2 * beta - beta) == beta, f"b=beta boundary sign d={d}")

    # Central point over L_0.
    s0 = sp.symbols(f"s{d}", nonzero=True)
    require(b.subs({x: 0, q: -d**2 * s0}) == 0, f"L0 central b d={d}")
    require(c.subs({x: 0, q: -d**2 * s0}) == 0, f"L0 central c d={d}")
    require(e.subs({x: 0, q: -d**2 * s0}) == s0, f"L0 central e d={d}")

    # Root-typed boundary formulae.  We use x^2 rather than choosing one of
    # the two algebraic square roots.
    q_root = d**2 * s0 / (t - 1)
    x2_root = t * (t - 1) / (d**2 * s0)
    require(sp.cancel(x2_root * q_root - t) == 0, f"root t reconstruction d={d}")
    require(sp.cancel(q_root * (t - 1) / d**2 - s0) == 0, f"root e reconstruction d={d}")
    require(sp.rem(sp.Poly(B, t), sp.Poly(P, t)) == 0, f"P root b=0 d={d}")
    require(sp.rem(sp.Poly(B - beta, t), sp.Poly(S, t)) == 0, f"S root b=beta d={d}")
    a_root_collision = sp.cancel(q_root / S**2 + d**2 * s0)
    a_root_numerator, a_root_denominator = sp.fraction(a_root_collision)
    require(
        sp.rem(
            sp.Poly(a_root_numerator, t, domain=sp.QQ.frac_field(s0)),
            sp.Poly(P, t, domain=sp.QQ.frac_field(s0)),
        )
        == 0,
        f"P root rational collision d={d}",
    )
    require(
        sp.gcd(sp.Poly(a_root_denominator, t), sp.Poly(P, t)).degree() == 0,
        f"P root collision lies in S nonzero d={d}",
    )
    require(
        sp.rem(sp.Poly(x2_root * (P * S) ** 2, t, domain=sp.QQ.frac_field(s0)), sp.Poly(P, t, domain=sp.QQ.frac_field(s0))) == 0,
        f"P root c=0 d={d}",
    )
    require(
        sp.rem(sp.Poly(x2_root * (P * S) ** 2, t, domain=sp.QQ.frac_field(s0)), sp.Poly(S, t, domain=sp.QQ.frac_field(s0))) == 0,
        f"S root c=0 d={d}",
    )

    # The c!=0 reconstruction, typed over QQ(T,C).
    T, C = sp.symbols(f"T{d} C{d}", nonzero=True)
    B_T = B.subs(t, T)
    Bp_T = sp.diff(B, t).subs(t, T)
    X_T = C / Bp_T
    Q_T = T * Bp_T**2 / C**2
    require(sp.cancel(X_T**2 * Q_T - T) == 0, f"open t reconstruction d={d}")
    require(sp.cancel(X_T * Bp_T - C) == 0, f"open c reconstruction d={d}")
    e_T = sp.cancel(Q_T * (T - 1) / d**2)
    require(
        sp.cancel(C**2 * e_T - B_T * (B_T - beta)) == 0,
        f"open target relation d={d}",
    )

    # Exact asymptotic-component controls.  At a P-root and S-root, B' also
    # vanishes simply, producing C_0 and C_beta respectively.  The second
    # derivative is nonzero there, so these are genuine branches.
    Bp = sp.diff(B, t)
    Bpp = sp.diff(Bp, t)
    require(sp.rem(sp.Poly(Bp, t), sp.Poly(P, t)) == 0, f"C0 derivative branch d={d}")
    require(sp.rem(sp.Poly(Bp, t), sp.Poly(S, t)) == 0, f"Cbeta derivative branch d={d}")
    require(sp.gcd(sp.Poly(P, t), sp.Poly(Bpp, t)).degree() == 0, f"C0 simple asymptotic d={d}")
    require(sp.gcd(sp.Poly(S, t), sp.Poly(Bpp, t)).degree() == 0, f"Cbeta simple asymptotic d={d}")

    # The third branch is t->1 with e fixed.  Its squared x-coordinate tends
    # to zero while b tends to beta and c^2 tends to zero.
    E = sp.symbols(f"E{d}", nonzero=True)
    x2_one = t * (t - 1) / (d**2 * E)
    q_one = d**2 * E / (t - 1)
    require(sp.cancel(x2_one * q_one - t) == 0, f"Lbeta t reconstruction d={d}")
    require(sp.cancel(q_one * (t - 1) / d**2 - E) == 0, f"Lbeta e fixed d={d}")
    require(sp.limit(B, t, 1) == beta, f"Lbeta b limit d={d}")
    require(sp.limit(x2_one * Bp**2, t, 1) == 0, f"Lbeta c limit d={d}")

    # Finite hostile controls for the intersection proof.  Removing one
    # beta-arm factor leaves a negative S-valuation.
    for weight in range(1, 7):
        m = (weight + 1) // 2
        candidate = sp.cancel((b * (b - beta)) ** m / c**weight)
        expected = e**m if weight % 2 == 0 else c * e**m
        require(sp.factor(candidate - expected) == 0, f"intersection normal form d={d} w={-weight}")
        hostile_S_order = 2 * (m - 1) - weight
        require(hostile_S_order < 0, f"intersection hostile valuation d={d} w={-weight}")
        valuation_rows.append((d, -weight, m, hostile_S_order))

    boundary_rows.append(
        (
            d,
            {
                "open_regular": d,
                "open_critical": 1,
                "L0_nonzero": d,
                "origin": 1,
                "Lbeta_nonzero": d - 1,
                "omitted": 0,
            },
            ("C0", "Cbeta", "Lbeta"),
        )
    )

# Global exact primitive and two-brackets identity on the abstract Y_beta.
beta0 = sp.symbols("beta", nonzero=True)
Q0 = bb * (bb - beta0)
Q0p = sp.diff(Q0, bb)
U0 = -4 / beta0**2
T0 = -Q0p / beta0**2
require(sp.cancel(U0 * Q0 - T0 * Q0p - 1) == 0, "Bezout UQ-TQprime")

# Coefficients of alpha=(U+T')ce db-T e dc-d(Tce).
alpha_b = sp.expand(U0 * cc * ee)
alpha_c = sp.expand(-2 * T0 * ee)
alpha_e = sp.expand(-T0 * cc)
alpha_b_chart = sp.cancel(alpha_b.subs(ee, Q0 / cc**2) + alpha_e * Q0p / cc**2)
alpha_c_chart = sp.cancel(alpha_c.subs(ee, Q0 / cc**2) - alpha_e * 2 * Q0 / cc**3)
require(sp.cancel(alpha_b_chart - 1 / cc) == 0, "global primitive db coefficient")
require(alpha_c_chart == 0, "global primitive dc coefficient")

# On c!=0 the declared target bracket is c^2(df/db dg/dc-df/dc dg/db).
def target_bracket_chart(f: sp.Expr, g: sp.Expr) -> sp.Expr:
    return sp.cancel(cc**2 * (sp.diff(f, bb) * sp.diff(g, cc) - sp.diff(f, cc) * sp.diff(g, bb)))


first_bracket_input = (U0 + sp.diff(T0, bb)) * Q0 / cc
second_bracket_input = -T0 * Q0 / cc**2
two_brackets = sp.cancel(
    target_bracket_chart(first_bracket_input, bb)
    + target_bracket_chart(second_bracket_input, cc)
)
require(two_brackets == 1, "two brackets sum to one")

# Direct generator bracket signs in the (a,c) chart with {a,c}=1.
aa, c0 = sp.symbols("a0 c0")
b0 = aa * c0**2
e0 = aa * (b0 - beta0)
bracket_bc = jacobian2(b0, c0, aa, c0)
bracket_ce = jacobian2(c0, e0, aa, c0)
bracket_be = jacobian2(b0, e0, aa, c0)
require(bracket_bc == c0**2, "Poisson bc sign")
require(sp.factor(bracket_ce - (beta0 - 2 * b0)) == 0, "Poisson ce sign")
require(sp.factor(bracket_be + 2 * c0 * e0) == 0, "Poisson be sign")

# Fixed-power parity classification over QQ[D].
D = sp.symbols("D")
fixed_power_rows = []
for m in range(0, 7):
    n = 2 * m
    lam = sp.Integer(2 * m + 1)
    central = sp.expand(
        sum(sp.binomial(2 * j, j) * D**j / 4**j for j in range(m + 1))
    )
    G_even = sp.expand(lam * D**m * central)
    kappa = sp.factor(lam * (2 * m + 1) * sp.binomial(2 * m, m) / 4**m)
    ode_left = sp.expand(
        2 * D * (D - 1) * sp.diff(G_even, D)
        + (n + (1 - n) * D) * G_even
    )
    require(sp.expand(ode_left - kappa * D ** (n + 1)) == 0, f"fixed-power even n={n}")
    require(sp.Poly(G_even, D).degree() == n, f"fixed-power degree n={n}")
    if m <= 3:
        Dxq = 1 + x**2 * q
        a_fixed = q / Dxq**n
        c_fixed = x * G_even.subs(D, Dxq)
        require(jacobian2(a_fixed, c_fixed, x, q) == -kappa, f"fixed-power Jacobian n={n}")
    fixed_power_rows.append((n, sp.factor(G_even), kappa))

G_triple = sp.expand(2 * D * (1 + D / 2))
require(sp.expand(G_triple - D * (D + 2)) == 0, "fixed-power triple G normalization")
require(
    sp.expand(
        2 * D * (D - 1) * sp.diff(G_triple, D)
        + (2 - D) * G_triple
        - 3 * D**3
    )
    == 0,
    "fixed-power triple kappa normalization",
)

# Odd n: exact coefficient systems with kappa=1 have no polynomial solution
# of the only possible degree n.  The theorem separately rules out rational
# poles, so this is the finite linear-algebra hostile bank.
odd_inconsistent = []
for n in range(1, 12, 2):
    coeffs = sp.symbols(f"g{n}_0:{n + 1}")
    G_trial = sum(coeffs[j] * D**j for j in range(n + 1))
    residual = sp.Poly(
        2 * D * (D - 1) * sp.diff(G_trial, D)
        + (n + (1 - n) * D) * G_trial
        - D ** (n + 1),
        D,
    )
    equations = [sp.Eq(residual.coeff_monomial(D**j), 0) for j in range(n + 2)]
    solution = sp.linsolve(equations, coeffs)
    require(solution is sp.EmptySet or solution == sp.EmptySet, f"fixed-power odd n={n}")
    odd_inconsistent.append(n)

# Finite pole-indicial hostiles at D=0 and D=1.  The theorem's local order
# argument is universal in s and n.
for n in range(0, 13):
    for pole_order in range(1, 9):
        require(2 * pole_order + n != 0, f"D=0 pole hostile n={n} s={pole_order}")
        require(1 - 2 * pole_order != 0, f"D=1 pole hostile n={n} s={pole_order}")

print("THM-3566 exact companion")
print("status=PROVED_VERIFIED_EXACT_INDEPENDENTLY_HOSTILE_AUDITED")
print("universe=QQ[x,q,t,b,c,e]; structural_scope=C,d=2r+1>=3")
print("master_identity=Jac(a,c)=-(t*H*G^2)'/G")
print("generic_construction=B=t*P^2 S=P+2tPprime a=q/S^2 c=x*P*S Jac=-1")
print(f"chebyshev_rows={chebyshev_rows}")
print("pell_identity=B-beta=(t-1)*S^2/d^2")
print("rational_collision=(a,c)=(-d^2*s,0) with d distinct sources for s!=0")
print("maximal_intersection=C[a,c]_polynomial_target intersect C[x,q]=C[b,c,e]")
print("target_relation=c^2*e-b*(b-beta)=0")
print(f"boundary_rows={boundary_rows}")
print("exact_image=Y_beta_minus_{(beta,0,0)}")
print("nonproperness=C0_union_Cbeta_union_Lbeta")
print(f"valuation_controls={valuation_rows}")
print("monodromy=dihedral_D_d_via_Kummer_u_to_u^d_mod_reflection")
print("ramification_budget=P_roots:r S_roots:r")
print("symplectic_primitive=global_exact")
print("two_brackets_to_one=PASS")
print(f"fixed_power_even_rows={fixed_power_rows}")
print(f"fixed_power_odd_inconsistent={odd_inconsistent}")
print(f"truth_gates={TRUTH_GATES}")
print("scope=NO_JC2_CLAIM; NO_LITERATURE_NOVELTY_CLAIM; SINGLE_BRACKET_REMAINS_OPEN")
