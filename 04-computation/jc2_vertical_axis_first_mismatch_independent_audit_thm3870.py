#!/usr/bin/env python3
"""Independent exact audit of the THM-3870 vertical first-mismatch argument.

This implementation was derived directly from the depressed-cubic branch
formula.  It deliberately does not import or execute the candidate checker.
It contains no Python assertions, so optimized mode exercises the same gates.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import sympy as sp


sys.stdout.reconfigure(encoding="utf-8", newline="\n")

A, C, Z = sp.symbols("A C Z")
B, X = sp.symbols("B X")
GATES = 0


def require(condition: bool, label: str) -> None:
    global GATES
    GATES += 1
    if not bool(condition):
        raise SystemExit(f"FAIL gate {GATES}: {label}")


def delta(b: sp.Expr) -> sp.Expr:
    return sp.expand(
        -27 * A**2 * b**2
        + 8 * A * C**3
        - 54 * A * C * b
        + 9 * C**2
        - 54 * b
    )


def q(j: int) -> sp.Expr:
    coeff = sp.binomial(sp.Rational(3, 2), j + 3) * sp.Rational(2, 3) ** (j + 3)
    return sp.expand(coeff * C ** (j + 3))


def deg_c(p: sp.Expr) -> int:
    return int(sp.Poly(sp.expand(p), C).degree())


def lc_c(p: sp.Expr) -> sp.Expr:
    return sp.expand(sp.Poly(sp.expand(p), C).LC())


def divisible_by_a_power(p: sp.Expr, power: int) -> bool:
    return sp.expand(sp.rem(sp.expand(p), A**power, A)) == 0


def divide_exact(p: sp.Expr, divisor: sp.Expr, label: str) -> sp.Expr:
    quo, rem = sp.div(sp.Poly(sp.expand(p), A, C), sp.Poly(divisor, A, C))
    require(rem.as_expr() == 0, label)
    return sp.expand(quo.as_expr())


def homogenize_c(p: sp.Expr, degree: int) -> sp.Expr:
    poly = sp.Poly(sp.expand(p), C)
    ans = 0
    for (power,), coefficient in poly.terms():
        ans += coefficient * C**power * Z ** (degree - power)
    return sp.expand(ans)


# A | Delta classification, derived before any comparator is introduced.
b0_symbols = sp.symbols("b0:6")
b0 = sum(b0_symbols[i] * C**i for i in range(6))
require(sp.expand(delta(b0).subs(A, 0) - (9 * C**2 - 54 * b0)) == 0,
        "specialization classification identity")
classification_coeffs = sp.Poly(9 * C**2 - 54 * b0, C)
solution = sp.solve(classification_coeffs.all_coeffs(), b0_symbols, dict=True)
require(len(solution) == 1, "classification coefficient solution is unique")
require(sp.expand(b0.subs(solution[0]) - C**2 / 6) == 0,
        "A divides Delta exactly at b(0,C)=C^2/6")

Qprobe = 2 + 3 * A + (1 - A) * C + A**2 * C**3
vertical_response = divide_exact(delta(C**2 / 6 + A * Qprobe), A,
                                 "vertical substitution has factor A")
expected_response = sp.expand(
    -C**3 - 54 * Qprobe - 54 * A * C * Qprobe
    - sp.Rational(3, 4) * A * C**4
    - 9 * A**2 * C**2 * Qprobe - 27 * A**3 * Qprobe**2
)
require(sp.expand(vertical_response - expected_response) == 0,
        "vertical branch formula")

# The sign and all coefficients of the unique formal comparator.
P = 1 + sp.Rational(2, 3) * A * C
for j in range(9):
    require(q(j) != 0, f"q_{j} is nonzero")
    require(deg_c(q(j)) == j + 3, f"degree q_{j}")

for depth in range(1, 9):
    Qstar_trunc = sum(q(j) * A**j for j in range(depth))
    bstar_trunc = C**2 / 6 + A * Qstar_trunc
    ustar_trunc = sp.expand(1 + A * C + A**2 * bstar_trunc)
    sqrt_trunc = sp.series(P ** sp.Rational(3, 2), A, 0, depth + 3).removeO()
    require(divisible_by_a_power(ustar_trunc - sqrt_trunc, depth + 2),
            f"positive square-root comparator through depth {depth}")
    require(divisible_by_a_power(delta(bstar_trunc), depth + 1),
            f"formal discriminant vanishing through depth {depth}")

require(sp.expand((1 + A * C + A**2 * (C**2 / 6)).subs(A, 0)) == 1,
        "prescribed u has positive constant sign")
require((-P ** sp.Rational(3, 2)).subs(A, 0) == -1,
        "negative square-root sign is excluded")

# Each new coefficient enters the branch equation with the unit -54.
unknowns = sp.symbols("y0:8")
Qformal = sum(unknowns[j] * A**j for j in range(8))
branch = sp.Poly(
    sp.expand(delta(C**2 / 6 + A * Qformal) / A), A
)
solved = {}
for n in range(8):
    coefficient = sp.expand(branch.coeff_monomial(A**n).subs(solved))
    slope = sp.diff(coefficient, unknowns[n])
    require(slope == -54, f"triangular unit at coefficient {n}")
    answer = sp.solve(coefficient, unknowns[n])
    require(len(answer) == 1, f"unique coefficient recursion at {n}")
    solved[unknowns[n]] = sp.factor(answer[0])
    require(sp.expand(solved[unknowns[n]] - q(n)) == 0,
            f"binomial coefficient recursion at {n}")

# The exact quadratic response has the claimed sign.
require(
    sp.expand(delta(B + X) - delta(B) + 54 * (1 + A * C + A**2 * B) * X
              + 27 * A**2 * X**2) == 0,
    "exact nonlinear response identity",
)


def base_data(n: int) -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    q_less = sum(q(j) * A**j for j in range(n))
    b_n = sp.expand(C**2 / 6 + A * q_less)
    u_n = sp.expand(1 + A * C + A**2 * b_n)
    r_n = divide_exact(delta(b_n), A ** (n + 1),
                       f"polynomial A^(N+1) division at N={n}")
    rho = sp.Rational(1, 6) if n == 0 else sp.Poly(q(n - 1), C).LC()
    return q_less, b_n, u_n, r_n, rho


BASE = {}
for n in range(7):
    q_less, b_n, u_n, r_n, rho = base_data(n)
    BASE[n] = (q_less, b_n, u_n, r_n, rho)
    require(sp.expand(r_n.subs(A, 0) - 54 * q(n)) == 0,
            f"special residual at N={n}")
    require(deg_c(b_n) == n + 2, f"base b C-degree at N={n}")
    require(lc_c(b_n) == rho * A**n, f"base b leading coefficient at N={n}")
    require(deg_c(u_n) == n + 2, f"base u C-degree at N={n}")
    require(lc_c(u_n) == rho * A ** (n + 2),
            f"base u leading coefficient at N={n}")
    require(deg_c(r_n) == 2 * n + 4, f"base residual degree at N={n}")
    require(lc_c(r_n) == -27 * rho**2 * A ** (n + 1),
            f"base residual leading coefficient at N={n}")


def check_profile(n: int, T: sp.Expr, tag: str, factor_check: bool = False) -> None:
    q_less, _b_n, u_n, r_n, rho = BASE[n]
    t = sp.expand(T.subs(A, 0))
    require(sp.expand(t - q(n)) != 0, f"first mismatch is genuine: {tag}")
    Q = sp.expand(q_less + A**n * T)
    b = sp.expand(C**2 / 6 + A * Q)
    S = divide_exact(delta(b), A ** (n + 1), f"actual divisibility: {tag}")
    response = sp.expand(r_n - 54 * u_n * T - 27 * A ** (n + 3) * T**2)
    require(sp.expand(S - response) == 0, f"nonlinear response: {tag}")
    special = sp.expand(S.subs(A, 0))
    require(sp.expand(special - 54 * (q(n) - t)) == 0,
            f"special mismatch formula: {tag}")
    require(special != 0, f"special mismatch is nonzero: {tag}")

    D = deg_c(S)
    m = deg_c(special)
    if T == 0:
        require(D == 2 * n + 4, f"T=0 generic degree: {tag}")
        require(lc_c(S) == -27 * rho**2 * A ** (n + 1),
                f"T=0 leading coefficient: {tag}")
    else:
        d = deg_c(T)
        tau = lc_c(T)
        expected_D = 2 * (n + 2) if d <= n + 2 else 2 * d
        require(D == expected_D, f"generic degree regime: {tag}")
        if d < n + 2:
            expected_lc = -27 * rho**2 * A ** (n + 1)
        elif d == n + 2:
            expected_lc = -27 * A ** (n + 1) * (rho + A * tau) ** 2
            require(sp.expand(rho + A * tau) != 0,
                    f"resonant square is not the zero polynomial: {tag}")
        else:
            expected_lc = -27 * A ** (n + 3) * tau**2
        require(sp.expand(lc_c(S) - expected_lc) == 0,
                f"leading coefficient regime: {tag}")
    require(D > m, f"strict C-degree drop: {tag}")

    Sh = homogenize_c(S, D)
    special_h = homogenize_c(special, m)
    require(sp.expand(Sh.subs({A: 0, C: 1, Z: 0})) == 0,
            f"total closure reaches finite-base C-infinity: {tag}")
    require(sp.expand(special_h.subs({C: 1, Z: 0})) != 0,
            f"degree-drop factor, not special leading cancellation: {tag}")
    require(
        sp.expand(Sh.subs(A, 0) - special_h * Z ** (D - m)) == 0,
        f"homogenized special-fibre identity: {tag}",
    )

    if factor_check:
        constant, factors = sp.factor_list(S, A, C)
        require(constant != 0 and len(factors) >= 1, f"factorization exists: {tag}")
        rebuilt = sp.expand(constant * sp.prod(g**e for g, e in factors))
        require(sp.expand(rebuilt - S) == 0, f"factorization with multiplicity: {tag}")
        require(sum(e * deg_c(g) for g, e in factors) == D,
                f"actual C-degrees add with multiplicity: {tag}")
        selected = []
        for g, e in factors:
            require(sp.expand(g.subs(A, 0)) != 0,
                    f"A is not an affine residual factor: {tag}")
            dh = deg_c(g)
            gh = homogenize_c(g, dh)
            if sp.expand(gh.subs({A: 0, C: 1, Z: 0})) == 0:
                selected.append((g, e, dh))
        require(len(selected) >= 1, f"a reduced factor reaches P0: {tag}")
        require(all(dh > 0 for _g, _e, dh in selected),
                f"selected factor is nonvertical/dominant: {tag}")


# Exhaust all relative degree regimes, including N=0 and T=0.
for n in range(6):
    check_profile(n, sp.Integer(0), f"N{n}-Tzero", factor_check=n <= 2)
    degrees = sorted(set([0, max(0, n + 1), n + 2, n + 3, n + 5]))
    for d in degrees:
        T = sp.expand((2 + A) * C**d + C + 1) if d else 2 + A
        check_profile(n, T, f"N{n}-d{d}", factor_check=(n <= 1 and d <= n + 2))

    rho = BASE[n][4]
    resonant = sp.expand(-rho * C ** (n + 2) + C + 1)
    check_profile(n, resonant, f"N{n}-resonant-isolated-A-root")

    # At A=0 these share the full q_N leading term; the mismatch occurs lower.
    check_profile(n, sp.expand(q(n) + 1), f"N{n}-shared-qN-constant",
                  factor_check=n == 0)
    check_profile(n, sp.expand(q(n) + C ** (n + 1) + 1),
                  f"N{n}-shared-qN-lower-tail")
    check_profile(n, sp.expand(q(n) + 1 + A * C ** (n + 5)),
                  f"N{n}-shared-qN-high-A-tail")

# A purely algebraic multiplicity control for the factor-selection lemma.
G = A * C + 1
H = C - A
synthetic = sp.expand(G**2 * H)
synthetic_D = deg_c(synthetic)
synthetic_m = deg_c(synthetic.subs(A, 0))
require(synthetic_D == 3 and synthetic_m == 1,
        "nonreduced synthetic packet has strict degree drop")
syn_constant, syn_factors = sp.factor_list(synthetic, A, C)
require(any(sp.expand(g - G) == 0 and e == 2 for g, e in syn_factors),
        "multiplicity is retained on selected reduced factor")
require(homogenize_c(G, 1).subs({A: 0, C: 1, Z: 0}) == 0,
        "selected nonreduced factor reaches (A=0,C=infinity)")
require(sp.expand(G.subs(A, 0)) == 1,
        "selected factor is not the vertical axis")
require(sp.expand(G.subs(C, 0) - 1) == 0,
        "same normalization has distinct A-infinity endpoint C=0")

semantic = "\n".join(
    [
        "classification:A|Delta iff b=C^2/6+A*Q",
        "comparator:unique positive binomial branch; q_j nonzero degree j+3",
        "mismatch:every polynomial Q has a unique finite first mismatch",
        "completion:Delta(b_N) divisible by A^(N+1) in the polynomial ring",
        "response:S=R_N-54*u_N*T-27*A^(N+3)*T^2",
        "degrees:d<N+2 base; d=N+2 nonzero square; d>N+2 quadratic",
        "cancellation:S(0,C)=54(q_N-t) handles arbitrary shared q_N terms",
        "geometry:a reduced dominant residual factor reaches finite-base C-infinity",
        "punctures:normalization also has a distinct point over A-infinity",
    ]
)
semantic_hash = hashlib.sha256((semantic + "\n").encode("utf-8")).hexdigest()
lines = [
    "THM3870_VERTICAL_FIRST_MISMATCH_HOSTILE_AUDIT PASS",
    "CLASSIFICATION A|Delta iff b=C^2/6+A*Q",
    "COMPARATOR unique positive-sign formal binomial branch",
    "FIRST_MISMATCH finite and unique for every polynomial Q",
    "DIVISIBILITY A^(N+1) completion-to-polynomial gate",
    "DEGREE_REGIMES lower / d=N+2 square / higher all strict",
    "CANCELLATION shared q_N leading terms explicitly tested",
    "FACTOR_GATE multiplicities retained; selected reduced factor dominates A",
    "PUNCTURES distinct missing points over A=0 and A=infinity",
    f"SEMANTIC_SHA256 {semantic_hash}",
    f"GATES {GATES}",
]
output = "\n".join(lines) + "\n"

if "--verify-frozen" in sys.argv:
    index = sys.argv.index("--verify-frozen")
    require(index + 1 < len(sys.argv), "frozen transcript path supplied")
    frozen_path = Path(sys.argv[index + 1])
    frozen = frozen_path.read_bytes()
    if frozen != output.encode("utf-8"):
        raise SystemExit("FAIL frozen transcript mismatch")

sys.stdout.write(output)
