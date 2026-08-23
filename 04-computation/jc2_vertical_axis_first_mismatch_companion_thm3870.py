#!/usr/bin/env python3
"""Exact controls for the THM-3870 vertical-axis first-mismatch theorem.

This promoted proof companion checks the universal axis classification, the
formal binomial comparator, the exact nonlinear response, and every C-degree
regime of the first mismatch.
"""

from __future__ import annotations

import ast
import hashlib
import json
import sys
from pathlib import Path

import sympy as sp


sys.stdout.reconfigure(newline="\n")

A, C = sp.symbols("A C")
tau, eta = sp.symbols("tau eta")
CHECKS = 0
MAX_N = 5


def require(condition: object, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(f"{label}: {condition}")


def zero(expression: sp.Expr, label: str) -> None:
    require(sp.factor(expression) == 0, label)


def equal(left: object, right: object, label: str) -> None:
    require(left == right, label)


def beta(j: int) -> sp.Rational:
    return sp.factor(
        sp.binomial(sp.Rational(3, 2), j + 2)
        * sp.Rational(2, 3) ** (j + 2)
    )


def gamma(j: int) -> sp.Rational:
    return beta(j + 1)


def branch(profile: sp.Expr) -> sp.Expr:
    return sp.expand(
        -27 * A**2 * profile**2
        + 8 * A * C**3
        - 54 * A * C * profile
        + 9 * C**2
        - 54 * profile
    )


b_special = sp.symbols("b_special")
zero(branch(b_special).subs(A, 0) - (9 * C**2 - 54 * b_special),
     "axis restriction")
zero(sp.solve(sp.Eq(branch(b_special).subs(A, 0), 0), b_special)[0]
     - C**2 / 6, "axis classification")

Q_free = sp.symbols("Q_free")
axis_profile = C**2 / 6 + A * Q_free
axis_quotient = sp.cancel(branch(axis_profile) / A)
require(sp.denom(axis_quotient) == 1, "axis factor is polynomial")

formal_order = MAX_N + 5
b_star = sp.series(
    ((1 + sp.Rational(2, 3) * A * C) ** sp.Rational(3, 2)
     - 1 - A * C) / A**2,
    A,
    0,
    formal_order,
).removeO().expand()
zero(b_star.coeff(A, 0) - C**2 / 6, "formal axis coefficient")
for j in range(MAX_N + 2):
    zero(b_star.coeff(A, j + 1) - gamma(j) * C ** (j + 3),
         f"formal vertical quotient coefficient {j}")
    require(gamma(j) != 0, f"gamma {j} nonzero")

b_symbol, delta_symbol = sp.symbols("b_symbol delta_symbol")
u_symbol = 1 + A * C + A**2 * b_symbol
zero(
    branch(b_symbol + delta_symbol)
    - branch(b_symbol)
    + 54 * u_symbol * delta_symbol
    + 27 * A**2 * delta_symbol**2,
    "universal profile perturbation",
)

rows: list[tuple[object, ...]] = []
for N in range(MAX_N + 1):
    q = [gamma(j) * C ** (j + 3) for j in range(N + 1)]
    QN = sp.expand(sum(A**j * q[j] for j in range(N)))
    bN = sp.expand(C**2 / 6 + A * QN)
    RN = sp.cancel(branch(bN) / A ** (N + 1))
    require(sp.denom(RN) == 1, f"N={N} base residual polynomial")
    RN = sp.expand(RN)
    uN = sp.expand(1 + A * C + A**2 * bN)
    lead_previous = sp.Rational(1, 6) if N == 0 else gamma(N - 1)

    zero(RN.subs(A, 0) - 54 * q[N], f"N={N} base special residual")
    equal(sp.degree(RN, C), 2 * N + 4, f"N={N} base generic degree")
    zero(
        sp.LC(sp.Poly(RN, C))
        + 27 * A ** (N + 1) * lead_previous**2,
        f"N={N} base leading coefficient",
    )
    equal(sp.degree(uN, C), N + 2, f"N={N} u degree")
    zero(
        sp.LC(sp.Poly(uN, C)) - A ** (N + 2) * lead_previous,
        f"N={N} u leading coefficient",
    )

    S_zero = RN
    equal(sp.degree(S_zero, C), 2 * N + 4, f"N={N} T=0 generic degree")
    equal(sp.degree(S_zero.subs(A, 0), C), N + 3,
          f"N={N} T=0 special degree")
    require(2 * N + 4 > N + 3, f"N={N} T=0 strict drop")

    regimes = sorted({0, max(0, N + 1), N + 2, N + 3, N + 4})
    for d in regimes:
        T = tau * C**d + eta
        b = sp.expand(bN + A ** (N + 1) * T)
        S = sp.expand(RN - 54 * uN * T - 27 * A ** (N + 3) * T**2)
        zero(branch(b) - A ** (N + 1) * S,
             f"N={N},d={d} first-mismatch factorization")
        zero(S.subs(A, 0) - 54 * (q[N] - T),
             f"N={N},d={d} special residual")

        generic_degree = 2 * N + 4 if d <= N + 2 else 2 * d
        special_degree = max(N + 3, d)
        equal(sp.degree(S, C), generic_degree,
              f"N={N},d={d} generic degree")
        equal(sp.degree(S.subs(A, 0), C), special_degree,
              f"N={N},d={d} special degree")
        require(generic_degree > special_degree,
                f"N={N},d={d} strict degree drop")

        leading = sp.LC(sp.Poly(S, C))
        if d < N + 2:
            expected = -27 * A ** (N + 1) * lead_previous**2
        elif d > N + 2:
            expected = -27 * A ** (N + 3) * tau**2
        else:
            expected = -27 * A ** (N + 1) * (
                lead_previous + A * tau
            ) ** 2
        zero(leading - expected, f"N={N},d={d} leading coefficient")
        if d == N + 2:
            require(expected != 0, f"N={N} resonant square nonzero")

        rows.append((N, d, generic_degree, special_degree,
                     sp.factor(leading)))

source = Path(__file__).read_text(encoding="utf-8")
require(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "checker contains no inactive assert",
)

semantic = {
    "status": "PROVED_VERIFIED_EXACT",
    "axis": "A divides Delta iff b=C^2/6+A Q",
    "formal": "Q_*=sum gamma_j A^j C^(j+3)",
    "mismatch": "Q=Q_<N+A^N T; t!=gamma_N C^(N+3)",
    "factor": "Delta=A^(N+1)S",
    "special": "S(0,C)=54(q_N-t) nonzero",
    "degrees": "d<N+2:2N+4;d>N+2:2d;d=N+2:2N+4",
    "resonance": "-27A^(N+1)(gamma_(N-1)+A tau)^2",
    "geometry": "strict finite-base C-degree drop selects nonvertical companion",
    "scope": "vertical axis in depressed-cubic carrier only; JC2 open",
}
semantic_sha = hashlib.sha256(
    json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
).hexdigest()

output = "\n".join(
    [
        "experiment=THM3870-vertical-axis-first-mismatch-companion",
        "status=PROVED;VERIFIED_EXACT",
        "classification=A_divides_Delta_iff_b_equals_C2_over_6_plus_AQ",
        "formal_comparator=Q_star_sum_gamma_j_Aj_C_to_j_plus_3",
        "first_mismatch=Q_equals_Q_below_N_plus_A_to_N_T",
        "factorization=Delta_equals_A_to_N_plus_1_times_S",
        "special_residual=S_at_A0_equals_54_times_qN_minus_t_nonzero",
        "degree_regimes=d_below_Nplus2:2Nplus4;d_above:2d;resonance:2Nplus4",
        "resonance=nonzero_perfect_square_leading_coefficient",
        f"finite_rows={len(rows)};N=0..{MAX_N}",
        f"GATES={CHECKS}",
        f"semantic_sha256={semantic_sha}",
        "scope=vertical_axis_depressed_cubic_carrier_only;no_JC2_conclusion",
        "RESULT PASS",
    ]
) + "\n"

if "--verify-frozen" in sys.argv:
    index = sys.argv.index("--verify-frozen")
    require(index + 1 < len(sys.argv), "frozen transcript path supplied")
    frozen = Path(sys.argv[index + 1]).read_bytes()
    if frozen != output.encode("utf-8"):
        raise SystemExit("FAIL frozen transcript mismatch")

sys.stdout.write(output)
