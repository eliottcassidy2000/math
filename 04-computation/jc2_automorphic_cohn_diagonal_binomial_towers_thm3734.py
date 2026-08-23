#!/usr/bin/env python3
"""Exact companion for THM-3734's paired diagonal binomial towers."""

from __future__ import annotations

import ast
import hashlib
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


X, T, Z = sp.symbols("X T Z")


def curl(one_form: tuple[sp.Expr, sp.Expr]) -> sp.Expr:
    return sp.expand(sp.diff(one_form[0], T) - sp.diff(one_form[1], X))


def jac(first: sp.Expr, second: sp.Expr) -> sp.Expr:
    return sp.expand(
        sp.diff(first, X) * sp.diff(second, T)
        - sp.diff(first, T) * sp.diff(second, X)
    )


M0 = sp.Matrix(((4 * T**2, 2 * X * T - 1), (1 + 2 * X * T, X**2)))
gate(M0.det() == 1, "base Cohn determinant")

# Derive the two full diagonal closure PDEs before specializing a charge or
# a resonance.  These are proof guards against checking only the displayed
# particular solutions.
a_symbol = sp.symbols("a_symbol", nonzero=True)
h_function = sp.Function("h_function")(X, T)
A_symbol = a_symbol**2
N_symbol = M0 * sp.diag(a_symbol, 1 / a_symbol)
alpha_symbol = (N_symbol[0, 0], N_symbol[0, 1])
beta_symbol = (N_symbol[1, 0], N_symbol[1, 1])
lower_symbol = (
    beta_symbol[0] + h_function * alpha_symbol[0],
    beta_symbol[1] + h_function * alpha_symbol[1],
)
upper_symbol = (
    alpha_symbol[0] + h_function * beta_symbol[0],
    alpha_symbol[1] + h_function * beta_symbol[1],
)
lower_pde = (
    4 * A_symbol * T**2 * sp.diff(h_function, T)
    - (2 * X * T - 1) * sp.diff(h_function, X)
    + 2 * (4 * A_symbol - 1) * T * h_function
    + 2 * (A_symbol - 1) * X
)
upper_pde = (
    A_symbol * (1 + 2 * X * T) * sp.diff(h_function, T)
    - X**2 * sp.diff(h_function, X)
    + 2 * (A_symbol - 1) * X * h_function
    + 2 * (4 * A_symbol - 1) * T
)
gate(sp.simplify(a_symbol * curl(lower_symbol) - lower_pde) == 0,
     "universal lower diagonal closure PDE")
gate(sp.simplify(a_symbol * curl(upper_symbol) - upper_pde) == 0,
     "universal upper diagonal closure PDE")

semantic_rows: list[str] = []
tested_depths = range(1, 13)

for r in tested_depths:
    # Lower-row tower.
    A = sp.Rational(r + 1, 2 * r)
    a = sp.sqrt(A)
    K = sp.expand((1 + 2 * Z / r) ** r)
    H = sp.cancel((K - 1 - 2 * Z) / (4 * Z**2))
    phi = sp.cancel(
        r * ((1 + 2 * Z / r) ** (r + 1) - 1)
        / (2 * (r + 1) * Z)
    )
    gate(not sp.denom(sp.together(H)).has(Z),
         f"lower H polynomial at r={r}")
    gate(not sp.denom(sp.together(phi)).has(Z),
         f"lower divided difference polynomial at r={r}")
    gate(sp.expand(phi + Z * sp.diff(phi, Z) - K) == 0,
         f"lower discrete antiderivative at r={r}")
    gate(sp.expand(1 + 2 * Z + 4 * Z**2 * H - K) == 0,
         f"lower closure numerator at r={r}")
    gate(sp.Poly(phi, Z).degree() == r, f"lower profile degree at r={r}")
    gate(sp.Poly(K, Z).degree() == r, f"lower kernel degree at r={r}")
    if r == 1:
        gate(H == 0, "lower depth-one constant closure")
    else:
        gate(sp.Poly(H, Z).degree() == r - 2,
             f"lower exposed-parameter degree at r={r}")
    lower_phi_lead = sp.Poly(phi, Z).LC()
    lower_K_lead = sp.Poly(K, Z).LC()
    gate(sp.expand(lower_K_lead - (r + 1) * lower_phi_lead) == 0,
         f"lower leading ratio at r={r}")

    R = sp.diag(a, 1 / a)
    N = sp.expand(M0 * R)
    alpha = (N[0, 0], N[0, 1])
    beta = (N[1, 0], N[1, 1])
    h = sp.expand(X**2 * H.subs(Z, X * T))
    q = sp.expand(X * phi.subs(Z, X * T))
    Q = sp.expand(a * q)
    lower = (
        sp.expand(beta[0] + h * alpha[0]),
        sp.expand(beta[1] + h * alpha[1]),
    )
    gate(sp.simplify(N.det() - 1) == 0, f"lower right determinant at r={r}")
    gate(curl(lower) == 0, f"lower row closure at r={r}")
    gate(sp.simplify(sp.diff(Q, X) - lower[0]) == 0,
         f"lower potential X derivative at r={r}")
    gate(sp.simplify(sp.diff(Q, T) - lower[1]) == 0,
         f"lower potential T derivative at r={r}")
    lower_target = sp.Rational(4 * (r + 2), r + 1) * T
    gate(sp.simplify(curl(alpha) / a - lower_target) == 0,
         f"lower normalized debt at r={r}")
    gate(sp.expand(phi.subs(Z, 0) - 1) == 0,
         f"lower zero-fibre flank at r={r}")
    gate(sp.gcd(sp.Poly(phi, Z), sp.Poly(sp.diff(phi, Z), Z)).degree() == 0,
         f"lower cyclotomic flanks square-free at r={r}")
    gate(sp.expand(
        2 * (r + 1) * Z * phi
        - r * ((1 + 2 * Z / r) ** (r + 1) - 1)
    ) == 0, f"lower roots-of-unity fibre identity at r={r}")
    if r >= 2:
        n = r - 2
        gate(4 * A * (n + 2) - 2 * n - 6 == 0,
             f"lower finite-termination resonance at r={r}")

    # A hostile perturbation of the diagonal slope destroys the same closure.
    bad_A = A + sp.Rational(1, 17)
    bad_a = sp.sqrt(bad_A)
    bad_N = sp.expand(M0 * sp.diag(bad_a, 1 / bad_a))
    bad_alpha = (bad_N[0, 0], bad_N[0, 1])
    bad_beta = (bad_N[1, 0], bad_N[1, 1])
    bad_lower = (
        sp.expand(bad_beta[0] + h * bad_alpha[0]),
        sp.expand(bad_beta[1] + h * bad_alpha[1]),
    )
    gate(sp.simplify(curl(bad_lower)) != 0,
         f"lower resonance hostile control at r={r}")

    # Charge -2 is the only component which can hit the target T.  Its
    # general polynomial is T^2 F(XT), and the Hamiltonian equation reduces
    # to one profile ODE.  Its top coefficient can never vanish.
    for n in range(0, 9):
        coefficients = sp.symbols(f"lo_{r}_{n}_0:{n + 1}")
        F = sum(coefficients[i] * Z**i for i in range(n + 1))
        f = sp.expand(T**2 * F.subs(Z, X * T))
        profile = sp.expand(Z * phi * sp.diff(F, Z) + 2 * K * F)
        gate(
            sp.expand(jac(f, q) + T * profile.subs(Z, X * T)) == 0,
            f"lower charge ODE at r={r}, degree={n}",
        )
        predicted = sp.expand(
            (n + 2 * (r + 1)) * lower_phi_lead * coefficients[n]
        )
        actual = sp.Poly(profile, Z).coeff_monomial(Z ** (n + r))
        gate(sp.expand(actual - predicted) == 0,
             f"lower unique top edge at r={r}, degree={n}")
        equation = sp.Poly(profile + sp.Rational(4 * (r + 2), r + 1), Z)
        system = [coefficient for _, coefficient in equation.terms()]
        gate(sp.linsolve(system, coefficients) == sp.EmptySet,
             f"lower finite hostile cokernel at r={r}, degree={n}")

    # Upper-row dual tower.
    A_top = sp.Rational(r, 2 * (r + 1))
    a_top = sp.sqrt(A_top)
    K_top = sp.expand(-(1 - 2 * Z / r) ** r)
    H_top = sp.cancel((K_top + 1 - 2 * Z) / Z**2)
    psi = sp.cancel(
        -r * (1 - (1 - 2 * Z / r) ** (r + 1))
        / (2 * (r + 1) * Z)
    )
    gate(not sp.denom(sp.together(H_top)).has(Z),
         f"upper H polynomial at r={r}")
    gate(not sp.denom(sp.together(psi)).has(Z),
         f"upper divided difference polynomial at r={r}")
    gate(sp.expand(psi + Z * sp.diff(psi, Z) - K_top) == 0,
         f"upper discrete antiderivative at r={r}")
    gate(sp.expand(2 * Z - 1 + Z**2 * H_top - K_top) == 0,
         f"upper closure numerator at r={r}")
    gate(sp.Poly(psi, Z).degree() == r, f"upper profile degree at r={r}")
    gate(sp.Poly(K_top, Z).degree() == r, f"upper kernel degree at r={r}")
    if r == 1:
        gate(H_top == 0, "upper depth-one constant closure")
    else:
        gate(sp.Poly(H_top, Z).degree() == r - 2,
             f"upper exposed-parameter degree at r={r}")
    upper_psi_lead = sp.Poly(psi, Z).LC()
    upper_K_lead = sp.Poly(K_top, Z).LC()
    gate(sp.expand(upper_K_lead - (r + 1) * upper_psi_lead) == 0,
         f"upper leading ratio at r={r}")

    R_top = sp.diag(a_top, 1 / a_top)
    N_top = sp.expand(M0 * R_top)
    alpha_top = (N_top[0, 0], N_top[0, 1])
    beta_top = (N_top[1, 0], N_top[1, 1])
    h_top = sp.expand(T**2 * H_top.subs(Z, X * T))
    q_top = sp.expand(T * psi.subs(Z, X * T))
    Q_top = sp.expand(q_top / a_top)
    upper = (
        sp.expand(alpha_top[0] + h_top * beta_top[0]),
        sp.expand(alpha_top[1] + h_top * beta_top[1]),
    )
    gate(sp.simplify(N_top.det() - 1) == 0,
         f"upper right determinant at r={r}")
    gate(curl(upper) == 0, f"upper row closure at r={r}")
    gate(sp.simplify(sp.diff(Q_top, X) - upper[0]) == 0,
         f"upper potential X derivative at r={r}")
    gate(sp.simplify(sp.diff(Q_top, T) - upper[1]) == 0,
         f"upper potential T derivative at r={r}")
    upper_target = -sp.Rational(r + 2, r + 1) * X
    gate(sp.simplify(a_top * curl(beta_top) - upper_target) == 0,
         f"upper normalized debt at r={r}")
    gate(sp.expand(psi.subs(Z, 0) + 1) == 0,
         f"upper zero-fibre flank at r={r}")
    gate(sp.gcd(sp.Poly(psi, Z), sp.Poly(sp.diff(psi, Z), Z)).degree() == 0,
         f"upper cyclotomic flanks square-free at r={r}")
    gate(sp.expand(
        2 * (r + 1) * Z * psi
        + r * (1 - (1 - 2 * Z / r) ** (r + 1))
    ) == 0, f"upper roots-of-unity fibre identity at r={r}")
    if r >= 2:
        n = r - 2
        gate(A_top * (2 * n + 6) - (n + 2) == 0,
             f"upper finite-termination resonance at r={r}")

    bad_A_top = A_top + sp.Rational(1, 17)
    bad_a_top = sp.sqrt(bad_A_top)
    bad_N_top = sp.expand(M0 * sp.diag(bad_a_top, 1 / bad_a_top))
    bad_alpha_top = (bad_N_top[0, 0], bad_N_top[0, 1])
    bad_beta_top = (bad_N_top[1, 0], bad_N_top[1, 1])
    bad_upper = (
        sp.expand(bad_alpha_top[0] + h_top * bad_beta_top[0]),
        sp.expand(bad_alpha_top[1] + h_top * bad_beta_top[1]),
    )
    gate(sp.simplify(curl(bad_upper)) != 0,
         f"upper resonance hostile control at r={r}")

    # Charge +2 is the only component which can hit the target X.  Its
    # profile equation has the same unique-top-edge coefficient.
    for n in range(0, 9):
        coefficients = sp.symbols(f"up_{r}_{n}_0:{n + 1}")
        F = sum(coefficients[i] * Z**i for i in range(n + 1))
        f = sp.expand(X**2 * F.subs(Z, X * T))
        profile = sp.expand(Z * psi * sp.diff(F, Z) + 2 * K_top * F)
        gate(
            sp.expand(jac(f, q_top) - X * profile.subs(Z, X * T)) == 0,
            f"upper charge ODE at r={r}, degree={n}",
        )
        predicted = sp.expand(
            (n + 2 * (r + 1)) * upper_psi_lead * coefficients[n]
        )
        actual = sp.Poly(profile, Z).coeff_monomial(Z ** (n + r))
        gate(sp.expand(actual - predicted) == 0,
             f"upper unique top edge at r={r}, degree={n}")
        equation = sp.Poly(profile + sp.Rational(r + 2, r + 1), Z)
        system = [coefficient for _, coefficient in equation.terms()]
        gate(sp.linsolve(system, coefficients) == sp.EmptySet,
             f"upper finite hostile cokernel at r={r}, degree={n}")

    semantic_rows.append(
        f"r={r}:"
        + hashlib.sha256(
            "|".join(
                sp.srepr(item)
                for item in (A, K, H, phi, A_top, K_top, H_top, psi)
            ).encode()
        ).hexdigest()
    )

# The formal r -> infinity coefficient limits are the exponential kernels.
# Check the first eight coefficients directly from the closed binomial forms.
for degree in range(0, 8):
    r_symbol = sp.symbols("r_symbol", positive=True, integer=True)
    lower_coefficient = sp.binomial(r_symbol, degree) * (2 / r_symbol) ** degree
    upper_coefficient = -sp.binomial(r_symbol, degree) * (-2 / r_symbol) ** degree
    gate(sp.limit(lower_coefficient, r_symbol, sp.oo) == sp.Rational(2**degree, sp.factorial(degree)),
         f"lower exponential coefficient limit degree={degree}")
    gate(sp.limit(upper_coefficient, r_symbol, sp.oo) == -sp.Rational((-2) ** degree, sp.factorial(degree)),
         f"upper exponential coefficient limit degree={degree}")

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()
print("theorem=THM-3734-automorphic-Cohn-diagonal-binomial-towers")
print("scope=all_integer_depths_r>=1;paired_lower_and_upper_diagonal_SL2_slices")
print("first_stage=closed_unimodular_gradients_with_divided_power_potentials")
print("second_stage=charge_reduction_and_unique_top_edge_exclude_polynomial_mates")
print("formal_limit=exp(2XT)_and_-exp(-2XT)_kernels")
print(f"tested_depths={tested_depths.start}..{tested_depths.stop - 1}")
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
