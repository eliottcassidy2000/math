#!/usr/bin/env python3
"""Independent hostile audit of the coupled THM-3831/THM-3832 packet."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import sympy as sp


GATES = 0


def gate(condition: object, label: str) -> None:
    global GATES
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(label)
    GATES += 1


def lf(data: bytes) -> bytes:
    return data.replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def sha256_raw(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def zero(expression: sp.Expr, label: str) -> None:
    gate(sp.cancel(sp.factor(expression)) == 0, label)


def zero_mod(
    expression: sp.Expr, modulus: sp.Expr, variable: sp.Symbol, label: str
) -> None:
    numerator = sp.expand(sp.cancel(expression).as_numer_denom()[0])
    remainder = sp.rem(
        sp.Poly(numerator, variable), sp.Poly(modulus, variable)
    ).as_expr()
    gate(sp.expand(remainder) == 0, label)


def repository_root() -> Path:
    """Locate the repository from either canonical or session-scratch paths."""
    for parent in Path(__file__).resolve().parents:
        if (parent / "04-computation").is_dir() and (parent / "01-canon").is_dir():
            return parent
    raise RuntimeError("repository root not found")


ROOT = repository_root()
PRIMARY = [
    (
        ROOT / "04-computation" / "jc2_intrinsic_spectral_pencil_fibre_atlas_thm3831.py",
        ROOT / "05-knowledge" / "results" / "jc2_intrinsic_spectral_pencil_fibre_atlas_thm3831.out",
        "4dc7892f2e36d467c68e8268bbce48125e41fa968f913781a7ddc61c0dfd7ae9",
        "79217dcd7e35d491266ff3a0d2d22b6502e068dbf4b98d4f1bd37854961633d9",
        b"CHECKS=34\n",
    ),
    (
        ROOT / "04-computation" / "jc2_nonlinear_cubic_root_ratio_chart_thm3832.py",
        ROOT / "05-knowledge" / "results" / "jc2_nonlinear_cubic_root_ratio_chart_thm3832.out",
        "3e21279f5d1deac6960d6d3a7b17d9b8d7accfb9428120b518627192737c5460",
        "c3aafc9a27ae4daea7c82ec9811a8dcfdeee6bec56093857d225092be7846c52",
        b"CHECKS=25\n",
    ),
]

for script, output, script_hash, output_hash, count_line in PRIMARY:
    gate(sha256_raw(script) == script_hash, f"primary script hash {script.name}")
    gate(sha256_raw(output) == output_hash, f"primary output hash {output.name}")
    normal = lf(subprocess.check_output([sys.executable, str(script)], cwd=ROOT))
    optimized = lf(
        subprocess.check_output([sys.executable, "-O", str(script)], cwd=ROOT)
    )
    frozen = lf(output.read_bytes())
    gate(normal == optimized == frozen, f"normal optimized frozen {script.name}")
    gate(count_line in normal, f"active primary gate count {script.name}")

# ---------------------------------------------------------------------------
# Spectral fibre quotient, components, and signs.
# ---------------------------------------------------------------------------

alpha, kappa, t = sp.symbols("alpha kappa t", nonzero=True)
Z = sp.symbols("Z")
Q = 7 * alpha**2 + 3
P = 3 * alpha**3 + 7 * alpha**2 + 1
b_alpha = (alpha + 1) * (2 * alpha + 1) * (3 * alpha - 1)

gate(sp.discriminant(Q, alpha) != 0, "quadratic roots simple")
gate(sp.discriminant(P, alpha) != 0, "cubic roots simple")
gate(sp.gcd(sp.Poly(Q, alpha), sp.Poly(P, alpha)).degree() == 0,
     "quadratic and cubic roots disjoint")
gate(sp.gcd(sp.Poly(Q, alpha), sp.Poly(alpha * (9 * alpha + 14) * b_alpha, alpha)).degree() == 0,
     "quadratic denominators are units")
gate(sp.gcd(sp.Poly(P, alpha), sp.Poly(alpha * Q * b_alpha, alpha)).degree() == 0,
     "cubic denominators are units")
zero_mod(b_alpha + Q, P, alpha, "b=-Q on cubic slopes")

# Starting from zeta=C*k, the two lift laws have a single compatibility
# numerator.  This derivation is independent of either primary parametrization.
N0 = alpha**2 * (9 * alpha + 14) + Q
N1 = 2 * alpha**2 * (9 * alpha + 14) - Q
zero(N0 - 3 * P, "compatibility constant is 3P")
zero(N1 - 3 * b_alpha, "compatibility linear term is 3b")
compatibility = sp.expand(
    kappa * (N0 + N1 * Z) - 3 * alpha**2 * Q * Z**2
)
zero_mod(
    compatibility + 3 * Q * Z * (kappa + alpha**2 * Z),
    P,
    alpha,
    "cubic compatibility factors as Z(k+a^2Z)",
)
zero(
    (Z * (kappa + alpha**2 * Z)).subs(Z, sp.Symbol("Cc") * kappa)
    - sp.Symbol("Cc") * kappa**2 * (1 + alpha**2 * sp.Symbol("Cc")),
    "C-form cubic factorization",
)
Cc = sp.Symbol("Cc")
zero((1 + alpha**2 * Cc) - alpha**2 * Cc - 1,
     "cubic factors are comaximal by an explicit Bezout identity")


def packet(A0: sp.Expr, C0: sp.Expr, k0: sp.Expr) -> dict[str, sp.Expr]:
    h0 = alpha * k0
    D0 = sp.cancel(A0 / h0)
    omega0 = sp.cancel(k0 * D0)
    m0 = sp.cancel((C0 * k0 - 1) / h0)
    theta0 = sp.cancel((m0 - 14 * A0) / 3)
    return {
        "A": A0,
        "C": C0,
        "k": k0,
        "h": h0,
        "D": D0,
        "omega": omega0,
        "m": m0,
        "theta": theta0,
    }


def relations(values: dict[str, sp.Expr]) -> list[sp.Expr]:
    A0, C0 = values["A"], values["C"]
    k0, h0, D0 = values["k"], values["h"], values["D"]
    omega0, m0, theta0 = values["omega"], values["m"], values["theta"]
    return [
        C0 * k0 - m0 * h0 - 1,
        D0 * (7 * h0**2 + 3 * k0**2) - 1 - 2 * C0 * k0,
        h0 * D0 * (9 * h0 + 14 * k0) - k0 * m0 - 3 * h0 * C0**2,
        A0 - h0 * D0,
        omega0 - k0 * D0,
        m0 - 3 * theta0 - 14 * A0,
        omega0**2 + 7 * A0**2 - C0 * omega0 + A0 * theta0,
        omega0 * theta0 - 3 * A0**2 + A0 * C0**2,
        theta0**2 - 3 * A0 * C0 + C0**3
        - (C0**2 - 3 * A0) * omega0 + 7 * A0 * theta0,
        D0 - C0 * omega0 + 3 * A0 * theta0 + 14 * A0**2,
    ]


def normalized_sign(values: dict[str, sp.Expr]) -> sp.Expr:
    h0, k0, D0 = values["h"], values["k"], values["D"]
    q0 = 7 * h0**2 + 3 * k0**2
    B0 = (
        k0**5 - 7 * h0**2 * k0**3 - 3 * h0**2 * k0**2
        - 6 * h0**3 * k0**2 - 7 * h0**4
    )
    W0 = sp.cancel((h0**2 * q0**2 * D0 + B0) / k0)
    B30 = (h0 + k0) * (2 * h0 + k0) * (3 * h0 - k0)
    return sp.cancel(W0 / (k0 * B30))


D_quadratic = (14 * t + 3) / (4 * (9 * alpha + 14) * t**3)
quadratic = packet(alpha * t * D_quadratic, -1 / (2 * t), t)
D_minus = 1 / (Q * t**2)
cubic_minus = packet(alpha * t * D_minus, sp.Integer(0), t)
D_plus = (1 - 2 * t / alpha**2) / (Q * t**2)
cubic_plus = packet(alpha * t * D_plus, -1 / alpha**2, t)

for name, values, modulus in [
    ("quadratic", quadratic, Q),
    ("cubic_minus", cubic_minus, P),
    ("cubic_plus", cubic_plus, P),
]:
    for index, relation in enumerate(relations(values)):
        zero_mod(relation, modulus, alpha, f"{name} relation {index}")

zero_mod(quadratic["D"] - D_quadratic, Q, alpha,
         "quadratic D formula")
zero_mod(cubic_minus["D"] - D_minus, P, alpha,
         "minus cubic D formula")
zero_mod(cubic_plus["D"] - D_plus, P, alpha,
         "plus cubic D formula")
zero_mod(normalized_sign(quadratic) + 1, Q, alpha,
         "quadratic sign is minus")
zero_mod(normalized_sign(cubic_minus) + 1, P, alpha,
         "C=0 cubic sign is minus")
zero_mod(normalized_sign(cubic_plus) - 1, P, alpha,
         "C=-1/a^2 cubic sign is plus")

# The candidate Laurent components contain real D=0 boundary points, but D is
# not identically zero on any component.  Thus D!=0 is dense on each reduced
# Laurent component and every generic inverse identity extends to those points.
gate(sp.factor(sp.cancel(quadratic["D"]).as_numer_denom()[0]) != 0,
     "quadratic D not identically zero")
gate(sp.factor(sp.cancel(cubic_minus["D"]).as_numer_denom()[0]) != 0,
     "minus cubic D not identically zero")
gate(sp.factor(sp.cancel(cubic_plus["D"]).as_numer_denom()[0]) != 0,
     "plus cubic D not identically zero")
zero_mod(quadratic["D"].subs(t, -sp.Rational(3, 14)), Q, alpha,
         "quadratic Laurent component includes D=0 point")
zero_mod(cubic_plus["D"].subs(t, alpha**2 / 2), P, alpha,
         "plus cubic component includes D=0 point")
gate(-sp.Rational(3, 14) != 0, "quadratic D=0 point retains Laurent parameter")
gate(sp.gcd(sp.Poly(P, alpha), sp.Poly(alpha, alpha)).degree() == 0,
     "plus D=0 parameter alpha^2/2 is nonzero")


# ---------------------------------------------------------------------------
# Birational root-ratio chart, density, and spectral-boundary coupling.
# ---------------------------------------------------------------------------

z, C = sp.symbols("z C")
r = 3 * z**3 + 7 * z**2 + 1
q = 7 * z**2 + 3
b = 6 * z**3 + 7 * z**2 - 1
s = z**2 * q * C - b
A_chart = z * C * (1 + z**2 * C) / r
omega_chart = C * (1 + z**2 * C) / r
theta_chart = -z * C * (7 * C * z**2 + C - 3 * z) / r
D_chart = C**2 * (1 + z**2 * C) * s / r**2
k_chart = r / (C * s)
h_chart = z * k_chart
m_chart = -z * C * (7 * C * z**2 + 3 * C - 9 * z - 14) / r

A_symbol = sp.Symbol("A_symbol")
characteristic_after_ratio = sp.expand(
    z**3 * (
        (A_symbol / z) ** 3
        - C * (A_symbol / z) ** 2
        + 7 * A_symbol**2 * (A_symbol / z)
        + 3 * A_symbol**3
        - A_symbol**2 * C**2
    )
)
zero(
    characteristic_after_ratio
    - A_symbol**2 * (A_symbol * r - z * C * (1 + z**2 * C)),
    "marked-root cubic gives triangular equation",
)

chart_values = {
    "A": A_chart,
    "C": C,
    "k": k_chart,
    "h": h_chart,
    "D": D_chart,
    "omega": omega_chart,
    "m": m_chart,
    "theta": theta_chart,
}
for index, relation in enumerate(relations(chart_values)):
    zero(relation, f"chart intrinsic relation {index}")

inverse_k = sp.cancel(A_chart * q / z - 2 * C)
zero(inverse_k - C * s / r, "birational inverse 1/k=C*s/r")
zero(inverse_k - 1 / k_chart, "birational inverse composes with k")
zero(h_chart / k_chart - z, "birational inverse composes with z")
zero(A_chart * k_chart - omega_chart * h_chart,
     "A/omega=h/k only by function-field cross multiplication")

target_density = sp.factor(sp.diff(A_chart, z))
zero(target_density - C * s / r**2,
     "triangular Jacobian density factorization")
zero(target_density - inverse_k / r,
     "density is 1/(k*r)")
zx, zy, Cx, Cy, lam = sp.symbols("z_x z_y C_x C_y lambda")
source_jacobian = zx * Cy - zy * Cx
zero(target_density * source_jacobian - source_jacobian / (k_chart * r),
     "source chain rule")
zero(
    (target_density * source_jacobian - lam).subs(
        source_jacobian, lam * k_chart * r
    ),
    "weighted area equation is equivalent to Keller",
)

zero(b + q - 2 * r, "b+q=2r")
gate(sp.discriminant(r, z) != 0, "cubic chart poles simple")
gate(sp.gcd(sp.Poly(r, z), sp.Poly(z * q * b, z)).degree() == 0,
     "cubic poles avoid every displayed unit denominator")

# At r(alpha)=0 the two numerator cancellations are exactly the two intrinsic
# components.  They also resolve two different denominator factors of k:
# C=0 on the minus arm, s=0 on the plus arm.
zero_mod((z * C * (1 + z**2 * C)).subs({z: alpha, C: 0}), P, alpha,
         "minus address cancels A numerator")
zero_mod(
    (z * C * (1 + z**2 * C)).subs({z: alpha, C: -1 / alpha**2}),
    P,
    alpha,
    "plus address cancels A numerator",
)
zero_mod(s.subs({z: alpha, C: 0}) - Q, P, alpha,
         "minus address has C=0 but s=Q unit")
zero_mod(s.subs({z: alpha, C: -1 / alpha**2}), P, alpha,
         "plus address has s=0 but C unit")

# The boundary polynomial identities compose with all three Laurent gauges.
triangular = A_symbol * r - z * C * (1 + z**2 * C)
invk_general = A_symbol * q / z - 2 * C
zero_mod(
    triangular.subs({z: alpha, C: quadratic["C"], A_symbol: quadratic["A"]}),
    Q,
    alpha,
    "quadratic Laurent gauge lies on triangular chart equation",
)
zero_mod(
    invk_general.subs({z: alpha, C: quadratic["C"], A_symbol: quadratic["A"]})
    - 1 / quadratic["k"],
    Q,
    alpha,
    "quadratic Laurent gauge matches inverse k",
)
for name, values in [("minus", cubic_minus), ("plus", cubic_plus)]:
    zero_mod(
        triangular.subs({z: alpha, C: values["C"], A_symbol: values["A"]}),
        P,
        alpha,
        f"cubic {name} gauge lies on triangular boundary equation",
    )
    zero_mod(
        invk_general.subs({z: alpha, C: values["C"], A_symbol: values["A"]})
        - 1 / values["k"],
        P,
        alpha,
        f"cubic {name} gauge matches inverse k",
    )

# On q=0 there is no two-address pole: s=-2r and k=-1/(2C).
zero_mod(s + 2 * r, q, z, "quadratic slope has s=-2r")
zero_mod(k_chart + 1 / (2 * C), q, z,
         "quadratic slope has the single Laurent C chart")

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "independent checker has no Python assert",
)

semantic = {
    "verdict": "PASS for promoted THM-3831 and THM-3832",
    "quotients": "Q roots: one Laurent Gm; P roots: C(1+a^2C)=0, two comaximal Laurent Gm components",
    "signs": "quadratic minus; cubic C=0 minus; cubic C=-1/a^2 plus",
    "boundary_extension": "each reduced Laurent component has dense D!=0; generic inverse identities extend across its isolated D=0 points",
    "birational": "A*r=z*C*(1+z^2C); 1/k=A*q/z-2C=C*s/r; all intrinsic generators reconstructed",
    "density": "partial_z A=C*s/r^2=1/(k*r); Jac(z,C)=lambda*k*r",
    "coupling": "at cubic r=0, minus resolves C=0 while plus resolves s=0; both match the intrinsic components",
    "scope": "rational chart and necessary irreducible-h incidence passport only; polynomial atlas, reducible-h branch, and JC2 remain open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("audit=THM-3831-THM-3832-coupled-independent-hostile")
print("verdict=PASS")
print("primary_replay=THM3831:34;THM3832:25;normal==optimized==frozen")
print("quotients=quadratic:Gm;cubic:Gm_minus_disjoint_union_Gm_plus")
print("signs=quadratic:minus;cubic_C0:minus;cubic_C_minus_a^-2:plus")
print("boundary=D_nonzero_dense_on_each_reduced_component;D0_points_included")
print("chart=A*r=z*C*(1+z^2*C);inverse_k=A*q/z-2*C=C*s/r")
print("density=partial_z_A=C*s/r^2=1/(k*r)")
print("spectral_coupling=minus_resolves_C0;plus_resolves_s0")
print("scope=polynomial_atlas_reducible_h_and_JC2_OPEN")
print(f"GATES={GATES}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
