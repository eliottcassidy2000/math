#!/usr/bin/env python3
"""Exact companion for THM-3257's tuned cubic-root infinity wall.

The inherited THM-3237 family is enlarged from

    B_t = 1 + t*x^9

to

    B_(s,t) = 1 + s*x^8 + t*x^9.

The degree-eight term cannot change the generic infinity face by itself.
On the degree-nine wall, however, it tunes the next strict-transform
coefficient.  At one exact value it cancels that coefficient and exposes a
nonzero cubic carry.  No assertion statement, floating literal, or random
choice is used.
"""

from __future__ import annotations

import ast
import importlib.util
import sys
from fractions import Fraction
from hashlib import sha256
from pathlib import Path

import sympy as sp


Q = Fraction
BASE_PATH = Path(__file__).with_name(
    "jc_heptic_degree_nine_infinity_wall_thm3237.py"
)
BASE_SHA256 = "92518f258afeca233e90790fa2f713fcfd375c295271ff85dc4c5c66c0057d81"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_hash(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load_base():
    require(lf_hash(BASE_PATH) == BASE_SHA256, "THM-3237 companion pin")
    spec = importlib.util.spec_from_file_location("thm3237_dependency", BASE_PATH)
    require(spec is not None and spec.loader is not None, "dependency import spec")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


base = load_base()


def rational_text(value: Q) -> str:
    return f"{value.numerator}/{value.denominator}"


def universal_top_jet_checks() -> None:
    x = sp.symbols("x")
    v, v15, v14, v13 = sp.symbols("v v15 v14 v13", nonzero=True)
    a, a7, a6, a5 = sp.symbols("a a7 a6 a5", nonzero=True)
    s, t = sp.symbols("s t")

    V = v * x**16 + v15 * x**15 + v14 * x**14 + v13 * x**13
    A = a * x**8 + a7 * x**7 + a6 * x**6 + a5 * x**5
    B = 1 + s * x**8 + t * x**9
    J = 2 * V * sp.diff(B, x) - B * sp.diff(V, x)
    K = (
        -A**3 * J**3
        + 12 * A**2 * J**2 * V**2
        - 4 * A * B**3 * J**2 * V
        + 4 * A * B**2 * J * V**2
        + 24 * A * B * J * V**3
        - 48 * A * J * V**4
        - 16 * A * V**4
        - 8 * B**4 * J * V**2
        + 8 * B**3 * J * V**3
        + 32 * B**2 * V**4
        - 96 * B * V**5
        + 64 * V**6
    )
    polynomial = sp.Poly(sp.expand(K), x)
    c99 = polynomial.coeff_monomial(x**99)
    c98 = polynomial.coeff_monomial(x**98)
    c97 = polynomial.coeff_monomial(x**97)
    c96 = polynomial.coeff_monomial(x**96)

    require(
        sp.factor(c99 + 16 * t**4 * v**3 * (a * t - v)) == 0,
        "universal x^99 face",
    )

    a7_value = a * v15 / (2 * v)
    a6_value = a * (4 * v * v14 - v15**2) / (8 * v**2)
    a5_value = (
        a * (8 * v**2 * v13 - 4 * v * v14 * v15 + v15**3)
        / (16 * v**3)
    )
    wall_substitution = {
        t: v / a,
        a7: a7_value,
        a6: a6_value,
        a5: a5_value,
    }
    require(
        sp.factor(c98.subs(wall_substitution)) == 0,
        "degree-eight term preserves second wall cancellation",
    )

    kappa = 4 * v * v14 - 3 * v15**2
    expected_c97 = -2 * v**6 * (kappa + 4 * a * s * v15) / a**4
    require(
        sp.factor(c97.subs(wall_substitution) - expected_c97) == 0,
        "linear tuned x^97 carry",
    )

    s_star = -kappa / (4 * a * v15)
    tau = (
        -8 * a**4 * v * v15
        + 16 * v**2 * v13 * v15
        - 16 * v**2 * v14**2
        + v15**4
    )
    expected_c96 = -v**5 * tau / (a**4 * v15)
    require(
        sp.factor(
            c96.subs(wall_substitution).subs({s: s_star}) - expected_c96
        )
        == 0,
        "cubic x^96 strict-transform carry",
    )

    response = sp.Poly(sp.expand(2 * V * sp.diff(A, x) - A * sp.diff(V, x)), x)
    require(
        sp.factor(response.coeff_monomial(x**22) - (a * v15 - 2 * a7 * v))
        == 0,
        "first response jet",
    )
    print(
        "universal_top_jet="
        "c99=-16*t^4*v^3*(a*t-v);"
        "wall_c98=0;"
        "wall_c97=-2*v^6*(kappa+4*a*s*v15)/a^4;"
        "tuned_c96=-v^5*tau/(a^4*v15)"
    )


EXPECTED = {
    "4111": {
        "r1_norm": Q(-1024, 175),
        "K_norm": Q(23328, 30625),
        "R_norm": Q(58451308942336, 937890625),
        "s_norm": Q(-250047, 32768000),
    },
    "3211": {
        "r1_norm": Q(-69632, 525),
        "K_norm": Q(-2147586048, 30625),
        "R_norm": Q(122694886932611072, 8441015625),
        "s_norm": Q(-173456171875, 35651584),
    },
}


def accessory_case(name: str):
    if name == "4111":
        base.configure(100, 244, 237, 44)
        u = base.U_A
        exponent_a, exponent_b = 4, 1
        accessory_v = (8 * u**2 + 9 * u + 8) / 7
        shift = 5 * (u + 1) / 7
        A0 = 80 * accessory_v**2 * (u + 1) / 343
        extras = (9, 0)
    else:
        base.configure(75, -89, -31, 61)
        u = base.U_A
        exponent_a, exponent_b = 3, 2
        accessory_v = (24 * u**2 - 16 * u - 16) / 21
        shift = (5 * u - 4) / 7
        A0 = 9 * accessory_v**2 * (5 * u - 4) / 343
        extras = (6, 3)

    gamma = -7 * A0
    q2 = base.X_P**2 - base.X_P.scalar(u) + accessory_v
    D = base.X_P**exponent_a * (base.X_P - 1) ** exponent_b * q2
    T = base.X_P * (base.X_P - 1) * q2
    E = (
        (base.X_P - 1) * q2 * exponent_a
        + base.X_P * q2 * exponent_b
        + base.X_P * (base.X_P - 1) * (2 * base.X_P - u)
    ).scalar(Q(1, 7))
    S = base.X_P + shift
    V = (S * D * T**2).scalar(4 / gamma**2)
    A = (S * E * T).scalar(2 / gamma)
    boundary = (
        S**3
        * T**8
        * base.X_P**extras[0]
        * (base.X_P - 1) ** extras[1]
    )
    require(V.degree == 16 and A.degree == 8, f"{name} response degrees")
    require(boundary.degree == 44, f"{name} boundary degree")
    require(2 * V * A.derivative() - A * V.derivative() == 2 * V, f"{name} response identity")

    v = V.c[16]
    a = A.c[8]
    require(v == a**2 and a == 2 / gamma, f"{name} normalized leading jets")
    r1 = V.c[15] / v
    r2 = V.c[14] / v
    r3 = V.c[13] / v
    K = 4 * r2 - 3 * r1**2
    R = r1**4 + 16 * r1 * r3 - 16 * r2**2 - 8 * r1
    s_star = -a * K / (4 * r1)
    expected = EXPECTED[name]
    require(base.norm(r1) == expected["r1_norm"], f"{name} r1 norm")
    require(base.norm(K) == expected["K_norm"], f"{name} K norm")
    require(base.norm(R) == expected["R_norm"], f"{name} R norm")
    require(base.norm(s_star) == expected["s_norm"], f"{name} s-star norm")
    require(
        expected["r1_norm"] != 0 and expected["K_norm"] != 0 and expected["R_norm"] != 0,
        f"{name} nonzero normalized carries",
    )

    B = (
        base.ONE_P
        + (base.X_P**8).scalar(s_star)
        + (base.X_P**9).scalar(a)
    )
    J = V * B.derivative().scalar(2) - B * V.derivative()
    H = base.K_general(V, A, B, J).exact_quotient(boundary)
    require(H.degree == 52, f"{name} tuned degree")
    h52 = -a**12 * R / r1
    require(H.c[52] == h52, f"{name} normalized cubic carry")

    print(
        f"accessory_case={name} "
        f"r1_norm={rational_text(expected['r1_norm'])} "
        f"K_norm={rational_text(expected['K_norm'])} "
        f"R_norm={rational_text(expected['R_norm'])} "
        f"s_star_norm={rational_text(expected['s_norm'])} "
        f"degree={H.degree} "
        f"H_sha256={sha256(base.serialize(H)).hexdigest()}"
    )
    return gamma, s_star, H


def finite_case(name: str) -> None:
    x = sp.symbols("x")
    if name == "4111":
        p, u, exponent_a, exponent_b = 113, 85, 4, 1
        q_coefficients = (100, 244, 237, 44)
        accessory_v = (8 * u**2 + 9 * u + 8) * pow(7, -1, p) % p
        shift = 5 * (u + 1) * pow(7, -1, p) % p
        A0 = 80 * accessory_v**2 * (u + 1) * pow(343, -1, p) % p
        extras = (9, 0)
        expected = (85, 65, 64)
    else:
        p, u, exponent_a, exponent_b = 101, 64, 3, 2
        q_coefficients = (75, -89, -31, 61)
        accessory_v = (24 * u**2 - 16 * u - 16) * pow(21, -1, p) % p
        shift = (5 * u - 4) * pow(7, -1, p) % p
        A0 = 9 * accessory_v**2 * (5 * u - 4) * pow(343, -1, p) % p
        extras = (6, 3)
        expected = (89, 73, 91)

    c3, c2, c1, c0 = q_coefficients
    require((c3 * u**3 + c2 * u**2 + c1 * u + c0) % p == 0, f"{name} accessory root")
    require((3 * c3 * u**2 + 2 * c2 * u + c1) % p != 0, f"{name} simple accessory root")
    gamma = -7 * A0 % p
    require(gamma != 0, f"{name} gamma unit")
    P = lambda expression: sp.Poly(expression, x, modulus=p)
    q2 = P(x**2 - u * x + accessory_v)
    D = P(x**exponent_a * (x - 1) ** exponent_b) * q2
    T = P(x * (x - 1)) * q2
    E = P(
        (
            exponent_a * (x - 1) * (x**2 - u * x + accessory_v)
            + exponent_b * x * (x**2 - u * x + accessory_v)
            + x * (x - 1) * (2 * x - u)
        )
        * pow(7, -1, p)
    )
    S = P(x + shift)
    V = (S * D * T**2).mul_ground(4 * pow(gamma**2, -1, p) % p)
    A = (S * E * T).mul_ground(2 * pow(gamma, -1, p) % p)
    boundary = S**3 * T**8 * P(x) ** extras[0] * P(x - 1) ** extras[1]
    owner_boundary = S * T

    v = int(V.nth(16)) % p
    r1 = int(V.nth(15)) * pow(v, -1, p) % p
    r2 = int(V.nth(14)) * pow(v, -1, p) % p
    K = (4 * r2 - 3 * r1**2) % p
    a = 2 * pow(gamma, -1, p) % p
    s_star = -a * K * pow(4 * r1, -1, p) % p
    t_star = a
    B = P(1 + s_star * x**8 + t_star * x**9)
    J = V.mul_ground(2) * B.diff() - B * V.diff()
    H = sp.exquo(base.K_general(V, A, B, J), boundary)
    controls = (
        H.degree(),
        sp.gcd(B, owner_boundary).degree(),
        sp.gcd(H, owner_boundary).degree(),
        sp.gcd(H, H.diff()).degree(),
        int(H.nth(52)) % p,
    )
    require((t_star, s_star, controls[-1]) == expected, f"{name} tuned finite values")
    require(controls[:4] == (52, 0, 0, 0), f"{name} tuned finite gcd gates")
    print(
        f"finite_case={name} good_reduction=(p={p},u={u}) "
        f"(t_star,s_star)=({t_star},{s_star}) "
        f"controls=(degree,gcd_B_boundary,gcd_H_boundary,gcd_H_Hprime,h52)={controls}"
    )


def source_audit() -> None:
    tree = ast.parse(Path(__file__).read_text())
    assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    float_nodes = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(tree)
    )
    require(assert_nodes == 0 and float_nodes == 0, "source AST gate")
    print(f"source_ast=(assert_nodes={assert_nodes},float_literals={float_nodes})")


def main() -> None:
    print(f"thm3237_dependency_sha256={lf_hash(BASE_PATH)}")
    universal_top_jet_checks()
    for name in ("4111", "3211"):
        accessory_case(name)
    for name in ("4111", "3211"):
        finite_case(name)
    print("cubic_law=w^3~-16*r1*delta/(a*R)=-8*r1*epsilon/R")
    print("local_inertia=one_3_cycle_on_the_three_escaping_critical_roots")
    source_audit()


if __name__ == "__main__":
    main()
