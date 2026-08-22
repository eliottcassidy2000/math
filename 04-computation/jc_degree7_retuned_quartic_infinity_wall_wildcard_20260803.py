#!/usr/bin/env python3
"""Exact wildcard probe after THM-3257's tuned cubic infinity wall.

The three-parameter clutch is

    B_(r,s,t) = 1 + r*x^7 + s*x^8 + t*x^9.

The direct hope that ``r`` first changes the cubic carry is false: with the
old THM-3257 value of ``s``, it reopens the preceding x^97 layer.  Retuning
``s`` along that layer produces a one-dimensional strict transform.  On it,
``r`` changes x^96 linearly, has one exact cancelling value, and exposes a
nonzero x^95 quartic carry in both accessory fields.

This is a critical-resultant computation.  It supplies no polynomial inverse
cover, branchwise Keller cofactor, or conclusion about JC(2).
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
ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "04-computation/jc_heptic_degree_nine_infinity_wall_thm3237.py"


DEPENDENCIES = {
    "04-computation/jc_heptic_degree_nine_infinity_wall_thm3237.py":
        "92518f258afeca233e90790fa2f713fcfd375c295271ff85dc4c5c66c0057d81",
    "05-knowledge/results/jc_heptic_degree_nine_infinity_wall_thm3237.out":
        "38599c5d2a9b527098274d0dfccd427b5f091528671df168ca3a7ffb31c3b9cb",
    "04-computation/jc_degree8_tuned_cubic_infinity_wall_thm3257.py":
        "f7c9b6d92204ab0af0271311bfbcfa11fb51014a1e3d875d271ff406898da06d",
    "05-knowledge/results/jc_degree8_tuned_cubic_infinity_wall_thm3257.out":
        "1ade6b505bf92f7cb1a17395a3fcff22ee9ac9692ea1ac449b2bee812831ab1b",
    "04-computation/quartic_tame_inertia_clutch_index_resonance_thm3057.py":
        "6fa2cfc98b36270d4615fad4645cb6985b0bd55d57135fa71953f8ee20b51573",
    "05-knowledge/results/quartic_tame_inertia_clutch_index_resonance_thm3057.out":
        "fb66c28e57c57939bdebf84e9b639a6da83e430ec48ce45b555002f5f6169d9d",
    "04-computation/quartic_twojet_even_jelonek_c3_escape_thm3059.py":
        "16a1d4813950bc623a48f99e0d140f5aff69b465ff872378cd14815084abb694",
    "05-knowledge/results/quartic_twojet_even_jelonek_c3_escape_thm3059.out":
        "af579b31fd427579cfcf6797e9d10a7c794646d8452e6a65d60ba7809fd50b26",
    "04-computation/pointed_cubic_norm_keller_decoder_thm3064.py":
        "fb82a95d8ccd17e4d509157353590e895e55a0fb90c14ba8d9f3b5fb82dcd696",
    "05-knowledge/results/pointed_cubic_norm_keller_decoder_thm3064.out":
        "9b51211a8c6f24a29f79d74db4311242719e15b1ad84789411eb80fcf6686458",
    "04-computation/quartic_c3_hafnian_cofactor_blindness_thm3066.py":
        "06a9d2d4329a82a05127c789139180141c6ea345e53b77954798b26b957b7211",
    "05-knowledge/results/quartic_c3_hafnian_cofactor_blindness_thm3066.out":
        "f4c36f36e710846025133559c4b7cf470e00829d577ef7010e1b61ded1bbb9ed",
    "04-computation/c3_escape_reciprocal_cofactor_shift_thm3068.py":
        "97072cec79f764752364557fb69fd4bbd7dd9b34e51592f3c66908736682dec5",
    "05-knowledge/results/c3_escape_reciprocal_cofactor_shift_thm3068.out":
        "884dec30edb43d68149038bf857280d35eea626e4913896d893c2bdb951bfa26",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_hash(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


for dependency, expected_hash in DEPENDENCIES.items():
    require(lf_hash(ROOT / dependency) == expected_hash, f"dependency drift: {dependency}")


def load_base():
    spec = importlib.util.spec_from_file_location("thm3237_dependency", BASE_PATH)
    require(spec is not None and spec.loader is not None, "dependency import spec")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


base = load_base()


def rational_text(value: Q) -> str:
    return f"{value.numerator}/{value.denominator}"


def normalized_invariants(r1, r2, r3, r4):
    kappa = 4 * r2 - 3 * r1**2
    lam = 4 * r2 - r1**2
    rho = r1**4 + 16 * r1 * r3 - 16 * r2**2 - 8 * r1
    sigma = (
        128
        - 128 * r3
        + 256 * r2 * r4
        - 256 * r3**2
        - 64 * r1**2 * r4
        + 128 * r1 * r2 * r3
        - 64 * r2**3
        - 32 * r1**3 * r3
        + 48 * r1**2 * r2**2
        - 12 * r1**4 * r2
        + r1**6
    )
    return kappa, lam, rho, sigma


def universal_top_jet_checks() -> None:
    w = sp.symbols("w")
    a, r1, r2, r3, r4 = sp.symbols("a r1 r2 r3 r4", nonzero=True)
    r, s, t = sp.symbols("r s t")

    # Reciprocal jets after using v=a^2 and the top four homogeneous
    # consequences of 2*V*A'-A*V'=2*V.
    vbar = a**2 * (1 + r1 * w + r2 * w**2 + r3 * w**3 + r4 * w**4)
    abar = a * (
        1
        + r1 * w / 2
        + (4 * r2 - r1**2) * w**2 / 8
        + (8 * r3 - 4 * r1 * r2 + r1**3) * w**3 / 16
        + (
            64 * r4
            - 32 * r1 * r3
            - 16 * r2**2
            + 24 * r1**2 * r2
            - 5 * r1**4
        )
        * w**4
        / 128
    )
    bbar = t + s * w + r * w**2
    jbar = (
        2 * vbar * bbar
        - 2 * w * vbar * sp.diff(bbar, w)
        + w * bbar * sp.diff(vbar, w)
    )

    # Only the degree-99 and degree-96 groups can reach x^95.  Thus this
    # truncated reciprocal computation is an exact top-jet computation, not
    # a numerical series approximation.
    degree99 = -4 * abar * bbar**3 * jbar**2 * vbar + 8 * bbar**3 * jbar * vbar**3
    degree96 = (
        -abar**3 * jbar**3
        + 12 * abar**2 * jbar**2 * vbar**2
        - 48 * abar * jbar * vbar**4
        + 64 * vbar**6
    )
    reciprocal = sp.series(degree99 + w**3 * degree96, w, 0, 5).removeO().expand()
    coefficients = tuple(sp.factor(reciprocal.coeff(w, index)) for index in range(5))

    kappa, lam, rho, sigma = normalized_invariants(r1, r2, r3, r4)
    wall = {t: a}
    c99, c98, c97, c96, c95 = tuple(sp.factor(value.subs(wall)) for value in coefficients)
    expected_c97 = -2 * a**11 * (a * kappa + 4 * r1 * s - 8 * r)
    require(c99 == 0 and c98 == 0, "primitive wall top cancellations")
    require(sp.factor(c97 - expected_c97) == 0, "degree-seven reopens x^97")

    old_s = -a * kappa / (4 * r1)
    require(sp.factor(c97.subs(s, old_s) - 16 * a**11 * r) == 0,
            "old cubic wall transverse obstruction")

    s_on_wall = (8 * r - a * kappa) / (4 * r1)
    expected_c96 = -a**12 * (rho + 8 * r * lam / a) / r1
    require(sp.factor(c96.subs(s, s_on_wall) - expected_c96) == 0,
            "retuned linear x^96 carry")

    r_star = -a * rho / (8 * lam)
    s_star = sp.factor(s_on_wall.subs(r, r_star))
    expected_c95 = -3 * a**12 * sigma / (8 * lam)
    require(sp.factor(c95.subs({r: r_star, s: s_star}) - expected_c95) == 0,
            "quartic x^95 carry")
    require(sp.factor(sp.diff(coefficients[0], t).subs(t, a) + 16 * a**11) == 0,
            "transverse leading derivative")

    print(
        "universal_top_jet="
        "wall_c99=wall_c98=0;"
        "wall_c97=-2*a^11*(a*Kappa+4*r1*s-8*r);"
        "old_s_c97=16*a^11*r;"
        "retuned_c96=-a^12*(Rho+8*r*Lambda/a)/r1;"
        "quartic_c95=-3*a^12*Sigma/(8*Lambda)"
    )


EXPECTED = {
    "4111": {
        "lambda_norm": Q(12448544, 30625),
        "sigma_norm": Q(-157881139404422797328384, 28722900390625),
        "r_star_norm": Q(2796986463061, 2489708800000),
        "s_star_norm": Q(-15118672709, 6224272000),
        "H_sha256": "3c2ac8e2cbd13ee2568f2565bd780e2b0b8bb11facad583cc4b626e1f2265626",
    },
    "3211": {
        "lambda_norm": Q(-42471424, 275625),
        "sigma_norm": Q(64559323245631953174528, 28722900390625),
        "r_star_norm": Q(8958651810434375, 82564448256),
        "s_star_norm": Q(-1836115611739609375, 14090999169024),
        "H_sha256": "2808b49be3b3ac590efdd39e1015da59e620d45b11bd103271e21e82bdf6d499",
    },
}


def accessory_case(name: str) -> None:
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
    boundary = S**3 * T**8 * base.X_P**extras[0] * (base.X_P - 1) ** extras[1]
    require(V.degree == 16 and A.degree == 8 and boundary.degree == 44,
            f"{name} inherited degrees")
    require(2 * V * A.derivative() - A * V.derivative() == 2 * V,
            f"{name} response identity")

    v = V.c[16]
    a = A.c[8]
    require(v == a**2 and a == 2 / gamma, f"{name} normalized leading jets")
    r1, r2, r3, r4 = tuple(V.c[16 - index] / v for index in range(1, 5))
    _, lam, rho, sigma = normalized_invariants(r1, r2, r3, r4)
    require(r1 != 0 and lam != 0 and rho != 0 and sigma != 0,
            f"{name} nonzero strict-transform gates")
    r_star = -a * rho / (8 * lam)
    kappa = 4 * r2 - 3 * r1**2
    s_star = (8 * r_star - a * kappa) / (4 * r1)
    expected = EXPECTED[name]
    require(base.norm(lam) == expected["lambda_norm"], f"{name} Lambda norm")
    require(base.norm(sigma) == expected["sigma_norm"], f"{name} Sigma norm")
    require(base.norm(r_star) == expected["r_star_norm"], f"{name} r-star norm")
    require(base.norm(s_star) == expected["s_star_norm"], f"{name} s-star norm")

    B = (
        base.ONE_P
        + (base.X_P**7).scalar(r_star)
        + (base.X_P**8).scalar(s_star)
        + (base.X_P**9).scalar(a)
    )
    J = V * B.derivative().scalar(2) - B * V.derivative()
    H = base.K_general(V, A, B, J).exact_quotient(boundary)
    h51 = -3 * a**12 * sigma / (8 * lam)
    require(H.degree == 51 and H.c[51] == h51, f"{name} quartic carry")
    H_hash = sha256(base.serialize(H)).hexdigest()
    require(H_hash == expected["H_sha256"], f"{name} H digest")

    print(
        f"accessory_case={name} "
        f"Lambda_norm={rational_text(expected['lambda_norm'])} "
        f"Sigma_norm={rational_text(expected['sigma_norm'])} "
        f"r_star_norm={rational_text(expected['r_star_norm'])} "
        f"s_star_norm={rational_text(expected['s_star_norm'])} "
        f"degree={H.degree} H_sha256={H_hash}"
    )


def finite_case(name: str) -> None:
    x = sp.symbols("x")
    if name == "4111":
        prime, u, exponent_a, exponent_b = 113, 85, 4, 1
        q_coefficients = (100, 244, 237, 44)
        accessory_v = (8 * u**2 + 9 * u + 8) * pow(7, -1, prime) % prime
        shift = 5 * (u + 1) * pow(7, -1, prime) % prime
        A0 = 80 * accessory_v**2 * (u + 1) * pow(343, -1, prime) % prime
        extras = (9, 0)
        expected = (85, 103, 91, 41)
    else:
        prime, u, exponent_a, exponent_b = 101, 64, 3, 2
        q_coefficients = (75, -89, -31, 61)
        accessory_v = (24 * u**2 - 16 * u - 16) * pow(21, -1, prime) % prime
        shift = (5 * u - 4) * pow(7, -1, prime) % prime
        A0 = 9 * accessory_v**2 * (5 * u - 4) * pow(343, -1, prime) % prime
        extras = (6, 3)
        expected = (89, 87, 55, 9)

    c3, c2, c1, c0 = q_coefficients
    require((c3 * u**3 + c2 * u**2 + c1 * u + c0) % prime == 0,
            f"{name} accessory root")
    require((3 * c3 * u**2 + 2 * c2 * u + c1) % prime != 0,
            f"{name} simple accessory root")
    gamma = -7 * A0 % prime
    require(gamma != 0, f"{name} gamma unit")
    P = lambda expression: sp.Poly(expression, x, modulus=prime)
    q2 = P(x**2 - u * x + accessory_v)
    D = P(x**exponent_a * (x - 1) ** exponent_b) * q2
    T = P(x * (x - 1)) * q2
    E = P(
        (
            exponent_a * (x - 1) * (x**2 - u * x + accessory_v)
            + exponent_b * x * (x**2 - u * x + accessory_v)
            + x * (x - 1) * (2 * x - u)
        )
        * pow(7, -1, prime)
    )
    S = P(x + shift)
    V = (S * D * T**2).mul_ground(4 * pow(gamma**2, -1, prime) % prime)
    A = (S * E * T).mul_ground(2 * pow(gamma, -1, prime) % prime)
    boundary = S**3 * T**8 * P(x) ** extras[0] * P(x - 1) ** extras[1]
    owner_boundary = S * T

    v = int(V.nth(16)) % prime
    inverse_v = pow(v, -1, prime)
    r1, r2, r3, r4 = tuple(
        int(V.nth(16 - index)) * inverse_v % prime for index in range(1, 5)
    )
    kappa, lam, rho, sigma = normalized_invariants(r1, r2, r3, r4)
    lam %= prime
    rho %= prime
    sigma %= prime
    require(r1 != 0 and lam != 0 and rho != 0 and sigma != 0,
            f"{name} finite strict-transform units")
    a = 2 * pow(gamma, -1, prime) % prime
    r_star = -a * rho * pow(8 * lam, -1, prime) % prime
    s_star = (8 * r_star - a * kappa) * pow(4 * r1, -1, prime) % prime
    t_star = a
    B = P(1 + r_star * x**7 + s_star * x**8 + t_star * x**9)
    J = V.mul_ground(2) * B.diff() - B * V.diff()
    H = sp.exquo(base.K_general(V, A, B, J), boundary)
    controls = (
        H.degree(),
        sp.gcd(B, owner_boundary).degree(),
        sp.gcd(H, owner_boundary).degree(),
        sp.gcd(H, H.diff()).degree(),
        int(H.nth(51)) % prime,
    )
    require((t_star, s_star, r_star, controls[-1]) == expected,
            f"{name} finite tuned values")
    require(controls[:4] == (51, 0, 0, 0), f"{name} finite gcd controls")
    print(
        f"finite_case={name} good_reduction=(p={prime},u={u}) "
        f"(t_star,s_star,r_star)=({t_star},{s_star},{r_star}) "
        "controls=(degree,gcd_B_boundary,gcd_H_boundary,gcd_H_Hprime,h51)="
        f"{controls}"
    )


def source_audit() -> None:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    float_literals = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(tree)
    )
    require(assert_nodes == 0 and float_literals == 0, "source AST gate")
    print(f"source_ast=(assert_nodes={assert_nodes},float_literals={float_literals})")


def main() -> None:
    print("degree-seven retuned quartic infinity-wall wildcard")
    print(f"dependency_hash_checks={len(DEPENDENCIES)}")
    universal_top_jet_checks()
    for name in ("4111", "3211"):
        accessory_case(name)
    for name in ("4111", "3211"):
        finite_case(name)
    print("quartic_law=w^4~-128*Lambda*delta/(3*a*Sigma)=-64*Lambda*epsilon/(3*Sigma)")
    print("local_inertia=one_4_cycle_on_four_escaping_critical_roots;sign=odd")
    print("remaining_critical_points=51;branchwise_Keller_cofactor=NOT_SUPPLIED")
    print("scope=critical-resultant-strict-transform-not-polynomial-inverse-cover-not-JC2")
    source_audit()


if __name__ == "__main__":
    main()
