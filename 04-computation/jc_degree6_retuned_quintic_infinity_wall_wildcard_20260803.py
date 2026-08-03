#!/usr/bin/env python3
"""Exact wildcard probe below THM-3263's quartic infinity wall.

The four-parameter clutch is

    B_(q,r,s,t) = 1 + q*x^6 + r*x^7 + s*x^8 + t*x^9.

On the primitive wall, q changes the x^96 carry after s is retuned.  Retuning
r as a function of q leaves one linear x^95 carry.  Its unique zero exposes
a nonzero quintic x^94 carry in both inherited accessory fields.

This remains a critical-resultant computation.  It supplies no polynomial
inverse cover, branchwise Keller cofactor, or conclusion about JC(2).
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
    "04-computation/jc_degree7_retuned_quartic_infinity_wall_wildcard_20260803.py":
        "aabdb06854ef9045d31c4fdbe76e2933282f44b49230ce574a77591af9a6df56",
    "05-knowledge/results/jc_degree7_retuned_quartic_infinity_wall_wildcard_20260803.out":
        "2ad12c392ff4c3228a7a242cffa58449a734339e102628d3e73386883a9bd0b2",
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
    require(lf_hash(ROOT / dependency) == expected_hash,
            f"dependency drift: {dependency}")


def load_base():
    spec = importlib.util.spec_from_file_location("thm3237_degree6_dependency", BASE_PATH)
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
    theta = r1**3 - 4 * r1 * r2 + 8 * r3 + 8
    return kappa, lam, rho, sigma, theta


def universal_top_jet_checks():
    w = sp.symbols("w")
    a = sp.symbols("a", nonzero=True)
    r1, r2, r3, r4, r5 = sp.symbols("r1 r2 r3 r4 r5")
    q, r, s, t = sp.symbols("q r s t")

    vbar = a**2 * (
        1 + r1 * w + r2 * w**2 + r3 * w**3 + r4 * w**4 + r5 * w**5
    )
    # The top response identity says A/a is the square-root jet of V/a^2.
    abar = a * sp.series(sp.sqrt(vbar / a**2), w, 0, 6).removeO()
    bbar = t + s * w + r * w**2 + q * w**3
    jbar = (
        2 * vbar * bbar
        - 2 * w * vbar * sp.diff(bbar, w)
        + w * bbar * sp.diff(vbar, w)
    )

    # Terms of maximal degrees 99 and 96 are the only K_general groups which
    # can reach x^94.  The w^3 shift aligns their reciprocal degrees.
    degree99 = -4 * abar * bbar**3 * jbar**2 * vbar + 8 * bbar**3 * jbar * vbar**3
    degree96 = (
        -abar**3 * jbar**3
        + 12 * abar**2 * jbar**2 * vbar**2
        - 48 * abar * jbar * vbar**4
        + 64 * vbar**6
    )
    reciprocal = sp.series(degree99 + w**3 * degree96, w, 0, 6).removeO().expand()
    coefficients = tuple(sp.factor(reciprocal.coeff(w, index)) for index in range(6))

    kappa, lam, rho, sigma, theta = normalized_invariants(r1, r2, r3, r4)
    wall_coefficients = tuple(sp.factor(value.subs(t, a)) for value in coefficients)
    c99, c98, c97, c96, c95, c94 = wall_coefficients
    require(c99 == 0 and c98 == 0, "primitive wall top cancellations")
    require(sp.factor(c97 + 2 * a**11 * (a * kappa + 4 * r1 * s - 8 * r)) == 0,
            "degree-six preserves x97 formula")

    s_line = (8 * r - a * kappa) / (4 * r1)
    expected_c96 = -a**11 * (a * rho + 8 * r * lam - 32 * q * r1) / r1
    require(sp.factor(c96.subs(s, s_line) - expected_c96) == 0,
            "degree-six linear x96 carry")

    r_line = (32 * q * r1 - a * rho) / (8 * lam)
    s_line_q = sp.factor(s_line.subs(r, r_line))
    expected_c95 = -3 * a**11 * (a * sigma + 64 * q * theta) / (8 * lam)
    require(sp.factor(c95.subs({r: r_line, s: s_line_q}) - expected_c95) == 0,
            "retuned linear x95 carry")

    q_star = -a * sigma / (64 * theta)
    r_star = sp.factor(r_line.subs(q, q_star))
    s_star = sp.factor(s_line_q.subs(q, q_star))
    c94_star = sp.factor(c94.subs({q: q_star, r: r_star, s: s_star}))
    psi = sp.factor(-16 * theta**2 * c94_star / a**12)
    psi_poly = sp.Poly(psi, r1, r2, r3, r4, r5, domain="ZZ")
    psi_text = str(psi_poly.as_expr())
    require(len(psi_poly.terms()) == 58 and psi_poly.total_degree() == 11,
            "quintic invariant shape")
    require(sha256(psi_text.encode("ascii")).hexdigest()
            == "6b1814756894c97cfe6a31a3bd4ae51bc1dc591879211ec5f67d65f2a960e860",
            "quintic invariant digest")
    require(sp.factor(c94_star + a**12 * psi / (16 * theta**2)) == 0,
            "quintic x94 carry")

    # q=0 must recover the promoted quartic wall exactly.
    require(sp.factor(r_line.subs(q, 0) + a * rho / (8 * lam)) == 0,
            "THM-3263 r-star recovery")
    require(sp.factor(expected_c95.subs(q, 0) + 3 * a**12 * sigma / (8 * lam)) == 0,
            "THM-3263 quartic carry recovery")
    require(sp.factor(sp.diff(coefficients[0], t).subs(t, a) + 16 * a**11) == 0,
            "transverse leading derivative")

    print(
        "universal_top_jet="
        "wall_c99=wall_c98=0;"
        "wall_c97=-2*a^11*(a*Kappa+4*r1*s-8*r);"
        "retuned_c96=-a^11*(a*Rho+8*r*Lambda-32*q*r1)/r1;"
        "retuned_c95=-3*a^11*(a*Sigma+64*q*Theta)/(8*Lambda);"
        "quintic_c94=-a^12*Psi/(16*Theta^2)"
    )
    print(
        "quintic_invariant="
        "Theta=r1^3-4*r1*r2+8*r3+8;"
        "Psi_terms=58;Psi_total_degree=11;"
        "Psi_sha256=6b1814756894c97cfe6a31a3bd4ae51bc1dc591879211ec5f67d65f2a960e860"
    )
    return psi_poly


EXPECTED = {
    "4111": {
        "theta_norm": Q(-37200002176, 5359375),
        "sigma_norm": Q(-157881139404422797328384, 28722900390625),
        "psi_norm": Q(-572298646716642082095739982919499776, 301716116790771484375),
        "q_norm": Q(2107940627729338801, 186000010880000000),
        "r_norm": Q(-3018470836755969, 1860000108800000),
        "s_norm": Q(1568249062443, 1162500068000),
        "H_sha256": "7926bec47cfbfad38a1de91dfeebf84b824da3c780f478240bf66bc6e9981d26",
    },
    "3211": {
        "theta_norm": Q(249238609408, 48234375),
        "sigma_norm": Q(64559323245631953174528, 28722900390625),
        "psi_norm": Q(
            -4402683590765317259172542063793800401125376,
            1697153156948089599609375,
        ),
        "q_norm": Q(-8978749567656192625, 9187932097216512),
        "r_norm": Q(8285100359018853125, 484519856689152),
        "s_norm": Q(-4924346007190859375, 161506618896384),
        "H_sha256": "631e99ab7332025f3a338fd15919afa90cd9e542e2ed3c314dc0a4b4867015ec",
    },
}


def evaluate_integer_poly(poly, values):
    result = 0
    for monomial, coefficient in poly.terms():
        term = int(coefficient)
        for value, exponent in zip(values, monomial):
            term = term * value**exponent
        result = result + term
    return result


def accessory_case(name: str, psi_poly) -> None:
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
    r1, r2, r3, r4, r5 = tuple(V.c[16 - index] / v for index in range(1, 6))
    kappa, lam, rho, sigma, theta = normalized_invariants(r1, r2, r3, r4)
    require(r1 != 0 and lam != 0 and sigma != 0 and theta != 0,
            f"{name} strict-transform denominators")
    q_star = -a * sigma / (64 * theta)
    r_star = (32 * q_star * r1 - a * rho) / (8 * lam)
    s_star = (8 * r_star - a * kappa) / (4 * r1)
    psi = evaluate_integer_poly(psi_poly, (r1, r2, r3, r4, r5))
    require(psi != 0, f"{name} quintic invariant")

    expected = EXPECTED[name]
    require(base.norm(theta) == expected["theta_norm"], f"{name} Theta norm")
    require(base.norm(sigma) == expected["sigma_norm"], f"{name} Sigma norm")
    require(base.norm(psi) == expected["psi_norm"], f"{name} Psi norm")
    require(base.norm(q_star) == expected["q_norm"], f"{name} q-star norm")
    require(base.norm(r_star) == expected["r_norm"], f"{name} r-star norm")
    require(base.norm(s_star) == expected["s_norm"], f"{name} s-star norm")

    B = (
        base.ONE_P
        + (base.X_P**6).scalar(q_star)
        + (base.X_P**7).scalar(r_star)
        + (base.X_P**8).scalar(s_star)
        + (base.X_P**9).scalar(a)
    )
    J = V * B.derivative().scalar(2) - B * V.derivative()
    H = base.K_general(V, A, B, J).exact_quotient(boundary)
    h50 = -a**12 * psi / (16 * theta**2)
    require(H.degree == 50 and H.c[50] == h50, f"{name} quintic carry")
    H_hash = sha256(base.serialize(H)).hexdigest()
    require(H_hash == expected["H_sha256"], f"{name} H digest")

    print(
        f"accessory_case={name} "
        f"Theta_norm={rational_text(expected['theta_norm'])} "
        f"Sigma_norm={rational_text(expected['sigma_norm'])} "
        f"Psi_norm={rational_text(expected['psi_norm'])} "
        f"q_star_norm={rational_text(expected['q_norm'])} "
        f"r_star_norm={rational_text(expected['r_norm'])} "
        f"s_star_norm={rational_text(expected['s_norm'])} "
        f"degree={H.degree} H_sha256={H_hash}"
    )


def evaluate_integer_poly_mod(poly, values, prime):
    result = 0
    for monomial, coefficient in poly.terms():
        term = int(coefficient) % prime
        for value, exponent in zip(values, monomial):
            term = term * pow(value, exponent, prime) % prime
        result = (result + term) % prime
    return result


def finite_case(name: str, psi_poly) -> None:
    x = sp.symbols("x")
    if name == "4111":
        prime, u, exponent_a, exponent_b = 113, 85, 4, 1
        q_coefficients = (100, 244, 237, 44)
        accessory_v = (8 * u**2 + 9 * u + 8) * pow(7, -1, prime) % prime
        shift = 5 * (u + 1) * pow(7, -1, prime) % prime
        A0 = 80 * accessory_v**2 * (u + 1) * pow(343, -1, prime) % prime
        extras = (9, 0)
        expected = (85, 1, 43, 54, 97)
    else:
        prime, u, exponent_a, exponent_b = 101, 64, 3, 2
        q_coefficients = (75, -89, -31, 61)
        accessory_v = (24 * u**2 - 16 * u - 16) * pow(21, -1, prime) % prime
        shift = (5 * u - 4) * pow(7, -1, prime) % prime
        A0 = 9 * accessory_v**2 * (5 * u - 4) * pow(343, -1, prime) % prime
        extras = (6, 3)
        expected = (89, 75, 80, 100, 56)

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
    r1, r2, r3, r4, r5 = tuple(
        int(V.nth(16 - index)) * inverse_v % prime for index in range(1, 6)
    )
    kappa, lam, rho, sigma, theta = normalized_invariants(r1, r2, r3, r4)
    kappa %= prime
    lam %= prime
    rho %= prime
    sigma %= prime
    theta %= prime
    require(r1 != 0 and lam != 0 and sigma != 0 and theta != 0,
            f"{name} finite strict-transform units")
    a = 2 * pow(gamma, -1, prime) % prime
    q_star = -a * sigma * pow(64 * theta, -1, prime) % prime
    r_star = (32 * q_star * r1 - a * rho) * pow(8 * lam, -1, prime) % prime
    s_star = (8 * r_star - a * kappa) * pow(4 * r1, -1, prime) % prime
    psi = evaluate_integer_poly_mod(psi_poly, (r1, r2, r3, r4, r5), prime)
    require(psi != 0, f"{name} finite quintic invariant")

    B = P(1 + q_star * x**6 + r_star * x**7 + s_star * x**8 + a * x**9)
    J = V.mul_ground(2) * B.diff() - B * V.diff()
    H = sp.exquo(base.K_general(V, A, B, J), boundary)
    h50 = -a**12 * psi * pow(16 * theta**2, -1, prime) % prime
    controls = (
        H.degree(),
        sp.gcd(B, owner_boundary).degree(),
        sp.gcd(H, owner_boundary).degree(),
        sp.gcd(H, H.diff()).degree(),
        int(H.nth(50)) % prime,
    )
    require((a, s_star, r_star, q_star, h50) == expected,
            f"{name} finite tuned values")
    require(controls == (50, 0, 0, 0, h50), f"{name} finite gcd controls")
    print(
        f"finite_case={name} good_reduction=(p={prime},u={u}) "
        f"(t_star,s_star,r_star,q_star)=({a},{s_star},{r_star},{q_star}) "
        "controls=(degree,gcd_B_boundary,gcd_H_boundary,gcd_H_Hprime,h50)="
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
    print("degree-six retuned quintic infinity-wall wildcard")
    print(f"dependency_hash_checks={len(DEPENDENCIES)}")
    psi_poly = universal_top_jet_checks()
    for name in ("4111", "3211"):
        accessory_case(name, psi_poly)
    for name in ("4111", "3211"):
        finite_case(name, psi_poly)
    print("quintic_law=w^5~-256*Theta^2*delta/(a*Psi)=-128*Theta^2*epsilon/Psi")
    print("local_inertia=one_5_cycle_on_five_escaping_critical_roots;sign=even")
    print("remaining_critical_points=50;branchwise_Keller_cofactor=NOT_SUPPLIED")
    print("scope=critical-resultant-strict-transform-not-polynomial-inverse-cover-not-JC2")
    source_audit()


if __name__ == "__main__":
    main()
