#!/usr/bin/env python3
"""Exact audit for affine transverse C_0 clutches in THM-3212.

For either cubic accessory response pair, consider

    P_C(x,z) = (V(x) z^2 + z + C(x))^2 + A(x) z + e,

where C is affine and e is constant.  The script derives the localized
gradient resultant from scratch, proves the universal owner-boundary ledger,
and checks the sharp two-wall S-jet obstruction in both exact cubic fields.

This is a critical-point obstruction for one explicit first-coordinate
family.  It is not an inverse-cover construction and proves no instance of
the planar Jacobian conjecture.
"""

from __future__ import annotations

import ast
from hashlib import sha256
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    "04-computation/jc_centered_heptic_source_morse_obstruction_thm3212.py":
        "03d121e57dd2edaece7cd8693d792349f03757c6e781eb5d9d0c897fcc889448",
    "05-knowledge/results/jc_centered_heptic_source_morse_obstruction_thm3212.out":
        "729e0c7b9fa51fa5c4ac5f18f50dc4413c8a8bb7bf5f0ebf2a7709650304bc85",
    "04-computation/jc_heptic_constant_parameter_discriminant_audit_thm3212.py":
        "adeb3c548d5fe3966eefc7ec4eeadfd1410a62356eca8a6c387e39cbe8fc6122",
    "05-knowledge/results/jc_heptic_constant_parameter_discriminant_audit_thm3212.out":
        "d170cf2212848ef76722579a40b65506bedf6e65a031012ca06c27dcd1ef77bb",
}
EXPECTED_DIGESTS = (
    "c07d66a389368e57ff32cdcc1c10f134951a0f2008ed65b707033d8ea8844e8e",
    "41228f21741350a603d05efd59a11ec6157d5d0fb17b6b77f03f360cd81e78ad",
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_hash(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def critical_polynomial(V, A, C, derivative_C, variable):
    """The factor K_C in Res_y(R_1,R_2)=V^3 K_C for B=1, E'=0."""

    derivative_V = V.diff(variable)
    d = derivative_C
    return (
        A**3 * derivative_V**3
        + 8 * A**2 * C * V**2 * derivative_V**2 * d
        + 2 * A**2 * C * derivative_V**3
        + 32 * A**2 * V**5 * d**3
        + 24 * A**2 * V**3 * derivative_V * d**2
        + 24 * A**2 * V**3 * derivative_V * d
        + 4 * A**2 * V * derivative_V**2 * d
        + 128 * A * C * V**5 * d**2
        + 32 * A * C * V**3 * derivative_V * d
        + 16 * A * C * V**3 * derivative_V
        - 32 * A * V**4 * d**2
        - 48 * A * V**4 * d
        - 16 * A * V**4
        - 8 * A * V**2 * derivative_V * d
        - 4 * A * V**2 * derivative_V
        + 128 * C**2 * V**5 * d
        + 32 * C**2 * V**3 * derivative_V
        - 32 * C * V**4 * d
        - 32 * C * V**4
        - 8 * C * V**2 * derivative_V
    )


def symbolic_critical_polynomial(A, C, V, derivative_V, derivative_C):
    """Symbol-only twin used to audit critical_polynomial without .diff()."""

    d = derivative_C
    return (
        A**3 * derivative_V**3
        + 8 * A**2 * C * V**2 * derivative_V**2 * d
        + 2 * A**2 * C * derivative_V**3
        + 32 * A**2 * V**5 * d**3
        + 24 * A**2 * V**3 * derivative_V * d**2
        + 24 * A**2 * V**3 * derivative_V * d
        + 4 * A**2 * V * derivative_V**2 * d
        + 128 * A * C * V**5 * d**2
        + 32 * A * C * V**3 * derivative_V * d
        + 16 * A * C * V**3 * derivative_V
        - 32 * A * V**4 * d**2
        - 48 * A * V**4 * d
        - 16 * A * V**4
        - 8 * A * V**2 * derivative_V * d
        - 4 * A * V**2 * derivative_V
        + 128 * C**2 * V**5 * d
        + 32 * C**2 * V**3 * derivative_V
        - 32 * C * V**4 * d
        - 32 * C * V**4
        - 8 * C * V**2 * derivative_V
    )


def universal_resultant_and_reduction() -> None:
    y, A, B, C, J, K, V, e = sp.symbols("y A B C J K V e")
    derivative_A, derivative_B, derivative_C, derivative_V = sp.symbols(
        "Ap Bp Cp Vp"
    )
    L = y**2 + B * y + C * V
    R1 = 2 * L * (2 * y + B) + V * A
    R2 = V**3 * e + V**2 * y + L * (J * y + K)

    resultant = sp.resultant(R1, R2, y)
    F = sp.cancel(resultant / V**2)
    require(sp.expand(resultant - V**2 * F) == 0, "universal V^2 factor")
    require(
        len(sp.Poly(F, A, B, C, J, K, V, e).terms()) == 40,
        "universal resultant term count",
    )

    raw = (
        2 * L * (
            derivative_V * y**2
            + derivative_B * V * y
            + derivative_C * V**2
        )
        + derivative_A * y * V**2
        + e * V**3
    )
    response_A_prime = (2 * V + A * derivative_V) / (2 * V)
    covariant_J = 2 * V * derivative_B - B * derivative_V
    reduced = R2.subs(
        {J: covariant_J, K: 2 * V**2 * derivative_C}
    )
    reduction_defect = sp.cancel(
        raw.subs(derivative_A, response_A_prime)
        - reduced
        - derivative_V * y * R1 / 2
    )
    require(reduction_defect == 0, "localized gradient reduction sign")

    specialized_F = sp.expand(
        F.subs({B: 1, J: -derivative_V, K: 2 * V**2 * derivative_C, e: 0})
    )
    specialized_K = sp.cancel(specialized_F / V)
    expected_K = symbolic_critical_polynomial(
        A, C, V, derivative_V, derivative_C
    )
    require(
        sp.expand(specialized_K - expected_K) == 0,
        "specialized critical factor",
    )
    require(
        len(sp.Poly(expected_K, A, C, V, derivative_V, derivative_C).terms())
        == 20,
        "specialized critical term count",
    )

    degree_bounds = (
        69, 79, 62, 96, 79, 79, 62, 89, 72, 72,
        72, 72, 72, 55, 55, 82, 65, 65, 65, 48,
    )
    require(max(degree_bounds) == 96, "affine C degree bound")
    require(
        tuple(i for i, degree in enumerate(degree_bounds) if degree == 96) == (3,),
        "unique affine C infinity term",
    )
    print(
        "universal=Res_y(R1,R2)=V^3*K_C;"
        "terms=(F40,K20);affine_nonconstant_degK=96;"
        "unique_top=32*A_lead^2*V_lead^5*d^3"
    )


def universal_boundary_checks() -> None:
    t, v, c, d = sp.symbols("t v c d", nonzero=True)
    leading_rows = []
    for multiplicity in (3, 4, 5, 6):
        V = v * t**multiplicity
        A = sp.Rational(2, 2 - multiplicity) * t
        C = c + d * t
        K = sp.Poly(sp.expand(critical_polynomial(V, A, C, d, t)), t)
        target = 3 * multiplicity - 1
        expected = (
            sp.Rational(
                32 * multiplicity * (multiplicity - 1),
                (multiplicity - 2) ** 2,
            )
            * c
            * v**3
        )
        require(K.nth(target) == expected, f"T leading row m={multiplicity}")
        require(
            all(K.nth(index) == 0 for index in range(target)),
            f"T order row m={multiplicity}",
        )
        leading_rows.append((multiplicity, target))

    v1, v2, v3, v4 = sp.symbols("v1 v2 v3 v4", nonzero=True)
    V = v1 * t + v2 * t**2 + v3 * t**3 + v4 * t**4
    A = (
        2 * t
        + 2 * v2 / (3 * v1) * t**2
        + 4 * (3 * v1 * v3 - v2**2) / (15 * v1**2) * t**3
        + 2 * (45 * v1**2 * v4 - 29 * v1 * v2 * v3 + 8 * v2**3)
        / (105 * v1**3)
        * t**4
    )
    response = sp.Poly(
        sp.cancel(2 * V * sp.diff(A, t) - A * sp.diff(V, t) - 2 * V), t
    )
    require(
        all(response.nth(index) == 0 for index in range(1, 5)),
        "simple-boundary response jet",
    )
    C = c + d * t
    K = sp.Poly(sp.expand(critical_polynomial(V, A, C, d, t)), t)
    q3 = sp.factor(K.nth(3))
    expected_q3 = sp.Rational(32, 3) * c * v1**2 * (3 * c * v1**2 + 2 * v2)
    require(q3 == expected_q3, "first S wall")
    wall_c = -2 * v2 / (3 * v1**2)
    q4_on_first_wall = sp.factor(K.nth(4).subs(c, wall_c))
    expected_q4 = (
        -sp.Rational(64, 135)
        * v2
        / v1
        * (90 * d * v1**3 - 45 * v1**3 + 54 * v1 * v3 - 28 * v2**2)
    )
    require(q4_on_first_wall == expected_q4, "second S wall")
    wall_d = (
        45 * v1**3 - 54 * v1 * v3 + 28 * v2**2
    ) / (90 * v1**3)
    obstruction = (
        315 * v1**3 * v2
        + 360 * v1**2 * v4
        - 414 * v1 * v2 * v3
        + 148 * v2**3
    )
    q5_on_both_walls = sp.factor(K.nth(5).subs({c: wall_c, d: wall_d}))
    require(
        q5_on_both_walls == -32 * v2 * obstruction / (315 * v1**2),
        "untunable S coefficient",
    )
    print(
        f"T_boundary_rows={tuple(leading_rows)};"
        "lead=32*m*(m-1)/(m-2)^2*c*v^3"
    )
    print(
        "S_boundary=ord(K_C)>=3;"
        "wall1:c=-2*v2/(3*v1^2);"
        "wall2:d=(45*v1^3-54*v1*v3+28*v2^2)/(90*v1^3)"
    )
    print("S_untunable=q5=-32*v2*R_S/(315*v1^2)")


def shifted_coefficients(poly, root, field, symbol, count):
    shifted_ring = field.poly_ring(symbol)
    t = shifted_ring.gens[0]
    shifted = shifted_ring.zero
    for (exponent,), coefficient in poly.to_dict().items():
        shifted += coefficient * (t + root) ** exponent
    entries = shifted.to_dict()
    return [entries.get((index,), field.zero) for index in range(count)]


def canonical_poly_text(poly) -> str:
    return "|".join(
        f"{degree}:{coefficient}"
        for (degree,), coefficient in sorted(poly.to_dict().items())
    )


def build_case(name: str) -> str:
    u, x, t = sp.symbols("u x t")
    if name == "4111":
        accessory = sp.Poly(
            100 * u**3 + 244 * u**2 + 237 * u + 44, u, domain=QQ
        )
        exponent_a, exponent_b = 4, 1
    else:
        accessory = sp.Poly(
            75 * u**3 - 89 * u**2 - 31 * u + 61, u, domain=QQ
        )
        exponent_a, exponent_b = 3, 2

    field = QQ.alg_field_from_poly(accessory, alias="u")
    alpha = field.ext
    ring = field.poly_ring(x)
    X = ring.gens[0]
    if name == "4111":
        accessory_v = (8 * alpha**2 + 9 * alpha + 8) / 7
        shift = 5 * (alpha + 1) / 7
        A0 = 80 * accessory_v**2 * (alpha + 1) / 343
        extras = (9, 0)
    else:
        accessory_v = (24 * alpha**2 - 16 * alpha - 16) / 21
        shift = (5 * alpha - 4) / 7
        A0 = 9 * accessory_v**2 * (5 * alpha - 4) / 343
        extras = (6, 3)

    gamma = -7 * A0
    quadratic = X**2 - alpha * X + accessory_v
    D = X**exponent_a * (X - 1) ** exponent_b * quadratic
    T = X * (X - 1) * quadratic
    E = (
        exponent_a * (X - 1) * quadratic
        + exponent_b * X * quadratic
        + X * (X - 1) * (2 * X - alpha)
    ) / 7
    S = X + shift
    V = 4 * S * D * T**2 / gamma**2
    A = 2 * S * E * T / gamma
    require(V.degree() == 16 and A.degree() == 8, f"{name} response degrees")
    require(
        2 * V * A.diff(X) - A * V.diff(X) == 2 * V,
        f"{name} response identity",
    )

    V_shift = shifted_coefficients(V, -shift, field, t, 7)
    v1, v2, v3, v4 = V_shift[1:5]
    require(V_shift[0] == field.zero and v1 != field.zero, f"{name} simple S")
    wall_c = -2 * v2 / (3 * v1**2)
    wall_d = (
        45 * v1**3 - 54 * v1 * v3 + 28 * v2**2
    ) / (90 * v1**3)
    obstruction = (
        315 * v1**3 * v2
        + 360 * v1**2 * v4
        - 414 * v1 * v2 * v3
        + 148 * v2**3
    )
    require(
        v1 != field.zero
        and v2 != field.zero
        and wall_c != field.zero
        and wall_d != field.zero
        and obstruction != field.zero,
        f"{name} exact S obstruction",
    )

    C = wall_c + wall_d * S
    require(C.gcd(S * T).degree() == 0, f"{name} wall control unit on ST")
    K = critical_polynomial(V, A, C, wall_d, X)
    require(K.degree() == 96, f"{name} critical degree")
    require(
        K.LC == 32 * A.LC**2 * V.LC**5 * wall_d**3,
        f"{name} unique top coefficient",
    )
    boundary = S**3 * T**8 * X ** extras[0] * (X - 1) ** extras[1]
    require(boundary.degree() == 44, f"{name} boundary degree")
    H = K.exquo(boundary)
    require(H.degree() == 52, f"{name} residual degree")
    require(H.gcd(T).degree() == 0, f"{name} no T escape")
    H_shift = shifted_coefficients(H, -shift, field, t, 4)
    require(
        H_shift[0] == field.zero
        and H_shift[1] == field.zero
        and H_shift[2] != field.zero,
        f"{name} sharp S order two",
    )
    K_shift = shifted_coefficients(K, -shift, field, t, 6)
    require(
        K_shift[5] == -32 * v2 * obstruction / (315 * v1**2),
        f"{name} exact untunable coefficient",
    )

    digest = sha256(canonical_poly_text(H.monic()).encode("ascii")).hexdigest()
    print(
        f"case={name} jets=(v1_nonzero,v2_nonzero,R_S_nonzero)=(1,1,1);"
        "wall=(c_nonzero,d_nonzero)=(1,1)"
    )
    print(
        f"case={name} wall_control=unit_mod_ST;degrees=(K,H)=(96,52);"
        f"ord_S(H)=2;gcd(H,T)=1;digest={digest}"
    )
    print(f"case={name} consequence=off_owner_resultant_multiplicity>=50")
    return digest


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
    print("affine transverse C0 clutch critical no-go")
    for dependency, expected_hash in DEPENDENCIES.items():
        require(
            lf_hash(ROOT / dependency) == expected_hash,
            f"dependency drift: {dependency}",
        )
    print(f"dependency_hash_checks={len(DEPENDENCIES)}")
    universal_resultant_and_reduction()
    universal_boundary_checks()
    digests = tuple(build_case(name) for name in ("4111", "3211"))
    require(digests == EXPECTED_DIGESTS, "case digest drift")
    print(f"case_digests={digests}")
    print(
        "scope=all_affine_C0_with_B=1_and_constant_E0;"
        "constant_C0_inherited_THM3212;not_inverse_cover_not_JC2"
    )
    source_audit()


if __name__ == "__main__":
    main()
