#!/usr/bin/env python3
"""Independent FLINT audit of the weighted odd-family m=3 coordinate atlas.

This companion imports neither the Sympy candidate nor any repository
computation.  It reconstructs the determinant-one map, inverse quintic,
branch discriminant, all three coordinate resultants, all three discriminant
index squares, and the target pullback in python-flint.  A separate direct
F_31 split-fibre search supplies an all-coordinate separator witness.

Scope: the explicit THM-3448 cyclic weighted subfamily at m=3.  The inherited
all-m Jelonek statement remains the structural theorem THM-3448; no identity
with the unstored historical THM-1605 family is asserted.
"""

from __future__ import annotations

import ast
from hashlib import sha256
from itertools import combinations
from pathlib import Path

from flint import fmpq, fmpq_mpoly_ctx


EXPECTED_SEMANTIC_SHA256 = "bd208dad9732439dfa14a794ef54dbfe57d66360f5f5a39b26353ffb82b6bba3"
EXPECTED_MAP_HASHES = (
    "01c6c9d962c66bab431ff42f55062a6c1d4f91762002ff73c0fe1fd14bea57a1",
    "e61a42ed722980bcf4e121a1a893bfe1fb67e3cc51409e23a5f64aac26b2f77f",
    "d2146188cbaaaac3d33b121d9800a49aebad5d84449f27a21fbd377d9b063e7e",
)
EXPECTED_BRANCH_HASH = "71dc7c97764eaf08b005c35ce60e4576e5edcc83a5d0f86328ce2d4d3a49f40e"
EXPECTED_COORDINATES = {
    "x": {
        "terms": 17,
        "total_degree": 10,
        "degrees": (5, 5, 4, 5),
        "eliminant": "f87f18f30f3d96d608b5c00b782320aa913008fc956247c269f5d64eac727d65",
        "core_terms": 17,
        "core_degree": 6,
        "core_degrees": (6, 6),
        "core": "9c793ac55f302a53713f9d6559eea0a9fcc2ac559dd539e595500207c63cde10",
        "c_exponent": 10,
        "discriminant": "580f96c990cec73828be9a4a1cfd969016e21ffacfcba7ed2bc79994ed140ebc",
    },
    "y": {
        "terms": 29,
        "total_degree": 10,
        "degrees": (5, 5, 4, 5),
        "eliminant": "36ab976ebf5830c9af1c51f8dab5aa7b71054e3b6670cc082c0a7514a2873436",
        "core_terms": 26,
        "core_degree": 6,
        "core_degrees": (6, 6),
        "core": "8e76afe528328bcca40f49bf5ef10ee64b33fb27c759046336bda31a26fe64c4",
        "c_exponent": 10,
        "discriminant": "882a5ad8fa76f575fd8c2a3ac374bfb74d25052cea425abe9f2f99cc79d5d26b",
    },
    "z": {
        "terms": 191,
        "total_degree": 15,
        "degrees": (5, 15, 12, 10),
        "eliminant": "b1ca9211c75046e9fab474e4fb842b8fdfe384d820b6efe73436add8af852b75",
        "core_terms": 268,
        "core_degree": 26,
        "core_degrees": (26, 22),
        "core": "1191484cf4a2a0c87ec7994f56ff368444ce72184f21f46d042ea16838e97e5b",
        "c_exponent": 20,
        "discriminant": "d2ea53ca388a5ee3b4b0f87fb9ab2bf76f6ab11b0486fc8d500779cc9d42fae7",
    },
}


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def digest(poly: object, indices: tuple[int, ...]) -> str:
    projected = [
        (
            tuple(monomial[index] for index in indices),
            int(coefficient.p),
            int(coefficient.q),
        )
        for monomial, coefficient in poly.terms()
    ]
    payload = tuple(sorted(projected, reverse=True))
    return sha256(repr(payload).encode("ascii")).hexdigest()


def term_count(poly: object) -> int:
    return sum(1 for _ in poly.terms())


def exact_quotient(numerator: object, denominator: object, label: str) -> object:
    quotient, remainder = divmod(numerator, denominator)
    require(remainder == 0, (label, "nonzero remainder", remainder))
    return quotient


SRC = fmpq_mpoly_ctx.get(["x", "y", "z"])
sx, sy, sz = SRC.gens()


def source_map() -> tuple[object, ...]:
    a_value = -fmpq(7, 6)
    unit = 1 + sx * sy
    gamma = 1 + a_value * sx * sy + sx**2 * sz
    alpha_numerator = unit + 3 * unit**4 * gamma**2 - 4 * unit**5 * gamma**3
    beta_numerator = 1 + 4 * unit**3 * gamma**2 - 5 * unit**4 * gamma**3
    alpha = exact_quotient(alpha_numerator, sx**2, "alpha polynomiality")
    beta = exact_quotient(beta_numerator, sx, "beta polynomiality")
    third = sx * gamma
    fmap = (alpha, beta, third)
    require(tuple(value.total_degree() for value in fmap) == (17, 16, 4),
            "map degrees")
    require(tuple(value.degrees()[2] for value in fmap) == (3, 3, 1),
            "z degrees")
    map_hashes = tuple(digest(value, (0, 1, 2)) for value in fmap)
    require(map_hashes == EXPECTED_MAP_HASHES, ("map hashes", map_hashes))

    jacobian = [[value.derivative(index) for index in range(3)] for value in fmap]
    determinant = (
        jacobian[0][0] * (
            jacobian[1][1] * jacobian[2][2]
            - jacobian[1][2] * jacobian[2][1]
        )
        - jacobian[0][1] * (
            jacobian[1][0] * jacobian[2][2]
            - jacobian[1][2] * jacobian[2][0]
        )
        + jacobian[0][2] * (
            jacobian[1][0] * jacobian[2][1]
            - jacobian[1][1] * jacobian[2][0]
        )
    )
    require(determinant == 1, ("Jacobian", determinant))
    source_w = unit * gamma
    source_p = beta * third
    source_q = alpha * third**2
    require(source_w**5 - source_w**4 + source_p * source_w - source_q == 0,
            "inverse quintic identity")
    return fmap, gamma, source_w, map_hashes


R = fmpq_mpoly_ctx.get(["w", "P", "Q", "C", "X", "Y", "Z"])
w, P, Q, C, X, Y, Z = R.gens()


def coordinate_row(
    label: str,
    relation: object,
    coordinate_index: int,
    coordinate_digest_indices: tuple[int, ...],
    inverse_equation: object,
    branch_discriminant: object,
) -> tuple[object, ...]:
    expected = EXPECTED_COORDINATES[label]
    eliminant = inverse_equation.resultant(relation, 0)
    scalar, factors = eliminant.factor()
    require(len(factors) == 1 and factors[0][1] == 1,
            (label, "irreducibility", scalar, factors))
    degrees = tuple(eliminant.degrees()[index] for index in coordinate_digest_indices)
    eliminant_digest = digest(eliminant, coordinate_digest_indices)
    require(
        (
            term_count(eliminant),
            eliminant.total_degree(),
            degrees,
            eliminant_digest,
        )
        == (
            expected["terms"],
            expected["total_degree"],
            expected["degrees"],
            expected["eliminant"],
        ),
        (label, "eliminant ledger"),
    )

    coordinate_discriminant = eliminant.discriminant(coordinate_index)
    require(coordinate_discriminant != 0, (label, "zero discriminant"))
    discriminant_digest = digest(coordinate_discriminant, (1, 2, 3))
    require(discriminant_digest == expected["discriminant"],
            (label, "discriminant hash", discriminant_digest))
    ratio = exact_quotient(
        coordinate_discriminant, branch_discriminant, label + " branch quotient"
    )
    index_root = ratio.sqrt()
    require(index_root * index_root == ratio, (label, "nonsquare index"))
    index_core = exact_quotient(
        index_root, C ** expected["c_exponent"], label + " C-index"
    )
    core_digest = digest(index_core, (1, 2))
    core_degrees = tuple(index_core.degrees()[index] for index in (1, 2))
    require(
        (
            term_count(index_core),
            index_core.total_degree(),
            core_degrees,
            core_digest,
        )
        == (
            expected["core_terms"],
            expected["core_degree"],
            expected["core_degrees"],
            expected["core"],
        ),
        (label, "index-core ledger"),
    )
    require(branch_discriminant.gcd(index_core) == 1,
            (label, "branch/index common factor"))
    return (
        label,
        str(scalar),
        term_count(eliminant),
        eliminant.total_degree(),
        degrees,
        eliminant_digest,
        term_count(index_core),
        index_core.total_degree(),
        core_degrees,
        core_digest,
        expected["c_exponent"],
        discriminant_digest,
    )


def inverse_atlas() -> tuple[object, ...]:
    inverse_equation = w**5 - w**4 + P * w - Q
    gamma = P - 4 * w**3 + 5 * w**4
    require(inverse_equation.derivative(0) == gamma, "T prime")
    branch_discriminant = inverse_equation.discriminant(0)
    scalar, factors = branch_discriminant.factor()
    require(scalar == 1 and len(factors) == 1 and factors[0][1] == 1,
            ("branch factor", scalar, factors))
    branch_digest = digest(branch_discriminant, (1, 2))
    require(branch_digest == EXPECTED_BRANCH_HASH,
            ("branch digest", branch_digest))
    a_value = -fmpq(7, 6)
    rows = (
        coordinate_row(
            "x", X * gamma - C, 4, (4, 1, 2, 3),
            inverse_equation, branch_discriminant,
        ),
        coordinate_row(
            "y", C * Y - w + gamma, 5, (5, 1, 2, 3),
            inverse_equation, branch_discriminant,
        ),
        coordinate_row(
            "z",
            C**2 * Z - gamma * (gamma * (gamma - 1 + a_value) - a_value * w),
            6,
            (6, 1, 2, 3),
            inverse_equation,
            branch_discriminant,
        ),
    )
    flat = inverse_equation.resultant(X - P, 0)
    require(flat == (X - P) ** 5 and flat.discriminant(4) == 0,
            "flat-coordinate hostile")
    return inverse_equation, branch_discriminant, branch_digest, rows


TGT = fmpq_mpoly_ctx.get(["A", "B", "C"])
tA, tB, tC = TGT.gens()


def pullback_atlas() -> tuple[object, ...]:
    target_p = tB * tC
    target_q = tA * tC**2
    pullback = (
        256 * target_p**5 - 27 * target_p**4
        - 36 * target_p**3 * target_q
        - 50 * target_p**2 * target_q**2
        - 2500 * target_p * target_q**3
        + 3125 * target_q**4 + 256 * target_q**3
    )
    l5 = exact_quotient(pullback, tC**4, "D5 pullback")
    expected_l5 = (
        3125 * tA**4 * tC**4 - 2500 * tA**3 * tB * tC**3
        + 256 * tA**3 * tC**2 - 50 * tA**2 * tB**2 * tC**2
        - 36 * tA * tB**3 * tC + 256 * tB**5 * tC - 27 * tB**4
    )
    require(l5 == expected_l5, "L5 formula")
    scalar, factors = l5.factor()
    require(scalar != 0 and len(factors) == 1 and factors[0][1] == 1,
            ("L5 irreducibility", scalar, factors))
    l5_digest = digest(l5, (0, 1, 2))
    require(l5_digest == "1557f5747f0df337384f864227850aafb7264236c1c5283e98f9eaca698d4e8e",
            ("L5 digest", l5_digest))
    return term_count(l5), l5.total_degree(), l5.degrees(), l5_digest


def evaluate_mod(poly: object, values: tuple[int, ...], prime: int) -> int:
    answer = 0
    for monomial, coefficient in poly.terms():
        term = int(coefficient.p) * pow(int(coefficient.q), prime - 2, prime)
        for value, exponent in zip(values, monomial):
            term = term * pow(value, exponent, prime) % prime
        answer = (answer + term) % prime
    return answer


def vandermonde(values: tuple[int, ...], prime: int) -> int:
    answer = 1
    for left, right in combinations(range(len(values)), 2):
        answer = answer * (values[right] - values[left]) % prime
    return answer


def finite_separator(fmap: tuple[object, ...]) -> tuple[object, ...]:
    prime = 31
    inverse_six = pow(6, prime - 2, prime)
    a_value = -7 * inverse_six % prime
    for p_value in range(prime):
        for q_value in range(prime):
            roots = tuple(
                value for value in range(prime)
                if (
                    value**5 - value**4 + p_value * value - q_value
                ) % prime == 0
            )
            if len(roots) != 5:
                continue
            rows = []
            lawful = True
            for w_value in roots:
                gamma = (
                    p_value - 4 * w_value**3 + 5 * w_value**4
                ) % prime
                if gamma == 0:
                    lawful = False
                    break
                x_value = pow(gamma, prime - 2, prime)
                y_value = (w_value - gamma) % prime
                z_value = (
                    gamma
                    * (gamma * (gamma - 1 + a_value) - a_value * w_value)
                ) % prime
                image = tuple(
                    evaluate_mod(component, (x_value, y_value, z_value), prime)
                    for component in fmap
                )
                if image != (q_value, p_value, 1):
                    lawful = False
                    break
                rows.append((x_value, y_value, z_value, w_value, gamma))
            if not lawful:
                continue
            coordinate_values = tuple(
                tuple(row[index] for row in rows) for index in range(3)
            )
            if all(len(set(values)) == 5 for values in coordinate_values):
                vandermondes = tuple(
                    vandermonde(values, prime) for values in coordinate_values
                )
                require(all(vandermondes), ("zero Vandermonde", vandermondes))
                result = (
                    prime,
                    (q_value, p_value, 1),
                    p_value,
                    q_value,
                    roots,
                    tuple(rows),
                    coordinate_values,
                    vandermondes,
                )
                require(
                    result
                    == (
                        31,
                        (26, 23, 1),
                        23,
                        26,
                        (8, 9, 12, 16, 18),
                        (
                            (28, 29, 19, 8, 10),
                            (10, 12, 15, 9, 28),
                            (7, 3, 13, 12, 9),
                            (11, 30, 11, 16, 17),
                            (6, 23, 0, 18, 26),
                        ),
                        (
                            (28, 10, 7, 11, 6),
                            (29, 12, 3, 30, 23),
                            (19, 15, 13, 11, 0),
                        ),
                        (28, 14, 1),
                    ),
                    ("finite witness drift", result),
                )
                return result
    raise RuntimeError("no all-coordinate split separator")


fmap, _source_gamma, _source_w, map_hashes = source_map()
inverse_equation, branch_discriminant, branch_digest, coordinate_rows = inverse_atlas()
pullback_row = pullback_atlas()
separator = finite_separator(fmap)

source = Path(__file__).resolve()
security_tree = ast.parse(source.read_text(encoding="utf-8"))
require(not any(isinstance(node, ast.Assert) for node in ast.walk(security_tree)),
        "assert node")
imported_modules = tuple(
    alias.name
    for node in ast.walk(security_tree)
    if isinstance(node, ast.Import)
    for alias in node.names
)
require("sympy" not in imported_modules, ("Sympy import", imported_modules))

semantic = (
    map_hashes,
    str(inverse_equation),
    str(branch_discriminant),
    branch_digest,
    coordinate_rows,
    pullback_row,
    separator,
    "all three m3 coordinate resultants irreducible quintics",
    "all three discriminants equal D5 times exact index squares",
    "F31 split fibre separates x,y,z simultaneously",
    "pullback D5=C^4 L5; sign class misses the inherited C3 Jelonek component",
)
semantic_hash = sha256(repr(semantic).encode("utf-8")).hexdigest()
if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
    require(semantic_hash == EXPECTED_SEMANTIC_SHA256,
            ("semantic hash", semantic_hash))

print("WEIGHTED ODD-FAMILY M3 THREE-COORDINATE INDEPENDENT AUDIT")
print("status=VERIFIED-EXACT independent FLINT reconstruction; explicit THM-3448 cyclic subfamily only")
print("candidate_imported=False; sympy_imported=False")
print(f"map=(term_counts={(28,25,3)},degrees={(17,16,4)},jacobian=1,hashes={map_hashes})")
print(f"inverse=(T={inverse_equation},branch_sha256={branch_digest})")
for row in coordinate_rows:
    print(f"coordinate_row={row}")
print(f"pullback=(D5=C^4*L5,L5_stats={pullback_row})")
print(f"finite_separator=(field=F_{separator[0]},target={separator[1]},P={separator[2]},Q={separator[3]},roots={separator[4]},rows={separator[5]},vandermonde_xyz={separator[7]})")
print("hostile=flat base-field view has eliminant (X-P)^5 and zero discriminant")
print("consequence=independent exact audit of the 191-term z quintic, 268-term index core, and common [D5]=[L5] coordinate class")
print("boundary=V(C) remains genuine sign-blind C3 Jelonek data inherited from THM-3448; no historical-family identity, arbitrary Keller, JC2, or LRC claim")
print(f"security_ast=({len(tuple(ast.walk(security_tree)))},assert_free=True)")
print(f"semantic_sha256={semantic_hash}")
print("STATUS=PASS")
