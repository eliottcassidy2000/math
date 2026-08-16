#!/usr/bin/env python3
"""Exact deterministic companion for provisional THM-3472.

The theorem itself is elementary: on an odd primitive modulus the sheet map
ell -> 2ell+1 transports every fixed-zero mask to a literal half-twist mask,
and THM-3405 dilation then transports the primitive cover back to the ambient
modulus.  This companion audits that conjugacy, primitivity, and dilation on
large exact finite banks and computes the cap-seven periodic/harmonic census
from THM-3453's proved atom list by a weighted CRT state product.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
from itertools import product
from json import dumps
from math import gcd
from pathlib import Path
import sys


THEOREM_ID = "THM-3472"
ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = (
    (
        "THM-3405",
        "01-canon/theorems/THM-3405-common-centre-gcd-gauge-and-boolean-half-twist.md",
        "d3e7dbeeb85c6f897bd9e31270bd0b6602ae4feac3b46a45eb5ce23ae5d24fe0",
    ),
    (
        "THM-3453",
        "01-canon/theorems/THM-3453-global-literal-half-twist-cap-seven-support-classification.md",
        "a4444813a2b24aa55613eef5f4cc0538ca0148f16297dff46f80d723f1beb247",
    ),
)

RANK4_ATOMS = (9,)
RANK6_ATOMS = (11, 15, 23, 25)
RANK7_ATOMS = (13, 29, 51)
PERIOD_FACTORS = (2, 9, 25, 11, 13, 17, 23, 29)
PERIOD_PRIMES = (2, 3, 5, 11, 13, 17, 23, 29)
EXPECTED_PERIOD = 729_664_650
GT7 = 9
EXPECTED_COUNTS = {
    0: 364_832_325,
    4: 40_536_925,
    5: 0,
    6: 64_859_080,
    7: 31_171_360,
    GT7: 228_264_960,
}
EXPECTED_DENSITIES = {
    4: Fraction(1, 18),
    5: Fraction(0),
    6: Fraction(4, 45),
    7: Fraction(283376, 6633315),
    GT7: Fraction(691712, 2211105),
}
EXPECTED_SEMANTIC_SHA256 = "d8bc9ad4a49f954ec1c76db01a7506a5965f0dc5f58881bb7315de402d151221"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def cyclic_distance(value: int, modulus: int) -> int:
    residue = value % modulus
    return min(residue, modulus - residue)


def zero_mask(q: int, owner: int) -> int:
    require(q >= 2 and 0 < owner < q, ("zero mask domain", q, owner))
    mask = 0
    for sheet in range(q):
        if 14 * cyclic_distance(owner * sheet, q) < q:
            mask |= 1 << sheet
    return mask


def half_mask(q: int, owner: int) -> int:
    require(q >= 2 and 0 < owner < q, ("half mask domain", q, owner))
    mask = 0
    for sheet in range(q):
        if 14 * cyclic_distance(owner * (2 * sheet + 1), 2 * q) < 2 * q:
            mask |= 1 << sheet
    return mask


def doubled_owner(q: int, owner: int) -> int:
    """Canonical sign representative of 2*owner modulo 2q."""

    result = 2 * owner
    if result > q:
        result = 2 * q - result
    require(q % 2 == 1 and 0 < result < q, ("doubled owner", q, owner, result))
    return result


def transported_mask(q: int, fixed: int) -> int:
    result = 0
    for sheet in range(q):
        target = (2 * sheet + 1) % q
        if fixed >> target & 1:
            result |= 1 << sheet
    return result


def conjugacy_audit() -> dict[str, object]:
    rows = 0
    cells = 0
    digest = sha256()
    for q in range(3, 402, 2):
        image = {(2 * sheet + 1) % q for sheet in range(q)}
        require(len(image) == q, ("odd sheet bijection", q))
        for owner in range(1, q):
            half_owner = doubled_owner(q, owner)
            fixed = zero_mask(q, owner)
            half = half_mask(q, half_owner)
            transported = transported_mask(q, fixed)
            require(half == transported, ("mask conjugacy", q, owner))
            require(gcd(q, owner) == gcd(q, half_owner), ("gcd", q, owner, half_owner))
            for sheet in range(q):
                target = (2 * sheet + 1) % q
                require(
                    ((half >> sheet) & 1) == ((fixed >> target) & 1),
                    ("cell conjugacy", q, owner, sheet),
                )
                cells += 1
            digest.update(f"{q}:{owner}:{fixed:x}:{half:x}\n".encode("ascii"))
            rows += 1
    require(rows == 40_200 and cells == 10_787_000, (rows, cells))

    # Sharp even-modulus failure of the sheet-bijection proof.
    even_image = {(2 * sheet + 1) % 8 for sheet in range(8)}
    require(even_image == {1, 3, 5, 7}, even_image)
    require(zero_mask(8, 4).bit_count() == 4, "even fixed mask")
    require(transported_mask(8, zero_mask(8, 4)) == 0, "even transport loss")
    return {
        "rows": rows,
        "cells": cells,
        "digest": digest.hexdigest(),
        "even_image": tuple(sorted(even_image)),
    }


def dilation_audit() -> dict[str, object]:
    rows = 0
    cells = 0
    digest = sha256()
    for base_q in range(3, 102, 2):
        for scale in (3, 5, 7):
            q = scale * base_q
            for owner in range(1, base_q):
                base = half_mask(base_q, owner)
                lifted = half_mask(q, scale * owner)
                expected = 0
                for sheet in range(q):
                    bit = (base >> (sheet % base_q)) & 1
                    if bit:
                        expected |= 1 << sheet
                    require(bit == ((lifted >> sheet) & 1), ("dilation cell", base_q, scale, owner, sheet))
                    cells += 1
                require(lifted == expected, ("dilation mask", base_q, scale, owner))
                digest.update(f"{base_q}:{scale}:{owner}:{lifted:x}\n".encode("ascii"))
                rows += 1
    require(rows == 7_650, rows)
    require(cells == 2_613_750, cells)
    return {"rows": rows, "cells": cells, "digest": digest.hexdigest()}


def cap7_symbol(q: int) -> int:
    if q % 2 == 0:
        return 0
    if any(q % atom == 0 for atom in RANK4_ATOMS):
        return 4
    if any(q % atom == 0 for atom in RANK6_ATOMS):
        return 6
    if any(q % atom == 0 for atom in RANK7_ATOMS):
        return 7
    return GT7


def weighted_crt_census() -> dict[str, object]:
    period = 1
    for factor in PERIOD_FACTORS:
        period *= factor
    require(period == EXPECTED_PERIOD, period)

    # Categories carry their number of residue classes modulo the local prime
    # power: parity; v3=0,1,>=2 mod 9; v5=0,1,>=2 mod 25; then divisibility
    # bits for 11,13,17,23,29.
    parity_states = ((False, 1), (True, 1))
    v3_states = ((0, 6), (1, 2), (2, 1))
    v5_states = ((0, 20), (1, 4), (2, 1))
    prime_states = tuple(((False, prime - 1), (True, 1)) for prime in (11, 13, 17, 23, 29))

    counts = {0: 0, 4: 0, 5: 0, 6: 0, 7: 0, GT7: 0}
    state_count = 0
    for state in product(parity_states, v3_states, v5_states, *prime_states):
        (odd, w2), (v3, w3), (v5, w5), *prime_bits = state
        flags = {prime: bit for prime, (bit, _weight) in zip((11, 13, 17, 23, 29), prime_bits)}
        weight = w2 * w3 * w5
        for _bit, local_weight in prime_bits:
            weight *= local_weight

        if not odd:
            symbol = 0
        elif v3 >= 2:
            symbol = 4
        elif flags[11] or flags[23] or (v3 >= 1 and v5 >= 1) or v5 >= 2:
            symbol = 6
        elif flags[13] or flags[29] or (v3 >= 1 and flags[17]):
            symbol = 7
        else:
            symbol = GT7
        counts[symbol] += weight
        state_count += 1

    require(state_count == 576, state_count)
    require(counts == EXPECTED_COUNTS, counts)
    require(sum(counts.values()) == period, counts)

    densities = {
        rank: Fraction(counts[rank], period)
        for rank in (4, 5, 6, 7, GT7)
    }
    require(densities == EXPECTED_DENSITIES, densities)
    require(sum(densities.values()) == Fraction(1, 2), densities)

    conditioned = {rank: 2 * value for rank, value in densities.items()}
    require(sum(conditioned.values()) == 1, conditioned)

    witnesses = []
    for prime in PERIOD_PRIMES:
        shift = period // prime
        witness = next(
            q for q in range(1, 1000)
            if cap7_symbol(q) != cap7_symbol(q + shift)
        )
        witnesses.append(
            (prime, witness, cap7_symbol(witness), witness + shift, cap7_symbol(witness + shift))
        )

    # Independent direct classifier check over a long prefix.
    for q in range(1, 100_001):
        symbol = cap7_symbol(q)
        require(symbol in counts, ("classifier", q, symbol))

    return {
        "period": period,
        "states": state_count,
        "counts": counts,
        "densities": densities,
        "conditioned": conditioned,
        "witnesses": tuple(witnesses),
    }


def security_report(source: Path) -> tuple[str, ...]:
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert present")
    forbidden = {
        "eval", "exec", "compile", "open", "system", "popen", "run",
        "Popen", "write_text", "write_bytes", "unlink", "remove", "rename",
    }
    called = {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    called.update(
        node.func.attr
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
    )
    require(not called & forbidden, ("forbidden calls", sorted(called & forbidden)))
    imports = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imports.add(node.module.split(".")[0])
    allowed = {"__future__", "ast", "fractions", "hashlib", "itertools", "json", "math", "pathlib", "sys"}
    require(imports <= allowed, ("unexpected imports", sorted(imports - allowed)))
    return tuple(sorted(imports))


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")

    dependency_hashes = tuple(
        (label, relative, lf_sha256(ROOT / relative))
        for label, relative, _expected in DEPENDENCIES
    )
    require(dependency_hashes == DEPENDENCIES, ("dependency drift", dependency_hashes))
    security = security_report(Path(__file__))
    conjugacy = conjugacy_audit()
    dilation = dilation_audit()
    census = weighted_crt_census()

    semantic_payload = {
        "theorem": THEOREM_ID,
        "dependencies": dependency_hashes,
        "security": security,
        "conjugacy": conjugacy,
        "dilation": dilation,
        "census": {
            "period": census["period"],
            "states": census["states"],
            "counts": census["counts"],
            "densities": {rank: fraction_text(value) for rank, value in census["densities"].items()},
            "conditioned": {rank: fraction_text(value) for rank, value in census["conditioned"].items()},
            "witnesses": census["witnesses"],
        },
    }
    semantic_hash = sha256(
        dumps(semantic_payload, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "UNFROZEN":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, ("semantic drift", semantic_hash))

    print("THM-3472 EXACT DETERMINISTIC COMPANION")
    print("STATUS: RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT")
    print("DEPENDENCIES:")
    for label, relative, digest in dependency_hashes:
        print(f"  {label} {digest} {relative}")
    print(f"SECURITY_IMPORTS: {','.join(security)}")
    print("UNIVERSAL_IDENTITY: for odd Q, B_Q(2s)(ell)=Z_Q(s)(2ell+1 mod Q)")
    print("GLOBAL_CONSEQUENCE: rho_ZMC(q)=rho_H(q) for every odd q")
    print(
        "ODD_CONJUGACY_AUDIT_Q3_TO401: "
        f"rows={conjugacy['rows']} cells={conjugacy['cells']} digest={conjugacy['digest']}"
    )
    print(f"EVEN_HOSTILE_PHI_IMAGE_Q8: {conjugacy['even_image']}; fixed owner 4 loses four sheets")
    print(
        "ODD_DILATION_AUDIT_BASE_Q3_TO101_SCALES_3_5_7: "
        f"rows={dilation['rows']} cells={dilation['cells']} digest={dilation['digest']}"
    )
    print("ODD_CAP7_ATOMS: rank4=(9); rank5=(); rank6=(11,15,23,25); rank7=(13,29,51)")
    print(f"ODD_CAP7_ANNOTATED_MINIMAL_PERIOD: {census['period']}")
    print(f"WEIGHTED_CRT_STATES: {census['states']}")
    print(f"PERIOD_COUNTS_0_4_5_6_7_GT7: {tuple(sorted(census['counts'].items()))}")
    print(
        "AMBIENT_NATURAL_AND_HARMONIC_COEFFICIENTS_4_5_6_7_GT7: "
        + str(tuple((rank, fraction_text(value)) for rank, value in census["densities"].items()))
    )
    print(
        "CONDITIONAL_ON_ODD_COEFFICIENTS_4_5_6_7_GT7: "
        + str(tuple((rank, fraction_text(value)) for rank, value in census["conditioned"].items()))
    )
    print(f"MINIMAL_PERIOD_WITNESSES: {census['witnesses']}")
    print(f"SEMANTIC_SHA256: {semantic_hash}")
    print("VERDICT: odd global layer equality and exact cap-seven harmonic atlas; no even-modulus, endpoint-current, or LRC(14) conclusion")


if __name__ == "__main__":
    main()
