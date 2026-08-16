#!/usr/bin/env python3
"""Exact deterministic companion for repaired and audited THM-3472.

The universal map is one-way at even modulus: a fixed-zero mask is sampled
on the image of ell -> 2ell+1 and becomes a literal doubled-owner half mask.
Odd moduli give a bijective conjugacy; even moduli still give cover transport
because the self-opposite fixed owner is empty on the odd image.  Together
with the explicitly audited active-gcd divisor interface and dilation, this
proves equality of the global zero-cochain and literal half ranks for every
q>=2.  The companion checks the identities and computes both the all-modulus
and odd cap-seven harmonic atlases from THM-3453's atom list.
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
        "THM-3415",
        "01-canon/theorems/THM-3415-zero-mode-cochain-global-rank-five-support.md",
        "de8d2f615695070dc75cad38ad4a9c22308d8bf900cd6f9a69cd2003f815a14d",
    ),
    (
        "THM-3453",
        "01-canon/theorems/THM-3453-global-literal-half-twist-cap-seven-support-classification.md",
        "a4444813a2b24aa55613eef5f4cc0538ca0148f16297dff46f80d723f1beb247",
    ),
)

RANK_ATOMS = {
    4: (8, 9),
    5: (10, 12),
    6: (11, 15, 23, 25),
    7: (13, 14, 29, 38, 51, 68, 148),
}
GT7 = 8
PRIME_EXPONENTS = ((2, 3), (3, 2), (5, 2), (7, 1), (11, 1), (13, 1),
                   (17, 1), (19, 1), (23, 1), (29, 1), (37, 1))
EXPECTED_ALL_PERIOD = 14_362_718_970_600
EXPECTED_ALL_COUNTS = {
    4: 3_191_715_326_800,
    5: 1_276_686_130_720,
    6: 1_734_627_895_000,
    7: 1_531_452_347_040,
    GT7: 6_628_237_271_040,
}
EXPECTED_ALL_DENSITIES = {
    4: Fraction(2, 9),
    5: Fraction(4, 45),
    6: Fraction(25, 207),
    7: Fraction(165741596, 1554406815),
    GT7: Fraction(717341696, 1554406815),
}
EXPECTED_ODD_PERIOD = 729_664_650
EXPECTED_ODD_COUNTS = {
    0: 364_832_325,
    4: 40_536_925,
    5: 0,
    6: 64_859_080,
    7: 31_171_360,
    GT7: 228_264_960,
}
EXPECTED_TRANSPORT_SHA256 = "3405151ea290a3724b8fc4ef419c3216789d3e8f1593264a49f630517dfb3ff4"
EXPECTED_DIVISOR_SHA256 = "98fe0dedb72963aab691506b5c4e65f309c43a8454ab61e47a1713d9785d5276"
EXPECTED_SEMANTIC_SHA256 = "116818fa2bbc5a0cada41b425f08c4b7afb9a3051e17d804463f002b1027d81a"


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


def doubled_owner(q: int, owner: int) -> int | None:
    """Sign-normalized 2*owner; None is the empty self-opposite class."""

    require(0 < owner <= q // 2, ("canonical fixed owner", q, owner))
    result = 2 * owner
    if result == q:
        return None
    require(0 < result < q, ("doubled owner", q, owner, result))
    return result


def transported_mask(q: int, fixed: int) -> int:
    result = 0
    for sheet in range(q):
        target = (2 * sheet + 1) % q
        if fixed >> target & 1:
            result |= 1 << sheet
    return result


def all_modulus_transport_audit() -> dict[str, object]:
    rows = 0
    cells = 0
    odd_rows = 0
    even_rows = 0
    empty_self_opposite = 0
    digest = sha256()
    for q in range(2, 402):
        image = {(2 * sheet + 1) % q for sheet in range(q)}
        require(len(image) == q // gcd(q, 2), ("affine image", q, image))
        for owner in range(1, q // 2 + 1):
            fixed = zero_mask(q, owner)
            transported = transported_mask(q, fixed)
            half_owner = doubled_owner(q, owner)
            if half_owner is None:
                require(q % 2 == 0 and owner == q // 2, (q, owner))
                require(transported == 0, ("self-opposite image", q, owner))
                half = 0
                empty_self_opposite += 1
            else:
                half = half_mask(q, half_owner)
                require(half == transported, ("one-way mask identity", q, owner))
            for sheet in range(q):
                target = (2 * sheet + 1) % q
                require(
                    ((half >> sheet) & 1) == ((fixed >> target) & 1),
                    ("transport cell", q, owner, sheet),
                )
                cells += 1
            digest.update(f"{q}:{owner}:{fixed:x}:{half:x}\n".encode("ascii"))
            rows += 1
            if q % 2:
                odd_rows += 1
            else:
                even_rows += 1
    require(rows == 40_200 and cells == 10_766_900, (rows, cells))
    require(odd_rows == 20_100 and even_rows == 20_100, (odd_rows, even_rows))
    require(empty_self_opposite == 200, empty_self_opposite)

    # MISTAKE-407: a valid literal cover need not be augmented-primitive.
    hostile_q = 15
    hostile_fixed = (1, 2, 3, 4, 5, 7)
    hostile_half = tuple(2 * owner for owner in hostile_fixed)
    full = (1 << hostile_q) - 1
    fixed_union = 0
    half_union = 0
    for owner in hostile_fixed:
        fixed_union |= zero_mask(hostile_q, owner)
    for owner in hostile_half:
        half_union |= half_mask(hostile_q, owner)
    require(fixed_union == full and half_union == full, "Q15 cover hostile")
    require(gcd(hostile_q, *hostile_fixed) == 1, "fixed primitive")
    require(gcd(2 * hostile_q, *hostile_half) == 2, "augmented half hostile")

    even_image = tuple(sorted((2 * sheet + 1) % 8 for sheet in range(8)))
    require(set(even_image) == {1, 3, 5, 7}, even_image)
    require(zero_mask(8, 4).bit_count() == 4, "even fixed mask")
    require(transported_mask(8, zero_mask(8, 4)) == 0, "empty on image")
    result_digest = digest.hexdigest()
    if EXPECTED_TRANSPORT_SHA256 != "UNFROZEN":
        require(result_digest == EXPECTED_TRANSPORT_SHA256, result_digest)
    return {
        "rows": rows,
        "odd_rows": odd_rows,
        "even_rows": even_rows,
        "cells": cells,
        "empty_self_opposite": empty_self_opposite,
        "digest": result_digest,
        "q15_augmented_hostile": (hostile_fixed, hostile_half, 1, 2),
        "q8_image": tuple(sorted(set(even_image))),
    }


def divisor_interface_audit() -> dict[str, object]:
    """Audit the THM-3405 reduction formula before either Boolean layer."""

    rows = 0
    cells = 0
    digest = sha256()
    for q in range(2, 81):
        for d in range(1, 2 * q + 1):
            g = gcd(q, d)
            Q = q // g
            d0 = d // g
            require(gcd(Q, d0) == 1, (q, d, Q, d0))
            for v in range(1, min(Q, 11)):
                u = d * v
                for epsilon in (0, 1):
                    centre = Fraction(epsilon * g, 2 * q * d)
                    for sheet in range(q):
                        ambient = Fraction(u) * (centre + Fraction(sheet, q))
                        if epsilon == 0:
                            reduced = Fraction(d0 * v * sheet, Q)
                        else:
                            reduced = Fraction(v * (2 * d0 * sheet + 1), 2 * Q)
                        require(ambient == reduced, ("divisor formula", q, d, v, epsilon, sheet))
                        left = 14 * cyclic_distance(ambient.numerator, ambient.denominator) < ambient.denominator
                        right = 14 * cyclic_distance(reduced.numerator, reduced.denominator) < reduced.denominator
                        require(left == right, ("divisor strict mask", q, d, v, epsilon, sheet))
                        cells += 1
                    digest.update(f"{q}:{d}:{v}:{epsilon}:{Q}:{d0}\n".encode("ascii"))
                    rows += 1
    result_digest = digest.hexdigest()
    if EXPECTED_DIVISOR_SHA256 != "UNFROZEN":
        require(result_digest == EXPECTED_DIVISOR_SHA256, result_digest)
    return {"q_max": 80, "rows": rows, "cells": cells, "digest": result_digest}


def dilation_audit() -> dict[str, object]:
    rows = 0
    cells = 0
    digest = sha256()
    for base_q in range(2, 102):
        for scale in (2, 3, 5):
            q = scale * base_q
            for owner in range(1, base_q):
                base = half_mask(base_q, owner)
                lifted = half_mask(q, scale * owner)
                expected = 0
                for sheet in range(q):
                    bit = (base >> (sheet % base_q)) & 1
                    if bit:
                        expected |= 1 << sheet
                    require(bit == ((lifted >> sheet) & 1),
                            ("dilation cell", base_q, scale, owner, sheet))
                    cells += 1
                require(lifted == expected, ("dilation mask", base_q, scale, owner))
                digest.update(f"{base_q}:{scale}:{owner}:{lifted:x}\n".encode("ascii"))
                rows += 1
    require(rows == 15_150 and cells == 3_434_000, (rows, cells))
    return {"rows": rows, "cells": cells, "digest": digest.hexdigest()}


def cap7_symbol(q: int) -> int:
    require(q >= 1, q)
    for rank in (4, 5, 6, 7):
        if any(q % atom == 0 for atom in RANK_ATOMS[rank]):
            return rank
    return GT7


def local_states(prime: int, exponent: int) -> tuple[tuple[int, int], ...]:
    states = []
    for valuation in range(exponent + 1):
        weight = 1 if valuation == exponent else (
            prime ** (exponent - valuation) - prime ** (exponent - valuation - 1)
        )
        states.append((valuation, weight))
    return tuple(states)


def atom_divides_state(atom: int, valuations: dict[int, int]) -> bool:
    value = atom
    for prime, cap in PRIME_EXPONENTS:
        needed = 0
        while value % prime == 0:
            value //= prime
            needed += 1
        if valuations[prime] < needed:
            return False
    require(value == 1, ("atom outside period", atom, value))
    return True


def all_modulus_crt_census() -> dict[str, object]:
    period = 1
    state_axes = []
    for prime, exponent in PRIME_EXPONENTS:
        period *= prime ** exponent
        state_axes.append(local_states(prime, exponent))
    require(period == EXPECTED_ALL_PERIOD, period)

    counts = {4: 0, 5: 0, 6: 0, 7: 0, GT7: 0}
    state_count = 0
    for state in product(*state_axes):
        valuations = {
            prime: valuation
            for (prime, _exponent), (valuation, _weight) in zip(PRIME_EXPONENTS, state)
        }
        weight = 1
        for _valuation, local_weight in state:
            weight *= local_weight
        symbol = GT7
        for rank in (4, 5, 6, 7):
            if any(atom_divides_state(atom, valuations) for atom in RANK_ATOMS[rank]):
                symbol = rank
                break
        counts[symbol] += weight
        state_count += 1
    require(state_count == 9_216, state_count)
    require(counts == EXPECTED_ALL_COUNTS, counts)
    require(sum(counts.values()) == period, counts)
    densities = {rank: Fraction(count, period) for rank, count in counts.items()}
    require(densities == EXPECTED_ALL_DENSITIES, densities)
    require(sum(densities.values()) == 1, densities)

    witnesses = []
    for prime, _exponent in PRIME_EXPONENTS:
        shift = period // prime
        q = next(q for q in range(1, 1000) if cap7_symbol(q) != cap7_symbol(q + shift))
        witnesses.append((prime, q, cap7_symbol(q), q + shift, cap7_symbol(q + shift)))
    return {
        "period": period,
        "states": state_count,
        "counts": counts,
        "densities": densities,
        "witnesses": tuple(witnesses),
    }


def odd_subatlas_census() -> dict[str, object]:
    period = EXPECTED_ODD_PERIOD
    parity_states = ((False, 1), (True, 1))
    v3_states = ((0, 6), (1, 2), (2, 1))
    v5_states = ((0, 20), (1, 4), (2, 1))
    prime_states = tuple(((False, prime - 1), (True, 1)) for prime in (11, 13, 17, 23, 29))
    counts = {0: 0, 4: 0, 5: 0, 6: 0, 7: 0, GT7: 0}
    states = 0
    for state in product(parity_states, v3_states, v5_states, *prime_states):
        (odd, w2), (v3, w3), (v5, w5), *bits = state
        flags = {prime: bit for prime, (bit, _weight) in zip((11, 13, 17, 23, 29), bits)}
        weight = w2 * w3 * w5
        for _bit, local_weight in bits:
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
        states += 1
    require(states == 576 and counts == EXPECTED_ODD_COUNTS, (states, counts))
    require(sum(counts.values()) == period, counts)
    return {
        "period": period,
        "states": states,
        "counts": counts,
        "ambient_densities": {
            rank: Fraction(count, period) for rank, count in counts.items()
        },
    }


def security_report(source: Path) -> tuple[str, ...]:
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert present")
    forbidden = {
        "eval", "exec", "compile", "open", "system", "popen", "run", "Popen",
        "write_text", "write_bytes", "unlink", "remove", "rename",
    }
    called = {
        node.func.id for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    called.update(
        node.func.attr for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
    )
    require(not called & forbidden, ("forbidden calls", sorted(called & forbidden)))
    imports = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imports.add(node.module.split(".")[0])
    allowed = {
        "__future__", "ast", "fractions", "hashlib", "itertools", "json", "math",
        "pathlib", "sys",
    }
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
    transport = all_modulus_transport_audit()
    divisor = divisor_interface_audit()
    dilation = dilation_audit()
    all_census = all_modulus_crt_census()
    odd_census = odd_subatlas_census()

    semantic_payload = {
        "theorem": THEOREM_ID,
        "dependencies": dependency_hashes,
        "security": security,
        "transport": transport,
        "divisor_interface": divisor,
        "dilation": dilation,
        "all_census": {
            "period": all_census["period"],
            "states": all_census["states"],
            "counts": all_census["counts"],
            "densities": {r: fraction_text(v) for r, v in all_census["densities"].items()},
            "witnesses": all_census["witnesses"],
        },
        "odd_subatlas": {
            "period": odd_census["period"],
            "states": odd_census["states"],
            "counts": odd_census["counts"],
            "ambient_densities": {
                r: fraction_text(v) for r, v in odd_census["ambient_densities"].items()
            },
        },
    }
    semantic_hash = sha256(
        dumps(semantic_payload, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "UNFROZEN":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, ("semantic drift", semantic_hash))

    print("THM-3472 EXACT DETERMINISTIC COMPANION")
    print("STATUS: PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED")
    for label, relative, digest in dependency_hashes:
        print(f"DEPENDENCY: {label} {digest} {relative}")
    print(f"SECURITY_IMPORTS: {','.join(security)}")
    print("UNIVERSAL_TRANSPORT: B_Q(2s)(ell)=Z_Q(s)(2ell+1 mod Q) for every Q; self-opposite even owner deleted")
    print("GLOBAL_CONSEQUENCE_CANDIDATE: rho_ZMC(q)=rho_H(q) for every q>=2")
    print(f"ALL_MODULUS_TRANSPORT_Q2_TO401: {transport}")
    print(f"DIVISOR_INTERFACE_Q2_TO80: {divisor}")
    print(f"DILATION_BASE_Q2_TO101_SCALES_2_3_5: {dilation}")
    print(f"ALL_CAP7_ATOMS: {RANK_ATOMS}")
    print(f"ALL_CAP7_MINIMAL_PERIOD: {all_census['period']}")
    print(f"ALL_WEIGHTED_CRT_STATES: {all_census['states']}")
    print(f"ALL_PERIOD_COUNTS_4_5_6_7_GT7: {tuple(sorted(all_census['counts'].items()))}")
    print("ALL_NATURAL_AND_HARMONIC_COEFFICIENTS: " + str(tuple(
        (rank, fraction_text(value)) for rank, value in all_census["densities"].items()
    )))
    print(f"ALL_MINIMAL_PERIOD_WITNESSES: {all_census['witnesses']}")
    print(f"ODD_SUBATLAS_PERIOD_STATES_COUNTS: {(odd_census['period'], odd_census['states'], tuple(sorted(odd_census['counts'].items())))}")
    print(f"SEMANTIC_SHA256: {semantic_hash}")
    print("VERDICT: repaired all-modulus rank equality candidate and exact cap-seven harmonic atlas; augmented primitivity, endpoint current, and LRC(14) do not follow")


if __name__ == "__main__":
    main()
