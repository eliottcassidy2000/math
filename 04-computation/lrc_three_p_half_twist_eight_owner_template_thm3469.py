#!/usr/bin/env python3
"""Exact deterministic companion for audited THM-3469.

The universal proof in the theorem reduces every odd phase to its signed
remainder from the nearest multiple of p.  This companion checks the symbolic
owner table and strict gap inequalities, audits the exact cover/failure
boundary for a large finite range, reattaches the common half-twist modes, and
checks the odd-modulus zero/half-layer conjugacy needed for the lower bound.
It then computes the exact periodic rank and harmonic-density census on
p=14k-1 and the annotated U-spine intersection.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
from json import dumps
from math import gcd
from pathlib import Path
import sys


THEOREM_ID = "THM-3469"
ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = (
    (
        "THM-3405",
        "01-canon/theorems/THM-3405-common-centre-gcd-gauge-and-boolean-half-twist.md",
        "d3e7dbeeb85c6f897bd9e31270bd0b6602ae4feac3b46a45eb5ce23ae5d24fe0",
    ),
    (
        "THM-3455",
        "01-canon/theorems/THM-3455-berggren-q-spine-cap-seven-atom-sieve-and-fibonacci-rank-spectrum.md",
        "e357aaf8036f485c4abbfa5969a0d1f9761372dcd7337b658daaadf4cd049c09",
    ),
)

RANK4_ATOMS = (8, 9)
RANK5_ATOMS = (10, 12)
RANK6_ATOMS = (11, 15, 23, 25)
RANK7_ATOMS = (13, 14, 29, 38, 51, 68, 148)
RANK6_K = ((5, 4), (11, 4), (23, 5))
RANK7_K = ((13, 1), (17, 11), (29, 27))
PERIOD_PRIMES = (3, 5, 11, 13, 17, 23, 29)
EXPECTED_PERIOD = 24_322_155
EXPECTED_COUNTS = {4: 8_107_385, 5: 0, 6: 4_934_930, 7: 1_818_080, 8: 9_461_760}
EXPECTED_U_SPINE_PERIOD = 11_781
EXPECTED_U_SPINE_COUNTS = {4: 748, 6: 272, 7: 144, 8: 1_080}
EXPECTED_SEMANTIC_SHA256 = "d8c6b1c2baf1df178ca09793b3f69bc5e537a12c77e9397d83fd3ce8c3e3e530"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def template(p: int) -> tuple[int, ...]:
    return (1, p - 1, p + 1, 2 * p - 1, 2 * p, 2 * p + 1, 3 * p - 6, 3 * p - 1)


def cyclic_distance(word: int, modulus: int) -> int:
    residue = word % modulus
    return min(residue, modulus - residue)


def direct_mask(q: int, owner: int) -> int:
    require(q >= 2 and 0 < owner < q, ("mask domain", q, owner))
    result = 0
    for sheet in range(q):
        word = owner * (2 * sheet + 1)
        if 14 * cyclic_distance(word, 2 * q) < 2 * q:
            result |= 1 << sheet
    return result


def zero_mask(q: int, owner: int) -> int:
    """Literal fixed-zero danger mask on Z/qZ."""

    require(q >= 2 and 0 < owner < q, ("zero-mask domain", q, owner))
    result = 0
    for sheet in range(q):
        if 14 * cyclic_distance(owner * sheet, q) < q:
            result |= 1 << sheet
    return result


def fraction_mask(q: int, owner: int) -> int:
    result = 0
    for sheet in range(q):
        value = Fraction(owner * (2 * sheet + 1), 2 * q)
        residue = value - value.numerator // value.denominator
        if min(residue, 1 - residue) < Fraction(1, 14):
            result |= 1 << sheet
    return result


def union_mask(q: int, owners: tuple[int, ...]) -> int:
    result = 0
    for owner in owners:
        result |= direct_mask(q, owner)
    return result


def predicted_cover(p: int) -> bool:
    require(p >= 13, p)
    return p % 42 not in (7, 35)


def threshold(p: int) -> int:
    """Largest integral numerator distance satisfying the strict danger test."""

    return (3 * p - 1) // 7


def mode_at_half(q: int, owner: int) -> dict[str, object]:
    divisor = gcd(q, owner)
    phase_modulus = q // divisor
    multiplier = owner // divisor % phase_modulus
    mask = direct_mask(q, owner)
    phase_block = {
        multiplier * sheet % phase_modulus
        for sheet in range(q)
        if mask >> sheet & 1
    }
    size = len(phase_block)
    require(size > 0, ("empty active mode", q, owner))
    starts = tuple(
        start
        for start in range(phase_modulus)
        if {(start + offset) % phase_modulus for offset in range(size)} == phase_block
    )
    require(len(starts) == 1, ("mode start", q, owner, starts))
    start = starts[0]
    h = (-divisor * (2 * start + size - 1)) % (2 * q)
    width = q - 7 * divisor * (size - 1)
    require(width > 0, ("mode width", q, owner, width))
    require(h == owner % (2 * q), ("half residue", q, owner, h))
    n = (owner - h) // (2 * q)
    centre = (Fraction(n) + Fraction(h, 2 * q)) / owner
    require(centre == Fraction(1, 2 * q), ("centre", q, owner, centre))
    return {
        "divisor": divisor,
        "order": phase_modulus,
        "size": size,
        "width": width,
        "centre": centre,
    }


def symbolic_owner_table() -> dict[str, object]:
    """Check the six nearest-multiple channels without choosing p or x."""

    # owner r=a*p+b; r*x == (a*sigma+b*m)*p+b*y mod 6p
    owner_parameters = {
        "1": (0, 1),
        "p-1": (1, -1),
        "p+1": (1, 1),
        "2p-1": (2, -1),
        "2p+1": (2, 1),
        "3p-1": (3, -1),
    }
    tables = {
        1: ("1", "p-1", "2p-1", "3p-1", "2p+1", "p+1"),
        -1: ("1", "p+1", "2p+1", "3p-1", "2p-1", "p-1"),
    }
    rows = []
    for sigma, names in tables.items():
        for multiple, name in enumerate(names):
            a, b = owner_parameters[name]
            coefficient = (a * sigma + b * multiple) % 6
            require(coefficient == 0, ("owner channel", sigma, multiple, name, coefficient))
            rows.append((sigma, multiple, name, b))

    residue_rows = []
    for residue in range(42):
        p = residue if residue >= 13 else residue + 42
        h = threshold(p)
        require(7 * h < 3 * p <= 7 * (h + 1), ("strict threshold", p, h))
        if p % 7:
            require(3 * p - 6 * (h + 1) <= h, ("generic gap", p, h))
        else:
            s = p // 7
            require(h == 3 * s - 1, ("septimal threshold", p, h, s))
            exceptional_phase_is_admissible = bool(s % 2 and s % 3)
            require(
                exceptional_phase_is_admissible == (p % 42 in (7, 35)),
                ("septimal boundary", p, s),
            )
        residue_rows.append((residue, predicted_cover(p), h))
    return {"channels": tuple(rows), "residues": tuple(residue_rows)}


def direct_boundary_audit() -> dict[str, object]:
    checked_sheets = 0
    reference_sheets = 0
    failures = []
    mode_rows = []
    for p in range(13, 601):
        q = 3 * p
        owners = template(p)
        require(len(set(owners)) == 8 and all(0 < owner < q for owner in owners), ("owners", p))
        masks = tuple(direct_mask(q, owner) for owner in owners)
        require(all(mask for mask in masks), ("nonempty modes", p))
        actual = union_mask(q, owners) == (1 << q) - 1
        require(actual == predicted_cover(p), ("cover boundary", p, actual))
        if not actual:
            failures.append(p)
        checked_sheets += 8 * q
        if p <= 60:
            for owner, mask in zip(owners, masks):
                require(mask == fraction_mask(q, owner), ("fraction mask", p, owner))
                reference_sheets += q

    for p in range(13, 55):
        q = 3 * p
        modes = tuple(mode_at_half(q, owner) for owner in template(p))
        mode_rows.append((p, tuple(mode["width"] for mode in modes)))

    # Sharp septimal controls.
    p = 35
    s = p // 7
    x = 17 * s
    sheet = (x - 1) // 2
    require(x % 2 == 1 and x % 3 != 0, ("hostile phase", x))
    require(
        all(not (direct_mask(3 * p, owner) >> sheet & 1) for owner in template(p)),
        ("p35 hostile", sheet),
    )
    require(union_mask(42, template(14)) == (1 << 42) - 1, "even septimal repair")
    require(union_mask(63, template(21)) == (1 << 63) - 1, "ternary septimal repair")

    expected_failures = tuple(p for p in range(13, 601) if p % 42 in (7, 35))
    require(tuple(failures) == expected_failures, ("failure list", failures))
    return {
        "checked_sheets": checked_sheets,
        "reference_sheets": reference_sheets,
        "failures": tuple(failures),
        "mode_rows": tuple(mode_rows),
        "hostile": (p, x, sheet),
    }


def odd_layer_conjugacy_audit() -> dict[str, int]:
    """Check B_Q(2s)=phi^*Z_Q(s) for every odd Q<=301.

    The canonical half owner is the sign representative of 2s modulo 2Q.
    The sheet permutation is phi(ell)=2ell+1 modulo Q.
    """

    checked_rows = 0
    checked_cells = 0
    for q in range(3, 302, 2):
        for owner in range(1, q):
            half_owner = 2 * owner
            if half_owner > q:
                half_owner = 2 * q - half_owner
            require(0 < half_owner < q, ("half owner", q, owner, half_owner))
            fixed = zero_mask(q, owner)
            half = direct_mask(q, half_owner)
            transported = 0
            for sheet in range(q):
                target = (2 * sheet + 1) % q
                if fixed >> target & 1:
                    transported |= 1 << sheet
                require(
                    ((half >> sheet) & 1) == ((fixed >> target) & 1),
                    ("odd zero/half conjugacy", q, owner, sheet),
                )
                checked_cells += 1
            require(half == transported, ("transported mask", q, owner))
            checked_rows += 1
    require(checked_rows == 22_650, checked_rows)
    require(checked_cells == 4_567_750, checked_cells)
    return {"rows": checked_rows, "cells": checked_cells}


def cap_seven_rank(q: int) -> int:
    if any(q % atom == 0 for atom in RANK4_ATOMS):
        return 4
    if any(q % atom == 0 for atom in RANK5_ATOMS):
        return 5
    if any(q % atom == 0 for atom in RANK6_ATOMS):
        return 6
    if any(q % atom == 0 for atom in RANK7_ATOMS):
        return 7
    return 8


def family_rank(k: int) -> int:
    require(k >= 1, k)
    if k % 3 == 2:
        return 4
    if any(k % modulus == residue for modulus, residue in RANK6_K):
        return 6
    if any(k % modulus == residue for modulus, residue in RANK7_K):
        return 7
    return 8


def crt(residues: tuple[tuple[int, int], ...]) -> tuple[int, int]:
    value = 0
    modulus = 1
    for next_modulus, next_residue in residues:
        require(gcd(modulus, next_modulus) == 1, ("CRT coprimality", modulus, next_modulus))
        inverse = pow(modulus, -1, next_modulus)
        step = ((next_residue - value) * inverse) % next_modulus
        value += modulus * step
        modulus *= next_modulus
        value %= modulus
    return value, modulus


def rank_census() -> dict[str, object]:
    period = 1
    for prime in PERIOD_PRIMES:
        period *= prime
    require(period == EXPECTED_PERIOD, period)

    for k in range(1, 100_001):
        q = 3 * (14 * k - 1)
        require(family_rank(k) == cap_seven_rank(q), ("rank classifier", k, q))

    count4 = period // 3
    rank6_union = 5 * 11 * 23 - 4 * 10 * 22
    count6 = 2 * rank6_union * 13 * 17 * 29
    rank7_union = 13 * 17 * 29 - 12 * 16 * 28
    count7 = 2 * 4 * 10 * 22 * rank7_union
    count8 = 2 * 4 * 10 * 12 * 16 * 22 * 28
    counts = {4: count4, 5: 0, 6: count6, 7: count7, 8: count8}
    require(counts == EXPECTED_COUNTS, counts)
    require(sum(counts.values()) == period, counts)

    trigger = {3: 2, 5: 4, 11: 4, 13: 1, 17: 11, 23: 5, 29: 27}
    minimality = []
    for prime in PERIOD_PRIMES:
        residues = tuple(
            (modulus, trigger[modulus] if modulus == prime else 0)
            for modulus in PERIOD_PRIMES
        )
        base, crt_modulus = crt(residues)
        require(crt_modulus == period and base > 0, ("CRT witness", prime, base))
        shifted = (base + period // prime) % period
        if shifted == 0:
            shifted = period
        require(family_rank(base) != family_rank(shifted), ("minimal period", prime, base, shifted))
        minimality.append((prime, base, family_rank(base), shifted, family_rank(shifted)))

    densities = {rank: Fraction(count, period) for rank, count in counts.items()}
    require(densities[4] == Fraction(1, 3), densities)
    require(densities[6] == Fraction(14, 69), densities)
    require(densities[7] == Fraction(33056, 442221), densities)
    require(densities[8] == Fraction(57344, 147407), densities)
    label_harmonic_rank8 = densities[8] / 42
    require(label_harmonic_rank8 == Fraction(4096, 442221), label_harmonic_rank8)

    return {
        "period": period,
        "counts": counts,
        "densities": densities,
        "minimality": tuple(minimality),
        "label_harmonic_rank8": label_harmonic_rank8,
    }


def u_spine_report() -> dict[str, object]:
    residues = tuple(
        t
        for t in range(21)
        if ((2 * t + 1) ** 2 + 2) % 42 == 39
    )
    require(residues == (5, 8, 12, 15), residues)
    examples = []
    polynomial_lanes = []
    for t in residues:
        q = (2 * t + 1) ** 2 + 2
        p = q // 3
        k = (p + 1) // 14
        require(q == 3 * (14 * k - 1), ("U-spine parametrization", t, q, k))
        require(union_mask(q, template(p)) == (1 << q) - 1, ("U-spine cover", t))
        examples.append((t, k, q, family_rank(k)))
        polynomial_lanes.append((t, 42, 4 * t + 2, k))
    require(examples[:2] == [(5, 3, 123, 8), (8, 7, 291, 8)], examples)

    def annotated_symbol(t: int) -> int:
        if t % 21 not in residues:
            return 0
        q = (2 * t + 1) ** 2 + 2
        k = (q + 3) // 42
        require(q == 42 * k - 3, ("annotated lane", t, q, k))
        rank = cap_seven_rank(q)
        require(rank == family_rank(k), ("annotated rank", t, k, rank))
        return rank

    period = EXPECTED_U_SPINE_PERIOD
    require(period == 3 * 3 * 7 * 11 * 17, period)
    word = tuple(annotated_symbol(t) for t in range(period))
    require(
        all(word[t] == annotated_symbol(t + period) for t in range(period)),
        "U-spine period",
    )
    counts = {
        rank: sum(symbol == rank for symbol in word)
        for rank in (4, 6, 7, 8)
    }
    require(counts == EXPECTED_U_SPINE_COUNTS, counts)
    require(sum(counts.values()) == period * 4 // 21, counts)

    period_witnesses = []
    for prime in (3, 7, 11, 17):
        shift = period // prime
        witness = next(
            t
            for t in range(1, period + 1)
            if annotated_symbol(t) != annotated_symbol(t + shift)
        )
        period_witnesses.append(
            (prime, witness, annotated_symbol(witness), witness + shift, annotated_symbol(witness + shift))
        )

    ambient_densities = {
        rank: Fraction(count, period)
        for rank, count in counts.items()
    }
    conditioned = {
        rank: Fraction(count, sum(counts.values()))
        for rank, count in counts.items()
    }
    require(ambient_densities == {
        4: Fraction(4, 63),
        6: Fraction(16, 693),
        7: Fraction(16, 1309),
        8: Fraction(120, 1309),
    }, ambient_densities)
    require(conditioned == {
        4: Fraction(1, 3),
        6: Fraction(4, 33),
        7: Fraction(12, 187),
        8: Fraction(90, 187),
    }, conditioned)
    return {
        "residues_mod21": residues,
        "density": Fraction(4, 21),
        "examples": tuple(examples),
        "polynomial_lanes": tuple(polynomial_lanes),
        "annotated_period": period,
        "annotated_counts": counts,
        "ambient_densities": ambient_densities,
        "conditioned": conditioned,
        "period_witnesses": tuple(period_witnesses),
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
    allowed = {"__future__", "ast", "fractions", "hashlib", "json", "math", "pathlib", "sys"}
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
    symbolic = symbolic_owner_table()
    boundary = direct_boundary_audit()
    layer_conjugacy = odd_layer_conjugacy_audit()
    ranks = rank_census()
    spine = u_spine_report()

    semantic_payload = {
        "theorem": THEOREM_ID,
        "dependencies": dependency_hashes,
        "symbolic": symbolic,
        "boundary": boundary,
        "odd_layer_conjugacy": layer_conjugacy,
        "ranks": {
            "period": ranks["period"],
            "counts": ranks["counts"],
            "densities": {rank: fraction_text(value) for rank, value in ranks["densities"].items()},
            "minimality": ranks["minimality"],
            "label_harmonic_rank8": fraction_text(ranks["label_harmonic_rank8"]),
        },
        "spine": {
            "residues": spine["residues_mod21"],
            "density": fraction_text(spine["density"]),
            "examples": spine["examples"],
            "polynomial_lanes": spine["polynomial_lanes"],
            "annotated_period": spine["annotated_period"],
            "annotated_counts": spine["annotated_counts"],
            "ambient_densities": {
                rank: fraction_text(value)
                for rank, value in spine["ambient_densities"].items()
            },
            "conditioned": {
                rank: fraction_text(value)
                for rank, value in spine["conditioned"].items()
            },
            "period_witnesses": spine["period_witnesses"],
        },
        "security": security,
    }
    semantic_hash = sha256(
        dumps(semantic_payload, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "UNFROZEN":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, ("semantic drift", semantic_hash))

    print("THM-3469 EXACT DETERMINISTIC COMPANION")
    print("STATUS: PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED")
    print("DEPENDENCIES:")
    for label, relative, digest in dependency_hashes:
        print(f"  {label} {digest} {relative}")
    print(f"SECURITY_IMPORTS: {','.join(security)}")
    print(f"SYMBOLIC_OWNER_CHANNELS: {symbolic['channels']}")
    print("STRICT_THRESHOLD: h=floor((3p-1)/7); generic gap closes; only p=7s tie remains")
    print("COVER_BOUNDARY_P_GE_13: cover iff p mod 42 not in {7,35}")
    print(f"DIRECT_SHEET_CHECKS_P13_TO600: {boundary['checked_sheets']}")
    print(f"FRACTION_REFERENCE_SHEETS_P13_TO60: {boundary['reference_sheets']}")
    print(
        "ODD_ZERO_HALF_CONJUGACY_Q3_TO301: "
        f"rows={layer_conjugacy['rows']} cells={layer_conjugacy['cells']}"
    )
    print(f"FAILURES_P13_TO600: {boundary['failures']}")
    print(f"SHARP_HOSTILE_P_X_SHEET: {boundary['hostile']}")
    print("SEPTIMAL_REPAIRS: p=14 parity removes tie; p=21 ternary backbone absorbs tie")
    print("SPECIAL_FAMILY: p=14k-1, q=3p, owners=(1,p-1,p+1,2p-1,2p,2p+1,3p-6,3p-1)")
    print(f"EXACT_RANK_PERIOD_K: {ranks['period']}")
    print(f"EXACT_RANK_COUNTS: {tuple(sorted(ranks['counts'].items()))}")
    print(
        "EXACT_RANK_DENSITIES: "
        + str(tuple((rank, fraction_text(value)) for rank, value in sorted(ranks["densities"].items())))
    )
    print(f"MINIMAL_PERIOD_CRT_WITNESSES: {ranks['minimality']}")
    print(
        "RANK8_LINEAR_LABEL_HARMONIC_COEFFICIENT: "
        f"{fraction_text(ranks['label_harmonic_rank8'])}"
    )
    print(f"U_SPINE_TEMPLATE_RESIDUES_MOD21: {spine['residues_mod21']}")
    print(f"U_SPINE_TEMPLATE_DENSITY_AND_HARMONIC_COEFFICIENT: {fraction_text(spine['density'])}")
    print(f"U_SPINE_FIRST_RESIDUE_EXAMPLES: {spine['examples']}")
    print(f"U_SPINE_LANE_K_POLYNOMIALS: {spine['polynomial_lanes']}")
    print(f"U_SPINE_ANNOTATED_MINIMAL_PERIOD: {spine['annotated_period']}")
    print(f"U_SPINE_ANNOTATED_COUNTS: {tuple(sorted(spine['annotated_counts'].items()))}")
    print(
        "U_SPINE_AMBIENT_DENSITIES: "
        + str(tuple(
            (rank, fraction_text(value))
            for rank, value in sorted(spine["ambient_densities"].items())
        ))
    )
    print(
        "U_SPINE_CONDITIONED_DENSITIES: "
        + str(tuple(
            (rank, fraction_text(value))
            for rank, value in sorted(spine["conditioned"].items())
        ))
    )
    print(f"U_SPINE_MINIMAL_PERIOD_WITNESSES: {spine['period_witnesses']}")
    print(f"SEMANTIC_SHA256: {semantic_hash}")
    print("VERDICT: universal eight-owner cover boundary and exact periodic rank-4/6/7/8 family; no endpoint current or LRC(14) conclusion")


if __name__ == "__main__":
    main()
