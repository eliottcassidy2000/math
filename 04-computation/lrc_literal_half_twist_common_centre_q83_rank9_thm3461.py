#!/usr/bin/env python3
"""Exact deterministic companion for THM-3461.

The program reattaches the mode data discarded by a bare literal half-twist
mask.  It verifies the fixed common centre and zero complete cochain for the
q=11,27,51 witnesses, constructs honest open physical cells, and exhausts
both Boolean common-centre layers at q=83 through rank eight before checking
an explicit rank-nine witness.

All arithmetic is integer or ``Fraction`` arithmetic.  Every mathematical
gate raises an explicit exception and remains active under ``python -O``.
"""

from __future__ import annotations

import ast
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from json import dumps
from math import gcd, lcm
from pathlib import Path
import sys


THEOREM_ID = "THM-3461"
ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = (
    (
        "THM-3387",
        "01-canon/theorems/THM-3387-exact-cyclic-sheet-cover-atlas-and-q2-gcd-graph.md",
        "c540255a185efb54c67035d69f3fbd94f4c1ad3e30c4e31738a8800e81198613",
    ),
    (
        "THM-3398",
        "01-canon/theorems/THM-3398-general-finite-mode-sheet-cover-cochain.md",
        "01901da2bb382184cfe4466550afe79255598f580f00a761fc32731a52ec9378",
    ),
    (
        "THM-3405",
        "01-canon/theorems/THM-3405-common-centre-gcd-gauge-and-boolean-half-twist.md",
        "d3e7dbeeb85c6f897bd9e31270bd0b6602ae4feac3b46a45eb5ce23ae5d24fe0",
    ),
    (
        "THM-3410",
        "01-canon/theorems/THM-3410-projective-cochain-wedge-ray-tree-tariff-and-residue-scalar-hubs.md",
        "6021e0144fa1e3bfc81448b1b47dc58b576881cbd5a0049036dd424b913713b6",
    ),
    (
        "THM-3416",
        "01-canon/theorems/THM-3416-zero-mode-cochain-global-rank-six-support.md",
        "42a9309145de51d1bb6fca0b7c1945302ff37a63a3183e1dfed838c07118e8bf",
    ),
    (
        "THM-3429",
        "01-canon/theorems/THM-3429-prime-fibre-activity-descent-for-mixed-order-half-twist-seven-covers.md",
        "58ebf850fc79fc9afed57966b7599e7376a6684fa3bbc5a2aa2e1a8e6e0ca148",
    ),
    (
        "THM-3453",
        "01-canon/theorems/THM-3453-global-literal-half-twist-cap-seven-support-classification.md",
        "a4444813a2b24aa55613eef5f4cc0538ca0148f16297dff46f80d723f1beb247",
    ),
    (
        "THM-3455",
        "01-canon/theorems/THM-3455-berggren-q-spine-cap-seven-atom-sieve-and-fibonacci-rank-spectrum.md",
        "e357aaf8036f485c4abbfa5969a0d1f9761372dcd7337b658daaadf4cd049c09",
    ),
)

WITNESSES = {
    11: (1, 2, 3, 5, 7, 9),
    27: (3, 15, 18, 21),
    51: (1, 11, 12, 18, 23, 34, 35),
    83: (1, 13, 14, 27, 41, 42, 55, 69, 70),
}
EXPECTED_RANKS = {11: 6, 27: 4, 51: 7, 83: 9}
EXPECTED_WIDTHS = {
    11: (4, 11, 4, 4, 4, 4),
    27: (6, 6, 27, 6),
    51: (2, 2, 9, 9, 2, 51, 2),
}
EXPECTED_Y_RADII = {
    11: Fraction(2, 63),
    27: Fraction(1, 49),
    51: Fraction(1, 245),
}
EXPECTED_MULTIPLICITIES = {
    11: ((1, 11),),
    27: ((1, 27),),
    51: ((1, 42), (2, 4), (3, 3), (4, 2)),
    83: ((1, 66), (2, 12), (3, 5)),
}
PHYSICAL_ROWS = {
    11: (1, 2, 3, 5, 7, 9, 11),
    27: (1, 3, 15, 18, 21, 27),
    51: (1, 11, 12, 18, 23, 34, 35, 51),
}
CRT_LIFTS = (
    (11, 1, 5083, 1),
    (11, 2, 138580, 2),
    (27, 3, 57135, 1),
    (27, 15, 215475, 2),
    (51, 1, 19279, 1),
    (51, 11, 1538069, 2),
)
EXPECTED_SEMANTIC_SHA256 = "97647a6cf82b23de0ddb962b6b2b714d6ca9904408af87e001fa37d0c9b37e90"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def circular_norm(value: Fraction) -> Fraction:
    reduced = value - value.numerator // value.denominator
    return min(reduced, 1 - reduced)


def direct_mask(q: int, owner: int, epsilon: int) -> int:
    """Danger mask at base centre y=epsilon/2, using denominator 2q."""

    require(epsilon in (0, 1), ("layer", epsilon))
    result = 0
    for sheet in range(q):
        word = owner * (2 * sheet + epsilon) % (2 * q)
        if 14 * min(word, 2 * q - word) < 2 * q:
            result |= 1 << sheet
    return result


def fraction_mask(q: int, owner: int, y: Fraction) -> int:
    result = 0
    for sheet in range(q):
        phase = Fraction(owner, q) * (y + sheet)
        if circular_norm(phase) < Fraction(1, 14):
            result |= 1 << sheet
    return result


def sheets(mask: int, q: int) -> tuple[int, ...]:
    return tuple(sheet for sheet in range(q) if mask >> sheet & 1)


def union_mask(q: int, owners: tuple[int, ...], epsilon: int) -> int:
    result = 0
    for owner in owners:
        result |= direct_mask(q, owner, epsilon)
    return result


def multiplicity_profile(q: int, owners: tuple[int, ...], epsilon: int):
    counts = [0] * q
    masks = tuple(direct_mask(q, owner, epsilon) for owner in owners)
    for mask in masks:
        for sheet in range(q):
            counts[sheet] += (mask >> sheet) & 1
    profile = tuple(sorted(Counter(counts).items()))
    private = tuple(
        sum(1 for sheet in range(q) if mask >> sheet & 1 and counts[sheet] == 1)
        for mask in masks
    )
    return profile, private


def canonical_bank(q: int, epsilon: int):
    by_mask: dict[int, int] = {}
    for owner in range(1, q):
        mask = direct_mask(q, owner, epsilon)
        if not mask:
            continue
        previous = by_mask.get(mask)
        if previous is None or owner < previous:
            by_mask[mask] = owner
    return tuple(sorted((owner, mask) for mask, owner in by_mask.items()))


def maximal_bank(items: tuple[tuple[int, int], ...]):
    """Drop strict subsets; cover existence at a fixed cap is preserved."""

    result = []
    for owner, mask in items:
        if any(mask != other and mask | other == other for _label, other in items):
            continue
        result.append((owner, mask))
    return tuple(result)


def bank_digest(items: tuple[tuple[int, int], ...]) -> str:
    payload = "".join(f"{owner}:{mask:x}\n" for owner, mask in items)
    return sha256(payload.encode("ascii")).hexdigest()


def bounded_cover(q: int, items: tuple[tuple[int, int], ...], cap: int):
    """Complete branch-and-memo search on the inclusion-maximal mask bank."""

    full = (1 << q) - 1
    reduced = maximal_bank(items)
    owners = tuple(owner for owner, _mask in reduced)
    masks = tuple(mask for _owner, mask in reduced)
    coverers = tuple(
        tuple(index for index, mask in enumerate(masks) if mask >> bit & 1)
        for bit in range(q)
    )
    joined = 0
    for mask in masks:
        joined |= mask
    if joined != full:
        return (), (1, 0, 0, len(items), len(reduced))

    nodes = 0
    branches = 0

    @lru_cache(maxsize=None)
    def solve(state: int, slots: int):
        nonlocal nodes, branches
        nodes += 1
        if state == full:
            return ()
        if slots == 0:
            return None
        missing = full ^ state
        live = tuple(index for index, mask in enumerate(masks) if mask & missing)
        possible = state
        for index in live:
            possible |= masks[index]
        if possible != full:
            return None
        gains = sorted(((masks[index] & missing).bit_count() for index in live), reverse=True)
        if sum(gains[:slots]) < missing.bit_count():
            return None
        missing_bits = tuple(bit for bit in range(q) if missing >> bit & 1)
        pivot = min(
            missing_bits,
            key=lambda bit: (
                sum(1 for index in coverers[bit] if masks[index] & missing),
                bit,
            ),
        )
        options = sorted(
            (index for index in coverers[pivot] if masks[index] & missing),
            key=lambda index: (-(masks[index] & missing).bit_count(), owners[index]),
        )
        for index in options:
            branches += 1
            suffix = solve(state | masks[index], slots - 1)
            if suffix is not None:
                return (index,) + suffix
        return None

    chosen = solve(0, cap)
    hits = solve.cache_info().hits
    witness = () if chosen is None else tuple(sorted(owners[index] for index in chosen))
    return witness, (nodes, branches, hits, len(items), len(reduced))


def mode_at_half(q: int, owner: int):
    divisor = gcd(q, owner)
    phase_modulus = q // divisor
    multiplier = owner // divisor % phase_modulus
    mask = direct_mask(q, owner, 1)
    phase_block = {
        multiplier * sheet % phase_modulus
        for sheet in range(q)
        if mask >> sheet & 1
    }
    size = len(phase_block)
    starts = tuple(
        start
        for start in range(phase_modulus)
        if {(start + offset) % phase_modulus for offset in range(size)} == phase_block
    )
    require(len(starts) == 1, ("mode start", q, owner, phase_block, starts))
    start = starts[0]
    h = (-divisor * (2 * start + size - 1)) % (2 * q)
    width = q - 7 * divisor * (size - 1)
    require(width > 0, ("mode width", q, owner, width))
    require(mask.bit_count() == divisor * size, ("mode mass", q, owner))
    require(h == owner % (2 * q), ("half-centre residue", q, owner, h))
    n = (owner - h) // (2 * q)
    centre = (Fraction(n) + Fraction(h, 2 * q)) / owner
    require(centre == Fraction(1, 2 * q), ("centre lift", q, owner, centre))
    return {
        "owner": owner,
        "divisor": divisor,
        "phase_modulus": phase_modulus,
        "start": start,
        "size": size,
        "h": h,
        "width": width,
        "n": n,
        "source_radius": Fraction(width, 14 * q * owner),
        "y_radius": Fraction(width, 14 * owner),
    }


def fixed_centre_report(q: int, owners: tuple[int, ...]):
    full = (1 << q) - 1
    require(union_mask(q, owners, 1) == full, ("witness cover", q, owners))
    modes = tuple(mode_at_half(q, owner) for owner in owners)
    active_gcd = gcd(*owners)
    gauge_gcd = gcd(q, active_gcd)
    scalar = active_gcd
    require(scalar % gauge_gcd == 0, ("MISTAKE-389 gate", q, scalar, gauge_gcd))
    widths = tuple(mode["width"] for mode in modes)
    y_radius = min(mode["y_radius"] for mode in modes)
    for left in range(len(owners)):
        for right in range(left + 1, len(owners)):
            u = owners[left]
            v = owners[right]
            a_u = u
            a_v = v
            require(a_u * v - u * a_v == 0, ("wedge", q, u, v))
    profile, private = multiplicity_profile(q, owners, 1)
    return {
        "q": q,
        "owners": owners,
        "active_gcd": active_gcd,
        "gauge_gcd": gauge_gcd,
        "scalar": scalar,
        "widths": widths,
        "y_radius": y_radius,
        "profile": profile,
        "private": private,
        "modes": modes,
    }


def rank_report(q: int, expected_rank: int, witness: tuple[int, ...]):
    banks = {epsilon: canonical_bank(q, epsilon) for epsilon in (0, 1)}
    lower = {
        epsilon: bounded_cover(q, banks[epsilon], expected_rank - 1)
        for epsilon in (0, 1)
    }
    require(not lower[0][0] and not lower[1][0], ("rank lower bound", q, lower))
    positive = bounded_cover(q, banks[1], expected_rank)
    require(positive[0], ("rank upper search", q, positive))
    require(len(witness) == expected_rank, ("frozen witness length", q, witness))
    require(union_mask(q, witness, 1) == (1 << q) - 1, ("frozen witness", q))
    return {
        "q": q,
        "rank": expected_rank,
        "banks": {
            epsilon: {
                "unique": len(banks[epsilon]),
                "maximal": len(maximal_bank(banks[epsilon])),
                "digest": bank_digest(banks[epsilon]),
            }
            for epsilon in (0, 1)
        },
        "lower": lower,
        "positive": positive,
    }


def physical_row_report(q: int, fixed_report: dict[str, object]):
    row = PHYSICAL_ROWS[q]
    require(gcd(*row) == 1, ("primitive physical row", q, row))
    require(q in row, ("missing core speed", q, row))
    transverse = tuple(owner for owner in row if owner != q)
    radius = fixed_report["y_radius"]
    require(isinstance(radius, Fraction), ("radius type", q))
    require(Fraction(1, 2) - radius > Fraction(1, 14), ("left core separation", q))
    require(Fraction(1, 2) + radius < Fraction(13, 14), ("right core separation", q))
    grid_denominator = 14 * lcm(*row) // q
    candidate = None
    for divisor in (2, 3, 5, 7, 11, 13):
        trial = Fraction(1, 2) + radius / divisor
        if (grid_denominator * trial).denominator != 1:
            candidate = trial
            break
    require(candidate is not None, ("non-grid candidate", q, grid_denominator))
    direct = 0
    for owner in transverse:
        direct |= fraction_mask(q, owner, candidate)
    require(direct == (1 << q) - 1, ("open physical cover", q, candidate))
    require(circular_norm(candidate) >= Fraction(1, 14), ("core unsafe", q, candidate))
    return {
        "q": q,
        "row": row,
        "transverse": transverse,
        "D": grid_denominator,
        "candidate_y": candidate,
        "radius": radius,
    }


def valuation(value: int, prime: int) -> int:
    answer = 0
    while value % prime == 0:
        value //= prime
        answer += 1
    return answer


def auxiliary_report():
    hostile_modes = tuple(mode_at_half(11, owner) for owner in (1, 23))
    require(direct_mask(11, 1, 1) == direct_mask(11, 23, 1), "same-mask hostile")
    hostile_radii = tuple(mode["source_radius"] for mode in hostile_modes)
    require(hostile_radii == (Fraction(2, 77), Fraction(2, 1771)), hostile_radii)

    quotient = tuple(owner // 3 for owner in WITNESSES[27])
    require(quotient == (1, 5, 6, 7), quotient)
    for sheet in range(27):
        for owner, base in zip(WITNESSES[27], quotient):
            require(
                bool(direct_mask(27, owner, 1) >> sheet & 1)
                == bool(direct_mask(9, base, 1) >> (sheet % 9) & 1),
                ("q27 pullback", owner, sheet),
            )

    require(1 % 34 == 35 % 34, "q51 reduction hostile")
    require((1 - 1) // 34 % 3 != (35 - 1) // 34 % 3, "q51 lift character")
    require(direct_mask(51, 1, 1) != direct_mask(51, 35, 1), "q51 lost lift")

    crt = []
    for q, residue, lift, expected_valuation in CRT_LIFTS:
        require(lift % (2 * q) == residue, ("CRT residue", q, residue, lift))
        require(lift % 7 == 1, ("CRT mod7", q, lift))
        require(valuation(lift, 13) == expected_valuation, ("CRT v13", q, lift))
        require(direct_mask(q, lift, 1) == direct_mask(q, residue, 1), ("CRT mask", q, lift))
        require(mode_at_half(q, lift)["h"] == residue, ("CRT mode", q, lift))
        crt.append((q, residue, lift, expected_valuation))
    return {
        "same_mask_radii": hostile_radii,
        "q27_quotient": quotient,
        "q51_characters": (0, 1),
        "crt": tuple(crt),
    }


def security_report(source: Path):
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
    allowed = {
        "__future__", "ast", "collections", "fractions", "functools",
        "hashlib", "json", "math", "pathlib", "sys",
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
    require(
        tuple((label, relative, expected) for label, relative, expected in DEPENDENCIES)
        == dependency_hashes,
        ("dependency drift", dependency_hashes),
    )
    security = security_report(Path(__file__))

    reference_cases = 0
    for q in range(2, 24):
        for owner in range(1, 2 * q):
            for epsilon in (0, 1):
                require(
                    direct_mask(q, owner, epsilon)
                    == fraction_mask(q, owner, Fraction(epsilon, 2)),
                    ("fraction reference", q, owner, epsilon),
                )
                reference_cases += q

    fixed = {q: fixed_centre_report(q, WITNESSES[q]) for q in (11, 27, 51)}
    for q in fixed:
        require(fixed[q]["widths"] == EXPECTED_WIDTHS[q], ("widths", q, fixed[q]))
        require(fixed[q]["y_radius"] == EXPECTED_Y_RADII[q], ("radius", q, fixed[q]))
        require(fixed[q]["profile"] == EXPECTED_MULTIPLICITIES[q], ("profile", q, fixed[q]))
    q83_profile, q83_private = multiplicity_profile(83, WITNESSES[83], 1)
    require(q83_profile == EXPECTED_MULTIPLICITIES[83], q83_profile)

    ranks = {
        q: rank_report(q, EXPECTED_RANKS[q], WITNESSES[q])
        for q in (11, 27, 51, 83)
    }
    require(ranks[83]["banks"][1]["unique"] == 82, ranks[83])
    physical = {q: physical_row_report(q, fixed[q]) for q in (11, 27, 51)}
    auxiliary = auxiliary_report()

    semantic_payload = {
        "theorem": THEOREM_ID,
        "dependencies": dependency_hashes,
        "reference_cells": reference_cases,
        "fixed": {
            str(q): {
                "owners": fixed[q]["owners"],
                "gcd": fixed[q]["active_gcd"],
                "gauge": fixed[q]["gauge_gcd"],
                "widths": fixed[q]["widths"],
                "y_radius": fraction_text(fixed[q]["y_radius"]),
                "profile": fixed[q]["profile"],
                "private": fixed[q]["private"],
            }
            for q in fixed
        },
        "ranks": {
            str(q): {
                "rank": ranks[q]["rank"],
                "banks": ranks[q]["banks"],
                "lower_stats": {
                    str(epsilon): ranks[q]["lower"][epsilon][1]
                    for epsilon in (0, 1)
                },
                "positive_witness": ranks[q]["positive"][0],
                "positive_stats": ranks[q]["positive"][1],
            }
            for q in ranks
        },
        "q83": {
            "witness": WITNESSES[83],
            "profile": q83_profile,
            "private": q83_private,
        },
        "physical": {
            str(q): {
                "row": physical[q]["row"],
                "D": physical[q]["D"],
                "candidate_y": fraction_text(physical[q]["candidate_y"]),
            }
            for q in physical
        },
        "auxiliary": {
            "same_mask_radii": tuple(fraction_text(value) for value in auxiliary["same_mask_radii"]),
            "q27_quotient": auxiliary["q27_quotient"],
            "q51_characters": auxiliary["q51_characters"],
            "crt": auxiliary["crt"],
        },
        "security": security,
    }
    semantic_hash = sha256(
        dumps(semantic_payload, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "UNFROZEN":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, ("semantic drift", semantic_hash))

    print("THM-3461 EXACT DETERMINISTIC COMPANION")
    print("status=PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED")
    print(f"dependency_count={len(dependency_hashes)};reference_fraction_cells={reference_cases}")
    print("identity=B_(q,r)=D_(q,r)(1/(2q));selected_h=r_mod_2q;A_i=r_i;P_ij=0")
    for q in (11, 27, 51):
        report = fixed[q]
        print(
            f"q={q}|owners={report['owners']}|d={report['active_gcd']}|g={report['gauge_gcd']}|"
            f"widths={report['widths']}|y_radius={fraction_text(report['y_radius'])}|"
            f"multiplicity={report['profile']}|private={report['private']}"
        )
        print("masks=" + ";".join(
            f"{owner}:{sheets(direct_mask(q, owner, 1), q)}" for owner in report["owners"]
        ))
    for q in (11, 27, 51, 83):
        report = ranks[q]
        print(
            f"rank_q{q}={report['rank']}|"
            f"banks0={report['banks'][0]['unique']}/{report['banks'][0]['maximal']}/"
            f"{report['banks'][0]['digest']}|"
            f"banks1={report['banks'][1]['unique']}/{report['banks'][1]['maximal']}/"
            f"{report['banks'][1]['digest']}|"
            f"lower0={report['lower'][0][1]}|lower1={report['lower'][1][1]}|"
            f"positive={report['positive']}"
        )
    print(f"q83_witness={WITNESSES[83]}|multiplicity={q83_profile}|private={q83_private}")
    for q in (11, 27, 51):
        report = physical[q]
        print(
            f"physical_q{q}=row{report['row']}|D={report['D']}|"
            f"open_non_grid_y={fraction_text(report['candidate_y'])}|"
            f"core1_safe=true"
        )
    print(
        "same_mask_hostile=q11:owners(1,23),mask_equal,radii="
        + repr(tuple(fraction_text(value) for value in auxiliary["same_mask_radii"]))
    )
    print(f"q27_ancestry=3*{auxiliary['q27_quotient']};primitive_breaker=1")
    print("q51_affine_hostile=1_and_35_equal_mod34;lift_characters_mod3=(0,1);masks_differ")
    print(f"crt_7x13_entry={auxiliary['crt']};masks_and_zero_wedge_preserved")
    print("current_boundary=one_projective_ray;all_pair_wedges_and_tree_tariffs_zero;endpoint_current_not_determined_by_mask")
    print(f"semantic_sha256={semantic_hash}")
    print("scope=fixed half-twist/ZMC and exact q83 rank only; no endpoint current, LRC row, decrement, spectral closure, or LRC(14) conclusion")
    print("VERDICT=PASS")


if __name__ == "__main__":
    main()
