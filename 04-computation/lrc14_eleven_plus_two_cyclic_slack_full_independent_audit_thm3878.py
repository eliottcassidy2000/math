#!/usr/bin/env python3
"""Independent full audit of the conditional t>=U 11+2 ledger.

This checker imports only the separate real-line interval implementation in
lrc14_eleven_plus_two_cyclic_slack_independent_audit_thm3878.py.  It subtracts periodic open danger intervals directly,
rather than sampling wall cells or importing the primary scale-one search.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import importlib.util
import json
from math import floor
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")
HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("interval_audit", HERE / "lrc14_eleven_plus_two_cyclic_slack_independent_audit_thm3878.py")
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load independent interval audit")
IA = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(IA)


Interval = tuple[Fraction, Fraction]


def subtract_open(base: Interval, removed: list[Interval]) -> list[Interval]:
    """Lengths of the components of base minus a union of open intervals."""
    left, right = base
    cursor = left
    pieces: list[Interval] = []
    for rleft, rright in removed:
        if rright <= cursor:
            continue
        if rleft >= right:
            break
        if rleft > cursor:
            pieces.append((cursor, min(rleft, right)))
        cursor = max(cursor, rright)
        if cursor >= right:
            break
    if cursor < right:
        pieces.append((cursor, right))
    return [(a, b) for a, b in pieces if a < b]


def danger_components(p: int, q: int) -> list[Interval]:
    full = IA.merge_open(IA.raw_danger(p) + IA.raw_danger(q))
    return IA.one_period(full, Fraction(1))


def residual_width(p: int, q: int, a: int) -> Fraction:
    base = danger_components(p, q)
    removed = IA.merge_open(IA.raw_danger(a))
    pieces = []
    for component in base:
        pieces.extend(subtract_open(component, removed))
    return max((right - left for left, right in pieces), default=Fraction(0))


def fmt(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def main() -> None:
    pairs = IA.build_atlas()
    IA.require(len(pairs) == 5855, "atlas")
    table = []
    for p, q in pairs:
        beta1, _, _, beta2, _, _ = IA.interval_widths(p, q)
        table.append((p, q, beta1, beta2))

    odd_residual = [(p, q, b2) for p, q, _, b2 in table if p % 2 and q % 2 and b2 > 0]
    scale2_closed = [(p, q) for p, q, b2 in odd_residual if 42 * b2 <= 1]
    scale2_open = [(p, q, b2) for p, q, b2 in odd_residual if 42 * b2 > 1]
    IA.require(len(scale2_closed) == 1649, "scale2 closed count")
    IA.require(scale2_open == [(1, 9, Fraction(2, 63))], "scale2 exception")

    base_closed = [(p, q) for p, q, b1, _ in table if 42 * b1 <= 1]
    base_open = [(p, q, b1) for p, q, b1, _ in table if 42 * b1 > 1]
    IA.require(len(base_closed) == 5445 and len(base_open) == 410, "scale1 base split")

    deletion_closed = []
    deletion_open = []
    for p, q, beta1 in base_open:
        gamma_p = residual_width(p, q, p)
        gamma_q = residual_width(p, q, q)
        if gamma_p < Fraction(1, 91 * p) or gamma_q < Fraction(1, 91 * q):
            deletion_closed.append((p, q))
        else:
            deletion_open.append((p, q, beta1, gamma_p, gamma_q))
    IA.require(len(deletion_closed) == 353 and len(deletion_open) == 57,
               "deletion split")

    final_open = []
    auxiliary_closures = []
    for p, q, beta1, _, _ in deletion_open:
        limit = floor(Fraction(13, 7) / beta1) + 1
        winners = []
        best = None
        for a in range(1, limit + 1):
            gamma = residual_width(p, q, a)
            record = (a * gamma, a, gamma)
            if best is None or record < best:
                best = record
            if gamma < Fraction(1, 91 * a):
                winners.append((a, gamma))
        if best is None:
            raise RuntimeError("empty auxiliary search")
        if winners:
            auxiliary_closures.append((p, q, tuple(winners)))
        else:
            final_open.append((p, q, beta1, limit, best))
    IA.require(not auxiliary_closures and len(final_open) == 57,
               "arbitrary auxiliary stopping boundary")

    combined = tuple((1, p, q) for p, q, *_ in final_open) + ((2, 1, 9),)
    survivor_digest = sha256(json.dumps(combined, separators=(",", ":")).encode("ascii")).hexdigest()
    classification = {
        "base_closed": tuple(base_closed),
        "deletion_closed": tuple(deletion_closed),
        "final_open": tuple(
            (p, q, fmt(42 * beta), limit, fmt(best[0]), best[1], fmt(best[2]))
            for p, q, beta, limit, best in final_open
        ),
        "scale2_open": tuple((p, q, fmt(42 * beta)) for p, q, beta in scale2_open),
    }
    semantic_digest = sha256(json.dumps(classification, sort_keys=True, separators=(",", ":")).encode("ascii")).hexdigest()

    print("LRC14_REMAINING_CYCLIC_FULL_INDEPENDENT_AUDIT_20260823")
    print("method=real_line_periodic_open_unions+direct_subtraction;no_primary_import")
    print(f"scale1_base={len(base_closed)}/5855;deletion_additional={len(deletion_closed)};final={len(final_open)}")
    print(f"scale2_base={len(scale2_closed)}/1650;exception={tuple((p,q,fmt(42*b)) for p,q,b in scale2_open)}")
    print(f"combined_t_ge_U_survivors={len(combined)}")
    print("combined_survivors=" + repr(combined))
    print(f"combined_survivor_sha256={survivor_digest}")
    print(f"classification_semantic_sha256={semantic_digest}")
    print(f"gates={IA.GATES}")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
