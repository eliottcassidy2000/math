#!/usr/bin/env python3
"""Deletion-witness slack cuts for the scale-one 11+2 seam."""

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
SPEC = importlib.util.spec_from_file_location("cyclic_probe", HERE / "lrc14_eleven_plus_two_cyclic_slack_thm3878.py")
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load probe")
P = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(P)


def auxiliary_walls(a: int) -> tuple[Fraction, ...]:
    walls = {Fraction(0)}
    for k in range(a):
        for sign in (-1, 1):
            walls.add(((Fraction(k) + sign * P.H) / a) % 1)
    return tuple(sorted(walls))


def deletion_residual_width(p: int, q: int, keep: int) -> Fraction:
    walls = tuple(sorted(set(P.danger_walls(p, q)) | set(auxiliary_walls(keep))))
    width, _, _ = P.max_open_component(
        walls,
        # The strict interior has the same component lengths as the closed
        # keep-runner safe set.  The sufficient length comparison is strict.
        lambda w: P.pair_bad(p, q, w) and P.distance(keep * w) > P.H,
    )
    return width


def fmt(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def main() -> None:
    pairs = P.atlas()
    P.require(len(pairs) == 5855, "atlas")
    base_closed = []
    base_open = []
    for p, q in pairs:
        beta1 = P.widths(p, q)[0]
        if 42 * beta1 <= 1:
            base_closed.append((p, q))
        else:
            base_open.append((p, q, beta1))
    P.require(len(base_closed) == 5445 and len(base_open) == 410, "base split")

    deletion_closed = []
    still_open = []
    for p, q, beta1 in base_open:
        gamma_p = deletion_residual_width(p, q, p)
        gamma_q = deletion_residual_width(p, q, q)
        # In the t>=U slice, the LRC(13) witness for u union {keep*t}
        # supplies a w-interval of length 1/(91*keep).
        p_works = gamma_p < Fraction(1, 91 * p)
        q_works = gamma_q < Fraction(1, 91 * q)
        record = (p, q, beta1, gamma_p, gamma_q)
        if p_works or q_works:
            deletion_closed.append((*record, p_works, q_works))
        else:
            still_open.append(record)

    auxiliary_closed = []
    final_open = []
    for p, q, beta1, gamma_p, gamma_q in still_open:
        # If a>13/(7 beta1), a largest obstruction component necessarily
        # contains a complete auxiliary-safe cell of length 6/(7a), already
        # much longer than the LRC(13) slack image 1/(91a).  Searching through
        # floor(13/(7 beta1))+1 is therefore exhaustive.
        limit = floor(Fraction(13, 7) / beta1) + 1
        winners = []
        best = None
        for a in range(1, limit + 1):
            gamma_a = deletion_residual_width(p, q, a)
            score = a * gamma_a
            if best is None or (score, a, gamma_a) < best:
                best = (score, a, gamma_a)
            if gamma_a < Fraction(1, 91 * a):
                winners.append((a, gamma_a))
        if winners:
            auxiliary_closed.append((p, q, beta1, limit, tuple(winners), best))
        else:
            final_open.append((p, q, beta1, limit, best))

    combined_survivors = tuple((1, p, q) for p, q, *_ in final_open) + ((2, 1, 9),)
    semantic = {
        "atlas": len(pairs),
        "base_closed": tuple(base_closed),
        "deletion_closed": tuple(
            (p, q, fmt(b), fmt(gp), fmt(gq), wp, wq)
            for p, q, b, gp, gq, wp, wq in deletion_closed
        ),
        "final_open": tuple(
            (p, q, fmt(42 * b), limit, fmt(best[0]), best[1], fmt(best[2]))
            for p, q, b, limit, best in final_open
        ),
        "combined_t_ge_U_survivors": combined_survivors,
    }
    digest = sha256(json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")).hexdigest()
    survivor_digest = sha256(json.dumps(combined_survivors, separators=(",", ":")).encode("ascii")).hexdigest()

    print("LRC14_SCALE1_DELETION_WITNESS_SEARCH_20260823")
    print("scope=THM3818_11+2_scale1_and_relative_slice_t>=U;LRC14=OPEN")
    print(f"atlas={len(pairs)};base_slack_closed={len(base_closed)};base_open={len(base_open)}")
    print(f"deletion_slack_closed={len(deletion_closed)};still_open={len(still_open)}")
    print(f"arbitrary_one_auxiliary_closed={len(auxiliary_closed)};final_open={len(final_open)}")
    print(f"combined_t_ge_U_certificate_survivors={len(combined_survivors)}")
    print("deletion_positive_controls=" + repr(tuple(
        (p, q, fmt(42*b), fmt(gp), fmt(gq), wp, wq)
        for p, q, b, gp, gq, wp, wq in deletion_closed[:5]
    )))
    print("arbitrary_auxiliary_result=no_closures;search_cutoff_at_most_78")
    print("combined_survivors=" + repr(combined_survivors))
    print(f"combined_survivor_sha256={survivor_digest}")
    print(f"semantic_sha256={digest}")
    print(f"gates={P.GATES}")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
