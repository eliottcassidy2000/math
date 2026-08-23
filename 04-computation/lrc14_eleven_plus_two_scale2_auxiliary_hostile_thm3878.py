#!/usr/bin/env python3
"""Search one-extra-runner LRC(13) cuts of the residual pair obstruction."""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import importlib.util
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


def gamma(p: int, q: int, a: int) -> tuple[Fraction, int, Fraction]:
    walls = tuple(sorted(set(P.doubled_two_lift_walls(p, q)) | set(auxiliary_walls(a))))
    return P.max_open_component(
        walls,
        # Strict interior has the same component lengths as the actual closed
        # auxiliary-safe set.  The final sufficient comparison is strict.
        lambda w: P.two_lift_bad(p, q, w / 2) and P.distance(a * w) > P.H,
    )


def fmt(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def main() -> None:
    print("AUXILIARY_RUNNER_SEARCH_20260823")
    # The sole t>=U exception in the main theorem is (1,9).  Its quotient
    # obstruction consists of two arcs of length 2/63.  For a>=59, every
    # such arc contains a complete auxiliary-safe cell: the worst distance
    # from an arbitrary start to the end of the next complete cell is
    # 13/(7a)<2/63.  That cell has length 6/(7a)>1/(91a), so no a>=59 can
    # satisfy the desired strict cut.  It is therefore exhaustive to test
    # a=1,...,58.
    pair = (1, 9)
    winners = []
    rows = []
    for a in range(1, 59):
        width, components, measure = gamma(*pair, a)
        rows.append((a * width, a, width, components, measure))
        if width < Fraction(1, 91 * a):
            winners.append((a, width, components, measure))
    best = min(rows)
    score, a, width, components, measure = best
    if winners:
        raise RuntimeError(("unexpected auxiliary winner", winners))
    print(
        f"pair={pair};tested_a=1..58;winners=();best_a={a};"
        f"best_gamma={fmt(width)};best_a_gamma={fmt(score)};"
        f"components={components};measure={fmt(measure)}"
    )
    print("a_ge_59_stopping=each_2/63_obstruction_arc_contains_full_6/(7a)_auxiliary_safe_cell")
    print("consequence=no_single_auxiliary_multiplier_closes_(1,9)_via_LRC13_slack_in_t>=U_slice")
    print(f"gates={P.GATES}")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
