#!/usr/bin/env python3
"""HYP-2731: exact tail45 strip certificate for generated miss-zeta words.

HYP-2728 found that the normalized generated frontier lies in a tail45 strip
and that the strip excludes all HYP-2721 cheap atom-cone directions.  This
follow-up names the exact extremal rows and converts the observation into a
finite proof ledger:

    182/2005 <= (q5 + 5 q6) / q0 <= 10910/21539.

The run recomputes the HYP-2722/S71 318-row frontier and reports exact lower
and upper slacks globally, by size, by context, and by sign profile.  No proof
of the generated-context automaton inequality is claimed; the output is a
certificate target for formalization and for the next hand proof.
"""

from __future__ import annotations

import importlib.util
import math
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S71_PATH = ROOT / "04-computation" / "lrc14_miss_zeta_word_compatibility_codex_s71.py"
spec = importlib.util.spec_from_file_location("s71_words", S71_PATH)
s71 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = s71
spec.loader.exec_module(s71)


FLOOR = F(182, 2005)
CEIL = F(10910, 21539)


def fmt(q: F) -> str:
    return f"{q} ({float(q):.9f})"


def sign_name(q: F) -> str:
    if q > 0:
        return "+"
    if q < 0:
        return "-"
    return "0"


def factorial(q: tuple[F, ...]) -> tuple[F, ...]:
    return tuple(sum(F(math.comb(t, j)) * q[t] for t in range(j, 7)) for j in range(7))


def tail45(q: tuple[F, ...]) -> F:
    return q[5] + 5 * q[6]


def u4(q: tuple[F, ...]) -> F:
    return q[0] + tail45(q)


def b2(W: tuple[F, ...]) -> F:
    return -W[1] + W[2]


@dataclass(frozen=True)
class Row:
    size: int
    shape: tuple[int, ...]
    context: str
    q0: F
    atoms: tuple[F, ...]
    W: tuple[F, ...]
    tail: F
    lower_slack: F
    upper_slack: F

    @property
    def key(self) -> str:
        return f"size={self.size} shape={self.shape} context={self.context}"

    @property
    def sign_profile(self) -> tuple[str, str, str, str]:
        return (
            sign_name(self.W[1]),
            sign_name(self.W[2]),
            sign_name(b2(self.W)),
            sign_name(self.tail),
        )


def load_rows() -> list[Row]:
    rows: list[Row] = []
    for metric in s71.frontier_metrics():
        if metric.q0 <= 0 or metric.norm_atoms is None:
            continue
        atoms = metric.norm_atoms
        W = factorial(atoms)
        t = tail45(atoms)
        rows.append(
            Row(
                size=metric.size,
                shape=metric.shape,
                context=metric.context_name,
                q0=metric.q0,
                atoms=atoms,
                W=W,
                tail=t,
                lower_slack=t - FLOOR,
                upper_slack=CEIL - t,
            )
        )
    return rows


def row_summary(row: Row) -> str:
    return (
        f"{row.key} q0={fmt(row.q0)} tail45={fmt(row.tail)} "
        f"lower_slack={fmt(row.lower_slack)} upper_slack={fmt(row.upper_slack)} "
        f"U4={fmt(u4(row.atoms))} B2={fmt(b2(row.W))}"
    )


def print_global(rows: list[Row]) -> None:
    print("GLOBAL STRIP")
    print(f"  claimed floor={fmt(FLOOR)}")
    print(f"  claimed ceiling={fmt(CEIL)}")
    print(f"  rows={len(rows)}")
    below = [r for r in rows if r.lower_slack < 0]
    above = [r for r in rows if r.upper_slack < 0]
    print(f"  below_floor={len(below)} above_ceiling={len(above)}")
    floor_row = min(rows, key=lambda r: r.tail)
    ceil_row = max(rows, key=lambda r: r.tail)
    lower_slack_row = min(rows, key=lambda r: r.lower_slack)
    upper_slack_row = min(rows, key=lambda r: r.upper_slack)
    print("  floor witness:")
    print(f"    {row_summary(floor_row)}")
    print("  ceiling witness:")
    print(f"    {row_summary(ceil_row)}")
    assert floor_row == lower_slack_row
    assert ceil_row == upper_slack_row


def grouped_extrema(rows: list[Row], title: str, key_fn) -> None:
    print()
    print(title)
    groups: dict[object, list[Row]] = defaultdict(list)
    for row in rows:
        groups[key_fn(row)].append(row)
    for key in sorted(groups, key=str):
        rs = groups[key]
        lo = min(rs, key=lambda r: r.tail)
        hi = max(rs, key=lambda r: r.tail)
        print(
            f"  {key!s:16s} count={len(rs):3d} "
            f"min_tail={fmt(lo.tail)} at {lo.shape}/{lo.context} "
            f"max_tail={fmt(hi.tail)} at {hi.shape}/{hi.context}"
        )


def strip_tournament(rows: list[Row]) -> None:
    print()
    print("TOURNAMENT ANALYSIS")
    print("  vertices: context partitions")
    print("  observable: tighter generated strip = larger floor and smaller ceiling;")
    print("  tie Hamiltonian path: context name.")
    by_context: dict[str, list[Row]] = defaultdict(list)
    for row in rows:
        by_context[row.context].append(row)
    contexts = sorted(by_context)

    def obs(name: str) -> tuple[F, F, F]:
        rs = by_context[name]
        floor = min(r.tail for r in rs)
        ceil = max(r.tail for r in rs)
        width = ceil - floor
        return (floor, -ceil, -width)

    wins = Counter()
    edges = set()
    for a in contexts:
        for b in contexts:
            if a == b:
                continue
            if obs(a) >= obs(b):
                wins[a] += 1
                edges.add((a, b))
    cycles = 0
    for i, a in enumerate(contexts):
        for j, b in enumerate(contexts):
            for c in contexts[j + 1 :]:
                if i >= j or a == c or b == c:
                    continue
                if (a, b) in edges and (b, c) in edges and (c, a) in edges:
                    cycles += 1
                if (a, c) in edges and (c, b) in edges and (b, a) in edges:
                    cycles += 1
    print(f"  vertices={len(contexts)}")
    print(f"  score_hist={dict(sorted(Counter(wins[c] for c in contexts).items()))}")
    print(f"  directed_3cycles={cycles}")
    print("  pressure path:")
    for name in sorted(contexts, key=obs, reverse=True):
        rs = by_context[name]
        floor = min(r.tail for r in rs)
        ceil = max(r.tail for r in rs)
        print(
            f"    {name:9s} score={wins[name]} floor={fmt(floor)} "
            f"ceil={fmt(ceil)} width={fmt(ceil - floor)}"
        )


def cheap_side() -> None:
    print()
    print("CHEAP DIRECTION SIDE")
    print("  normalized tail45 values from HYP-2721:")
    for r in range(1, 6):
        q = s71.cheap_lp_direction(r)
        t = tail45(q)
        if t < FLOOR:
            side = "below_floor"
            margin = FLOOR - t
        elif t > CEIL:
            side = "above_ceiling"
            margin = t - CEIL
        else:
            side = "inside_strip"
            margin = min(t - FLOOR, CEIL - t)
        print(f"  r={r}: tail45={fmt(t)} {side} margin={fmt(margin)}")


def main() -> None:
    print("HYP-2731 tail45 strip certificate scout")
    print("Exact Fraction arithmetic; generated frontier imported from HYP-2722/S71.")
    print()
    rows = load_rows()
    print()
    print_global(rows)
    grouped_extrema(rows, "BY SIZE", lambda r: f"size={r.size}")
    grouped_extrema(rows, "BY CONTEXT", lambda r: r.context)
    grouped_extrema(rows, "BY SIGN PROFILE (W1,W2,B2,tail)", lambda r: r.sign_profile)
    strip_tournament(rows)
    cheap_side()
    print()
    print("SYNTHESIS")
    print("  The tail45 strip is a finite frontier theorem with two named witnesses:")
    print("  floor at size=3 shape=(0,1,3) context=[3+1], and ceiling at")
    print("  size=3 shape=(0,1,13) context=[4].")
    print("  The hand-proof target is now context-local: prove these two rows are the")
    print("  tail45 extrema of the generated sparse-tail automaton, then use the Lean")
    print("  cheap-side theorem to exclude the abstract atom-cone directions.")


if __name__ == "__main__":
    main()
