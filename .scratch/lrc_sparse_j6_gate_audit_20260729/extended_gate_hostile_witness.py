#!/usr/bin/env python3
"""Extend the failed top-thirty j=6 gate on its first exact hostile root."""

from __future__ import annotations

import hashlib
import importlib.util
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
KERNEL = (
    ROOT
    / "04-computation/lrc14_thm2885_eight_body_top15_hitting_gate_codex_20260729.py"
)
KERNEL_SHA = "dff97f67b1104c25589802a6a2f216b6e7bfedd58eebfa1bcce615d59c1e872f"
BODY = (1, 8, 10, 11, 12, 13, 14)
FIRST = 15
BASE = 1600
COUNT = 40
S2 = F(99, 70)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load():
    require(
        hashlib.sha256(KERNEL.read_bytes()).hexdigest() == KERNEL_SHA,
        "single-comb kernel changed",
    )
    spec = importlib.util.spec_from_file_location("j6_extended_gate", KERNEL)
    require(spec is not None and spec.loader is not None, "cannot load kernel")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


T = load()


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def main() -> None:
    good, components, mass = T.CORE.good_norm(BODY)
    rows = T.coverages_many(good, range(FIRST, BASE + 1))
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    q = ranked[COUNT - 1][0]
    require(q > mass / 7, "rank forty misses discrepancy limit")
    threshold = S2 * components / (7 * (q - mass / 7))
    tail_first = max(BASE + 1, ceiling(threshold))
    rows.extend(
        T.coverages_many(good, range(BASE + 1, tail_first))
    )
    require(
        mass / 7 + S2 * components / (7 * tail_first) <= q,
        "rank-forty tail did not seal",
    )
    top = tuple(sorted(rows, key=lambda item: (-item[0], item[1]))[:COUNT])
    margins = tuple(
        mass
        - sum((value for value, _ in top[k : k + 6]), F(0))
        for k in range(COUNT - 5)
    )
    least = next((k for k, margin in enumerate(margins) if margin > 0), None)
    require(least is not None, "top forty still has no adaptive gate")
    print("J6 EXTENDED GATE HOSTILE WITNESS")
    print(f"E={BODY};h={ftext(mass)};r={components}")
    print(
        f"K24_margin={ftext(margins[24])};"
        f"least_K={least};least_margin={ftext(margins[least])}"
    )
    print(
        "least_complement="
        + ",".join(
            f"{speed}:{ftext(value)}"
            for value, speed in top[least : least + 6]
        )
    )
    print(
        f"rank40_threshold={ftext(threshold)};"
        f"tail_first={tail_first};"
        f"max_top40_speed={max(speed for _, speed in top)}"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
