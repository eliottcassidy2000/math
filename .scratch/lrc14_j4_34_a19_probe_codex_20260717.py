#!/usr/bin/env python3
"""Scratch-first exact census for the THM-741 `(3,4)`, a=19 branch."""

from fractions import Fraction as F
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PREVIOUS = ROOT / "04-computation/lrc14_j4_34_a18_pure_tail_exact_codex_20260717.py"
E = tuple(sorted((*range(8, 15), 3, 4)))
A = 19


def load_previous():
    spec = importlib.util.spec_from_file_location("lrc14_j4_34_a18_probe_helpers", PREVIOUS)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    previous = load_previous()
    core = previous.load_core()
    good0, r0, m0 = core.good_norm(E)
    r1, m1, good1 = core.subtract(good0, A)
    _, _, full1 = core.good_norm(tuple(sorted(E + (A,))))
    assert m1 == full1 == core.subtract_sparse(good0, A) > 0

    r18, m18, good18 = core.subtract(good0, 18)
    common = previous.intersection_measure(good18, good1)
    only18, only19 = m18 - common, m1 - common
    V2 = core.minV(3, *(F(4) * m1 / (core.S2 * r1)).as_integer_ratio())

    counts = dict(c=0, preclosed=0, exact=0, closed=0, fallback=0,
                  integral=0, candidates=0, sweeps=0)
    minimum = None
    minimum_family = None
    active = 0
    first_failure = None

    for b in range(A + 1, V2):
        r2, m2, good2 = core.subtract(good1, b)
        assert m2 == core.subtract_sparse(good1, b) > 0
        V3 = core.minV(2, *(F(5) * m2 / (core.S2 * r2)).as_integer_ratio())
        values = tuple(previous.p2_slack(core, m2, r2, c) for c in range(b + 1, V3))
        first = next((i for i, value in enumerate(values) if value > 0), len(values))
        assert all(value <= 0 for value in values[:first])
        assert all(value > 0 for value in values[first:])
        cutoff = b + 1 + first
        counts["c"] += max(0, V3 - b - 1)
        counts["preclosed"] += max(0, V3 - cutoff)
        row_sweeps = 0

        for c in range(b + 1, cutoff):
            counts["exact"] += 1
            sparse_m3 = core.subtract_sparse(good2, c)
            assert sparse_m3 > 0
            denominator = sparse_m3 - m2 / 7
            if denominator > 0:
                cap = core.S2 * r2 / (7 * denominator)
            else:
                r3, m3, _ = core.subtract(good2, c)
                assert m3 == sparse_m3
                cap = core.S2 * r3 / (6 * m3)
                counts["fallback"] += 1
            counts["integral"] += cap.denominator == 1
            d_max = cap.numerator // cap.denominator
            if d_max <= c:
                counts["closed"] += 1
                continue
            counts["candidates"] += 1
            r3, m3, good3 = core.subtract(good2, c)
            assert m3 == sparse_m3 and r3 == len(good3)
            for d in range(c + 1, d_max + 1):
                clearance = core.subtract_sparse(good3, d)
                if clearance <= 0 and first_failure is None:
                    first_failure = (b, c, d, clearance)
                counts["sweeps"] += 1
                row_sweeps += 1
                family = tuple(sorted(E + (A, b, c, d)))
                if minimum is None or clearance < minimum:
                    minimum, minimum_family = clearance, family
        active += row_sweeps > 0

    print(f"root r={r0} m={m0}")
    print(f"a={A} r1={r1} m1={m1} V2={V2}")
    print(f"a18 r={r18} m={m18} common={common} only18={only18} only19={only19}")
    print(f"b={V2-A-1} active={active} counts={counts}")
    print(f"minimum={minimum} family={minimum_family}")
    print(f"first_failure={first_failure}")


if __name__ == "__main__":
    main()
