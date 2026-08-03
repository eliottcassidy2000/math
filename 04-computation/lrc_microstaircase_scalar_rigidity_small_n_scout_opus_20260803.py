#!/usr/bin/env python3
"""Small-n scout for the all-n scalar-ramp rigidity question.

This is deliberately a dependent scout: it imports the direct one-hot
classifier, rebuilds each exact open-cell system for 3 <= n <= 12, and asks
CaDiCaL whether a normalized nonzero full blocker exists.  It does not claim
n=13 and does not turn the finite cell statement into LRC(n).
"""

from __future__ import annotations

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

from pysat.solvers import Solver


ROOT = Path(__file__).resolve().parents[1]
CLASSIFIER = ROOT / (
    "04-computation/"
    "lrc14_microstaircase_scalar_ramp_sat_classification_opus_20260803.py"
)
EXPECTED_PATTERNS = {
    3: 6,
    4: 16,
    5: 30,
    6: 60,
    7: 84,
    8: 144,
    9: 198,
    10: 280,
    11: 352,
    12: 504,
}


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def load_classifier():
    spec = spec_from_file_location("microstaircase_direct_classifier", CLASSIFIER)
    require(spec is not None and spec.loader is not None, CLASSIFIER)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    classifier = load_classifier()
    print("Small-n micro-staircase scalar-rigidity scout")
    print("scope=FINITE-EXACT n=3..12 residue-cell systems; dependent encoding")
    print("n=2;normalized_nonzero=TRIVIALLY_EMPTY")
    for n, expected in EXPECTED_PATTERNS.items():
        patterns = classifier.cell_patterns(n)
        require(len(patterns) == expected, (n, len(patterns), expected))
        clauses = classifier.build_cnf(
            n, patterns, fix_gauge=True, forbid_zero=True
        )
        with Solver(name="cadical195", bootstrap_with=clauses) as solver:
            verdict = solver.solve()
        require(not verdict, f"n={n} has a normalized nonzero full blocker")
        print(
            f"n={n};patterns={len(patterns)};candidates={n * len(patterns)};"
            f"variables={n * (n - 1)};clauses={len(clauses)};"
            "normalized_nonzero=UNSAT"
        )
    print("n=13;status=OPEN_IN_THIS_SCOUT")
    print("all_claimed_dimensions=PASS")


if __name__ == "__main__":
    main()
