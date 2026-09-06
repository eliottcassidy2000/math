#!/usr/bin/env python3
"""Exact coefficient box for THM-4437's arbitrary-parity reduction.

This removes the odd-speed parity restriction from THM-4434's coefficient
universe.  The three small full-support circuits are retained as explicit
residuals rather than silently discarded.
"""

import importlib.util
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations_with_replacement, product
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PRIMARY_PATH = ROOT / "04-computation/lrc14_coefficient_box_empty_core_three_ray_sep06.py"
EDGE_PATH = ROOT / "04-computation/lrc14_empty_core_certificate_sep06.py"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


primary = load("coefficient_clipper", PRIMARY_PATH)
edge = load("cube_edge_compiler", EDGE_PATH)
CHECKS = 0


def need(condition, payload):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(payload)


LOW = {(1, 1, 1), (1, 1, 2), (1, 2, 2)}
IMPOSSIBLE = (0, 1, 1)
TARGET = F(15, 98)


def eligible(pattern):
    return (
        sum(x != 0 for x in pattern) >= 2
        and gcd(gcd(pattern[0], pattern[1]), pattern[2]) == 1
        and sum(x % 3 == 0 for x in pattern) <= 1
        and pattern != IMPOSSIBLE
    )


patterns = tuple(
    p for p in combinations_with_replacement(range(19), 3) if eligible(p)
)
second = {
    tuple(sorted(p))
    for p in product(range(19), repeat=3)
    if eligible(tuple(sorted(p)))
}
need(set(patterns) == second, "independent universe construction")
need(len(patterns) == 750, ("all-parity cardinality", len(patterns)))
need(
    Counter(sum(x != 0 for x in p) for p in patterns) == {2: 49, 3: 701},
    "support split",
)

digest = sha256()
leaders = []
low_rows = {}
for pattern in patterns:
    defects, slope, intercept, witness = primary.compile_by_clipping(pattern)
    if pattern[0]:
        alternate = edge.compile_pattern(pattern)
        need(
            (defects, slope, intercept)
            == (alternate["deltas"], alternate["slope"], alternate["intercept"]),
            ("edge disagreement", pattern),
        )
    else:
        need(
            slope == primary.support_two_formula(pattern, defects),
            ("support-two disagreement", pattern),
        )
    digest.update(repr((pattern, defects, slope, intercept)).encode())
    if pattern in LOW:
        low_rows[pattern] = (slope, defects, intercept, witness)
    else:
        need(slope <= TARGET, ("expanded slope bound", pattern, slope, witness))
        leaders.append((slope, pattern, defects, intercept, witness))

leaders.sort(reverse=True)
need(
    leaders[0][0] == TARGET and leaders[0][1] == (1, 7, 8),
    ("sharp leader", leaders[0]),
)
need(set(low_rows) == LOW, ("low residual", low_rows))

print("ALL_PARITY_COEFFICIENT_BOX")
print("UNIVERSE", len(patterns), "support", Counter(sum(x != 0 for x in p) for p in patterns))
print("RESIDUAL_LOW")
for pattern in sorted(low_rows):
    print(pattern, low_rows[pattern])
print("TOP_NONLOW")
for row in leaders[:20]:
    print(row)
print("MAX_NONLOW", leaders[0][0], "unique_pattern", leaders[0][1])
print("SEMANTIC_SHA256", digest.hexdigest())
print("CHECKS", CHECKS, "CLIPPER_CHECKS", primary.CHECKS, "EDGE_CHECKS", edge.CHECKS)
print("PASS: every max-coefficient<=18 pattern outside the three low circuits has slope<=15/98")
print("SCOPE: coefficient box only; the low circuits and finite speed head remain separate obligations")
