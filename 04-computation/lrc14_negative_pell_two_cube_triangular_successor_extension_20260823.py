#!/usr/bin/env python3
"""Finite-exact negative-Pell extension of the two-cube triangular tag."""

from __future__ import annotations

import ast
import hashlib
import json
from math import isqrt
from pathlib import Path


PELL_CAP = 1_000_000
DIRECT_CAP = 5_000
CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def triangular(n: int) -> int:
    return n * (n + 1) // 2


def cube_root_floor(n: int) -> int:
    """Return floor(cuberoot(n)) using integer Newton iteration."""
    if n < 0:
        raise ValueError("cube_root_floor expects n>=0")
    if n < 2:
        return n
    x = 1 << ((n.bit_length() + 2) // 3)
    while True:
        y = (2 * x + n // (x * x)) // 3
        if y >= x:
            while (x + 1) ** 3 <= n:
                x += 1
            while x**3 > n:
                x -= 1
            return x
        x = y


def pell_nodes(cap: int) -> tuple[list[tuple[int, int, int, int, int]], int]:
    """Enumerate all relevant positive odd X^2-3Y^2=-2 solutions."""
    max_value = (cap - 1) ** 3 + 1
    nodes: list[tuple[int, int, int, int, int]] = []
    X, Y = 1, 1
    while True:
        check(X * X - 3 * Y * Y == -2, "negative-Pell recurrence stays on shell")
        check(X % 2 == 1 and Y % 2 == 1, "negative-Pell coordinates stay odd")
        N, M = (Y - 1) // 2, (X - 1) // 2
        value = triangular(N)
        check(3 * value == triangular(M), "Pell node decodes triangular successor")
        if value > max_value:
            break
        nodes.append((X, Y, N, M, value))
        X, Y = 2 * X + 3 * Y, X + 2 * Y
    check(nodes[-1][4] <= max_value, "last retained Pell node is in value range")
    check(value > max_value, "first omitted Pell node exceeds the value range")
    return nodes, max_value


def pell_first_scan(cap: int) -> tuple[
    list[tuple[int, int, int, int]],
    list[tuple[int, int, int, int, int]],
    int,
    int,
]:
    """Search every cube pair on every relevant Pell target."""
    nodes, max_value = pell_nodes(cap)
    hits: list[tuple[int, int, int, int]] = []
    probes = 0
    for _X, _Y, N, M, value in nodes:
        top = min((cap - 1) // 2, cube_root_floor(value // 2))
        for a in range(1, top + 1):
            probes += 1
            remainder = value - a**3
            b = cube_root_floor(remainder)
            check(
                b**3 <= remainder < (b + 1) ** 3,
                "integer cube-root probe is exact",
            )
            if b**3 == remainder and a < b and a + b <= cap:
                hits.append((a, b, N, M))
    return hits, nodes, probes, max_value


def direct_pair_scan(cap: int) -> tuple[list[tuple[int, int, int, int]], int]:
    """Independent pair-first scan, without using the Pell recurrence."""
    hits: list[tuple[int, int, int, int]] = []
    pairs = 0
    for a in range(1, (cap - 1) // 2 + 1):
        a3 = a**3
        for b in range(a + 1, cap - a + 1):
            pairs += 1
            value = a3 + b**3
            disc = 8 * value + 1
            y = isqrt(disc)
            if y * y != disc:
                continue
            triple_disc = 24 * value + 1
            x = isqrt(triple_disc)
            if x * x != triple_disc:
                continue
            N, M = (y - 1) // 2, (x - 1) // 2
            check(triangular(N) == value, "direct triangular discriminant decodes")
            check(triangular(M) == 3 * value, "direct successor discriminant decodes")
            hits.append((a, b, N, M))
    return hits, pairs


pell_hits, nodes, pell_probes, max_value = pell_first_scan(PELL_CAP)
direct_hits, direct_pairs = direct_pair_scan(DIRECT_CAP)

expected_hit = [(9, 13, 76, 132)]
check(pell_hits == expected_hit, "million-shell Pell-first scan has one hit")
check(direct_hits == expected_hit, "independent direct scan has the same one hit")
check(
    [hit for hit in pell_hits if hit[0] + hit[1] <= DIRECT_CAP] == direct_hits,
    "Pell-first and pair-first searches agree on their common universe",
)
check(len(nodes) == 17, "million-shell has seventeen relevant Pell nodes")
check(pell_probes == 731_088, "million-shell candidate cube-probe count")
check(direct_pairs == 6_247_500, "direct-control pair count")
check(nodes[-1][:4] == (1_934_726_305, 1_117_014_753, 558_507_376, 967_363_152),
      "last relevant Pell address")

# Classical square tag: this is a genuine two-cube square, but it does not
# lie in the triangular-successor fibre.
square_value = 56**3 + 65**3
check(square_value == 671**2, "classical square-tag positive control")
square_disc = 8 * square_value + 1
check(isqrt(square_disc) ** 2 != square_disc,
      "square tag is not also a triangular tag")
check(9**3 + 13**3 == triangular(76), "triangular-tag positive control")
check(3 * triangular(76) == triangular(132), "triangular successor positive control")
check(isqrt(triangular(76)) ** 2 != triangular(76),
      "triangular-successor hit is not a square tag")

source = Path(__file__).read_text(encoding="utf-8")
check(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
      "no inactive Python assert")

semantic = {
    "classification": "FINITE-EXACT",
    "universe": "positive distinct a<b with a+b<=1000000",
    "predicate": "a^3+b^3=T_N and 3*T_N=T_M",
    "pell_map": "X=2M+1,Y=2N+1; X^2-3Y^2=-2",
    "hits": expected_hit,
    "independent_control": "direct pair-first scan through a+b<=5000",
    "square_hostile": "56^3+65^3=671^2 is not triangular",
    "scope": "finite address tag only; no LRC(14) implication",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

last_X, last_Y, last_N, last_M, _last_value = nodes[-1]
print("experiment=LRC14-negative-Pell-two-cube-triangular-successor-extension")
print("classification=FINITE-EXACT")
print("universe=positive_distinct_a_lt_b;a_plus_b_le_1000000")
print("predicate=a3_plus_b3_equals_T_N;three_T_N_equals_T_M")
print("pell_equation=X2_minus_3Y2_equals_minus_2;X_equals_2M_plus_1;Y_equals_2N_plus_1")
print(f"pell_nodes={len(nodes)};candidate_cube_probes={pell_probes}")
print(f"largest_retained_node=X:{last_X},Y:{last_Y},N:{last_N},M:{last_M}")
print("pell_hits=(a,b,N,M):(9,13,76,132)")
print(f"direct_control=a_plus_b_le_{DIRECT_CAP};pairs={direct_pairs};same_unique_hit=1")
print("square_hostile=56^3+65^3=671^2;not_triangular")
print("tag_separation=triangular_successor_hit_is_not_square")
print("scope=finite_address_tag_only;no_LRC14_consequence")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print("RESULT PASS")
