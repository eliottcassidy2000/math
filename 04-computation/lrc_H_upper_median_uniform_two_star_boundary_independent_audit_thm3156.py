#!/usr/bin/env python3
r"""Independent audit of the THM-3156 uniform two-star boundary companion.

The main referee uses the proved arc engine and a top-four matching lemma.
This audit does not import either the main referee or that engine.  It solves
the modular interval inequality directly and uses a full 16-state level-scan
matching DP, plus complete brute-force controls at D=6,7.
"""

from __future__ import annotations

import argparse
import hashlib
import re
from ast import literal_eval
from fractions import Fraction as F
from itertools import permutations
from math import gcd
from pathlib import Path


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "00-navigation").is_dir())
MAIN = ROOT / "04-computation/lrc_H_upper_median_uniform_two_star_boundary_thm3156.py"
MAIN_OUTPUT = ROOT / "05-knowledge/results/lrc_H_upper_median_uniform_two_star_boundary_thm3156.out"
OUTPUT = ROOT / "05-knowledge/results/lrc_H_upper_median_uniform_two_star_boundary_independent_audit_thm3156.out"

EXPECTED_MAIN_SHA256 = "0f8e87a4f00b1af78aeee2a8d2160e31830a9e464b01e6abef28a7be02395454"
EXPECTED_MAIN_OUTPUT_SHA256 = "c477a863d8cd211984343e891a4924d1173040b40e81c7db2f75820a443b788d"
EXPECTED_DATA_SEMANTIC = "3d86198c361d1d62322489b37be3b759bcc7a8ac6efdff22858a1f619cb7fbae"
EXPECTED_GLOBAL_MINIMUM = F(2942257539652583, 320491868062929255)
EXPECTED_CONSERVATIVE_MINIMUM = 85295963189446

H = (1, 2, 3, 4, 6, 12)
L = 168
CELL = 90
RADIUS_NUMERATOR = L // 14
SCALE = 10**15
PIVOTS = (0, 1)
LEAVES = (2, 3, 4, 5)
EDGES = ((0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
         (1, 2), (1, 3), (1, 4), (1, 5))
MASK_CONTROLS = (6, 7, 12, 13, 37, 62, 100)
PHASE_BLIND_WITNESS = (136, 119, 102, 103, 101, 100)
RAY_CORRECTIONS = (0, -336, 84)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def reflected_arcs(e: int, q: int):
    """Solve ||((qL-e)(CELL+u))/L||<1/14 directly on 0<=u<=1."""
    z = q * L - e
    rows = []
    for n in range(-10, q + 10):
        centre = e * CELL + n * L
        left = max(F(0), F(centre - RADIUS_NUMERATOR, z))
        right = min(F(1), F(centre + RADIUS_NUMERATOR, z))
        if left < right:
            rows.append((left, right))
    rows.sort()
    merged = []
    for left, right in rows:
        if merged and left <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], right))
        else:
            merged.append((left, right))
    return tuple(merged)


def mass(rows) -> F:
    return sum((right - left for left, right in rows), F(0))


def intersection_mass(first, second) -> F:
    i = 0
    j = 0
    total = F(0)
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            total += right - left
        if first[i][1] <= second[j][1]:
            i += 1
        else:
            j += 1
    return total


def strip_intersection(e: int, f: int, p: int, q: int) -> F:
    """Independent integer-tent formula (25), with p<=q<=2p."""
    require(1 <= p <= q <= 2 * p and e != f, (e, f, p, q))
    z = L * p - e
    w = L * q - f
    r = (CELL * e) % L
    s = (CELL * f) % L
    cross = r * w - s * z
    require(cross % L == 0, ("nonintegral strip offset", e, f, p, q, cross))
    B = cross // L
    numerator = 0
    for k in range(p):
        shifted = B + k * w
        lower = shifted // z
        for ell in (lower, lower + 1):
            if not 0 <= ell < q:
                continue
            distance = abs(shifted - ell * z)
            numerator += (
                max(0, 12 * (z + w) - L * distance)
                - max(0, 12 * abs(z - w) - L * distance)
            )
    return F(numerator, z * w)


def ray_formula(g: int) -> F:
    require(g >= 1, g)
    numerator = 4284 * g * g - 2520 * g + RAY_CORRECTIONS[g % 3]
    return F(numerator, (840 * g - 1) * (504 * g - 6))


def polynomial_add(first, second, sign: int = 1):
    size = max(len(first), len(second))
    return tuple(
        (first[i] if i < len(first) else 0)
        + sign * (second[i] if i < len(second) else 0)
        for i in range(size)
    )


def polynomial_multiply(first, second):
    result = [0] * (len(first) + len(second) - 1)
    for i, a in enumerate(first):
        for j, b in enumerate(second):
            result[i + j] += a * b
    return tuple(result)


def ray_difference_polynomial(residue: int):
    """Cleared numerator of I_(g+1)-I_g as a polynomial in g."""
    current_numerator = (RAY_CORRECTIONS[residue], -2520, 4284)
    next_numerator = (1764 + RAY_CORRECTIONS[(residue + 1) % 3], 6048, 4284)
    current_denominator = (6, -5544, 423360)
    next_denominator = (417822, 841176, 423360)
    return polynomial_add(
        polynomial_multiply(next_numerator, current_denominator),
        polynomial_multiply(current_numerator, next_denominator),
        sign=-1,
    )


def phase_blind_floor(p: int, q: int) -> F:
    divisor = gcd(p, q)
    P, Q = sorted((p // divisor, q // divisor))
    if P + Q <= 7:
        return F(0)
    return max(F(1, 105), F(1, 49) - F(1, 2 * P * Q))


def phase_blind_edge_floor(a: int, p: int, b: int, q: int) -> F:
    """The standard oriented one-edge transport floor used after THM-2941."""
    phase = phase_blind_floor(p, q)
    if phase == 0:
        return F(0)
    rows = []
    for aa, pp, bb, qq in ((a, p, b, q), (b, q, a, p)):
        prefactor = F(pp * L - aa, pp * L)
        eta = F(qq * aa - pp * bb, pp * L - aa)
        rows.append(max(F(0), (phase - 2 * abs(eta)) / prefactor))
    return max(rows)


def floor_scaled(value: F) -> int:
    return value.numerator * SCALE // value.denominator


def ceil_scaled(value: F) -> int:
    return -((-value.numerator * SCALE) // value.denominator)


def banks(D: int):
    levels = tuple(range(D, 2 * D + 1))
    arcs = {(i, q): reflected_arcs(H[i], q) for i in range(6) for q in levels}
    gains = {
        (i, j, p, q): intersection_mass(arcs[i, p], arcs[j, q])
        for i, j in EDGES
        for p in levels
        for q in levels
        if p != q
    }
    debts = {
        (i, q): F(H[i], 7 * (q * L - H[i]))
        for i in range(6)
        for q in levels
    }
    return levels, gains, debts


def objective(assignment, gains, debts, conservative: bool):
    if conservative:
        return (
            sum(floor_scaled(gains[i, j, assignment[i], assignment[j]]) for i, j in EDGES)
            - 9 * sum(ceil_scaled(debts[i, assignment[i]]) for i in range(6))
        )
    return (
        sum((gains[i, j, assignment[i], assignment[j]] for i, j in EDGES), F(0)) / 9
        - sum((debts[i, assignment[i]] for i in range(6)), F(0))
    )


def full_mask_minimum(D: int, conservative: bool):
    """Full level-scan subset DP, without the main referee's top-four lemma."""
    levels, gains, debts = banks(D)
    best = None
    for q0, q1 in permutations(levels, 2):
        base = (
            floor_scaled(gains[0, 1, q0, q1])
            - 9 * (ceil_scaled(debts[0, q0]) + ceil_scaled(debts[1, q1]))
            if conservative
            else gains[0, 1, q0, q1] / 9 - debts[0, q0] - debts[1, q1]
        )
        costs = {
            (leaf - 2, q): (
                floor_scaled(gains[0, leaf, q0, q])
                + floor_scaled(gains[1, leaf, q1, q])
                - 9 * ceil_scaled(debts[leaf, q])
                if conservative
                else (gains[0, leaf, q0, q] + gains[1, leaf, q1, q]) / 9
                - debts[leaf, q]
            )
            for leaf in LEAVES
            for q in levels
            if q not in (q0, q1)
        }
        dp = [None] * 16
        dp[0] = 0
        for q in levels:
            if q in (q0, q1):
                continue
            nxt = list(dp)
            for mask, value in enumerate(dp):
                if value is None:
                    continue
                for leaf_index in range(4):
                    bit = 1 << leaf_index
                    if mask & bit:
                        continue
                    candidate = value + costs[leaf_index, q]
                    old = nxt[mask | bit]
                    if old is None or candidate < old:
                        nxt[mask | bit] = candidate
            dp = nxt
        require(dp[15] is not None, ("mask empty", D, q0, q1))
        candidate = base + dp[15]
        if best is None or candidate < best:
            best = candidate
    require(best is not None, ("global empty", D))
    return best


def brute_minimum(D: int):
    levels, gains, debts = banks(D)
    return min(objective(assignment, gains, debts, False)
               for assignment in permutations(levels, 6))


def parse_main_output():
    lines = MAIN_OUTPUT.read_text(encoding="utf-8").splitlines()
    require(len(lines) == 110, ("line count", len(lines)))
    require(lines[0] == "LRC H UNIFORM TWO-STAR BOUNDARY -- THM-3156 EXACT COMPANION", lines[0])
    data_lines = tuple(line for line in lines if line.startswith("D="))
    require(len(data_lines) == 95, len(data_lines))
    pattern = re.compile(
        r"^D=(\d+);m=(\d+);pivot_checks=(\d+);"
        r"(exact_minimum|conservative_floor_for_9S_margin)=([^;]+);"
        r"representative_margin=([^;]+);levels=(\([^)]*\));"
        r"nine_edge_gain_sum=([^;]+);debt=(.+)$"
    )
    digest = hashlib.sha256()
    rows = {}
    for expected_D, line in zip(range(6, 101), data_lines):
        digest.update((line + "\n").encode())
        match = pattern.match(line)
        require(match is not None, ("parse", line))
        D, m, checks = map(int, match.group(1, 2, 3))
        kind = match.group(4)
        score = F(match.group(5)) if kind == "exact_minimum" else int(match.group(5))
        margin = F(match.group(6))
        assignment = literal_eval(match.group(7))
        gain_sum = F(match.group(8))
        debt_sum = F(match.group(9))
        require((D, m, checks) == (expected_D, D, D * (D + 1)), (D, m, checks))
        require(kind == ("exact_minimum" if D <= 12 else
                         "conservative_floor_for_9S_margin"), (D, kind))
        require(len(assignment) == len(set(assignment)) == 6, (D, assignment))
        require(all(D <= q <= 2 * D for q in assignment), (D, assignment))
        direct_gains = sum(
            (
                intersection_mass(
                    reflected_arcs(H[i], assignment[i]),
                    reflected_arcs(H[j], assignment[j]),
                )
                for i, j in EDGES
            ),
            F(0),
        )
        direct_debt = sum(
            (F(H[i], 7 * (assignment[i] * L - H[i])) for i in range(6)),
            F(0),
        )
        require((gain_sum, debt_sum) == (direct_gains, direct_debt),
                ("witness components", D))
        require(margin == gain_sum / 9 - debt_sum > 0, ("witness margin", D))
        if D <= 12:
            require(score == margin, ("exact witness", D, score, margin))
        else:
            rounded = (
                sum(
                    floor_scaled(
                        intersection_mass(
                            reflected_arcs(H[i], assignment[i]),
                            reflected_arcs(H[j], assignment[j]),
                        )
                    )
                    for i, j in EDGES
                )
                - 9 * sum(
                    ceil_scaled(F(H[i], 7 * (assignment[i] * L - H[i])))
                    for i in range(6)
                )
            )
            require(score == rounded <= 9 * SCALE * margin,
                    ("rounding", D, score, rounded, margin))
        rows[D] = (kind, score, margin, assignment)
    require(digest.hexdigest() == EXPECTED_DATA_SEMANTIC,
            (digest.hexdigest(), EXPECTED_DATA_SEMANTIC))
    require(f"data_semantic_sha256={EXPECTED_DATA_SEMANTIC}" in lines, "semantic token")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    require(sha256(MAIN) == EXPECTED_MAIN_SHA256,
            ("main changed", sha256(MAIN), EXPECTED_MAIN_SHA256))
    require(sha256(MAIN_OUTPUT) == EXPECTED_MAIN_OUTPUT_SHA256,
            ("main output changed", sha256(MAIN_OUTPUT), EXPECTED_MAIN_OUTPUT_SHA256))
    main_text = MAIN.read_text(encoding="utf-8")
    require("lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.py" in main_text,
            "proved engine pin missing")
    require("reflected_one_cone" not in main_text and
            "lrc_H_upper_median_two_star_matching_scout" not in main_text,
            "unproved dependency leaked into main")

    rows = parse_main_output()
    for e in H:
        for q in (6, 13, 62, 100, 200):
            require(mass(reflected_arcs(e, q)) ==
                    F(1, 7) + F(e, 7 * (q * L - e)),
                    ("singleton mass", e, q))

    brute_rows = tuple((D, brute_minimum(D)) for D in (6, 7))
    require(all(value == rows[D][1] for D, value in brute_rows), brute_rows)
    mask_rows = tuple(
        (D, full_mask_minimum(D, D >= 13))
        for D in MASK_CONTROLS
    )
    require(all(value == rows[D][1] for D, value in mask_rows), mask_rows)

    # Independently audit the exact digital-strip formula in theorem section 5.
    strip_checks = 0
    for e in H:
        for f in H:
            if e == f:
                continue
            for p in range(1, 21):
                for q in range(p, 2 * p + 1):
                    direct = intersection_mass(reflected_arcs(e, p), reflected_arcs(f, q))
                    stripped = strip_intersection(e, f, p, q)
                    require(direct == stripped,
                            ("integer-strip formula", e, f, p, q, direct, stripped))
                    strip_checks += 1
    require(strip_checks == 6900, strip_checks)

    ray_rows = []
    previous = None
    for g in range(1, 101):
        direct = intersection_mass(reflected_arcs(6, 3 * g), reflected_arcs(1, 5 * g))
        formula = ray_formula(g)
        require(direct == formula, ("periodic ray", g, direct, formula))
        if previous is not None:
            require(formula > previous, ("ray monotonicity", g, previous, formula))
        previous = formula
        ray_rows.append((g, formula))
    expected_difference_polynomials = (
        (504 * 17, 504 * 2073474, 504 * 1787436, 0, 0),
        (1008 * 139285, 1008 * 1314819, 1008 * 1211238, 0, 0),
        (1008 * -34808, 1008 * 964791, 1008 * 999558, 0, 0),
    )
    actual_difference_polynomials = tuple(ray_difference_polynomial(r) for r in range(3))
    require(actual_difference_polynomials == expected_difference_polynomials,
            (actual_difference_polynomials, expected_difference_polynomials))
    # Each quadratic is positive for every real g>=1: the linear and quadratic
    # coefficients are positive, and the value at one is positive.
    require(all(poly[2] > 0 and poly[1] > 0 and sum(poly[:3]) > 0
                for poly in actual_difference_polynomials),
            actual_difference_polynomials)
    require(ray_rows[1][1] == F(2030, 280393) > F(1, 140), ray_rows[1])

    phase_blind_rows = tuple(
        (
            edge,
            phase_blind_floor(PHASE_BLIND_WITNESS[edge[0]],
                              PHASE_BLIND_WITNESS[edge[1]]),
            phase_blind_edge_floor(
                H[edge[0]], PHASE_BLIND_WITNESS[edge[0]],
                H[edge[1]], PHASE_BLIND_WITNESS[edge[1]],
            ),
        )
        for edge in EDGES
    )
    require(all(row[2] == 0 for row in phase_blind_rows), phase_blind_rows)
    phase_blind_debt = sum(
        (F(e, 7 * (q * L - e)) for e, q in zip(H, PHASE_BLIND_WITNESS)),
        F(0),
    )
    require(phase_blind_debt ==
            F(656518999793763784, 2839155148366624716475) > 0,
            phase_blind_debt)

    exact_global = min((rows[D][1], D, rows[D][3]) for D in range(6, 13))
    conservative_global = min((rows[D][1], D) for D in range(13, 101))
    require(exact_global == (EXPECTED_GLOBAL_MINIMUM, 6, (6, 9, 12, 11, 8, 7)),
            exact_global)
    require(conservative_global == (EXPECTED_CONSERVATIVE_MINIMUM, 24),
            conservative_global)
    require(F(EXPECTED_CONSERVATIVE_MINIMUM, 9 * SCALE) > EXPECTED_GLOBAL_MINIMUM,
            (conservative_global, exact_global))

    lines = [
        "LRC H UNIFORM TWO-STAR BOUNDARY -- THM-3156 INDEPENDENT AUDIT",
        "status=INDEPENDENTLY HOSTILE-AUDITED;direct modular arcs;full 16-state matching DP",
        "main_witness_rows=95;all exact gain/debt/margin rows reconstructed",
        f"brute_controls={tuple(D for D, _ in brute_rows)};mask_controls={MASK_CONTROLS}",
        f"integer_strip_controls={strip_checks};ordered_distinct_labels=30;p=1..20;p<=q<=2p",
        "periodic_ray_controls=g=1..100;direct_equals_closed_form;strictly_increasing",
        f"ray_difference_polynomials={actual_difference_polynomials};positive_for_all_g>=1",
        f"ray_floor=I2:{qtext(ray_rows[1][1])}>1/140",
        f"phase_blind_zero_witness={PHASE_BLIND_WITNESS};edge_floors={tuple(row[2] for row in phase_blind_rows)};debt={qtext(phase_blind_debt)}",
        f"global_exact_minimum={qtext(EXPECTED_GLOBAL_MINIMUM)}@D=6,levels=(6, 9, 12, 11, 8, 7)",
        f"smallest_conservative_9S_floor={EXPECTED_CONSERVATIVE_MINIMUM}@D=24;rounding_direction=PASS",
        "dependency_audit=main imports only proved reflected arc engine;unproved one-cone/scout absent",
        "scope=fixed reflected H packet,D=6..100,m=D;not arbitrary D;not physical LRC14",
        "reproduction=python3 04-computation/lrc_H_upper_median_uniform_two_star_boundary_independent_audit_thm3156.py",
        f"main_sha256={EXPECTED_MAIN_SHA256}",
        f"main_output_sha256={EXPECTED_MAIN_OUTPUT_SHA256}",
        f"main_data_semantic_sha256={EXPECTED_DATA_SEMANTIC}",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"source_sha256={sha256(HERE)}",
        "all_independent_controls=PASS",
    ]
    text = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
