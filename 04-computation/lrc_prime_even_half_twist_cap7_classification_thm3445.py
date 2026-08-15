#!/usr/bin/env python3
"""Exact companion for the provisional prime-even cap-seven classification.

The all-prime argument has an infinite elementary part and three finite gates.
This companion checks the finite gates without using floating point:

* the parity-constrained line-section lower bounds for every reduced relation
  with 1 <= a < 100 and 1 <= b <= 13, plus the uniform a >= 100 plateau;
* the A/A orientation branches and the three exact height-eight tangent counts;
* the typed height-seven zero graph in every unit class modulo 240, including
  its clique-six equality packet and direct-mask hostile controls;
* the cover-profile overlap budgets, the p=601 negative boundary, and the two
  positive witnesses Q=14,38.

The script pins, but does not silently replace, the two predecessor bounded
companions.  THM-3445 remains a provisional proof candidate until an
independent full-package audit promotes it.
"""

from __future__ import annotations

import ast
from collections import Counter
from hashlib import sha256
from math import gcd, isqrt
from pathlib import Path
import runpy


ROOT = Path(__file__).resolve().parents[1]
FINITE_PATH = ROOT / "04-computation/lrc_prime_even_half_twist_cap7_finite_atlas_20260815.py"
FINITE_OUTPUT = ROOT / "05-knowledge/results/lrc_prime_even_half_twist_cap7_finite_atlas_20260815.out"
TAIL_PATH = ROOT / "04-computation/lrc_prime_even_half_twist_low_weight_tail_probe_20260815.py"
TAIL_OUTPUT = ROOT / "05-knowledge/results/lrc_prime_even_half_twist_low_weight_tail_probe_20260815.out"
EXPECTED_FINITE_SHA256 = "ed54b01f9bf155643b6407c6af8ee6f15c7a099f8ac3b60399dd49b496fb1d12"
EXPECTED_FINITE_OUTPUT_SHA256 = "9751da7e0464f16f92b675ca9c95b118fd3fe95fc11cc99c9a7015d8425adb65"
EXPECTED_TAIL_SHA256 = "4bb0c853fa3bebabf1414d75aa2b43d4b00c22a2b49cbaaac5e06bde5986c4ba"
EXPECTED_TAIL_OUTPUT_SHA256 = "a9012981cef620a0eeac6883f6e6cb080272c3d99cb73703a0b45212f9cb9b80"
EXPECTED_SEMANTIC_SHA256 = "45e7e843b1767ef7254be729add71ea6ef501160503b46a03cbb43e3e53517c3"

ALPHA = {"A": 4, "B": 2, "E": 1}
DIVISOR = {"A": 1, "B": 2, "E": 4}
CROSS_SECTORS = ("AB", "AE", "BA", "BB", "BE", "EA", "EB")
POSITIVE_CONTROLS = {
    7: (1, 3, 4, 5, 9, 11, 13),
    19: (1, 9, 17, 20, 21, 29, 37),
}


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def is_prime(value: int) -> bool:
    return value >= 2 and all(
        value % divisor for divisor in range(2, isqrt(value) + 1)
    )


def first_primes_in_class(residue: int, modulus: int, lower: int,
                          count: int = 2) -> tuple[int, ...]:
    values = []
    candidate = lower + (residue - lower) % modulus
    while len(values) < count:
        if is_prime(candidate):
            values.append(candidate)
        candidate += modulus
    return tuple(values)


def coefficient_type(coefficient: int) -> str:
    if coefficient % 2:
        return "A"
    if coefficient % 4:
        return "B"
    return "E"


def danger_mask(prime: int, coefficient: int) -> int:
    period = 4 * prime
    mask = 0
    for sheet in range(2 * prime):
        phase = coefficient * (2 * sheet + 1) % period
        if 14 * min(phase, period - phase) < period:
            mask |= 1 << sheet
    return mask


def literal_weight(prime: int, first: int, second: int,
                   masks: dict[int, int] | None = None) -> int:
    if masks is None:
        masks = {
            first: danger_mask(prime, first),
            second: danger_mask(prime, second),
        }
    return (masks[first] & masks[second]).bit_count()


def section_shape(left: int, right: int, line: int) -> int:
    """Scaled strict x-section length for two centered open intervals."""
    return max(
        0,
        min(left, 14 * line + right) - max(-left, 14 * line - right),
    )


def allowed_line(kind: str, line: int) -> bool:
    if kind == "all":
        return True
    if kind == "odd":
        return line % 2 == 1
    if kind == "even":
        return line % 2 == 0
    if kind == "mod4-0":
        return line % 4 == 0
    if kind == "mod4-2":
        return line % 4 == 2
    raise RuntimeError(("line kind", kind))


def cross_rule(first_type: str, second_type: str, a: int, b: int):
    """Return mismatch, admissible-line class, AP multiplier, deck factor."""
    if first_type in "AB" and second_type in "AB":
        mismatch = (a + b) % 2 == 1
        kind = "odd" if mismatch else "even"
        step_multiplier = 2
    elif first_type in "AB" and second_type == "E":
        mismatch = b % 2 == 0
        kind = "odd" if mismatch else "all"
        step_multiplier = 1 if mismatch else 2
    elif first_type == "E" and second_type in "AB":
        mismatch = a % 2 == 0
        kind = "odd" if mismatch else "all"
        step_multiplier = 1 if mismatch else 2
    else:
        raise RuntimeError(("not a cross sector", first_type, second_type))
    literal_factor = 1 if "A" in (first_type, second_type) else 2
    return mismatch, kind, step_multiplier, literal_factor


def cross_zero(first_type: str, second_type: str, a: int, b: int) -> bool:
    mismatch, _kind, _step, _factor = cross_rule(
        first_type, second_type, a, b
    )
    return mismatch and ALPHA[first_type] * a + ALPHA[second_type] * b <= 14


def line_lower_bound(prime_floor: int, first_type: str, second_type: str,
                     a: int, b: int, kind: str, step_multiplier: int,
                     literal_factor: int) -> int:
    left = ALPHA[first_type] * a
    right = ALPHA[second_type] * b
    diameter = left + right
    denominator = 14 * a * step_multiplier * b
    total = 0
    radius = diameter // 14 + 2
    for line in range(-radius, radius + 1):
        shape = section_shape(left, right, line)
        if shape and allowed_line(kind, line):
            # ceil(prime_floor*shape/denominator)-1, with strict endpoints.
            total += (prime_floor * shape - 1) // denominator
    return literal_factor * total


def cross_lower_bound(first_type: str, second_type: str, a: int, b: int) -> int:
    mismatch, kind, step, factor = cross_rule(first_type, second_type, a, b)
    del mismatch
    prime_floor = max(607, 14 * a + 1)
    return line_lower_bound(
        prime_floor, first_type, second_type, a, b, kind, step, factor
    )


EXPECTED_CROSS_ZEROS = {
    "AB": ((1, 2), (1, 4), (2, 1), (2, 3)),
    "AE": ((1, 2), (1, 4), (1, 6), (1, 8), (1, 10), (3, 2)),
    "BA": ((1, 2), (2, 1), (3, 2), (4, 1)),
    "BB": (
        (1, 2), (1, 4), (1, 6), (2, 1), (2, 3), (2, 5),
        (3, 2), (3, 4), (4, 1), (4, 3), (5, 2), (6, 1),
    ),
    "BE": (
        (1, 2), (1, 4), (1, 6), (1, 8), (1, 10), (1, 12),
        (3, 2), (3, 4), (3, 8), (5, 2), (5, 4),
    ),
    "EA": ((2, 1), (2, 3), (4, 1), (6, 1), (8, 1), (10, 1)),
    "EB": (
        (2, 1), (2, 3), (2, 5), (4, 1), (4, 3), (4, 5),
        (6, 1), (8, 1), (8, 3), (10, 1), (12, 1),
    ),
}


def finite_cross_certificate():
    expected_minima = {
        "AB": 14, "AE": 13, "BA": 14, "BB": 12,
        "BE": 8, "EA": 14, "EB": 8,
    }
    expected_low = {
        "AB": (), "AE": (), "BA": (), "BB": (),
        "BE": ((8, 3, 10), (8, 5, 6)),
        "EA": (),
        "EB": ((8, 6, 5), (8, 10, 3)),
    }
    reports = []
    for sector in CROSS_SECTORS:
        first_type, second_type = sector
        rows = []
        zeros = []
        for a in range(1, 100):
            for b in range(1, 14):
                if gcd(a, b) != 1:
                    continue
                bound = cross_lower_bound(first_type, second_type, a, b)
                if cross_zero(first_type, second_type, a, b):
                    require(bound == 0, ("cross zero bound", sector, a, b, bound))
                    zeros.append((a, b))
                else:
                    rows.append((bound, a, b))
        require(tuple(zeros) == EXPECTED_CROSS_ZEROS[sector],
                ("cross zero atlas", sector, zeros))
        minimum = min(bound for bound, _a, _b in rows)
        low_rows = tuple(row for row in rows if row[0] <= 8)
        require(minimum == expected_minima[sector],
                ("cross minimum", sector, minimum))
        require(low_rows == expected_low[sector],
                ("cross low rows", sector, low_rows))
        reports.append((sector, len(rows) + len(zeros), minimum, low_rows))
    return tuple(reports)


def aa_rule(branch: str):
    if branch == "opposite":
        return "odd", 4
    if branch == "kappa0":
        return "mod4-0", 2
    if branch == "kappa2":
        return "mod4-2", 2
    raise RuntimeError(("A/A branch", branch))


def aa_zero(branch: str, a: int, b: int) -> bool:
    if branch == "opposite":
        return a + b <= 3
    if branch == "kappa0":
        return False
    if branch == "kappa2":
        return a + b <= 7
    raise RuntimeError(("A/A branch", branch))


def aa_lower_bound(branch: str, a: int, b: int) -> int:
    kind, step = aa_rule(branch)
    return line_lower_bound(
        max(607, 14 * a + 1), "A", "A", a, b, kind, step, 1
    )


def finite_aa_certificate():
    expected_zeros = {
        "opposite": ((1, 2), (2, 1)),
        "kappa0": (),
        "kappa2": ((1, 1), (1, 3), (1, 5), (3, 1), (5, 1)),
    }
    expected_minima = {"opposite": 16, "kappa0": 13, "kappa2": 10}
    expected_low12 = {
        "opposite": (),
        "kappa0": (),
        "kappa2": ((10, 3, 5), (10, 5, 3)),
    }
    reports = []
    for branch in ("opposite", "kappa0", "kappa2"):
        rows = []
        zeros = []
        for a in range(1, 100):
            for b in range(1, 14):
                if gcd(a, b) != 1:
                    continue
                if branch == "opposite" and (a + b) % 2 != 1:
                    continue
                if branch != "opposite" and not (a % 2 and b % 2):
                    continue
                bound = aa_lower_bound(branch, a, b)
                if aa_zero(branch, a, b):
                    require(bound == 0, ("A/A zero bound", branch, a, b, bound))
                    zeros.append((a, b))
                else:
                    rows.append((bound, a, b))
        require(tuple(zeros) == expected_zeros[branch],
                ("A/A zero atlas", branch, zeros))
        minimum = min(bound for bound, _a, _b in rows)
        low12 = tuple(row for row in rows if row[0] < 12)
        require(minimum == expected_minima[branch],
                ("A/A minimum", branch, minimum))
        require(low12 == expected_low12[branch],
                ("A/A low rows", branch, low12))
        reports.append((branch, len(rows) + len(zeros), minimum, low12))
    return tuple(reports)


def plateau_line_count(radius: int, kind: str) -> int:
    return sum(allowed_line(kind, line) for line in range(-radius, radius + 1))


def plateau_cross_bound(first_type: str, second_type: str,
                        a: int, b: int) -> int:
    _mismatch, kind, step, factor = cross_rule(first_type, second_type, a, b)
    left = ALPHA[first_type] * a
    right = ALPHA[second_type] * b
    require(left >= right, ("plateau ordering", first_type, second_type, a, b))
    radius = (left - right) // 14
    per_line = (2 * right) // (step * b)
    return factor * plateau_line_count(radius, kind) * per_line


def plateau_certificate():
    cross_reports = []
    for sector in CROSS_SECTORS:
        first_type, second_type = sector
        rows = tuple(
            (plateau_cross_bound(first_type, second_type, a, b), a, b)
            for a in (100, 101)
            for b in range(1, 14)
        )
        minimum = min(row[0] for row in rows)
        require(minimum > 8, ("large-a cross plateau", sector, minimum))
        cross_reports.append((sector, minimum))

    aa_reports = []
    for branch in ("opposite", "kappa0", "kappa2"):
        kind, step = aa_rule(branch)
        rows = []
        for a in (100, 101):
            for b in range(1, 14):
                if branch == "opposite" and (a + b) % 2 != 1:
                    continue
                if branch != "opposite" and not (a % 2 and b % 2):
                    continue
                radius = (4 * a - 4 * b) // 14
                per_line = (8 * b) // (step * b)
                rows.append((plateau_line_count(radius, kind) * per_line, a, b))
        minimum = min(row[0] for row in rows)
        require(minimum > 8, ("large-a A/A plateau", branch, minimum))
        aa_reports.append((branch, minimum))
    return tuple(cross_reports), tuple(aa_reports)


def count_residue(lower: int, upper: int, residue: int, modulus: int) -> int:
    if lower > upper:
        return 0
    residue %= modulus
    return (upper - residue) // modulus - ((lower - 1) - residue) // modulus


def aa_tangent_count(prime: int) -> int:
    lower = 4 * prime // 21 + 1
    upper = (2 * prime - 1) // 7
    residues = tuple(
        residue for residue in range(10)
        if residue % 2 and (3 * residue - 2 * prime) % 5 == 0
    )
    require(len(residues) == 1, ("A/A tangent residue", prime, residues))
    return count_residue(lower, upper, residues[0], 10)


def be_35_tangent_count(prime: int) -> int:
    lower = 2 * prime // 21 + 1
    upper = (prime - 1) // 7
    residue = 7 * prime % 10
    return count_residue(lower, upper, residue, 10)


def be_53_tangent_count(prime: int) -> int:
    lower = 4 * prime // 35 + 1
    upper = (prime - 1) // 7
    residues = tuple(
        residue for residue in range(6)
        if residue % 2 and (5 * residue - prime) % 6 == 0
    )
    require(len(residues) == 1, ("B/E 5/3 tangent residue", prime, residues))
    return count_residue(lower, upper, residues[0], 6)


def tangent_certificate():
    for odd in range(1, 421, 2):
        require(aa_tangent_count(odd + 210) == aa_tangent_count(odd) + 2,
                ("A/A recurrence", odd))
        require(be_35_tangent_count(odd + 210) == be_35_tangent_count(odd) + 1,
                ("B/E 3/5 recurrence", odd))
        require(be_53_tangent_count(odd + 210) == be_53_tangent_count(odd) + 1,
                ("B/E 5/3 recurrence", odd))

    last_aa_below_six = max(
        odd for odd in range(1, 816, 2) if aa_tangent_count(odd) < 6
    )
    last_35_below_three = max(
        odd for odd in range(1, 816, 2) if be_35_tangent_count(odd) < 3
    )
    last_53_below_three = max(
        odd for odd in range(1, 816, 2) if be_53_tangent_count(odd) < 3
    )
    require((last_aa_below_six, last_35_below_three, last_53_below_three)
            == (605, 601, 609), "tangent odd boundaries")
    require(min(aa_tangent_count(odd) for odd in range(607, 816, 2)) == 6,
            "A/A one-period floor")
    require(min(be_35_tangent_count(odd) for odd in range(607, 816, 2)) == 3,
            "B/E 3/5 one-period floor")
    require(min(
        be_53_tangent_count(odd)
        for odd in range(607, 816, 2) if gcd(odd, 210) == 1
    ) == 3, "B/E 5/3 prime-class floor")

    controls = []
    for prime in (521, 593, 601, 607, 1009):
        controls.append((
            prime,
            2 * aa_tangent_count(prime),
            4 * be_35_tangent_count(prime),
            4 * be_53_tangent_count(prime),
        ))
    return (
        (last_aa_below_six, last_35_below_three, last_53_below_three),
        tuple(controls),
    )


def height_seven_fractions() -> tuple[tuple[int, int], ...]:
    values = []
    for numerator in range(1, 7):
        for denominator in range(1, 8 - numerator):
            if gcd(numerator, denominator) == 1:
                values.extend(((numerator, denominator), (-numerator, denominator)))
    return tuple(values)


HEIGHT_SEVEN = height_seven_fractions()


def coefficient_from_label(prime: int, root: int,
                           label: tuple[int, int, int]) -> int:
    numerator, denominator, lift = label
    residue = numerator * root * pow(denominator, -1, prime) % prime
    require(residue, ("zero labelled coefficient", prime, root, label))
    return residue + lift * prime


VERTEX_LABELS = tuple(
    (numerator, denominator, lift)
    for numerator, denominator in HEIGHT_SEVEN
    for lift in (0, 1)
)
ROOT_LABEL = (1, 1, 0)


def aa_formula_zero(prime: int, first: int, second: int) -> bool:
    modulus = 4 * prime
    normalized = second * pow(first, -1, modulus) % modulus
    for a in range(1, 7):
        for b in range(1, 7):
            if gcd(a, b) != 1 or a + b > 7:
                continue
            for sign in (1, -1):
                difference = (b * normalized - sign * a) % modulus
                if difference % prime:
                    continue
                kappa = difference // prime
                if (a + b) % 2 == 1 and a + b <= 3:
                    return True
                if a % 2 and b % 2 and kappa == 2:
                    return True
    return False


def cross_formula_zero(prime: int, first: int, second: int,
                       first_type: str, second_type: str) -> bool:
    first_effective = first // DIVISOR[first_type]
    second_effective = second // DIVISOR[second_type]
    for a in range(1, 15):
        for b in range(1, 14):
            if gcd(a, b) != 1:
                continue
            if not cross_zero(first_type, second_type, a, b):
                continue
            for sign in (1, -1):
                if (b * second_effective - sign * a * first_effective) % prime == 0:
                    return True
    return False


def formula_zero(prime: int, first: int, second: int) -> bool:
    if first == second:
        return False
    first_type = coefficient_type(first)
    second_type = coefficient_type(second)
    if first_type == second_type == "E":
        return False
    if first_type == second_type == "A":
        return aa_formula_zero(prime, first, second)
    return cross_formula_zero(
        prime, first, second, first_type, second_type
    )


def maximum_rooted_cliques(prime: int, root: int):
    coefficients = {
        label: coefficient_from_label(prime, root, label)
        for label in VERTEX_LABELS
    }
    require(len(set(coefficients.values())) == len(coefficients),
            ("label collision", prime, root))
    require(coefficients[ROOT_LABEL] == root, ("root label", prime, root))

    def adjacent(first_label, second_label):
        return formula_zero(
            prime, coefficients[first_label], coefficients[second_label]
        )

    neighbors = tuple(
        label for label in VERTEX_LABELS
        if label != ROOT_LABEL and adjacent(ROOT_LABEL, label)
    )
    best_size = 1
    best = {tuple([ROOT_LABEL])}

    def extend(chosen, candidates):
        nonlocal best_size, best
        if len(chosen) + len(candidates) < best_size:
            return
        if len(chosen) > best_size:
            best_size = len(chosen)
            best = {tuple(sorted(chosen))}
        elif len(chosen) == best_size:
            best.add(tuple(sorted(chosen)))
        remaining = candidates
        while remaining:
            if len(chosen) + len(remaining) < best_size:
                return
            vertex = remaining[0]
            tail = remaining[1:]
            next_candidates = tuple(
                other for other in tail if adjacent(vertex, other)
            )
            extend(chosen + (vertex,), next_candidates)
            remaining = tail

    extend((ROOT_LABEL,), neighbors)
    coefficient_cliques = {
        tuple(sorted(coefficients[label] for label in clique))
        for clique in best if len(clique) == best_size
    }
    return len(neighbors), best_size, coefficient_cliques, coefficients


def graph_signature(prime: int, root: int):
    coefficients = tuple(
        coefficient_from_label(prime, root, label) for label in VERTEX_LABELS
    )
    types = tuple(coefficient_type(coefficient) for coefficient in coefficients)
    edges = tuple(
        formula_zero(prime, coefficients[i], coefficients[j])
        for i in range(len(coefficients))
        for j in range(i + 1, len(coefficients))
    )
    return types, edges


def mod240_graph_certificate():
    residues = tuple(residue for residue in range(240) if gcd(residue, 240) == 1)
    require(len(residues) == 64, ("mod-240 units", len(residues)))
    reports = []
    for residue in residues:
        first_prime, second_prime = first_primes_in_class(residue, 240, 607)
        root_reports = []
        all_cliques = set()
        for root in (1, 2, 4):
            require(graph_signature(first_prime, root)
                    == graph_signature(second_prime, root),
                    ("mod-240 signature", residue, root, first_prime, second_prime))
            degree, clique_size, cliques, _coefficients = maximum_rooted_cliques(
                first_prime, root
            )
            require(clique_size == 6, ("clique size", residue, root, clique_size))
            root_reports.append((root, degree, clique_size, len(cliques)))
            all_cliques.update(cliques)
        require(tuple(row[1] for row in root_reports) == (19, 31, 23),
                ("zero degrees", residue, root_reports))
        require(tuple(row[3] for row in root_reports) == (2, 1, 1),
                ("rooted equality count", residue, root_reports))
        require(len(all_cliques) == 4, ("equality presentations", residue, all_cliques))
        for clique in all_cliques:
            require(all(2 * first_prime - vertex in clique for vertex in clique),
                    ("complement packets", residue, clique))
            require(tuple(sorted(coefficient_type(vertex) for vertex in clique))
                    == ("A", "A", "A", "A", "B", "E"),
                    ("equality types", residue, clique))
        reports.append((residue, first_prime, second_prime, tuple(root_reports)))
    return tuple(reports)


def direct_graph_control(prime: int):
    all_coefficients = tuple(
        coefficient for coefficient in range(1, 2 * prime)
        if coefficient != prime
    )
    masks = {
        coefficient: danger_mask(prime, coefficient)
        for coefficient in all_coefficients
    }
    root_reports = []
    for root in (1, 2, 4):
        weights = {
            coefficient: literal_weight(prime, root, coefficient, masks)
            for coefficient in all_coefficients if coefficient != root
        }
        for coefficient, weight in weights.items():
            require(formula_zero(prime, root, coefficient) == (weight == 0),
                    ("full-root formula", prime, root, coefficient, weight))
        degree = sum(weight == 0 for weight in weights.values())
        minimum = min(weight for weight in weights.values() if weight)
        root_reports.append((root, degree, minimum))

        _degree, _size, _cliques, labelled = maximum_rooted_cliques(prime, root)
        candidate_values = tuple(labelled.values())
        for index, first in enumerate(candidate_values):
            for second in candidate_values[index + 1:]:
                require(formula_zero(prime, first, second)
                        == (literal_weight(prime, first, second, masks) == 0),
                        ("candidate graph formula", prime, root, first, second))
    return prime, tuple(root_reports)


def profile_budget_certificate():
    offsets = {
        1: {"A": 0, "B": 0, "E": 2},
        3: {"A": 0, "B": 0, "E": 2},
        5: {"A": 2, "B": 0, "E": 2},
        9: {"A": 2, "B": 4, "E": 2},
        11: {"A": 4, "B": 4, "E": 2},
        13: {"A": 4, "B": 4, "E": 2},
    }
    reports = []
    for residue in (1, 3, 5, 9, 11, 13):
        surviving = []
        for a_count in range(2, 7):
            for b_count in range(0, 8 - a_count):
                e_count = 7 - a_count - b_count
                if e_count < 1:
                    continue
                # The 4k terms sum to 28k=2(14k); offsets decide Omega.
                omega = (
                    a_count * offsets[residue]["A"]
                    + b_count * offsets[residue]["B"]
                    + e_count * offsets[residue]["E"]
                    - 2 * residue
                )
                fixed_cost = 2 * (e_count - 1)
                if omega >= 0 and fixed_cost <= omega:
                    surviving.append((a_count, b_count, e_count, omega))
        maximum = max((row[3] for row in surviving), default=None)
        reports.append((residue, maximum, tuple(surviving)))
    require(tuple((row[0], row[1]) for row in reports)
            == ((1, 8), (3, None), (5, 4), (9, 4), (11, 4), (13, 0)),
            ("profile budgets", reports))
    return tuple(reports)


def witness_control(prime: int, witness: tuple[int, ...]):
    masks = {coefficient: danger_mask(prime, coefficient) for coefficient in witness}
    multiplicities = []
    for sheet in range(2 * prime):
        multiplicities.append(sum(bool(mask & (1 << sheet)) for mask in masks.values()))
    require(all(multiplicities), ("positive witness misses", prime, witness))
    return (
        prime,
        witness,
        tuple(sorted(Counter(coefficient_type(c) for c in witness).items())),
        tuple(masks[c].bit_count() for c in witness),
        tuple(sorted(Counter(multiplicities).items())),
    )


def main() -> None:
    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
            "assert found")
    dependency_hashes = (
        lf_sha256(FINITE_PATH), lf_sha256(FINITE_OUTPUT),
        lf_sha256(TAIL_PATH), lf_sha256(TAIL_OUTPUT),
    )
    require(dependency_hashes == (
        EXPECTED_FINITE_SHA256, EXPECTED_FINITE_OUTPUT_SHA256,
        EXPECTED_TAIL_SHA256, EXPECTED_TAIL_OUTPUT_SHA256,
    ), ("dependency hash drift", dependency_hashes))

    cross_report = finite_cross_certificate()
    aa_report = finite_aa_certificate()
    plateau_report = plateau_certificate()
    tangent_report = tangent_certificate()

    # h=floor((p-1)/14); h^2>=2p is the only finite gate in the E/E
    # successive-minima argument, and the polynomial margin then increases.
    h607 = (607 - 1) // 14
    require(h607 * h607 >= 2 * 607, ("E/E boundary", h607))
    require(all(
        ((prime - 1) // 14) ** 2 >= 2 * prime
        for prime in range(607, 900)
    ), "E/E monotonic boundary control")

    profile_report = profile_budget_certificate()
    graph_report = mod240_graph_certificate()
    direct_controls = tuple(direct_graph_control(prime) for prime in (607, 1009))

    finite = runpy.run_path(str(FINITE_PATH))
    p601 = finite["decide_prime"](601)
    require(p601[0] is None, ("p=601 cover", p601[0]))
    p601_report = (p601[3].profiles, p601[3].nodes, p601[3].branches)
    require(p601_report == (15, 65, 50), ("p=601 tally", p601_report))
    positive_report = tuple(
        witness_control(prime, witness)
        for prime, witness in POSITIVE_CONTROLS.items()
    )

    semantic_payload = (
        dependency_hashes,
        cross_report,
        aa_report,
        plateau_report,
        tangent_report,
        h607,
        profile_report,
        graph_report,
        direct_controls,
        p601_report,
        positive_report,
    )
    semantic_sha256 = sha256(repr(semantic_payload).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
                ("semantic digest", semantic_sha256))

    print("THM-3445 prime-even literal half-twist cap-seven classification companion")
    print("status=VERIFIED-EXACT companion for RESERVED/PROVISIONAL theorem; audit required")
    print(f"dependency_hashes={dependency_hashes}")
    print(f"cross_finite_certificate={cross_report}")
    print(f"aa_orientation_finite_certificate={aa_report}")
    print(f"large_a_plateau_certificate={plateau_report}")
    print(f"height8_tangent_certificate={tangent_report}")
    print(f"EE_successive_minima_floor: p0=607 h={h607} h_squared={h607*h607}")
    print("EE_literal_intersection_lower_bound=10 for every prime p>=607")
    print(f"profile_overlap_budgets={tuple((row[0], row[1]) for row in profile_report)}")
    print("mod240_zero_graph: classes=64 degrees=(19,31,23) clique=6 "
          "rooted_presentations=(2,1,1) equality=three_complementary_pairs types=A4BE")
    print(f"mod240_first_prime_range=({min(row[1] for row in graph_report)},"
          f"{max(row[1] for row in graph_report)})")
    print(f"direct_graph_controls={direct_controls}")
    print(f"p601_exact_negative_profiles_nodes_branches={p601_report}")
    print(f"positive_controls={positive_report}")
    print("finite_stitch=p<=599 pinned FINITE-EXACT atlas; p=601 direct negative; tail p>=607")
    print(f"semantic_sha256={semantic_sha256}")
    print(f"script_lf_sha256={lf_sha256(source)}")


if __name__ == "__main__":
    main()
