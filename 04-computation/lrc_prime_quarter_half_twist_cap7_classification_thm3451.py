#!/usr/bin/env python3
"""Exact companion for the provisional prime-quarter cap-seven theorem.

The all-prime proof for Q=4p has an elementary infinite tail and three exact
finite gates.  This companion uses integer arithmetic only and checks:

* the four owner types, exact size offsets, central-fibre invoice, and every
  feasible seven-owner profile in the six prime classes modulo 14;
* the complete strengthened-target-free finite stitch through p=443, with a
  literal-mask DFS rooted at the mandatory E_e owner and no search-node cap;
* the ordered non-E_e/E_e line-and-coset compiler, all small and middle
  positive-weight bounds, the large-a plateau, and the E_e/E_e boundary;
* the raw-height-seven zero bank and every unit class modulo 420, including
  direct literal masks at two primes per class and the rooted clique-six
  equality orbit;
* all four positive profiles at p=17, the unique positive profile at p=37,
  and independent root-E_e hostile controls.

The theorem remains RESERVED / PROVISIONAL until a separate full-package
audit checks this proof and promotes its status.
"""

from __future__ import annotations

import ast
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from math import gcd, isqrt
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
THM3434_PATH = ROOT / "01-canon/theorems/THM-3434-seventeen-fibre-two-sided-mass-closure.md"
THM3435_PATH = ROOT / "01-canon/theorems/THM-3435-dyadic-fibre-grid-decomposition-for-literal-half-twists.md"
THM3445_PATH = ROOT / "01-canon/theorems/THM-3445-prime-even-half-twist-cap-seven-classification.md"
EXPECTED_DEPENDENCY_SHA256 = (
    "52a5a5caed75ad48ab35ce287c15f6cece074c88dddda2b87329792a3df70af7",
    "e2d202f97b1f39af032597767633c55f1abc4da261662824058b05b70724bb52",
    "06b28dd920ab525da087a6eba7bc3ecd25b51a878ce7e81f07725fb00e40b6bc",
)
EXPECTED_SEMANTIC_SHA256 = "b455ca48ce8454d754b1a7d4f816ecd76bb431a5b9202757152badd5e0993482"

TYPES = ("A", "B", "E_o", "E_e")
F = {"A": 0, "B": 1, "E_o": 2, "E_e": 2}
ALPHA = {"A": 8, "B": 4, "E_o": 2, "E_e": 2}
ETA = {"A": 1, "B": 1, "E_o": 1, "E_e": 0}
COEFFICIENT_RANGE = {"A": 4, "B": 2, "E_o": 1, "E_e": 1}
ROOTS = {"A": 1, "B": 2, "E_o": 4, "E_e": 8}

A0_BASES = (8, 9, 10, 11, 12, 13, 15, 23, 25, 29, 51)
POST_THM3445_BASES = A0_BASES + (14, 38)
LEGACY_EXCLUDED_PRIMES = (3, 5, 11, 13, 23, 29)
STRENGTHENED_EXCLUDED_PRIMES = (3, 5, 7, 11, 13, 19, 23, 29)

OFFSETS = {
    1: (0, 0, 0, 4),
    3: (2, 0, 0, 4),
    5: (2, 4, 0, 4),
    9: (6, 4, 8, 4),
    11: (6, 8, 8, 4),
    13: (8, 8, 8, 4),
}

SECTORS = (
    ("AA", "A", "A"),
    ("AB", "A", "B"),
    ("AE_o", "A", "E_o"),
    ("AE_e", "A", "E_e"),
    ("BB", "B", "B"),
    ("BE_o", "B", "E_o"),
    ("BE_e", "B", "E_e"),
    ("E_oE_o", "E_o", "E_o"),
    ("E_oE_e", "E_o", "E_e"),
)

EXPECTED_SMALL_ROWS = (
    ("AA", 2030, 2023, 7, 17668, 45040, 16),
    ("AB", 1020, 1015, 5, 7148, 17276, 16),
    ("AE_o", 515, 514, 1, 3694, 7388, 16),
    ("AE_e", 515, 512, 3, 3683, 7366, 16),
    ("BB", 1020, 1011, 9, 4646, 11228, 16),
    ("BE_o", 515, 510, 5, 2170, 4340, 16),
    ("BE_e", 515, 507, 8, 2165, 4330, 16),
    ("E_oE_o", 515, 499, 16, 1368, 2736, 16),
    ("E_oE_e", 515, 500, 15, 1415, 2830, 16),
)

EXPECTED_MIDDLE_ROWS = (
    ("AA", 9768, 16),
    ("AB", 4909, 20),
    ("AE_o", 2479, 22),
    ("AE_e", 2479, 22),
    ("BB", 4909, 16),
    ("BE_o", 2479, 20),
    ("BE_e", 2479, 20),
    ("E_oE_o", 2479, 16),
    ("E_oE_e", 2479, 16),
)

ROOT_A_ZERO_BANK = {
    "A": (
        (1, 1, 2), (1, 1, 4), (1, 3, 4), (1, 5, 4),
        (1, 5, 12), (3, 1, 4), (5, 1, 4),
    ),
    "B": ((1, 1, 2), (1, 3, 2), (1, 5, 2), (1, 5, 6), (3, 1, 2)),
    "E_o": ((1, 2, 1),),
    "E_e": ((1, 1, 1), (1, 2, 1), (1, 3, 1)),
}

POSITIVE_CONTROLS = (
    (17, (8, 11, 23, 24, 45, 56, 57), (4, 0, 0, 3), 8),
    (17, (1, 24, 32, 33, 35, 36, 67), (4, 0, 1, 2), 4),
    (17, (8, 11, 12, 23, 28, 45, 57), (4, 0, 2, 1), 0),
    (17, (8, 13, 21, 30, 38, 47, 55), (4, 2, 0, 1), 0),
    (37, (8, 33, 41, 100, 107, 115, 140), (4, 0, 2, 1), 8),
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def is_prime(value: int) -> bool:
    return value >= 2 and all(
        value % divisor for divisor in range(2, isqrt(value) + 1)
    )


def primes_through(limit: int) -> tuple[int, ...]:
    return tuple(value for value in range(3, limit + 1, 2) if is_prime(value))


def first_primes_in_class(residue: int, modulus: int, lower: int,
                          count: int = 2) -> tuple[int, ...]:
    values = []
    candidate = lower + (residue - lower) % modulus
    while len(values) < count:
        if is_prime(candidate):
            values.append(candidate)
        candidate += modulus
    return tuple(values)


def in_arc(word: int, ambient: int) -> bool:
    word %= ambient
    return 14 * min(word, ambient - word) < ambient


def danger_mask(prime: int, coefficient: int) -> int:
    period = 4 * prime
    ambient = 2 * period
    mask = 0
    for sheet in range(period):
        if in_arc(coefficient * (2 * sheet + 1), ambient):
            mask |= 1 << sheet
    return mask


def owner_type(coefficient: int) -> str:
    if coefficient % 2:
        return "A"
    if coefficient % 4:
        return "B"
    if (coefficient // 4) % 2:
        return "E_o"
    return "E_e"


def target_free_from_bases(prime: int, bases: tuple[int, ...]) -> bool:
    period = 4 * prime
    return not any(period % base == 0 for base in bases)


def exact_size(prime: int, type_name: str) -> int:
    quotient, residue = divmod(prime, 14)
    return 8 * quotient + OFFSETS[residue][TYPES.index(type_name)]


def profile_rows(residue: int, central_invoice: bool = True):
    rows = []
    offsets = OFFSETS[residue]
    for a_count in range(8):
        for b_count in range(8 - a_count):
            for eo_count in range(8 - a_count - b_count):
                ee_count = 7 - a_count - b_count - eo_count
                profile = (a_count, b_count, eo_count, ee_count)
                if a_count < 2 or ee_count < 1 or a_count + 2 * b_count < 4:
                    continue
                omega = sum(profile[index] * offsets[index] for index in range(4))
                omega -= 4 * residue
                if omega < 0:
                    continue
                if central_invoice and 4 * (ee_count - 1) > omega:
                    continue
                rows.append((profile, omega))
    return tuple(rows)


def profile_certificate():
    expected = {
        1: (26, 12, (((2, 1, 0, 4), 12),)),
        3: (10, 8, (((4, 0, 0, 3), 8),)),
        5: (5, 4, (((2, 3, 0, 2), 4), ((2, 4, 0, 1), 4))),
        9: (13, 8, (((2, 1, 3, 1), 8), ((4, 0, 2, 1), 8))),
        11: (
            10, 4,
            (((2, 1, 3, 1), 4), ((2, 2, 2, 1), 4),
             ((2, 3, 1, 1), 4), ((2, 4, 0, 1), 4)),
        ),
        13: (13, 0, None),
    }
    report = []
    for residue in (1, 3, 5, 9, 11, 13):
        rows = profile_rows(residue)
        raw_rows = profile_rows(residue, central_invoice=False)
        maximum = max(omega for _profile, omega in rows)
        maximizers = tuple(row for row in rows if row[1] == maximum)
        want_count, want_maximum, want_maximizers = expected[residue]
        require(len(rows) == want_count, ("profile count", residue, len(rows)))
        require(maximum == want_maximum, ("profile maximum", residue, maximum))
        if want_maximizers is not None:
            require(maximizers == want_maximizers,
                    ("profile maximizers", residue, maximizers))
        expected_raw_count = {1: 26, 3: 19, 5: 13, 9: 19, 11: 13, 13: 13}[residue]
        require(len(raw_rows) == expected_raw_count,
                ("raw profile count", residue, len(raw_rows)))
        report.append((
            residue, OFFSETS[residue], len(raw_rows), len(rows), maximum, maximizers,
        ))
    return tuple(report)


@dataclass
class SearchStats:
    states: int = 0
    branches: int = 0
    budget_prunes: int = 0
    capacity_prunes: int = 0
    reach_prunes: int = 0


def search_profile(period: int, profile: tuple[int, int, int, int],
                   omega: int, residues: tuple[int, ...],
                   types: tuple[str, ...], masks: tuple[int, ...],
                   stats: SearchStats):
    """Complete typed cover search, normalized at the E_e owner r=8."""
    full = (1 << period) - 1
    remaining = dict(zip(TYPES, profile))
    remaining["E_e"] -= 1
    if remaining["E_e"] < 0:
        return None

    by_type = {
        type_name: tuple(
            index for index in range(1, len(residues))
            if types[index] == type_name
        )
        for type_name in TYPES
    }
    if any(len(by_type[type_name]) < remaining[type_name] for type_name in TYPES):
        return None

    positions = {
        type_name: {index: position for position, index in enumerate(by_type[type_name])}
        for type_name in TYPES
    }
    lower = {type_name: 0 for type_name in TYPES}
    type_sizes = {type_name: exact_size(period // 4, type_name) for type_name in TYPES}
    chosen = [0]

    def recurse(current_mask: int, overlap_used: int):
        stats.states += 1
        slots = sum(remaining.values())
        if slots == 0:
            if current_mask == full:
                require(overlap_used == omega,
                        ("cover excess identity", period, profile, overlap_used, omega))
                return tuple(residues[index] for index in chosen)
            return None

        budget_left = omega - overlap_used
        if budget_left < 0:
            stats.budget_prunes += 1
            return None
        remaining_mass = sum(
            remaining[type_name] * type_sizes[type_name]
            for type_name in TYPES
        )
        if current_mask.bit_count() + remaining_mass < period:
            stats.capacity_prunes += 1
            return None

        reachable = current_mask
        options_by_type = {}
        for type_name in TYPES:
            if not remaining[type_name]:
                continue
            options = []
            for index in by_type[type_name]:
                if positions[type_name][index] < lower[type_name]:
                    continue
                increment = (masks[index] & current_mask).bit_count()
                if increment <= budget_left:
                    options.append((index, increment))
                    reachable |= masks[index]
            if len(options) < remaining[type_name]:
                stats.capacity_prunes += 1
                return None
            options_by_type[type_name] = tuple(options)
        if reachable != full:
            stats.reach_prunes += 1
            return None

        type_name = min(
            options_by_type,
            key=lambda name: (len(options_by_type[name]), -remaining[name], name),
        )
        options = options_by_type[type_name]
        need = remaining[type_name]
        for option_position, (index, increment) in enumerate(options):
            if len(options) - option_position < need:
                break
            stats.branches += 1
            chosen.append(index)
            remaining[type_name] -= 1
            old_lower = lower[type_name]
            lower[type_name] = positions[type_name][index] + 1
            result = recurse(current_mask | masks[index], overlap_used + increment)
            lower[type_name] = old_lower
            remaining[type_name] += 1
            chosen.pop()
            if result is not None:
                return result
        return None

    return recurse(masks[0], 0)


def finite_prime_row(prime: int):
    period = 4 * prime
    root = ROOTS["E_e"]
    root_mask = danger_mask(prime, root)
    sizes = {type_name: exact_size(prime, type_name) for type_name in TYPES}
    for type_name in TYPES:
        direct_size = danger_mask(prime, ROOTS[type_name]).bit_count()
        require(direct_size == sizes[type_name],
                ("type size", prime, type_name, direct_size, sizes[type_name]))

    profiles = profile_rows(prime % 14)
    maximum_omega = max(omega for _profile, omega in profiles)
    root_sheets = tuple(
        sheet for sheet in range(period) if root_mask & (1 << sheet)
    )
    low = []
    for coefficient in range(1, period):
        if coefficient in (prime, 2 * prime, 3 * prime, root):
            continue
        weight = sum(
            in_arc(coefficient * (2 * sheet + 1), 8 * prime)
            for sheet in root_sheets
        )
        if weight <= maximum_omega:
            low.append(coefficient)
    residues = (root,) + tuple(low)
    types = tuple(owner_type(coefficient) for coefficient in residues)
    masks = tuple(danger_mask(prime, coefficient) for coefficient in residues)
    for coefficient, type_name, mask in zip(residues, types, masks):
        require(mask.bit_count() == sizes[type_name],
                ("candidate size", prime, coefficient, type_name, mask.bit_count()))

    stats = SearchStats()
    supported_profiles = []
    witnesses = []
    for profile, omega in sorted(profiles, key=lambda row: (row[1], row[0])):
        witness = search_profile(
            period, profile, omega, residues, types, masks, stats
        )
        if witness is not None:
            supported_profiles.append((profile, omega))
            witnesses.append(witness)
    return (
        len(profiles),
        tuple(supported_profiles),
        tuple(witnesses),
        len(residues),
        stats,
    )


def finite_stitch_certificate():
    all_primes = primes_through(443)
    legacy_excluded = tuple(
        prime for prime in all_primes
        if not target_free_from_bases(prime, A0_BASES)
    )
    strengthened_excluded = tuple(
        prime for prime in all_primes
        if not target_free_from_bases(prime, POST_THM3445_BASES)
    )
    require(legacy_excluded == LEGACY_EXCLUDED_PRIMES,
            ("legacy exclusions", legacy_excluded))
    require(strengthened_excluded == STRENGTHENED_EXCLUDED_PRIMES,
            ("strengthened exclusions", strengthened_excluded))
    primes = tuple(
        prime for prime in all_primes
        if target_free_from_bases(prime, POST_THM3445_BASES)
    )

    total_raw_profiles = 0
    total_profiles = 0
    total_stats = SearchStats()
    support = []
    root_bank_range = []
    for prime in primes:
        total_raw_profiles += len(profile_rows(prime % 14, central_invoice=False))
        profile_count, supported_profiles, witnesses, root_bank, stats = finite_prime_row(prime)
        total_profiles += profile_count
        root_bank_range.append(root_bank)
        for field in SearchStats.__dataclass_fields__:
            setattr(total_stats, field, getattr(total_stats, field) + getattr(stats, field))
        if supported_profiles:
            support.append((prime, supported_profiles, witnesses))

    support_profile_summary = tuple(
        (prime, supported_profiles) for prime, supported_profiles, _witnesses in support
    )
    require(len(primes) == 77, ("finite prime count", len(primes)))
    require(total_raw_profiles == 1312,
            ("finite raw profile count", total_raw_profiles))
    require(total_profiles == 959, ("finite profile count", total_profiles))
    stats_tuple = tuple(
        getattr(total_stats, field) for field in SearchStats.__dataclass_fields__
    )
    require(stats_tuple == (14839, 13966, 0, 12054, 1220),
            ("finite DFS tally", stats_tuple))
    require((min(root_bank_range), max(root_bank_range)) == (36, 433),
            ("finite root-bank range", min(root_bank_range), max(root_bank_range)))
    require(support_profile_summary == (
        (17, (
            ((4, 0, 2, 1), 0),
            ((4, 2, 0, 1), 0),
            ((4, 0, 1, 2), 4),
            ((4, 0, 0, 3), 8),
        )),
        (37, (((4, 0, 2, 1), 8),)),
    ), ("finite support profiles", support_profile_summary))
    require(primes[-1] == 443 and not is_prime(445) and not is_prime(447)
            and is_prime(449), "finite/tail boundary")
    return (
        len(primes), total_raw_profiles, total_profiles,
        stats_tuple,
        (min(root_bank_range), max(root_bank_range)),
        support_profile_summary,
        primes[-1], 449,
    )


def section_shape(first_type: str, second_type: str,
                  a_value: int, b_value: int, line: int) -> int:
    left = ALPHA[first_type] * a_value
    right = ALPHA[second_type] * b_value
    return max(
        0,
        min(left, 14 * line + right) - max(-left, 14 * line - right),
    )


@lru_cache(maxsize=None)
def odd_class_count(second_type: str, b_value: int,
                    lift: int, line_residue: int) -> int:
    congruence_modulus = b_value * ALPHA[second_type]
    word_period = 2 * congruence_modulus
    return sum(
        word % 2 == 1
        and (line_residue - lift * word) % congruence_modulus == 0
        for word in range(word_period)
    )


def admissible_lifts(second_type: str, a_value: int, b_value: int):
    return tuple(
        lift
        for lift in range(COEFFICIENT_RANGE[second_type] * b_value + 1)
        if gcd(lift, b_value) == 1
        and (lift - (b_value * ETA[second_type] - a_value)) % 2 == 0
    )


def compiler_row(prime_floor: int, first_type: str, second_type: str,
                 a_value: int, b_value: int, lift: int):
    diameter = ALPHA[first_type] * a_value + ALPHA[second_type] * b_value
    line_count = 0
    residue_class_count = 0
    lower_bound = 0
    for line in range(-(diameter // 14 + 2), diameter // 14 + 3):
        shape = section_shape(first_type, second_type, a_value, b_value, line)
        if not shape:
            continue
        classes = odd_class_count(
            second_type, b_value, lift,
            line % (b_value * ALPHA[second_type]),
        )
        if not classes:
            continue
        line_count += 1
        residue_class_count += classes
        denominator = 28 * a_value * b_value * ALPHA[second_type]
        lower_bound += classes * max(0, (prime_floor * shape - 1) // denominator)
    lower_bound *= 1 << F[first_type]
    return line_count, residue_class_count, lower_bound


def line_compiler_certificate():
    small_report = []
    middle_report = []
    zero_rows = {}
    for sector, first_type, second_type in SECTORS:
        total = 0
        active = 0
        no_line = 0
        lines = 0
        classes = 0
        bounds = []
        zeros = []
        for a_value in range(1, 21):
            for b_value in range(1, 14):
                if gcd(a_value, b_value) != 1:
                    continue
                for lift in admissible_lifts(second_type, a_value, b_value):
                    total += 1
                    line_count, class_count, lower_bound = compiler_row(
                        449, first_type, second_type, a_value, b_value, lift
                    )
                    if not line_count:
                        no_line += 1
                        zeros.append((a_value, b_value, lift))
                    else:
                        active += 1
                        lines += line_count
                        classes += class_count
                        bounds.append(lower_bound)
        zero_rows[sector] = tuple(zeros)
        small_report.append(
            (sector, total, active, no_line, lines, classes, min(bounds))
        )

        middle_bounds = []
        middle_total = 0
        for a_value in range(21, 120):
            for b_value in range(1, 14):
                if gcd(a_value, b_value) != 1:
                    continue
                for lift in admissible_lifts(second_type, a_value, b_value):
                    middle_total += 1
                    line_count, _class_count, lower_bound = compiler_row(
                        max(383, 14 * a_value + 1),
                        first_type, second_type, a_value, b_value, lift,
                    )
                    require(line_count > 0,
                            ("middle no-line", sector, a_value, b_value, lift))
                    middle_bounds.append(lower_bound)
        middle_report.append((sector, middle_total, min(middle_bounds)))

    small_report = tuple(small_report)
    middle_report = tuple(middle_report)
    require(small_report == EXPECTED_SMALL_ROWS,
            ("small compiler table", small_report))
    require(middle_report == EXPECTED_MIDDLE_ROWS,
            ("middle compiler table", middle_report))
    require(tuple(zero_rows[sector] for sector in ("AA", "AB", "AE_o", "AE_e"))
            == tuple(ROOT_A_ZERO_BANK[type_name] for type_name in TYPES),
            ("root A zero bank", zero_rows))
    require(max(
        a_value + b_value
        for rows in zero_rows.values()
        for a_value, b_value, _lift in rows
    ) <= 7, ("raw zero height", zero_rows))
    return small_report, middle_report, tuple(
        (sector, zero_rows[sector]) for sector, _first, _second in SECTORS
    )


def admissible_line_residues(second_type: str, b_value: int, lift: int):
    modulus = b_value * ALPHA[second_type]
    return tuple(sorted({
        lift * word % modulus
        for word in range(2 * modulus) if word % 2
    }))


def cyclic_gap(residues: tuple[int, ...], modulus: int) -> int:
    require(residues, ("empty residue cycle", modulus))
    if len(residues) == 1:
        return modulus
    return max(
        (residues[(index + 1) % len(residues)] - residues[index]) % modulus
        for index in range(len(residues))
    )


def plateau_certificate():
    reports = []
    for sector, first_type, second_type in SECTORS:
        maximum_gap = 0
        minimum_lines = None
        for parity in (0, 1):
            a_value = 120 + parity
            for b_value in range(1, 14):
                for lift in admissible_lifts(second_type, a_value, b_value):
                    modulus = b_value * ALPHA[second_type]
                    residues = admissible_line_residues(second_type, b_value, lift)
                    gap = cyclic_gap(residues, modulus)
                    maximum_gap = max(maximum_gap, gap)
                    require(gap <= ALPHA[second_type],
                            ("plateau residue gap", sector, parity, b_value, lift, gap))
                    radius = (
                        ALPHA[first_type] * a_value
                        - ALPHA[second_type] * b_value
                    ) // 14
                    count = sum(
                        line % modulus in residues
                        for line in range(-radius, radius + 1)
                    )
                    minimum_lines = count if minimum_lines is None else min(minimum_lines, count)
        require(minimum_lines is not None and minimum_lines >= 15,
                ("plateau line floor", sector, minimum_lines))
        reports.append((sector, maximum_gap, minimum_lines,
                        minimum_lines * (1 << F[first_type])))
    return tuple(reports)


def ee_successive_minima_certificate():
    boundary = 409
    h_boundary = (boundary - 1) // 14
    require(h_boundary * h_boundary >= 2 * boundary,
            ("E_e/E_e boundary", boundary, h_boundary))
    require(all(
        ((prime - 1) // 14) ** 2 >= 2 * prime
        for prime in range(boundary, 1000)
    ), "E_e/E_e monotone boundary control")
    return boundary, h_boundary, h_boundary * h_boundary, 20


def candidate_coefficients(prime: int):
    labelled = [(('A', 0, 0, 0, 0), 1)]
    for second_type in TYPES:
        for a_value, b_value, lift in ROOT_A_ZERO_BANK[second_type]:
            for sign in (-1, 1):
                numerator = sign * a_value + lift * prime
                if numerator % b_value:
                    continue
                effective = numerator // b_value
                if not 0 < effective < COEFFICIENT_RANGE[second_type] * prime:
                    continue
                if effective % 2 != ETA[second_type]:
                    continue
                coefficient = (1 << F[second_type]) * effective
                labelled.append(
                    ((second_type, a_value, b_value, lift, sign), coefficient)
                )
    require(len(labelled) == 16 and len({value for _label, value in labelled}) == 16,
            ("candidate bank", prime, labelled))
    require(Counter(owner_type(value) for _label, value in labelled[1:])
            == Counter({"A": 7, "B": 4, "E_o": 1, "E_e": 3}),
            ("candidate types", prime, labelled))
    return tuple(labelled)


def candidate_graph(prime: int):
    labelled = candidate_coefficients(prime)
    masks = {label: danger_mask(prime, coefficient) for label, coefficient in labelled}
    labels = tuple(label for label, _coefficient in labelled)
    coefficients = dict(labelled)
    root = labels[0]
    require(all(not (masks[root] & masks[label]) for label in labels[1:]),
            ("candidate is not a root zero", prime))
    edges = {
        (first, second): not (masks[first] & masks[second])
        for first, second in combinations(labels, 2)
    }

    best_size = 1
    best = set()

    def adjacent(first, second):
        if labels.index(first) < labels.index(second):
            return edges[(first, second)]
        return edges[(second, first)]

    def extend(chosen, candidates):
        nonlocal best_size, best
        if len(chosen) + len(candidates) < best_size:
            return
        if len(chosen) > best_size:
            best_size = len(chosen)
            best = {tuple(chosen)}
        elif len(chosen) == best_size:
            best.add(tuple(chosen))
        remaining = candidates
        while remaining:
            if len(chosen) + len(remaining) < best_size:
                return
            vertex = remaining[0]
            tail = remaining[1:]
            extend(
                chosen + (vertex,),
                tuple(other for other in tail if adjacent(vertex, other)),
            )
            remaining = tail

    extend((root,), labels[1:])
    maximum_coefficients = tuple(sorted({
        tuple(sorted(coefficients[label] for label in clique))
        for clique in best if len(clique) == best_size
    }))
    signature = (
        labels,
        tuple(owner_type(coefficients[label]) for label in labels),
        tuple(edges[pair] for pair in combinations(labels, 2)),
    )
    return signature, best_size, maximum_coefficients


def mod420_graph_certificate():
    residues = tuple(residue for residue in range(420) if gcd(residue, 420) == 1)
    require(len(residues) == 96, ("mod-420 units", len(residues)))
    reports = []
    signatures = Counter()
    for residue in residues:
        first_prime, second_prime = first_primes_in_class(residue, 420, 449)
        first_signature, clique_size, cliques = candidate_graph(first_prime)
        second_signature, second_clique_size, second_cliques = candidate_graph(second_prime)
        require(first_signature == second_signature,
                ("mod-420 graph drift", residue, first_prime, second_prime))
        require(clique_size == second_clique_size == 6,
                ("rooted clique", residue, clique_size, second_clique_size))
        expected = tuple(sorted((
            1, 4 * first_prime - 1,
            2 * first_prime - 2, 2 * first_prime - 1,
            2 * first_prime + 1, 2 * first_prime + 2,
        )))
        expected_second = tuple(sorted((
            1, 4 * second_prime - 1,
            2 * second_prime - 2, 2 * second_prime - 1,
            2 * second_prime + 1, 2 * second_prime + 2,
        )))
        require(cliques == (expected,),
                ("rooted equality", residue, first_prime, cliques))
        require(second_cliques == (expected_second,),
                ("second rooted equality", residue, second_prime, second_cliques))
        signatures[first_signature] += 1
        reports.append((residue, first_prime, second_prime))
    require(len(signatures) == 16 and set(signatures.values()) == {6},
            ("graph signatures", len(signatures), signatures.values()))
    return (
        len(residues), len(signatures), tuple(sorted(signatures.values())),
        (min(row[1] for row in reports), max(row[1] for row in reports)),
        (min(row[2] for row in reports), max(row[2] for row in reports)),
        15, (7, 4, 1, 3), 6,
    )


def full_root_control(prime: int, root: int):
    period = 4 * prime
    root_mask = danger_mask(prime, root)
    root_sheets = tuple(
        sheet for sheet in range(period) if root_mask & (1 << sheet)
    )
    weights = {}
    for coefficient in range(1, period):
        if coefficient in (prime, 2 * prime, 3 * prime, root):
            continue
        weights[coefficient] = sum(
            in_arc(coefficient * (2 * sheet + 1), 8 * prime)
            for sheet in root_sheets
        )
    zeros = tuple(sorted(coefficient for coefficient, weight in weights.items() if not weight))
    minimum_positive = min(weight for weight in weights.values() if weight)
    if root == 1:
        expected = tuple(sorted(value for _label, value in candidate_coefficients(prime)[1:]))
        require(zeros == expected, ("full root-A bank", prime, zeros, expected))
        return prime, root, len(zeros), tuple(
            Counter(owner_type(value) for value in zeros)[type_name] for type_name in TYPES
        ), minimum_positive

    masks = {coefficient: danger_mask(prime, coefficient) for coefficient in (root,) + zeros}
    neighbor_edges = sum(
        not (masks[first] & masks[second])
        for first, second in combinations(zeros, 2)
    )
    nodes = (root,) + zeros
    adjacency = {
        (first, second): not (masks[first] & masks[second])
        for first, second in combinations(nodes, 2)
    }

    def adjacent(first, second):
        return adjacency[(first, second)] if first < second else adjacency[(second, first)]

    best_size = 1
    best = set()

    def extend(chosen, candidates):
        nonlocal best_size, best
        if len(chosen) + len(candidates) < best_size:
            return
        if len(chosen) > best_size:
            best_size = len(chosen)
            best = {tuple(chosen)}
        elif len(chosen) == best_size:
            best.add(tuple(chosen))
        remaining = candidates
        while remaining:
            vertex = remaining[0]
            tail = remaining[1:]
            extend(chosen + (vertex,), tuple(
                other for other in tail if adjacent(vertex, other)
            ))
            remaining = tail

    extend((root,), zeros)
    cliques = tuple(sorted({tuple(sorted(clique)) for clique in best if len(clique) == best_size}))
    expected_cliques = tuple(sorted((
        tuple(sorted((root, prime - 4, prime + 4,
                      3 * prime - 4, 3 * prime + 4, 4 * prime - 8))),
        tuple(sorted((root, 2 * prime - 16, 2 * prime - 8,
                      2 * prime + 8, 2 * prime + 16, 4 * prime - 8))),
    )))
    require(len(zeros) == 35, ("root-E_e degree", prime, len(zeros)))
    require(tuple(Counter(owner_type(value) for value in zeros)[type_name] for type_name in TYPES)
            == (12, 12, 11, 0), ("root-E_e types", prime, zeros))
    require(neighbor_edges == 109, ("root-E_e edges", prime, neighbor_edges))
    require(best_size == 6 and cliques == expected_cliques,
            ("root-E_e cliques", prime, best_size, cliques))
    return (
        prime, root, len(zeros), (12, 12, 11, 0), minimum_positive,
        neighbor_edges, best_size, cliques,
    )


def witness_control(prime: int, witness: tuple[int, ...],
                    expected_profile: tuple[int, int, int, int],
                    expected_omega: int):
    period = 4 * prime
    require(len(witness) == len(set(witness)) == 7,
            ("witness owners", prime, witness))
    masks = tuple(danger_mask(prime, coefficient) for coefficient in witness)
    union = 0
    for mask in masks:
        union |= mask
    require(union == (1 << period) - 1, ("witness misses", prime, witness))
    profile = tuple(Counter(owner_type(value) for value in witness)[name] for name in TYPES)
    omega = sum(mask.bit_count() for mask in masks) - period
    require(profile == expected_profile and omega == expected_omega,
            ("witness profile", prime, witness, profile, omega))
    multiplicities = Counter(
        sum(bool(mask & (1 << sheet)) for mask in masks)
        for sheet in range(period)
    )
    positive_pairs = []
    for first_index, second_index in combinations(range(7), 2):
        intersection = masks[first_index] & masks[second_index]
        if intersection:
            sheets = tuple(
                sheet for sheet in range(period) if intersection & (1 << sheet)
            )
            positive_pairs.append((
                witness[first_index], witness[second_index], len(sheets), sheets,
            ))
    return (
        prime, witness, profile, tuple(mask.bit_count() for mask in masks), omega,
        tuple(sorted(multiplicities.items())), tuple(positive_pairs),
    )


def tiny_period_certificate():
    reports = []
    for period in (1, 2, 4):
        union = 0
        sizes = []
        for coefficient in range(1, period):
            mask = 0
            for sheet in range(period):
                if in_arc(coefficient * (2 * sheet + 1), 2 * period):
                    mask |= 1 << sheet
            union |= mask
            sizes.append(mask.bit_count())
        require(union != (1 << period) - 1,
                ("tiny proper-period cover", period, sizes))
        reports.append((period, tuple(sizes), union.bit_count()))
    return tuple(reports)


def security_certificate(source: Path):
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
            "assert found")
    forbidden_calls = {"eval", "exec", "compile", "__import__"}
    called = {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    require(not (called & forbidden_calls), ("forbidden calls", called & forbidden_calls))
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split(".")[0])
    allowed = {
        "__future__", "ast", "collections", "dataclasses", "functools",
        "hashlib", "itertools", "math", "pathlib",
    }
    require(imported <= allowed, ("unexpected imports", imported - allowed))
    return tuple(sorted(imported)), tuple(sorted(forbidden_calls))


def main() -> None:
    source = Path(__file__)
    security_report = security_certificate(source)
    dependency_hashes = tuple(
        lf_sha256(path) for path in (THM3434_PATH, THM3435_PATH, THM3445_PATH)
    )
    require(dependency_hashes == EXPECTED_DEPENDENCY_SHA256,
            ("dependency hash drift", dependency_hashes))

    profile_report = profile_certificate()
    finite_report = finite_stitch_certificate()
    small_report, middle_report, zero_report = line_compiler_certificate()
    plateau_report = plateau_certificate()
    ee_report = ee_successive_minima_certificate()
    graph_report = mod420_graph_certificate()
    direct_root_a = tuple(full_root_control(prime, 1) for prime in (449, 1009))
    direct_root_ee = tuple(full_root_control(prime, 8) for prime in (409, 449, 1009))
    positive_report = tuple(
        witness_control(prime, witness, profile, omega)
        for prime, witness, profile, omega in POSITIVE_CONTROLS
    )
    tiny_report = tiny_period_certificate()

    semantic_payload = (
        dependency_hashes, security_report, profile_report, finite_report,
        small_report, middle_report, zero_report, plateau_report, ee_report,
        graph_report, direct_root_a, direct_root_ee, positive_report, tiny_report,
    )
    semantic_sha256 = sha256(repr(semantic_payload).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
                ("semantic digest", semantic_sha256))

    print("THM-3451 prime-quarter literal half-twist cap-seven classification companion")
    print("status=VERIFIED-EXACT companion for RESERVED/PROVISIONAL theorem; independent audit required")
    print(f"dependency_hashes={dependency_hashes}")
    print(f"security_certificate={security_report}")
    print("target_free_universes="
          f"legacy_A0_excludes={LEGACY_EXCLUDED_PRIMES}; "
          f"post_THM3445_excludes={STRENGTHENED_EXCLUDED_PRIMES}")
    print("type_size_offsets="
          f"{tuple((residue, OFFSETS[residue]) for residue in OFFSETS)}")
    print("profile_budgets=(pmod14,raw_mass,central_admissible,Omega_max)="
          f"{tuple((row[0], row[2], row[3], row[4]) for row in profile_report)}")
    print("finite_stitch=(prime_count,raw_profiles,central_DFS_profiles,"
          f"DFS_tally,root_bank_range,support,last_prime,tail_prime)={finite_report}")
    print(f"line_compiler_small={small_report}")
    print(f"line_compiler_middle={middle_report}")
    print("line_compiler_table="
          "AA:(8b,1),AB:(4b,1),AE:(2b,1),BB:(4b,2),"
          "BE:(2b,2),EoE:(2b,4); entries=(J modulus,literal factor)")
    print("raw_zero_atlas="
          f"rows={sum(len(rows) for _sector, rows in zero_report)} "
          f"sector_counts={tuple((sector, len(rows)) for sector, rows in zero_report)} "
          "maximum_raw_height=7")
    print(f"large_a_plateau={plateau_report}")
    print(f"EeEe_successive_minima={ee_report}")
    print(f"mod420_zero_graph={graph_report}")
    print(f"direct_root_A_controls={direct_root_a}")
    print(f"direct_root_Ee_controls={direct_root_ee}")
    print(f"positive_profile_controls={positive_report}")
    print(f"tiny_proper_period_controls={tiny_report}")
    print("tail=p>=449 every pair is zero or has literal weight >=15; rooted zero clique number=6")
    print("classification=post-THM3445 target-free support iff p in (17,37)")
    print(f"semantic_sha256={semantic_sha256}")
    print(f"script_lf_sha256={lf_sha256(source)}")


if __name__ == "__main__":
    main()
