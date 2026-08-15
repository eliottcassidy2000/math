#!/usr/bin/env python3
"""Exact companion for THM-3434's seventeen-fibre mass closure.

The all-tower exclusions are elementary integer inequalities.  This program
checks every local size identity used by those inequalities, gives a separate
complete set-cover proof at the exceptional modulus 85, and verifies the
period descent and positive pullbacks used in the global classification.
"""

from functools import lru_cache
from hashlib import sha256
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PINNED = (
    (
        "THM-3416",
        ROOT / "01-canon/theorems/THM-3416-zero-mode-cochain-global-rank-six-support.md",
        "42a9309145de51d1bb6fca0b7c1945302ff37a63a3183e1dfed838c07118e8bf",
    ),
    (
        "THM-3421",
        ROOT / "01-canon/theorems/THM-3421-prime-half-twist-rank-seven-classification.md",
        "2f577354a06628660d90f70aca34378a186b7e67a55795e6c4a5c32a255d9736",
    ),
    (
        "THM-3429",
        ROOT / "01-canon/theorems/THM-3429-prime-fibre-activity-descent-for-mixed-order-half-twist-seven-covers.md",
        "58ebf850fc79fc9afed57966b7599e7376a6684fa3bbc5a2aa2e1a8e6e0ca148",
    ),
    (
        "THM-3432",
        ROOT / "01-canon/theorems/THM-3432-order-two-fixed-half-parity-transplant.md",
        "2d3d0e6b59d9f8f7d2d4da6d7e23c78660634b54b5f6b807540ea6c82490135b",
    ),
)
EXPECTED_SEMANTIC_DIGEST = "13984be8cee76a72b5d1ff98011538e56eb9c5f86ff3db41ed80f7c536189f2d"
SUPPORT_DIVISORS = (9, 11, 13, 15, 23, 25, 29, 51)


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def half_mask(modulus, residue):
    answer = 0
    double = 2 * modulus
    for sheet in range(modulus):
        word = residue * (2 * sheet + 1) % double
        if 14 * min(word, double - word) < double:
            answer |= 1 << sheet
    return answer


def quotient_order(modulus, residue):
    return modulus // gcd(modulus, residue)


def ceil_seventh(value):
    return (value + 6) // 7


def pullback_mask(mask, small, large):
    require(large % small == 0, (small, large))
    return sum(1 << sheet for sheet in range(large) if mask >> (sheet % small) & 1)


def exact_size_and_density_audit():
    rows = []
    cells = 0
    for exponent in range(1, 9):
        order = 17 ** exponent
        expected_bound = 3 * order // 17
        # Exact parity size formulas depend only on order modulo 14.
        k, residue_class = divmod(order, 14)
        even_size = 2 * k + 1
        odd_size = 2 * k + (2 if residue_class in (9, 11, 13) else 0)
        maximum = max(even_size, odd_size)
        require(maximum == ceil_seventh(order), (order, maximum))
        require(maximum <= expected_bound, (order, maximum, expected_bound))
        rows.append(("pure", exponent, order, even_size, odd_size, expected_bound))
        cells += 2

    for exponent in range(0, 8):
        order = 5 * (17 ** exponent)
        k, residue_class = divmod(order, 14)
        even_size = 2 * k + 1
        odd_size = 2 * k + (2 if residue_class in (9, 11, 13) else 0)
        maximum = max(even_size, odd_size)
        require(maximum <= order // 5, (order, maximum, order // 5))
        rows.append(("five", exponent, order, even_size, odd_size, order // 5))
        cells += 2
    return tuple(rows), cells


def fibre_section_audit():
    rows = []
    cells = 0
    for base_modulus in (3, 5, 17, 19):
        modulus = 17 * base_modulus
        minimum = 17
        maximum = 0
        active_count = 0
        for residue in range(1, modulus + 1):
            if residue % 17 == 0:
                continue
            mask = half_mask(modulus, residue)
            active_count += 1
            for base_sheet in range(base_modulus):
                size = sum(
                    bool(mask >> (base_sheet + base_modulus * point) & 1)
                    for point in range(17)
                )
                require(size in (2, 3),
                        ("seventeen-section", modulus, residue, base_sheet, size))
                minimum = min(minimum, size)
                maximum = max(maximum, size)
                cells += 17
        require((minimum, maximum) == (2, 3), (modulus, minimum, maximum))
        rows.append((modulus, active_count, minimum, maximum))

    # The exceptional Q=85 proof also uses the sharp five-fibre upper bound.
    modulus = 85
    maximum = 0
    for residue in range(1, modulus + 1):
        if residue % 5 == 0:
            continue
        mask = half_mask(modulus, residue)
        for base_sheet in range(17):
            size = sum(
                bool(mask >> (base_sheet + 17 * point) & 1)
                for point in range(5)
            )
            require(size <= 1, ("five-section", residue, base_sheet, size))
            maximum = max(maximum, size)
            cells += 5
    require(maximum == 1, maximum)
    rows.append((85, 68, 0, maximum))
    return tuple(rows), cells


def tower_inequality_audit():
    pure_rows = []
    five_rows = []
    for exponent in range(2, 11):
        base = 17 ** (exponent - 1)
        inactive_upper = 3 * base // 17
        lower = 17 * base - 5 * inactive_upper
        upper = 6 * ceil_seventh(17 * base)
        require(lower > upper, ("pure tower", exponent, lower, upper))
        require(184 * base > 612, ("pure symbolic", exponent, base))
        pure_rows.append((exponent, base, inactive_upper, lower, upper, lower - upper))

        multiplier = base
        inactive_upper = multiplier
        lower = 85 * multiplier - 5 * inactive_upper
        upper = 30 * ceil_seventh(17 * multiplier)
        require(lower > upper, ("five tower", exponent, lower, upper))
        require(50 * multiplier > 180, ("five symbolic", exponent, multiplier))
        five_rows.append(
            (exponent, multiplier, inactive_upper, lower, upper, lower - upper)
        )
    return tuple(pure_rows), tuple(five_rows)


def q85_invoice():
    modulus = 85
    prime = 17
    base = 5
    inactive = tuple(
        residue
        for residue in range(1, modulus + 1)
        if residue % prime == 0 and residue % 2 == 0
    )
    require(inactive == (34, 68), inactive)
    inactive_masks = {half_mask(modulus, residue) for residue in inactive}
    require(len(inactive_masks) == 1, inactive_masks)
    inactive_mask = next(iter(inactive_masks))
    covered_bases = tuple(
        sheet
        for sheet in range(base)
        if all(inactive_mask >> (sheet + base * point) & 1 for point in range(prime))
    )
    require(covered_bases == (2,), covered_bases)

    contribution_rows = []
    for residue in range(1, modulus + 1):
        if residue % prime == 0:
            continue
        mask = half_mask(modulus, residue)
        fixed_section = sum(
            bool(mask >> (2 + base * point) & 1) for point in range(prime)
        )
        missed_contribution = mask.bit_count() - fixed_section
        if residue % 5:
            require(missed_contribution == 10,
                    ("Q85 full-order invoice", residue, missed_contribution))
            kind = "full"
        elif residue % 2 == 0:
            require(missed_contribution == 12,
                    ("Q85 even order17 invoice", residue, missed_contribution))
            kind = "order17-even"
        else:
            require(missed_contribution == 8,
                    ("Q85 odd order17 invoice", residue, missed_contribution))
            kind = "order17-odd"
        contribution_rows.append(
            (residue, kind, mask.bit_count(), fixed_section, missed_contribution)
        )

    # If b of the six 17-active owners are divisible by five, the missed
    # five-fibre argument gives b<=2.  The four missed seventeen-fibres then
    # receive at most 60+2b incidences, below the required 68.
    inequalities = []
    for count_order17 in range(3):
        available = 60 + 2 * count_order17
        require(available < 68, (count_order17, available))
        inequalities.append((count_order17, available, 68 - available))
    return inactive, covered_bases, tuple(contribution_rows), tuple(inequalities)


def exact_cover_search(modulus, cap):
    raw = tuple(half_mask(modulus, residue) for residue in range(1, modulus + 1))
    masks = tuple(sorted(set(raw) - {0}, key=lambda item: (-item.bit_count(), item)))
    full = (1 << modulus) - 1
    coverers = tuple(
        tuple(index for index, mask in enumerate(masks) if mask >> sheet & 1)
        for sheet in range(modulus)
    )
    states = 0
    branches = 0

    @lru_cache(maxsize=None)
    def visit(uncovered, remaining):
        nonlocal states, branches
        states += 1
        if uncovered == 0:
            return ()
        if remaining == 0:
            return None
        gains = sorted(
            ((mask & uncovered).bit_count() for mask in masks), reverse=True
        )[:remaining]
        if sum(gains) < uncovered.bit_count():
            return None
        sheets = tuple(
            sheet for sheet in range(modulus) if uncovered >> sheet & 1
        )
        pivot = min(
            sheets,
            key=lambda sheet: sum(
                bool(masks[index] & uncovered) for index in coverers[sheet]
            ),
        )
        options = sorted(
            coverers[pivot],
            key=lambda index: -(masks[index] & uncovered).bit_count(),
        )
        for index in options:
            reduced = uncovered & ~masks[index]
            if reduced == uncovered:
                continue
            branches += 1
            suffix = visit(reduced, remaining - 1)
            if suffix is not None:
                return (index,) + suffix
        return None

    answer = visit(full, cap)
    return answer, len(raw), len(masks), states, branches, visit.cache_info().hits


def period_descent_audit():
    cells = 0
    rows = []
    for modulus in range(3, 102, 2):
        for residue in range(1, modulus + 1):
            order = quotient_order(modulus, residue)
            reduced = residue // (modulus // order)
            observed = half_mask(modulus, residue)
            expected = pullback_mask(half_mask(order, reduced), order, modulus)
            require(observed == expected,
                    ("period descent", modulus, residue, order, reduced))
            cells += modulus
        rows.append((modulus, modulus, cells))
    return tuple(rows), cells


def positive_pullback_audit():
    atoms = (
        (13, (1, 2, 3, 5, 7, 9, 11)),
        (29, (1, 5, 7, 8, 12, 13, 22)),
        (51, (1, 11, 12, 18, 23, 34, 35)),
    )
    multipliers = (1, 3, 5, 7)
    rows = []
    for base, residues in atoms:
        for multiplier in multipliers:
            modulus = base * multiplier
            scaled = tuple(multiplier * residue for residue in residues)
            union = 0
            for residue in scaled:
                union |= half_mask(modulus, residue)
            require(union == (1 << modulus) - 1,
                    ("positive pullback", base, multiplier, scaled))
            rows.append((modulus, base, scaled, len(residues)))
    return tuple(rows)


def finite_classification_audit():
    rows = []
    total_states = 0
    total_branches = 0
    support = []
    for modulus in range(3, 102, 2):
        result = exact_cover_search(modulus, 7)
        observed = result[0] is not None
        expected = any(modulus % divisor == 0 for divisor in SUPPORT_DIVISORS)
        require(observed == expected,
                ("finite classification", modulus, observed, expected, result))
        if observed:
            support.append(modulus)
        total_states += result[3]
        total_branches += result[4]
        rows.append((modulus, observed, result[2], result[3], result[4], result[5]))
    return tuple(rows), tuple(support), total_states, total_branches


def main():
    dependencies = tuple((label, lf_hash(path)) for label, path, _ in PINNED)
    for label, path, expected in PINNED:
        if expected is not None:
            require(lf_hash(path) == expected, (label, lf_hash(path), expected))

    density = exact_size_and_density_audit()
    fibres = fibre_section_audit()
    towers = tower_inequality_audit()
    q85 = q85_invoice()
    q85_search = exact_cover_search(85, 7)
    require(q85_search[0] is None, q85_search)
    period = period_descent_audit()
    positives = positive_pullback_audit()
    finite = finite_classification_audit()

    semantic_surface = (
        dependencies,
        density,
        fibres,
        towers,
        q85,
        q85_search,
        period,
        positives,
        finite,
    )
    semantic_digest = sha256(repr(semantic_surface).encode("ascii")).hexdigest()
    if EXPECTED_SEMANTIC_DIGEST is not None:
        require(semantic_digest == EXPECTED_SEMANTIC_DIGEST,
                (semantic_digest, EXPECTED_SEMANTIC_DIGEST))

    density_rows, density_cells = density
    pure_towers, five_towers = towers
    q85_inactive, q85_bases, q85_rows, q85_inequalities = q85
    q85_kinds = tuple(
        (kind, sum(row[1] == kind for row in q85_rows),
         tuple(sorted({row[4] for row in q85_rows if row[1] == kind})))
        for kind in ("full", "order17-even", "order17-odd")
    )
    period_rows, period_cells = period
    finite_rows, finite_support, finite_states, finite_branches = finite

    print("THM-3434 seventeen-fibre two-sided mass closure")
    print(f"dependency_sha256_lf={dependencies}")
    print(
        "quotient_density_summary="
        f"(rows,cells,first_pure,last_pure,first_five,last_five)="
        f"{(len(density_rows), density_cells, density_rows[0], density_rows[7], density_rows[8], density_rows[-1])}"
    )
    print(f"prime_fibre_section_rows_cells={fibres}")
    print(
        "tower_integer_inequalities="
        f"(pure_rows,pure_min_margin,five_rows,five_min_margin)="
        f"{(len(pure_towers), min(row[-1] for row in pure_towers), len(five_towers), min(row[-1] for row in five_towers))}"
    )
    print(
        "Q85_two_prime_invoice="
        f"(inactive,covered_bases,kinds,budget_rows)="
        f"{(q85_inactive, q85_bases, q85_kinds, q85_inequalities)}"
    )
    print(f"Q85_exact_rank7_search=(answer,raw,unique,states,branches,hits)={q85_search}")
    print(
        "period_descent=(row_count,cells,last_row)="
        f"{(len(period_rows), period_cells, period_rows[-1])}"
    )
    print(f"positive_pullbacks=(Q,base,residues,rank)={positives}")
    print(
        "finite_odd_Q<=101_classification="
        f"(rows,support,states,branches)="
        f"{(len(finite_rows), finite_support, finite_states, finite_branches)}"
    )
    print(f"semantic_sha256={semantic_digest}")
    print("status=VERIFIED_EXACT_COMPANION_FOR_PROVED_THM3434;ordinary_transverse_literal_half_rank7_only;no_LRC14_decrement")


if __name__ == "__main__":
    main()
