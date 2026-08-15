#!/usr/bin/env python3
"""Exact companion for THM-3429 prime-fibre descent.

This verifies the elementary projection/fibre identities and the Q=51 lift
cocycle hostile.  The all-q proof is analytic and depends on the pinned
THM-3416 and THM-3421 classifications; this finite bank does not replace it.
"""

from fractions import Fraction
from hashlib import sha256
from math import gcd, isqrt, lcm
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
)
CONTROL_MODULI = (21, 35, 49, 51, 77, 85, 87, 119, 121, 143, 169, 203, 289, 493, 529, 841)
EXPECTED_PRIME_SOLUTIONS = (2, 3, 5, 11, 17, 23, 29)
TRANSVERSE_UNIVERSE = "Q_does_not_divide_each_selected_residue"
EXPECTED_SEMANTIC_DIGEST = "c46b24a98d32f83e1933276372e8bc1688ac9d3c3efa54c154ec19fcd1094c2d"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def prime(value):
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    return all(value % divisor for divisor in range(3, isqrt(value) + 1, 2))


def prime_divisors(value):
    answer = []
    residual = value
    divisor = 2
    while divisor * divisor <= residual:
        if residual % divisor == 0:
            answer.append(divisor)
            while residual % divisor == 0:
                residual //= divisor
        divisor = 3 if divisor == 2 else divisor + 2
    if residual > 1:
        answer.append(residual)
    return tuple(answer)


def half_mask(modulus, residue):
    result = 0
    double = 2 * modulus
    for sheet in range(modulus):
        word = residue * (2 * sheet + 1) % double
        if 14 * min(word, double - word) < double:
            result |= 1 << sheet
    return result


def ceil_seventh(value):
    return (value + 6) // 7


def lifted_mask(mask, small_modulus, large_modulus):
    return sum(
        1 << sheet
        for sheet in range(large_modulus)
        if mask >> (sheet % small_modulus) & 1
    )


def prime_capacity_solutions():
    # For p>36, 6*ceil(p/7) <= 6*(p+6)/7 < p.
    require(6 * (37 + 6) < 7 * 37, "capacity cutoff")
    solutions = tuple(
        value
        for value in range(2, 37)
        if prime(value) and value <= 6 * ceil_seventh(value)
    )
    require(solutions == EXPECTED_PRIME_SOLUTIONS, solutions)
    require(7 not in solutions and ceil_seventh(7) == 1, "strict p=7 boundary")
    return solutions


def projection_and_fibre_audit():
    projection_cells = 0
    active_fibres = 0
    active_points = 0
    strict_seven_maxima = []
    records = []
    for modulus in CONTROL_MODULI:
        modulus_masks = tuple(half_mask(modulus, residue) for residue in range(2 * modulus))
        for divisor in prime_divisors(modulus):
            base_modulus = modulus // divisor
            fixed_base = (base_modulus - 1) // 2
            maximum = 0
            inactive_even = 0
            inactive_odd = 0
            for residue in range(1, 2 * modulus):
                if residue % modulus == 0:
                    continue
                mask = modulus_masks[residue]
                if residue % divisor == 0:
                    base = half_mask(base_modulus, residue // divisor)
                    require(
                        mask == lifted_mask(base, base_modulus, modulus),
                        ("inactive projection", modulus, divisor, residue),
                    )
                    projection_cells += modulus
                    fixed_section = tuple(
                        bool(mask >> ((fixed_base + base_modulus * t) % modulus) & 1)
                        for t in range(divisor)
                    )
                    expected = (residue // divisor) % 2 == 0
                    require(all(value == expected for value in fixed_section),
                            ("inactive fixed parity", modulus, divisor, residue, fixed_section))
                    if expected:
                        inactive_even += 1
                    else:
                        inactive_odd += 1
                    continue

                local_fixed = half_mask(divisor, residue)
                observed_fixed = 0
                for t in range(divisor):
                    sheet = (fixed_base + base_modulus * t) % modulus
                    if mask >> sheet & 1:
                        observed_fixed |= 1 << t
                require(observed_fixed == local_fixed,
                        ("active fixed fibre", modulus, divisor, residue))

                for base_sheet in range(base_modulus):
                    section_size = sum(
                        bool(mask >> ((base_sheet + base_modulus * t) % modulus) & 1)
                        for t in range(divisor)
                    )
                    require(section_size <= ceil_seventh(divisor),
                            ("active fibre capacity", modulus, divisor, residue,
                             base_sheet, section_size))
                    maximum = max(maximum, section_size)
                    active_fibres += 1
                    active_points += section_size
            if divisor == 7:
                require(maximum == 1, ("strict seven maximum", modulus, maximum))
                strict_seven_maxima.append((modulus, maximum))
            records.append(
                (modulus, divisor, base_modulus, maximum, inactive_even, inactive_odd)
            )
    require(tuple(strict_seven_maxima) == ((21, 1), (35, 1), (49, 1),
                                           (77, 1), (119, 1), (203, 1)),
            strict_seven_maxima)
    return tuple(records), projection_cells, active_fibres, active_points, tuple(strict_seven_maxima)


def fraction_is_dangerous(value):
    denominator = value.denominator
    residue = value.numerator % denominator
    return 14 * min(residue, denominator - residue) < denominator


def q51_cocycle_hostile():
    modulus = 51
    prime_coordinate = 17
    base_modulus = 3
    residues = (1, 11, 12, 18, 23, 34, 35)
    orders = tuple(modulus // gcd(modulus, residue) for residue in residues)
    require(orders == (51, 51, 17, 17, 51, 3, 51), orders)
    require(lcm(*orders) == modulus, orders)

    masks = {residue: half_mask(modulus, residue) for residue in residues}
    union = 0
    for mask in masks.values():
        union |= mask
    require(union == (1 << modulus) - 1, "Q51 literal cover")

    expected = (
        (
            0,
            (),
            ((1, (0, 1, 16)), (11, (3, 6, 9)), (12, (4, 7, 14)),
             (18, (13, 14, 15)), (23, (2, 5, 8)), (35, (10, 11, 12))),
            17,
            1,
        ),
        (
            1,
            (34,),
            ((1, (0, 16)), (11, (1, 15)), (12, (1, 8, 15)),
             (18, (7, 8, 9)), (23, (1, 15)), (35, (0, 16))),
            7,
            7,
        ),
        (
            2,
            (),
            ((1, (0, 15, 16)), (11, (7, 10, 13)), (12, (2, 9, 12)),
             (18, (1, 2, 3)), (23, (8, 11, 14)), (35, (4, 5, 6))),
            17,
            1,
        ),
    )
    rows = []
    cocycle_cells = 0
    for base_sheet in range(base_modulus):
        inactive_whole = tuple(
            residue
            for residue in residues
            if residue % prime_coordinate == 0
            and all(
                masks[residue] >> ((base_sheet + base_modulus * t) % modulus) & 1
                for t in range(prime_coordinate)
            )
        )
        sections = []
        for residue in residues:
            if residue % prime_coordinate == 0:
                continue
            section = tuple(
                t
                for t in range(prime_coordinate)
                if masks[residue] >> ((base_sheet + base_modulus * t) % modulus) & 1
            )
            sections.append((residue, section))

            reduced = residue % (2 * prime_coordinate)
            lift = (residue - reduced) // (2 * prime_coordinate)
            for t in range(prime_coordinate):
                phase = (
                    Fraction(reduced * t, prime_coordinate)
                    + Fraction(reduced * (2 * base_sheet + 1),
                               2 * prime_coordinate * base_modulus)
                    + Fraction(lift * (2 * base_sheet + 1), base_modulus)
                )
                observed = t in section
                require(observed == fraction_is_dangerous(phase),
                        ("affine cocycle", base_sheet, residue, t, phase, observed))
                cocycle_cells += 1

        active_union = set().union(*(set(section) for _, section in sections))
        active_mass = sum(len(section) for _, section in sections)
        rows.append(
            (base_sheet, inactive_whole, tuple(sections), len(active_union),
             active_mass - len(active_union))
        )
    require(tuple(rows) == expected, rows)

    by_row = {row[0]: dict(row[2]) for row in rows}
    require(set(by_row[0][1]).isdisjoint(by_row[0][35]), "Q51 y0 lift separation")
    require(set(by_row[2][1]).isdisjoint(by_row[2][35]), "Q51 y2 lift separation")
    require(by_row[1][1] == by_row[1][35], "Q51 fixed-fibre cocycle cancellation")
    return residues, orders, tuple(rows), cocycle_cells


def prime_positive_controls():
    witnesses = (
        (11, (1, 2, 3, 5, 7, 9)),
        (13, (1, 2, 3, 5, 7, 9, 11)),
        (23, (1, 4, 5, 7, 9, 11)),
        (29, (1, 5, 7, 8, 12, 13, 22)),
    )
    rows = []
    for modulus, residues in witnesses:
        union = 0
        for residue in residues:
            require(gcd(modulus, residue) == 1, (modulus, residue))
            union |= half_mask(modulus, residue)
        require(union == (1 << modulus) - 1, ("prime control", modulus))
        rows.append((modulus, len(residues)))
    return tuple(rows)


def transverse_boundary_controls():
    rows = tuple(
        (modulus, half_mask(modulus, 0).bit_count(), half_mask(modulus, modulus).bit_count())
        for modulus in (13, 31, 51)
    )
    require(rows == ((13, 13, 0), (31, 31, 0), (51, 51, 0)), rows)
    return rows


def main():
    dependencies = tuple((label, lf_hash(path)) for label, path, _ in PINNED)
    for label, path, expected in PINNED:
        require(lf_hash(path) == expected, (label, lf_hash(path), expected))

    prime_solutions = prime_capacity_solutions()
    fibre_audit = projection_and_fibre_audit()
    q51 = q51_cocycle_hostile()
    prime_controls = prime_positive_controls()
    transverse_controls = transverse_boundary_controls()
    semantic_surface = (
        TRANSVERSE_UNIVERSE,
        dependencies,
        prime_solutions,
        fibre_audit,
        q51,
        prime_controls,
        transverse_controls,
    )
    semantic_digest = sha256(repr(semantic_surface).encode("ascii")).hexdigest()
    if EXPECTED_SEMANTIC_DIGEST is not None:
        require(semantic_digest == EXPECTED_SEMANTIC_DIGEST,
                (semantic_digest, EXPECTED_SEMANTIC_DIGEST))

    print("THM-3429 prime-fibre activity descent")
    print(f"typed_universe={TRANSVERSE_UNIVERSE};r=0_universal_and_r=Q_empty_are_excluded")
    print(f"dependency_sha256_lf={dependencies}")
    print(f"capacity_prime_solutions={prime_solutions}")
    print("mixed_target_free_odd_primes=(3,5,17,29);all_active_target_free_primes=(13,29)")
    print(f"projection_fibre_audit=(records,projection_cells,active_fibres,active_points,p7_controls)={fibre_audit}")
    print(f"Q51_affine_lift_cocycle=(residues,orders,fibres,cells)={q51}")
    print(f"prime_positive_controls=(Q,rank)={prime_controls}")
    print(f"transverse_boundary_controls=(Q,r0_size,rQ_size)={transverse_controls}")
    print(f"semantic_sha256={semantic_digest}")
    print("status=PROVED_EXACT_COMPANION;analytic_reduction_not_finite_classification;residual_17adic_towers_open;no_LRC14_decrement")


if __name__ == "__main__":
    main()
