#!/usr/bin/env python3
"""Exact companion for THM-3435's dyadic half-twist decomposition.

The all-modulus theorem is the elementary fibre identity proved in the
theorem file.  This standard-library companion checks that identity on a
large finite universe, freezes the four even rank-seven controls, and runs
an inductive finite-exact ordinary-cover census only through Q=362.  It does
not extrapolate the census past that boundary.
"""

import ast
from functools import lru_cache
from hashlib import sha256
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
IDENTITY_LIMIT = 160
CENSUS_LIMIT = 362
KNOWN_BASES = (8, 9, 10, 11, 12, 13, 15, 23, 25, 29, 51)
EXPECTED_EVEN_ATOMS = (14, 38, 68, 148)
WITNESSES = (
    (14, (1, 3, 4, 5, 9, 11, 13)),
    (38, (1, 9, 17, 20, 21, 29, 37)),
    (68, (8, 11, 23, 24, 45, 56, 57)),
    (148, (8, 33, 41, 100, 107, 115, 140)),
)
PINNED = (
    (
        "THM-3416",
        ROOT / "01-canon/theorems/THM-3416-zero-mode-cochain-global-rank-six-support.md",
        "42a9309145de51d1bb6fca0b7c1945302ff37a63a3183e1dfed838c07118e8bf",
    ),
    (
        "THM-3434",
        ROOT / "01-canon/theorems/THM-3434-seventeen-fibre-two-sided-mass-closure.md",
        "52a5a5caed75ad48ab35ce287c15f6cece074c88dddda2b87329792a3df70af7",
    ),
)
EXPECTED_SEMANTIC_DIGEST = "35dfbba192b4b17b72cda5a019692a9af5ee1ce406f3ca4a7d7ff551fe349096"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def in_arc(word, modulus):
    word %= modulus
    return 14 * min(word, modulus - word) < modulus


def half_mask(modulus, residue):
    ambient = 2 * modulus
    answer = 0
    for sheet in range(modulus):
        if in_arc(residue * (2 * sheet + 1), ambient):
            answer |= 1 << sheet
    return answer


def factor_two(modulus):
    exponent = 0
    odd_part = modulus
    while odd_part % 2 == 0:
        exponent += 1
        odd_part //= 2
    return exponent, odd_part


def truncated_valuation(residue, exponent):
    value = residue
    answer = 0
    while answer < exponent and value % 2 == 0:
        answer += 1
        value //= 2
    return answer, value


def fibre_grid_mask(modulus, residue):
    exponent, odd_part = factor_two(modulus)
    valuation, normalized = truncated_valuation(residue, exponent)
    grid_order = 1 << (exponent - valuation)
    denominator = 2 * grid_order * odd_part
    answer = 0
    for base_sheet in range(odd_part):
        base_word = normalized * (2 * base_sheet + 1)
        for fibre_sheet in range(1 << exponent):
            word = base_word + 2 * odd_part * normalized * fibre_sheet
            if in_arc(word, denominator):
                answer |= 1 << (base_sheet + odd_part * fibre_sheet)
    return answer


def low_grid_support(odd_part, normalized, grid_order, base_sheet):
    ambient = 2 * odd_part
    word = normalized * (2 * base_sheet + 1) % ambient
    return 14 * min(word, ambient - word) < grid_order * ambient


def direct_union(modulus, residues):
    answer = 0
    for residue in residues:
        answer |= half_mask(modulus, residue)
    return answer


def translated_mask(mask, modulus, shift):
    answer = 0
    for sheet in range(modulus):
        if mask >> sheet & 1:
            answer |= 1 << ((sheet + shift) % modulus)
    return answer


def identity_audit():
    mask_rows = 0
    fibre_rows = 0
    low_grid_rows = 0
    dense_grid_rows = 0
    deepest_rows = 0
    for modulus in range(2, IDENTITY_LIMIT + 1):
        exponent, odd_part = factor_two(modulus)
        fibre_size = 1 << exponent
        for residue in range(1, 2 * modulus):
            direct = half_mask(modulus, residue)
            reconstructed = fibre_grid_mask(modulus, residue)
            require(direct == reconstructed,
                    (modulus, residue, "fibre reconstruction"))
            mask_rows += 1

            valuation, normalized = truncated_valuation(residue, exponent)
            grid_order = 1 << (exponent - valuation)
            repeat = 1 << valuation
            for base_sheet in range(odd_part):
                selected = tuple(
                    fibre_sheet
                    for fibre_sheet in range(fibre_size)
                    if direct >> (base_sheet + odd_part * fibre_sheet) & 1
                )
                require(len(selected) % repeat == 0,
                        (modulus, residue, base_sheet, selected, repeat))
                reduced = len(selected) // repeat
                if grid_order < 7:
                    expected = low_grid_support(
                        odd_part, normalized, grid_order, base_sheet
                    )
                    require(reduced == int(expected),
                            (modulus, residue, base_sheet, reduced, expected))
                    if selected:
                        classes = {value % grid_order for value in selected}
                        require(len(classes) == 1,
                                (modulus, residue, base_sheet, classes))
                    low_grid_rows += 1
                else:
                    require(
                        reduced in (grid_order // 7, (grid_order + 6) // 7),
                        (modulus, residue, base_sheet, grid_order, reduced),
                    )
                    dense_grid_rows += 1
                fibre_rows += 1

            if valuation == exponent:
                base = half_mask(odd_part, normalized)
                expected_pullback = sum(
                    1 << sheet
                    for sheet in range(modulus)
                    if base >> (sheet % odd_part) & 1
                )
                require(direct == expected_pullback,
                        (modulus, residue, "deepest pullback"))
                # Q does not divide r iff R does not divide s in this branch.
                require((residue % modulus != 0) == (normalized % odd_part != 0),
                        (modulus, residue, normalized, odd_part, "transverse"))
                deepest_rows += 1

    return (
        mask_rows,
        fibre_rows,
        low_grid_rows,
        dense_grid_rows,
        deepest_rows,
    )


def two_sheet_pair_audit():
    rows = 0
    cells = 0
    for odd_part in range(3, 152, 2):
        modulus = 2 * odd_part
        for residue in range(1, modulus, 2):
            if residue == odd_part:
                continue
            first = half_mask(modulus, residue)
            partner = half_mask(modulus, modulus - residue)
            require(first & partner == 0,
                    (modulus, residue, "partner overlap"))
            require(translated_mask(first, modulus, odd_part) == partner,
                    (modulus, residue, "deck partner"))
            widened = 0
            for base_sheet in range(odd_part):
                if low_grid_support(odd_part, residue, 2, base_sheet):
                    widened |= 1 << base_sheet
            pullback = sum(
                1 << sheet
                for sheet in range(modulus)
                if widened >> (sheet % odd_part) & 1
            )
            require(first | partner == pullback,
                    (modulus, residue, "widened pullback"))
            rows += 1
            cells += modulus
    return rows, cells


def witness_audit():
    rows = []
    for modulus, residues in WITNESSES:
        require(len(residues) == 7 and len(set(residues)) == 7,
                (modulus, residues, "owner count"))
        require(all(residue % modulus for residue in residues),
                (modulus, residues, "transverse"))
        masks = tuple(half_mask(modulus, residue) for residue in residues)
        require(direct_union(modulus, residues) == (1 << modulus) - 1,
                (modulus, residues, "cover"))
        multiplicities = tuple(
            sum(mask >> sheet & 1 for mask in masks)
            for sheet in range(modulus)
        )
        orders = tuple(sorted(modulus // gcd(modulus, residue)
                              for residue in residues))
        rows.append(
            (
                modulus,
                residues,
                orders,
                tuple(mask.bit_count() for mask in masks),
                tuple(sorted((value, multiplicities.count(value))
                             for value in set(multiplicities))),
            )
        )

    q7_union = 0
    for residue in range(1, 14):
        if residue != 7:
            q7_union |= half_mask(7, residue)
    require(q7_union != (1 << 7) - 1, "Q7 hostile unexpectedly covers")
    require(rows[0][4] == ((1, 14),), rows[0])
    require(rows[2][4] == ((1, 64), (3, 4)), rows[2])
    return tuple(rows), tuple(
        sheet for sheet in range(7) if q7_union >> sheet & 1
    )


def prime_factors(value):
    answer = []
    divisor = 2
    residual = value
    while divisor * divisor <= residual:
        if residual % divisor == 0:
            answer.append(divisor)
            while residual % divisor == 0:
                residual //= divisor
        divisor = 3 if divisor == 2 else divisor + 2
    if residual > 1:
        answer.append(residual)
    return tuple(answer)


def joint_period_bank(modulus):
    primes = prime_factors(modulus)
    grouped = {}
    raw = 0
    for residue in range(1, 2 * modulus):
        if residue == modulus:
            continue
        sheet_mask = half_mask(modulus, residue)
        augmented = sheet_mask
        for offset, prime in enumerate(primes):
            if residue % prime:
                augmented |= 1 << (modulus + offset)
        if augmented == 0:
            continue
        raw += 1
        old = grouped.get(augmented)
        if old is None or residue < old:
            grouped[augmented] = residue
    unique = tuple(sorted((mask, residue) for mask, residue in grouped.items()))
    maximal = tuple(
        item
        for item in unique
        if not any(
            item[0] != other[0] and item[0] | other[0] == other[0]
            for other in unique
        )
    )
    full = (1 << (modulus + len(primes))) - 1
    return raw, unique, maximal, full, primes


def exact_cover(full, items, cap):
    masks = tuple(mask for mask, _ in items)
    residues = tuple(residue for _, residue in items)
    width = full.bit_length()
    coverers = tuple(
        tuple(index for index, mask in enumerate(masks) if mask >> bit & 1)
        for bit in range(width)
    )
    joined = 0
    for mask in masks:
        joined |= mask
    if joined != full:
        return (), 1, 0, 0

    nodes = 0
    branches = 0

    @lru_cache(maxsize=None)
    def solve(state, slots):
        nonlocal nodes, branches
        nodes += 1
        if state == full:
            return ()
        if slots == 0:
            return None
        missing = full ^ state
        gains = sorted(
            ((mask & missing).bit_count() for mask in masks), reverse=True
        )
        if sum(gains[:slots]) < missing.bit_count():
            return None
        missing_bits = tuple(
            bit for bit in range(width) if missing >> bit & 1
        )
        pivot = min(
            missing_bits,
            key=lambda bit: (
                sum(masks[index] | state != state for index in coverers[bit]),
                bit,
            ),
        )
        options = sorted(
            (
                index
                for index in coverers[pivot]
                if masks[index] | state != state
            ),
            key=lambda index: (
                -(masks[index] & missing).bit_count(), residues[index]
            ),
        )
        for index in options:
            branches += 1
            suffix = solve(state | masks[index], slots - 1)
            if suffix is not None:
                return (index,) + suffix
        return None

    chosen = solve(0, cap)
    if chosen is None:
        return (), nodes, branches, solve.cache_info().hits
    witness = tuple(sorted(residues[index] for index in chosen))
    return witness, nodes, branches, solve.cache_info().hits


def census():
    atoms = []
    searched = []
    pullbacks = []
    skipped_known = []
    for modulus in range(2, CENSUS_LIMIT + 1, 2):
        if any(modulus % base == 0 for base in KNOWN_BASES):
            skipped_known.append(modulus)
            continue
        if any(modulus % atom == 0 for atom in atoms):
            pullbacks.append(modulus)
            continue
        raw, unique, maximal, full, primes = joint_period_bank(modulus)
        witness, nodes, branches, hits = exact_cover(full, maximal, 7)
        if witness:
            sheet_union = direct_union(modulus, witness)
            require(sheet_union == (1 << modulus) - 1,
                    (modulus, witness, "census sheet union"))
            require(gcd(modulus, *witness) == 1,
                    (modulus, witness, "joint period"))
            atoms.append(modulus)
        searched.append(
            (
                modulus,
                raw,
                len(unique),
                len(maximal),
                witness,
                nodes,
                branches,
                hits,
                primes,
            )
        )
    require(tuple(atoms) == EXPECTED_EVEN_ATOMS, atoms)
    expected_witnesses = dict(WITNESSES)
    for row in searched:
        if row[4]:
            require(row[4] == expected_witnesses[row[0]],
                    (row[0], row[4], expected_witnesses[row[0]]))
    positives = tuple(
        modulus
        for modulus in range(2, CENSUS_LIMIT + 1, 2)
        if not any(modulus % base == 0 for base in KNOWN_BASES)
        and any(modulus % atom == 0 for atom in atoms)
    )
    negatives = tuple(row[0] for row in searched if not row[4])
    require(
        all(
            any(value % atom == 0 for atom in atoms)
            for value in positives
        ),
        positives,
    )
    next_boundary = next(
        modulus
        for modulus in range(CENSUS_LIMIT + 1, 1000)
        if modulus % 2 == 0
        and not any(modulus % base == 0 for base in KNOWN_BASES)
        and not any(modulus % atom == 0 for atom in atoms)
    )
    require(next_boundary == 366, next_boundary)
    counts = (
        len(searched),
        sum(row[1] for row in searched),
        sum(row[2] for row in searched),
        sum(row[3] for row in searched),
        sum(row[5] for row in searched),
        sum(row[6] for row in searched),
        sum(row[7] for row in searched),
    )
    events = tuple(
        (row[0], row[4], row[5], row[6], row[7]) for row in searched
    )
    return (
        tuple(atoms),
        positives,
        negatives,
        tuple(pullbacks),
        tuple(skipped_known),
        counts,
        events,
        next_boundary,
    )


def main():
    dependencies = tuple((label, lf_hash(path)) for label, path, _ in PINNED)
    for label, path, expected in PINNED:
        require(lf_hash(path) == expected, (label, lf_hash(path), expected))

    identity = identity_audit()
    pairs = two_sheet_pair_audit()
    witnesses, q7_visible = witness_audit()
    finite = census()
    event_digest = sha256(repr(finite[6]).encode("ascii")).hexdigest()
    semantic_surface = (identity, pairs, witnesses, q7_visible, finite)
    semantic_digest = sha256(repr(semantic_surface).encode("ascii")).hexdigest()
    if EXPECTED_SEMANTIC_DIGEST is not None:
        require(semantic_digest == EXPECTED_SEMANTIC_DIGEST,
                (semantic_digest, EXPECTED_SEMANTIC_DIGEST))

    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
            "assert node")
    require(
        not any(
            isinstance(node, ast.Constant) and isinstance(node.value, float)
            for node in ast.walk(tree)
        ),
        "float literal",
    )

    print("THM-3435 dyadic fibre-grid decomposition -- exact companion")
    print("dependency_sha256_lf=" + repr(dependencies))
    print("status=VERIFIED_EXACT_CONTROLS_FOR_ELEMENTARY_ALL_MODULUS_FIBRE_IDENTITY;FINITE_EXACT_even_target_free_Q2_Q362;no_all_q_even_rank7_classification;no_LRC14_cut")
    print("identity_universe=(Q_max,mask_rows,fibre_rows,low_grid_rows,dense_grid_rows,deepest_rows)="
          + repr((IDENTITY_LIMIT,) + identity))
    print("two_sheet_partner_audit=(rows,cells)=" + repr(pairs))
    print("Q7_naive_descent_hostile_visible_sheets=" + repr(q7_visible))
    print("even_atom_controls=(Q,residues,orders,sizes,multiplicity_hist)="
          + repr(witnesses))
    finite_summary = (
        finite[0], finite[1], finite[2], finite[3], len(finite[4]), finite[5],
        event_digest, finite[7],
    )
    print("finite_census=(atoms,global_positive_target_free,primitive_negative_rows,pullback_rows,known_base_row_count,counts,event_sha256,next_unresolved_boundary)="
          + repr(finite_summary))
    print("semantic_sha256=" + semantic_digest)
    print("script_sha256_lf=" + lf_hash(source))
    print("scope=literal_transverse_half_twist_masks;dyadic_base_projection_retains_grid_support_and_fibre_count_but_loses_selected_cosets;Q362_finite_only;Q366_not_run;LRC14_open")


if __name__ == "__main__":
    main()
