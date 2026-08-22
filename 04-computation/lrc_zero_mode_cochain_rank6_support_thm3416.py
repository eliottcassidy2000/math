#!/usr/bin/env python3
"""Exact companion for THM-3416, global zero-mode-cochain rank six.

The all-q converse is analytic.  This companion freezes every finite cutoff
used by that proof, the reflection/CRT hostile at orders 17 and 29, the four
positive atoms, and an independent primitive cap-six census through Q=200.
The census is finite-only and is not used to extrapolate the theorem.

Only the Python standard library is used.  Every gate remains active under
``python -O``.
"""

from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations_with_replacement
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCAN_LIMIT = 200
LOWER_BASES = (8, 9, 10, 12)
NEW_BASES = (11, 15, 23, 25)
ALL_BASES = LOWER_BASES + NEW_BASES
PERIOD = 455_400

PINNED = (
    (
        "THM-3405-two-gauge",
        ROOT / "01-canon/theorems/THM-3405-common-centre-gcd-gauge-and-boolean-half-twist.md",
        "d3e7dbeeb85c6f897bd9e31270bd0b6602ae4feac3b46a45eb5ce23ae5d24fe0",
    ),
    (
        "THM-3414-fixed-zero-six",
        ROOT / "01-canon/theorems/THM-3414-fixed-zero-six-owner-base-classification.md",
        "5568a4e93bc4478566335e2722c488c999797462eeb7c95af364b20dba41e998",
    ),
    (
        "THM-3415-rank-five",
        ROOT / "01-canon/theorems/THM-3415-zero-mode-cochain-global-rank-five-support.md",
        "de8d2f615695070dc75cad38ad4a9c22308d8bf900cd6f9a69cd2003f815a14d",
    ),
)

EXPECTED_SEMANTIC_DIGEST = "99892baf39b3d2b1b6a802bf21e0fe4164f155030d8ad051bc7ae26513b01ca3"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def prime_factors(value):
    factors = []
    trial = 2
    while trial * trial <= value:
        if value % trial == 0:
            factors.append(trial)
            while value % trial == 0:
                value //= trial
        trial += 1
    if value > 1:
        factors.append(value)
    return tuple(factors)


def danger_mask(q, residue, epsilon):
    modulus = 2 * q
    mask = 0
    for sheet in range(q):
        word = residue * (2 * sheet + epsilon) % modulus
        centred = min(word, modulus - word)
        if 14 * centred < modulus:
            mask |= 1 << sheet
    return mask


def zero_count(order):
    return 1 + 2 * ((order - 1) // 14)


def odd_coset_count(order):
    return 2 * ((((order - 1) // 7) + 1) // 2)


def half_count(order):
    odd = odd_coset_count(order)
    return odd if order % 2 == 0 else max(odd, zero_count(order))


def exact_order_masks(q, order, epsilon):
    modulus = q if epsilon == 0 else 2 * q
    return tuple(
        sorted(
            {
                danger_mask(q, residue, epsilon)
                for residue in range(1, modulus)
                if residue % q
                and q // gcd(q, residue) == order
                and danger_mask(q, residue, epsilon)
            }
        )
    )


def formula_audit():
    # Direct masks audit the closed count formula on a finite hostile range.
    for order in range(2, 201):
        observed = max(
            (mask.bit_count() for mask in exact_order_masks(order, order, 1)),
            default=0,
        )
        require(observed == half_count(order), (order, observed, half_count(order)))
        require(7 * half_count(order) <= order + 6, (order, half_count(order)))

    h6 = tuple(order for order in range(2, 37) if 6 * half_count(order) >= order)
    expected_h6 = (3, 5, 8, 9, 10, 11, 12, 15, 17, 22, 23, 24, 29, 36)
    require(h6 == expected_h6, h6)
    require(6 * (37 + 6) < 7 * 37, "strict one-sixth tail")

    small_low = tuple(
        (Fraction(half_count(order), order), order)
        for order in range(2, 43)
        if not any(base and order % base == 0 for base in ALL_BASES)
        and order not in (3, 5, 17, 29)
    )
    require(max(small_low) == (Fraction(6, 37), 37), max(small_low))
    require(half_count(43) == 7, half_count(43))
    require(Fraction(43 + 6, 7 * 43) == Fraction(7, 43), "tail equality")
    return h6, max(small_low)


def complement_classifier(base):
    # The analytic ceiling bound in the theorem reduces the divisible case
    # to m<=32; the coprime case is already reduced to m<=15.
    winners = []
    records = []
    for order in range(2, 33):
        q = lcm(base, order)
        base_masks = exact_order_masks(q, base, 1)
        require(len(base_masks) == 1, (base, order, q, len(base_masks)))
        complement = ((1 << q) - 1) ^ base_masks[0]
        contribution = max(
            ((mask & complement).bit_count() for mask in exact_order_masks(q, order, 1)),
            default=0,
        )
        if 5 * contribution >= complement.bit_count():
            winners.append(order)
            records.append(
                (
                    order,
                    q,
                    complement.bit_count(),
                    contribution,
                    "tie" if 5 * contribution == complement.bit_count() else "strict",
                )
            )
    expected = {
        3: (5, 8, 9, 10, 12, 15),
        5: (3, 8, 9, 10, 25),
    }[base]
    require(tuple(winners) == expected, (base, winners))
    return tuple(records)


def union_maxima(masks, cap):
    maxima = []
    for size in range(1, cap + 1):
        best = 0
        for packet in combinations_with_replacement(masks, size):
            joined = 0
            for mask in packet:
                joined |= mask
            best = max(best, joined.bit_count())
        maxima.append(best)
    return tuple(maxima)


def reflection_crt_audit():
    records = []
    for order, expected_sizes, expected_unions in (
        (17, (2, 3), (3, 5, 7, 9, 11, 13)),
        (29, (4, 5), (5, 9, 13, 17, 21, 25)),
    ):
        masks = exact_order_masks(order, order, 1)
        require(
            all(
                all(
                    ((mask >> sheet) & 1)
                    == ((mask >> ((-1 - sheet) % order)) & 1)
                    for sheet in range(order)
                )
                for mask in masks
            ),
            (order, "reflection"),
        )
        sizes = tuple(sorted({mask.bit_count() for mask in masks}))
        maxima = union_maxima(masks, 6)
        require(sizes == expected_sizes, (order, sizes))
        require(maxima == expected_unions, (order, maxima))
        records.append((order, len(masks), sizes, maxima))

    q = 17 * 29
    masks17 = exact_order_masks(q, 17, 1)
    masks29 = exact_order_masks(q, 29, 1)
    for order, masks in ((17, masks17), (29, masks29)):
        require(
            all(
                all(
                    ((mask >> sheet) & 1) == ((mask >> (sheet % order)) & 1)
                    for sheet in range(q)
                )
                for mask in masks
            ),
            (q, order, "cylinder"),
        )
    full = (1 << q) - 1
    for left in masks17:
        for right in masks29:
            missed = full ^ (left | right)
            expected = (17 - left.bit_count() // 29) * (29 - right.bit_count() // 17)
            require(missed.bit_count() == expected, (left.bit_count(), right.bit_count()))

    deficits = []
    for count17 in range(7):
        for count29 in range(7 - count17):
            if count17 + count29 == 0:
                continue
            low = 6 - count17 - count29
            core = Fraction(1)
            if count17:
                core *= Fraction(16 - 2 * count17, 17)
            if count29:
                core *= Fraction(28 - 4 * count29, 29)
            deficits.append((core - Fraction(7 * low, 43), count17, count29, low))
    require(min(deficits) == (Fraction(7, 731), 1, 0, 5), min(deficits))
    mixed_min = min(row for row in deficits if row[1] and row[2])
    require(mixed_min == (Fraction(644, 21_199), 1, 1, 4), mixed_min)
    return tuple(records), min(deficits), mixed_min


def witness_record(q, epsilon, residues):
    masks = tuple(danger_mask(q, residue, epsilon) for residue in residues)
    full = (1 << q) - 1
    joined = 0
    multiplicities = []
    for mask in masks:
        joined |= mask
    require(joined == full, (q, epsilon, residues, joined))
    for sheet in range(q):
        multiplicities.append(sum((mask >> sheet) & 1 for mask in masks))
    require(gcd(q if epsilon == 0 else 2 * q, *residues) == 1, (q, epsilon, residues))
    return (
        q,
        epsilon,
        residues,
        tuple(mask.bit_count() for mask in masks),
        tuple(sorted((value, multiplicities.count(value)) for value in set(multiplicities))),
        sha256(repr(masks).encode("ascii")).hexdigest(),
    )


def positive_atoms():
    atoms = (
        witness_record(11, 1, (1, 2, 3, 5, 7, 9)),
        witness_record(15, 0, (1, 2, 3, 4, 5, 7)),
        witness_record(15, 1, (1, 4, 6, 7, 8, 10)),
        witness_record(23, 1, (1, 4, 5, 7, 9, 11)),
        witness_record(25, 1, (1, 9, 10, 11, 19, 21)),
    )
    require(atoms[0][4] == ((1, 11),), atoms[0])
    require(atoms[1][4] == ((1, 14), (6, 1)), atoms[1])
    require(atoms[2][4] == ((1, 14), (4, 1)), atoms[2])
    require(atoms[3][4] == ((1, 23),), atoms[3])
    require(atoms[4][4] == ((1, 25),), atoms[4])
    return atoms


def augmented_bank(q, epsilon):
    modulus = q if epsilon == 0 else 2 * q
    primes = prime_factors(modulus)
    grouped = {}
    for residue in range(1, modulus):
        if residue % q == 0:
            continue
        sheet_mask = danger_mask(q, residue, epsilon)
        if not sheet_mask:
            continue
        mask = sheet_mask
        for offset, prime in enumerate(primes):
            if residue % prime:
                mask |= 1 << (q + offset)
        old = grouped.get(mask)
        if old is None or residue < old:
            grouped[mask] = residue
    unique = tuple(sorted(((mask, residue) for mask, residue in grouped.items()), key=lambda row: row[1]))
    maximal = tuple(
        row
        for row in unique
        if not any(row[0] != other[0] and row[0] | other[0] == other[0] for other in unique)
    )
    full = (1 << (q + len(primes))) - 1
    return maximal, full


def rare_coordinate_cover(items, full, cap):
    masks = tuple(mask for mask, _ in items)
    residues = tuple(residue for _, residue in items)
    width = full.bit_length()
    by_bit = tuple(
        tuple(index for index, mask in enumerate(masks) if mask & (1 << bit))
        for bit in range(width)
    )
    joined = 0
    for mask in masks:
        joined |= mask
    if joined != full:
        return ()

    @lru_cache(maxsize=None)
    def solve(state, slots):
        if state == full:
            return ()
        if slots == 0:
            return None
        missing = full ^ state
        gains = sorted((mask & missing).bit_count() for mask in masks)[::-1]
        if sum(gains[:slots]) < missing.bit_count():
            return None
        missing_bits = tuple(bit for bit in range(width) if missing & (1 << bit))
        pivot = min(
            missing_bits,
            key=lambda bit: (
                sum(1 for index in by_bit[bit] if masks[index] | state != state),
                bit,
            ),
        )
        candidates = sorted(
            (index for index in by_bit[pivot] if masks[index] | state != state),
            key=lambda index: (-(masks[index] & missing).bit_count(), residues[index]),
        )
        for index in candidates:
            answer = solve(state | masks[index], slots - 1)
            if answer is not None:
                return (index,) + answer
        return None

    selected = solve(0, cap)
    if selected is None:
        return ()
    witness = tuple(sorted(residues[index] for index in selected))
    require(len(witness) <= cap and len(set(witness)) == len(witness), witness)
    return witness


def primitive_census():
    rows = []
    keys = []
    for q in range(2, SCAN_LIMIT + 1):
        for epsilon in (0, 1):
            bank, full = augmented_bank(q, epsilon)
            witness = rare_coordinate_cover(bank, full, 6)
            if witness:
                modulus = q if epsilon == 0 else 2 * q
                require(gcd(modulus, *witness) == 1, (q, epsilon, witness))
                keys.append((q, epsilon))
                rows.append((q, epsilon, witness))

    actual_support = tuple(
        q
        for q in range(2, SCAN_LIMIT + 1)
        if any(q % key_q == 0 for key_q, _ in keys)
    )
    predicted_support = tuple(
        q
        for q in range(2, SCAN_LIMIT + 1)
        if any(q % base == 0 for base in ALL_BASES)
    )
    require(actual_support == predicted_support, (actual_support, predicted_support))
    actual_rank6 = tuple(
        q
        for q in actual_support
        if not any(q % base == 0 for base in LOWER_BASES)
    )
    predicted_rank6 = tuple(
        q
        for q in range(2, SCAN_LIMIT + 1)
        if any(q % base == 0 for base in NEW_BASES)
        and not any(q % base == 0 for base in LOWER_BASES)
    )
    require(actual_rank6 == predicted_rank6, (actual_rank6, predicted_rank6))
    epsilon_counts = tuple(sum(epsilon == value for _, epsilon in keys) for value in (0, 1))
    return (
        tuple(rows),
        epsilon_counts,
        actual_support,
        actual_rank6,
        sha256(repr(tuple(rows)).encode("ascii")).hexdigest(),
    )


def rank6_condition(q):
    return (
        any(q % base == 0 for base in NEW_BASES)
        and not any(q % base == 0 for base in LOWER_BASES)
    )


def fibonacci_rank(modulus):
    previous, current = 0, 1
    for index in range(1, 10 * modulus + 1):
        previous, current = current, (previous + current) % modulus
        if previous == 0:
            return index
    raise RuntimeError((modulus, "rank not found"))


def arithmetic_corollaries():
    rank6_count = sum(rank6_condition(q) for q in range(1, PERIOD + 1))
    support_count = sum(
        any(q % base == 0 for base in ALL_BASES)
        for q in range(1, PERIOD + 1)
    )
    require(rank6_count == 55_000, rank6_count)
    require(Fraction(rank6_count, PERIOD) == Fraction(25, 207), rank6_count)
    require(Fraction(support_count, PERIOD) == Fraction(149, 345), support_count)

    fib_ranks = tuple((base, fibonacci_rank(base)) for base in ALL_BASES)
    require(
        fib_ranks == ((8, 6), (9, 12), (10, 15), (12, 12), (11, 10), (15, 20), (23, 24), (25, 25)),
        fib_ranks,
    )
    fib = [0, 1]
    for _ in range(2, 601):
        fib.append((fib[-1] + fib[-2]) % PERIOD)
    fib_residues = tuple(
        n % 600
        for n in range(1, 601)
        if rank6_condition(fib[n])
    )
    formula_residues = tuple(
        n % 600
        for n in range(1, 601)
        if (n % 10 == 0 or n % 24 == 0 or n % 25 == 0)
        and n % 6
        and n % 15
    )
    require(fib_residues == formula_residues, (fib_residues, formula_residues))
    require(len(fib_residues) == 48, len(fib_residues))

    spine = tuple(
        n
        for n in range(99)
        if rank6_condition(4 * n * n + 12 * n + 11)
    )
    spine_formula = tuple(
        n
        for n in range(99)
        if n % 11 in (0, 8) and n % 9 not in (1, 5)
    )
    require(spine == spine_formula, (spine, spine_formula))
    require(len(spine) == 14, len(spine))
    return rank6_count, support_count, fib_ranks, fib_residues, spine


def main():
    pin_rows = []
    for label, path, expected in PINNED:
        observed = lf_hash(path)
        require(observed == expected, (label, observed, expected))
        pin_rows.append((label, observed))

    formulas = formula_audit()
    complement3 = complement_classifier(3)
    complement5 = complement_classifier(5)
    reflection = reflection_crt_audit()
    atoms = positive_atoms()
    census = primitive_census()
    arithmetic = arithmetic_corollaries()

    semantic_payload = (
        formulas,
        complement3,
        complement5,
        reflection,
        atoms,
        census[1:4],
        census[4],
        arithmetic,
    )
    semantic = sha256(repr(semantic_payload).encode("ascii")).hexdigest()
    if EXPECTED_SEMANTIC_DIGEST is not None:
        require(semantic == EXPECTED_SEMANTIC_DIGEST, (semantic, EXPECTED_SEMANTIC_DIGEST))

    print("THM-3416 exact zero-mode-cochain global rank-six companion")
    print(f"pinned_dependencies={tuple(pin_rows)}")
    print(f"half_density=(H6,small_low_max,q43)={(formulas[0], formulas[1], (43, half_count(43), Fraction(7, 43)))}")
    print(f"order3_complement_classifier={complement3}")
    print(f"order5_complement_classifier={complement5}")
    print(f"reflection_union_and_CRT_deficits={reflection}")
    print(f"positive_atoms={atoms}")
    print(
        "primitive_cap6_Q2_Q200="
        f"(epsilon_counts,support_count,rank6_count,row_sha256)="
        f"{(census[1], len(census[2]), len(census[3]), census[4])}"
    )
    print(f"rank6_Q2_Q200={census[3]}")
    print(
        "periodic_arithmetic="
        f"(period,rank6_count,density,support_count,density)="
        f"{(PERIOD, arithmetic[0], Fraction(25, 207), arithmetic[1], Fraction(149, 345))}"
    )
    print(
        "fibonacci_period600=(ranks,count,residues,density)="
        f"{(arithmetic[2], len(arithmetic[3]), arithmetic[3], Fraction(2, 25))}"
    )
    print(f"berggren_U_period99=(count,residues,density)={(len(arithmetic[4]), arithmetic[4], Fraction(14, 99))}")
    print("scope=all-q proof is analytic; primitive cap-six census is finite-exact only through Q=200; no physical-time or LRC(14) consequence")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
