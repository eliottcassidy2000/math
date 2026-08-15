#!/usr/bin/env python3
"""Exact rank-at-most-five scout for zero-mode-cochain divisor ancestry.

This is deliberately a different search route from the union-state BFS and
combination census in lrc_unrestricted_zero_mode_cochain_rank_probe_20260815.
It branches on a rare uncovered sheet/prime-breaker coordinate, memoizes only
the covered state and remaining budget, and tracks no quotient-order profile
assumption.  The output separates primitive covers on Q sheets from the
global divisor minimum on q sheets.

Status while exploratory: FINITE-EXACT only.  No LRC(14) consequence.
"""

from functools import lru_cache
from hashlib import sha256
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCAN_LIMIT = 200


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


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


def divisors(value):
    return tuple(divisor for divisor in range(2, value + 1) if value % divisor == 0)


def danger_mask(q, residue, epsilon):
    modulus = 2 * q
    mask = 0
    for sheet in range(q):
        phase_word = residue * (2 * sheet + epsilon) % modulus
        centred = min(phase_word, modulus - phase_word)
        if 14 * centred < modulus:
            mask |= 1 << sheet
    return mask


def augmented_bank(q, epsilon):
    modulus = q if epsilon == 0 else 2 * q
    primes = prime_factors(modulus)
    grouped = {}
    raw = 0
    for residue in range(1, modulus):
        if residue % q == 0:
            continue
        sheet_mask = danger_mask(q, residue, epsilon)
        if not sheet_mask:
            continue
        raw += 1
        mask = sheet_mask
        for offset, prime in enumerate(primes):
            if residue % prime:
                mask |= 1 << (q + offset)
        old = grouped.get(mask)
        if old is None or residue < old:
            grouped[mask] = residue

    unique = tuple(sorted(((mask, residue) for mask, residue in grouped.items()), key=lambda x: x[1]))
    maximal = tuple(
        item
        for item in unique
        if not any(item[0] != other[0] and item[0] | other[0] == other[0] for other in unique)
    )
    full = (1 << (q + len(primes))) - 1
    return modulus, primes, raw, unique, maximal, full


def rare_coordinate_cover(full, items, cap):
    masks = tuple(mask for mask, _ in items)
    residues = tuple(residue for _, residue in items)
    width = full.bit_length()
    by_bit = tuple(
        tuple(index for index, mask in enumerate(masks) if mask & (1 << bit))
        for bit in range(width)
    )
    all_union = 0
    for mask in masks:
        all_union |= mask
    if all_union != full:
        return (), 0, 0

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
            ((mask & missing).bit_count() for mask in masks),
            reverse=True,
        )
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
        candidates = tuple(
            sorted(
                (index for index in by_bit[pivot] if masks[index] | state != state),
                key=lambda index: (-(masks[index] & missing).bit_count(), residues[index]),
            )
        )
        for index in candidates:
            branches += 1
            answer = solve(state | masks[index], slots - 1)
            if answer is not None:
                return (index,) + answer
        return None

    chosen = solve(0, cap)
    if chosen is None:
        return (), nodes, branches
    witness = tuple(sorted(residues[index] for index in chosen))
    joined = 0
    for index in chosen:
        joined |= masks[index]
    require(joined == full, (cap, witness, joined, full))
    require(len(witness) <= cap, (cap, witness))
    return witness, nodes, branches


def record(q, epsilon):
    modulus, primes, raw, unique, maximal, full = augmented_bank(q, epsilon)
    witness4, nodes4, branches4 = rare_coordinate_cover(full, maximal, 4)
    if witness4:
        rank = len(witness4)
        witness = witness4
        nodes5 = branches5 = 0
    else:
        witness5, nodes5, branches5 = rare_coordinate_cover(full, maximal, 5)
        rank = len(witness5) if witness5 else None
        witness = witness5
    orders = tuple(q // gcd(q, residue) for residue in witness)
    augmented_gcd = gcd(modulus, *witness) if witness else None
    if witness:
        require(augmented_gcd == 1, (q, epsilon, witness, augmented_gcd))
        require(len(set(witness)) == len(witness), (q, epsilon, witness))
    return (
        q,
        epsilon,
        raw,
        len(unique),
        len(maximal),
        rank,
        witness,
        tuple(sorted(orders)),
        nodes4 + nodes5,
        branches4 + branches5,
    )


def main():
    records = tuple(record(q, epsilon) for q in range(2, SCAN_LIMIT + 1) for epsilon in (0, 1))
    primitive = tuple(row for row in records if row[5] is not None)
    primitive_rank4 = tuple((row[0], row[1], row[6], row[7]) for row in primitive if row[5] <= 4)
    primitive_rank5 = tuple((row[0], row[1], row[6], row[7]) for row in primitive if row[5] == 5)
    require(
        primitive_rank4
        == (
            (8, 1, (1, 3, 5, 7), (8, 8, 8, 8)),
            (9, 1, (1, 5, 6, 7), (3, 9, 9, 9)),
        ),
        primitive_rank4,
    )

    by_key = {(row[0], row[1]): row for row in records}
    global_rows = []
    for q in range(2, SCAN_LIMIT + 1):
        candidates = []
        for quotient in divisors(q):
            for epsilon in (0, 1):
                rank = by_key[(quotient, epsilon)][5]
                if rank is not None:
                    candidates.append((rank, quotient, epsilon))
        best = min((candidate[0] for candidate in candidates), default=None)
        minimizers = tuple(candidate[1:] for candidate in candidates if candidate[0] == best)
        global_rows.append((q, best, minimizers))
    global_rows = tuple(global_rows)
    global_rank5 = tuple(row for row in global_rows if row[1] == 5)
    predicted_global_rank5 = tuple(
        q
        for q in range(2, SCAN_LIMIT + 1)
        if (q % 10 == 0 or q % 12 == 0) and q % 8 and q % 9
    )
    require(tuple(row[0] for row in global_rank5) == predicted_global_rank5, global_rank5)

    event_digest = sha256(repr(records).encode("ascii")).hexdigest()
    print("status=FINITE-EXACT exploratory_rank_at_most_five_scan;no_LRC14_consequence")
    print(f"scan_limit={SCAN_LIMIT}")
    print(f"primitive_rank4=(Q,epsilon,residues,orders)={primitive_rank4}")
    print(f"primitive_rank5=(Q,epsilon,residues,orders)={primitive_rank5}")
    print(f"global_rank5=(q,rank,minimizing_(Q,epsilon))={global_rank5}")
    print("global_rank5_conjecture=(10|q_or_12|q)_and_8_not|q_and_9_not|q")
    print(
        "search_totals=(raw,unique,maximal,nodes,branches)="
        f"{(sum(row[2] for row in records), sum(row[3] for row in records), sum(row[4] for row in records), sum(row[8] for row in records), sum(row[9] for row in records))}"
    )
    print(f"event_sha256={event_digest}")


if __name__ == "__main__":
    main()
