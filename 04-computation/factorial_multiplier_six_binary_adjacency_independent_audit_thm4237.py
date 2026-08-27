#!/usr/bin/env python3
"""Independent set-of-bit-positions audit for THM-4237.

This implementation does not import the primary companion and does not use
its bitwise submask predicate.  Binary expansions are converted to explicit
sets of occupied positions and the THM-3474 candidate intersection is rebuilt
literally.
"""

from hashlib import sha256


ODD_H_BOUND = 1 << 14
PRIME_BOUND = 350
EXPONENT_BOUND = 19


def bit_positions(integer: int):
    answer = set()
    position = 0
    while integer:
        if integer % 2:
            answer.add(position)
        integer //= 2
        position += 1
    return frozenset(answer)


def adjacent_by_positions(integer: int) -> bool:
    positions = bit_positions(integer)
    return any(position + 1 in positions for position in positions)


def contained(left: int, right: int) -> bool:
    return bit_positions(left) <= bit_positions(right)


def prime(candidate: int) -> bool:
    if candidate < 2:
        return False
    for divisor in range(2, candidate):
        if divisor * divisor > candidate:
            break
        if candidate % divisor == 0:
            return False
    return True


def literal_candidates(prime_value: int, exponent: int):
    height = prime_value**exponent
    limit = min(5, (prime_value - 1) // 2)
    return tuple(
        t for t in range(1, limit + 1) if contained(t * height, 6 * height)
    )


def main():
    exhaustive = True
    odd_overlap = True
    even_reduction = True
    cells = 0
    for height in range(1, ODD_H_BOUND, 2):
        positions = bit_positions(height)
        adjacent = any(position + 1 in positions for position in positions)
        for limit in (3, 5):
            actual = tuple(
                t for t in range(1, limit + 1) if contained(t * height, 6 * height)
            )
            predicted = () if adjacent else ((2,) if limit == 3 else (2, 4))
            exhaustive &= actual == predicted
            for t in range(1, limit + 1, 2):
                odd_overlap &= bool(
                    bit_positions(t * height) & bit_positions((6 - t) * height)
                )
            even_reduction &= (
                bool(bit_positions(2 * height) & bit_positions(4 * height))
                == adjacent
            )
            cells += 1
    print(
        f"independent_odd_height_cells={cells} set_containment={exhaustive} "
        f"odd_low_bit_overlap={odd_overlap} even_shift_adjacency={even_reduction}"
    )

    rows = []
    congruence = True
    classification = True
    for p in range(7, PRIME_BOUND):
        if not prime(p):
            continue
        for k in range(1, EXPONENT_BOUND + 1):
            height = p**k
            candidates = literal_candidates(p, k)
            classification &= (not candidates) == adjacent_by_positions(height)
            if p % 4 == 3 and k % 2 == 1:
                congruence &= tuple(sorted(bit_positions(height) & {0, 1})) == (0, 1)
                congruence &= not candidates
            rows.append((p, k, adjacent_by_positions(height), candidates))
    print(
        f"independent_prime_power_cells={len(rows)} "
        f"classification={classification} congruence_family={congruence}"
    )
    print(f"independent_semantic_sha256={sha256(repr(rows).encode()).hexdigest()}")

    hostile_height = 23**2
    hostile = literal_candidates(23, 2)
    exact_difference = bit_positions(6 * hostile_height) - bit_positions(
        2 * hostile_height
    )
    print(
        f"hostile_rebuilt=H:{hostile_height},bits:{sorted(bit_positions(hostile_height))},"
        f"candidates:{hostile},difference_for_t2:{sorted(exact_difference)}"
    )


if __name__ == "__main__":
    main()
