#!/usr/bin/env python3
"""Exact companion for THM-3421, the prime half-twist cap-seven theorem.

The proof is elementary.  This standard-library companion freezes the
24-ratio complement, its sharp dyadic three-colouring, the modular stability
bound, every normalized scalar-feasible profile below the analytic threshold,
and the four positive literal atoms.  Every truth gate remains active under
``python -O``.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from math import gcd, isqrt


EXPECTED_SEMANTIC_DIGEST = "f0c7ceb122e90fd3792d72e80a205a6c4bd823cfb5b04064764b95ddc9d9471d"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    for divisor in range(2, isqrt(value) + 1):
        if value % divisor == 0:
            return False
    return True


def valuation_two(value: int) -> int:
    require(value != 0, "zero has no two-adic valuation")
    value = abs(value)
    exponent = 0
    while value % 2 == 0:
        exponent += 1
        value //= 2
    return exponent


def fraction_valuation_two(value: Fraction) -> int:
    return valuation_two(value.numerator) - valuation_two(value.denominator)


def short_ratio_set() -> tuple[Fraction, ...]:
    values = {
        sign * Fraction(numerator, denominator)
        for numerator in range(1, 8)
        for denominator in range(1, 8)
        if numerator + denominator <= 7
        and gcd(numerator, denominator) == 1
        and (numerator - denominator) % 2
        for sign in (1, -1)
    }
    return tuple(sorted(values))


def ratio_graph_audit():
    short = short_ratio_set()
    require(len(short) == 24, short)
    vertices = tuple(sorted({Fraction(1), *short}))
    short_set = set(short)
    edges = tuple(
        (left, right)
        for left_index, left in enumerate(vertices)
        for right in vertices[left_index + 1 :]
        if left / right in short_set
    )
    adjacency = {vertex: set() for vertex in vertices}
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)

    # The colour is v_2(x) modulo three.  Every ratio in D has valuation
    # +-1 or +-2, so adjacent vertices always receive different colours.
    colours = {
        vertex: fraction_valuation_two(vertex) % 3
        for vertex in vertices
    }
    require(
        all(colours[left] != colours[right] for left, right in edges),
        "dyadic colouring",
    )
    colour_classes = tuple(
        tuple(vertex for vertex in vertices if colours[vertex] == colour)
        for colour in range(3)
    )
    require(tuple(len(part) for part in colour_classes) == (1, 12, 12), colour_classes)

    triangle_count = sum(
        1
        for first_index, first in enumerate(vertices)
        for second_index, second in enumerate(vertices[first_index + 1 :], first_index + 1)
        if second in adjacency[first]
        for third in vertices[second_index + 1 :]
        if third in adjacency[first] and third in adjacency[second]
    )
    require(len(edges) == 84 and triangle_count == 60, (len(edges), triangle_count))
    require(
        all(
            right in adjacency[left]
            for left, right in ((Fraction(1), Fraction(2)), (Fraction(1), Fraction(4)), (Fraction(2), Fraction(4)))
        ),
        "sharp triangle {1,2,4}",
    )

    # If u/v=d modulo p without the rational equality, p divides the reduced
    # numerator of u/v-d.  The maximum below is therefore a complete modular
    # stability threshold, not a sample over primes.
    false_numerators = tuple(
        abs((left / right - target).numerator)
        for left in vertices
        for right in vertices
        for target in short
        if left / right != target
    )
    maximum_false_numerator = max(false_numerators)
    require(maximum_false_numerator == 217, maximum_false_numerator)
    return (
        short,
        len(vertices),
        len(edges),
        triangle_count,
        colour_classes,
        maximum_false_numerator,
    )


def ratio_complement(prime: int, length: int) -> set[int]:
    odd_interval = {
        value % prime
        for value in range(-length, length + 1, 2)
    }
    ratios = {
        left * pow(right, -1, prime) % prime
        for left in odd_interval
        for right in odd_interval
    }
    return set(range(1, prime)) - ratios


def ratio_hostile_controls(short: tuple[Fraction, ...]):
    samples = ((547, 1, 8), (541, 9, 2), (557, 11, 4), (587, 13, 6))
    rows = []
    for prime, residue, offset in samples:
        require(is_prime(prime) and prime % 14 == residue, (prime, residue))
        length = (prime - offset) // 7
        require(length % 2 == 1 and prime == 7 * length + offset, (prime, length, offset))
        observed = ratio_complement(prime, length)
        expected = {
            value.numerator * pow(value.denominator, -1, prime) % prime
            for value in short
        }
        require(observed == expected and len(observed) == 24, (prime, observed ^ expected))
        rows.append((prime, residue, length, len(observed)))
    thresholds = tuple((offset, 448 + 8 * offset) for offset in (2, 4, 6, 8))
    require(max(value for _, value in thresholds) == 512, thresholds)
    return tuple(rows), thresholds


def finite_field_mask(base: set[int], coefficient: int, prime: int) -> int:
    inverse = pow(coefficient, -1, prime)
    result = 0
    for value in base:
        result |= 1 << ((inverse * value) % prime)
    return result


def block_bank(prime: int):
    k = (prime - 1) // 14
    even_base = {0, *range(1, k + 1), *range(prime - k, prime)}
    odd_ceiling = max(
        (value for value in range(1, prime, 2) if 7 * value < prime),
        default=-1,
    )
    odd_base = (
        {value % prime for value in range(-odd_ceiling, odd_ceiling + 1, 2)}
        if odd_ceiling >= 1
        else set()
    )
    representatives = range(1, (prime + 1) // 2)
    even = [
        (coefficient, finite_field_mask(even_base, coefficient, prime))
        for coefficient in representatives
    ]
    odd = [
        (coefficient, finite_field_mask(odd_base, coefficient, prime))
        for coefficient in representatives
    ]
    return even, odd, len(even_base), len(odd_base)


def find_seven_cover(prime: int, even_count: int):
    even, odd, even_size, odd_size = block_bank(prime)
    full = (1 << prime) - 1
    overlap_budget = even_count * even_size + (7 - even_count) * odd_size - prime
    if overlap_budget < 0:
        return None

    def choose_odd(start, left, union, overlap, chosen_even, chosen_odd):
        if left == 0:
            if union == full and overlap == overlap_budget:
                return chosen_even, chosen_odd
            return None
        if union.bit_count() + left * odd_size < prime:
            return None
        for index in range(start, len(odd) - left + 1):
            coefficient, candidate = odd[index]
            added_overlap = (candidate & union).bit_count()
            if overlap + added_overlap > overlap_budget:
                continue
            answer = choose_odd(
                index + 1,
                left - 1,
                union | candidate,
                overlap + added_overlap,
                chosen_even,
                chosen_odd + (coefficient,),
            )
            if answer is not None:
                return answer
        return None

    def choose_even(start, left, union, overlap, chosen):
        if left == 0:
            return choose_odd(0, 7 - even_count, union, overlap, chosen, ())
        for index in range(start, len(even) - left + 1):
            coefficient, candidate = even[index]
            added_overlap = (candidate & union).bit_count()
            if overlap + added_overlap > overlap_budget:
                continue
            answer = choose_even(
                index + 1,
                left - 1,
                union | candidate,
                overlap + added_overlap,
                chosen + (coefficient,),
            )
            if answer is not None:
                return answer
        return None

    # Every cover contains an even block because only those blocks contain the
    # reflection-fixed sheet.  Common multiplicative scaling normalizes it to
    # coefficient one.
    return choose_even(1, even_count - 1, even[0][1], 0, (1,))


def finite_bank():
    positives = []
    tested_profiles = 0
    tested_primes = 0
    for prime in range(3, 512):
        if not is_prime(prime):
            continue
        tested_primes += 1
        _, _, even_size, odd_size = block_bank(prime)
        for even_count in range(1, 8):
            overlap = even_count * even_size + (7 - even_count) * odd_size - prime
            if overlap < even_count - 1 or (overlap - (even_count - 1)) % 2:
                continue
            tested_profiles += 1
            answer = find_seven_cover(prime, even_count)
            if answer is not None:
                positives.append((prime, even_count, answer[0], answer[1]))
    support = tuple(sorted({prime for prime, _, _, _ in positives}))
    require(tested_primes == 96 and tested_profiles == 197, (tested_primes, tested_profiles))
    require(support == (11, 13, 23, 29), support)
    require(len(positives) == 12, positives)
    return tested_primes, tested_profiles, tuple(positives), support


def danger_mask(prime: int, residue: int) -> int:
    modulus = 2 * prime
    result = 0
    for sheet in range(prime):
        word = residue * (2 * sheet + 1) % modulus
        if 14 * min(word, modulus - word) < modulus:
            result |= 1 << sheet
    return result


def atom_record(prime: int, residues: tuple[int, ...]):
    masks = tuple(danger_mask(prime, residue) for residue in residues)
    full = (1 << prime) - 1
    union = 0
    for mask in masks:
        union |= mask
    require(union == full, (prime, residues))
    multiplicities = tuple(
        sum((mask >> sheet) & 1 for mask in masks)
        for sheet in range(prime)
    )
    return (
        prime,
        residues,
        tuple(mask.bit_count() for mask in masks),
        tuple((value, multiplicities.count(value)) for value in sorted(set(multiplicities))),
        sum(mask.bit_count() for mask in masks) - prime,
    )


def positive_atoms():
    rows = (
        atom_record(11, (1, 2, 3, 5, 7, 9)),
        atom_record(13, (1, 2, 3, 5, 7, 9, 11)),
        atom_record(23, (1, 4, 5, 7, 9, 11)),
        atom_record(29, (1, 5, 7, 8, 12, 13, 22)),
    )
    require(tuple(row[4] for row in rows) == (0, 0, 0, 2), rows)
    return rows


def even_prime_control():
    masks = tuple(danger_mask(2, residue) for residue in range(1, 4))
    require(masks == (0, 0, 0), masks)
    return 2, masks


def main() -> None:
    graph = ratio_graph_audit()
    hostiles = ratio_hostile_controls(graph[0])
    bank = finite_bank()
    atoms = positive_atoms()
    even_prime = even_prime_control()
    semantic_payload = (graph, hostiles, bank, atoms, even_prime)
    semantic = sha256(repr(semantic_payload).encode("ascii")).hexdigest()
    if EXPECTED_SEMANTIC_DIGEST is not None:
        require(semantic == EXPECTED_SEMANTIC_DIGEST, (semantic, EXPECTED_SEMANTIC_DIGEST))

    print("THM-3421 exact prime half-twist rank-seven companion")
    print(
        "ratio_graph="
        f"(D_size={len(graph[0])},vertices={graph[1]},edges={graph[2]},"
        f"triangles={graph[3]},colour_sizes={tuple(len(part) for part in graph[4])},"
        f"max_false_numerator={graph[5]})"
    )
    print(f"ratio_complement_hostiles_and_thresholds={hostiles}")
    print(
        "finite_bank="
        f"(primes={bank[0]},scalar_profiles={bank[1]},"
        f"positive_profiles={len(bank[2])},support={bank[3]})"
    )
    print(f"positive_scalar_profiles={bank[2]}")
    print(f"literal_positive_atoms={atoms}")
    print(f"even_prime_control={even_prime}")
    print(
        "scope=literal prime half-twist cap7 iff p in {11,13,23,29}; "
        "exact half ranks 6,7,6,7 respectively; composite rank7 and LRC14 remain open"
    )
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
