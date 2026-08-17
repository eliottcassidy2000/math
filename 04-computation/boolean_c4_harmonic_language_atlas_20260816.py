#!/usr/bin/env python3
"""Exact Boolean-C4, regular-language, and harmonic-height atlas.

This companion separates three structures that are easy to conflate:

* XOR on two membership bits is the four-state group (Z/2)^2;
* a tournament chooses one orientation on every unordered pair;
* reciprocal mass depends on the height used to inject a language into N.

Only integer and rational arithmetic is used.  The asymptotic conclusions in
the companion reflection are proved from the displayed counting identities
and elementary comparison tests, not inferred from floating-point data.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json


EXPECTED_SEMANTIC_SHA256 = "78ee6a60805603c034a58e93b04a99f8b927bc76d4de5b1715f11f4719ce5459"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def fibonacci(index: int) -> int:
    require(index >= 0, ("fibonacci index", index))
    left, right = 0, 1
    for _ in range(index):
        left, right = right, left + right
    return left


def binary_words(length: int) -> tuple[tuple[int, ...], ...]:
    require(length >= 0, ("word length", length))
    return tuple(
        tuple((mask >> position) & 1 for position in range(length))
        for mask in range(1 << length)
    )


def no_adjacent_ones(word: tuple[int, ...]) -> bool:
    return all(not (left == right == 1) for left, right in zip(word, word[1:]))


def zeckendorf_value(word: tuple[int, ...]) -> int:
    """Weights low-to-high bits by F_2,F_3,... ."""
    return sum(bit * fibonacci(position + 2) for position, bit in enumerate(word))


def radix_value(word: tuple[int, ...], base: int) -> int:
    value = 0
    for digit in reversed(word):
        require(0 <= digit < base, ("radix digit", word, base))
        value = base * value + digit
    return value


def leading_marker_code(word: tuple[int, ...], base: int) -> int:
    """Inject length-d words into [b^d,2b^d-1]."""
    require(base >= 2, ("base", base))
    return base ** len(word) + radix_value(word, base)


def xor_translation(vertex: int, displacement: int) -> int:
    return vertex ^ displacement


def cayley_arc_ledger(generators: tuple[int, ...]) -> dict[str, int]:
    vertices = range(4)
    arcs = {
        (vertex, xor_translation(vertex, generator))
        for vertex in vertices
        for generator in generators
    }
    require(all(left != right for left, right in arcs), ("Cayley loops", arcs))
    unordered = {tuple(sorted(arc)) for arc in arcs}
    bidirected = sum(
        1 for left, right in unordered
        if (left, right) in arcs and (right, left) in arcs
    )
    return {
        "directed_arcs": len(arcs),
        "unordered_pairs": len(unordered),
        "bidirected_pairs": bidirected,
        "missing_pairs": 6 - len(unordered),
    }


def tournament_arcs(size: int, mask: int) -> frozenset[tuple[int, int]]:
    pairs = tuple(
        (left, right)
        for left in range(size)
        for right in range(left + 1, size)
    )
    return frozenset(
        (right, left) if (mask >> index) & 1 else (left, right)
        for index, (left, right) in enumerate(pairs)
    )


def invariant_tournament_count(size: int, translations: tuple[tuple[int, ...], ...]) -> int:
    pair_count = size * (size - 1) // 2
    count = 0
    for mask in range(1 << pair_count):
        arcs = tournament_arcs(size, mask)
        if all(
            frozenset((permutation[left], permutation[right]) for left, right in arcs)
            == arcs
            for permutation in translations
        ):
            count += 1
    return count


def cyclic_translations(size: int) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple((vertex + shift) % size for vertex in range(size))
        for shift in range(size)
    )


def xor_translations() -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(vertex ^ shift for vertex in range(4))
        for shift in range(4)
    )


def harmonic_prefix(values: set[int], cutoff: int) -> Fraction:
    return sum((Fraction(1, value) for value in values if 1 <= value <= cutoff),
               Fraction(0))


def digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def main() -> None:
    # XOR is a four-cell membership partition.  Its natural Cayley graph has
    # bidirected involution edges, not tournament orientations.
    xor_square = cayley_arc_ledger((1, 2))
    xor_complete = cayley_arc_ledger((1, 2, 3))
    require(xor_square == {
        "directed_arcs": 8,
        "unordered_pairs": 4,
        "bidirected_pairs": 4,
        "missing_pairs": 2,
    }, xor_square)
    require(xor_complete == {
        "directed_arcs": 12,
        "unordered_pairs": 6,
        "bidirected_pairs": 6,
        "missing_pairs": 0,
    }, xor_complete)

    invariant_counts = {
        "V4_xor": invariant_tournament_count(4, xor_translations()),
        "C3": invariant_tournament_count(3, cyclic_translations(3)),
        "C4": invariant_tournament_count(4, cyclic_translations(4)),
        "C5": invariant_tournament_count(5, cyclic_translations(5)),
        "C6": invariant_tournament_count(6, cyclic_translations(6)),
    }
    require(invariant_counts == {
        "V4_xor": 0,
        "C3": 2,
        "C4": 0,
        "C5": 4,
        "C6": 0,
    }, invariant_counts)

    # The four Boolean atoms make the exact XOR harmonic identity pointwise.
    cutoff = 120
    set_a = {value for value in range(1, cutoff + 1) if value % 2 == 0}
    set_b = {value for value in range(1, cutoff + 1) if value % 3 == 0}
    cells = {
        "00": set(range(1, cutoff + 1)) - (set_a | set_b),
        "10": set_a - set_b,
        "01": set_b - set_a,
        "11": set_a & set_b,
    }
    require(set().union(*cells.values()) == set(range(1, cutoff + 1)),
            "Boolean cells cover prefix")
    require(sum(len(cell) for cell in cells.values()) == cutoff,
            "Boolean cells are disjoint")
    mass_a = harmonic_prefix(set_a, cutoff)
    mass_b = harmonic_prefix(set_b, cutoff)
    mass_intersection = harmonic_prefix(set_a & set_b, cutoff)
    mass_xor = harmonic_prefix(set_a ^ set_b, cutoff)
    require(mass_xor == mass_a + mass_b - 2 * mass_intersection,
            (mass_xor, mass_a, mass_b, mass_intersection))
    require(mass_xor
            == harmonic_prefix(cells["10"], cutoff)
            + harmonic_prefix(cells["01"], cutoff), "XOR cells")

    # The regular languages that produce Fibonacci and XOR counts.
    count_rows = []
    for length in range(13):
        words = binary_words(length)
        no_11 = tuple(word for word in words if no_adjacent_ones(word))
        even_xor = tuple(word for word in words if sum(word) % 2 == 0)
        require(len(no_11) == fibonacci(length + 2),
                ("Fibonacci language", length, len(no_11)))
        expected_even = 1 if length == 0 else 2 ** (length - 1)
        require(len(even_xor) == expected_even,
                ("XOR language", length, len(even_xor)))

        # Under Fibonacci weights the same no-11 language is a complete
        # interval, whereas its ordinary binary shell has sub-full growth.
        zeckendorf_values = {zeckendorf_value(word) for word in no_11}
        require(zeckendorf_values == set(range(fibonacci(length + 2))),
                ("Zeckendorf interval", length, zeckendorf_values))

        binary_codes = tuple(leading_marker_code(word, 2) for word in no_11)
        require(len(set(binary_codes)) == len(no_11),
                ("binary marker injective", length))
        count_rows.append({
            "length": length,
            "no11": len(no_11),
            "even_xor": len(even_xor),
            "zeckendorf_interval_size": len(zeckendorf_values),
        })

    # Exact normalized count series.  These are upper bounds for the
    # corresponding marker-coded reciprocal masses, and the shell lemma gives
    # the reverse comparison up to a factor two.
    no11_normalized_partial = sum(
        (Fraction(fibonacci(length + 2), 2 ** length) for length in range(13)),
        Fraction(0),
    )
    require(no11_normalized_partial < 6, no11_normalized_partial)
    require(
        sum((Fraction(fibonacci(length + 2), 2 ** length)
             for length in range(80)), Fraction(0)) < 6,
        "finite no11 partials lie below generating-function value 6",
    )
    for length in range(1, 13):
        require(Fraction(2 ** (length - 1), 2 ** length) == Fraction(1, 2),
                ("XOR full entropy", length))

    # The shell comparison is checked exactly on all words through these
    # depths; its proof is simply b^d <= code < 2b^d.
    shell_rows = []
    for base, max_depth in ((2, 10), (3, 7), (4, 5), (6, 4)):
        for depth in range(max_depth + 1):
            word_count = base ** depth
            mass = sum(
                (Fraction(1, base ** depth + value) for value in range(word_count)),
                Fraction(0),
            )
            normalized_count = Fraction(word_count, base ** depth)
            require(normalized_count / 2 <= mass <= normalized_count,
                    ("shell comparison", base, depth, mass))
        shell_rows.append({
            "base": base,
            "checked_through_depth": max_depth,
            "full_language_normalized_count_each_shell": "1",
        })

    # Height dependence on the same Berggren U-spine and on the Fibonacci ray.
    height_rows = []
    for depth in range(12):
        farey_denominator = depth + 2
        hypotenuse = 2 * depth * depth + 6 * depth + 5
        q_height = (2 * depth + 3) ** 2 + 2
        ternary_marker_floor = 3 ** depth
        require(hypotenuse >= 2 * (depth + 1) ** 2,
                ("quadratic hypotenuse bound", depth, hypotenuse))
        require(q_height >= 4 * (depth + 1) ** 2,
                ("quadratic Q bound", depth, q_height))
        height_rows.append({
            "depth": depth,
            "farey_denominator": farey_denominator,
            "hypotenuse": hypotenuse,
            "Q": q_height,
            "ternary_marker_floor": ternary_marker_floor,
        })
    for depth in range(1, 11):
        denominator = fibonacci(2 * depth + 1)
        next_denominator = fibonacci(2 * depth + 3)
        require(next_denominator >= 2 * denominator,
                ("Fibonacci denominator growth", depth))

    # Every loopless m-state walk language is harmonic-thin under m-radix
    # coding: a_d <= m(m-1)^d.  C4/C6 successor cycles have the sharper
    # exact series displayed here.
    walk_bounds = {
        "loopless_C4_upper_normalized_sum": "4",
        "loopless_C6_upper_normalized_sum": "6",
        "successor_C4_exact_normalized_sum": "4/3",
        "successor_C6_exact_normalized_sum": "6/5",
    }
    require(sum((Fraction(4, 4 ** length) for length in range(1, 40)),
                Fraction(0)) < Fraction(4, 3), "C4 finite geometric partial")
    require(sum((Fraction(6, 6 ** length) for length in range(1, 40)),
                Fraction(0)) < Fraction(6, 5), "C6 finite geometric partial")

    ledger = {
        "xor_cayley_square": xor_square,
        "xor_cayley_complete": xor_complete,
        "invariant_tournaments": invariant_counts,
        "harmonic_xor_prefix": {
            "cutoff": cutoff,
            "cell_sizes": {key: len(value) for key, value in cells.items()},
            "identity": "H(A xor B)=H(A)+H(B)-2H(A intersect B)",
            "mass_xor": str(mass_xor),
        },
        "language_counts": count_rows,
        "no11_normalized_limit": "6",
        "xor_normalized_shell_mass_for_positive_lengths": "1/2",
        "shell_rows": shell_rows,
        "height_rows": height_rows,
        "walk_bounds": walk_bounds,
    }
    semantic_sha256 = digest(ledger)
    require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest", semantic_sha256, EXPECTED_SEMANTIC_SHA256))

    print("== Boolean C4 / harmonic regular-language atlas ==")
    print(f"xor_two_generators={xor_square}")
    print(f"xor_three_generators={xor_complete}")
    print(f"translation_invariant_tournaments={invariant_counts}")
    print("harmonic_xor_identity=H_N(A xor B)=H_N(A)+H_N(B)-"
          f"2H_N(A intersect B);cutoff={cutoff};mass={mass_xor}")
    print("regular_counts="
          + str(tuple((row["length"], row["no11"], row["even_xor"])
                      for row in count_rows)))
    print("radix_height=no11 has rho=phi<2 and normalized sum 6: convergent;"
          "even-XOR has rho=2 and shell ratio 1/2: divergent")
    print("fibonacci_height=same no11 words fill [0,F_(d+2));union=N: divergent")
    print("berggren_height=Farey denominators n+2 divergent;"
          "hypotenuse,Q,and ternary marker heights convergent")
    print(f"walk_language_bounds={walk_bounds}")
    print(f"semantic_sha256={semantic_sha256}")
    print("status=VERIFIED-EXACT FINITE IDENTITIES + PROVED COMPARISON CRITERIA;"
          "no tournament equivalence, Keller bridge, current, or LRC14 claim")


if __name__ == "__main__":
    main()
