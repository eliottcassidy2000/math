#!/usr/bin/env python3
"""Exact companion for THM-3382.

The same flipped Fibonacci/Berggren ray is embedded in the natural numbers in
two canonical ways.  Fibonacci time gives an arithmetic progression; ternary
breadth-first heap address gives an affine-base-nine recurrence.  This script
checks the address formulas, the harmonic bifurcation, and the full-tree
automatic hostile without floating point or optimization-sensitive asserts.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
from pathlib import Path


EXPECTED_SEMANTIC_DIGEST = "67e8bf54060f0269d1d47e9b22c0928da7d0c1377ee0deba55f95d98be7a0b6a"
RAY_CHECK_MAX = 800
TREE_MAX_DEPTH = 11
LAMBERT_TERMS = 12
RAY_MAHLER_TERMS = 20

PINNED_THEOREMS = (
    (
        "THM-3359",
        "01-canon/theorems/THM-3359-modular-c-finite-supports-harmonic-density-and-periodic-scar.md",
        "0ff08a37801c5f3b96fa3ee74b0d7a897b3f87c5db4397b4abdadbe99f38562a",
    ),
    (
        "THM-3379",
        "01-canon/theorems/THM-3379-fibonacci-ray-local-t4-bit-is-mod3-boolean-projection.md",
        "d1c9cb25ca02ef77d5e1ff681972554856247ad9326fa2b9d0a480586e3d499b",
    ),
)


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(Path(path).read_bytes().replace(b"\r\n", b"\n")).hexdigest()


class ExactDigest:
    def __init__(self):
        self._hash = sha256()

    def add(self, item):
        self._hash.update(repr(item).encode("ascii"))
        self._hash.update(b"\n")

    def hexdigest(self):
        return self._hash.hexdigest()


def heap_child(address, letter):
    offsets = {"A": 1, "B": 2, "C": 3}
    require(letter in offsets, ("unknown branch", letter))
    return 3 * address + offsets[letter]


def heap_word(word):
    address = 0
    for letter in word:
        address = heap_child(address, letter)
    return address


def epsilon_word(word):
    bit = 0
    for letter in word:
        if letter == "A":
            bit = 0
        elif letter == "B":
            bit = 1 - bit
        elif letter == "C":
            bit = 1
        else:
            raise RuntimeError(("unknown branch", letter))
    return bit


def ray_records(index):
    power = 9**index
    return (
        (
            "R2",
            "BA" * index,
            3 * index + 2,
            7 * (power - 1) // 8,
            0,
        ),
        (
            "R0",
            "A" + "BA" * index,
            3 * index + 3,
            (15 * power - 7) // 8,
            0,
        ),
        (
            "R1",
            "C" + "BC" * index,
            3 * index + 4,
            (33 * power - 9) // 8,
            1,
        ),
    )


def ray_address(ray, index):
    return {record[0]: record[3] for record in ray_records(index)}[ray]


def flipped_heap_address(index):
    return ray_address("R1", index)


def level_interval(depth):
    return (3**depth - 1) // 2, (3 ** (depth + 1) - 3) // 2


def floor_log_three(number):
    require(number >= 1, ("positive logarithm input", number))
    exponent = 0
    power = 1
    while 3 * power <= number:
        power *= 3
        exponent += 1
    return exponent


def base_three(number):
    require(number >= 0, ("base-three input", number))
    if number == 0:
        return "0"
    digits = []
    while number:
        digits.append(str(number % 3))
        number //= 3
    return "".join(reversed(digits))


def fraction_key(value):
    return value.numerator, value.denominator


def main():
    repo_root = Path(__file__).resolve().parents[1]
    source_path = Path(__file__).resolve()
    source_hash = lf_hash(source_path)
    theorem_hashes = []
    for theorem_id, relative_path, expected_hash in PINNED_THEOREMS:
        actual_hash = lf_hash(repo_root / relative_path)
        require(
            actual_hash == expected_hash,
            ("theorem dependency drift", theorem_id, actual_hash, expected_hash),
        )
        theorem_hashes.append((theorem_id, actual_hash))

    syntax = ast.parse(source_path.read_text(encoding="utf-8"))
    assertion_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax))
    floating_literals = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(syntax)
    )
    require(assertion_nodes == 0, ("assertion nodes", assertion_nodes))
    require(floating_literals == 0, ("floating literals", floating_literals))

    digest = ExactDigest()
    digest.add(tuple(theorem_hashes))

    # The heap map is a bijective shortlex enumeration of all finite words.
    level = (0,)
    level_table = []
    heap_bijection_checks = 0
    for depth in range(TREE_MAX_DEPTH + 1):
        left, right = level_interval(depth)
        require(level == tuple(range(left, right + 1)), (
            "heap level interval", depth, left, right, level[:5], level[-5:]
        ))
        level_table.append((depth, left, right, len(level)))
        heap_bijection_checks += len(level)
        if depth < TREE_MAX_DEPTH:
            level = tuple(
                heap_child(address, letter)
                for address in level
                for letter in "ABC"
            )

    # All three exact ray formulas, affine recurrences, and C-finite lifts.
    first_ray_records = []
    last_ray_records = []
    ray_formula_checks = 0
    previous = {}
    before_previous = {}
    address_language_checks = 0
    constants = {"R2": 7, "R0": 7, "R1": 9}
    initials = {"R2": 0, "R0": 1, "R1": 3}
    for index in range(RAY_CHECK_MAX + 1):
        records = ray_records(index)
        addresses = []
        for ray, word, time_index, address, expected_bit in records:
            require(heap_word(word) == address, (
                "ray address formula", index, ray, word, heap_word(word), address
            ))
            require(epsilon_word(word) == expected_bit, (
                "ray epsilon", index, ray, word, epsilon_word(word), expected_bit
            ))
            require(time_index % 3 == {"R2": 2, "R0": 0, "R1": 1}[ray], (
                "time residue", index, ray, time_index
            ))
            require(address >= 0, ("negative heap address", index, ray, address))
            if index == 0:
                require(address == initials[ray], ("ray initial", ray, address))
            else:
                require(address == 9 * previous[ray] + constants[ray], (
                    "affine-nine recurrence", index, ray, address, previous[ray]
                ))
            if index >= 2:
                require(address == 10 * previous[ray] - 9 * before_previous[ray], (
                    "C-finite recurrence", index, ray, address
                ))
            before_previous[ray] = previous.get(ray, address)
            previous[ray] = address
            addresses.append(address)
            ray_formula_checks += 1
        require(len(set(addresses)) == 3, ("same-index ray collision", index, addresses))
        require(
            addresses[0] < addresses[1] < addresses[2] < ray_address("R2", index + 1),
            ("strict ray interlacing", index, addresses, ray_address("R2", index + 1)),
        )
        if index < 5:
            first_ray_records.append(tuple((r[0], r[2], r[3], r[4]) for r in records))
        if index > RAY_CHECK_MAX - 3:
            last_ray_records.append(tuple((r[0], r[2], r[3], r[4]) for r in records))
        flipped_language = "10" if index == 0 else "11" + "01" * (index - 1) + "00"
        require(base_three(addresses[2]) == flipped_language, (
            "flipped base-three language", index, base_three(addresses[2]), flipped_language
        ))
        address_language_checks += 1

    # The flipped fibre has arithmetic time support but lacunary heap support.
    flipped_pairs = tuple(
        (3 * index + 4, flipped_heap_address(index)) for index in range(20)
    )
    require(len({pair[0] for pair in flipped_pairs}) == len(flipped_pairs), "time injectivity")
    require(len({pair[1] for pair in flipped_pairs}) == len(flipped_pairs), "heap injectivity")
    require(all(
        flipped_pairs[index + 1][1] == 9 * flipped_pairs[index][1] + 9
        for index in range(len(flipped_pairs) - 1)
    ), "flipped affine-nine law")
    require(Fraction(1, flipped_pairs[0][1]) > Fraction(1, flipped_pairs[0][0]), (
        "finite-prefix hostile", flipped_pairs[0]
    ))
    require(heap_word("BA") == 7 and heap_word("AB") == 5, "word-direction hostile")

    # Rational OGFs: denominator (1-t)(1-9t)=1-10t+9t^2.
    ogf_numerators = {"R2": (0, 7), "R0": (1, 6), "R1": (3, 6)}
    ogf_checks = 0
    for ray in ("R2", "R0", "R1"):
        values = [ray_address(ray, index) for index in range(RAY_CHECK_MAX + 1)]
        for index, value in enumerate(values):
            coefficient = value
            if index >= 1:
                coefficient -= 10 * values[index - 1]
            if index >= 2:
                coefficient += 9 * values[index - 2]
            expected = ogf_numerators[ray][index] if index < 2 else 0
            require(coefficient == expected, (
                "rational OGF coefficient", ray, index, coefficient, expected
            ))
            ogf_checks += 1

    # Sparse ray indicators obey exact base-nine Mahler set equations.
    ray_mahler_addresses = {
        ray: tuple(ray_address(ray, index) for index in range(RAY_MAHLER_TERMS + 1))
        for ray in ("R2", "R0", "R1")
    }
    require(ray_mahler_addresses["R2"][1:] == tuple(
        [7] + [7 + 9 * value for value in ray_mahler_addresses["R2"][1:-1]]
    ), "R2 sparse Mahler equation")
    require(ray_mahler_addresses["R0"] == tuple(
        [1] + [7 + 9 * value for value in ray_mahler_addresses["R0"][:-1]]
    ), "R0 sparse Mahler equation")
    require(ray_mahler_addresses["R1"] == tuple(
        [3] + [9 + 9 * value for value in ray_mahler_addresses["R1"][:-1]]
    ), "R1 sparse Mahler equation")

    # Exact Lambert-series identity and rational tail enclosure.
    q = Fraction(1, 9)
    a = Fraction(3, 11)
    heap_harmonic_partial = Fraction(0)
    lambert_partial = Fraction(0)
    for index in range(LAMBERT_TERMS):
        address = flipped_heap_address(index)
        direct_term = Fraction(1, address)
        lambert_term = Fraction(8, 33) * q**index / (1 - a * q**index)
        require(direct_term == lambert_term, (
            "Lambert term", index, direct_term, lambert_term
        ))
        heap_harmonic_partial += direct_term
        lambert_partial += lambert_term
    require(heap_harmonic_partial == lambert_partial, "Lambert partial sum")
    lower_tail = Fraction(3, 11 * 9**LAMBERT_TERMS)
    upper_tail = lower_tail / (1 - a * q**LAMBERT_TERMS)
    require(lower_tail < upper_tail, "Lambert tail ordering")
    for index in range(LAMBERT_TERMS, LAMBERT_TERMS + 20):
        term = Fraction(1, flipped_heap_address(index))
        require(Fraction(8, 33 * 9**index) < term, ("tail lower term", index))
        require(term <= Fraction(8, 33) * q**index / (1 - a * q**LAMBERT_TERMS), (
            "tail upper term", index
        ))
    harmonic_bracket = (
        heap_harmonic_partial + lower_tail,
        heap_harmonic_partial + upper_tail,
    )
    require(harmonic_bracket[0] < harmonic_bracket[1], "harmonic bracket")
    dirichlet_term_checks = 0
    for exponent in range(1, 5):
        for index in range(20):
            direct = Fraction(1, flipped_heap_address(index) ** exponent)
            factored = (
                Fraction(8, 33) ** exponent
                * q ** (index * exponent)
                / (1 - a * q**index) ** exponent
            )
            require(direct == factored, (
                "Dirichlet factored term", exponent, index, direct, factored
            ))
            dirichlet_term_checks += 1

    # Full-tree hostile: the heap-indexed local-T4 bit is automatic and dense.
    tree_max_address = level_interval(TREE_MAX_DEPTH)[1]
    epsilon = [0] * (tree_max_address + 1)
    automatic_checks = 0
    for address in range(1, tree_max_address + 1):
        residue = address % 3
        if residue == 1:
            bit = 0
        elif residue == 2:
            bit = 1 - epsilon[(address - 2) // 3]
        else:
            bit = 1
        epsilon[address] = bit
        automatic_checks += 1

    level_bit_counts = []
    for depth in range(TREE_MAX_DEPTH + 1):
        left, right = level_interval(depth)
        ones = sum(epsilon[left : right + 1])
        expected_ones = (3**depth - (-1) ** depth) // 2
        require(ones == expected_ones, (
            "level bit count", depth, ones, expected_ones
        ))
        level_bit_counts.append((depth, ones, right - left + 1 - ones))

    discrepancy = 0
    prefix_discrepancies = [0]
    minimum_discrepancy = (0, 0)
    maximum_discrepancy = (0, 0)
    for address in range(1, tree_max_address + 1):
        discrepancy += 2 * epsilon[address] - 1
        prefix_discrepancies.append(discrepancy)
        require(abs(discrepancy) <= 1 + floor_log_three(address), (
            "prefix discrepancy bound", address, discrepancy
        ))
        if discrepancy < minimum_discrepancy[1]:
            minimum_discrepancy = (address, discrepancy)
        if discrepancy > maximum_discrepancy[1]:
            maximum_discrepancy = (address, discrepancy)

    prefix_recurrence_checks = 0
    for multiplier in range(1, tree_max_address // 3 + 1):
        require(
            prefix_discrepancies[3 * multiplier]
            == 1 - prefix_discrepancies[multiplier - 1],
            ("prefix recurrence 3m", multiplier),
        )
        prefix_recurrence_checks += 1
    for multiplier in range(1, (tree_max_address - 1) // 3 + 1):
        require(
            prefix_discrepancies[3 * multiplier + 1]
            == -prefix_discrepancies[multiplier - 1],
            ("prefix recurrence 3m+1", multiplier),
        )
        prefix_recurrence_checks += 1
    for multiplier in range((tree_max_address - 2) // 3 + 1):
        require(
            prefix_discrepancies[3 * multiplier + 2]
            == -prefix_discrepancies[multiplier],
            ("prefix recurrence 3m+2", multiplier),
        )
        prefix_recurrence_checks += 1

    # Coefficientwise check of E(z)=z^2(1+z)/(1-z^3)-z^2 E(z^3).
    mahler_checks = 0
    for exponent in range(tree_max_address + 1):
        rational_part = int(
            exponent >= 2 and exponent % 3 in (0, 2)
        )
        recursive_part = (
            epsilon[(exponent - 2) // 3]
            if exponent >= 2 and exponent % 3 == 2
            else 0
        )
        require(epsilon[exponent] == rational_part - recursive_part, (
            "Mahler coefficient", exponent, epsilon[exponent],
            rational_part, recursive_part
        ))
        mahler_checks += 1

    # The full 3-kernel closes on exactly six explicit states.
    def state_value(state, index):
        if state == "E":
            return epsilon[index]
        if state == "C":
            return 1 - epsilon[index]
        if state == "U":
            return int(index != 0)
        if state == "D":
            return int(index == 0)
        if state == "Z":
            return 0
        if state == "O":
            return 1
        raise RuntimeError(("unknown kernel state", state))

    kernel_transitions = {
        "E": ("U", "Z", "C"),
        "C": ("D", "O", "E"),
        "U": ("U", "O", "O"),
        "D": ("D", "Z", "Z"),
        "Z": ("Z", "Z", "Z"),
        "O": ("O", "O", "O"),
    }
    kernel_checks = 0
    for state, targets in kernel_transitions.items():
        for residue, target in enumerate(targets):
            for index in range((tree_max_address - residue) // 3 + 1):
                require(
                    state_value(state, 3 * index + residue) == state_value(target, index),
                    ("3-kernel transition", state, residue, target, index),
                )
                kernel_checks += 1
    require(len({
        tuple(state_value(state, index) for index in range(20))
        for state in kernel_transitions
    }) == 6, "six distinct 3-kernel states")

    depth_support = tuple(2 * index + 1 for index in range(20))
    require(depth_support == tuple(len("C" + "BC" * index) for index in range(20)), (
        "R1 depth support", depth_support
    ))
    require(depth_support == tuple(len("A" + "BA" * index) for index in range(20)), (
        "R0 depth collision", depth_support
    ))

    digest.add(tuple(level_table))
    digest.add(tuple(first_ray_records))
    digest.add(tuple(last_ray_records))
    digest.add(tuple(flipped_pairs))
    digest.add(tuple((ray, ogf_numerators[ray]) for ray in ("R2", "R0", "R1")))
    digest.add(tuple((ray, ray_mahler_addresses[ray]) for ray in ("R2", "R0", "R1")))
    digest.add(fraction_key(heap_harmonic_partial))
    digest.add(tuple(fraction_key(value) for value in harmonic_bracket))
    digest.add(tuple(level_bit_counts))
    digest.add((minimum_discrepancy, maximum_discrepancy))
    digest.add(tuple((state, kernel_transitions[state]) for state in kernel_transitions))
    digest.add(depth_support)
    digest.add((
        heap_bijection_checks, ray_formula_checks, address_language_checks,
        ogf_checks, dirichlet_term_checks, automatic_checks,
        prefix_recurrence_checks, mahler_checks, kernel_checks,
    ))
    semantic_digest = digest.hexdigest()
    if EXPECTED_SEMANTIC_DIGEST != "TO_BE_FROZEN":
        require(semantic_digest == EXPECTED_SEMANTIC_DIGEST, (
            "semantic digest drift", semantic_digest, EXPECTED_SEMANTIC_DIGEST
        ))

    print("THM-3382 FIBONACCI-RAY DUAL-INDEX HARMONIC BIFURCATION")
    print("source_sha256_lf", source_hash)
    print("theorem_source_sha256_lf", tuple(theorem_hashes))
    print("assertion_nodes", assertion_nodes, "floating_literals", floating_literals)
    print("heap_definition h(empty)=0; h(wA)=3h(w)+1; h(wB)=3h(w)+2; h(wC)=3h(w)+3")
    print("heap_bijection_depths", (0, TREE_MAX_DEPTH), "checks", heap_bijection_checks)
    print("heap_level_boundary", level_table[-1])
    print("ray_address_formulas", (
        "R2=7(9^r-1)/8",
        "R0=(15*9^r-7)/8",
        "R1=(33*9^r-9)/8",
    ))
    print("ray_affine_recursions", (
        "R2(r+1)=9R2(r)+7",
        "R0(r+1)=9R0(r)+7",
        "R1(r+1)=9R1(r)+9",
    ))
    print("common_C_finite_recurrence x(r+2)=10x(r+1)-9x(r)")
    print("strict_interlacing", "R2(r)<R0(r)<R1(r)<R2(r+1)")
    print("ray_rational_OGFs", (
        "R2:7t/((1-t)(1-9t))",
        "R0:(1+6t)/((1-t)(1-9t))",
        "R1:(3+6t)/((1-t)(1-9t))",
    ))
    print("ray_sparse_Mahler_equations", (
        "G2(z)=z^7+z^7G2(z^9) (root excluded)",
        "G0(z)=z+z^7G0(z^9)",
        "G1(z)=z^3+z^9G1(z^9)",
    ))
    print("ray_formula_range", (0, RAY_CHECK_MAX), "checks", ray_formula_checks)
    print("flipped_base_three_language", "{10} union 11(01)*00", "checks", address_language_checks)
    print("ray_initial_records", tuple(first_ray_records))
    print("flipped_dual_index_first20", flipped_pairs)
    print("time_support", "{3r+4:r>=0}", "density=1/3", "harmonic=DIVERGENT")
    print("heap_support", "{(33*9^r-9)/8:r>=0}", "density=0", "harmonic=CONVERGENT")
    print("heap_harmonic_Lambert", "(8/33) sum_(r>=0) q^r/(1-aq^r)", "q=1/9", "a=3/11")
    print("heap_harmonic_terms", LAMBERT_TERMS)
    print("heap_harmonic_partial", fraction_key(heap_harmonic_partial))
    print("heap_harmonic_tail_bounds", fraction_key(lower_tail), fraction_key(upper_tail))
    print("heap_harmonic_bracket", tuple(fraction_key(value) for value in harmonic_bracket))
    print("heap_Dirichlet_profile", "(8/33)^s sum_(m>=0) (s)_m/m! (3/11)^m/(1-9^(-(s+m))) for Re(s)>0")
    print("heap_Dirichlet_factored_term_checks", dirichlet_term_checks)
    print("full_tree_automatic_rules", (
        "e(3n+1)=0", "e(3n+2)=1-e(n)", "e(3n+3)=1"
    ))
    print("full_tree_level_counts_epsilon0_epsilon1", tuple(
        (depth, zeros, ones) for depth, ones, zeros in level_bit_counts
    ))
    print("full_tree_prefix_discrepancy_extrema", minimum_discrepancy, maximum_discrepancy)
    print("full_tree_prefix_recurrence_checks", prefix_recurrence_checks)
    print("full_tree_natural_density", "1/2")
    print("full_tree_harmonic_profile", "(1/2)log(x)+C_tree+O(log(x)/x)")
    print("full_tree_Mahler_equation", "E(z)=z^2(1+z)/(1-z^3)-z^2E(z^3)")
    print("full_tree_3_kernel", tuple(kernel_transitions.items()), "states=6")
    print("depth_carrier_hostile", "R0 and R1 both map to {2r+1}; depth loses the tournament bit")
    print("automatic_checks", automatic_checks, "Mahler_checks", mahler_checks, "kernel_checks", kernel_checks)
    print("hostile_conclusion", "heap indexing does not make the full tree sparse; the bifurcation is ray-and-encoding specific")
    print("typing_conclusion", "C-finite address values are not modular index support; harmonic class needs the index map sidecar")
    print("semantic_sha256", semantic_digest)
    print("status=ALL EXACT CHECKS PASSED; no encoding-invariant density claim")


if __name__ == "__main__":
    main()
