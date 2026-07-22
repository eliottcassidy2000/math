#!/usr/bin/env python3
"""Optimization-safe exact referee for THM-2069 and the THM-211 bridge."""

from collections import Counter
from itertools import combinations, permutations, product
from math import comb, gcd, lcm


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def gcd_many(values):
    answer = 0
    for value in values:
        answer = gcd(answer, abs(value))
    return answer


def determinant(matrix):
    size = len(matrix)
    require(all(len(row) == size for row in matrix), ("nonsquare", matrix))
    answer = 0
    for permutation in permutations(range(size)):
        inversions = sum(
            permutation[left] > permutation[right]
            for left in range(size)
            for right in range(left + 1, size)
        )
        term = 1
        for row, column in enumerate(permutation):
            term *= matrix[row][column]
        answer += (-1 if inversions % 2 else 1) * term
    return answer


def rank_mod(rows, prime):
    if not rows:
        return 0
    matrix = [[entry % prime for entry in row] for row in rows]
    row_count = len(matrix)
    column_count = len(matrix[0])
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, row_count) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        inverse = pow(matrix[pivot_row][column], -1, prime)
        matrix[pivot_row] = [(inverse * value) % prime for value in matrix[pivot_row]]
        for row in range(row_count):
            if row == pivot_row or not matrix[row][column]:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                (left - factor * right) % prime
                for left, right in zip(matrix[row], matrix[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def lattice_index(rows, rank):
    return gcd_many(
        determinant(tuple(rows[index] for index in chosen))
        for chosen in combinations(range(len(rows)), rank)
    )


def prime_factors(value):
    value = abs(value)
    answer = set()
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            answer.add(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 1
    if value > 1:
        answer.add(value)
    return answer


def radical(values):
    answer = 1
    for value in values:
        for prime in prime_factors(value):
            answer *= prime if answer % prime else 1
    return answer


def evaluate(rows, covector, prime):
    return tuple(
        sum(left * right for left, right in zip(row, covector)) % prime
        for row in rows
    )


def normalize_projective(vector, prime):
    first = next(entry for entry in vector if entry % prime)
    inverse = pow(first, -1, prime)
    return tuple((inverse * entry) % prime for entry in vector)


def code_data(rows, prime):
    rank = len(rows[0])
    parameters = list(product(range(prime), repeat=rank))
    words = {parameter: evaluate(rows, parameter, prime) for parameter in parameters}
    require(len(set(words.values())) == prime ** rank,
            ("evaluation not injective", rows, prime))
    distribution = Counter(sum(entry != 0 for entry in word) for word in words.values())
    return words, distribution


def direct_deletion_bad(word, deletion_size):
    length = len(word)
    return any(
        all(word[index] == 0 for index in range(length) if index not in deleted)
        for deleted in combinations(range(length), deletion_size)
    )


def cocircuit_supports(rows, prime):
    rank = len(rows[0])
    length = len(rows)
    answer = set()
    for size in range(1, length + 1):
        for support_tuple in combinations(range(length), size):
            support = frozenset(support_tuple)
            complement = [rows[index] for index in range(length) if index not in support]
            if rank_mod(complement, prime) == rank:
                continue
            minimal = True
            for removed in support:
                smaller = support - {removed}
                smaller_complement = [
                    rows[index] for index in range(length) if index not in smaller
                ]
                if rank_mod(smaller_complement, prime) < rank:
                    minimal = False
                    break
            if minimal:
                answer.add(support)
    return answer


def template_bank():
    answer = []
    for rank, target in ((2, 15), (3, 35)):
        basis = [tuple(int(i == j) for i in range(rank)) for j in range(rank)]
        candidates = [
            vector
            for vector in product((-1, 0, 1), repeat=rank)
            if any(vector) and vector not in basis
        ]
        for extras in combinations(candidates, 2):
            rows = tuple(basis + list(extras))
            if rows not in answer:
                answer.append(rows)
            if sum(len(template[0]) == rank for template in answer) == target:
                break
    return answer


def audit_local_code_theorem():
    templates = template_bank()
    require(len(templates) == 50, ("template count", len(templates)))
    local_k_checks = 0
    projective_checks = 0
    cocircuit_checks = 0
    crt_prime_checks = 0
    crt_residue_checks = 0

    for rows in templates:
        rank = len(rows[0])
        length = len(rows)
        require(lattice_index(rows, rank) == 1, ("not integral spanning", rows))
        for prime in (2, 3, 5):
            words, distribution = code_data(rows, prime)
            nonzero_supports = {
                frozenset(index for index, entry in enumerate(word) if entry)
                for parameter, word in words.items()
                if any(parameter)
            }
            minimal_supports = {
                support
                for support in nonzero_supports
                if not any(other < support for other in nonzero_supports)
            }
            cocircuits = cocircuit_supports(rows, prime)
            require(minimal_supports == cocircuits,
                    ("cocircuit mismatch", rows, prime, minimal_supports, cocircuits))
            minimum_weight = min(map(len, cocircuits))
            require(min(weight for weight in distribution if weight) == minimum_weight,
                    ("cogirth mismatch", rows, prime))
            coloop_count = sum(
                rank_mod([row for index, row in enumerate(rows) if index != deleted], prime)
                < rank
                for deleted in range(length)
            )
            require(distribution.get(1, 0) == (prime - 1) * coloop_count,
                    ("coloop weight-one mismatch", rows, prime))
            require(
                (distribution.get(1, 0) + distribution.get(2, 0) == 0)
                == (not any(len(cocircuit) <= 2 for cocircuit in cocircuits)),
                ("cosimplicity mismatch", rows, prime),
            )
            cocircuit_checks += 1

            for deletion_size in range(length + 1):
                bad_parameters = {
                    parameter
                    for parameter, word in words.items()
                    if any(parameter) and direct_deletion_bad(word, deletion_size)
                }
                predicted_bad = sum(
                    distribution.get(weight, 0)
                    for weight in range(1, deletion_size + 1)
                )
                require(len(bad_parameters) == predicted_bad,
                        ("low-weight count", rows, prime, deletion_size))
                require(all(
                    direct_deletion_bad(words[parameter], deletion_size)
                    == (sum(entry != 0 for entry in words[parameter]) <= deletion_size)
                    for parameter in words if any(parameter)
                ), ("deletion equivalence", rows, prime, deletion_size))
                bad_projective = {
                    normalize_projective(parameter, prime)
                    for parameter in bad_parameters
                }
                require(len(bad_projective) * (prime - 1) == predicted_bad,
                        ("projective quotient", rows, prime, deletion_size))
                require((predicted_bad == 0) == (deletion_size < minimum_weight),
                        ("first failure radius", rows, prime, deletion_size))
                local_k_checks += 1
                projective_checks += 1

        for deletion_size in range(length):
            indices = []
            for deleted in combinations(range(length), deletion_size):
                retained = [row for index, row in enumerate(rows) if index not in deleted]
                indices.append(lattice_index(retained, rank))
            if not all(indices):
                for prime in (2, 3, 5, 7):
                    _, distribution = code_data(rows, prime)
                    require(sum(distribution.get(w, 0) for w in range(1, deletion_size + 1)) > 0,
                            ("rank-deficient wheel should be bad", rows, deletion_size, prime))
                continue

            modulus = radical(indices)
            test_primes = {2, 3, 5, 7}
            test_primes.update(prime_factors(modulus))
            local_good = {}
            for prime in sorted(test_primes):
                _, distribution = code_data(rows, prime)
                bad = sum(distribution.get(w, 0) for w in range(1, deletion_size + 1))
                require((bad > 0) == (modulus % prime == 0),
                        ("index-prime equivalence", rows, deletion_size, prime, modulus))
                if modulus % prime == 0:
                    local_good[prime] = prime ** rank - 1 - bad
                crt_prime_checks += 1

            if modulus <= 15:
                direct = 0
                for residue in product(range(modulus), repeat=rank):
                    allowed = True
                    for prime, good_count in local_good.items():
                        parameter = tuple(entry % prime for entry in residue)
                        if not any(parameter):
                            allowed = False
                            break
                        word = evaluate(rows, parameter, prime)
                        if sum(entry != 0 for entry in word) <= deletion_size:
                            allowed = False
                            break
                        require(good_count >= 0, ("negative local factor", rows, prime))
                    direct += allowed
                predicted = 1
                for good_count in local_good.values():
                    predicted *= good_count
                require(direct == predicted,
                        ("CRT vector product", rows, deletion_size, modulus, direct, predicted))
                crt_residue_checks += modulus ** rank

    return {
        "templates": len(templates),
        "local_k_checks": local_k_checks,
        "projective_checks": projective_checks,
        "cocircuit_checks": cocircuit_checks,
        "crt_prime_checks": crt_prime_checks,
        "crt_residue_checks": crt_residue_checks,
    }


def audit_exact_prime_family():
    rows_checked = 0
    primitive_parameters = 0
    summaries = []
    for prime in (2, 3, 5, 7):
        rows = ((1, 0), (1, prime), (1, 2 * prime), (0, 1))
        indices = []
        for deleted in range(4):
            retained = [row for index, row in enumerate(rows) if index != deleted]
            indices.append(lattice_index(retained, 2))
        require(indices == [1, 1, 1, prime], ("prime-family indices", prime, indices))
        words, distribution = code_data(rows, prime)
        expected = Counter({0: 1, 1: prime - 1, 3: prime - 1, 4: (prime - 1) ** 2})
        require(distribution == expected, ("prime-family enumerator", prime, distribution))
        bad = distribution[1]
        good = prime ** 2 - 1 - bad
        require(good == prime * (prime - 1), ("prime-family local good", prime, good))
        bad_lines = {
            normalize_projective(parameter, prime)
            for parameter, word in words.items()
            if any(parameter) and sum(entry != 0 for entry in word) <= 1
        }
        require(len(bad_lines) == 1, ("one projective bad line", prime, bad_lines))
        for first in range(-20, 21):
            for second in range(-20, 21):
                if gcd(abs(first), abs(second)) != 1:
                    continue
                values = (
                    first,
                    first + prime * second,
                    first + 2 * prime * second,
                    second,
                )
                deletion_gcds = [
                    gcd_many(value for index, value in enumerate(values) if index != deleted)
                    for deleted in range(4)
                ]
                require(deletion_gcds[:3] == [1, 1, 1],
                        ("prime-family trivial deletions", prime, first, second, deletion_gcds))
                require(deletion_gcds[3] == gcd(abs(first), prime),
                        ("prime-family nontrivial deletion", prime, first, second, deletion_gcds))
                primitive_parameters += 1
        rows_checked += 1
        summaries.append((prime, tuple(sorted(distribution.items())), good, len(bad_lines)))
    return rows_checked, primitive_parameters, summaries


def is_directed_triangle(block, out_neighbors):
    internal_outdegrees = sorted(
        sum(other in out_neighbors[vertex] for other in block if other != vertex)
        for vertex in block
    )
    return internal_outdegrees == [1, 1, 1]


def audit_paley_e8_bridge():
    vertices = tuple(range(7))
    residues = frozenset((1, 2, 4))
    nonresidues = frozenset((3, 5, 6))
    out_neighbors = {
        vertex: frozenset((vertex + step) % 7 for step in residues)
        for vertex in vertices
    }
    in_neighbors = {
        vertex: frozenset((vertex + step) % 7 for step in nonresidues)
        for vertex in vertices
    }
    score_histogram = Counter(len(out_neighbors[vertex]) for vertex in vertices)
    triangles = {
        frozenset(block)
        for block in combinations(vertices, 3)
        if is_directed_triangle(block, out_neighbors)
    }
    plus_blocks = {out_neighbors[vertex] for vertex in vertices}
    minus_blocks = {in_neighbors[vertex] for vertex in vertices}
    require(len(triangles) == 14, ("Paley triangle count", triangles))
    require(plus_blocks.isdisjoint(minus_blocks), "two Fano planes must be block-disjoint")
    require(triangles == plus_blocks | minus_blocks, "triangle orbit exhaustion")
    for family in (plus_blocks, minus_blocks):
        pair_counts = Counter()
        for block in family:
            for pair in combinations(sorted(block), 2):
                pair_counts[pair] += 1
        require(set(pair_counts.values()) == {1} and len(pair_counts) == comb(7, 2),
                ("Fano difference design", family, pair_counts))
    union_pair_counts = Counter()
    for block in triangles:
        for pair in combinations(sorted(block), 2):
            union_pair_counts[pair] += 1
    require(set(union_pair_counts.values()) == {2}, ("2-(7,3,2)", union_pair_counts))
    disjoint_pairs = []
    for left, right in combinations(sorted(triangles, key=lambda block: tuple(sorted(block))), 2):
        if left.isdisjoint(right):
            disjoint_pairs.append((left, right))
    require(len(disjoint_pairs) == 7, ("near-resolution size", disjoint_pairs))
    leftovers = {
        next(iter(set(vertices) - set(left) - set(right)))
        for left, right in disjoint_pairs
    }
    require(leftovers == set(vertices), ("near-resolution leftovers", leftovers))

    infinity = 7
    gauge_supports = {
        frozenset({infinity}) | in_neighbors[vertex] for vertex in vertices
    }
    complement_supports = {
        frozenset({vertex}) | out_neighbors[vertex] for vertex in vertices
    }
    require(len(gauge_supports) == len(complement_supports) == 7, "seven marked rows")
    require(gauge_supports.isdisjoint(complement_supports), "marked support families")
    masks = [sum(1 << point for point in support) for support in gauge_supports]
    code = {0}
    for mask in masks:
        code |= {word ^ mask for word in tuple(code)}
    weight_distribution = Counter(word.bit_count() for word in code)
    require(weight_distribution == Counter({0: 1, 4: 14, 8: 1}),
            ("e8 enumerator", weight_distribution))
    weight_four = {
        frozenset(index for index in range(8) if word & (1 << index))
        for word in code if word.bit_count() == 4
    }
    require(weight_four == gauge_supports | complement_supports,
            ("marked e8 support layer", weight_four))
    require(all((support - {infinity}) in minus_blocks for support in gauge_supports),
            "fixed infinity marker")
    require(all((support - {vertex}) == out_neighbors[vertex]
                for vertex, support in zip(vertices, [frozenset({v}) | out_neighbors[v] for v in vertices])),
            "varying row marker")

    wheel_good = {}
    for deletion_size in range(9):
        bad = sum(
            count for weight, count in weight_distribution.items()
            if 1 <= weight <= deletion_size
        )
        wheel_good[deletion_size] = 15 - bad
    require([wheel_good[k] for k in range(4)] == [15, 15, 15, 15], wheel_good)
    require(wheel_good[4] == 1 and wheel_good[8] == 0, wheel_good)

    hamiltonian_paths = 0
    for ordering in permutations(vertices):
        if all(ordering[index + 1] in out_neighbors[ordering[index]] for index in range(6)):
            hamiltonian_paths += 1
    require(score_histogram == Counter({3: 7}), ("Paley scores", score_histogram))
    require(hamiltonian_paths == 189, ("Paley Hamiltonian paths", hamiltonian_paths))
    return {
        "score_histogram": tuple(sorted(score_histogram.items())),
        "triangles": len(triangles),
        "near_classes": len(disjoint_pairs),
        "weight_distribution": tuple(sorted(weight_distribution.items())),
        "wheel_good": tuple(wheel_good[index] for index in range(9)),
        "hamiltonian_paths": hamiltonian_paths,
    }


def audit_code72_gate():
    blocks = 78 * comb(72, 5) // comb(16, 5)
    require(78 * comb(72, 5) % comb(16, 5) == 0, "design block integrality")
    require(blocks == 249849, ("weight-16 design count", blocks))
    total_directions = 2 ** 36 - 1
    good_at_15 = total_directions
    good_at_16 = total_directions - blocks
    require(good_at_16 == 68719226886, ("code72 wheel count", good_at_16))
    return blocks, total_directions, good_at_15, good_at_16


def audit_guardrails():
    rational_only = ((2, 0), (0, 1))
    bad_word = evaluate(rational_only, (1, 0), 2)
    require(bad_word == (0, 0), ("rational-span hostile control", bad_word))

    rows = ((1, 0), (1, 3), (1, 6), (0, 1))
    for prime in (2, 3, 5, 7):
        _, distribution = code_data(rows, prime)
        require(sum(distribution.get(weight, 0) for weight in range(1, 4)) > 0,
                ("rank-deficient k=3 should be bad at every prime", prime))
    return "Q-span is insufficient; rank-deficient complements have no finite wheel"


def main():
    local = audit_local_code_theorem()
    family_rows, primitive_parameters, family_summaries = audit_exact_prime_family()
    paley = audit_paley_e8_bridge()
    code72 = audit_code72_gate()
    guardrail = audit_guardrails()

    print("THM-2069 K-DELETION CODE/COGIRTH CRT-WHEEL AUDIT")
    print("integrally spanning templates:", local["templates"])
    print("local deletion/weight checks:", local["local_k_checks"])
    print("projective orbit checks:", local["projective_checks"])
    print("cocircuit/cogirth checks:", local["cocircuit_checks"])
    print("index-prime checks:", local["crt_prime_checks"])
    print("direct CRT residues checked:", local["crt_residue_checks"])
    print("exact prime-family rows / primitive parameters:", family_rows, primitive_parameters)
    print("prime-family (p, enumerator, good vectors, bad lines):", family_summaries)
    print("Paley P7 score histogram:", paley["score_histogram"])
    print("Paley directed triangles / near-classes:", paley["triangles"], paley["near_classes"])
    print("Paley Hamiltonian paths:", paley["hamiltonian_paths"])
    print("e8 weight enumerator:", paley["weight_distribution"])
    print("e8 k=0..8 good directions:", paley["wheel_good"])
    print("[72,36,16] A16 / all binary directions:", code72[0], code72[1])
    print("[72,36,16] good directions at k=15 / k=16:", code72[2], code72[3])
    print("guardrail:", guardrail)
    print("carrier: prime-labelled evaluation code with deletion radius")
    print("tournament vertices: Paley vertices; gauge: quadratic-residue orientation")
    print("PASS")


if __name__ == "__main__":
    main()
