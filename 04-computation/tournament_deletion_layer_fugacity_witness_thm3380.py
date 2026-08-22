#!/usr/bin/env python3
"""Generator-free exact witness companion for proposed THM-3380.

This standard-library-only program verifies three claims from frozen raw
upper-triangle tournament encodings.

1. At order eight, the deletion-layer semiring M strictly refines the pair
   (spectral marked-deletion Gram Gamma, Hamiltonian transform Dham).
2. The first Hamiltonian skew current and spectral skew current both miss the
   palindromic SCC-word collision ABBA versus BAAB for A=K1, B=C3.
3. The formal odd-cycle-packing fugacity F(x,y) first contains information
   beyond Dham=F(2,y) at order nine.

The algebraic laws used by the theorem are also checked on the frozen factors:
M and F multiply under ordered join.  All checks use explicit RuntimeError
guards, so optimized Python executes the same validation.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import math
import multiprocessing
from collections import Counter, defaultdict
from functools import lru_cache


N8_GRAM_FIBER = (
    "1110100111011111011111111111",
    "1110100111011111111101111110",
    "1101000111111111011110111111",
    "1100011110100111111111101110",
)
N8_GRAM_LEX = (
    "0000001000011001000001000100",
    "0000011000001001000001010100",
    "0000011000101000010000110101",
    "0000011000001001000110010011",
)
N8_PALINDROMES = (
    "1111111101111111111111101111",
    "1011111111111111111111111101",
)
N8_PALINDROME_LEX = (
    "0000000000001000000011001111",
    "0000001000000000010001011011",
)
N9_FUGACITY = (
    "111111001111001110101110101111111111",
    "111001011111110111001111011010111111",
)
N9_FUGACITY_LEX = (
    "000000110001101000110010000000010100",
    "000000110001101000101010100000001001",
)
C3_OF_C3_LEX = "000011110000111011000010001000111010"

EXPECTED_COMMON_S = (
    (1, 8, 56, 224, 704, 1472, 2176, 1920, 768),
    (0, 8, 56, 352, 1200, 3200, 5248, 5248, 2304),
)
EXPECTED_COMMON_DHAM = (52, 110, 80, 108, 0, 32, 0, 0, 1)
EXPECTED_GAMMA_SHA256 = (
    "70c32ef0e0f4f5b53ae4f9b0dfa93a258cba051019fafd495877e52a59a22458"
)
EXPECTED_MARKED_SPECTRAL_DECK_SHA256 = (
    "f707600d2804a502d7d91a2e4db6109f1aad26007f1679716194c451353b28c8"
)
EXPECTED_N9_DHAM = (394, 656, 664, 304, 262, 0, 50, 0, 0, 1)
EXPECTED_C3_OF_C3_DHAM = (710, 648, 972, 444, 324, 0, 60, 0, 0, 1)
EXPECTED_CENSUS = {
    2: {
        "classes": 1,
        "input_sha256": "4355a46b19d348dc2f57c046f8ef63d4538ebb936000f3c9ee954a27460dd865",
        "distinct": (1, 1, 1, 1),
        "gd_splits": (),
        "mxi_fibres": (),
    },
    3: {
        "classes": 2,
        "input_sha256": "b8bd9f18dbaa6e48bf08aea7368585a4e3e001c4914a4012b138c7aa8b1bb6a0",
        "distinct": (2, 2, 2, 2),
        "gd_splits": (),
        "mxi_fibres": (),
    },
    4: {
        "classes": 4,
        "input_sha256": "91dde84a6b6286a5e0b7b6295ccb1d16955dc4d74b8b93163d50d8a9cd7c7921",
        "distinct": (3, 3, 4, 4),
        "gd_splits": (),
        "mxi_fibres": (),
    },
    5: {
        "classes": 12,
        "input_sha256": "5d16cab0f65c58f402ef49c3442e3e8fedbb5498a8cb09f58a5d8d733784f6e4",
        "distinct": (9, 9, 12, 12),
        "gd_splits": (),
        "mxi_fibres": (),
    },
    6: {
        "classes": 56,
        "input_sha256": "814e0ed10e5d3aaae92d809b9eb915f1eb5542a65a33fe4d47463d22159a525b",
        "distinct": (32, 32, 56, 56),
        "gd_splits": (),
        "mxi_fibres": (),
    },
    7: {
        "classes": 456,
        "input_sha256": "164260b94960af0cc63faf3f178ceb95f4dd23bbca376ec23872c33b30d94261",
        "distinct": (247, 247, 456, 456),
        "gd_splits": (),
        "mxi_fibres": (),
    },
    8: {
        "classes": 6880,
        "input_sha256": "fc96c6997724e54ccea3bd166f4117d9e27925d85f568e31b0623527e5139dad",
        "distinct": (3328, 3325, 6879, 6859),
        "gd_splits": (
            (1049, 1133, 2620, 5234),
            (3667, 3673, 4339, 4563),
            (5468, 6021, 6691, 6718),
        ),
        "mxi_fibres": ((5314, 5791),),
    },
}


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def trim(polynomial):
    result = list(polynomial)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return tuple(result)


def polynomial_add(*polynomials):
    size = max((len(polynomial) for polynomial in polynomials), default=1)
    return trim(tuple(
        sum(polynomial[degree] if degree < len(polynomial) else 0
            for polynomial in polynomials)
        for degree in range(size)
    ))


def polynomial_scale(polynomial, scalar):
    return trim(tuple(scalar * coefficient for coefficient in polynomial))


def polynomial_multiply(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for left_degree, left_value in enumerate(left):
        for right_degree, right_value in enumerate(right):
            result[left_degree + right_degree] += left_value * right_value
    return trim(result)


def shifted_monomial_power(power):
    return tuple(
        math.comb(power, degree) * (-1) ** (power - degree)
        for degree in range(power + 1)
    )


def decode(bits: str):
    length = len(bits)
    order = (1 + math.isqrt(1 + 8 * length)) // 2
    require(order * (order - 1) // 2 == length, (bits, "bad length"))
    matrix = [[0] * order for _ in range(order)]
    edge = 0
    for left, right in itertools.combinations(range(order), 2):
        require(bits[edge] in "01", (bits, edge))
        value = int(bits[edge])
        edge += 1
        matrix[left][right] = value
        matrix[right][left] = 1 - value
    return tuple(tuple(row) for row in matrix)


def encode(matrix):
    return "".join(
        str(matrix[left][right])
        for left, right in itertools.combinations(range(len(matrix)), 2)
    )


def out_masks(matrix):
    return tuple(
        sum(1 << column for column, value in enumerate(row) if value)
        for row in matrix
    )


def delete_vertex(matrix, vertex):
    return tuple(
        tuple(value for column, value in enumerate(row) if column != vertex)
        for row_index, row in enumerate(matrix)
        if row_index != vertex
    )


def converse(matrix):
    return tuple(
        tuple(0 if row == column else matrix[column][row]
              for column in range(len(matrix)))
        for row in range(len(matrix))
    )


def ordered_join(left, right):
    left_order = len(left)
    right_order = len(right)
    order = left_order + right_order
    return tuple(
        tuple(
            left[row][column]
            if row < left_order and column < left_order
            else right[row - left_order][column - left_order]
            if row >= left_order and column >= left_order
            else 1 if row < left_order else 0
            for column in range(order)
        )
        for row in range(order)
    )


def ordered_join_many(*factors):
    result = ()
    for factor in factors:
        result = ordered_join(result, factor)
    return result


@lru_cache(maxsize=None)
def canonical_code(matrix):
    order = len(matrix)
    return min(
        tuple(
            matrix[permutation[left]][permutation[right]]
            for left, right in itertools.combinations(range(order), 2)
        )
        for permutation in itertools.permutations(range(order))
    )


def canonical_string(matrix):
    return "".join(map(str, canonical_code(matrix)))


def score_sequence(matrix):
    return tuple(sorted(sum(row) for row in matrix))


def ordered_scc_sizes(matrix):
    order = len(matrix)
    reach = [[row == column or bool(matrix[row][column])
              for column in range(order)] for row in range(order)]
    for pivot in range(order):
        for left in range(order):
            if reach[left][pivot]:
                for right in range(order):
                    reach[left][right] = reach[left][right] or reach[pivot][right]
    unseen = set(range(order))
    components = []
    while unseen:
        first = min(unseen)
        component = frozenset(
            vertex for vertex in unseen
            if reach[first][vertex] and reach[vertex][first]
        )
        unseen -= component
        components.append(component)
    beat_counts = {}
    for component in components:
        beat_counts[component] = sum(
            all(matrix[left][right] for left in component for right in other)
            for other in components if other != component
        )
    components.sort(key=lambda component: beat_counts[component], reverse=True)
    return tuple(len(component) for component in components)


def hamiltonian_counts(matrix):
    out = out_masks(matrix)
    order = len(out)
    full = (1 << order) - 1
    paths = [[0] * order for _mask in range(full + 1)]
    for vertex in range(order):
        paths[1 << vertex][vertex] = 1
    for mask in range(1, full + 1):
        for endpoint, count in enumerate(paths[mask]):
            if not count:
                continue
            available = out[endpoint] & (full ^ mask)
            while available:
                bit = available & -available
                available ^= bit
                successor = bit.bit_length() - 1
                paths[mask | bit][successor] += count
    counts = [1]
    counts.extend(sum(paths[mask]) for mask in range(1, full + 1))
    return tuple(counts)


def deletion_layers(matrix, hamiltonian=None):
    order = len(matrix)
    hamiltonian = hamiltonian or hamiltonian_counts(matrix)
    layers = [Counter() for _ in range(order + 1)]
    for remaining, value in enumerate(hamiltonian):
        layers[order - remaining.bit_count()][value] += 1
    return tuple(tuple(sorted(layer.items())) for layer in layers)


def m_multiply(left, right):
    result = [Counter() for _ in range(len(left) + len(right) - 1)]
    for left_degree, left_layer in enumerate(left):
        for right_degree, right_layer in enumerate(right):
            for left_h, left_count in left_layer:
                for right_h, right_count in right_layer:
                    result[left_degree + right_degree][left_h * right_h] += (
                        left_count * right_count
                    )
    return tuple(tuple(sorted(layer.items())) for layer in result)


def diagonal_transform(matrix, hamiltonian=None, permitted=None):
    order = len(matrix)
    full = (1 << order) - 1
    hamiltonian = hamiltonian or hamiltonian_counts(matrix)
    permitted = full if permitted is None else permitted
    result = (0,)
    deleted = permitted
    while True:
        remaining = full ^ deleted
        result = polynomial_add(
            result,
            polynomial_scale(
                shifted_monomial_power(deleted.bit_count()),
                hamiltonian[remaining],
            ),
        )
        if deleted == 0:
            break
        deleted = (deleted - 1) & permitted
    return result


def marked_diagonal_deletions(matrix, hamiltonian=None):
    order = len(matrix)
    full = (1 << order) - 1
    hamiltonian = hamiltonian or hamiltonian_counts(matrix)
    return tuple(
        diagonal_transform(matrix, hamiltonian, full ^ (1 << vertex))
        for vertex in range(order)
    )


def sparse_outer(left, right, scalar=1):
    return {
        (left_degree, right_degree): scalar * left_value * right_value
        for left_degree, left_value in enumerate(left)
        for right_degree, right_value in enumerate(right)
        if scalar * left_value * right_value
    }


def sparse_add(*terms):
    result = defaultdict(int)
    for term in terms:
        for monomial, coefficient in term.items():
            result[monomial] += coefficient
    return tuple(sorted(
        (monomial, coefficient)
        for monomial, coefficient in result.items() if coefficient
    ))


def xi_current(matrix, marked):
    terms = []
    for left in range(len(matrix)):
        for right in range(len(matrix)):
            if left == right:
                continue
            sign = matrix[left][right] - matrix[right][left]
            terms.append(sparse_outer(marked[left], marked[right], sign))
    return sparse_add(*terms)


def identity(order):
    return tuple(
        tuple(int(row == column) for column in range(order))
        for row in range(order)
    )


def matrix_product(left, right):
    rows = len(left)
    middle = len(right)
    columns = len(right[0]) if middle else 0
    return tuple(
        tuple(sum(left[row][index] * right[index][column]
                  for index in range(middle))
              for column in range(columns))
        for row in range(rows)
    )


def matrix_add_scalar_identity(matrix, scalar):
    return tuple(
        tuple(value + (scalar if row == column else 0)
              for column, value in enumerate(matrix[row]))
        for row in range(len(matrix))
    )


def matrix_scale(matrix, scalar):
    return tuple(tuple(scalar * value for value in row) for row in matrix)


def char_polynomial_i_minus(matrix):
    order = len(matrix)
    if not order:
        return (1,)
    power = identity(order)
    traces = []
    for _degree in range(1, order + 1):
        power = matrix_product(power, matrix)
        traces.append(sum(power[index][index] for index in range(order)))
    coefficients = [1]
    for degree in range(1, order + 1):
        numerator = -sum(
            traces[index - 1] * coefficients[degree - index]
            for index in range(1, degree + 1)
        )
        require(numerator % degree == 0, (matrix, degree, numerator))
        coefficients.append(numerator // degree)
    return tuple(coefficients)


def centered(matrix):
    return tuple(
        tuple(2 * matrix[row][column] - 1 for column in range(len(matrix)))
        for row in range(len(matrix))
    )


def adjugate_coefficients(matrix, denominator):
    order = len(matrix)
    if not order:
        return ()
    coefficients = [identity(order)]
    for degree in range(1, order):
        coefficients.append(matrix_add_scalar_identity(
            matrix_product(matrix, coefficients[-1]), denominator[degree]
        ))
    require(
        matrix_product(matrix, coefficients[-1])
        == matrix_scale(identity(order), -denominator[-1]),
        (matrix, denominator, "Cayley-Hamilton"),
    )
    return tuple(coefficients)


def response(matrix):
    order = len(matrix)
    if not order:
        return (1,), (0,)
    core = centered(matrix)
    denominator = char_polynomial_i_minus(core)
    adjugate = adjugate_coefficients(core, denominator)
    numerator = tuple(
        sum(coefficient[row][column]
            for row in range(order) for column in range(order))
        for coefficient in adjugate
    )
    return denominator, (0,) + numerator


def pad(polynomial, size):
    require(len(polynomial) <= size, (polynomial, size))
    return tuple(polynomial) + (0,) * (size - len(polynomial))


def dense_outer(left, right):
    return tuple(tuple(a * b for b in right) for a in left)


def dense_kernel_sum(kernels, size):
    return tuple(
        tuple(sum(kernel[row][column] for kernel in kernels)
              for column in range(size))
        for row in range(size)
    )


def spectral_packet(matrix):
    size = len(matrix) + 1
    parent = tuple(pad(polynomial, size) for polynomial in response(matrix))
    marked = tuple(
        tuple(pad(polynomial, size) for polynomial in response(delete_vertex(matrix, vertex)))
        for vertex in range(len(matrix))
    )
    zero = tuple(tuple(0 for _column in range(size)) for _row in range(size))
    gram = tuple(
        tuple(
            dense_kernel_sum(
                tuple(dense_outer(card[left], card[right]) for card in marked),
                size,
            ) if marked else zero
            for right in range(2)
        )
        for left in range(2)
    )
    omega_rows = []
    for response_left in range(2):
        omega_row = []
        for response_right in range(2):
            kernels = []
            for left in range(len(matrix)):
                for right in range(len(matrix)):
                    sign = matrix[left][right] - matrix[right][left]
                    if sign:
                        kernels.append(matrix_scale(
                            dense_outer(marked[left][response_left], marked[right][response_right]),
                            sign,
                        ))
            omega_row.append(dense_kernel_sum(tuple(kernels), size) if kernels else zero)
        omega_rows.append(tuple(omega_row))
    return parent, marked, gram, tuple(omega_rows)


def freeze(value):
    if isinstance(value, dict):
        return tuple(sorted((key, freeze(item)) for key, item in value.items()))
    if isinstance(value, (tuple, list)):
        return tuple(freeze(item) for item in value)
    return value


def sha(value):
    return hashlib.sha256(repr(freeze(value)).encode("ascii")).hexdigest()


def is_zero_dense(value):
    if isinstance(value, int):
        return value == 0
    return all(is_zero_dense(item) for item in value)


def directed_odd_cycles(matrix):
    out = out_masks(matrix)
    order = len(matrix)
    cycles = []
    for length in range(3, order + 1, 2):
        for support in itertools.combinations(range(order), length):
            anchor = support[0]
            for tail in itertools.permutations(support[1:]):
                cycle = (anchor,) + tail
                if all(out[cycle[index]] & (1 << cycle[(index + 1) % length])
                       for index in range(length)):
                    cycles.append(sum(1 << vertex for vertex in support))
    return tuple(cycles)


def packing_fugacity(matrix):
    order = len(matrix)
    full = (1 << order) - 1
    states = {(0, 0): 1}
    for cycle_mask in directed_odd_cycles(matrix):
        additions = defaultdict(int)
        for (used, number), count in tuple(states.items()):
            if not used & cycle_mask:
                additions[(used | cycle_mask, number + 1)] += count
        for state, count in additions.items():
            states[state] = states.get(state, 0) + count
    result = defaultdict(int)
    for (used, number), count in states.items():
        result[(number, order - used.bit_count())] += count
    return tuple(sorted(
        (x_degree, y_degree, coefficient)
        for (x_degree, y_degree), coefficient in result.items() if coefficient
    ))


def fugacity_multiply(left, right):
    result = defaultdict(int)
    for left_x, left_y, left_value in left:
        for right_x, right_y, right_value in right:
            result[(left_x + right_x, left_y + right_y)] += left_value * right_value
    return tuple(sorted(
        (x_degree, y_degree, coefficient)
        for (x_degree, y_degree), coefficient in result.items() if coefficient
    ))


def evaluate_fugacity_at_two(fugacity, order):
    result = [0] * (order + 1)
    for x_degree, y_degree, coefficient in fugacity:
        result[y_degree] += coefficient * 2 ** x_degree
    return trim(result)


def anchored_hamiltonian_cycle_count(matrix):
    out = out_masks(matrix)
    order = len(matrix)
    full = (1 << order) - 1
    paths = [[0] * order for _mask in range(full + 1)]
    paths[1][0] = 1
    for mask in range(full + 1):
        if not mask & 1:
            continue
        for endpoint, count in enumerate(paths[mask]):
            if not count:
                continue
            available = out[endpoint] & (full ^ mask)
            while available:
                bit = available & -available
                available ^= bit
                successor = bit.bit_length() - 1
                paths[mask | bit][successor] += count
    return sum(paths[full][endpoint]
               for endpoint in range(1, order) if out[endpoint] & 1)


def cyclic_triangle(matrix, vertices):
    return all(sum(matrix[vertex][other] for other in vertices) == 1
               for vertex in vertices)


def triangle_factor_count(matrix):
    order = len(matrix)
    require(order == 9, order)
    vertices = frozenset(range(order))
    count = 0
    first = 0
    for pair in itertools.combinations(range(1, order), 2):
        block1 = frozenset((first,) + pair)
        remaining = vertices - block1
        second_first = min(remaining)
        for second_pair in itertools.combinations(sorted(remaining - {second_first}), 2):
            block2 = frozenset((second_first,) + second_pair)
            block3 = remaining - block2
            if (cyclic_triangle(matrix, block1)
                    and cyclic_triangle(matrix, block2)
                    and cyclic_triangle(matrix, block3)):
                count += 1
    require(count <= 280, count)
    return count


def c3_matrix():
    return (
        (0, 1, 0),
        (0, 0, 1),
        (1, 0, 0),
    )


def c3_of_c3():
    inner = c3_matrix()
    matrix = [[0] * 9 for _row in range(9)]
    for part in range(3):
        for left in range(3):
            for right in range(3):
                matrix[3 * part + left][3 * part + right] = inner[left][right]
    for source, target in ((0, 1), (1, 2), (2, 0)):
        for left in range(3 * source, 3 * source + 3):
            for right in range(3 * target, 3 * target + 3):
                matrix[left][right] = 1
                matrix[right][left] = 0
    return tuple(tuple(row) for row in matrix)


def hamiltonian_packet(matrix):
    hamiltonian = hamiltonian_counts(matrix)
    diagonal = diagonal_transform(matrix, hamiltonian)
    marked = marked_diagonal_deletions(matrix, hamiltonian)
    return (
        deletion_layers(matrix, hamiltonian),
        diagonal,
        marked,
        xi_current(matrix, marked),
    )


def census_packet(bits):
    matrix = decode(bits)
    hamiltonian = hamiltonian_packet(matrix)
    spectral = spectral_packet(matrix)
    return (
        hamiltonian[0], hamiltonian[1], hamiltonian[3],
        spectral[0], spectral[2], spectral[3],
    )


def nontrivial_fibres(keys):
    buckets = defaultdict(list)
    for index, key in enumerate(keys):
        buckets[key].append(index)
    return tuple(sorted(
        (tuple(indices) for indices in buckets.values() if len(indices) > 1),
        key=lambda fibre: (len(fibre), fibre),
    ))


def fibres_split_by(left_keys, right_keys):
    parts = defaultdict(lambda: defaultdict(list))
    for index, (left, right) in enumerate(zip(left_keys, right_keys)):
        parts[left][right].append(index)
    return tuple(sorted(
        tuple(sorted(index for fibre in subdivision.values() for index in fibre))
        for subdivision in parts.values() if len(subdivision) > 1
    ))


def run_census(path, workers):
    with open(path, encoding="ascii") as source:
        lines = tuple(line.strip() for line in source if line.strip())
    require(lines, (path, "empty representative universe"))
    order = (1 + math.isqrt(1 + 8 * len(lines[0]))) // 2
    require(order in EXPECTED_CENSUS, (path, order))
    expected = EXPECTED_CENSUS[order]
    payload = ("\n".join(lines) + "\n").encode("ascii")
    input_sha256 = hashlib.sha256(payload).hexdigest()
    require(len(lines) == expected["classes"], (path, len(lines)))
    require(input_sha256 == expected["input_sha256"], (path, input_sha256))
    require(all(len(line) == order * (order - 1) // 2 for line in lines),
            (path, "line length"))
    if workers == 1:
        packets = tuple(map(census_packet, lines))
    else:
        with multiprocessing.Pool(workers) as pool:
            packets = tuple(pool.map(census_packet, lines, chunksize=16))
    m_keys = tuple(packet[0] for packet in packets)
    d_keys = tuple(packet[1] for packet in packets)
    xi_keys = tuple(packet[2] for packet in packets)
    s_keys = tuple(packet[3] for packet in packets)
    gamma_keys = tuple(packet[4] for packet in packets)
    omega_keys = tuple(packet[5] for packet in packets)
    gd_keys = tuple(zip(gamma_keys, d_keys))
    mxi_keys = tuple(zip(m_keys, xi_keys))
    somega_keys = tuple(zip(s_keys, omega_keys))
    distinct = tuple(map(lambda keys: len(set(keys)),
                         (m_keys, gd_keys, mxi_keys, somega_keys)))
    gd_splits = fibres_split_by(gd_keys, m_keys)
    m_splits = fibres_split_by(m_keys, gd_keys)
    mxi_fibres = nontrivial_fibres(mxi_keys)
    require(distinct == expected["distinct"], (order, distinct))
    require(gd_splits == expected["gd_splits"], (order, gd_splits))
    require(not m_splits, (order, m_splits))
    require(mxi_fibres == expected["mxi_fibres"], (order, mxi_fibres))
    print(
        f"n{order}_census=PASS classes={len(lines)} input_sha256={input_sha256} "
        f"distinct_M_GammaD_Mxi_sOmega={distinct} "
        f"GammaD_split_fibres={gd_splits} Mxi_fibres={mxi_fibres}"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--census", action="append", default=[],
                        help="pinned default gentourng order-7 or order-8 file")
    parser.add_argument("--workers", type=int, default=1)
    arguments = parser.parse_args()
    require(arguments.workers >= 1, arguments.workers)
    gram_matrices = tuple(map(decode, N8_GRAM_FIBER))
    require(
        tuple(map(canonical_string, gram_matrices)) == N8_GRAM_LEX,
        "order-eight Gram-fibre canonical codes",
    )
    gram_hamiltonian = tuple(map(hamiltonian_packet, gram_matrices))
    gram_spectral = tuple(map(spectral_packet, gram_matrices))
    require(all(packet[1] == EXPECTED_COMMON_DHAM for packet in gram_hamiltonian),
            "common Dham")
    require(all(packet[0] == EXPECTED_COMMON_S for packet in gram_spectral),
            "common s")
    require(len({freeze(packet[2]) for packet in gram_spectral}) == 1,
            "common Gamma")
    require(sha(gram_spectral[0][2]) == EXPECTED_GAMMA_SHA256,
            "Gamma digest")
    marked_spectral_multisets = tuple(
        tuple(sorted(packet[1])) for packet in gram_spectral
    )
    require(len(set(marked_spectral_multisets)) == 1,
            "common full marked spectral deck")
    require(sha(marked_spectral_multisets[0]) == EXPECTED_MARKED_SPECTRAL_DECK_SHA256,
            "marked spectral deck digest")
    require(gram_hamiltonian[0][0] == gram_hamiltonian[2][0], "first M pair")
    require(gram_hamiltonian[1][0] == gram_hamiltonian[3][0], "second M pair")
    require(gram_hamiltonian[0][0] != gram_hamiltonian[1][0], "M split")
    require(
        canonical_code(converse(gram_matrices[0])) == canonical_code(gram_matrices[2])
        and canonical_code(converse(gram_matrices[1])) == canonical_code(gram_matrices[3]),
        "converse pairing",
    )
    differing_layers = tuple(
        degree for degree, (left, right) in enumerate(zip(
            gram_hamiltonian[0][0], gram_hamiltonian[1][0]
        )) if left != right
    )
    require(differing_layers == (1, 2), differing_layers)
    require(
        sum(value * multiplicity for value, multiplicity in gram_hamiltonian[0][0][1])
        == sum(value * multiplicity for value, multiplicity in gram_hamiltonian[1][0][1])
        == 762,
        "first deletion augmentation",
    )
    require(
        sum(value * multiplicity for value, multiplicity in gram_hamiltonian[0][0][2])
        == sum(value * multiplicity for value, multiplicity in gram_hamiltonian[1][0][2])
        == 752,
        "second deletion augmentation",
    )

    palindrome_matrices = tuple(map(decode, N8_PALINDROMES))
    require(
        tuple(map(canonical_string, palindrome_matrices)) == N8_PALINDROME_LEX,
        "palindrome canonical codes",
    )
    singleton = ((0,),)
    triangle = c3_matrix()
    constructed = (
        ordered_join_many(singleton, triangle, triangle, singleton),
        ordered_join_many(triangle, singleton, singleton, triangle),
    )
    require(
        tuple(map(canonical_code, palindrome_matrices))
        == tuple(map(canonical_code, constructed)),
        "ABBA/BAAB identification",
    )
    require(canonical_code(palindrome_matrices[0]) != canonical_code(palindrome_matrices[1]),
            "palindromes nonisomorphic")
    require(tuple(map(ordered_scc_sizes, palindrome_matrices))
            == ((1, 3, 3, 1), (3, 1, 1, 3)), "SCC words")
    require(tuple(map(score_sequence, palindrome_matrices)) == (
        (0, 2, 2, 2, 5, 5, 5, 7),
        (1, 1, 1, 3, 4, 6, 6, 6),
    ), "score controls")
    require(all(canonical_code(matrix) == canonical_code(converse(matrix))
                for matrix in palindrome_matrices), "self-converse palindromes")
    palindrome_hamiltonian = tuple(map(hamiltonian_packet, palindrome_matrices))
    palindrome_spectral = tuple(map(spectral_packet, palindrome_matrices))
    require(palindrome_hamiltonian[0][0] == palindrome_hamiltonian[1][0], "common M")
    require(palindrome_hamiltonian[0][1] == palindrome_hamiltonian[1][1], "common Dham")
    require(palindrome_spectral[0][0] == palindrome_spectral[1][0], "common s")
    require(palindrome_spectral[0][2] == palindrome_spectral[1][2], "common Gamma")
    require(palindrome_hamiltonian[0][3] == palindrome_hamiltonian[1][3] == (),
            "xi vanishes")
    require(all(is_zero_dense(packet[3]) for packet in palindrome_spectral),
            "Omega vanishes")
    factor_m = (deletion_layers(singleton), deletion_layers(triangle))
    expected_m = m_multiply(m_multiply(factor_m[0], factor_m[1]),
                            m_multiply(factor_m[1], factor_m[0]))
    require(palindrome_hamiltonian[0][0] == palindrome_hamiltonian[1][0] == expected_m,
            "M ordered-join product")
    factor_f = (packing_fugacity(singleton), packing_fugacity(triangle))
    expected_f = fugacity_multiply(fugacity_multiply(factor_f[0], factor_f[1]),
                                   fugacity_multiply(factor_f[1], factor_f[0]))
    require(packing_fugacity(palindrome_matrices[0]) == expected_f
            == packing_fugacity(palindrome_matrices[1]), "F ordered-join product")

    n9_matrices = tuple(map(decode, N9_FUGACITY))
    require(tuple(map(canonical_string, n9_matrices)) == N9_FUGACITY_LEX,
            "order-nine canonical codes")
    n9_diagonal = tuple(diagonal_transform(matrix) for matrix in n9_matrices)
    require(n9_diagonal == (EXPECTED_N9_DHAM, EXPECTED_N9_DHAM), "common n9 Dham")
    n9_fugacity = tuple(map(packing_fugacity, n9_matrices))
    require(tuple(evaluate_fugacity_at_two(value, 9) for value in n9_fugacity)
            == n9_diagonal, "n9 fugacity specialization")
    n9_full = tuple(
        (next(coefficient for x_degree, y_degree, coefficient in fugacity
              if (x_degree, y_degree) == (1, 0)),
         next(coefficient for x_degree, y_degree, coefficient in fugacity
              if (x_degree, y_degree) == (3, 0)))
        for fugacity in n9_fugacity
    )
    require(n9_full == ((161, 9), (157, 10)), n9_full)
    require(tuple(anchored_hamiltonian_cycle_count(matrix) for matrix in n9_matrices)
            == (161, 157), "independent Hamiltonian cycles")
    require(tuple(triangle_factor_count(matrix) for matrix in n9_matrices)
            == (9, 10), "independent triangle factors")
    require(n9_fugacity[0] != n9_fugacity[1], "fugacity split")

    blowup = c3_of_c3()
    require(canonical_string(blowup) == C3_OF_C3_LEX, "C3[C3] canonical code")
    blowup_diagonal = diagonal_transform(blowup)
    blowup_fugacity = packing_fugacity(blowup)
    require(blowup_diagonal == EXPECTED_C3_OF_C3_DHAM, "C3[C3] Dham")
    require(evaluate_fugacity_at_two(blowup_fugacity, 9) == blowup_diagonal,
            "C3[C3] specialization")
    blowup_full = (
        next(coefficient for x_degree, y_degree, coefficient in blowup_fugacity
             if (x_degree, y_degree) == (1, 0)),
        next(coefficient for x_degree, y_degree, coefficient in blowup_fugacity
             if (x_degree, y_degree) == (3, 0)),
    )
    require(blowup_full == (207, 37), blowup_full)
    require(anchored_hamiltonian_cycle_count(blowup) == 207,
            "C3[C3] Hamiltonian cycles")
    require(triangle_factor_count(blowup) == 37 == 1 + math.factorial(3) ** 2,
            "C3[C3] triangle factors")

    semantic = (
        tuple(map(canonical_code, gram_matrices)),
        tuple(packet[0] for packet in gram_hamiltonian),
        gram_spectral[0][0], gram_spectral[0][1], gram_spectral[0][2],
        tuple(map(canonical_code, palindrome_matrices)),
        tuple(packet[0] for packet in palindrome_hamiltonian),
        palindrome_spectral[0][0], palindrome_spectral[0][2],
        tuple(map(canonical_code, n9_matrices)), n9_diagonal, n9_fugacity,
        canonical_code(blowup), blowup_diagonal, blowup_fugacity,
    )
    print("THM-3380 GENERATOR-FREE WITNESS AUDIT")
    print("n8_symmetric_fibre=PASS four_classes two_M_converse_pairs")
    print(f"common_s={EXPECTED_COMMON_S}")
    print(f"common_Dham={EXPECTED_COMMON_DHAM}")
    print(f"M_differing_deletion_sizes={differing_layers} sums=(762,752)")
    print("n8_palindrome=PASS ABBA_vs_BAAB common_M_s_Gamma xi=Omega=0")
    print(f"palindrome_SCC_words={tuple(map(ordered_scc_sizes, palindrome_matrices))}")
    print(f"n9_fugacity=PASS common_Dham={EXPECTED_N9_DHAM} full={n9_full}")
    print(f"C3_of_C3=PASS full={blowup_full} Dham={blowup_diagonal}")
    print(f"semantic_sha256={sha(semantic)}")
    for census_path in arguments.census:
        run_census(census_path, arguments.workers)


if __name__ == "__main__":
    main()
