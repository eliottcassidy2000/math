#!/usr/bin/env python3
"""Primewise extension of the consecutive Hasse-jet Smith scout.

This is standalone: it imports no producer and no prior scratch module.  It
uses exact Fraction arithmetic to perform Smith elimination over Z_(p), plus
an exhaustive determinantal-divisor path on the two smallest hostiles.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from math import comb, factorial, gcd
import sys


sys.stdout.reconfigure(newline="\n")

EXPECTED_SEMANTIC_SHA256 = "a9e30379316861df868043a67f3306f8a343e474aec052d22889be3eb122a989"
GATES = 0


def gate(condition, label):
    global GATES
    GATES += 1
    if not bool(condition):
        raise RuntimeError(label)


def primes_through(limit):
    result = []
    for candidate in range(2, limit + 1):
        if all(candidate % prime for prime in result if prime * prime <= candidate):
            result.append(candidate)
    return tuple(result)


def vp_integer(value, prime):
    gate(value != 0, "integer valuation nonzero")
    value = abs(value)
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def vp_fraction(value, prime):
    gate(value != 0, "fraction valuation nonzero")
    value = Fraction(value)
    return vp_integer(value.numerator, prime) - vp_integer(value.denominator, prime)


def hasse_matrix(nodes, k):
    size = len(nodes) * k
    return [
        [comb(degree, order) * node ** (degree - order)
         if degree >= order else 0
         for degree in range(size)]
        for node in nodes
        for order in range(k)
    ]


def p_smith_exponents(nodes, k, prime):
    """Exact Smith exponents over the DVR Z_(p)."""
    raw = hasse_matrix(nodes, k)
    work = [[Fraction(entry) for entry in row] for row in raw]
    size = len(work)
    exponents = []

    for pivot_index in range(size):
        row, column = min(
            ((row, column)
             for row in range(pivot_index, size)
             for column in range(pivot_index, size)
             if work[row][column]),
            key=lambda pair: vp_fraction(work[pair[0]][pair[1]], prime),
        )
        work[pivot_index], work[row] = work[row], work[pivot_index]
        for current_row in work:
            current_row[pivot_index], current_row[column] = (
                current_row[column], current_row[pivot_index]
            )

        pivot = work[pivot_index][pivot_index]
        pivot_exponent = vp_fraction(pivot, prime)
        exponents.append(pivot_exponent)

        for row in range(pivot_index + 1, size):
            if work[row][pivot_index] == 0:
                continue
            multiplier = work[row][pivot_index] / pivot
            gate(vp_fraction(multiplier, prime) >= 0,
                 "DVR-integral row multiplier")
            work[row] = [left - multiplier * right
                         for left, right in zip(work[row], work[pivot_index])]
            gate(work[row][pivot_index] == 0, "DVR row clear")

        for column in range(pivot_index + 1, size):
            if work[pivot_index][column] == 0:
                continue
            multiplier = work[pivot_index][column] / pivot
            gate(vp_fraction(multiplier, prime) >= 0,
                 "DVR-integral column multiplier")
            for row in range(pivot_index, size):
                work[row][column] -= multiplier * work[row][pivot_index]
            gate(work[pivot_index][column] == 0, "DVR column clear")

    gate(exponents == sorted(exponents), "DVR exponent chain")
    expected_total = k * k * sum(
        vp_integer(nodes[j] - nodes[i], prime)
        for i in range(len(nodes)) for j in range(i + 1, len(nodes))
    )
    gate(sum(exponents) == expected_total, "DVR determinant valuation")
    return tuple(exponents)


def bareiss_det(matrix):
    size = len(matrix)
    gate(size > 0 and all(len(row) == size for row in matrix), "square Bareiss")
    work = [list(row) for row in matrix]
    if size == 1:
        return work[0][0]
    sign = 1
    previous = 1
    for pivot_index in range(size - 1):
        pivot_row = next((row for row in range(pivot_index, size)
                          if work[row][pivot_index]), None)
        if pivot_row is None:
            return 0
        if pivot_row != pivot_index:
            work[pivot_index], work[pivot_row] = work[pivot_row], work[pivot_index]
            sign = -sign
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = (work[row][column] * pivot
                             - work[row][pivot_index] * work[pivot_index][column])
                gate(numerator % previous == 0, "Bareiss exact division")
                work[row][column] = numerator // previous
        for row in range(pivot_index + 1, size):
            work[row][pivot_index] = 0
        previous = pivot
    return sign * work[-1][-1]


def submatrix(matrix, rows, columns):
    return [[matrix[row][column] for column in columns] for row in rows]


def determinantal_divisors(matrix):
    size = len(matrix)
    answer = [1]
    for rank in range(1, size + 1):
        divisor = 0
        for rows in combinations(range(size), rank):
            for columns in combinations(range(size), rank):
                divisor = gcd(divisor, abs(bareiss_det(submatrix(matrix, rows, columns))))
                if divisor == 1:
                    break
            if divisor == 1:
                break
        answer.append(divisor)
    return tuple(answer)


def global_smith_from_divisors(divisors):
    diagonal = []
    for left, right in zip(divisors, divisors[1:]):
        gate(left and right % left == 0, "determinantal chain")
        diagonal.append(right // left)
    return tuple(diagonal)


def rle(exponents):
    result = []
    for exponent in exponents:
        if result and result[-1][0] == exponent:
            result[-1] = (exponent, result[-1][1] + 1)
        else:
            result.append((exponent, 1))
    return tuple(result)


def atlas_pairs():
    pairs = {(m, k) for m in range(2, 11) for k in range(1, 6)}
    pairs.update((m, 2) for m in tuple(range(11, 17)) + (20, 25, 30))
    pairs.update({
        (12, 3), (15, 3), (12, 4),
        (8, 6), (10, 6), (6, 7), (6, 8),
    })
    return tuple(sorted(pairs))


def main():
    gate(sys.argv[1:] in ([], ["--full"]), "supported command line")
    full_output = sys.argv[1:] == ["--full"]
    cache = {}

    def cached_exponents(nodes, k, prime):
        key = (tuple(nodes), k, prime)
        if key not in cache:
            cache[key] = p_smith_exponents(tuple(nodes), k, prime)
        return cache[key]

    profiles = []
    cluster_controls = []
    for m, k in atlas_pairs():
        nodes = tuple(range(m + 1))
        n = m + 1
        prime_profiles = []
        for prime in primes_through(m):
            exponents = cached_exponents(nodes, k, prime)
            zero_count = exponents.count(0)
            gate(zero_count == k * prime, ("zero layer", m, k, prime))

            first_positive = exponents[zero_count]
            expected_first = 1 + vp_integer(k, prime) if k % prime == 0 else 1
            gate(first_positive == expected_first,
                 ("first positive exponent", m, k, prime))

            active_residues = min(prime, n - prime)
            expected_ones = 0 if k % prime == 0 else active_residues
            gate(exponents.count(1) == expected_ones,
                 ("exponent-one multiplicity", m, k, prime))

            quotient, remainder = divmod(n, prime)
            local_profiles = []
            for residue in range(prime):
                cluster_size = quotient + (1 if residue < remainder else 0)
                local_nodes = tuple(prime * h for h in range(cluster_size))
                local_profiles.extend(cached_exponents(local_nodes, k, prime))
            local_profiles.sort()
            gate(tuple(local_profiles) == exponents,
                 ("p-local CRT cluster decomposition", m, k, prime))
            cluster_controls.append((m, k, prime, quotient, remainder))

            prime_profiles.append((prime, rle(exponents)))
        profiles.append((m, k, (m + 1) * k, tuple(prime_profiles)))

    # Independent all-minor crosschecks on the two smallest false Smith corners.
    exhaustive_controls = []
    for m, k, expected in (
        (3, 2, (1, 1, 1, 1, 4, 4, 12, 108)),
        (2, 3, (1, 1, 1, 1, 1, 1, 2, 16, 16)),
    ):
        matrix = hasse_matrix(tuple(range(m + 1)), k)
        divisors = determinantal_divisors(matrix)
        diagonal = global_smith_from_divisors(divisors)
        gate(diagonal == expected, ("all-minor hostile", m, k))
        for prime in primes_through(m):
            p_exponents = tuple(vp_integer(entry, prime) if entry != 1 else 0
                                for entry in diagonal)
            gate(p_exponents == cached_exponents(tuple(range(m + 1)), k, prime),
                 ("DVR/all-minor agreement", m, k, prime))
        exhaustive_controls.append((m, k, divisors, diagonal))

    # Proven k=2 pair band, obtained by p-local CRT plus the arbitrary two-node law.
    pair_band_controls = []
    for prime in primes_through(13):
        for n in range(prime + 1, 2 * prime + 1):
            exponents = cached_exponents(tuple(range(n)), 2, prime)
            collisions = n - prime
            if prime == 2:
                expected = tuple(sorted((0,) * (2 * prime) + (2,) * (2 * collisions)))
            else:
                expected = tuple(sorted((0,) * (2 * prime)
                                        + (1,) * collisions + (3,) * collisions))
            gate(exponents == expected, ("k2 pair band", prime, n))
            pair_band_controls.append((prime, n, rle(exponents)))

    # Two natural but false extensions, with their first confluent hostiles.
    even_pair = cached_exponents((0, 1, 2), 2, 2)
    gate(even_pair == (0, 0, 0, 0, 2, 2), "p divides k hostile partition")
    gate(even_pair.count(1) == 0 and min(2, 3 - 2) == 1,
         "unconditioned first-layer formula refuted")

    triple_hostile = cached_exponents(tuple(range(7)), 2, 3)
    gate(triple_hostile == (0, 0, 0, 0, 0, 0, 1, 1, 1, 3, 3, 3, 4, 4),
         "first triple-cluster hostile partition")
    gate(triple_hostile.count(1) == 3 and 7 - 3 == 4,
         "one-layer-per-extra-node formula refuted")

    # Exact bridge to the arbitrary two-node k=2 formula, including p^e gaps.
    two_node_controls = []
    for prime in primes_through(11):
        for gap_exponent in range(1, 5):
            distance = prime ** gap_exponent
            exponents = cached_exponents((0, distance), 2, prime)
            if prime == 2:
                expected = (0, 0, gap_exponent + 1, 3 * gap_exponent - 1)
            else:
                expected = (0, 0, gap_exponent, 3 * gap_exponent)
            gate(exponents == expected,
                 ("arbitrary two-node bridge", prime, gap_exponent))
            two_node_controls.append((prime, gap_exponent, exponents))

    # A wider finite signal for k=2 local clusters; deliberately not promoted.
    stable_local_controls = []
    for prime in primes_through(13):
        for cluster_size in range(2, prime):
            exponents = cached_exponents(
                tuple(prime * h for h in range(cluster_size)), 2, prime
            )
            expected = tuple(
                sorted((0, 0)
                       + tuple(range(1, cluster_size))
                       + tuple(range(cluster_size + 1, 2 * cluster_size)))
            )
            gate(exponents == expected,
                 ("finite local k2 stable pattern", prime, cluster_size))
            stable_local_controls.append((prime, cluster_size, rle(exponents)))

    semantic = {
        "profiles": [
            [m, k, size,
             [[prime, [[exponent, count] for exponent, count in compressed]]
              for prime, compressed in prime_profiles]]
            for m, k, size, prime_profiles in profiles
        ],
        "cluster_controls": [list(row) for row in cluster_controls],
        "exhaustive_controls": [
            [m, k, list(divisors), list(diagonal)]
            for m, k, divisors, diagonal in exhaustive_controls
        ],
        "pair_band": [
            [prime, n, [[exponent, count] for exponent, count in compressed]]
            for prime, n, compressed in pair_band_controls
        ],
        "hostiles": {
            "p_divides_k": list(even_pair),
            "triple_cluster": list(triple_hostile),
        },
        "two_node": [[prime, exponent, list(profile)]
                     for prime, exponent, profile in two_node_controls],
        "stable_local_finite": [
            [prime, size, [[exponent, count] for exponent, count in compressed]]
            for prime, size, compressed in stable_local_controls
        ],
    }
    semantic_digest = sha256(json.dumps(
        semantic, sort_keys=True, separators=(",", ":")
    ).encode("ascii")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_FROZEN":
        gate(semantic_digest == EXPECTED_SEMANTIC_SHA256, "semantic freeze")

    print("CONFLUENT_SAMPLER_PRIMEWISE_EXTENSION_20260824")
    print("status=PROVED_FIRST_POSITIVE_LAYER;FINITE_EXACT_WIDE_ATLAS;FULL_SNF_OPEN")
    print("source=mod_p_kernel_G_power_k_Q;target=duplicate_node_divided_Hasse_jets")
    print("map_first_term=Q(c)_to_k*h*Gprime(c)^k*Q(c)_in_jet_r_equals_k_minus_1")
    print("preserved=zero_layer_first_positive_exponent_and_if_p_not_divide_k_its_multiplicity")
    print("lost=higher_Bocksteins_and_positive_exponent_partition_beyond_first_layer")
    print("theorem=n_gt_p:first_kp_exponents_zero;next_exponent=1+v_p(k)")
    print("theorem=number_exponent_one=0_if_p_divides_k_else_min(p,n-p)")
    print("theorem=p_local_CRT_partition_is_union_of_residue_cluster_partitions")
    print("theorem=k2_pair_band:p_lt_n_le_2p;odd_p=(0^(2p),1^(n-p),3^(n-p));p2=(0^4,2^(2n-4))")
    print(f"atlas_pairs={len(profiles)};prime_partitions={sum(len(row[3]) for row in profiles)};max_rank={max(row[2] for row in profiles)}")
    displayed_pairs = {
        (2, 2), (2, 3), (3, 2), (3, 3), (6, 2), (6, 4), (6, 8),
        (10, 2), (10, 6), (15, 2), (15, 3), (20, 2), (25, 2), (30, 2),
    }
    for m, k, size, prime_profiles in profiles:
        if full_output or (m, k) in displayed_pairs:
            print(f"atlas(m={m},k={k},N={size})=" + repr(prime_profiles).replace(" ", ""))
    print(f"atlas_rows_printed={'all' if full_output else len(displayed_pairs)};full_command=append_--full")
    print("hostile_ignore_p_divides_k=(p=2,m=2,k=2);actual=(0,0,0,0,2,2);false_ones=1")
    print("hostile_one_per_extra_node=(p=3,m=6,k=2);actual=(0^6,1^3,3^3,4^2);false_ones=4")
    print(f"pair_band_controls={len(pair_band_controls)};two_node_controls={len(two_node_controls)}")
    print(f"finite_local_k2_stable_controls={len(stable_local_controls)};status=FINITE_EXACT_NOT_PROVED_GENERAL")
    print(f"semantic_sha256={semantic_digest}")
    print(f"gates={GATES}")
    print("RESULT=PASS;HIGHER_POSITIVE_PARTITION=OPEN")


if __name__ == "__main__":
    main()
