#!/usr/bin/env python3
"""Independent referee for the primewise confluent Hasse sampler extension.

No producer or repository computation is imported.  The main exact path is a
modular Z/p^E elimination, deliberately different from the producer's
Fraction/Z_(p) elimination.  A division-free subset-DP determinant and an
all-minor determinantal-divisor calculation provide a second path at the
minimal dyadic hostile.
"""

from hashlib import sha256
from itertools import combinations
import json
from math import comb, gcd
import sys


sys.stdout.reconfigure(newline="\n")

EXPECTED_SEMANTIC_SHA256 = "0f8dc0d0d56792c2526b4395230393963e77df72c9b102d2b6086fae3ffd0f26"
GATES = 0


def gate(condition, label):
    global GATES
    GATES += 1
    if not bool(condition):
        raise RuntimeError(label)


def primes_through(limit):
    sieve = [True] * (limit + 1)
    if limit >= 0:
        sieve[0] = False
    if limit >= 1:
        sieve[1] = False
    for value in range(2, int(limit ** 0.5) + 1):
        if sieve[value]:
            for multiple in range(value * value, limit + 1, value):
                sieve[multiple] = False
    return tuple(value for value, is_prime in enumerate(sieve) if is_prime)


def vp_nonzero(value, prime):
    gate(value != 0, "valuation input nonzero")
    value = abs(value)
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def vp_mod(value, prime, precision):
    """Valuation of a residue, with zero assigned the precision cap."""
    if value == 0:
        return precision
    exponent = 0
    while exponent < precision and value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def vp_factorial(length, prime):
    exponent = 0
    quotient = length
    while quotient:
        quotient //= prime
        exponent += quotient
    return exponent


def hasse_matrix(nodes, jet_order):
    size = len(nodes) * jet_order
    return tuple(
        tuple(
            comb(degree, derivative_order) * node ** (degree - derivative_order)
            if degree >= derivative_order else 0
            for degree in range(size)
        )
        for node in nodes
        for derivative_order in range(jet_order)
    )


def determinant_valuation(nodes, jet_order, prime):
    """Confluent Vandermonde valuation, derived from pairwise differences."""
    return jet_order * jet_order * sum(
        vp_nonzero(nodes[right] - nodes[left], prime)
        for left in range(len(nodes))
        for right in range(left + 1, len(nodes))
    )


def p_smith_exponents_modular(nodes, jet_order, prime):
    """Smith exponents via elimination in Z/(p^(v(det)+1)).

    A minimum-valuation pivot divides every remaining entry in the DVR.  Its
    unit part is inverted only modulo the precision required to clear a row or
    column.  Since the precision exceeds the total determinant valuation, no
    true invariant exponent is hidden by reduction modulo p^E.
    """
    raw = hasse_matrix(nodes, jet_order)
    size = len(raw)
    total = determinant_valuation(nodes, jet_order, prime)
    precision = total + 1
    modulus = prime ** precision
    work = [[entry % modulus for entry in row] for row in raw]
    exponents = []

    for pivot_index in range(size):
        candidates = (
            (vp_mod(work[row][column], prime, precision), row, column)
            for row in range(pivot_index, size)
            for column in range(pivot_index, size)
            if work[row][column] % modulus
        )
        try:
            exponent, pivot_row, pivot_column = min(candidates)
        except ValueError as exc:
            raise RuntimeError(("modular reducer lost full rank", nodes,
                                jet_order, prime, pivot_index)) from exc
        gate(exponent < precision, "visible modular pivot")

        work[pivot_index], work[pivot_row] = work[pivot_row], work[pivot_index]
        for row in work:
            row[pivot_index], row[pivot_column] = (
                row[pivot_column], row[pivot_index]
            )

        pivot = work[pivot_index][pivot_index]
        prime_power = prime ** exponent
        unit_modulus = prime ** (precision - exponent)
        gate(pivot % prime_power == 0, "pivot valuation divisibility")
        unit = (pivot // prime_power) % unit_modulus
        gate(gcd(unit, prime) == 1, "pivot unit part")
        unit_inverse = pow(unit, -1, unit_modulus)

        for row_index in range(pivot_index + 1, size):
            entry = work[row_index][pivot_index]
            if entry == 0:
                continue
            gate(entry % prime_power == 0, "row entry divisible by pivot power")
            multiplier = ((entry // prime_power) * unit_inverse) % unit_modulus
            work[row_index] = [
                (left - multiplier * right) % modulus
                for left, right in zip(work[row_index], work[pivot_index])
            ]
            gate(work[row_index][pivot_index] == 0, "modular row clear")

        for column_index in range(pivot_index + 1, size):
            entry = work[pivot_index][column_index]
            if entry == 0:
                continue
            gate(entry % prime_power == 0,
                 "column entry divisible by pivot power")
            multiplier = ((entry // prime_power) * unit_inverse) % unit_modulus
            for row_index in range(pivot_index, size):
                work[row_index][column_index] = (
                    work[row_index][column_index]
                    - multiplier * work[row_index][pivot_index]
                ) % modulus
            gate(work[pivot_index][column_index] == 0, "modular column clear")

        exponents.append(exponent)

    gate(exponents == sorted(exponents), "modular DVR exponent chain")
    gate(sum(exponents) == total, "modular DVR determinant mass")
    return tuple(exponents)


def determinant_subset_dp(matrix):
    """Division-free determinant, summing permutations by subset DP."""
    size = len(matrix)
    gate(size > 0 and all(len(row) == size for row in matrix),
         "subset determinant square")
    states = {0: 1}
    for row_index in range(size):
        next_states = {}
        for mask, coefficient in states.items():
            for column in range(size):
                if mask & (1 << column):
                    continue
                entry = matrix[row_index][column]
                if entry == 0:
                    continue
                inversions_added = (mask >> (column + 1)).bit_count()
                signed = -entry if inversions_added % 2 else entry
                new_mask = mask | (1 << column)
                next_states[new_mask] = (
                    next_states.get(new_mask, 0) + coefficient * signed
                )
        states = next_states
    return states.get((1 << size) - 1, 0)


def submatrix(matrix, rows, columns):
    return tuple(tuple(matrix[row][column] for column in columns) for row in rows)


def smith_from_all_minors(matrix):
    """Global Smith form from every determinantal divisor (small ranks only)."""
    size = len(matrix)
    divisors = [1]
    for rank in range(1, size + 1):
        divisor = 0
        for rows in combinations(range(size), rank):
            for columns in combinations(range(size), rank):
                minor = determinant_subset_dp(submatrix(matrix, rows, columns))
                divisor = gcd(divisor, abs(minor))
                if divisor == 1:
                    break
            if divisor == 1:
                break
        divisors.append(divisor)
    diagonal = []
    for left, right in zip(divisors, divisors[1:]):
        gate(left != 0 and right % left == 0, "determinantal divisor chain")
        diagonal.append(right // left)
    return tuple(divisors), tuple(diagonal)


def rle(values):
    compressed = []
    for value in values:
        if compressed and compressed[-1][0] == value:
            compressed[-1] = (value, compressed[-1][1] + 1)
        else:
            compressed.append((value, 1))
    return tuple(compressed)


def local_cluster_sizes(n, prime):
    quotient, remainder = divmod(n, prime)
    return tuple(quotient + (1 if residue < remainder else 0)
                 for residue in range(prime))


def expected_pair_band(prime, n):
    collisions = n - prime
    if prime == 2:
        return tuple(sorted((0,) * 4 + (2,) * (2 * collisions)))
    return tuple(sorted((0,) * (2 * prime)
                        + (1,) * collisions + (3,) * collisions))


def atlas_cases():
    specifications = (
        (2, range(1, 7), range(3, 10)),
        (3, range(1, 7), range(4, 12)),
        (5, range(1, 6), range(6, 14)),
        (7, range(1, 5), range(8, 16)),
        (11, range(1, 4), range(12, 16)),
    )
    cases = {
        (prime, n, jet_order)
        for prime, jet_orders, node_counts in specifications
        for jet_order in jet_orders
        for n in node_counts
    }
    cases.update({
        (2, 3, 8), (2, 5, 8), (2, 3, 16), (2, 5, 16),
        (3, 4, 9), (3, 7, 9),
    })
    return tuple(sorted(cases))


def main():
    gate(sys.argv[1:] == [], "no command-line arguments")
    profile_cache = {}

    def profile(nodes, jet_order, prime):
        key = (tuple(nodes), jet_order, prime)
        if key not in profile_cache:
            profile_cache[key] = p_smith_exponents_modular(
                tuple(nodes), jet_order, prime
            )
        return profile_cache[key]

    # Arithmetic audit of the divided-Hasse bound.  t=k-r is in [1,k],
    # s>=0, and L=s+t.  The numerator interval of C(k+s,L) contains k.
    inequality_controls = []
    equality_risks = []
    for prime in primes_through(19):
        for jet_order in range(1, 65):
            for s in range(0, 65):
                for t in range(1, jet_order + 1):
                    length = s + t
                    binomial = comb(jet_order + s, length)
                    left = length + vp_nonzero(binomial, prime)
                    right = 1 + vp_nonzero(jet_order, prime)
                    gate(left >= right,
                         ("divided-Hasse valuation", prime, jet_order, s, t))
                    gate(vp_factorial(length, prime) <= length - 1,
                         ("factorial valuation", prime, length))
                    if length > 1 and left == right:
                        equality_risks.append((prime, jet_order, s, t, length))
        inequality_controls.append(prime)
    gate((2, 2, 1, 1, 2) in equality_risks,
         "small-characteristic higher-L equality risk exposed")

    # The cancellation-free attainment repair: after local CRT and clearing
    # the representative-node identity block, the residual entry at duplicate
    # h=1, row r=k-1, column q=k is literally p*k.  Every other residual
    # entry has L=q-r>=1 and obeys the same lower bound.
    attainment_controls = []
    for prime in primes_through(19):
        for jet_order in range(1, 65):
            target = 1 + vp_nonzero(jet_order, prime)
            minimum = None
            minimizers = []
            for h in range(1, 7):
                for derivative_order in range(jet_order):
                    for degree in range(jet_order, 2 * jet_order + 3):
                        entry = (comb(degree, derivative_order)
                                 * (prime * h) ** (degree - derivative_order))
                        valuation = vp_nonzero(entry, prime)
                        if minimum is None or valuation < minimum:
                            minimum = valuation
                            minimizers = [(h, derivative_order, degree)]
                        elif valuation == minimum:
                            minimizers.append((h, derivative_order, degree))
            gate(minimum == target,
                 ("cancellation-free attainment", prime, jet_order))
            gate((1, jet_order - 1, jet_order) in minimizers,
                 ("literal pk pivot", prime, jet_order))
            attainment_controls.append((prime, jet_order, minimum))

    # Independent finite atlas.  Each global profile is also compared with a
    # separately assembled sorted union of local p-dilated cluster profiles.
    atlas = []
    maximum_rank = 0
    for prime, n, jet_order in atlas_cases():
        nodes = tuple(range(n))
        exponents = profile(nodes, jet_order, prime)
        maximum_rank = max(maximum_rank, len(exponents))
        gate(exponents[:jet_order * prime] == (0,) * (jet_order * prime),
             ("zero layer", prime, n, jet_order))
        first_positive = exponents[jet_order * prime]
        gate(first_positive == 1 + vp_nonzero(jet_order, prime),
             ("first positive", prime, n, jet_order))
        expected_ones = (0 if jet_order % prime == 0
                         else min(prime, n - prime))
        gate(exponents.count(1) == expected_ones,
             ("exponent-one count", prime, n, jet_order))

        local_union = []
        cluster_profiles = []
        for cluster_size in local_cluster_sizes(n, prime):
            local = profile(
                tuple(prime * h for h in range(cluster_size)),
                jet_order,
                prime,
            )
            local_union.extend(local)
            cluster_profiles.append(rle(local))
        local_union.sort()
        gate(tuple(local_union) == exponents,
             ("CRT cluster union", prime, n, jet_order))
        atlas.append((prime, n, jet_order, rle(exponents), tuple(cluster_profiles)))

    # Complete k=2 pair band, extending the producer's p<=13 atlas to p<=17.
    pair_band = []
    for prime in primes_through(17):
        for n in range(prime + 1, 2 * prime + 1):
            exponents = profile(tuple(range(n)), 2, prime)
            expected = expected_pair_band(prime, n)
            gate(exponents == expected, ("complete pair band", prime, n))
            pair_band.append((prime, n, rle(exponents)))

    # Minimal hostile (2,2,2): a genuinely separate all-minor path obtains
    # the complete global Smith diagonal, then agrees primewise.
    dyadic_matrix = hasse_matrix((0, 1, 2), 2)
    dyadic_divisors, dyadic_smith = smith_from_all_minors(dyadic_matrix)
    gate(dyadic_smith == (1, 1, 1, 1, 4, 4),
         "minimal dyadic hostile all-minor Smith")
    dyadic_profile = profile((0, 1, 2), 2, 2)
    gate(dyadic_profile == (0, 0, 0, 0, 2, 2),
         "minimal dyadic hostile p-profile")
    gate(tuple(vp_nonzero(value, 2) if value != 1 else 0
               for value in dyadic_smith) == dyadic_profile,
         "all-minor/modular dyadic agreement")

    # Minimal triple-cluster hostile under the explicit genuinely-confluent
    # constraints k>=2 and p not dividing k.
    triple_profile = profile(tuple(range(7)), 2, 3)
    gate(triple_profile == ((0,) * 6 + (1,) * 3 + (3,) * 3 + (4,) * 2),
         "triple-cluster hostile profile")
    earlier_hostiles = []
    candidate_hostiles = []
    for prime in primes_through(11):
        for jet_order in range(2, 8):
            if jet_order % prime == 0:
                continue
            for n in range(2 * prime + 1, 3 * prime + 2):
                rank = n * jet_order
                if rank > 14:
                    continue
                candidate_profile = profile(tuple(range(n)), jet_order, prime)
                if candidate_profile.count(1) != n - prime:
                    candidate_hostiles.append((rank, prime, n, jet_order,
                                               candidate_profile.count(1),
                                               n - prime))
                    if rank < 14:
                        earlier_hostiles.append(candidate_hostiles[-1])
    gate(not earlier_hostiles, "no smaller genuinely confluent triple hostile")
    gate((14, 3, 7, 2, 3, 4) in candidate_hostiles,
         "rank-14 triple hostile found")

    # Boundary n<=p: determinant is a p-unit, hence every exponent is zero.
    boundary_controls = []
    for prime in primes_through(13):
        for n in range(1, prime + 1):
            for jet_order in (1, 2, 3, 5):
                boundary_profile = profile(tuple(range(n)), jet_order, prime)
                gate(boundary_profile == (0,) * (n * jet_order),
                     ("n<=p boundary", prime, n, jet_order))
                boundary_controls.append((prime, n, jet_order))

    semantic = {
        "inequality_primes": inequality_controls,
        "higher_L_equality_risk_first": list(min(equality_risks)),
        "higher_L_equality_risk_count": len(equality_risks),
        "attainment": [list(row) for row in attainment_controls],
        "atlas": [
            [prime, n, jet_order,
             [[value, count] for value, count in profile],
             [[[value, count] for value, count in local]
              for local in cluster_profiles]]
            for prime, n, jet_order, profile, cluster_profiles in atlas
        ],
        "pair_band": [
            [prime, n, [[value, count] for value, count in profile]]
            for prime, n, profile in pair_band
        ],
        "dyadic": {
            "divisors": list(dyadic_divisors),
            "smith": list(dyadic_smith),
            "profile": list(dyadic_profile),
        },
        "triple_profile": list(triple_profile),
        "triple_candidates": [
            [rank, prime, n, jet_order, actual, predicted]
            for rank, prime, n, jet_order, actual, predicted
            in sorted(candidate_hostiles)[:12]
        ],
        "boundary_count": len(boundary_controls),
    }
    semantic_digest = sha256(json.dumps(
        semantic, sort_keys=True, separators=(",", ":")
    ).encode("ascii")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_FROZEN":
        gate(semantic_digest == EXPECTED_SEMANTIC_SHA256, "semantic freeze")

    print("CONFLUENT_SAMPLER_PRIMEWISE_EXTENSION_INDEPENDENT_REFEREE_20260824")
    print("status=PASS_WITH_PROOF_REPAIR;NO_IMPORT;EXACT_MODULAR_DVR")
    print("producer_algorithm_not_used=Fraction_Z_(p)_elimination")
    print("independent_algorithm=Z_mod_p_power_elimination_plus_subset_DP_all_minors")
    print("theorem=n_gt_p:alpha_1_through_alpha_kp_zero;alpha_kp_plus_1=1+v_p(k)")
    print("theorem=ones_zero_if_p_divides_k_else_min(p,n-p)")
    print("theorem=p_local_partition_sorted_union_of_residue_clusters")
    print("theorem=k2_pair_band_odd=(0^(2p),1^(n-p),3^(n-p));p2=(0^4,2^(2n-4))")
    print("proof_audit=L_factorial_bound_valid_but_summand_attainment_not_cancellation_safe")
    print("proof_repair=local_CRT_clear_identity_block;literal_residual_entry_(h=1,r=k-1,q=k)=p*k")
    print(f"higher_L_same_valuation_risks={len(equality_risks)};first={min(equality_risks)}")
    print(f"inequality_primes={len(inequality_controls)};k_max=64;s_max=64")
    print(f"attainment_controls={len(attainment_controls)};k_max=64")
    print(f"atlas_cases={len(atlas)};max_rank={maximum_rank};p_max=11;k_max=16;n_max=15")
    print(f"pair_band_cases={len(pair_band)};p_max=17")
    print("hostile_(p,m,k)=(2,2,2):profile=(0^4,2^2);global_smith=(1,1,1,1,4,4)")
    print("hostile_(p,m,k)=(3,6,2):profile=(0^6,1^3,3^3,4^2);false_ones=4")
    print(f"dyadic_determinantal_divisors={dyadic_divisors}")
    print(f"boundary_n_le_p_controls={len(boundary_controls)}")
    print(f"semantic_sha256={semantic_digest}")
    print(f"gates={GATES}")
    print("RESULT=PASS;WRITTEN_ATTAINMENT_STEP_NEEDS_REPAIR;THEOREM_SURVIVES")


if __name__ == "__main__":
    main()
