#!/usr/bin/env python3
"""Independent finite boundary for THM-3428's rough maximal-order sector.

Universe: every odd 3 <= Q < 512 with least prime factor greater than seven;
every half-twist owner has quotient order exactly Q.  Coefficients are reduced
only by the lawful common sign symmetry.  Distinct coefficients whose masks
coincide accidentally remain distinct owner choices.  The search asks for
at most seven owners by checking every scalar-feasible exact size 1..7.
"""

from fractions import Fraction
from hashlib import sha256
from math import gcd, isqrt


EXPECTED_SEMANTIC_DIGEST = "b61db8412c1a45419cb4be595f75c802461472f24525b017b3d50d01af9f35ef"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def least_prime_factor(value):
    require(value >= 2, ("least-prime-factor domain", value))
    if value % 2 == 0:
        return 2
    for divisor in range(3, isqrt(value) + 1, 2):
        if value % divisor == 0:
            return divisor
    return value


def danger_mask(q, residue):
    modulus = 2 * q
    result = 0
    for sheet in range(q):
        word = residue * (2 * sheet + 1) % modulus
        if 14 * min(word, modulus - word) < modulus:
            result |= 1 << sheet
    return result


def field_mask(base, coefficient, q):
    inverse = pow(coefficient, -1, q)
    result = 0
    for value in base:
        result |= 1 << ((inverse * value) % q)
    return result


def odd_lift(coefficient, q):
    return coefficient if coefficient % 2 else coefficient + q


def block_bank(q):
    k = (q - 1) // 14
    even_base = {0, *range(1, k + 1), *range(q - k, q)}
    odd_radius = max(
        (value for value in range(1, q, 2) if 7 * value < q),
        default=-1,
    )
    odd_base = (
        {value % q for value in range(-odd_radius, odd_radius + 1, 2)}
        if odd_radius >= 1
        else set()
    )
    representatives = tuple(
        coefficient
        for coefficient in range(1, (q + 1) // 2)
        if gcd(coefficient, q) == 1
    )
    even = tuple(
        (coefficient, field_mask(even_base, coefficient, q))
        for coefficient in representatives
    )
    odd = tuple(
        (coefficient, field_mask(odd_base, coefficient, q))
        for coefficient in representatives
    )

    # Audit the coefficient-coordinate formulas against literal sheet masks.
    # The odd-sheet coordinate 2*ell+1 is a bijection modulo odd q.
    coordinate_to_sheet = {
        (2 * sheet + 1) % q: sheet
        for sheet in range(q)
    }
    require(len(coordinate_to_sheet) == q, (q, "odd coordinate not bijective"))
    for parity, bank in ((0, even), (1, odd)):
        for coefficient, coordinate_mask in bank:
            residue = 2 * coefficient if parity == 0 else odd_lift(coefficient, q)
            literal = danger_mask(q, residue)
            transported = 0
            for coordinate in range(q):
                if coordinate_mask >> coordinate & 1:
                    transported |= 1 << coordinate_to_sheet[coordinate]
            require(literal == transported,
                    (q, parity, coefficient, residue, literal ^ transported))
            require(gcd(q, residue) == 1, (q, residue, "not full order"))
    return even, odd, len(even_base), len(odd_base), odd_radius


def find_cover(q, even, odd, even_count, total_count, even_size, odd_size):
    require(1 <= even_count <= total_count <= 7, (q, even_count, total_count))
    odd_count = total_count - even_count
    full = (1 << q) - 1
    overlap_budget = even_count * even_size + odd_count * odd_size - q
    counters = [0, 0]

    def choose_odd(start, left, union, overlap, chosen_even, chosen_odd):
        counters[0] += 1
        if left == 0:
            return (chosen_even, chosen_odd) if union == full else None
        if union.bit_count() + left * odd_size < q:
            return None
        for index in range(start, len(odd) - left + 1):
            counters[1] += 1
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
        counters[0] += 1
        if left == 0:
            return choose_odd(0, odd_count, union, overlap, chosen, ())
        if union.bit_count() + left * even_size + odd_count * odd_size < q:
            return None
        for index in range(start, len(even) - left + 1):
            counters[1] += 1
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

    # Every cover contains an even block because only even residues contain
    # the reflection-fixed sheet.  Multiplicative normalization fixes its
    # coefficient to one.  Coefficient representatives, not masks, are chosen;
    # accidental duplicate masks therefore remain lawful distinct owners.
    answer = choose_even(1, even_count - 1, even[0][1], 0, (1,))
    return answer, counters[0], counters[1]


def literal_witness(q, even_coefficients, odd_coefficients):
    residues = tuple(
        sorted(
            tuple(2 * coefficient for coefficient in even_coefficients)
            + tuple(odd_lift(coefficient, q) for coefficient in odd_coefficients)
        )
    )
    require(len(residues) == len(set(residues)), (q, residues, "owner collision"))
    union = 0
    for residue in residues:
        require(gcd(q, residue) == 1, (q, residue, "not full quotient order"))
        union |= danger_mask(q, residue)
    require(union == (1 << q) - 1, (q, residues, "literal union failure"))
    return residues


def census():
    moduli = tuple(
        q
        for q in range(3, 512, 2)
        if least_prime_factor(q) > 7
    )
    rows = []
    positives = []
    feasible_profiles = 0
    nodes = 0
    branches = 0
    coefficient_cells = 0
    accidental_mask_duplicates = 0

    for q in moduli:
        even, odd, even_size, odd_size, odd_radius = block_bank(q)
        coefficient_cells += len(even) + len(odd)
        accidental_mask_duplicates += (
            len(even) - len({mask for _, mask in even})
            + len(odd) - len({mask for _, mask in odd})
        )
        q_positives = []
        for total_count in range(1, 8):
            for even_count in range(1, total_count + 1):
                odd_count = total_count - even_count
                overlap = even_count * even_size + odd_count * odd_size - q
                if overlap < even_count - 1:
                    continue
                if (overlap - (even_count - 1)) % 2:
                    continue
                feasible_profiles += 1
                answer, profile_nodes, profile_branches = find_cover(
                    q, even, odd, even_count, total_count, even_size, odd_size
                )
                nodes += profile_nodes
                branches += profile_branches
                if answer is not None:
                    residues = literal_witness(q, answer[0], answer[1])
                    q_positives.append((total_count, even_count, residues))
        minimum = min((row[0] for row in q_positives), default=None)
        if minimum is not None:
            positives.append((q, minimum, tuple(row for row in q_positives if row[0] == minimum)))
        rows.append(
            (
                q,
                least_prime_factor(q),
                len(even),
                len({mask for _, mask in even}),
                len(odd),
                len({mask for _, mask in odd}),
                even_size,
                odd_size,
                odd_radius,
                minimum,
            )
        )

    support = tuple(q for q, _, _ in positives)
    require(support == (11, 13, 23, 29), support)
    require(tuple((q, rank) for q, rank, _ in positives) ==
            ((11, 6), (13, 7), (23, 6), (29, 7)), positives)
    require(all(least_prime_factor(q) == q for q in support), support)
    composite_moduli = tuple(q for q in moduli if least_prime_factor(q) < q)
    require(not any(q in support for q in composite_moduli), (support, composite_moduli))
    hostiles = tuple(
        row
        for row in rows
        if row[0] in (31, 121, 143, 169, 323, 493, 509)
    )
    require(tuple(row[0] for row in hostiles) == (31, 121, 143, 169, 323, 493, 509), hostiles)
    return (
        moduli,
        composite_moduli,
        tuple(rows),
        tuple(positives),
        hostiles,
        feasible_profiles,
        nodes,
        branches,
        coefficient_cells,
        accidental_mask_duplicates,
    )


def main():
    result = census()
    semantic_surface = (
        result[0],
        result[1],
        result[3],
        result[4],
        result[5:],
    )
    semantic_digest = sha256(repr(semantic_surface).encode("ascii")).hexdigest()
    if EXPECTED_SEMANTIC_DIGEST is not None:
        require(semantic_digest == EXPECTED_SEMANTIC_DIGEST,
                (semantic_digest, EXPECTED_SEMANTIC_DIGEST))

    print("THM-3428 rough maximal-order finite boundary")
    print("status=FINITE_EXACT_INDEPENDENT_BOUNDARY;all_odd_Q_lt512_spf_gt7;all_full_order_coeff_reps;at_most7;duplicate_masks_retained_by_coefficient;no_mixed_orders;no_LRC14_decrement")
    print(f"universe=(moduli,composites,first,last)={(len(result[0]),len(result[1]),result[0][0],result[0][-1])}")
    print(f"search_counts=(scalar_feasible_profiles,nodes,branches,coefficient_cells,accidental_mask_duplicates)={result[5:]}")
    print(f"positive_support_and_minimal_witnesses={result[3]}")
    print(f"selected_hostiles=(Q,spf,even_coeffs,even_masks,odd_coeffs,odd_masks,even_size,odd_size,L,minrank)={result[4]}")
    print(f"semantic_sha256={semantic_digest}")
    print("scope=finite_boundary_for_THM3428_only;analytic_Q_ge512_theorem_separate;support_11,13,23,29_has_full_order_minranks_6,7,6,7;no_composite_positive;composite_mixed_order_rank7_open")


if __name__ == "__main__":
    main()
