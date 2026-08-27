#!/usr/bin/env python3
"""Dependency-free finite controls for THM-4274's observer transport lemma."""

from fractions import Fraction
from itertools import product
import sys


def require(condition, message="audit failure"):
    if not condition:
        raise RuntimeError(message)


def trim(poly, p):
    out = [x % p for x in poly]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def inv_mod(a, p):
    return pow(a % p, -1, p)


def divmod_poly(a, b, p):
    a = trim(a, p)
    b = trim(b, p)
    if b == [0]:
        raise ZeroDivisionError
    q = [0] * max(1, len(a) - len(b) + 1)
    while a != [0] and len(a) >= len(b):
        k = len(a) - len(b)
        c = a[-1] * inv_mod(b[-1], p) % p
        q[k] = c
        for j, bj in enumerate(b):
            a[j + k] = (a[j + k] - c * bj) % p
        a = trim(a, p)
    return trim(q, p), a


def gcd_poly(a, b, p):
    a, b = trim(a, p), trim(b, p)
    while b != [0]:
        _, remainder = divmod_poly(a, b, p)
        a, b = b, remainder
    lead_inv = inv_mod(a[-1], p)
    return trim([(lead_inv * x) % p for x in a], p)


def apply_poly(poly, word, p):
    """Apply sum poly[k] S^k to a cyclic word, Sx[j]=x[j+1]."""
    length = len(word)
    return tuple(
        sum(poly[k] * word[(j + k) % length] for k in range(len(poly))) % p
        for j in range(length)
    )


def monic_polys(p, max_degree):
    answer = []
    for degree in range(1, max_degree + 1):
        for coeffs in product(range(p), repeat=degree):
            answer.append(tuple(coeffs) + (1,))
    return answer


def audit_field_gcd_recurrence():
    pair_cells = 0
    kernel_equalities = 0
    common_words = 0
    prefix_pairs = 0
    for p, length, max_degree in ((2, 5, 3), (3, 5, 2)):
        words = list(product(range(p), repeat=length))
        polynomials = monic_polys(p, max_degree)
        for poly_p in polynomials:
            for poly_q in polynomials:
                poly_d = gcd_poly(list(poly_p), list(poly_q), p)
                degree = len(poly_d) - 1
                common = [
                    word for word in words
                    if not any(apply_poly(poly_p, word, p))
                    and not any(apply_poly(poly_q, word, p))
                ]
                d_kernel = [
                    word for word in words
                    if not any(apply_poly(poly_d, word, p))
                ]
                require(common == d_kernel)
                kernel_equalities += 1
                seen = {}
                for word in common:
                    require(not any(apply_poly(poly_d, word, p)))
                    prefix = word[:degree]
                    require(prefix not in seen or seen[prefix] == word)
                    seen[prefix] = word
                pair_cells += 1
                common_words += len(common)
                prefix_pairs += len(common) * (len(common) - 1) // 2
    return pair_cells, kernel_equalities, common_words, prefix_pairs


def audit_integral_and_monic_hostiles():
    # On a length-one word S=1. In Z/2, P(S)=0 and Q(S)=2=0.
    word = (1,)
    require(apply_poly((-1, 1), word, 2) == (0,))
    require(apply_poly((1, 1), word, 2) == (0,))
    poly_p = (Fraction(-1), Fraction(1))
    poly_q = (Fraction(1), Fraction(1))
    rational_bezout = tuple(
        -Fraction(1, 2) * poly_p[i] + Fraction(1, 2) * poly_q[i]
        for i in range(2)
    )
    require(rational_bezout == (Fraction(1), Fraction(0)))
    require(2 * word[0] % 2 == 0 and word != (0,))

    # Monicity is separate. Over Z/4, D(T)=2T has two kernel words with the
    # same one-entry prefix.
    zero = (0, 0)
    spike = (0, 2)
    require(zero[:1] == spike[:1] and zero != spike)
    require(not any(apply_poly((0, 2), zero, 4)))
    require(not any(apply_poly((0, 2), spike, 4)))

    # A resultant prime permits but need not cause degeneration. P=T and
    # Q=T+5 have resultant 5; modulo 5 both are the invertible cyclic shift.
    words = list(product(range(5), repeat=3))
    common = [
        candidate for candidate in words
        if not any(apply_poly((0, 1), candidate, 5))
        and not any(apply_poly((5, 1), candidate, 5))
    ]
    require(common == [(0, 0, 0)])


def audit_empty_candidate_edge():
    # For C=empty the difference shell is empty, not {0}; THM-4274 assumes
    # C nonempty.
    observed = set()
    differences = {left - right for left in observed for right in observed}
    require(differences == set())
    require(differences != {0})
    require(differences <= {0})


def audit_integer_box_threshold():
    cells = 0
    injective_cells = 0
    hostile_cells = 0
    for height in range(1, 13):
        values = list(range(-height, height + 1))
        for modulus in range(2, 31):
            residues = [value % modulus for value in values]
            injective = len(set(residues)) == len(values)
            require(injective == (modulus > 2 * height))
            cells += 1
            injective_cells += int(injective)
            hostile_cells += int(not injective)
    return cells, injective_cells, hostile_cells


def fibonacci_word(first, second, length):
    word = [first, second]
    while len(word) < length:
        word.append(word[-1] + word[-2])
    return tuple(word[:length])


def audit_confined_fibres_and_density():
    # A complete order-two recurrence observer, followed by a separating
    # reduction. Candidate duplicates model full-certificate multiplicity.
    height, modulus, length, raw_bound = 1, 5, 8, 3
    candidates = []
    for first, second in product(range(-height, height + 1), repeat=2):
        word = fibonacci_word(first, second, length)
        multiplicity = 1 + ((first + 2 * second) % raw_bound)
        for duplicate in range(multiplicity):
            candidates.append((word, duplicate))

    def full(candidate):
        return candidate[0]

    def reduced_prefix(candidate):
        return tuple(value % modulus for value in candidate[0][:2])

    full_fibres = {}
    quotient_fibres = {}
    for candidate in candidates:
        full_fibres.setdefault(full(candidate), []).append(candidate)
        quotient_fibres.setdefault(reduced_prefix(candidate), []).append(candidate)

    partition_pair_checks = 0
    for candidate in candidates:
        for other in candidates:
            require(
                (full(candidate) == full(other))
                == (reduced_prefix(candidate) == reduced_prefix(other))
            )
            partition_pair_checks += 1
    require(max(map(len, quotient_fibres.values())) == raw_bound)

    targets = sorted(quotient_fibres)
    subset_checks = 0
    weights = {target: len(quotient_fibres[target]) for target in targets}
    rho = {
        target: Fraction(len(targets) * weights[target], len(candidates))
        for target in targets
    }
    sharp = max(rho.values())
    require(
        sharp
        <= Fraction(raw_bound * len(targets), len(candidates))
        <= raw_bound
    )
    for mask in range(1 << len(targets)):
        chosen = {
            targets[index]
            for index in range(len(targets))
            if (mask >> index) & 1
        }
        source_mass = Fraction(
            sum(weights[target] for target in chosen), len(candidates)
        )
        target_mass = Fraction(len(chosen), len(targets))
        exact_weight_average = Fraction(1, len(targets)) * sum(
            (rho[target] for target in chosen), start=Fraction(0, 1)
        )
        require(source_mass == exact_weight_average)
        require(source_mass <= sharp * target_mass <= raw_bound * target_mass)
        subset_checks += 1
    return len(candidates), len(targets), partition_pair_checks, subset_checks


def audit_specialization_and_unbounded_order():
    # Store a*u+b*f. Evaluation u=f maps (a,b) to (a+b)f; the first normal
    # Hasse jet stores (value coefficient, d/du coefficient).
    u_term = (1, 0)
    f_term = (0, 1)
    evaluate = lambda value: value[0] + value[1]
    jet = lambda value: (evaluate(value), value[0])
    require(evaluate(u_term) == evaluate(f_term))
    require(jet(u_term) != jet(f_term))

    cases = 0
    for prefix_length in range(1, 10):
        length = prefix_length + 1
        zero = (0,) * length
        spike = (0,) * prefix_length + (1,)
        require(zero[:prefix_length] == spike[:prefix_length])
        annihilator = (-1,) + (0,) * (length - 1) + (1,)
        require(not any(apply_poly(annihilator, zero, 2)))
        require(not any(apply_poly(annihilator, spike, 2)))
        cases += 1
    return cases


def audit_hazards_and_heavy_fibres():
    targets = tuple(range(8))
    event = {0}
    survivors = set(targets)
    hazards = []
    for _ in range(6):
        hit = len(survivors & event)
        hazards.append(Fraction(hit, len(survivors)))
        survivors -= event
    require(hazards == [Fraction(1, 8)] + [Fraction(0, 1)] * 5)
    product_survival = Fraction(1, 1)
    for hazard in hazards:
        product_survival *= 1 - hazard
    require(product_survival == Fraction(7, 8))
    require(product_survival != Fraction(7, 8) ** 6)

    heavy_rows = []
    for size in range(4, 21):
        source_mass = Fraction(2**size, 2**size + size - 1)
        target_mass = Fraction(1, size)
        require(source_mass > Fraction(4, 5) if size >= 6 else source_mass > 0)
        heavy_rows.append((size, target_mass, source_mass))
    require(heavy_rows[-1][1] < Fraction(1, 10))
    require(heavy_rows[-1][2] > Fraction(999, 1000))
    return len(hazards), len(heavy_rows)


def audit_capacity_and_balanced_overload():
    collision_cells = 0
    balanced_cells = 0
    for alphabet_size in range(2, 6):
        for samples in range(0, 6):
            codes = list(product(range(alphabet_size), repeat=samples))
            capacity = alphabet_size**samples
            require(len(codes) == capacity)
            overloaded = [codes[index % capacity] for index in range(capacity + 1)]
            require(len(set(overloaded)) < len(overloaded))
            collision_cells += 1

            multiplicity = alphabet_size + samples + 1
            balanced = [code for code in codes for _ in range(multiplicity)]
            fibre_counts = {code: balanced.count(code) for code in codes}
            require(set(fibre_counts.values()) == {multiplicity})
            require(
                all(capacity * weight == len(balanced)
                    for weight in fibre_counts.values())
            )
            balanced_cells += 1
    return collision_cells, balanced_cells


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")
    field = audit_field_gcd_recurrence()
    audit_integral_and_monic_hostiles()
    audit_empty_candidate_edge()
    boxes = audit_integer_box_threshold()
    density = audit_confined_fibres_and_density()
    unbounded = audit_specialization_and_unbounded_order()
    hazards = audit_hazards_and_heavy_fibres()
    capacity = audit_capacity_and_balanced_overload()
    print("confined_confluent_observer_transport_audit_thm4274")
    print(f"field_gcd_pair_cells={field[0]}")
    print(f"field_gcd_kernel_equalities={field[1]}")
    print(f"field_gcd_common_words={field[2]}")
    print(f"field_gcd_prefix_pair_controls={field[3]}")
    print(f"integer_box_cells={boxes[0]}")
    print(f"integer_box_injective={boxes[1]}")
    print(f"integer_box_hostile={boxes[2]}")
    print(f"bounded_fibre_candidates={density[0]}")
    print(f"bounded_fibre_target_words={density[1]}")
    print(f"fibre_partition_pair_checks={density[2]}")
    print(f"bounded_fibre_subset_checks={density[3]}")
    print(f"unbounded_order_hostiles={unbounded}")
    print(f"merged_hazard_steps={hazards[0]}")
    print(f"heavy_fibre_rows={hazards[1]}")
    print(f"capacity_collision_cells={capacity[0]}")
    print(f"balanced_overload_cells={capacity[1]}")
    print("integral_torsion_hostile=PASS")
    print("nonmonic_prefix_hostile=PASS")
    print("resultant_bad_prime_no_degeneration=PASS")
    print("empty_candidate_edge=PASS")
    print("graph_specialization_hasse_repair=PASS")
    print("status=PASS")


if __name__ == "__main__":
    main()
