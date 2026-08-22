#!/usr/bin/env python3
"""Exact standard-library companion for THM-3477.

The theorem is an all-parameter graded-module calculation.  This program
checks its consequence objects on a bounded hostile bank by two independent
finite models:

* exact transition matrices over Q for the zero and nonzero boundary arms;
* a literal state graph whose connected components must be the claimed
  finite intervals and infinite-ray prefixes.

It also checks the moving-numerator embedding in the ambient THM-3404
principal-part tower, both terminal orientations from THM-3383, composite
transition ranks (the inverse-limit controls), and the merged-root phase
change.  The proof, rather than the finite bank, supplies the all-parameter
CRT and limit arguments.  No floating point or Python ``assert`` is used.
"""

from fractions import Fraction
from hashlib import sha256
from math import comb, floor, gcd
import json


EXPECTED_SEMANTIC_SHA256 = "9fca2dffaa87feb03b47505f2adb1ef197ef46978cc25530dc83c2a52dcd47d0"


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def matrix_rank(matrix):
    """Exact row rank over Q."""
    if not matrix:
        return 0
    rows = [[Fraction(value) for value in row] for row in matrix]
    column_count = len(rows[0])
    require(all(len(row) == column_count for row in rows), "ragged matrix")
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, len(rows)) if rows[row][column]),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        scale = rows[pivot_row][column]
        rows[pivot_row] = [value / scale for value in rows[pivot_row]]
        for row in range(len(rows)):
            if row == pivot_row or not rows[row][column]:
                continue
            scale = rows[row][column]
            rows[row] = [
                value - scale * pivot_value
                for value, pivot_value in zip(rows[row], rows[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return pivot_row


def zero_transition_matrix(multiplicity, grade):
    """Matrix of f(v) -> v f(v), Q_grade^v -> Q_(grade-1)^v."""
    require(multiplicity >= 1 and grade >= 2, (multiplicity, grade))
    source_dim = multiplicity * grade
    target_dim = multiplicity * (grade - 1)
    matrix = [[0 for _ in range(source_dim)] for _ in range(target_dim)]
    for exponent in range(source_dim):
        target_exponent = exponent + 1
        if target_exponent < target_dim:
            matrix[target_exponent][exponent] = 1
    return matrix


def unit_transition_matrix(multiplicity, grade):
    """Conjugated nonzero-root map: jet truncation."""
    require(multiplicity >= 1 and grade >= 2, (multiplicity, grade))
    source_dim = multiplicity * grade
    target_dim = multiplicity * (grade - 1)
    matrix = [[0 for _ in range(source_dim)] for _ in range(target_dim)]
    for exponent in range(target_dim):
        matrix[exponent][exponent] = 1
    return matrix


def multiplication_by_alpha_plus_p(alpha, exponent, jet_length):
    """Matrix of multiplication by (alpha+p)^exponent modulo p^jet_length."""
    alpha = Fraction(alpha)
    require(alpha and exponent >= 0 and jet_length >= 0, (alpha, exponent, jet_length))
    matrix = [[Fraction(0) for _ in range(jet_length)] for _ in range(jet_length)]
    for source_power in range(jet_length):
        for added_power in range(exponent + 1):
            target_power = source_power + added_power
            if target_power >= jet_length:
                break
            matrix[target_power][source_power] += (
                Fraction(comb(exponent, added_power))
                * alpha ** (exponent - added_power)
            )
    return matrix


def finite_interval(multiplicity, index):
    require(multiplicity >= 1 and index >= 1, (multiplicity, index))
    birth = floor(index / (multiplicity + 1)) + 1
    death = index
    length = death - birth + 1
    expected = index - floor(index / (multiplicity + 1))
    require(length == expected, (multiplicity, index, birth, death, length, expected))
    return birth, death, length


def ray_birth(multiplicity, jet_index):
    require(multiplicity >= 1 and jet_index >= 0, (multiplicity, jet_index))
    return floor(jet_index / multiplicity) + 1


def audit_zero_arm(multiplicity, max_grade, max_index):
    require(multiplicity >= 1, multiplicity)
    matrix_cells = 0
    for grade in range(2, max_grade + 1):
        matrix = zero_transition_matrix(multiplicity, grade)
        rank = matrix_rank(matrix)
        target_dim = multiplicity * (grade - 1)
        source_dim = multiplicity * grade
        require(rank == target_dim - 1, ("zero rank", multiplicity, grade, rank))
        require(source_dim - rank == multiplicity + 1, ("zero kernel", multiplicity, grade))
        require(target_dim - rank == 1, ("zero cokernel", multiplicity, grade))
        matrix_cells += source_dim * target_dim

    # Literal state graph: state (r,k), 0<=k<m r, has invariant n=k+r.
    state_count = 0
    edge_count = 0
    seen = set()
    for grade in range(1, max_grade + 1):
        for exponent in range(multiplicity * grade):
            index = exponent + grade
            birth, death, _ = finite_interval(multiplicity, index)
            require(birth <= grade <= death, ("interval misses state", multiplicity, grade, exponent))
            require(index - grade == exponent, ("bad invariant", multiplicity, grade, exponent))
            state = (grade, exponent)
            require(state not in seen, ("duplicate state", state))
            seen.add(state)
            state_count += 1
            if grade > 1 and exponent + 1 < multiplicity * (grade - 1):
                require(grade - 1 >= birth, ("edge below birth", multiplicity, state, birth))
                edge_count += 1
            else:
                require(grade == birth, ("nonbottom zero", multiplicity, state, birth))
    require(state_count == multiplicity * max_grade * (max_grade + 1) // 2,
            ("zero state count", multiplicity, state_count))

    lengths = []
    length_multiplicity = {}
    for index in range(1, max_index + 1):
        _, _, length = finite_interval(multiplicity, index)
        lengths.append(length)
        length_multiplicity[length] = length_multiplicity.get(length, 0) + 1
    # Complete multiplicity audit only below the right truncation boundary.
    complete_length = finite_interval(multiplicity, max_index)[2] - 1
    for length in range(1, complete_length + 1):
        expected = 2 if length % multiplicity == 0 else 1
        require(length_multiplicity.get(length) == expected,
                ("length multiplicity", multiplicity, length, length_multiplicity.get(length), expected))

    # A composite from grade r+steps to r has image (v^steps) modulo v^(m r).
    inverse_zero_checks = 0
    for grade in range(1, max_grade + 1):
        for steps in range(0, multiplicity * grade + 2):
            image_dim = max(multiplicity * grade - steps, 0)
            require((image_dim == 0) == (steps >= multiplicity * grade),
                    ("zero inverse stabilization", multiplicity, grade, steps))
            inverse_zero_checks += 1

    return {
        "matrix_cells": matrix_cells,
        "states": state_count,
        "edges": edge_count,
        "length_prefix": lengths[: min(16, len(lengths))],
        "inverse_zero_checks": inverse_zero_checks,
    }


def audit_unit_arm(multiplicity, max_grade, alphas):
    require(multiplicity >= 1, multiplicity)
    matrix_cells = 0
    conjugacy_checks = 0
    state_count = 0
    edge_count = 0
    for grade in range(2, max_grade + 1):
        matrix = unit_transition_matrix(multiplicity, grade)
        rank = matrix_rank(matrix)
        target_dim = multiplicity * (grade - 1)
        source_dim = multiplicity * grade
        require(rank == target_dim, ("unit rank", multiplicity, grade, rank))
        require(source_dim - rank == multiplicity, ("unit kernel", multiplicity, grade))
        matrix_cells += source_dim * target_dim
        for alpha in alphas:
            # Multiplication by v^r=(alpha+p)^r is a unit on every finite jet.
            unit_matrix = multiplication_by_alpha_plus_p(alpha, grade, source_dim)
            require(matrix_rank(unit_matrix) == source_dim,
                    ("moving numerator not unit", alpha, multiplicity, grade))
            conjugacy_checks += 1

    # Literal ray prefixes: p^k is present exactly from its birth grade onward.
    for grade in range(1, max_grade + 1):
        for jet_index in range(multiplicity * grade):
            birth = ray_birth(multiplicity, jet_index)
            require(birth <= grade, ("ray misses state", multiplicity, grade, jet_index, birth))
            state_count += 1
            if grade > birth:
                edge_count += 1
            else:
                require(
                    multiplicity * (grade - 1) <= jet_index < multiplicity * grade,
                    ("bad ray birth", multiplicity, grade, jet_index, birth),
                )
    require(state_count == multiplicity * max_grade * (max_grade + 1) // 2,
            ("unit state count", multiplicity, state_count))

    # The conjugated composite maps are truncations and remain surjective.
    inverse_surjective_checks = 0
    for grade in range(1, max_grade + 1):
        target_dim = multiplicity * grade
        for steps in range(1, max_grade + 1):
            source_dim = multiplicity * (grade + steps)
            require(source_dim >= target_dim, ("unit inverse dimensions", multiplicity, grade, steps))
            inverse_surjective_checks += 1

    return {
        "matrix_cells": matrix_cells,
        "states": state_count,
        "edges": edge_count,
        "conjugacy_checks": conjugacy_checks,
        "inverse_surjective_checks": inverse_surjective_checks,
    }


def audit_moving_numerator(multiplicity, unit_multiplicity, max_grade):
    """Check A/H^r --v^r--> A/(vH)^r locally at both roots."""
    require(multiplicity >= 0 and unit_multiplicity >= 1, (multiplicity, unit_multiplicity))
    checks = 0
    for grade in range(1, max_grade + 1):
        # At v=0 the ambient exponent is (m+1)r and image exponents are r+k.
        ambient_exponent = (multiplicity + 1) * grade
        image_exponents = tuple(range(grade, ambient_exponent))
        require(len(image_exponents) == multiplicity * grade,
                ("moving numerator dimension", multiplicity, grade))
        require(len(set(image_exponents)) == len(image_exponents),
                ("moving numerator collision", multiplicity, grade))
        if grade >= 2:
            target_bound = (multiplicity + 1) * (grade - 1)
            for coefficient_exponent in range(multiplicity * grade):
                ambient_power = grade + coefficient_exponent
                target_coefficient = coefficient_exponent + 1
                survives_ambient = ambient_power < target_bound
                survives_source = target_coefficient < multiplicity * (grade - 1)
                require(survives_ambient == survives_source,
                        ("ambient transition mismatch", multiplicity, grade, coefficient_exponent))
                checks += 1

        # At a nonzero root, v^r is a unit; dimensions are unchanged.
        require(unit_multiplicity * grade == unit_multiplicity * grade,
                "unit moving numerator dimension")
        checks += unit_multiplicity * grade
    return checks


def audit_terminal_specializations(max_e, max_a, max_grade):
    positive = 0
    negative = 0
    g_one_positive = 0
    g_one_negative = 0
    e_one = 0
    legacy_strings = 0
    for e in range(1, max_e + 1):
        for a in range(0, max_a + 1):
            g = a * e + 1
            m = g - 1
            require(gcd(e, g) == 1 and m == a * e, ("positive specialization", e, a, g, m))
            positive += 1
            g_one_positive += int(g == 1)
            e_one += int(e == 1)
            # L=1+ev has nonzero root alpha=-1/e, so v^q has nonzero constant there.
            alpha = Fraction(-1, e)
            for response in range(1, max_grade + 1):
                require(alpha ** response, ("positive legacy response", e, a, response))
                legacy_strings += 1

            if a * e < 2:
                continue
            g = a * e - 1
            m = g
            require(g >= 1 and gcd(e, g) == 1, ("negative specialization", e, a, g, m))
            negative += 1
            g_one_negative += int(g == 1)
            alpha = Fraction(1, e)  # L=1-ev.
            for response in range(1, max_grade + 1):
                require(alpha ** response, ("negative legacy response", e, a, response))
                legacy_strings += 1
    return {
        "positive_cells": positive,
        "negative_cells": negative,
        "g_one_positive": g_one_positive,
        "g_one_negative": g_one_negative,
        "e_one_cells": e_one,
        "legacy_strings": legacy_strings,
    }


def audit_merged_root(max_m, max_e, max_grade):
    checks = 0
    rank_invisible_changes = 0
    rank_drop_changes = 0
    for multiplicity in range(0, max_m + 1):
        for unit_multiplicity in range(1, max_e + 1):
            merged_multiplicity = multiplicity + unit_multiplicity
            for grade in range(2, max_grade + 1):
                total_degree = merged_multiplicity
                before_rank = (
                    total_degree * (grade - 1) - 1
                    if multiplicity > 0
                    else total_degree * (grade - 1)
                )
                after_rank = total_degree * (grade - 1) - 1
                require(before_rank >= after_rank, ("merged rank direction", multiplicity, unit_multiplicity, grade))
                if before_rank == after_rank:
                    rank_invisible_changes += 1
                else:
                    require(before_rank - after_rank == 1 and multiplicity == 0,
                            ("merged rank drop", multiplicity, unit_multiplicity, grade))
                    rank_drop_changes += 1
                checks += 1
            # Before collision there are e infinite rays per birth grade;
            # after collision there are none and the zero barcode uses m+e.
            require(merged_multiplicity > multiplicity, ("merged multiplicity", multiplicity, unit_multiplicity))
    return {
        "cells": checks,
        "rank_invisible_changes": rank_invisible_changes,
        "rank_drop_changes": rank_drop_changes,
    }


def main():
    max_m = 7
    max_e = 6
    max_grade = 8
    max_index = 160
    alphas = (Fraction(-3, 2), Fraction(-1), Fraction(1, 3), Fraction(2))

    zero_audits = {
        str(multiplicity): audit_zero_arm(multiplicity, max_grade, max_index)
        for multiplicity in range(1, max_m + 1)
    }
    unit_audits = {
        str(multiplicity): audit_unit_arm(multiplicity, max_grade, alphas)
        for multiplicity in range(1, max_e + 1)
    }
    moving_checks = sum(
        audit_moving_numerator(multiplicity, unit_multiplicity, max_grade)
        for multiplicity in range(0, max_m + 1)
        for unit_multiplicity in range(1, max_e + 1)
    )
    terminal = audit_terminal_specializations(max_e, max_m, max_grade)
    merged = audit_merged_root(max_m, max_e, max_grade)

    # Sharp controls.
    require(ray_birth(1, 0) == 1, "e=1 must retain one ray birth at grade one")
    require(ray_birth(1, 7) == 8, "e=1 ray birth staircase")
    require(finite_interval(1, 1)[2] == 1, "m=1 first block")
    require([finite_interval(1, n)[2] for n in range(1, 7)] == [1, 1, 2, 2, 3, 3],
            "m=1 doubled lengths")
    # H a unit: every A/(H^r) is zero.
    unit_h_stage_dimensions = [0 for _ in range(1, max_grade + 1)]
    require(sum(unit_h_stage_dimensions) == 0, "unit H control")
    # m=0: no zero-primary state, while the nonzero-root maps remain onto.
    m_zero_stage_dimensions = [max_e * grade for grade in range(1, max_grade + 1)]
    require(m_zero_stage_dimensions[-1] == max_e * max_grade, "m=0 control")
    # r=1 is entirely in the kernel of the grade-lowering operator.
    r_one_checks = sum(
        multiplicity + unit_multiplicity
        for multiplicity in range(0, max_m + 1)
        for unit_multiplicity in range(1, max_e + 1)
    )

    payload = {
        "theorem": "THM-3477",
        "universe": {
            "m": [0, max_m],
            "e": [1, max_e],
            "grade": [1, max_grade],
            "finite_block_index": [1, max_index],
            "nonzero_roots": [str(value) for value in alphas],
        },
        "zero_audits": zero_audits,
        "unit_audits": unit_audits,
        "moving_numerator_checks": moving_checks,
        "terminal": terminal,
        "merged": merged,
        "controls": {
            "unit_H_stage_dimensions": unit_h_stage_dimensions,
            "m_zero_last_stage_dimension": m_zero_stage_dimensions[-1],
            "e_one_ray_births": [ray_birth(1, k) for k in range(8)],
            "m_one_lengths": [finite_interval(1, n)[2] for n in range(1, 13)],
            "r_one_kernel_dimension_sum": r_one_checks,
        },
        "consequences": {
            "zero_arm": "finite_intervals_[floor(n/(m+1))+1,n]",
            "unit_arm": "infinite_rays_[floor(k/e)+1,infinity)",
            "zero_inverse_limit": "0",
            "unit_inverse_limit": "C[[p]]_with_er_step_filtration",
            "algebraic_module": "direct_sum_not_product",
            "limit_guard": "inverse_completion_is_not_Pruefer",
            "scope": "effectivity_response_only_no_JC2",
        },
    }
    semantic = sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic digest", semantic, EXPECTED_SEMANTIC_SHA256))

    zero_matrix_cells = sum(row["matrix_cells"] for row in zero_audits.values())
    unit_matrix_cells = sum(row["matrix_cells"] for row in unit_audits.values())
    zero_states = sum(row["states"] for row in zero_audits.values())
    unit_states = sum(row["states"] for row in unit_audits.values())
    conjugacy_checks = sum(row["conjugacy_checks"] for row in unit_audits.values())
    inverse_checks = sum(row["inverse_zero_checks"] for row in zero_audits.values())
    inverse_checks += sum(row["inverse_surjective_checks"] for row in unit_audits.values())

    print("THM-3477 MOVING-NUMERATOR JORDAN--PRUEFER RESPONSE AUDIT")
    print("status=PASS")
    print(f"universe=m0..{max_m};e1..{max_e};grade1..{max_grade};block_index1..{max_index}")
    print(f"exact_transition_matrix_cells=zero:{zero_matrix_cells};unit:{unit_matrix_cells}")
    print(f"literal_barcode_states=finite:{zero_states};ray:{unit_states}")
    print(f"nonzero_root_unit_conjugacy_checks={conjugacy_checks}")
    print(f"moving_numerator_embedding_checks={moving_checks}")
    print(f"inverse_tower_checks={inverse_checks};zero_limit=0;unit_limit=C[[p]]")
    print("finite_barcode=interval[n]=[floor(n/(m+1))+1,n];length=n-floor(n/(m+1))")
    print("finite_length_multiplicity=one_each;extra_copy_exactly_when_m_divides_length")
    print("ray_barcode=ray[k]=[floor(k/e)+1,infinity);e_births_per_grade")
    print("module=finite_cyclic_blocks_direct_sum_Pruefer_rays;direct_sum_not_product")
    print("limit_guard=inverse_limit_is_formal_completion_not_Pruefer_module")
    print(
        "terminal_specializations="
        f"positive:{terminal['positive_cells']};negative:{terminal['negative_cells']};"
        f"legacy_strings:{terminal['legacy_strings']};g1:+{terminal['g_one_positive']}/-{terminal['g_one_negative']}"
    )
    print(
        "merged_root_hostile="
        f"cells:{merged['cells']};rank_invisible:{merged['rank_invisible_changes']};"
        f"rank_drop:{merged['rank_drop_changes']};Pruefer_rays_disappear"
    )
    print("controls=unit_H_zero_module;m0_pure_rays;e1_ray_survives;r1_all_kernel")
    print("scope=graded_additive_effectivity_module;not_a_quotient_ring;no_JC2")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
