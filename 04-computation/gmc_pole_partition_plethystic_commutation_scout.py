#!/usr/bin/env python3
"""Exact pole x partition commutation scout.

Multiplication of a row OGF by (1-Mt) is virtual-letter subtraction X -> X-M,
not literal deletion only from alphabets which contain M.  This companion
checks the factorial-normalized labelled-deletion square termwise in both
THM-3110 banks on three small supports and records the sharp literal-deletion
hostile at support (1,2).
"""

import ast
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from math import factorial, prod
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
UPSTREAM = HERE / "gmc_product_gamma_arbitrary_anchored_schur_thm3110.py"


def load_bank_prefix():
    tree = ast.parse(UPSTREAM.read_text(encoding="utf-8"))
    prefix = []
    for node in tree.body:
        prefix.append(node)
        if (isinstance(node, ast.Assign)
                and any(isinstance(target, ast.Name) and target.id == "BANKS"
                        for target in node.targets)):
            break
    namespace = {}
    exec(compile(ast.Module(body=prefix, type_ignores=[]),
                 str(UPSTREAM), "exec"), namespace)
    return namespace["BANKS"]


BANKS = load_bank_prefix()
SUPPORTS = ((1, 2), (1, 3), (2, 3))
MAXIMUM_DEGREE = 10


@lru_cache(maxsize=None)
def partitions(total, largest=None):
    if total == 0:
        return ((),)
    if largest is None or largest > total:
        largest = total
    return tuple(
        (first,) + tail
        for first in range(largest, 0, -1)
        for tail in partitions(total - first, first)
    )


ALL_SHAPES = tuple(sorted(
    (shape for degree in range(1, MAXIMUM_DEGREE + 1)
     for shape in partitions(degree)),
    key=lambda shape: (sum(shape), len(shape), shape),
))


def weight(shape):
    return prod(factorial(part) for part in shape)


def residual_roots(invariant, row, a_value, b_value):
    middle = max(2 * a_value, b_value)
    counts = Counter()
    for alpha, beta in row:
        counts.update(range(alpha * a_value + beta * b_value))
    counts.subtract(3 * list(range(a_value)))
    counts.subtract((invariant + 1) * list(range(middle)))
    require(all(value >= 0 for value in counts.values()),
            "residual-root subtraction became negative")
    return tuple(sorted(counts.elements()))


def all_monomial_values(power_sums):
    values = {(): 1}
    for shape in ALL_SHAPES:
        exponent = shape[-1]
        remainder = shape[:-1]
        value = power_sums[exponent - 1] * values[remainder]
        for old_exponent in set(remainder):
            merged = list(remainder)
            merged.remove(old_exponent)
            merged.append(old_exponent + exponent)
            merged = tuple(sorted(merged, reverse=True))
            value -= Counter(merged)[old_exponent + exponent] * values[merged]
        multiplicity = Counter(shape)[exponent]
        require(value % multiplicity == 0,
                "virtual monomial recurrence lost integrality")
        values[shape] = value // multiplicity
    return values


def virtual_values(roots, removed):
    power_sums = tuple(
        sum(root**degree for root in roots)
        - sum(root**degree for root in removed)
        for degree in range(1, MAXIMUM_DEGREE + 1)
    )
    return all_monomial_values(power_sums)


def complete(values, degree):
    return sum(values[shape] for shape in partitions(degree))


def labelled_down(values, source_degree):
    answer = {shape: Fraction(0) for shape in partitions(source_degree - 1)}
    for shape in partitions(source_degree):
        normalized = Fraction(values[shape], weight(shape))
        for part, multiplicity in Counter(shape).items():
            target = list(shape)
            target.remove(part)
            if part > 1:
                target.append(part - 1)
            target = tuple(sorted(target, reverse=True))
            answer[target] += part * multiplicity * normalized
    return answer


termwise_commutations = row_strip_checks = bank_commutations = 0
literal_absence = {0: 0, 1: 0}
literal_presence = {0: 0, 1: 0}

for invariant, bank in enumerate(BANKS):
    for a_value, b_value in SUPPORTS:
        atom_data = []
        for coefficient, row in bank:
            roots = residual_roots(invariant, row, a_value, b_value)
            # The largest possible residual root is enough for the universal
            # commutation check.  At (1,2), THM-3120's reduced top pole is 5.
            pole = max(root for _, current in bank
                       for root in residual_roots(
                           invariant, current, a_value, b_value))
            plain = virtual_values(roots, ())
            stripped = virtual_values(roots, (pole,))
            atom_data.append((coefficient, roots, pole, plain, stripped))

            for degree in range(5, MAXIMUM_DEGREE + 1):
                require(
                    complete(stripped, degree)
                    == complete(plain, degree)
                    - pole * complete(plain, degree - 1),
                    "row pole strip is not virtual-letter subtraction",
                )
                row_strip_checks += 1

            for target_degree in range(5, MAXIMUM_DEGREE):
                left = labelled_down(stripped, target_degree + 1)
                for shape in partitions(target_degree):
                    right = Fraction(
                        (sum(roots) - pole) * stripped[shape], weight(shape))
                    require(left[shape] == right,
                            "labelled deletion and pole subtraction do not commute")
                termwise_commutations += 1

            if (a_value, b_value) == (1, 2):
                if 5 in roots:
                    literal_presence[invariant] += 1
                else:
                    literal_absence[invariant] += 1
                    require(complete(plain, 4) > 0,
                            "literal-deletion hostile has zero lower response")

        for target_degree in range(5, MAXIMUM_DEGREE):
            shapes = partitions(target_degree)
            aggregate_left = {shape: Fraction(0) for shape in shapes}
            aggregate_right = {shape: Fraction(0) for shape in shapes}
            for coefficient, roots, pole, _, stripped in atom_data:
                down = labelled_down(stripped, target_degree + 1)
                for shape in shapes:
                    aggregate_left[shape] += coefficient * down[shape]
                    aggregate_right[shape] += Fraction(
                        coefficient * (sum(roots) - pole) * stripped[shape],
                        weight(shape),
                    )
            require(aggregate_left == aggregate_right,
                    "signed bank commutation failed")
            bank_commutations += 1

require(literal_absence == {0: 21, 1: 21},
        "literal top-pole absence census drift")
require(literal_presence == {0: 3, 1: 4},
        "literal top-pole presence census drift")
require((termwise_commutations, row_strip_checks, bank_commutations)
        == (735, 882, 30), "commutation census drift")

print("pole_operator=(1-Mt)=plethystic_virtual_subtraction:X->X-M")
print("factorial_vector=Mtilde_N[X]=sum_lambda m_lambda[X]/w_lambda e_lambda")
print("mixed_square=A*P_M=P_M*A:exact_termwise_before_bank_sum")
print("common_value=(p1[X]-M)*Mtilde_N[X-M]")
print(f"termwise_commutations={termwise_commutations}")
print(f"row_strip_checks={row_strip_checks}")
print(f"bank_commutations={bank_commutations}")
print("literal_deletion_hostile_support=(1,2):M=5")
print(
    f"I1_M_absent={literal_absence[0]}/{len(BANKS[0])}:"
    f"I2_M_absent={literal_absence[1]}/{len(BANKS[1])}"
)
print("scope=algebraic_commutation_only;virtual_subtraction_is_signed")
print("all_exact_checks=PASS")
