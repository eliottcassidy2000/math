#!/usr/bin/env python3
"""Exact companion for THM-3154's sharp complete partition-profile jet."""

from collections import Counter, defaultdict
from functools import lru_cache
from hashlib import sha256
from itertools import product
from math import comb, factorial


MAXIMUM_PROFILE_DEGREE = 30
MAXIMUM_JET_DEGREE = 16
MAXIMUM_KERNEL_DEGREE = 12


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


@lru_cache(maxsize=None)
def partitions_bounded(total, maximum):
    if total == 0:
        return ((),)
    answer = []
    for part in range(min(total, maximum), 0, -1):
        for tail in partitions_bounded(total - part, part):
            answer.append((part,) + tail)
    return tuple(answer)


def partitions(total):
    return partitions_bounded(total, total)


def cutoff(degree):
    return (degree - 2) // 2


def profile_code(shape, maximum_marked):
    return ((len(shape),)
            + tuple(shape.count(part)
                    for part in range(1, maximum_marked + 1)))


def reconstruct_from_code(degree, code, maximum_marked):
    length = code[0]
    counts = code[1:]
    small = []
    for part, multiplicity in enumerate(counts, start=1):
        small.extend([part] * multiplicity)
    residual_count = length - len(small)
    residual_sum = degree - sum(small)
    require(residual_count >= 0 and residual_sum >= 0,
            "profile code has negative residual data")
    if residual_count == 0:
        require(residual_sum == 0, "zero-count residual has positive mass")
        residual = []
    else:
        base = maximum_marked + 1
        surplus = residual_sum - residual_count * base
        require(surplus >= 0, "unmarked residual fell below the cutoff")
        if residual_count >= 2:
            require(surplus <= 1,
                    "faithful code has an ambiguous unmarked residual")
        residual = [base] * residual_count
        residual[0] += surplus
    return tuple(sorted(small + residual, reverse=True))


def faithful_code_audit():
    digest = sha256()
    code_checks = cutoff_checks = 0
    for degree in range(2, MAXIMUM_PROFILE_DEGREE + 1):
        maximum_marked = cutoff(degree)
        seen = {}
        for shape in partitions(degree):
            code = profile_code(shape, maximum_marked)
            require(code not in seen, "sharp profile code lost injectivity")
            require(reconstruct_from_code(
                degree, code, maximum_marked) == shape,
                "explicit profile reconstruction drift")
            seen[code] = shape
            digest.update(repr((degree, shape, code)).encode())
            code_checks += 1

        if degree < 4:
            require(maximum_marked == 0,
                    "low-degree length-only boundary drift")
            continue
        if degree % 2 == 0:
            half = degree // 2
            left = tuple(sorted((half - 1, half + 1), reverse=True))
            right = (half, half)
        else:
            half = degree // 2
            left = tuple(sorted((half - 1, half + 2), reverse=True))
            right = tuple(sorted((half, half + 1), reverse=True))
        require(left != right
                and profile_code(left, maximum_marked - 1)
                == profile_code(right, maximum_marked - 1)
                and profile_code(left, maximum_marked)
                != profile_code(right, maximum_marked),
                "sharp initial-cutoff witness drift")
        cutoff_checks += 1
    require(code_checks == 28627 and cutoff_checks == 27,
            "sharp profile census drift")
    return code_checks, cutoff_checks, digest.hexdigest()


def all_monomial_values(roots, removed, maximum_degree):
    """Evaluate every m_lambda on the virtual alphabet roots-removed."""

    power_sums = tuple(
        sum(root ** degree for root in roots)
        - sum(root ** degree for root in removed)
        for degree in range(1, maximum_degree + 1)
    )
    shapes = tuple(sorted(
        (shape
         for degree in range(1, maximum_degree + 1)
         for shape in partitions(degree)),
        key=lambda shape: (sum(shape), len(shape), shape),
    ))
    values = {(): 1}
    for shape in shapes:
        exponent = shape[-1]
        remainder = shape[:-1]
        value = power_sums[exponent - 1] * values[remainder]
        for old_exponent in set(remainder):
            merged = list(remainder)
            merged.remove(old_exponent)
            merged.append(old_exponent + exponent)
            merged = tuple(sorted(merged, reverse=True))
            value -= (
                Counter(merged)[old_exponent + exponent] * values[merged]
            )
        multiplicity = Counter(shape)[exponent]
        require(value % multiplicity == 0,
                "virtual monomial recurrence lost integrality")
        values[shape] = value // multiplicity
    return values


def profile_series(values, maximum_degree, maximum_marked):
    series = defaultdict(int)
    series[(0, 0) + (0,) * maximum_marked] = values[()]
    for degree in range(1, maximum_degree + 1):
        for shape in partitions(degree):
            key = ((degree, len(shape))
                   + tuple(shape.count(part)
                           for part in range(1, maximum_marked + 1)))
            series[key] += values[shape]
    return {key: value for key, value in series.items() if value}


def add_term(polynomial, key, value):
    polynomial[key] = polynomial.get(key, 0) + value
    if polynomial[key] == 0:
        del polynomial[key]


def multiply(left, right, maximum_degree):
    answer = {}
    for left_key, left_value in left.items():
        for right_key, right_value in right.items():
            if left_key[0] + right_key[0] > maximum_degree:
                continue
            key = tuple(a + b for a, b in zip(left_key, right_key))
            add_term(answer, key, left_value * right_value)
    return answer


def pole_denominator(pole, maximum_marked):
    """Raw-variable denominator in the marked virtual-pole kernel."""

    zero = (0, 0) + (0,) * maximum_marked
    answer = {zero: 1}

    def key(degree, length=0, marked_part=None):
        out = [0] * (maximum_marked + 2)
        out[0] = degree
        out[1] = length
        if marked_part is not None:
            out[marked_part + 1] = 1
        return tuple(out)

    add_term(answer, key(1, length=1), pole)
    add_term(answer, key(1), -pole)
    for part in range(1, maximum_marked + 1):
        add_term(answer, key(part, 1, part), pole ** part)
        add_term(answer, key(part, 1), -(pole ** part))
        add_term(answer, key(part + 1, 1, part), -(pole ** (part + 1)))
        add_term(answer, key(part + 1, 1), pole ** (part + 1))
    return answer


def virtual_pole_kernel_audit():
    roots = (1, 2, 4, 7)
    pole = 3
    before_values = all_monomial_values(
        roots, (), MAXIMUM_KERNEL_DEGREE)
    after_values = all_monomial_values(
        roots, (pole,), MAXIMUM_KERNEL_DEGREE)
    digest = sha256()
    coefficient_checks = 0
    for maximum_marked in range(0, 6):
        before = profile_series(
            before_values, MAXIMUM_KERNEL_DEGREE, maximum_marked)
        after = profile_series(
            after_values, MAXIMUM_KERNEL_DEGREE, maximum_marked)
        denominator = pole_denominator(pole, maximum_marked)
        left = multiply(after, denominator, MAXIMUM_KERNEL_DEGREE)
        right = dict(before)
        for key, value in tuple(before.items()):
            if key[0] + 1 <= MAXIMUM_KERNEL_DEGREE:
                add_term(right, (key[0] + 1,) + key[1:], -pole * value)
        require(left == right, "multivariate virtual-pole kernel drift")
        coefficient_checks += len(right)
        digest.update(repr((maximum_marked, sorted(right.items()))).encode())

    require(pole_denominator(pole, 0) == {
        (0, 0): 1, (1, 1): pole, (1, 0): -pole,
    }, "length-only kernel specialization drift")
    require(coefficient_checks == 902,
            "virtual-pole coefficient census drift")
    return coefficient_checks, digest.hexdigest()


def coarsenings(shape):
    answer = set()
    for first in range(len(shape)):
        for second in range(first + 1, len(shape)):
            merged = [shape[index] for index in range(len(shape))
                      if index not in (first, second)]
            merged.append(shape[first] + shape[second])
            answer.add(tuple(sorted(merged, reverse=True)))
    return answer


def endpoint_jet_audit():
    reconstruction_digest = sha256()
    reconstruction_checks = upset_checks = 0
    expected_upsets = {2: 3, 3: 4, 4: 7, 5: 10, 6: 27, 7: 47}
    for degree in range(2, MAXIMUM_JET_DEGREE + 1):
        maximum_marked = cutoff(degree)
        shapes = partitions(degree)
        codes = {shape: profile_code(shape, maximum_marked)
                 for shape in shapes}
        raw = {shape: (-1) ** index * (index + 1) ** 2
               for index, shape in enumerate(shapes)}
        raw[shapes[0]] -= sum(raw.values())
        require(sum(raw.values()) == 0, "test current lost zero mass")

        jets = defaultdict(int)
        for shape in shapes:
            value = raw[shape]
            code = codes[shape]
            for beta in product(*(range(entry + 1) for entry in code)):
                weight = 1
                for entry, order in zip(code, beta):
                    weight *= comb(entry, order)
                jets[beta] += value * weight

        recovered = {}
        for shape in shapes:
            code = codes[shape]
            value = 0
            for beta, jet in jets.items():
                if not all(order >= entry
                           for order, entry in zip(beta, code)):
                    continue
                weight = (-1) ** sum(
                    order - entry for order, entry in zip(beta, code))
                for order, entry in zip(beta, code):
                    weight *= comb(order, entry)
                value += weight * jet
            require(value == raw[shape],
                    "multivariate endpoint inversion drift")
            recovered[shape] = value
            reconstruction_digest.update(
                repr((degree, shape, code, value)).encode())
            reconstruction_checks += 1

        if degree in expected_upsets:
            index = {shape: position
                     for position, shape in enumerate(shapes)}
            edges = tuple(
                (index[fine], index[coarse])
                for fine in shapes for coarse in coarsenings(fine)
            )
            count = 0
            for mask in range(1 << len(shapes)):
                if not all(not (mask >> fine) & 1
                           or (mask >> coarse) & 1
                           for fine, coarse in edges):
                    continue
                direct = sum(raw[shape] for position, shape in enumerate(shapes)
                             if (mask >> position) & 1)
                rebuilt = sum(
                    recovered[shape]
                    for position, shape in enumerate(shapes)
                    if (mask >> position) & 1
                )
                require(direct == rebuilt, "upset mass reconstruction drift")
                count += 1
            require(count == expected_upsets[degree],
                    "partition-upset census drift")
            upset_checks += count

    require(reconstruction_checks == 913 and upset_checks == 98,
            "endpoint-jet census drift")
    return reconstruction_checks, upset_checks, \
        reconstruction_digest.hexdigest()


def central_class_audit():
    class_checks = factorial_checks = 0
    digest = sha256()
    alphabet = (1, 2, 4)
    for degree in range(2, MAXIMUM_JET_DEGREE + 1):
        class_mass = 0
        maximum_marked = cutoff(degree)
        for shape in partitions(degree):
            code = profile_code(shape, maximum_marked)
            rebuilt = reconstruct_from_code(degree, code, maximum_marked)
            require(rebuilt == shape, "cycle-type reconstruction drift")
            multiplicities = Counter(shape)
            centralizer = 1
            for cycle_length, multiplicity in multiplicities.items():
                centralizer *= (cycle_length ** multiplicity
                                * factorial(multiplicity))
            class_size = factorial(degree) // centralizer
            class_mass += class_size
            coefficient = 1
            for cycle_length in shape:
                coefficient *= sum(
                    value ** cycle_length for value in alphabet)
            rebuilt_coefficient = 1
            for cycle_length in rebuilt:
                rebuilt_coefficient *= sum(
                    value ** cycle_length for value in alphabet)
            require(coefficient == rebuilt_coefficient,
                    "central cycle-weight reconstruction drift")
            digest.update(repr(
                (degree, code, class_size, coefficient)).encode())
            class_checks += 1
        require(class_mass == factorial(degree),
                "conjugacy-class mass formula drift")
        factorial_checks += 1

    left = (4, 1, 1)
    right = (2, 2, 1, 1)
    require(sum(left) - left.count(1) == sum(right) - right.count(1) == 4
            and profile_code(left, cutoff(6))
            != profile_code(right, cutoff(6)),
            "equal-moved-support central hostile drift")
    require(class_checks == 913 and factorial_checks == 15,
            "central-class census drift")
    return class_checks, factorial_checks, digest.hexdigest()


def main():
    code_checks, cutoff_checks, code_hash = faithful_code_audit()
    kernel_checks, kernel_hash = virtual_pole_kernel_audit()
    jet_checks, upset_checks, jet_hash = endpoint_jet_audit()
    class_checks, factorial_checks, class_hash = central_class_audit()
    print(f"profile_degrees=2..{MAXIMUM_PROFILE_DEGREE}")
    print(f"faithful_partition_codes={code_checks}")
    print(f"sharp_initial_cutoff_witnesses={cutoff_checks}")
    print(f"profile_code_sha256={code_hash}")
    print(f"virtual_kernel_degrees=0..{MAXIMUM_KERNEL_DEGREE}")
    print("virtual_kernel_marked_cutoffs=0..5")
    print(f"virtual_kernel_coefficients={kernel_checks}")
    print(f"virtual_kernel_sha256={kernel_hash}")
    print(f"endpoint_jet_reconstructions={jet_checks}")
    print(f"partition_upset_mass_checks={upset_checks}")
    print(f"endpoint_jet_sha256={jet_hash}")
    print(f"central_class_reconstructions={class_checks}")
    print(f"symmetric_group_mass_checks={factorial_checks}")
    print(f"central_class_sha256={class_hash}")
    print("status=PASS")


if __name__ == "__main__":
    main()
