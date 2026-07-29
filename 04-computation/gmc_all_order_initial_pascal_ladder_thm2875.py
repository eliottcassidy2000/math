#!/usr/bin/env python3
"""Exact companion for THM-2875.

The theorem proves an initial, support-sensitive part of the experimental
all-order Pascal-response ladder

    n |-> L(d_n W^m) / L(d_n W),              m >= 2,

for nonzero nonnegative adjacent-difference cone elements W.  This companion
checks the exact product expansion, the rook/marked-cycle tensor, the
private-label insertion inequality, the fully polarized common-denominator
identity, the antitone-weight certificate in the proved range, and explicit
hostile boundaries.  It does not test or claim the still-open high-row
multi-ray regime.

Every truth-bearing gate raises explicitly, so normal and optimized Python
executions run identical checks.
"""

from functools import lru_cache
from fractions import Fraction
from itertools import combinations_with_replacement, permutations
from math import comb, factorial


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def elementary_symmetric(values):
    result = [1] + [0] * len(values)
    for value in values:
        for degree in range(len(values), 0, -1):
            result[degree] += value * result[degree - 1]
    return result


@lru_cache(maxsize=None)
def avoiders(blocks):
    """Cleared factorial tensor: forbidden-board/marked-cycle count."""
    blocks = tuple(blocks)
    total = sum(blocks)
    elementary = elementary_symmetric(blocks)
    return sum(
        (-1) ** degree * elementary[degree] * factorial(total - degree)
        for degree in range(len(blocks) + 1)
    )


def tensor(indices):
    blocks = tuple(index + 1 for index in indices)
    denominator = 1
    for block in blocks:
        denominator *= factorial(block)
    return Fraction(avoiders(blocks), denominator)


def h(row, column):
    return Fraction(comb(row + column, row))


def poly_add(left, right):
    result = [Fraction(0)] * max(len(left), len(right))
    for index in range(len(result)):
        result[index] = (
            (left[index] if index < len(left) else 0)
            + (right[index] if index < len(right) else 0)
        )
    return result


def poly_scale(scale, poly):
    return [scale * coefficient for coefficient in poly]


def poly_mul(left, right):
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for first_index, first in enumerate(left):
        for second_index, second in enumerate(right):
            result[first_index + second_index] += first * second
    return result


def f_poly(index):
    result = [Fraction(0)] * (index + 1)
    result[index] = Fraction(1, factorial(index))
    return result


def d_poly(index):
    return poly_add(f_poly(index + 1), poly_scale(-1, f_poly(index)))


def moment(poly):
    return sum(
        coefficient * factorial(index)
        for index, coefficient in enumerate(poly)
    )


def canonical_product(blocks):
    """Return the constant and d-basis coefficients in the exact product."""
    blocks = tuple(blocks)
    total = sum(blocks)
    elementary = elementary_symmetric(blocks)
    denominator = 1
    for block in blocks:
        denominator *= factorial(block)

    f_coefficients = {
        total - degree: Fraction(
            (-1) ** degree
            * elementary[degree]
            * factorial(total - degree),
            denominator,
        )
        for degree in range(len(blocks) + 1)
    }
    constant = sum(f_coefficients.values())
    d_coefficients = tuple(
        sum(
            coefficient
            for degree, coefficient in f_coefficients.items()
            if degree > index
        )
        for index in range(total)
    )
    return constant, d_coefficients


def literal_admissible_count(blocks):
    """Enumerate the marked-cycle predecessor model on small universes."""
    marks = []
    private_sets = []
    cursor = 0
    for block in blocks:
        mark = cursor
        marks.append(mark)
        private_sets.append(set(range(cursor + 1, cursor + block)))
        cursor += block

    count = 0
    for permutation in permutations(range(cursor)):
        inverse = [0] * cursor
        for source, target in enumerate(permutation):
            inverse[target] = source
        admitted = True
        for mark, private in zip(marks, private_sets):
            if permutation[mark] == mark or inverse[mark] in private:
                admitted = False
                break
        count += int(admitted)
    return count


def kernel(row, values):
    total = Fraction(0)
    for singled, column in enumerate(values):
        complement = values[:singled] + values[singled + 1 :]
        total += (
            tensor((row + 1,) + complement) * h(row, column)
            - tensor((row,) + complement) * h(row + 1, column)
        )
    return total


def cleared_kernel(a_value, blocks):
    total = 0
    for singled, block in enumerate(blocks):
        complement = blocks[:singled] + blocks[singled + 1 :]
        total += (
            avoiders((a_value + 1,) + complement)
            * avoiders((a_value, block))
            - avoiders((a_value,) + complement)
            * avoiders((a_value + 1, block))
        )
    return total


def kernel_weights(a_value, blocks):
    return tuple(
        block
        * factorial(a_value + block - 2)
        * avoiders((a_value,) + blocks[:index] + blocks[index + 1 :])
        for index, block in enumerate(blocks)
    )


def insertion_lower_bound(a_value, blocks):
    block_sum = sum(blocks)
    weights = kernel_weights(a_value, blocks)
    return sum(
        weight * (1 + a_value * block_sum - (2 * a_value + 1) * block)
        for weight, block in zip(weights, blocks)
    )


def tuple_gate(row, blocks):
    block_sum = sum(blocks)
    for first, first_block in enumerate(blocks):
        for second, second_block in enumerate(blocks):
            if first_block < second_block:
                remaining = block_sum - first_block - second_block
                if first_block * remaining < row:
                    return False
    return True


def derangement(index):
    if index == 0:
        return 1
    if index == 1:
        return 0
    previous_previous, previous = 1, 0
    for value in range(2, index + 1):
        current = (value - 1) * (previous + previous_previous)
        previous_previous, previous = previous, current
    return previous


def direct_response(row, poly):
    return moment(poly_mul(d_poly(row), poly))


def main():
    semiring_profiles = 0
    for arity in range(2, 7):
        for blocks in combinations_with_replacement(range(1, 5), arity):
            direct = [Fraction(1)]
            for block in blocks:
                direct = poly_mul(direct, d_poly(block - 1))
            constant, coefficients = canonical_product(blocks)
            require(constant == moment(direct), "canonical constant")
            require(all(coefficient >= 0 for coefficient in coefficients), "positive d coefficient")
            rebuilt = [constant]
            for index, coefficient in enumerate(coefficients):
                rebuilt = poly_add(rebuilt, poly_scale(coefficient, d_poly(index)))
            require(rebuilt == direct, "canonical product reconstruction")
            require(
                tensor(tuple(block - 1 for block in blocks)) == moment(direct),
                "rook tensor versus direct polynomial",
            )
            semiring_profiles += 1

    literal_rook_profiles = 0
    for arity in range(2, 5):
        for blocks in combinations_with_replacement(range(1, 7), arity):
            if sum(blocks) <= 8:
                require(
                    literal_admissible_count(blocks) == avoiders(blocks),
                    "literal marked-cycle count",
                )
                literal_rook_profiles += 1

    insertion_profiles = 0
    for arity in range(2, 7):
        for first in range(1, 7):
            for rest in combinations_with_replacement(range(1, 6), arity - 1):
                old_blocks = (first,) + rest
                old_total = sum(old_blocks)
                old_count = avoiders(old_blocks)
                new_count = avoiders((first + 1,) + rest)
                require(new_count >= old_total * old_count, "private-label insertion")
                if arity >= 3:
                    require(new_count > old_total * old_count, "strict shielding residue")
                else:
                    require(
                        (new_count == old_total * old_count) == (rest[0] == 1),
                        "two-block insertion equality",
                    )
                insertion_profiles += 1

    initial_kernel_cells = 0
    all_zero_residues = []
    finite_minima = []
    for order in range(2, 8):
        arity = order + 1
        zero_value = kernel(0, (0,) * arity)
        expected_zero = Fraction(
            arity
            * (
                (order - 1) * derangement(order + 1)
                + order * derangement(order)
            ),
            2,
        )
        require(zero_value == expected_zero, "all-zero derangement formula")
        all_zero_residues.append(int(zero_value))
        minimum = None
        minimum_values = None

        for row in range(order):
            for values in combinations_with_replacement(range(6), arity):
                blocks = tuple(value + 1 for value in values)
                a_value = row + 1
                common_denominator = factorial(a_value) * factorial(a_value + 1)
                for block in blocks:
                    common_denominator *= factorial(block)

                original = kernel(row, values)
                cleared = cleared_kernel(a_value, blocks)
                lower = insertion_lower_bound(a_value, blocks)
                require(original * common_denominator == cleared, "common-denominator kernel")
                require(cleared > lower > 0, "strict insertion/rearrangement lower bound")
                require(tuple_gate(row, blocks), "initial-range tuple gate")

                weights = kernel_weights(a_value, blocks)
                for first in range(arity):
                    for second in range(arity):
                        if blocks[first] < blocks[second]:
                            require(weights[first] > weights[second], "antitone kernel weights")
                require(
                    arity * sum(weight * block for weight, block in zip(weights, blocks))
                    <= sum(blocks) * sum(weights),
                    "weighted rearrangement",
                )
                if row == 0 and (minimum is None or original < minimum):
                    minimum = original
                    minimum_values = values
                initial_kernel_cells += 1

        require(minimum == zero_value, "finite all-zero minimum value")
        require(minimum_values == (0,) * arity, "finite all-zero minimum address")
        finite_minima.append(int(minimum))

    support_boundary_cells = 0
    for order in range(2, 7):
        arity = order + 1
        for blocks in combinations_with_replacement(range(1, 6), arity):
            row = (order - 1) * min(blocks) ** 2
            a_value = row + 1
            require(tuple_gate(row, blocks), "support-square tuple gate")
            weights = kernel_weights(a_value, blocks)
            for first in range(arity):
                for second in range(arity):
                    if blocks[first] < blocks[second]:
                        require(weights[first] > weights[second], "support-boundary weight order")
            require(insertion_lower_bound(a_value, blocks) > 0, "support-boundary lower bound")
            require(cleared_kernel(a_value, blocks) > 0, "support-boundary kernel")
            support_boundary_cells += 1

    # Equal-label tuples satisfy the polarized gate at every row.
    equal_ray_cells = 0
    for order in range(2, 8):
        for block in range(1, 8):
            for row in (0, order, 10, 25):
                blocks = (block,) * (order + 1)
                require(insertion_lower_bound(row + 1, blocks) > 0, "equal-ray all-row lower bound")
                require(cleared_kernel(row + 1, blocks) > 0, "equal-ray all-row kernel")
                equal_ray_cells += 1

    # The sufficient antitone mechanism genuinely fails outside its range.
    hostile_a = 4
    hostile_blocks = (1, 1, 2)
    hostile_weights = kernel_weights(hostile_a, hostile_blocks)
    require(hostile_weights == (8928, 8928, 9216), "outside-gate weight hostile")
    require(not tuple_gate(hostile_a - 1, hostile_blocks), "outside-gate hostile entered theorem")
    require(insertion_lower_bound(hostile_a, hostile_blocks) == 133632, "hostile lower value")
    require(cleared_kernel(hostile_a, hostile_blocks) == 172800, "hostile cleared value")
    require(kernel(hostile_a - 1, (0, 0, 1)) == 30, "hostile positive survivor")

    # The order-one boundary is exactly constant, not strict.
    require(kernel(0, (0, 0)) == 0, "order-one equality boundary")

    # Signed coefficients can reverse the ladder while both denominators stay positive.
    signed_w = poly_add(d_poly(0), poly_scale(Fraction(-14, 29), d_poly(1)))
    signed_cube = poly_mul(poly_mul(signed_w, signed_w), signed_w)
    signed_d0 = direct_response(0, signed_w)
    signed_d1 = direct_response(1, signed_w)
    signed_r0 = direct_response(0, signed_cube)
    signed_r1 = direct_response(1, signed_cube)
    signed_cross = signed_r1 * signed_d0 - signed_r0 * signed_d1
    require((signed_d0, signed_d1) == (Fraction(15, 29), Fraction(1, 29)), "signed denominators")
    require(signed_cross == Fraction(-5422944, 707281), "signed reversal")

    expected_residues = [6, 48, 420, 3840, 38010, 407904]
    require(all_zero_residues == expected_residues, "all-zero residue sequence")
    require(finite_minima == expected_residues, "bounded minimum controls")

    print("THM-2875 exact companion")
    print(f"semiring_profiles={semiring_profiles}")
    print(f"literal_rook_profiles={literal_rook_profiles}")
    print(f"insertion_profiles={insertion_profiles}")
    print(f"initial_kernel_cells={initial_kernel_cells}")
    print(f"support_boundary_cells={support_boundary_cells}")
    print(f"equal_ray_cells={equal_ray_cells}")
    print("all_zero_residues=" + ",".join(str(value) for value in all_zero_residues))
    print("finite_all_zero_minima=PASS(universe:m=2..7,n=0,labels=0..5)")
    print("outside_gate_weight_reversal=a4,b=(1,1,2),w=(8928,8928,9216),K=30")
    print("signed_hostile=W=d0-(14/29)d1,m=3,n=0,cross=-5422944/707281")
    print("high_row_multi_ray_regime=OPEN_NOT_TESTED_AS_THEOREM")
    print("status=PASS")


if __name__ == "__main__":
    main()
