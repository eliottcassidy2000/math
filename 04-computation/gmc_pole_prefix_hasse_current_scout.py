#!/usr/bin/env python3
"""Exact pole-prefix versus Young-refinement-current scout.

For one reduced row-response pole prefix R_j, form the virtual alphabets

    Phi^j(f) = sum_S c_S f[S-R_j],       Q^j = Q-R_j,

and the zero-mass normalized coefficient vector

    G^j_mu = Phi^j(h_N)m_mu[Q^j] - Phi^j(m_mu)h_N[Q^j].

The j=0 vector is the pointwise THM-3115 control.  This script asks whether
each prefixed vector is still a nonnegative fine-to-coarse partition-lattice
boundary.  It stops at the lexicographically first exact obstruction.
"""

import ast
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import chain
from math import factorial, prod
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
UPSTREAM = HERE / "gmc_product_gamma_arbitrary_anchored_schur_thm3110.py"
MAXIMUM_DEGREE = 9


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

# The signed bank has zero one-row marginal: every row vector occurs with
# total signed multiplicity zero.  Since power sums are additive over rows,
# this forces Phi^j(p_N)=0 for every common virtual prefix R_j.
for bank in BANKS:
    marginal = Counter()
    for coefficient, row in bank:
        for row_vector in row:
            marginal[row_vector] += coefficient
    require(not any(marginal.values()), "signed row marginal did not cancel")


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
    chain.from_iterable(partitions(degree)
                        for degree in range(1, MAXIMUM_DEGREE + 1)),
    key=lambda shape: (sum(shape), len(shape), shape),
))


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


def dominant_row(invariant):
    if invariant == 0:
        return ((0, 3), (2, 0), (2, 0), (1, 0))
    return ((0, 3), (1, 1), (2, 0), (1, 0), (1, 0))


def poly_multiply(left, right):
    answer = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, first in enumerate(left):
        for j, second in enumerate(right):
            answer[i + j] += first * second
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def factor_polynomial(roots):
    answer = [Fraction(1)]
    for root in roots:
        if root:
            answer = poly_multiply(answer, (Fraction(1), Fraction(-root)))
    return answer


def poly_add_scaled(target, source, scalar):
    if len(target) < len(source):
        target.extend([Fraction(0)] * (len(source) - len(target)))
    for index, value in enumerate(source):
        target[index] += scalar * value
    while len(target) > 1 and target[-1] == 0:
        target.pop()


def poly_value(poly, point):
    answer = Fraction(0)
    for coefficient in reversed(poly):
        answer = answer * point + coefficient
    return answer


def divide_linear(poly, root):
    """Divide exactly by 1-root*t."""

    require(root != 0 and len(poly) >= 2, "invalid linear division")
    quotient = [poly[0]]
    for index in range(1, len(poly) - 1):
        quotient.append(poly[index] + root * quotient[-1])
    require(poly[-1] == -root * quotient[-1],
            "claimed pole factor did not divide numerator")
    while len(quotient) > 1 and quotient[-1] == 0:
        quotient.pop()
    return quotient


def reduced_poles(invariant, bank, a_value, b_value):
    """Return the exact reduced common denominator pole multiset."""

    atoms = tuple(
        (coefficient, residual_roots(invariant, row, a_value, b_value))
        for coefficient, row in bank
    )
    common = Counter()
    for _, roots in atoms:
        counts = Counter(root for root in roots if root)
        for root, multiplicity in counts.items():
            common[root] = max(common[root], multiplicity)

    numerator = [Fraction(0)]
    for coefficient, roots in atoms:
        missing = common.copy()
        missing.subtract(Counter(root for root in roots if root))
        require(all(value >= 0 for value in missing.values()),
                "denominator LCM subtraction failed")
        quotient_roots = tuple(
            root for root, multiplicity in missing.items()
            for _ in range(multiplicity)
        )
        poly_add_scaled(numerator, factor_polynomial(quotient_roots),
                        Fraction(coefficient))

    for root in sorted(common):
        while common[root] and poly_value(numerator, Fraction(1, root)) == 0:
            numerator = divide_linear(numerator, root)
            common[root] -= 1

    leading_zeros = 0
    while leading_zeros < len(numerator) and numerator[leading_zeros] == 0:
        leading_zeros += 1
    require(leading_zeros == 5,
            f"forced numerator zero order drift: {leading_zeros}")
    poles = tuple(sorted(
        (root for root, multiplicity in common.items()
         for _ in range(multiplicity)), reverse=True))
    return poles, tuple(numerator[5:])


def all_monomial_values(roots, removed):
    power_sums = tuple(
        sum(root**degree for root in roots)
        - sum(root**degree for root in removed)
        for degree in range(1, MAXIMUM_DEGREE + 1)
    )
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


@lru_cache(maxsize=None)
def coarsens(fine, coarse):
    """Return whether coarse is obtained by merging parts of fine."""

    if sum(fine) != sum(coarse) or len(coarse) > len(fine):
        return False
    fine = tuple(sorted(fine, reverse=True))
    coarse = tuple(sorted(coarse, reverse=True))

    @lru_cache(maxsize=None)
    def recurse(targets, remaining):
        if not targets:
            return not remaining
        target = targets[0]
        remainders = set()

        def choose(start, subtotal, selected):
            if subtotal == target:
                remainders.add(tuple(sorted(
                    (remaining[index] for index in range(len(remaining))
                     if index not in selected), reverse=True)))
                return
            if subtotal > target:
                return
            previous = None
            for index in range(start, len(remaining)):
                if remaining[index] == previous:
                    continue
                previous = remaining[index]
                choose(index + 1, subtotal + remaining[index],
                       selected + (index,))

        choose(0, 0, ())
        return any(recurse(targets[1:], rest) for rest in remainders)

    return recurse(coarse, fine)


def transport_diagnostic(coefficients):
    """Return exact max flow and the first unresolved negative type."""

    degree = sum(next(iter(coefficients)))
    singleton = (1,) * degree
    positive = tuple(sorted(
        (shape, value) for shape, value in coefficients.items() if value > 0))
    negative = tuple(sorted(
        (shape, -value) for shape, value in coefficients.items()
        if value < 0 and shape != singleton))
    demand = sum((value for _, value in negative), Fraction(0))
    source = ("source",)
    sink = ("sink",)
    adjacency = defaultdict(list)
    capacity = {}

    def add_edge(left, right, value):
        adjacency[left].append(right)
        adjacency[right].append(left)
        capacity[left, right] = value
        capacity[right, left] = Fraction(0)

    for shape, value in positive:
        add_edge(source, ("positive", shape), value)
    for shape, value in negative:
        add_edge(("negative", shape), sink, value)
    infinite = demand or Fraction(1)
    for high, _ in positive:
        for low, _ in negative:
            if coarsens(low, high):
                add_edge(("positive", high), ("negative", low), infinite)

    flow = Fraction(0)
    while True:
        parent = {source: None}
        queue = [source]
        for vertex in queue:
            for neighbour in adjacency[vertex]:
                if neighbour not in parent and capacity[vertex, neighbour] > 0:
                    parent[neighbour] = vertex
                    queue.append(neighbour)
            if sink in parent:
                break
        if sink not in parent:
            break
        increment = infinite
        vertex = sink
        while parent[vertex] is not None:
            increment = min(increment, capacity[parent[vertex], vertex])
            vertex = parent[vertex]
        vertex = sink
        while parent[vertex] is not None:
            previous = parent[vertex]
            capacity[previous, vertex] -= increment
            capacity[vertex, previous] += increment
            vertex = previous
        flow += increment

    unresolved = tuple(
        (shape, capacity[("negative", shape), sink])
        for shape, _ in negative
        if capacity[("negative", shape), sink] > 0
    )
    singleton_ok = coefficients[singleton] <= 0
    passed = flow == demand and singleton_ok
    witness = (unresolved[0] if unresolved else
               ((singleton, coefficients[singleton]) if not singleton_ok
                else None))
    return passed, flow, demand, witness


def labelled_delete(vector):
    """THM-3119 same-label deletion on partition coordinates."""

    answer = defaultdict(Fraction)
    for shape, coefficient in vector.items():
        for part, multiplicity in Counter(shape).items():
            target = list(shape)
            target.remove(part)
            if part > 1:
                target.append(part - 1)
            target = tuple(sorted(target, reverse=True))
            answer[target] += coefficient * part * multiplicity
    return {shape: value for shape, value in answer.items() if value}


def rational_rank(matrix):
    work = [[Fraction(value) for value in row] for row in matrix]
    rank = 0
    for column in range(len(work[0]) if work else 0):
        pivot = next((row for row in range(rank, len(work))
                      if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        scale = work[rank][column]
        work[rank] = [value / scale for value in work[rank]]
        for row in range(len(work)):
            if row == rank or not work[row][column]:
                continue
            scale = work[row][column]
            work[row] = [left - scale * right
                         for left, right in zip(work[row], work[rank])]
        rank += 1
    return rank


def coefficient_vectors(invariant, bank, a_value, b_value, removed):
    atoms = tuple(
        (coefficient, all_monomial_values(
            residual_roots(invariant, row, a_value, b_value), removed))
        for coefficient, row in bank
    )
    quotient = all_monomial_values(
        residual_roots(invariant, dominant_row(invariant), a_value, b_value),
        removed,
    )
    answer = {}
    for degree in range(5, MAXIMUM_DEGREE + 1):
        shapes = partitions(degree)
        phi_monomial = {
            shape: sum(coefficient * values[shape]
                       for coefficient, values in atoms)
            for shape in shapes
        }
        phi_h = sum(phi_monomial.values())
        quotient_h = sum(quotient[shape] for shape in shapes)
        require(phi_monomial[(degree,)] == 0,
                "virtual-prefixed power-sum response did not vanish")
        vector = {
            shape: Fraction(
                phi_h * quotient[shape]
                - phi_monomial[shape] * quotient_h)
            for shape in shapes
        }
        require(sum(vector.values()) == 0,
                "prefixed coefficient mass did not cancel")
        answer[degree] = vector
    return answer


def coarsest_components(invariant, bank, a_value, b_value, removed, degree):
    """Return Phi(h), h(Q), Phi(p_N), and p_N(Q)."""

    atoms = tuple(
        (coefficient, all_monomial_values(
            residual_roots(invariant, row, a_value, b_value), removed))
        for coefficient, row in bank
    )
    quotient = all_monomial_values(
        residual_roots(invariant, dominant_row(invariant), a_value, b_value),
        removed,
    )
    shapes = partitions(degree)
    phi_h = sum(
        coefficient * sum(values[shape] for shape in shapes)
        for coefficient, values in atoms
    )
    quotient_h = sum(quotient[shape] for shape in shapes)
    phi_power = sum(coefficient * values[(degree,)]
                    for coefficient, values in atoms)
    quotient_power = quotient[(degree,)]
    return phi_h, quotient_h, phi_power, quotient_power


UNIVERSE = tuple(
    (a_value, b_value)
    for a_value in range(1, 11)
    for b_value in range(a_value + 1, min(3 * a_value + 4, 21) + 1)
)
require(len(UNIVERSE) == 115, "THM-3120 support universe drift")

checks = controls = active_prefix_cases = 0
first_active_failure = None
degree_range = [MAXIMUM_DEGREE + 1, -1]
for a_value, b_value in UNIVERSE:
    for invariant, bank in enumerate(BANKS):
        poles, numerator = reduced_poles(invariant, bank, a_value, b_value)
        numerator_degree = len(numerator) - 1
        require(numerator_degree <= len(poles),
                "flag prefix exceeds reduced pole list")
        degree_range[0] = min(degree_range[0], numerator_degree)
        degree_range[1] = max(degree_range[1], numerator_degree)
        for prefix_length in range(numerator_degree + 1):
            active_prefix_cases += 1
            removed = poles[:prefix_length]
            for degree, vector in coefficient_vectors(
                    invariant, bank, a_value, b_value, removed).items():
                passed, flow, demand, witness = transport_diagnostic(vector)
                checks += 1
                if prefix_length == 0:
                    controls += 1
                    require(passed, "THM-3115 pointwise control failed")
                if not passed:
                    record = (
                        degree, a_value, b_value, invariant,
                        prefix_length, removed, vector,
                        flow, demand, witness, poles, numerator,
                    )
                    record_key = (
                        degree, a_value, b_value, prefix_length,
                        len(witness[0]), witness[0], invariant,
                    )
                    if (first_active_failure is None
                            or record_key < first_active_failure[0]):
                        first_active_failure = (record_key, record)

require(active_prefix_cases == 8241,
        "THM-3120 active flag-prefix census drift")

# Reordering cannot repair the support-(1,3), I2 flag.  Every possible first
# physical reduced pole already destroys Hasse positivity by N=9.
reordering_records = []
for first_pole in range(1, 9):
    vectors = coefficient_vectors(1, BANKS[1], 1, 3, (first_pole,))
    first = None
    for degree in range(5, MAXIMUM_DEGREE + 1):
        passed, flow, demand, witness = transport_diagnostic(vectors[degree])
        if not passed:
            singleton_value = vectors[degree][(1,) * degree]
            require(singleton_value > 0 and demand - flow == singleton_value,
                    "first-pole upset/flow deficit mismatch")
            first = (degree, demand - flow)
            break
    require(first is not None, "pole reordering unexpectedly repaired the flag")
    reordering_records.append((first_pole,) + first)
require(tuple(record[1] for record in reordering_records)
        == (9, 9, 8, 8, 8, 8, 8, 8),
        "first-pole failure-degree census drift")
require(tuple(record[2] for record in reordering_records) == (
    230041249606656,
    86086680574464,
    164016803808,
    842945556672,
    4210741742976,
    12707727217920,
    24539565851424,
    33677263528128,
), "first-pole upset deficit census drift")

# The first prefix beyond the active quadratic flag at support (1,2) is a
# sharp hostile: algebraic commutation continues, but Hasse positivity fails.
boundary_invariant = 0
boundary_a, boundary_b, boundary_j, boundary_degree = 1, 2, 4, 5
boundary_poles, boundary_numerator = reduced_poles(
    boundary_invariant, BANKS[boundary_invariant], boundary_a, boundary_b)
boundary_removed = boundary_poles[:boundary_j]
boundary_vector = coefficient_vectors(
    boundary_invariant, BANKS[boundary_invariant],
    boundary_a, boundary_b, boundary_removed,
)[boundary_degree]
boundary_passed, boundary_flow, boundary_demand, boundary_witness = (
    transport_diagnostic(boundary_vector)
)
require(not boundary_passed, "declared beyond-flag hostile disappeared")

nonzero_boundary = ";".join(
    f"{shape}:{value}"
    for shape, value in sorted(
        boundary_vector.items(), key=lambda item: (len(item[0]), item[0]))
    if value
)
print("prefix_current=Phi^j(hN)m_mu(Q-Rj)-Phi^j(m_mu)hN(Q-Rj)")
print("prefix_semantics=virtual_subtraction_by_reduced_row_poles")
print("row_marginal_cancellation=2/2;Phi^j(p_N)=0")
print(f"support_universe={len(UNIVERSE)}:bank_cases={2*len(UNIVERSE)}")
print(f"active_P_degree_range={degree_range[0]}..{degree_range[1]}")
print(f"active_prefix_cases={active_prefix_cases}:N5..9_vectors={checks}")
print(f"controls_j0={controls}:all=PASS")
for first_pole, degree, deficit in reordering_records:
    print(
        f"first_pole={first_pole}:failure=N{degree}:"
        f"upset=P_N-minus-singleton:mass={-deficit}"
    )
if first_active_failure:
    first_failure = first_active_failure[1]
    (degree, a_value, b_value, invariant, prefix_length, removed, vector,
     flow, demand, witness, poles, numerator) = first_failure
    nonzero_active = ";".join(
        f"{shape}:{value}"
        for shape, value in sorted(
            vector.items(), key=lambda item: (len(item[0]), item[0]))
        if value
    )
    components = coarsest_components(
        invariant, BANKS[invariant], a_value, b_value, removed, degree)
    require(components[0] * components[3] - components[2] * components[1]
            == vector[(degree,)], "coarsest component identity drift")
    derangement_ghost = {shape: Fraction(0) for shape in partitions(degree)}
    for j in range(degree - 1):
        derangement_ghost[(degree - j,) + (1,) * j] = Fraction(
            (-1) ** j * factorial(degree), factorial(j))
    singleton = (1,) * degree
    derangement_ghost[singleton] = -sum(derangement_ghost.values())
    blind_upset = ((degree,), (degree - 1, 1))
    blind_current_mass = sum(vector[shape] for shape in blind_upset)
    blind_ghost_mass = sum(derangement_ghost[shape] for shape in blind_upset)
    require(blind_current_mass < 0 and blind_ghost_mass == 0,
            "derangement-ghost blind upset boundary drift")
    require(degree == 5, "first failure degree changed")
    ordered_shapes = (
        (5,), (4, 1), (3, 2), (3, 1, 1),
        (2, 2, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1),
    )
    kernel_1 = dict(zip(
        ordered_shapes,
        map(Fraction, (-1, 5, -2, -8, 6, 0, 0)),
    ))
    kernel_2 = dict(zip(
        ordered_shapes,
        map(Fraction, (1, -5, 0, 10, 0, -10, 4)),
    ))
    require(labelled_delete(kernel_1) == {}
            and labelled_delete(kernel_2) == {},
            "declared full labelled-deletion kernel basis failed")
    require(sum(kernel_1.values()) == sum(kernel_2.values()) == 0,
            "kernel basis lost zero mass")
    target_shapes = partitions(4)
    deletion_matrix = []
    for target in target_shapes:
        row = []
        for source_shape in ordered_shapes:
            image = labelled_delete({source_shape: Fraction(1)})
            row.append(image.get(target, Fraction(0)))
        deletion_matrix.append(row)
    require(rational_rank(deletion_matrix) == 5,
            "labelled deletion did not have the claimed two-dimensional kernel")
    weighted_kernel_1 = {
        shape: kernel_1[shape] * prod(factorial(part) for part in shape)
        for shape in ordered_shapes
    }
    weighted_kernel_2 = {
        shape: kernel_2[shape] * prod(factorial(part) for part in shape)
        for shape in ordered_shapes
    }
    weighted_blind_masses = tuple(
        sum(kernel[shape] for shape in blind_upset)
        for kernel in (weighted_kernel_1, weighted_kernel_2)
    )
    require(weighted_blind_masses == (0, 0),
            "factorial kernel changed the blind upset")
    print(
        f"active_first_failure=N{degree}:support=({a_value},{b_value}):"
        f"I{invariant + 1}:j={prefix_length}:Rj={removed}:"
        f"deficit={demand-flow}:witness={witness}"
    )
    print(f"active_reduced_poles={poles}:P_coefficients={numerator}")
    print(
        "active_coarsest_components="
        f"Phi_h:{components[0]},Q_h:{components[1]},"
        f"Phi_p:{components[2]},Q_p:{components[3]}"
    )
    print(f"active_coefficient_vector={nonzero_active}")
    print(
        f"derangement_ghost_blind_upset={blind_upset}:"
        f"Gmass={blind_current_mass}:O_N_mass={blind_ghost_mass}"
    )
    print(
        "full_factorial_kernel_blind_upset="
        f"WK1_mass:{weighted_blind_masses[0]},"
        f"WK2_mass:{weighted_blind_masses[1]}:no_positive_preimage"
    )
else:
    print("active_flag_hasse_transport=ALL_PASS")
print(
    "beyond_flag_hostile=N5:support=(1,2):I1:j=4:"
    f"Rj={boundary_removed}"
)
print(
    f"boundary_flow={boundary_flow}:demand={boundary_demand}:"
    f"deficit={boundary_demand-boundary_flow}:witness={boundary_witness}"
)
print(f"boundary_coefficient_vector={nonzero_boundary}")
print("scope=active_flag_test_finite_exact;all_prefix_preservation_false")
print("all_exact_checks=PASS")
