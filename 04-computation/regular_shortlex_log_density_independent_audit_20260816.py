#!/usr/bin/env python3
"""Independent hostile audit for THM-3499.

This companion does not import the submitted control program.  It rebuilds
finite-chain Cesaro projections from closed SCCs, exact stationary measures,
and exact absorption probabilities; exhausts every complete three-state
binary DFA and every complete two-state ternary DFA; and checks the q-adic
cylinder formula, martingale error interval, vector functional equation,
periodic/reducible/transient hostiles, alphabet reversal, and the four
submitted finite controls.
"""

from __future__ import annotations

import hashlib
import itertools
import math
from collections import Counter, deque
from fractions import Fraction


EXPECTED_EXACT_LEDGER_SHA256 = "046420bde716c6f83c7314c0010e1e4b19feffae59de7c9684fb4d0bdc4b1459"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def solve_linear(matrix: list[list[Fraction]], rhs: list[Fraction]) -> list[Fraction]:
    """Solve a nonsingular square rational system by full elimination."""

    size = len(rhs)
    if size == 0:
        return []
    augmented = [list(matrix[row]) + [rhs[row]] for row in range(size)]
    for column in range(size):
        pivot = next((row for row in range(column, size) if augmented[row][column]), None)
        require(pivot is not None, "singular exact system in finite-chain audit")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [entry / scale for entry in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            scale = augmented[row][column]
            if scale:
                augmented[row] = [
                    augmented[row][index] - scale * augmented[column][index]
                    for index in range(size + 1)
                ]
    return [augmented[row][-1] for row in range(size)]


def reachable_states(transitions: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    reached = {0}
    queue = deque([0])
    while queue:
        state = queue.popleft()
        for target in transitions[state]:
            if target not in reached:
                reached.add(target)
                queue.append(target)
    return tuple(sorted(reached))


def strongly_connected_components(
    transitions: tuple[tuple[int, ...], ...], reachable: tuple[int, ...]
) -> tuple[tuple[int, ...], ...]:
    """Kosaraju SCCs on the reachable directed graph."""

    allowed = set(reachable)
    seen: set[int] = set()
    order: list[int] = []

    def forward(state: int) -> None:
        seen.add(state)
        for target in transitions[state]:
            if target in allowed and target not in seen:
                forward(target)
        order.append(state)

    for state in reachable:
        if state not in seen:
            forward(state)

    reverse = {state: set() for state in reachable}
    for state in reachable:
        for target in transitions[state]:
            if target in allowed:
                reverse[target].add(state)

    seen.clear()
    components: list[tuple[int, ...]] = []

    def backward(state: int, component: list[int]) -> None:
        seen.add(state)
        component.append(state)
        for source in reverse[state]:
            if source not in seen:
                backward(source, component)

    for state in reversed(order):
        if state not in seen:
            component: list[int] = []
            backward(state, component)
            components.append(tuple(sorted(component)))
    return tuple(sorted(components))


def class_period(
    transitions: tuple[tuple[int, ...], ...], component: tuple[int, ...]
) -> int:
    """Return the period (gcd of directed cycle lengths) of one closed SCC."""

    allowed = set(component)
    root = component[0]
    distance = {root: 0}
    queue = deque([root])
    while queue:
        state = queue.popleft()
        for target in transitions[state]:
            if target in allowed and target not in distance:
                distance[target] = distance[state] + 1
                queue.append(target)
    period = 0
    for state in component:
        for target in transitions[state]:
            if target in allowed:
                period = math.gcd(period, distance[state] + 1 - distance[target])
    return abs(period)


def stationary_distribution(
    transitions: tuple[tuple[int, ...], ...], q: int, component: tuple[int, ...]
) -> tuple[Fraction, ...]:
    """Exact stationary row vector on one closed irreducible component."""

    size = len(component)
    location = {state: index for index, state in enumerate(component)}
    matrix: list[list[Fraction]] = []
    rhs: list[Fraction] = []
    for target_index in range(size - 1):
        target = component[target_index]
        row = []
        for source_index, source in enumerate(component):
            count = transitions[source].count(target)
            entry = Fraction(count, q)
            if source_index == target_index:
                entry -= 1
            row.append(entry)
        matrix.append(row)
        rhs.append(Fraction(0))
    matrix.append([Fraction(1) for _ in component])
    rhs.append(Fraction(1))
    answer = tuple(solve_linear(matrix, rhs))
    require(sum(answer) == 1 and all(value > 0 for value in answer), "bad stationary law")
    for target in component:
        target_index = location[target]
        incoming = sum(
            answer[source_index] * Fraction(transitions[source].count(target), q)
            for source_index, source in enumerate(component)
        )
        require(incoming == answer[target_index], "stationarity failed")
    return answer


class ChainData:
    def __init__(self, transitions: tuple[tuple[int, ...], ...], q: int):
        self.transitions = transitions
        self.q = q
        self.reachable = reachable_states(transitions)
        components = strongly_connected_components(transitions, self.reachable)
        self.closed = tuple(
            component
            for component in components
            if all(target in component for state in component for target in transitions[state])
        )
        require(self.closed, "finite reachable chain has no closed class")
        self.stationary = tuple(
            stationary_distribution(transitions, q, component) for component in self.closed
        )
        closed_states = {state for component in self.closed for state in component}
        self.transient = tuple(state for state in self.reachable if state not in closed_states)
        self.absorption = self._absorption_probabilities()
        self.periods = tuple(class_period(transitions, component) for component in self.closed)

    def _absorption_probabilities(self) -> tuple[tuple[Fraction, ...], ...]:
        class_of = {
            state: class_index
            for class_index, component in enumerate(self.closed)
            for state in component
        }
        transient_location = {state: index for index, state in enumerate(self.transient)}
        answers = [[Fraction(0) for _ in self.closed] for _ in self.transitions]
        for state, class_index in class_of.items():
            answers[state][class_index] = Fraction(1)
        if self.transient:
            matrix: list[list[Fraction]] = []
            for state in self.transient:
                row = [Fraction(int(state == target)) for target in self.transient]
                for target in self.transitions[state]:
                    if target in transient_location:
                        row[transient_location[target]] -= Fraction(1, self.q)
                matrix.append(row)
            for class_index, component in enumerate(self.closed):
                component_set = set(component)
                rhs = [
                    Fraction(sum(target in component_set for target in self.transitions[state]), self.q)
                    for state in self.transient
                ]
                values = solve_linear(matrix, rhs)
                for state, value in zip(self.transient, values):
                    answers[state][class_index] = value
        for state in self.reachable:
            require(sum(answers[state]) == 1, "absorption probabilities do not sum to one")
            require(all(0 <= value <= 1 for value in answers[state]), "bad absorption probability")
        return tuple(tuple(row) for row in answers)

    def harmonic_data(
        self, accepting: tuple[int, ...]
    ) -> tuple[tuple[Fraction, ...], tuple[Fraction, ...], tuple[Fraction, ...]]:
        means = []
        for component, stationary in zip(self.closed, self.stationary):
            means.append(
                sum(
                    stationary[index] * accepting[state]
                    for index, state in enumerate(component)
                )
            )
        h = [Fraction(0) for _ in self.transitions]
        deviation = [Fraction(0) for _ in self.transitions]
        for state in self.reachable:
            h[state] = sum(
                self.absorption[state][index] * means[index] for index in range(len(means))
            )
            deviation[state] = sum(
                self.absorption[state][index] * abs(means[index] - h[state])
                for index in range(len(means))
            )
        for state in self.reachable:
            next_mean = sum(h[target] for target in self.transitions[state]) / self.q
            require(
                next_mean == h[state],
                f"Ph=h failed exactly at state {state}: h={h[state]}, Ph={next_mean}, "
                f"accepting={accepting}, transitions={self.transitions}",
            )
        if len(self.closed) == 1:
            require(all(h[state] == means[0] for state in self.reachable), "one-class h is not constant")
            require(all(deviation[state] == 0 for state in self.reachable), "one-class martingale moved")
        return tuple(means), tuple(h), tuple(deviation)


def transition_tuple(flat: tuple[int, ...], states: int, q: int) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(flat[q * state : q * (state + 1)]) for state in range(states))


def word_end_states(
    transitions: tuple[tuple[int, ...], ...], q: int, depth: int, start: int = 0
) -> list[int]:
    states = [start]
    for _ in range(depth):
        states = [transitions[state][digit] for state in states for digit in range(q)]
    return states


def cylinder_state_weights(
    transitions: tuple[tuple[int, ...], ...], q: int, depth: int, t: float, start: int = 0
) -> tuple[float, ...]:
    states = word_end_states(transitions, q, depth, start)
    scale = q**depth
    weights = [0.0 for _ in transitions]
    for rank, state in enumerate(states):
        # Integral over the rank cylinder of dx/(t+x).
        weights[state] += math.log1p(1.0 / (t * scale + rank))
    return tuple(weights)


def state_distribution(
    transitions: tuple[tuple[int, ...], ...], q: int, depth: int
) -> tuple[Fraction, ...]:
    counts = [0 for _ in transitions]
    for state in word_end_states(transitions, q, depth):
        counts[state] += 1
    denominator = q**depth
    return tuple(Fraction(count, denominator) for count in counts)


def cylinder_coefficient(h: tuple[Fraction, ...], weights: tuple[float, ...], q: int) -> float:
    return sum(float(h[state]) * weights[state] for state in range(len(h))) / math.log(q)


def martingale_bound(
    deviation: tuple[Fraction, ...], distribution: tuple[Fraction, ...], q: int
) -> float:
    expectation = sum(distribution[state] * deviation[state] for state in range(len(deviation)))
    return (q - 1) * float(expectation) / math.log(q)


def cesaro_projection_error(
    transitions: tuple[tuple[int, ...], ...], q: int, accepting: tuple[int, ...],
    h: tuple[Fraction, ...], steps: int = 512,
) -> float:
    values = [float(value) for value in accepting]
    sums = [0.0 for _ in transitions]
    for _ in range(steps):
        for state in range(len(values)):
            sums[state] += values[state]
        values = [
            sum(values[target] for target in transitions[state]) / q
            for state in range(len(values))
        ]
    return max(
        abs(sums[state] / steps - float(h[state])) for state in reachable_states(transitions)
    )


def finite_functional_equation(
    chain: ChainData, h: tuple[Fraction, ...], depth: int, t: float
) -> float:
    left_errors = []
    for state in chain.reachable:
        left_weights = cylinder_state_weights(chain.transitions, chain.q, depth + 1, t, state)
        left = sum(float(h[index]) * left_weights[index] for index in range(len(h)))
        right = 0.0
        for digit, target in enumerate(chain.transitions[state]):
            right_weights = cylinder_state_weights(
                chain.transitions, chain.q, depth, chain.q * t + digit, target
            )
            right += sum(float(h[index]) * right_weights[index] for index in range(len(h)))
        left_errors.append(abs(left - right))
    return max(left_errors, default=0.0)


def endpoint_density(
    transitions: tuple[tuple[int, ...], ...], q: int, accepting: tuple[int, ...],
    final_depth: int, partial_fraction: Fraction | None = None,
) -> float:
    """Direct harmonic sum through a complete or partial final shortlex level."""

    harmonic_mass = 0.0
    first_index = 1
    last_index = 1
    states = [0]
    for depth in range(final_depth + 1):
        count = len(states)
        take = count
        if depth == final_depth and partial_fraction is not None:
            take = int(count * partial_fraction)
        harmonic_mass += math.fsum(
            1.0 / (first_index + rank)
            for rank, state in enumerate(states[:take])
            if accepting[state]
        )
        last_index = first_index + take - 1
        if depth == final_depth:
            break
        first_index += count
        states = [transitions[state][digit] for state in states for digit in range(q)]
    return harmonic_mass / math.log(last_index)


def trim_has_full_component(
    transitions: tuple[tuple[int, ...], ...], accepting: tuple[int, ...]
) -> bool:
    """Graph form of the trimmed transition-count boundary rho=q."""

    reachable = set(reachable_states(transitions))
    reverse = {state: set() for state in range(len(transitions))}
    for state, row in enumerate(transitions):
        for target in row:
            reverse[target].add(state)
    coaccessible = {state for state in reachable if accepting[state]}
    queue = deque(coaccessible)
    while queue:
        state = queue.popleft()
        for source in reverse[state]:
            if source in reachable and source not in coaccessible:
                coaccessible.add(source)
                queue.append(source)
    if not coaccessible:
        return False
    components = strongly_connected_components(transitions, tuple(sorted(coaccessible)))
    return any(
        all(target in component for state in component for target in transitions[state])
        for component in components
    )


def cycle_type(permutation: tuple[int, ...]) -> tuple[int, ...]:
    seen: set[int] = set()
    lengths = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        state = start
        length = 0
        while state not in seen:
            seen.add(state)
            state = permutation[state]
            length += 1
        lengths.append(length)
    return tuple(sorted(lengths))


def permutation_after(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[right[index]] for index in range(len(left)))


def cayley_automaton(identity, generators, multiply):
    states = [identity]
    location = {identity: 0}
    queue = deque([identity])
    while queue:
        state = queue.popleft()
        for generator in generators:
            target = multiply(generator, state)
            if target not in location:
                location[target] = len(states)
                states.append(target)
                queue.append(target)
    transitions = tuple(
        tuple(location[multiply(generator, state)] for generator in generators) for state in states
    )
    return tuple(states), transitions


def exhaustive_small_dfas(states: int, q: int, depth_small: int, depth_large: int):
    table_count = states ** (states * q)
    language_count = table_count * 2**states
    stats = Counter()
    worst_cesaro = (0.0, None)
    worst_bound = (0.0, None)
    strongest_unweighted = (-1.0, None)
    strongest_reversal = (-1.0, None)
    ledger_rows = []

    for flat in itertools.product(range(states), repeat=states * q):
        transitions = transition_tuple(flat, states, q)
        chain = ChainData(transitions, q)
        stats["tables"] += 1
        stats["transient_tables"] += bool(chain.transient)
        stats["multi_closed_tables"] += len(chain.closed) > 1
        stats["periodic_tables"] += any(period > 1 for period in chain.periods)
        stats[f"closed_{len(chain.closed)}"] += 1
        weights_small = cylinder_state_weights(transitions, q, depth_small, 1 / (q - 1))
        weights_large = cylinder_state_weights(transitions, q, depth_large, 1 / (q - 1))
        distribution_small = state_distribution(transitions, q, depth_small)
        distribution_large = state_distribution(transitions, q, depth_large)
        level_profiles = []
        level_states = [0]
        first_index = 1
        for depth in range(depth_small + 1):
            counts = [0 for _ in transitions]
            harmonic_weights = [0.0 for _ in transitions]
            for rank, state in enumerate(level_states):
                counts[state] += 1
                harmonic_weights[state] += 1 / (first_index + rank)
            level_profiles.append((tuple(counts), tuple(harmonic_weights)))
            first_index += len(level_states)
            level_states = [
                transitions[state][digit] for state in level_states for digit in range(q)
            ]

        reversed_transitions = tuple(tuple(reversed(row)) for row in transitions)
        reversed_chain = ChainData(reversed_transitions, q)
        reversed_weights = cylinder_state_weights(
            reversed_transitions, q, depth_large, 1 / (q - 1)
        )
        reversed_distribution = state_distribution(reversed_transitions, q, depth_large)

        for mask in range(2**states):
            accepting = tuple((mask >> state) & 1 for state in range(states))
            means, h, deviation = chain.harmonic_data(accepting)
            accepting_closed = any(
                any(accepting[state] for state in component) for component in chain.closed
            )
            require(
                trim_has_full_component(transitions, accepting) == accepting_closed,
                "trimmed rho=q graph boundary disagrees with accepting recurrent class",
            )
            stats["divergent_languages"] += accepting_closed
            stats["convergent_languages"] += not accepting_closed
            for depth, (counts, harmonic_weights) in enumerate(level_profiles):
                accepted_count = sum(counts[state] * accepting[state] for state in range(states))
                level_mass = sum(
                    harmonic_weights[state] * accepting[state] for state in range(states)
                )
                lower = Fraction((q - 1) * accepted_count, q ** (depth + 1) - 1)
                upper = Fraction((q - 1) * accepted_count, q**depth + q - 2)
                require(float(lower) <= level_mass + 2e-15, "Kraft lower constant failed")
                require(level_mass <= float(upper) + 2e-15, "Kraft upper constant failed")
                if depth == 0:
                    require(level_mass == accepted_count, "n=0 Kraft boundary is not exact")
            coefficient_small = cylinder_coefficient(h, weights_small, q)
            coefficient_large = cylinder_coefficient(h, weights_large, q)
            bound_small = martingale_bound(deviation, distribution_small, q)
            bound_large = martingale_bound(deviation, distribution_large, q)
            require(
                abs(coefficient_small - coefficient_large) <= bound_small + bound_large + 2e-13,
                "nested cylinder approximants escaped certified martingale intervals",
            )
            require(-bound_large - 1e-13 <= coefficient_large <= 1 + bound_large + 1e-13,
                    "coefficient interval left [0,1]")
            if len(chain.closed) == 1:
                require(bound_large == 0, "single recurrent class has address uncertainty")
                require(abs(coefficient_large - float(means[0])) < 3e-13,
                        "single-class coefficient is not stationary mass")

            cesaro_error = cesaro_projection_error(transitions, q, accepting, h)
            if cesaro_error > worst_cesaro[0]:
                worst_cesaro = (cesaro_error, (flat, mask))
            if bound_large > worst_bound[0]:
                worst_bound = (bound_large, (flat, mask))

            certified_unweighted_gap = abs(coefficient_large - float(h[0])) - bound_large
            if certified_unweighted_gap > strongest_unweighted[0]:
                strongest_unweighted = (certified_unweighted_gap, (flat, mask, coefficient_large, h[0], bound_large))

            reverse_means, reverse_h, reverse_deviation = reversed_chain.harmonic_data(accepting)
            reverse_coefficient = cylinder_coefficient(reverse_h, reversed_weights, q)
            reverse_bound = martingale_bound(reverse_deviation, reversed_distribution, q)
            certified_reversal_gap = abs(coefficient_large - reverse_coefficient) - bound_large - reverse_bound
            if certified_reversal_gap > strongest_reversal[0]:
                strongest_reversal = (
                    certified_reversal_gap,
                    (flat, mask, coefficient_large, reverse_coefficient, bound_large, reverse_bound),
                )

            complement = tuple(1 - bit for bit in accepting)
            _, complement_h, _ = chain.harmonic_data(complement)
            complement_coefficient = cylinder_coefficient(complement_h, weights_large, q)
            require(abs(coefficient_large + complement_coefficient - 1) < 4e-13,
                    "language/complement coefficients failed to sum to one")

            ledger_rows.append(
                "|".join(
                    [
                        f"q={q}",
                        f"n={states}",
                        f"T={','.join(map(str, flat))}",
                        f"A={mask}",
                        "C=" + ";".join(
                            ",".join(map(str, component)) for component in chain.closed
                        ),
                        "p=" + ",".join(map(str, chain.periods)),
                        "m=" + ",".join(map(str, means)),
                        "h=" + ",".join(map(str, (h[state] for state in chain.reachable))),
                    ]
                )
            )

    require(stats["tables"] == table_count, "transition-table census changed")
    require(len(ledger_rows) == language_count, "language census changed")
    require(worst_cesaro[0] < 0.012, "finite Cesaro audit is unexpectedly far from exact projection")
    return stats, worst_cesaro, worst_bound, strongest_unweighted, strongest_reversal, ledger_rows


def main() -> None:
    # Shortlex normalization and explicit partial-level uniform bound.
    for q in range(2, 8):
        previous_count = 0
        for depth in range(9):
            first = (q**depth - 1) // (q - 1) + 1
            require(first == previous_count + 1, "shortlex first-index normalization failed")
            c_n = Fraction(first, q**depth)
            expected = Fraction(1, q - 1) + Fraction(q - 2, (q - 1) * q**depth)
            require(c_n == expected, "scaled shortlex offset failed")
            full_mass = math.fsum(1 / index for index in range(first, first + q**depth))
            require(full_mass <= 1 + math.log(q), "partial-level O(1) bound failed")
            previous_count += q**depth

    binary = exhaustive_small_dfas(states=3, q=2, depth_small=7, depth_large=12)
    ternary = exhaustive_small_dfas(states=2, q=3, depth_small=5, depth_large=8)
    exact_rows = binary[-1] + ternary[-1]

    # Submitted control 1: the S3 x C2 variable Berggren language.
    identity3 = (0, 1, 2)
    rotation3 = (1, 2, 0)
    reflection3 = (1, 0, 2)
    variable_generators = ((rotation3, 1), (reflection3, 1), (rotation3, 0))

    def variable_multiply(left, right):
        return permutation_after(left[0], right[0]), left[1] ^ right[1]

    variable_states, variable_transitions = cayley_automaton(
        (identity3, 0), variable_generators, variable_multiply
    )
    variable_accepting = tuple(
        int(
            (parity == 0 and cycle_type(permutation) == (1, 1, 1))
            or (parity == 1 and cycle_type(permutation) == (1, 2))
        )
        for permutation, parity in variable_states
    )
    variable_chain = ChainData(variable_transitions, 3)
    variable_means, variable_h, _ = variable_chain.harmonic_data(variable_accepting)
    require(len(variable_states) == 12 and sum(variable_accepting) == 4, "variable control census")
    require(len(variable_chain.closed) == 1 and variable_means == (Fraction(1, 3),),
            "variable control stationary mass")

    # Submitted control 2: the S4 x D4 fixed-drift language.
    identity4 = (0, 1, 2, 3)
    A = (1, 3, 2, 0)
    B = (1, 3, 0, 2)
    C = (3, 1, 0, 2)
    G = B
    T = (1, 0, 3, 2)
    fixed_generators = ((A, G), (B, G), (C, T))

    def fixed_multiply(left, right):
        return permutation_after(left[0], right[0]), permutation_after(left[1], right[1])

    fixed_states, fixed_transitions = cayley_automaton(
        (identity4, identity4), fixed_generators, fixed_multiply
    )
    fixed_accepting = tuple(
        int(cycle_type(first) == cycle_type(second)) for first, second in fixed_states
    )
    fixed_chain = ChainData(fixed_transitions, 3)
    fixed_means, fixed_h, _ = fixed_chain.harmonic_data(fixed_accepting)
    require(len(fixed_states) == 192 and sum(fixed_accepting) == 34, "fixed control census")
    require(fixed_chain.periods == (2,) and fixed_means == (Fraction(17, 96),),
            "fixed control period/stationary mass")
    distances = {0: 0}
    queue = deque([0])
    while queue:
        state = queue.popleft()
        for target in fixed_transitions[state]:
            if target not in distances:
                distances[target] = distances[state] + 1
                queue.append(target)
    require(
        all((distances[target] - distances[state]) % 2 == 1
            for state in range(len(fixed_states)) for target in fixed_transitions[state]),
        "fixed shortest-distance parity is not a lawful period character",
    )
    fixed_parity_accept = tuple(
        sum(fixed_accepting[state] and distances[state] % 2 == parity for state in range(len(fixed_states)))
        for parity in (0, 1)
    )
    require(fixed_parity_accept == (16, 18), "fixed parity chamber census")

    # Submitted controls 3 and 4: periodic even length and two prefix basins.
    even_transitions = ((1, 1), (0, 0))
    even_accepting = (1, 0)
    even_chain = ChainData(even_transitions, 2)
    even_means, even_h, _ = even_chain.harmonic_data(even_accepting)
    require(even_chain.periods == (2,) and even_means == (Fraction(1, 2),), "even-length control")

    prefix_zero_transitions = ((1, 2), (1, 1), (2, 2))
    prefix_one_transitions = ((2, 1), (1, 1), (2, 2))
    prefix_accepting = (0, 1, 0)
    prefix_zero_chain = ChainData(prefix_zero_transitions, 2)
    prefix_one_chain = ChainData(prefix_one_transitions, 2)
    _, prefix_zero_h, prefix_zero_dev = prefix_zero_chain.harmonic_data(prefix_accepting)
    _, prefix_one_h, prefix_one_dev = prefix_one_chain.harmonic_data(prefix_accepting)
    prefix_depth = 14
    prefix_zero_value = cylinder_coefficient(
        prefix_zero_h, cylinder_state_weights(prefix_zero_transitions, 2, prefix_depth, 1.0), 2
    )
    prefix_one_value = cylinder_coefficient(
        prefix_one_h, cylinder_state_weights(prefix_one_transitions, 2, prefix_depth, 1.0), 2
    )
    require(martingale_bound(prefix_zero_dev, state_distribution(prefix_zero_transitions, 2, prefix_depth), 2) == 0,
            "prefix basin should be decided after the first digit")
    require(abs(prefix_zero_value - math.log(3 / 2) / math.log(2)) < 3e-13,
            "prefix-zero coefficient")
    require(abs(prefix_one_value - math.log(4 / 3) / math.log(2)) < 3e-13,
            "prefix-one coefficient")
    require(abs(prefix_zero_value + prefix_one_value - 1) < 3e-13, "prefix complement")

    # Delayed reducible hostile: two transient states leave uncertainty at
    # every finite prefix, so the martingale certificate radius is positive
    # and decays rather than vanishing after the first digit.
    delayed_split_transitions = ((0, 1), (2, 3), (2, 2), (3, 3))
    delayed_split_accepting = (0, 0, 1, 0)
    delayed_chain = ChainData(delayed_split_transitions, 2)
    delayed_means, delayed_h, delayed_dev = delayed_chain.harmonic_data(delayed_split_accepting)
    delayed_bounds = []
    delayed_values = []
    for depth in (4, 12):
        delayed_values.append(
            cylinder_coefficient(
                delayed_h, cylinder_state_weights(delayed_split_transitions, 2, depth, 1.0), 2
            )
        )
        delayed_bounds.append(
            martingale_bound(
                delayed_dev, state_distribution(delayed_split_transitions, 2, depth), 2
            )
        )
    require(delayed_bounds[0] > delayed_bounds[1] > 0, "delayed martingale radius did not shrink")
    require(
        abs(delayed_values[0] - delayed_values[1]) <= sum(delayed_bounds) + 2e-13,
        "delayed cylinder approximants escaped their certified intervals",
    )

    # A reducible transient split into two period-two recurrent classes.
    periodic_split_transitions = (
        (1, 3),
        (2, 2),
        (1, 1),
        (4, 4),
        (3, 3),
    )
    periodic_split_accepting = (0, 1, 0, 1, 1)
    periodic_split_chain = ChainData(periodic_split_transitions, 2)
    split_means, split_h, split_dev = periodic_split_chain.harmonic_data(periodic_split_accepting)
    require(periodic_split_chain.periods == (2, 2), "split recurrent periods")
    require(sorted(split_means) == [Fraction(1, 2), Fraction(1)], "split class means")
    split_weights = cylinder_state_weights(periodic_split_transitions, 2, 12, 1.0)
    split_value = cylinder_coefficient(split_h, split_weights, 2)
    split_bound = martingale_bound(
        split_dev, state_distribution(periodic_split_transitions, 2, 12), 2
    )
    require(split_bound == 0, "first digit should decide the periodic split")

    # The finite-cylinder vector functional equation is exact at every depth.
    functional_errors = [
        finite_functional_equation(prefix_zero_chain, prefix_zero_h, 7, Fraction(7, 5)),
        finite_functional_equation(delayed_chain, delayed_h, 7, Fraction(9, 8)),
        finite_functional_equation(periodic_split_chain, split_h, 7, Fraction(11, 7)),
        finite_functional_equation(variable_chain, variable_h, 5, Fraction(5, 4)),
    ]
    require(max(functional_errors) < 2e-13, "finite vector functional equation failed")

    # Direct arbitrary-cutoff hostiles.  Mid-level and endpoint values must
    # approach the same coefficient; no complete-level-only theorem is used.
    direct_even_endpoint = endpoint_density(even_transitions, 2, even_accepting, 18)
    direct_even_midpoint = endpoint_density(
        even_transitions, 2, even_accepting, 18, Fraction(1, 2)
    )
    direct_prefix_endpoint = endpoint_density(prefix_zero_transitions, 2, prefix_accepting, 18)
    direct_prefix_midpoint = endpoint_density(
        prefix_zero_transitions, 2, prefix_accepting, 18, Fraction(1, 2)
    )
    require(abs(direct_even_endpoint - 0.5) < 0.06, "even endpoint convergence hostile")
    require(abs(direct_even_midpoint - 0.5) < 0.06, "even partial-level convergence hostile")
    require(abs(direct_prefix_endpoint - prefix_zero_value) < 0.04, "prefix endpoint convergence")
    require(abs(direct_prefix_midpoint - prefix_zero_value) < 0.04, "prefix partial convergence")

    exact_rows.extend(
        [
            f"variable:{len(variable_states)}:{sum(variable_accepting)}:{variable_means[0]}:{variable_chain.periods}",
            f"fixed:{len(fixed_states)}:{sum(fixed_accepting)}:{fixed_means[0]}:{fixed_chain.periods}:{fixed_parity_accept}",
            f"even:{even_means[0]}:{even_chain.periods}",
            f"prefix:{prefix_zero_h}:{prefix_one_h}",
            f"delayed:{delayed_means}:{delayed_h}",
            f"split:{split_means}:{periodic_split_chain.periods}",
        ]
    )
    exact_digest = hashlib.sha256("\n".join(exact_rows).encode("ascii")).hexdigest()
    require(
        exact_digest == EXPECTED_EXACT_LEDGER_SHA256,
        f"independent exact ledger changed: {exact_digest}",
    )

    print("== THM-3499 independent hostile audit ==")
    print("shortlex normalization: exact for q=2..7, depths 0..8")
    print("partial final level: explicit uniform bound <= 1+log(q)")
    for label, result in (("binary-3-state", binary), ("ternary-2-state", ternary)):
        stats, worst_cesaro, worst_bound, unweighted, reversal, _ = result
        print(
            f"{label}: tables={stats['tables']} languages="
            f"{stats['tables'] * (8 if label.startswith('binary') else 4)} "
            f"transient={stats['transient_tables']} multi-closed={stats['multi_closed_tables']} "
            f"periodic={stats['periodic_tables']} closed-class-hist="
            f"{dict(sorted((key, value) for key, value in stats.items() if key.startswith('closed_')))}"
        )
        print(
            f"  Kraft/spectral dichotomy: convergent={stats['convergent_languages']} "
            f"positive-log={stats['divergent_languages']}"
        )
        print(
            f"  worst Cesaro(512) error={worst_cesaro[0]:.12g}; "
            f"largest depth certificate radius={worst_bound[0]:.12g}"
        )
        print(
            f"  strongest certified unweighted-basin failure={unweighted[0]:.12g}; "
            f"strongest certified alphabet-reversal change={reversal[0]:.12g}"
        )
    print("variable Berggren control: 12 states, 4 accept, coefficient=1/3, one irreducible class")
    print("fixed Berggren control: 192 states, 34 accept, period=2, chambers=16/18, coefficient=17/96")
    print("even-length control: period=2, natural endpoint oscillation survives, logarithmic coefficient=1/2")
    print(
        "prefix controls: zero=log(3/2)/log(2), one=log(4/3)/log(2), "
        f"values={prefix_zero_value:.15g}/{prefix_one_value:.15g}"
    )
    print(
        "reducible periodic split: class means=1/2,1; address-weighted coefficient="
        f"{split_value:.15g}"
    )
    print(
        "delayed two-basin transient: depth-4/12 radii="
        f"{delayed_bounds[0]:.12g}/{delayed_bounds[1]:.12g}; "
        f"approximants={delayed_values[0]:.12g}/{delayed_values[1]:.12g}"
    )
    print(f"finite functional-equation max error={max(functional_errors):.12g}")
    print(
        "direct depth-18 arbitrary-cutoff controls: "
        f"even(endpoint/mid)={direct_even_endpoint:.12g}/{direct_even_midpoint:.12g}; "
        f"prefix(endpoint/mid)={direct_prefix_endpoint:.12g}/{direct_prefix_midpoint:.12g}"
    )
    print(f"exact exhaustive ledger sha256={exact_digest}")
    print("verdict: no counterexample; all exact and hostile gates passed")
    print("scope: complete finite automata in shortlex order, not arbitrary subsets of N")


if __name__ == "__main__":
    main()
