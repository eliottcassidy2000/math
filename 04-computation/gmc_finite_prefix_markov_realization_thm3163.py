#!/usr/bin/env python3
"""Exact controls for THM-3163's finite-prefix Markov realization.

The theorem is an elementary finite identity.  This companion independently
propagates the posterior Markov kernel for a dense four-label law and for the
labelled lift of THM-3158's seven-state product-Gamma selector law.
"""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, product
from math import comb


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def submasks(mask):
    answer = []
    current = mask
    while True:
        answer.append(current)
        if current == 0:
            return tuple(answer)
        current = (current - 1) & mask


def posterior_realization(number_of_labels, terminal_law):
    """Build and independently propagate the posterior stopping kernel."""

    require(all(weight > 0 for weight in terminal_law.values())
            and sum(terminal_law.values(), Fraction(0)) == 1,
            "terminal law is not a probability law")
    reachable = set()
    for target in terminal_law:
        require(target >> number_of_labels == 0,
                "terminal mask left the labelled universe")
        reachable.update(submasks(target))
    reachable = tuple(sorted(reachable, key=lambda mask: (
        mask.bit_count(), mask)))

    reach = {}
    stop = {}
    transition = {}
    normalization_checks = 0
    for state in reachable:
        size = state.bit_count()
        reach[state] = sum(
            weight / comb(target.bit_count(), size)
            for target, weight in terminal_law.items()
            if state & target == state
        )
        require(reach[state] > 0, "reachable state acquired zero mass")
        stop[state] = terminal_law.get(state, Fraction(0)) / reach[state]
        outgoing = {}
        for label in range(number_of_labels):
            bit = 1 << label
            if state & bit:
                continue
            numerator = sum(
                weight
                / (comb(target.bit_count(), size)
                   * (target.bit_count() - size))
                for target, weight in terminal_law.items()
                if target & (state | bit) == state | bit
            )
            outgoing[label] = numerator / reach[state]
        require(stop[state] + sum(outgoing.values(), Fraction(0)) == 1,
                "posterior kernel lost stochastic normalization")
        require(all(value >= 0 for value in outgoing.values()),
                "posterior transition became negative")
        transition[state] = outgoing
        normalization_checks += 1

    arrival = defaultdict(Fraction)
    terminal = defaultdict(Fraction)
    arrival[0] = Fraction(1)
    for state in reachable:
        require(arrival[state] == reach[state],
                "independent propagation disagreed with reach formula")
        terminal[state] += arrival[state] * stop[state]
        for label, probability in transition[state].items():
            arrival[state | (1 << label)] += arrival[state] * probability

    require({state: mass for state, mass in terminal.items() if mass}
            == terminal_law,
            "propagated stopping distribution drifted from target law")
    return {
        "reachable": len(reachable),
        "normalization_checks": normalization_checks,
        "transition_count": sum(
            value > 0 for outgoing in transition.values()
            for value in outgoing.values()),
        "reach": reach,
        "stop": stop,
        "transition": transition,
    }


# Dense control, including positive mass at the empty state.
DENSE_DENOMINATOR = sum(range(1, 17))
DENSE_LAW = {
    mask: Fraction(mask + 1, DENSE_DENOMINATOR) for mask in range(16)
}
DENSE = posterior_realization(4, DENSE_LAW)
require(DENSE["reachable"] == 16,
        "dense four-label fixture lost a state")


# The exact THM-3158 law, lifted uniformly to the 16 labelled pole copies.
POLE_VALUES = (1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 6, 7, 8)
VALUE_LABELS = {
    value: tuple(index for index, old in enumerate(POLE_VALUES)
                 if old == value)
    for value in sorted(set(POLE_VALUES))
}
SELECTOR_NUMERATORS = {
    (1,): 97856,
    (1, 1, 1): 56951,
    (1, 1, 1, 1): 140643,
    (1, 1, 2, 2, 2): 194498,
    (1, 2, 2, 3, 3): 7398,
    (2, 2, 2, 3, 3): 118572,
    (5, 5, 6, 7, 8): 384082,
}


def labelled_lifts(state):
    counts = Counter(state)
    choices = tuple(
        tuple(combinations(VALUE_LABELS[value], multiplicity))
        for value, multiplicity in sorted(counts.items())
    )
    answer = []
    for blocks in product(*choices):
        mask = 0
        for block in blocks:
            for label in block:
                mask |= 1 << label
        answer.append(mask)
    return tuple(answer)


SELECTOR_LAW = {}
for state, numerator in SELECTOR_NUMERATORS.items():
    lifts = labelled_lifts(state)
    weight = Fraction(numerator, 10**6 * len(lifts))
    for mask in lifts:
        require(mask not in SELECTOR_LAW,
                "distinct selector states acquired the same labelled lift")
        SELECTOR_LAW[mask] = weight

require(len(SELECTOR_LAW) == 29
        and sum(SELECTOR_LAW.values(), Fraction(0)) == 1,
        "labelled selector lift census/normalization drift")
SELECTOR = posterior_realization(len(POLE_VALUES), SELECTOR_LAW)

# Project terminal mass back to value multisets.
projected = defaultdict(Fraction)
for mask, mass in SELECTOR_LAW.items():
    state = tuple(sorted(
        POLE_VALUES[label] for label in range(len(POLE_VALUES))
        if mask & (1 << label)))
    projected[state] += mass
require(projected == {
    state: Fraction(numerator, 10**6)
    for state, numerator in SELECTOR_NUMERATORS.items()
}, "labelled lift failed to recover the unlabeled selector law")

# The same support cannot come from independent Bernoulli deletion, even after
# conditioning on a nonempty outcome.  All four singleton 1-label sets and the
# set of all four 1-labels have positive mass, whereas every two-1 set is zero.
ONE_LABELS = VALUE_LABELS[1]
ONE_SINGLETON_MASKS = tuple(1 << label for label in ONE_LABELS)
ALL_ONES_MASK = sum(ONE_SINGLETON_MASKS)
TWO_ONES_MASK = ONE_SINGLETON_MASKS[0] | ONE_SINGLETON_MASKS[1]
require(all(mask in SELECTOR_LAW for mask in ONE_SINGLETON_MASKS)
        and ALL_ONES_MASK in SELECTOR_LAW
        and TWO_ONES_MASK not in SELECTOR_LAW,
        "independent-deletion support hostile drift")


# Nonuniqueness: delta_{0,1} is realized by either deterministic order.
def deterministic_order_kernel(order):
    state = 0
    kernel = {}
    for label in order:
        kernel[state] = (Fraction(0), {label: Fraction(1)})
        state |= 1 << label
    kernel[state] = (Fraction(1), {})
    return kernel


def propagate_explicit_kernel(kernel):
    arrival = defaultdict(Fraction)
    terminal = defaultdict(Fraction)
    arrival[0] = Fraction(1)
    for state in sorted(kernel, key=lambda mask: (mask.bit_count(), mask)):
        stop_probability, outgoing = kernel[state]
        terminal[state] += arrival[state] * stop_probability
        for label, probability in outgoing.items():
            arrival[state | (1 << label)] += arrival[state] * probability
    return {state: mass for state, mass in terminal.items() if mass}


FIRST_ORDER_KERNEL = deterministic_order_kernel((0, 1))
SECOND_ORDER_KERNEL = deterministic_order_kernel((1, 0))
require(FIRST_ORDER_KERNEL != SECOND_ORDER_KERNEL
        and propagate_explicit_kernel(FIRST_ORDER_KERNEL)
        == propagate_explicit_kernel(SECOND_ORDER_KERNEL)
        == {0b11: Fraction(1)},
        "deterministic-order nonuniqueness control failed")


print("THM-3163 universal finite-prefix Markov realization")
print("dense_n4_terminal_states=16")
print(f"dense_n4_reachable_states={DENSE['reachable']}")
print(f"dense_n4_kernel_normalization_checks={DENSE['normalization_checks']}")
print(f"dense_n4_positive_transitions={DENSE['transition_count']}")
print("selector_pole_multiset=" + repr(POLE_VALUES))
print("selector_unlabeled_terminal_states=7")
print(f"selector_labelled_terminal_states={len(SELECTOR_LAW)}")
print(f"selector_reachable_prefix_states={SELECTOR['reachable']}")
print(f"selector_kernel_normalization_checks={SELECTOR['normalization_checks']}")
print(f"selector_positive_transitions={SELECTOR['transition_count']}")
print("selector_unlabeled_projection=EXACT")
print("selector_independent_deletion_support=IMPOSSIBLE")
print("same_terminal_law_distinct_order_kernels=2")
print("scope=abstract_state_dependent_chain_not_response_transport")
print("all_exact_checks=PASS")
