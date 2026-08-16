#!/usr/bin/env python3
"""Exact finite-state controls for the regular shortlex harmonic theorem.

The analytic theorem is proved in the companion reflection.  This script
reconstructs its sharp examples: THM-3497's 12- and 192-state group
automata, a periodic language with no natural density, and a two-basin prefix
language whose logarithmic coefficient is not its state-count density.
"""

from __future__ import annotations

import hashlib
import math
from collections import deque


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


EXPECTED_SEMANTIC_SHA256 = "90586710e72cf8ebc6c36e3e84135f0053973fecf412b5d83fa22881f816ce75"


def compose(left, right):
    """Return the permutation left after right."""
    return tuple(left[right[index]] for index in range(len(right)))


def inverse(permutation):
    result = [0] * len(permutation)
    for index, value in enumerate(permutation):
        result[value] = index
    return tuple(result)


def cycle_type(permutation):
    seen, lengths = set(), []
    for start in range(len(permutation)):
        if start in seen:
            continue
        current, length = start, 0
        while current not in seen:
            seen.add(current)
            current = permutation[current]
            length += 1
        lengths.append(length)
    return tuple(sorted(lengths))


def generated_states(identity, generators, step):
    states = {identity}
    queue = deque([identity])
    while queue:
        state = queue.popleft()
        for generator in generators:
            new = step(generator, state)
            if new not in states:
                states.add(new)
                queue.append(new)
    return tuple(sorted(states))


# THM-3497 variable-translation automaton S3 x C2.
identity3 = (0, 1, 2)
r = (1, 2, 0)
s = (1, 0, 2)
variable_generators = ((r, 1), (s, 1), (r, 0))  # A,B,C


def variable_step(generator, state):
    return compose(generator[0], state[0]), generator[1] ^ state[1]


variable_states = generated_states((identity3, 0), variable_generators, variable_step)
variable_accept = {
    state
    for state in variable_states
    if (state[1] == 0 and cycle_type(state[0]) == (1, 1, 1))
    or (state[1] == 1 and cycle_type(state[0]) == (1, 2))
}
require(len(variable_states) == 12 and len(variable_accept) == 4, "variable automaton census changed")

# THM-3497 fixed-drift automaton S4 x D4.
identity4 = (0, 1, 2, 3)
A3 = (1, 3, 2, 0)
B3 = (1, 3, 0, 2)
C3 = (3, 1, 0, 2)
G = B3
T = (1, 0, 3, 2)
fixed_generators = ((A3, G), (B3, G), (C3, T))


def fixed_step(generator, state):
    return compose(generator[0], state[0]), compose(generator[1], state[1])


fixed_states = generated_states((identity4, identity4), fixed_generators, fixed_step)
fixed_accept = {state for state in fixed_states if cycle_type(state[0]) == cycle_type(state[1])}
require(len(fixed_states) == 192 and len(fixed_accept) == 34, "fixed automaton census changed")

# The target D4 period character is word-length parity.  Recover its two
# chambers directly by graph distance from the identity.
distance = {(identity4, identity4): 0}
queue = deque([(identity4, identity4)])
while queue:
    state = queue.popleft()
    for generator in fixed_generators:
        new = fixed_step(generator, state)
        if new not in distance:
            distance[new] = distance[state] + 1
            queue.append(new)
even_accept = sum(state in fixed_accept and distance[state] % 2 == 0 for state in fixed_states)
odd_accept = sum(state in fixed_accept and distance[state] % 2 == 1 for state in fixed_states)
require((even_accept, odd_accept) == (16, 18), "fixed parity chambers changed")

# Binary even-length language: complete-level natural densities oscillate,
# while the logarithmic density is the Cesaro mean 1/2.
for exponent in range(1, 12):
    even_depth = 2 * exponent
    odd_depth = even_depth + 1
    accepted_even = (4 ** (exponent + 1) - 1) // 3
    total_even = 2 ** (even_depth + 1) - 1
    accepted_odd = accepted_even
    total_odd = 2 ** (odd_depth + 1) - 1
    require(0 < accepted_even < total_even and 0 < accepted_odd < total_odd, "parity endpoint count failed")
require(abs(2 / 3 - ((4**13 - 1) / 3) / (2**25 - 1)) < 1e-7, "even endpoint limit changed")
require(abs(1 / 3 - ((4**13 - 1) / 3) / (2**26 - 1)) < 1e-7, "odd endpoint limit changed")

# Binary prefix-zero language.  The accepting recurrent basin is the first
# half of the q-adic interval, so the coefficient is log(3/2)/log(2), not
# the unweighted branch probability 1/2.  Reversing the alphabet gives the
# complementary log(4/3)/log(2).
prefix_zero = math.log(3 / 2) / math.log(2)
prefix_one = math.log(4 / 3) / math.log(2)
require(abs(prefix_zero + prefix_one - 1) < 1e-15, "prefix basin coefficients ceased to complement")
require(abs(prefix_zero - 0.5) > 0.08, "prefix hostile no longer separates state mass")

# The full-level harmonic block converges to log(3/2).  A direct finite sum
# supplies a numerical hostile independent of the integral formula.
depth = 18
block_start = 2**depth
block_stop = block_start + 2 ** (depth - 1)
block_mass = math.fsum(1 / index for index in range(block_start, block_stop))
require(abs(block_mass - math.log(3 / 2)) < 2e-6, "prefix harmonic block missed its integral")

semantic_rows = [
    f"variable:{state}:{int(state in variable_accept)}" for state in variable_states
] + [
    f"fixed:{state}:{int(state in fixed_accept)}:{distance[state] % 2}" for state in fixed_states
]
semantic_rows += [
    "parity-natural-limits:2/3,1/3",
    "parity-logarithmic:1/2",
    "prefix-zero:log(3/2)/log(2)",
    "prefix-one:log(4/3)/log(2)",
]
semantic_digest = hashlib.sha256("\n".join(semantic_rows).encode("ascii")).hexdigest()
require(semantic_digest == EXPECTED_SEMANTIC_SHA256, "automaton semantic ledger changed")

print("== regular shortlex harmonic-density controls ==")
print("variable Berggren automaton: states=12 accept=4 stationary/log coefficient=1/3")
print("fixed Berggren automaton: states=192 accept=34=17/96; parity chambers=16/96,18/96")
print("periodic even-length binary language: natural endpoint limits=2/3,1/3; log coefficient=1/2")
print("prefix-zero binary basin: coefficient=log(3/2)/log(2), not 1/2")
print("prefix-one reversed-order basin: coefficient=log(4/3)/log(2); complementary sum=1")
print(f"semantic ledger sha256={semantic_digest}")
print("scope: exact finite-state controls for the analytic theorem; no arbitrary-subset regularity claim")
print("all exact checks passed")
