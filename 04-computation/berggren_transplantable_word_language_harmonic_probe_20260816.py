#!/usr/bin/env python3
"""Exact automaton for wordwise Berggren-to-T4 affine calibratability.

A branch word is accepted when its action on P1(F3) is conjugate to an
affine permutation of V4 whose linear part is the word's true mod-two
direction action.  The 48-state raw action factors through a 12-state
matching-action/linear-bit automaton.  This script proves the level recurrence
and records the transfer spectrum used for density and harmonic deductions.
"""

from __future__ import annotations

from collections import Counter, deque
from hashlib import sha256
from itertools import permutations, product
from json import dumps
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "04-computation/berggren_transplantable_word_language_harmonic_probe_20260816.py"
OUTPUT = "05-knowledge/results/berggren_transplantable_word_language_harmonic_probe_20260816.out"
PIN = (
    ROOT / "05-knowledge/results/berggren_full_branch_t4_static_calibration_no_go_20260816.out",
    "77b7a038cf4a6fdb7f1b609729f5e387af8ea8290ef2a976c31d2141e6e85dc3",
)

LETTERS = "ABC"
P1_ACTION = {
    "A": (1, 3, 2, 0),
    "B": (1, 3, 0, 2),
    "C": (3, 1, 0, 2),
}
I4 = (0, 1, 2, 3)
J4 = (0, 2, 1, 3)
LINEAR_ACTION = {"A": J4, "B": J4, "C": I4}
LINEAR_BIT = {"A": 1, "B": 1, "C": 0}

MATCHINGS = (
    ((0, 1), (2, 3)),
    ((0, 2), (1, 3)),
    ((0, 3), (1, 2)),
)
I3 = (0, 1, 2)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def compose(left, right):
    return tuple(left[right[index]] for index in range(len(left)))


def cycle_lengths(permutation):
    unseen = set(range(len(permutation)))
    lengths = []
    while unseen:
        start = min(unseen)
        current = start
        length = 0
        while current in unseen:
            unseen.remove(current)
            current = permutation[current]
            length += 1
        require(current == start, (permutation, start, current))
        lengths.append(length)
    return tuple(sorted(lengths))


def canonical_matching(matching):
    return tuple(sorted(tuple(sorted(edge)) for edge in matching))


def matching_action(permutation):
    result = []
    for matching in MATCHINGS:
        image = canonical_matching(
            tuple((permutation[left], permutation[right]) for left, right in matching)
        )
        result.append(MATCHINGS.index(image))
    return tuple(result)


MATCHING_ACTION = {
    letter: matching_action(permutation)
    for letter, permutation in P1_ACTION.items()
}


def state_step(letter, state):
    matching_permutation, linear_bit = state
    return (
        compose(MATCHING_ACTION[letter], matching_permutation),
        LINEAR_BIT[letter] ^ linear_bit,
    )


def is_transposition(permutation):
    return cycle_lengths(permutation) == (1, 2)


def accepted(state):
    matching_permutation, linear_bit = state
    return (
        (linear_bit == 0 and matching_permutation == I3)
        or (linear_bit == 1 and is_transposition(matching_permutation))
    )


def word_state(word):
    state = (I3, 0)
    for letter in word:
        state = state_step(letter, state)
    return state


def affine_cycle_types(linear):
    # All translations of the specified linear permutation on V4.  XOR in
    # the fixed point order 0,p,q,r is simply integer XOR.
    return {
        cycle_lengths(tuple(linear[index] ^ translation for index in range(4)))
        for translation in range(4)
    }


def raw_affine_compatible(projective, linear):
    return cycle_lengths(projective) in affine_cycle_types(linear)


def enumerate_states():
    start = (I3, 0)
    states = [start]
    index = {start: 0}
    queue = deque((start,))
    while queue:
        state = queue.popleft()
        for letter in LETTERS:
            candidate = state_step(letter, state)
            if candidate not in index:
                index[candidate] = len(states)
                states.append(candidate)
                queue.append(candidate)
    return tuple(states), index


def transfer_matrix(states, index):
    matrix = sp.zeros(len(states))
    for state in states:
        for letter in LETTERS:
            matrix[index[state_step(letter, state)], index[state]] += 1
    return matrix


def level_counts(limit):
    distribution = {(I3, 0): 1}
    counts = []
    supports = []
    for _ in range(limit + 1):
        counts.append(sum(weight for state, weight in distribution.items() if accepted(state)))
        supports.append(len(distribution))
        next_distribution = Counter()
        for state, weight in distribution.items():
            for letter in LETTERS:
                next_distribution[state_step(letter, state)] += weight
        distribution = dict(next_distribution)
    return tuple(counts), tuple(supports)


def minimal_partition(states, index):
    accepting = {index[state] for state in states if accepted(state)}
    partition = [accepting, set(range(len(states))) - accepting]
    while True:
        block_of = {
            state_index: block_index
            for block_index, block in enumerate(partition)
            for state_index in block
        }
        refined = []
        for block in partition:
            fibers = {}
            for state_index in block:
                signature = tuple(
                    block_of[index[state_step(letter, states[state_index])]]
                    for letter in LETTERS
                )
                fibers.setdefault(signature, set()).add(state_index)
            refined.extend(fibers.values())
        if len(refined) == len(partition):
            return tuple(frozenset(block) for block in partition)
        partition = refined


def oscillatory_integer_sequence(limit):
    values = [1]
    if limit:
        values.append(-1)
    while len(values) <= limit:
        values.append(-2 * values[-1] - 3 * values[-2])
    return tuple(values)


def main() -> None:
    pin_path, pin_expected = PIN
    pin_actual = lf_sha256(pin_path)
    require(pin_actual == pin_expected, (pin_expected, pin_actual))

    require(MATCHING_ACTION == {
        "A": (1, 2, 0),
        "B": (1, 0, 2),
        "C": (1, 2, 0),
    }, MATCHING_ACTION)

    states, index = enumerate_states()
    require(len(states) == 12, len(states))
    accepted_states = tuple(state for state in states if accepted(state))
    require(len(accepted_states) == 4, accepted_states)
    minimal_blocks = minimal_partition(states, index)
    require(len(minimal_blocks) == 12, minimal_blocks)

    # The quotient criterion is equivalent to the raw four-point cycle-type
    # criterion on every element of S4 x <J>, not merely on sampled words.
    raw_pairs = tuple(product(permutations(range(4)), (I4, J4)))
    for projective, linear in raw_pairs:
        quotient_state = (matching_action(projective), int(linear == J4))
        require(
            accepted(quotient_state) == raw_affine_compatible(projective, linear),
            (projective, linear, quotient_state),
        )
    require(sum(raw_affine_compatible(*pair) for pair in raw_pairs) == 16, "48-state count")

    matrix = transfer_matrix(states, index)
    variable = sp.symbols("lambda")
    characteristic = sp.factor(matrix.charpoly(variable).as_expr())
    expected_characteristic = (
        (variable - 3) * (variable - 1) ** 4 * (variable + 1) ** 3
        * (variable ** 2 + 2 * variable + 3) ** 2
    )
    require(sp.expand(characteristic - expected_characteristic) == 0, characteristic)

    identity_matrix = sp.eye(12)
    squarefree_annihilator = (
        (matrix - 3 * identity_matrix)
        * (matrix - identity_matrix)
        * (matrix + identity_matrix)
        * (matrix ** 2 + 2 * matrix + 3 * identity_matrix)
    )
    require(squarefree_annihilator == sp.zeros(12), "transfer is not squarefree-annihilated")

    counts, supports = level_counts(40)
    require(counts[:12] == (1, 1, 3, 11, 25, 81, 251, 715, 2193, 6593, 19603, 59115), counts[:12])
    require(supports[6:] == (12,) * (len(supports) - 6), supports)

    recurrence_ok = all(
        counts[n] == 2 * counts[n - 1] + 2 * counts[n - 2]
        + 6 * counts[n - 3] - 9 * counts[n - 4]
        for n in range(4, len(counts))
    )
    require(recurrence_ok, counts)

    # Exact all-n recurrence gate: the scalar recurrence polynomial vanishes
    # on the full reachable Krylov space.
    start = sp.zeros(12, 1)
    start[index[(I3, 0)]] = 1
    accept_row = sp.Matrix([[int(accepted(state)) for state in states]])
    krylov = sp.Matrix.hstack(*(matrix ** exponent * start for exponent in range(12)))
    scalar_polynomial = (
        matrix ** 4 - 2 * matrix ** 3 - 2 * matrix ** 2
        - 6 * matrix + 9 * identity_matrix
    )
    require(krylov.rank() == 5, krylov.rank())
    require(accept_row * scalar_polynomial * krylov == sp.zeros(1, 12), "recurrence gate")

    oscillation = oscillatory_integer_sequence(len(counts) - 1)
    closed_form = tuple((3 ** n + 1 + oscillation[n]) // 3 for n in range(len(counts)))
    require(all((3 ** n + 1 + oscillation[n]) % 3 == 0 for n in range(len(counts))), oscillation)
    require(closed_form == counts, (closed_form, counts))

    accepted_small_words = {}
    for length in range(1, 4):
        words = []
        for digits in product(LETTERS, repeat=length):
            word = "".join(digits)
            if accepted(word_state(word)):
                words.append(word)
        accepted_small_words[length] = tuple(words)
    require(accepted_small_words == {
        1: ("B",),
        2: ("BB", "BC", "CB"),
        3: ("AAB", "AAC", "ABA", "ACA", "BAA", "BBB", "BCC", "CAA", "CBC", "CCB", "CCC"),
    }, accepted_small_words)
    require(
        accepted(word_state("B"))
        and accepted(word_state("BC"))
        and not accepted(word_state("BBC")),
        "submonoid hostile",
    )

    # THM-3339's three interleaved Fibonacci addresses.  The right-address
    # strings are read root-to-child, while state_step left-multiplies the new
    # branch action, matching the theorem's "matrix products act rightmost
    # first" convention.
    fibonacci_ray_gate = []
    for ray_index in range(24):
        ray_words = (
            "BA" * ray_index,
            "A" + "BA" * ray_index,
            "C" + "BC" * ray_index,
        )
        observed = tuple(accepted(word_state(word)) for word in ray_words)
        expected = (
            ray_index % 2 == 0,
            ray_index % 2 == 1,
            ray_index % 2 == 1,
        )
        require(observed == expected, (ray_index, ray_words, observed, expected))
        fibonacci_ray_gate.append(observed)

    generating_numerator = "1-x-x^2-3x^3"
    generating_denominator = "(1-x)(1-3x)(1+2x+3x^2)"
    semantic = {
        "matching_action": MATCHING_ACTION,
        "accepted_states": accepted_states,
        "characteristic": str(characteristic),
        "counts": counts[:20],
        "recurrence": (2, 2, 6, -9),
        "accepted_small_words": accepted_small_words,
        "minimal_dfa_states": len(minimal_blocks),
        "fibonacci_ray_gate": fibonacci_ray_gate[:4],
        "generating_function": (generating_numerator, generating_denominator),
    }
    semantic_sha256 = sha256(
        dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

    print("BERGGREN TRANSPLANTABLE-WORD LANGUAGE / HARMONIC SUBSET EXACT PROBE")
    print("STATUS: VERIFIED-EXACT FINITE AUTOMATON; DENSITY/HARMONIC CONSEQUENCES ARE ELEMENTARY DEDUCTIONS")
    print(f"SCRIPT: {SCRIPT}")
    print(f"OUTPUT: {OUTPUT}")
    print(f"PIN_FULL_BRANCH_OUTPUT_SHA256: {pin_actual}")
    print(f"MATCHING_QUOTIENT_GENERATORS_A_B_C: {MATCHING_ACTION}")
    print("LINEAR_BITS_A_B_C: (1,1,0), recording whether the true mod-two direction action is J or I")
    print("TWELVE_STATE_IMAGE: S3(perfect-matchings) x C2(linear bit); all 12 states reached by depth 6")
    print("MINIMAL_DFA_AND_SYNTACTIC_MONOID: partition refinement leaves all 12 states distinct; the transition monoid is the regular left action of S3 x C2, hence is the same 12-element group")
    print(f"ACCEPTED_STATES: {accepted_states}")
    print("ACCEPTANCE_RULE: (mu,eps) is wordwise affine-calibratable iff eps=0 and mu=id, or eps=1 and mu is a transposition")
    print("RAW_EQUIVALENCE_GATE: checked on all 48 pairs S4 x <J>; exactly 16 pass, so the accepted-state fraction is 4/12=1/3")
    print(f"SMALL_ACCEPTED_WORDS: {accepted_small_words}")
    print("NOT_A_SUBMONOID_HOSTILE: B and BC pass, but BBC fails")
    print("FIBONACCI_THREE_RAY_RESTRICTION: (BA)^r passes iff r is even; A(BA)^r and C(BC)^r pass iff r is odd; equivalently the THM-3339 Fibonacci index k passes exactly for k mod 6 in {0,1,2}")
    print(f"LEVEL_COUNTS_A_0_TO_19: {counts[:20]}")
    print("RECURRENCE: a_n=2a_(n-1)+2a_(n-2)+6a_(n-3)-9a_(n-4), n>=4; initials (1,1,3,11)")
    print(f"GENERATING_FUNCTION: ({generating_numerator})/({generating_denominator})")
    print("CLOSED_FORM: a_n=(3^n+1+c_n)/3, where c_0=1, c_1=-1, c_n=-2c_(n-1)-3c_(n-2)=Re((-1+i*sqrt(2))^n)")
    print("LEVEL_DENSITY: a_n/3^n=1/3+O(3^(-n/2)); in particular exactly one third asymptotically, not every branch word")
    print(f"TRANSFER_CHARACTERISTIC: {characteristic}")
    print("TRANSFER_SPECTRAL_GATE: squarefree annihilator (lambda-3)(lambda-1)(lambda+1)(lambda^2+2lambda+3); every non-Perron eigenvalue has modulus <=sqrt(3)")
    print("SHORTLEX_SUBSET: index the empty word and successive ternary levels consecutively; the accepted indices have A(N)=N/3+O(sqrt(N))")
    print("HARMONIC_SUBSERIES: partial summation gives sum_{m<=N, accepted} 1/m=(1/3)log N+C+O(N^(-1/2)); the coefficient remembers automaton density, not ancestry labels")
    print("TYPING: calibration may depend on the composite word; this is a regular language of individually calibratable addresses, not one monoid-equivariant calibration and not a physical current")
    print(f"SEMANTIC_SHA256: {semantic_sha256}")
    print("VERDICT: the ternary branch tree contains an exact density-one-third, harmonically divergent subset of wordwise T4-calibratable addresses, recognized by the matching-action/linear-bit automaton")


if __name__ == "__main__":
    main()
