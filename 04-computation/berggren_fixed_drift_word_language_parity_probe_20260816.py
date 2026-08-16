#!/usr/bin/env python3
"""Exact word language for THM-3339's fixed affine branch cocycle.

Unlike the variable-translation language, the target permutation of a word is
fixed by A~=B~=(u -> Ju+p) and C~=(u -> u+p).  A word is accepted precisely
when its projective P1(F3) permutation is conjugate to that fixed affine
permutation.  The resulting 192-state group automaton is bipartite.
"""

from __future__ import annotations

from collections import Counter, deque
from hashlib import sha256
from itertools import product
from json import dumps
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "04-computation/berggren_fixed_drift_word_language_parity_probe_20260816.py"
OUTPUT = "05-knowledge/results/berggren_fixed_drift_word_language_parity_probe_20260816.out"
PIN = (
    ROOT / "05-knowledge/results/berggren_transplantable_word_language_harmonic_probe_20260816.out",
    "8690f9254744ac2fae5ad700e1932f45a7b996cb88be5153d864e6e9f1a02a9a",
)

LETTERS = "ABC"
I4 = (0, 1, 2, 3)
SOURCE = {
    "A": (1, 3, 2, 0),
    "B": (1, 3, 0, 2),
    "C": (3, 1, 0, 2),
}
G = SOURCE["B"]                  # u -> Ju+p in 0,p,q,r order
T = (1, 0, 3, 2)                # u -> u+p
TARGET = {"A": G, "B": G, "C": T}


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def compose(left, right):
    return tuple(left[right[index]] for index in range(len(left)))


def inverse(permutation):
    return tuple(permutation.index(index) for index in range(len(permutation)))


def power(permutation, exponent):
    result = tuple(range(len(permutation)))
    for _ in range(exponent):
        result = compose(permutation, result)
    return result


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


def step(letter, state):
    source, target = state
    return compose(SOURCE[letter], source), compose(TARGET[letter], target)


def word_state(word):
    state = (I4, I4)
    for letter in word:
        state = step(letter, state)
    return state


def accepted(state):
    return cycle_lengths(state[0]) == cycle_lengths(state[1])


def enumerate_group():
    start = (I4, I4)
    states = [start]
    index = {start: 0}
    paths = {start: ""}
    parity = {start: 0}
    queue = deque((start,))
    while queue:
        state = queue.popleft()
        for letter in LETTERS:
            candidate = step(letter, state)
            candidate_parity = parity[state] ^ 1
            if candidate not in index:
                index[candidate] = len(states)
                states.append(candidate)
                paths[candidate] = paths[state] + letter
                parity[candidate] = candidate_parity
                queue.append(candidate)
            else:
                require(parity[candidate] == candidate_parity, (state, letter, candidate))
    return tuple(states), index, paths, parity


def minimal_partition(states, index):
    accepting = {index[state] for state in states if accepted(state)}
    partition = [accepting, set(range(len(states))) - accepting]
    rounds = []
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
                    block_of[index[step(letter, states[state_index])]]
                    for letter in LETTERS
                )
                fibers.setdefault(signature, set()).add(state_index)
            refined.extend(fibers.values())
        rounds.append(len(refined))
        if len(refined) == len(partition):
            return tuple(frozenset(block) for block in partition), tuple(rounds)
        partition = refined


def level_counts(limit):
    distribution = {(I4, I4): 1}
    counts = []
    supports = []
    for _ in range(limit + 1):
        counts.append(sum(weight for state, weight in distribution.items() if accepted(state)))
        supports.append(len(distribution))
        next_distribution = Counter()
        for state, weight in distribution.items():
            for letter in LETTERS:
                next_distribution[step(letter, state)] += weight
        distribution = dict(next_distribution)
    return tuple(counts), tuple(supports)


def generating_function(counts, recurrence):
    x = sp.symbols("x")
    order = len(recurrence)
    denominator = 1 - sum(
        recurrence[index] * x ** (index + 1)
        for index in range(order)
    )
    numerator = sum(counts[index] * x ** index for index in range(order))
    for lag, coefficient in enumerate(recurrence, start=1):
        numerator -= coefficient * x ** lag * sum(
            counts[index] * x ** index for index in range(order - lag)
        )
    return sp.factor(sp.expand(numerator)), sp.factor(sp.expand(denominator))


def main() -> None:
    pin_path, pin_expected = PIN
    pin_actual = lf_sha256(pin_path)
    require(pin_actual == pin_expected, (pin_expected, pin_actual))

    require(cycle_lengths(G) == (4,), cycle_lengths(G))
    require(cycle_lengths(T) == (2, 2), cycle_lengths(T))
    require(power(G, 4) == I4 and power(T, 2) == I4, (G, T))
    require(compose(compose(T, G), T) == inverse(G), "D4 relation")

    states, index, paths, parity = enumerate_group()
    require(len(states) == 192, len(states))
    source_projection = {state[0] for state in states}
    target_projection = {state[1] for state in states}
    pure_source = {state[0] for state in states if state[1] == I4}
    pure_target = {state[1] for state in states if state[0] == I4}
    require((len(source_projection), len(target_projection)) == (24, 8), (len(source_projection), len(target_projection)))
    require((len(pure_source), len(pure_target)) == (24, 8), (len(pure_source), len(pure_target)))

    color_sizes = Counter(parity.values())
    accepting_by_color = Counter(parity[state] for state in states if accepted(state))
    require(color_sizes == Counter({0: 96, 1: 96}), color_sizes)
    require(accepting_by_color == Counter({0: 16, 1: 18}), accepting_by_color)

    target_normal_forms = {}
    for rotation in range(4):
        target_normal_forms[power(G, rotation)] = (rotation, 0)
        target_normal_forms[compose(power(G, rotation), T)] = (rotation, 1)
    require(len(target_normal_forms) == 8, target_normal_forms)
    target_character = {
        target: (rotation + reflection) % 2
        for target, (rotation, reflection) in target_normal_forms.items()
    }
    require(target_character[G] == 1 and target_character[T] == 1, target_character)
    require(all(target_character[state[1]] == parity[state] for state in states), "period character")
    require(word_state("B" * 4) == (I4, I4), "length-four return")
    require(word_state("C" * 6) == (I4, I4), "length-six return")

    minimal_blocks, refinement_rounds = minimal_partition(states, index)
    require(len(minimal_blocks) == 192, len(minimal_blocks))
    require(refinement_rounds == (9, 85, 172, 192, 192), refinement_rounds)

    recurrence = (0, 7, 0, 7, 0, 71, 0, 213, 0, 189, 0, 1701, 0, -2187)
    counts, supports = level_counts(220)
    require(counts[:20] == (
        1, 1, 1, 4, 13, 46, 113, 421, 1121, 3667,
        9801, 33166, 88893, 299248, 796489, 2690239,
        7173601, 24211333, 64565297, 217930552,
    ), counts[:20])
    require(supports[11:] == (96,) * (len(supports) - 11), supports[8:16])

    residuals = tuple(
        counts[n] - sum(
            recurrence[lag - 1] * counts[n - lag]
            for lag in range(1, len(recurrence) + 1)
        )
        for n in range(len(recurrence), len(counts))
    )
    require(residuals == (0,) * len(residuals), residuals[:10])
    require(len(residuals) >= len(states), len(residuals))

    numerator, denominator = generating_function(counts, recurrence)
    x = sp.symbols("x")
    expected_denominator = (
        (x - 1) * (x + 1) * (3 * x - 1) * (3 * x + 1)
        * (3 * x ** 2 + 1)
        * (3 * x ** 2 - 2 * x + 1)
        * (3 * x ** 2 + 2 * x + 1)
        * (9 * x ** 4 - 2 * x ** 2 + 1)
    )
    require(sp.expand(denominator - expected_denominator) == 0, denominator)
    generating_series = numerator / denominator
    coefficient_plus = sp.simplify(sp.limit((1 - 3 * x) * generating_series, x, sp.Rational(1, 3)))
    coefficient_minus = sp.simplify(sp.limit((1 + 3 * x) * generating_series, x, -sp.Rational(1, 3)))
    require(coefficient_plus == sp.Rational(17, 96), coefficient_plus)
    require(coefficient_minus == -sp.Rational(1, 96), coefficient_minus)

    small_words = {}
    for length in range(1, 5):
        small_words[length] = tuple(
            "".join(digits)
            for digits in product(LETTERS, repeat=length)
            if accepted(word_state("".join(digits)))
        )
    require(small_words[1] == ("B",) and small_words[2] == ("BB",), small_words)
    require(all(accepted(word_state("B" * exponent)) for exponent in range(25)), "B spine")

    fibonacci_gate = []
    for ray_index in range(48):
        words = (
            "BA" * ray_index,
            "A" + "BA" * ray_index,
            "C" + "BC" * ray_index,
        )
        observed = tuple(accepted(word_state(word)) for word in words)
        expected = (
            ray_index % 4 == 0,
            ray_index % 4 == 3,
            ray_index % 4 == 3,
        )
        require(observed == expected, (ray_index, words, observed, expected))
        fibonacci_gate.append(observed)

    semantic = {
        "group_order": len(states),
        "projection_orders": (len(source_projection), len(target_projection)),
        "color_sizes": dict(color_sizes),
        "accepting_by_color": dict(accepting_by_color),
        "minimal_states": len(minimal_blocks),
        "refinement_rounds": refinement_rounds,
        "counts": counts[:30],
        "recurrence": recurrence,
        "gf": (str(numerator), str(denominator)),
        "small_words": small_words,
        "fibonacci_gate": fibonacci_gate[:8],
    }
    semantic_sha256 = sha256(
        dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

    print("BERGGREN FIXED-DRIFT WORD LANGUAGE / PARITY-SPLIT EXACT PROBE")
    print("STATUS: VERIFIED-EXACT FINITE AUTOMATON; ASYMPTOTIC CONSEQUENCES ARE ELEMENTARY DEDUCTIONS")
    print(f"SCRIPT: {SCRIPT}")
    print(f"OUTPUT: {OUTPUT}")
    print(f"PIN_VARIABLE_TRANSLATION_OUTPUT_SHA256: {pin_actual}")
    print("FIXED_LIFTS: A~=B~=G=(u |-> Ju+p), cycle type 4; C~=T=(u |-> u+p), cycle type 2^2; TGT=G^(-1)")
    print("ACCEPTANCE: projective rho_3(w) and fixed affine rho_fix(w) have the same S4 cycle type, equivalently are point-conjugate")
    print("PAIR_IMAGE: all S4 x D4, order 24*8=192; pure source and pure target subgroups have orders 24 and 8")
    print("PERIOD_CHARACTER: the D4 abelianization character chi(G)=chi(T)=1 equals word-length parity; every letter crosses the 96+96 bipartition")
    print("TWO_STEP_MIXING_GATE: B^4 and C^6 are returns, of two-step lengths 2 and 3; together with reachability of each 96-state class this gives period one for the two-step walk")
    print("ACCEPTING_STATES_BY_PARITY: even=16 of 96, odd=18 of 96")
    print(f"MINIMAL_DFA: 192 states; partition-refinement block counts {refinement_rounds}; syntactic monoid is the regular S4 x D4 action")
    print(f"SMALL_ACCEPTED_WORDS: {small_words}")
    print("B_SPINE_POSITIVE: every B^n passes under the identity point calibration")
    print(f"LEVEL_COUNTS_A_0_TO_29: {counts[:30]}")
    print("ORDER_14_RECURRENCE_LAGS_1_TO_14: (0,7,0,7,0,71,0,213,0,189,0,1701,0,-2187); 192+ consecutive residual gates vanish, licensing all n by Cayley-Hamilton")
    print(f"GENERATING_FUNCTION_NUMERATOR: {numerator}")
    print(f"GENERATING_FUNCTION_DENOMINATOR: {denominator}")
    print("PARITY_ASYMPTOTIC: a_n=(17/96)3^n-(1/96)(-3)^n+O(3^(n/2)); even levels tend to 1/6, odd levels to 3/16")
    print("SHORTLEX_NATURAL_DENSITY_NO_GO: at complete level boundaries, density tends to 11/64 along even levels and 35/192 along odd levels")
    print("SHORTLEX_HARMONIC_DENSITY: two-level logarithmic averaging gives sum_{m<=N,accepted}1/m=(17/96)log N+o(log N), even though natural density fails")
    print("FIBONACCI_FIXED_DRIFT_GATE: (BA)^r passes iff r=0 mod 4; A(BA)^r and C(BC)^r pass iff r=3 mod 4; equivalently k mod 12 in {0,1,2}")
    print("FIBONACCI_INDEX_HARMONIC: the fixed-drift passing indices have coefficient 1/4 in the harmonic series; actual Fibonacci/triple-coordinate reciprocal weights still converge")
    print("TYPING: the point gauge may still depend on the composite word, but the affine target drift is fixed globally; this is not one simultaneous branch calibration and has no LRC/Jacobian consequence")
    print(f"SEMANTIC_SHA256: {semantic_sha256}")
    print("VERDICT: retaining THM-3339's drift doubles the Fibonacci gate from mod 6 to mod 12 and turns the clean one-third tree density into an exact parity-split language whose harmonic density survives at 17/96")


if __name__ == "__main__":
    main()
