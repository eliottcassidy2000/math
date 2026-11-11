#!/usr/bin/env python3
"""Exact finite controls for the local C1--C3 proofs; no external packages."""

from collections import defaultdict
from hashlib import sha256
import json
from math import gcd
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def step(x, sign):
    return 3 * x + sign if x % 2 else x // 2


def modulus_step(m, sign):
    require(m % 8 in (sign % 8, (4 + sign) % 8), "modulus domain")
    return 3 * m + 2 * sign if m % 8 == (4 + sign) % 8 else (m + sign) // 2


def cylinders(length, sign):
    """Generate legal words without enumerating starting integers.

    State tuple is (residue, exponent, odd_count, even_count, positive_carry).
    """
    states = {"": (0, 0, 0, 0, 0)}
    for _ in range(length):
        next_states = {}
        for word, (residue, exponent, a, b, carry) in states.items():
            choices = "E" if word.endswith("O") else "EO"
            for letter in choices:
                forced = word.endswith("O")
                new_exponent = exponent if forced else exponent + 1
                candidates = (residue,) if forced else (residue, residue + 2**exponent)
                selected = []
                for candidate in candidates:
                    numerator = 3**a * candidate + sign * carry
                    require(numerator % 2**b == 0, "source certificate divisibility")
                    value = numerator // 2**b
                    if value % 2 == (letter == "O"):
                        selected.append(candidate)
                require(len(selected) == 1, "unique cylinder lift")
                if letter == "O":
                    next_state = (selected[0], new_exponent, a + 1, b, 3 * carry + 2**b)
                else:
                    next_state = (selected[0], new_exponent, a, b + 1, carry)
                next_states[word + letter] = next_state
        states = next_states
    return states


def fib(index):
    a, b = 0, 1
    for _ in range(index):
        a, b = b, a + b
    return a


def main():
    central_count = 0
    for m in range(3, 502, 2):
        k = (m - 1) // 2
        actual = [[i * j % m for j in (k, k + 1)] for i in (k, k + 1)]
        d = pow(4, -1, m)
        require(actual == [[d, m - d], [m - d, d]], f"central inverse at {m}")
        if m % 4 == 1:
            x = (m - 1) // 4
            expected = [[3 * x + 1, x], [x, 3 * x + 1]]
        else:
            x = (m + 1) // 4
            expected = [[x, 3 * x - 1], [3 * x - 1, x]]
        require(actual == expected, f"central affine at {m}")
        require(all(sum(row) == m for row in actual), f"row sums at {m}")
        require(all(actual[0][j] + actual[1][j] == m for j in (0, 1)), f"column sums at {m}")
        require(all(gcd(value, m) == 1 for row in actual for value in row), f"entry units at {m}")
        require(gcd(k, m) == gcd(k + 1, m) == 1, f"label units at {m}")
        central_count += 1

    codec_count = 0
    for sign in (-1, 1):
        for x in range(1, 1001):
            m = 4 * x + sign
            y = 3 * x + sign
            if x % 2:
                require(m % 8 == (4 + sign) % 8, "odd mod8 branch")
                require((y - sign) % 6 != 0, "odd division guard")
                decoded = y
            else:
                require(m % 8 == sign % 8, "even mod8 branch")
                require((y - sign) % 6 == 0, "even division guard")
                decoded = (y - sign) // 6
            require(decoded == step(x, sign), "decoded edge")
            require(modulus_step(m, sign) == 4 * step(x, sign) + sign, "graph conjugacy")
            codec_count += 1

    direct_inputs = 0
    direct_step_count = 0
    word_counts = []
    for length in range(13):
        per_sign_count = []
        for sign in (-1, 1):
            generated = cylinders(length, sign)
            observed = defaultdict(list)
            for x in range(1, 2**length + 1):
                value, word = x, ""
                for _ in range(length):
                    word += "O" if value % 2 else "E"
                    value = step(value, sign)
                    direct_step_count += 1
                require("OO" not in word, "no consecutive ordinary odd steps")
                observed[word].append(x)
                require(word in generated, "direct word generated")
                residue, exponent, a, b, carry = generated[word]
                require(exponent == b + int(word.endswith("O")), "final-O precision")
                require(x % 2**exponent == residue, "source congruence")
                require(3**a * x + sign * carry == 2**b * value, "ordered affine endpoint")
                direct_inputs += 1
            require(set(observed) == set(generated), "all generated words directly realized")
            for word, sources in observed.items():
                residue, exponent, _, _, _ = generated[word]
                expected = [x for x in range(1, 2**length + 1) if x % 2**exponent == residue]
                require(sources == expected, "complete source fibre")
            require(len(generated) == fib(length + 2), "Fibonacci word count")
            per_sign_count.append(len(generated))
        require(per_sign_count[0] == per_sign_count[1], "same two-sign language size")
        word_counts.append(per_sign_count[0])

    require(step(2, 1) == 1 and step(2, 1) != 7, "raw-even-pair hostile")
    require(step(7, 1) == 22 and step(7, 1) != 2, "literal-reversal hostile")
    require((5 + 1) // 6 == step(2, -1) == 1, "minus decoding control")
    require(cylinders(1, 1)["O"][:2] == (1, 1), "final-O bit hostile")
    plus = cylinders(3, 1)
    require(plus["OEE"] == (1, 2, 1, 2, 1), "OEE carry")
    require(plus["EOE"] == (2, 2, 1, 2, 2), "EOE carry")
    plus_cycle = [5]
    for _ in range(3):
        plus_cycle.append(modulus_step(plus_cycle[-1], 1))
    require(plus_cycle == [5, 17, 9, 5], "prime-colour cycle")
    minus_cycle = [3]
    for _ in range(2):
        minus_cycle.append(modulus_step(minus_cycle[-1], -1))
    require(minus_cycle == [3, 7, 3], "minus base cycle")

    root = Path(__file__).resolve().parents[4]
    paths = [
        "05-knowledge/results/central_modular_chanpak_20261110/central/check.py",
        "05-knowledge/results/tiling_modular_atoms_20261003.md",
        "01-canon/theorems/THM-4792-two-sheets-mod-odd-m-merged-centre-reflections-unit-shift-and-alternation-classification.md",
        "01-canon/theorems/THM-4528-collatz-parity-in-base-phi-golden-beta-map-and-holonomy.md",
        "05-knowledge/results/odd_eighths_uniform_families_20261030/COLLATZ.md",
    ]
    report = {
        "status": "FINITE-EXACT: all checks passed",
        "central_universe": {"odd_moduli_min": 3, "odd_moduli_max": 501, "count": central_count},
        "codec_universe": {"signs": [-1, 1], "positive_x_min": 1, "positive_x_max": 1000, "count": codec_count},
        "word_universe": {"signs": [-1, 1], "lengths": [0, 12], "representatives": "1 through 2^length", "direct_inputs": direct_inputs, "direct_steps": direct_step_count},
        "word_counts_by_length": word_counts,
        "hostile_controls": ["raw even pair is not an edge", "literal reversal is not an edge", "final-O bit is necessary", "operation counts do not determine carry"],
        "plus_modulus_cycle": plus_cycle,
        "minus_modulus_cycle": minus_cycle,
        "source_sha256": {path: sha256((root / path).read_bytes()).hexdigest() for path in paths},
        "limits": ["No Collatz convergence conclusion", "No tournament count-realization theorem", "Prime-only moduli do not form a closed graph"],
    }
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
