"""Exact signed Collatz word replay; stdlib only and safe under python -O.

Run with --output PATH to save the deterministic JSON; otherwise print it.
The complete main universe is every positive ordered exponent word of length
1 <= L <= 10 and L <= K <= 2L. The note proves this covers every signed
integer cycle at b=+/-1 whose least odd period is at most ten, at any height.
This is an independent replay, not a stronger census than prior repo work.
"""

import argparse
from fractions import Fraction
from math import comb, gcd
import json
from pathlib import Path


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def compositions(total, length):
    if length == 1:
        yield (total,)
        return
    for first in range(1, total - length + 2):
        for rest in compositions(total - first, length - 1):
            yield (first,) + rest


def word_data(word):
    carry, total = 0, 0
    for exponent in word:
        check(exponent >= 1, "nonpositive exponent")
        carry = 3 * carry + 2**total
        total += exponent
    gap = 2**total - 3**len(word)
    check(gap != 0 and carry % 2 == 1, "bad nonempty word")
    return total, carry, gap, abs(gap) // gcd(carry, abs(gap))


def step(n, b):
    check(n % 2 != 0 and b % 2 != 0, "step domain")
    numerator = 3 * n + b
    check(numerator != 0, "zero numerator is outside odd dynamics")
    exponent = 0
    while numerator % 2 == 0:
        numerator //= 2
        exponent += 1
    return numerator, exponent


def canonical_cycle(nodes):
    # First return removes repetition. Rotations, never reversals, are identified.
    first = nodes[0]
    period = next((j for j in range(1, len(nodes)) if nodes[j] == first), len(nodes))
    primitive = tuple(nodes[:period])
    check(len(set(primitive)) == period, "unexpected non-simple cycle")
    index = min(range(period), key=lambda j: (abs(primitive[j]), primitive[j]))
    return primitive[index:] + primitive[:index]


def reconstruct(word, b):
    _, carry, gap, parameter = word_data(word)
    check(b % parameter == 0, "not an integral supporting parameter")
    check((b * carry) % gap == 0, "gate iff failed")
    n0 = b * carry // gap
    n = n0
    nodes = []
    for exponent in word:
        check(n % 2 != 0 and n != 0, "nonzero odd node required")
        nodes.append(n)
        n, actual = step(n, b)
        check(actual == exponent, "formal word does not equal actual valuation")
    check(n == n0, "cycle did not close")
    content = 0
    for node in nodes:
        content = gcd(content, node)
    check(content == abs(b) // parameter, "content theorem failed")
    return canonical_cycle(nodes)


def main():
    total_words = accepted_words = 0
    cycles = set()
    negative_parameter_cycles = set()
    by_length = []
    for length in range(1, 11):
        count = accepted = 0
        for total in range(length, 2 * length + 1):
            for word in compositions(total, length):
                count += 1
                k, carry, gap, parameter = word_data(word)
                check(k == total, "composition total mismatch")
                if parameter == 1:
                    accepted += 1
                    cycle = reconstruct(word, 1)
                    reflected = reconstruct(word, -1)
                    check(reflected == tuple(-n for n in cycle), "conjugacy failure")
                    cycles.add(cycle)
                    negative_parameter_cycles.add(reflected)
        check(count == comb(2 * length, length), "universe count mismatch")
        total_words += count
        accepted_words += accepted
        by_length.append({"length": length, "word_count": count, "accepted_word_count": accepted})

    expected = {(1,), (-1,), (-5, -7), (-17, -25, -37, -55, -41, -61, -91)}
    check(cycles == expected, "unexpected complete bounded-period census")
    check(total_words == 250952 and accepted_words == 37, "headline counts differ")
    cycle_rows = []
    for cycle in sorted(cycles, key=lambda c: (len(c), abs(c[0]), c)):
        word = tuple(step(n, 1)[1] for n in cycle)
        total, carry, gap, parameter = word_data(word)
        check(reconstruct(word, 1) == cycle, "least-period replay")
        cycle_rows.append({"nodes": cycle, "least_period": len(cycle), "word": word,
                           "K": total, "B": carry, "Delta": gap, "q": parameter})

    # Exhaustive independent parameter/rotation checks in a smaller word universe.
    parameter_checks = rotation_checks = 0
    for length in range(1, 7):
        for total in range(length, 2 * length + 1):
            for word in compositions(total, length):
                _, carry, gap, parameter = word_data(word)
                rotated = word[1:] + word[:1]
                _, next_carry, next_gap, next_parameter = word_data(rotated)
                check(next_carry * 2**word[0] == 3 * carry + gap, "rotation carry")
                check((next_gap, next_parameter) == (gap, parameter), "rotation denominator")
                rotation_checks += 1
                for b in range(-15, 16, 2):
                    check((b * carry) % gap == 0 if b % parameter == 0 else (b * carry) % gap != 0,
                          "parameter divisibility gate")
                    if b % parameter == 0:
                        reconstruct(word, b)
                    parameter_checks += 1

    # All unit-gap clocks in this finite box supplement the all-exponent proof.
    unit_gaps = [(k, length) for k in range(1, 65) for length in range(1, min(k, 64) + 1)
                 if abs(2**k - 3**length) == 1]
    check(unit_gaps == [(1, 1), (2, 1), (3, 2)], "unit-gap hostile")
    unit_rows = []
    for k, length in unit_gaps:
        for word in compositions(k, length):
            _, carry, gap, parameter = word_data(word)
            check(parameter == 1, "unit-gap denominator")
            unit_rows.append({"word": word, "K": k, "L": length, "B": carry,
                              "Delta": gap, "b1_cycle": reconstruct(word, 1)})

    repetition_rows = []
    for word in [(1,), (2,), (1, 2), (1, 3), (1, 1, 1, 2, 1, 1, 4)]:
        k, carry, gap, parameter = word_data(word)
        for repeat in range(2, 6):
            rk, rb, rd, rq = word_data(word * repeat)
            factor = sum(2**(k * (repeat - 1 - j)) * 3**(len(word) * j) for j in range(repeat))
            check((rk, rb, rd, rq) == (repeat * k, factor * carry, factor * gap, parameter),
                  "repeated word factorization")
            check(reconstruct(word * repeat, parameter) == reconstruct(word, parameter),
                  "repeat changes cycle")
        repetition_rows.append({"word": word, "q": parameter, "repeats_checked": [2, 3, 4, 5]})

    clock_rows = []
    for word in compositions(4, 2):
        _, carry, gap, parameter = word_data(word)
        clock_rows.append({"word": word, "B": carry, "Delta": gap, "q": parameter,
                           "rational_source_b1": str(Fraction(carry, gap))})
    check([row["rational_source_b1"] for row in clock_rows] == ["5/7", "1", "11/7"], "same-clock control")
    rational = Fraction(5, 7)
    rational_nodes = [rational]
    for k in (1, 3):
        rational = (3 * rational + 1) / 2**k
        check(rational.numerator % 2 and rational.denominator % 2, "rational odd valuation guard")
        rational_nodes.append(rational)
    check(rational_nodes == [Fraction(5, 7), Fraction(11, 7), Fraction(5, 7)], "rational cycle")

    sign_scaling_checks = 0
    for b in range(-9, 10, 2):
        for n in range(-101, 102, 2):
            if 3 * n + b == 0:
                continue
            value, k = step(n, b)
            check(step(-n, -b) == (-value, k), "sign conjugacy")
            for d in (-5, -3, -1, 1, 3, 5):
                check(step(d * n, d * b) == (d * value, k), "odd scaling")
                sign_scaling_checks += 1
    check(abs(step(1, 1)[0]) == abs(step(-1, 1)[0]) == 1, "magnitude boundary")
    check((abs(step(3, 1)[0]), abs(step(-3, 1)[0])) == (5, 1), "minimal sign-forgetting hostile")

    return {
        "status": "FINITE-EXACT independent replay; no claim about longer cycles or convergence",
        "universe": {"parameters": [-1, 1], "least_odd_period_max": 10,
                     "word_lengths": [1, 10], "total_exponent_bounds": "L <= K <= 2L",
                     "node_height_cap": None, "ordered_words": total_words,
                     "integer_words_each_parameter": accepted_words,
                     "distinct_cycles_each_parameter": len(cycles)},
        "by_length": by_length,
        "b1_cycles": cycle_rows,
        "bminus1_cycles": sorted(negative_parameter_cycles, key=lambda c: (len(c), abs(c[0]), c)),
        "unit_gap_controls": {"finite_clock_box": "1 <= L <= K <= 64", "rows": unit_rows},
        "same_clock_carry_hostile": clock_rows,
        "rational_cycle": [str(n) for n in rational_nodes],
        "repetition_controls": repetition_rows,
        "other_controls": {"rotation_checks": rotation_checks, "parameter_checks": parameter_checks,
                           "parameter_word_universe": "1 <= L <= 6, L <= K <= 2L, odd -15 <= b <= 15",
                           "sign_scaling_checks": sign_scaling_checks,
                           "sign_scaling_universe": "odd -9 <= b <= 9; odd -101 <= n <= 101; 3n+b != 0; d in +/-{1,3,5}",
                           "minimal_sign_forgetting_hostile": {"magnitude": 3, "positive_image_magnitude": 5,
                                                               "negative_image_magnitude": 1}},
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(main(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
