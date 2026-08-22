#!/usr/bin/env python3
"""Independent audit of THM-3166.

This file does not import the main subset-DP engine.  It obtains path-cover
profiles from permutation backward-adjacency cuts, checks joins directly, and
audits the falling-factorial product in the ordinary power basis.
"""

from __future__ import annotations

import argparse
import hashlib
from itertools import permutations
from math import comb, factorial
from pathlib import Path


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "00-navigation").is_dir())
MAIN = ROOT / "04-computation/tournament_order_join_falling_factorial_transform_thm3166.py"
MAIN_OUTPUT = ROOT / "05-knowledge/results/tournament_order_join_falling_factorial_transform_thm3166.out"
OUTPUT = ROOT / "05-knowledge/results/tournament_order_join_falling_factorial_independent_audit_thm3166.out"

EXPECTED_MAIN_SHA256 = "5d38326311fc8fea8976f105348b131a106b7d0f2b912b50937f969152d310fd"
EXPECTED_MAIN_OUTPUT_SHA256 = "55f31dbf60b23184d69167a26d39fad192a3230bb2163cc046c8b12df3885eff"


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def edge_index(n: int, i: int, j: int) -> int:
    return i * (2 * n - i - 1) // 2 + (j - i - 1)


def beats(n: int, code: int, i: int, j: int) -> bool:
    if i < j:
        return bool((code >> edge_index(n, i, j)) & 1)
    return not bool((code >> edge_index(n, j, i)) & 1)


def cut_profile(n: int, code: int):
    """Permutation-cut formula, independent of subset path enumeration."""
    backward = [0] * n
    for order in permutations(range(n)):
        b = sum(not beats(n, code, order[i], order[i + 1]) for i in range(n - 1))
        backward[b] += 1
    profile = [0] * (n + 1)
    for d in range(1, n + 1):
        ordered = sum(
            comb(n - 1 - b, d - 1 - b) * backward[b]
            for b in range(d)
        )
        require(ordered % factorial(d) == 0, (n, code, d, ordered))
        profile[d] = ordered // factorial(d)
    return tuple(profile)


def backward_distribution(n: int, code: int):
    counts = [0] * n
    for order in permutations(range(n)):
        backward = sum(
            not beats(n, code, order[i], order[i + 1])
            for i in range(n - 1)
        )
        counts[backward] += 1
    return tuple(counts)


def negative_binomial_sum(backward, m: int) -> int:
    n = len(backward)
    return sum(count * comb(m + b, n) for b, count in enumerate(backward))


def recover_backward_from_negative(profile):
    n = len(profile) - 1
    recovered = [0] * n
    for m in range(1, n + 1):
        b = n - m
        positive = (-1) ** n * qvalue(profile, -m)
        recovered[b] = positive - sum(
            recovered[c] * comb(m + c, n)
            for c in range(b + 1, n)
        )
    return tuple(recovered)


def join_code(n: int, first: int, m: int, second: int) -> int:
    total = n + m
    out = 0
    for i in range(total):
        for j in range(i + 1, total):
            value = (
                beats(n, first, i, j)
                if j < n
                else beats(m, second, i - n, j - n)
                if i >= n
                else True
            )
            if value:
                out |= 1 << edge_index(total, i, j)
    return out


def convolution(first, second):
    out = [0] * (len(first) + len(second) - 1)
    for a, first_count in enumerate(first):
        for b, second_count in enumerate(second):
            for k in range(min(a, b) + 1):
                out[a + b - k] += (
                    first_count * second_count * comb(a, k) * comb(b, k) * factorial(k)
                )
    return tuple(out)


def falling_value(t: int, d: int) -> int:
    out = 1
    for j in range(d):
        out *= t - j
    return out


def qvalue(profile, t: int) -> int:
    return sum(count * falling_value(t, d) for d, count in enumerate(profile))


def inverse_from_values(values):
    profile = []
    for d in range(len(values)):
        numerator = sum(
            (-1) ** (d - j) * comb(d, j) * values[j]
            for j in range(d + 1)
        )
        require(numerator % factorial(d) == 0, (d, numerator))
        profile.append(numerator // factorial(d))
    return tuple(profile)


def polynomial_product(first, second):
    out = [0] * (len(first) + len(second) - 1)
    for i, a in enumerate(first):
        for j, b in enumerate(second):
            out[i + j] += a * b
    return tuple(out)


def falling_polynomial(d: int):
    out = (1,)
    for root in range(d):
        out = polynomial_product(out, (-root, 1))
    return out


def polynomial_add_scaled(target, source, scale: int):
    if len(target) < len(source):
        target.extend([0] * (len(source) - len(target)))
    for i, value in enumerate(source):
        target[i] += scale * value


def falling_product_rhs(a: int, b: int):
    out = []
    for k in range(min(a, b) + 1):
        polynomial_add_scaled(
            out,
            falling_polynomial(a + b - k),
            comb(a, k) * comb(b, k) * factorial(k),
        )
    return tuple(out)


def score_sequence(n: int, code: int):
    return tuple(sorted(sum(beats(n, code, i, j) for j in range(n) if i != j)
                        for i in range(n)))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    require(sha256(MAIN) == EXPECTED_MAIN_SHA256,
            ("main changed", sha256(MAIN), EXPECTED_MAIN_SHA256))
    require(sha256(MAIN_OUTPUT) == EXPECTED_MAIN_OUTPUT_SHA256,
            ("main output changed", sha256(MAIN_OUTPUT), EXPECTED_MAIN_OUTPUT_SHA256))
    main_output = MAIN_OUTPUT.read_text(encoding="utf-8")
    for token in ("labelled_tournament_controls=1099",
                  "negative_colour_reciprocity_controls=7693",
                  "negative_colour_inverse_controls=1099",
                  "ordered_join_controls=5625",
                  "repeated_join_closed_form_controls=4275", "all_exact_controls=PASS"):
        require(token in main_output, ("missing main token", token))

    profiles = {}
    tournament_checks = 0
    negative_reciprocity_checks = 0
    negative_inverse_checks = 0
    for n in range(1, 6):
        for code in range(1 << comb(n, 2)):
            profile = cut_profile(n, code)
            backward = backward_distribution(n, code)
            values = tuple(qvalue(profile, j) for j in range(n + 1))
            require(inverse_from_values(values) == profile, ("inverse", n, code))
            for m in range(1, 8):
                require(
                    (-1) ** n * qvalue(profile, -m)
                    == negative_binomial_sum(backward, m),
                    ("negative reciprocity", n, code, m),
                )
                negative_reciprocity_checks += 1
            require(recover_backward_from_negative(profile) == backward,
                    ("negative inverse", n, code))
            negative_inverse_checks += 1
            profiles[n, code] = profile
            tournament_checks += 1
    require(tournament_checks == 1099, tournament_checks)
    require(negative_reciprocity_checks == 7693, negative_reciprocity_checks)
    require(negative_inverse_checks == 1099, negative_inverse_checks)

    # Complete direct joins of total order at most six, plus a hostile 4+4 grid.
    join_pairs = []
    for n in range(1, 5):
        for first in range(1 << comb(n, 2)):
            for m in range(1, 5):
                if n + m > 6:
                    continue
                for second in range(1 << comb(m, 2)):
                    join_pairs.append((n, first, m, second))
    selected_four = tuple(range(0, 64, 4))
    join_pairs.extend((4, first, 4, second)
                      for first in selected_four for second in selected_four)
    require(len(join_pairs) == 761, len(join_pairs))

    for n, first, m, second in join_pairs:
        left = profiles[n, first]
        right = profiles[m, second]
        direct = cut_profile(n + m, join_code(n, first, m, second))
        require(direct == convolution(left, right),
                ("direct convolution", n, first, m, second))
        values = tuple(qvalue(direct, j) for j in range(n + m + 1))
        product_values = tuple(qvalue(left, j) * qvalue(right, j)
                               for j in range(n + m + 1))
        require(values == product_values, ("transform product", n, first, m, second))
        require(inverse_from_values(values) == direct, ("direct inverse", n, m))

    basis_checks = 0
    for a in range(1, 33):
        for b in range(a, 33):
            require(polynomial_product(falling_polynomial(a), falling_polynomial(b)) ==
                    falling_product_rhs(a, b), ("falling product", a, b))
            basis_checks += 1
    require(basis_checks == 528, basis_checks)

    # Fixed-depth exponential formulas, reconstructed only from Q(0),...,Q(d).
    repeated_checks = 0
    seeds = ((1, 0, 8), (2, 0, 4), (2, 1, 4), (3, 5, 2), (3, 7, 2))
    for n, code, max_repetitions in seeds:
        profile = profiles[n, code]
        running = (1,)
        for repetitions in range(1, max_repetitions + 1):
            running = profile if repetitions == 1 else convolution(running, profile)
            for d in range(1, len(running)):
                numerator = sum(
                    (-1) ** (d - j) * comb(d, j) * qvalue(profile, j) ** repetitions
                    for j in range(d + 1)
                )
                require(numerator == factorial(d) * running[d],
                        ("exponential formula", n, code, repetitions, d))
                repeated_checks += 1

    k1 = profiles[1, 0]
    c3 = profiles[3, 5]
    transitive3 = profiles[3, 7]
    source_code = join_code(1, 0, 3, 5)
    sink_code = join_code(3, 5, 1, 0)
    source_profile = cut_profile(4, source_code)
    sink_profile = cut_profile(4, sink_code)
    require(source_profile == sink_profile == convolution(k1, c3),
            (source_profile, sink_profile))
    source_scores = score_sequence(4, source_code)
    sink_scores = score_sequence(4, sink_code)
    require(source_scores == (1, 1, 1, 3) and sink_scores == (0, 2, 2, 2),
            (source_scores, sink_scores))
    require(tuple(qvalue(transitive3, j) for j in range(4)) == (0, 1, 8, 27),
            transitive3)
    require(tuple(qvalue(c3, j) for j in range(4)) == (0, 3, 12, 33), c3)

    lines = [
        "TOURNAMENT ORDER-JOIN FALLING-FACTORIAL TRANSFORM -- THM-3166 INDEPENDENT AUDIT",
        "status=INDEPENDENTLY HOSTILE-AUDITED;permutation-cut profile engine;no main import",
        f"permutation_profile_controls={tournament_checks};all labelled tournaments orders 1..5",
        f"negative_colour_reciprocity_controls={negative_reciprocity_checks};m=1..7;direct_backward_distribution=PASS",
        f"negative_colour_inverse_controls={negative_inverse_checks};triangular_recovery=PASS",
        f"direct_join_controls={len(join_pairs)};all total-order<=6 plus 256 selected order-8 joins",
        f"falling_basis_controls={basis_checks};ordinary-power identity through a<=b<=32",
        f"fixed_depth_exponential_controls={repeated_checks};five seeds;direct convolution=PASS",
        f"order_loss_hostile=profiles_equal:{source_profile};scores={source_scores}!={sink_scores}",
        "cyclic_hostile=Q_T3(j)=(0,1,8,27);Q_C3(j)=(0,3,12,33);product law restricted to order-join",
        "main_token_audit=1099 tournaments;7693 negative values;1099 negative inverses;5625 joins;4275 repeated controls;PASS",
        f"main_sha256={EXPECTED_MAIN_SHA256}",
        f"main_output_sha256={EXPECTED_MAIN_OUTPUT_SHA256}",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"source_sha256={sha256(HERE)}",
        "all_independent_controls=PASS",
    ]
    text = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
