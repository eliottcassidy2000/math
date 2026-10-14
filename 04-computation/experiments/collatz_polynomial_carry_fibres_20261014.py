#!/usr/bin/env python3
"""Exact slope-admissible carry frontier through L=16; no integer-orbit census."""
from collections import Counter
from fractions import Fraction
from math import ceil, log2
import hashlib
import json

checks = 0


def require(ok, message):
    global checks
    checks += 1
    if not ok:
        raise ArithmeticError(message)


def decode(L, A, B):
    word = []
    for length in range(L, 1, -1):
        remainder = B - 3 ** (length - 1)
        require(remainder > 0, "decoder positive remainder")
        a = (remainder & -remainder).bit_length() - 1
        require(a >= 1, "decoder positive first valuation")
        word.append(a)
        A, B = A - a, remainder // (1 << a)
    require(B == 1 and A >= 1, "decoder terminal letter")
    return tuple(word + [A])


# For each carry and cost, retain its exact multiplicity and one exact B witness.
frontier, minimum, cost_counts = {0: {0: (1, 0)}}, {0: 0}, Counter({0: 1})
rows = []
for L in range(1, 17):
    q, budget = 3 ** L, (3 ** L).bit_length() - 1
    inv = [0] + [pow(1 << a, -1, q) for a in range(1, budget - L + 2)]
    next_frontier = {}
    for old_r, costs in frontier.items():
        for A, (count, B) in costs.items():
            Bnext = 3 * B + (1 << A)
            for a in range(1, budget - A + 1):
                r, total = (3 * old_r + 1) * inv[a] % q, A + a
                cell = next_frontier.setdefault(r, {})
                if total in cell:
                    old_count, old_B = cell[total]
                    cell[total] = old_count + count, min(old_B, Bnext)
                else:
                    cell[total] = count, Bnext
    frontier = next_frontier

    # A separate minimum-cost compiler, and an independent cost-only word counter.
    next_minimum, next_counts = {}, Counter()
    for old_r, A in minimum.items():
        for a in range(1, budget - A + 1):
            r, total = (3 * old_r + 1) * inv[a] % q, A + a
            next_minimum[r] = min(next_minimum.get(r, budget + 1), total)
    for A, count in cost_counts.items():
        for a in range(1, budget - A + 1):
            next_counts[A + a] += count
    minimum, cost_counts = next_minimum, next_counts
    require(minimum == {r: min(costs) for r, costs in frontier.items()}, "minimum-cost compiler")
    actual_cost_counts = Counter()
    for costs in frontier.values():
        for A, (count, _) in costs.items():
            actual_cost_counts[A] += count
    require(actual_cost_counts == cost_counts, "independent cost-only counter")

    # B witnesses and their exact inverse decoder: a bounded sample per level.
    for r, costs in list(frontier.items())[:16]:
        for A, (_, B) in costs.items():
            require(0 < 3 * B <= L * q, "polynomial fibre height bound")
            require(B * pow(1 << A, -1, q) % q == r, "carry witness")
            word = decode(L, A, B)
            running, BB = 0, 0
            for i, a in enumerate(word, 1):
                BB, running = 3 * BB + (1 << running), running + a
                require(1 << running < 3 ** i, "decoded prefix slope budget")
            require(running == A and BB == B, "decoder round trip")

    histogram = Counter(sum(count for count, _ in costs.values()) for costs in frontier.values())
    W, R = sum(cost_counts.values()), len(frontier)
    bound = (budget - L + 1) * ceil(L / 3)
    require(max(histogram) <= bound, "all finite fibres within proved polynomial bound")
    total_mass_units = sum(count * (1 << (budget - A)) for A, count in cost_counts.items())
    minimum_mass_units = sum(1 << (budget - A) for A in minimum.values())
    require(minimum_mass_units <= total_mass_units <= bound * minimum_mass_units,
            "weighted polynomial comparison")

    # Descriptive finite entropy only; exact theorem supplies the asymptotic result.
    fibre_mass_log_sum = 0.0
    for costs in frontier.values():
        units = sum(count * (1 << (budget - A)) for A, (count, _) in costs.items())
        fibre_mass_log_sum += units * log2(units)
    word_log_weight_sum = sum(count * (1 << (budget - A)) * (budget - A)
                              for A, count in cost_counts.items())
    conditional_bits = (fibre_mass_log_sum - word_log_weight_sum) / total_mass_units
    rows.append({"length": L, "budget": budget, "budget_jump": budget - ((3 ** (L - 1)).bit_length() - 1),
                 "word_count": W, "carry_count": R, "carry_word_ratio": str(Fraction(R, W)),
                 "maximum_fibre": max(histogram), "fibre_histogram": dict(sorted(histogram.items())),
                 "maximum_fixed_cost_multiplicity_observed": max(count for c in frontier.values() for count, _ in c.values()),
                 "distinct_carry_cost_pairs": sum(map(len, frontier.values())),
                 "proved_fibre_bound": bound,
                 "source_mass": str(Fraction(total_mass_units, 1 << budget)),
                 "minimum_representative_mass": str(Fraction(minimum_mass_units, 1 << budget)),
                 "source_weighted_word_uncertainty_given_carry_bits_approx": round(conditional_bits, 12)})

alpha = log2(3)
beta = alpha ** alpha / (alpha - 1) ** (alpha - 1)
result = {"status": "FINITE-EXACT frontier counts; polynomial-fibre and rate theorems proved in note",
          "checks": checks, "maximum_length": 16, "rows": rows,
          "beta_approx": beta, "weighted_rate_beta_over_3_approx": beta / 3,
          "conclusion": "The strict unweighted carry-entropy gap is refuted by a general polynomial fibre bound.",
          "limitations": ["Fixed-cost injectivity is observed only through length 16, not used as a theorem.",
                         "These are slope-admissible native words, not all finite-source-only no-descent prefixes.",
                         "No claim about convergence of one fixed positive integer orbit."]}
payload = json.dumps(result, sort_keys=True, separators=(",", ":")).encode()
result["payload_sha256"] = hashlib.sha256(payload).hexdigest()
print(json.dumps(result, sort_keys=True, indent=2))
