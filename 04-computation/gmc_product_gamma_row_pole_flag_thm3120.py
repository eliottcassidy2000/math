#!/usr/bin/env python3
"""Exact pole-prefix/Newton-flag companion for THM-3120.

For each fixed integer support ``0<a<b`` in the stated 115-support bank, this
script forms the complete-homogeneous row generating function

    F_i(t) = sum_R c_R prod_(r in S_R) (1-r*t)^(-1)

from the exact THM-3110 residual alphabets.  It clears the least common pole
denominator, cancels every common linear factor exactly, and writes the
reduced numerator as ``t^5 P(t)``.  If the remaining poles are
``r_1 >= ... >= r_E`` and ``d=deg(P)``, it then certifies the unique flag

    P(t) = sum_(k=0)^d c_k t^k prod_(ell=1)^(d-k) (1-r_ell*t)

with every ``c_k`` strictly positive.  Dividing by the denominator gives a
literal positive sum of suffix-denominator series and proves every row
coefficient in every degree ``n>=5`` for each certified support.

The reciprocal polynomial is independently checked in its Newton form,
including confluent divided differences at repeated pole nodes.  The script
also records the exceptional top-pole cancellation family and exact hostile
controls against stronger but false coefficientwise mechanisms.

Only Python integers and exact fractions are used.
"""

from collections import Counter
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from math import comb
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


COMPANION = Path(__file__).with_name(
    "gmc_product_gamma_arbitrary_anchored_schur_thm3110.py"
)
SPEC = spec_from_file_location("thm3110_companion", COMPANION)
g = module_from_spec(SPEC)
SPEC.loader.exec_module(g)


def support_universe():
    """The same 115 exact supports as the THM-3110 Young-cover scout."""

    return tuple(
        (a, b)
        for a in range(1, 11)
        for b in range(a + 1, min(3 * a + 4, 21) + 1)
    )


def trim(poly):
    out = list(poly)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def add_scaled(target, source, scalar):
    if len(target) < len(source):
        target.extend([0] * (len(source) - len(target)))
    for index, coefficient in enumerate(source):
        target[index] += scalar * coefficient
    return trim(target)


def multiply(left, right):
    out = [0] * (len(left) + len(right) - 1)
    for i, coefficient in enumerate(left):
        if coefficient:
            for j, other in enumerate(right):
                out[i + j] += coefficient * other
    return trim(out)


def multiply_linear(poly, root):
    """Multiply a low-to-high coefficient list by 1-root*t."""

    out = [0] * (len(poly) + 1)
    for index, coefficient in enumerate(poly):
        out[index] += coefficient
        out[index + 1] -= root * coefficient
    return trim(out)


def product_from_counter(counter):
    out = [1]
    for root in sorted(counter):
        for _ in range(counter[root]):
            out = multiply_linear(out, root)
    return out


def divide_linear_if_exact(poly, root):
    """Return the quotient by 1-root*t, or None when it is not a factor."""

    if len(poly) <= 1:
        return None
    quotient = [0] * (len(poly) - 1)
    quotient[0] = poly[0]
    for index in range(1, len(quotient)):
        quotient[index] = poly[index] + root * quotient[index - 1]
    if poly[-1] != -root * quotient[-1]:
        return None
    return trim(quotient)


def atom_counters(invariant, a, b):
    return tuple(
        Counter(root for root in g.residual_roots(invariant, rows, a, b)
                if root)
        for _, rows in g.BANKS[invariant]
    )


def reduced_row_fraction(invariant, a, b):
    """Return reduced N,D and exact pole metadata for one row OGF."""

    counters = atom_counters(invariant, a, b)
    envelope = Counter()
    for counter in counters:
        for root, multiplicity in counter.items():
            envelope[root] = max(envelope[root], multiplicity)

    numerator = [0]
    for (coefficient, _), counter in zip(g.BANKS[invariant], counters):
        missing = envelope.copy()
        missing.subtract(counter)
        require(all(value >= 0 for value in missing.values()),
                "atom denominator escaped its envelope")
        numerator = add_scaled(
            numerator,
            product_from_counter(+missing),
            coefficient,
        )
    numerator = trim(numerator)
    require(numerator != [0], "row numerator vanished")

    remaining = envelope.copy()
    cancelled = Counter()
    for root in sorted(envelope, reverse=True):
        while remaining[root]:
            quotient = divide_linear_if_exact(numerator, root)
            if quotient is None:
                break
            numerator = quotient
            remaining[root] -= 1
            cancelled[root] += 1
        if remaining[root] == 0:
            del remaining[root]

    denominator = product_from_counter(remaining)
    # Reconstructing with the cancelled factors is an independent guard on
    # exact division and on the least-common-denominator construction.
    cancelled_poly = product_from_counter(cancelled)
    raw_numerator = multiply(numerator, cancelled_poly)
    raw_denominator = multiply(denominator, cancelled_poly)
    require(raw_denominator == product_from_counter(envelope),
            "denominator cancellation reconstruction failed")

    direct_raw = [0]
    for (coefficient, _), counter in zip(g.BANKS[invariant], counters):
        missing = envelope.copy()
        missing.subtract(counter)
        direct_raw = add_scaled(
            direct_raw, product_from_counter(+missing), coefficient)
    require(trim(direct_raw) == raw_numerator,
            "numerator cancellation reconstruction failed")
    require(all(divide_linear_if_exact(numerator, root) is None
                for root in remaining),
            "reduced numerator and denominator retain a common pole factor")
    return numerator, denominator, remaining, cancelled, envelope, counters


def reciprocal(poly):
    """Coefficient list of x^d poly(1/x)."""

    return list(reversed(poly))


def derivative_over_factorial_at(poly, order, point):
    return sum(
        coefficient * comb(degree, order) * point ** (degree - order)
        for degree, coefficient in enumerate(poly)
        if degree >= order
    )


def confluent_newton_coefficients(poly, nodes):
    """Newton coefficients on consecutive, possibly repeated, nodes."""

    degree = len(poly) - 1
    require(len(nodes) >= degree + 1, "insufficient pole nodes for flag")
    nodes = tuple(nodes[:degree + 1])
    table = [[Fraction(0) for _ in range(degree + 1)]
             for _ in range(degree + 1)]
    for index, node in enumerate(nodes):
        table[index][0] = Fraction(sum(
            coefficient * node ** power
            for power, coefficient in enumerate(poly)
        ))
    for order in range(1, degree + 1):
        for start in range(degree + 1 - order):
            if nodes[start] == nodes[start + order]:
                # Sorted nodes make the whole block confluent.
                table[start][order] = Fraction(
                    derivative_over_factorial_at(
                        poly, order, nodes[start]))
            else:
                table[start][order] = Fraction(
                    table[start + 1][order - 1]
                    - table[start][order - 1],
                    nodes[start + order] - nodes[start],
                )
    return tuple(table[0][order] for order in range(degree + 1))


def pole_flag(poly, roots):
    """Return c_k in P=sum c_k t^k q_(d-k), checking Newton duality."""

    degree = len(poly) - 1
    require(len(roots) >= degree + 1,
            "pole count is too small for the reciprocal Newton flag")
    prefixes = [[1]]
    for root in roots[:degree]:
        prefixes.append(multiply_linear(prefixes[-1], root))

    residual = list(poly)
    coefficients = []
    for k in range(degree + 1):
        coefficient = residual[k]
        coefficients.append(coefficient)
        basis = [0] * k + prefixes[degree - k]
        residual = add_scaled(residual, basis, -coefficient)
        if len(residual) < degree + 1:
            residual.extend([0] * (degree + 1 - len(residual)))
    require(all(value == 0 for value in residual),
            "pole-prefix triangular reconstruction failed")

    newton = confluent_newton_coefficients(reciprocal(poly), roots)
    require(newton == tuple(reversed(coefficients)),
            "reciprocal confluent-Newton coefficients disagree with flag")
    return tuple(coefficients), tuple(tuple(poly) for poly in prefixes)


def series_coefficients(numerator, denominator, maximum_degree):
    require(denominator[0] == 1, "series denominator is not normalized")
    out = [0] * (maximum_degree + 1)
    for degree in range(maximum_degree + 1):
        value = numerator[degree] if degree < len(numerator) else 0
        for lag in range(1, min(degree, len(denominator) - 1) + 1):
            value -= denominator[lag] * out[degree - lag]
        out[degree] = value
    return tuple(out)


def cycle_123_coefficient(invariant, a, b):
    """Coefficient of z1*z2*z3 in sum c exp(sum p_j(S)z_j)."""

    total = 0
    for coefficient, rows in g.BANKS[invariant]:
        roots = g.residual_roots(invariant, rows, a, b)
        powers = [sum(root ** degree for root in roots)
                  for degree in (1, 2, 3)]
        total += coefficient * powers[0] * powers[1] * powers[2]
    return total


def local_top_constant_ratios(invariant, k):
    """Local numerator constants at t=1/(4k) when (a,b)=(2k,3k)."""

    a, b, pole = 2 * k, 3 * k, 4 * k
    counters = atom_counters(invariant, a, b)
    envelope = Counter()
    for counter in counters:
        for root, multiplicity in counter.items():
            envelope[root] = max(envelope[root], multiplicity)
    maximum = envelope[pole]
    records = []
    for index, ((coefficient, _), counter) in enumerate(
            zip(g.BANKS[invariant], counters)):
        if counter[pole] != maximum:
            continue
        missing = envelope.copy()
        missing.subtract(counter)
        local = Fraction(1)
        for root, multiplicity in (+missing).items():
            if root != pole:
                local *= Fraction(pole - root, pole) ** multiplicity
        records.append((index, coefficient, local))
    base = records[0][2]
    return tuple((index, coefficient, local / base)
                 for index, coefficient, local in records)


def cancellation_audit():
    """Check the exact 4|a, b=3a/2 top-pole cancellation family."""

    expected_indices = ((20, 21), (19, 21, 22))
    family = []
    for k in range(1, 11):
        expected_ratios = (
            ((20, 2, Fraction(1)),
             (21, -4, Fraction((-1) ** k, 2))),
            ((19, 1, Fraction(1)),
             (21, -4, Fraction((-1) ** k, 2)),
             (22, 4, Fraction(1, 4))),
        )
        for invariant in (0, 1):
            ratios = local_top_constant_ratios(invariant, k)
            require(tuple(index for index, _, _ in ratios)
                    == expected_indices[invariant],
                    "top-pole contributing-atom set drift")
            require(ratios == expected_ratios[invariant],
                    "top-pole local-constant ratio drift")
            residue = sum(coefficient * ratio
                          for _, coefficient, ratio in ratios)
            require(residue == 2 - 2 * (-1) ** k,
                    "top-pole residue formula drift")
            _, _, _, cancelled, _, _ = reduced_row_fraction(
                invariant, 2 * k, 3 * k)
            require((cancelled[4 * k] > 0) == (k % 2 == 0),
                    "top-pole cancellation parity drift")
            if cancelled[4 * k]:
                family.append((invariant + 1, 2 * k, 3 * k, 4 * k))
    require(family == [
        (invariant, a, b, 2 * a)
        for a, b in ((4, 6), (8, 12), (12, 18), (16, 24), (20, 30))
        for invariant in (1, 2)
    ], "exceptional cancellation-family census drift")
    return tuple(family)


def exact_census():
    supports = support_universe()
    require(len(supports) == 115, "support-universe count drift")

    stats = [
        {"cases": 0, "flags": 0, "degree_min": None, "degree_max": None,
         "order_min": None, "order_max": None, "gap_min": None,
         "minimum": None, "minimum_record": None}
        for _ in range(2)
    ]
    cancellation_records = []
    all_positive = 0

    for invariant in (0, 1):
        for a, b in supports:
            numerator, denominator, remaining, cancelled, _, _ = (
                reduced_row_fraction(invariant, a, b))
            valuation = next(index for index, coefficient
                             in enumerate(numerator) if coefficient)
            require(valuation == 5, "row numerator valuation drift")
            poly = numerator[5:]
            roots = sorted(remaining.elements(), reverse=True)
            degree = len(poly) - 1
            order = len(roots)
            coefficients, prefixes = pole_flag(poly, roots)
            require(all(coefficient > 0 for coefficient in coefficients),
                    "nonpositive pole-flag coefficient")

            # Exact positive rational decomposition, recombined over the
            # common reduced denominator.
            reconstructed = [0]
            for k, coefficient in enumerate(coefficients):
                term = [0] * (5 + k) + list(prefixes[degree - k])
                reconstructed = add_scaled(reconstructed, term, coefficient)
            require(trim(reconstructed) == numerator,
                    "positive suffix-denominator decomposition drift")
            require(order >= degree + 1,
                    "strict-tail denominator vanished")

            # The reduced rational function has minimal recurrence order E.
            require(len(denominator) - 1 == order,
                    "denominator order/root count drift")
            require(all(divide_linear_if_exact(numerator, root) is None
                        for root in remaining),
                    "minimal recurrence gcd guard failed")

            for root, multiplicity in cancelled.items():
                cancellation_records.append(
                    (invariant + 1, a, b, root, multiplicity))

            record = stats[invariant]
            record["cases"] += 1
            record["flags"] += len(coefficients)
            all_positive += len(coefficients)
            for key, value in (("degree_min", degree),
                               ("order_min", order),
                               ("gap_min", order - degree)):
                if record[key] is None or value < record[key]:
                    record[key] = value
            for key, value in (("degree_max", degree),
                               ("order_max", order)):
                if record[key] is None or value > record[key]:
                    record[key] = value
            for k, coefficient in enumerate(coefficients):
                candidate = (coefficient, a, b, degree, k)
                if (record["minimum"] is None
                        or candidate < record["minimum_record"]):
                    record["minimum"] = coefficient
                    record["minimum_record"] = candidate

    # In the canonical finite universe these are the only denominator
    # cancellations.  The wider all-k family is checked separately above.
    require(cancellation_records == [
        (1, 4, 6, 8, 1),
        (1, 8, 12, 16, 1),
        (2, 4, 6, 8, 1),
        (2, 8, 12, 16, 1),
    ], "115-support cancellation census drift")
    expected_stats = (
        (115, 3953, 2, 67, 11, 132, 9,
         (36, 1, 2, 2, 0)),
        (115, 4288, 2, 68, 11, 133, 9,
         (32, 1, 2, 2, 0)),
    )
    for record, expected in zip(stats, expected_stats):
        actual = (
            record["cases"], record["flags"],
            record["degree_min"], record["degree_max"],
            record["order_min"], record["order_max"],
            record["gap_min"], record["minimum_record"],
        )
        require(actual == expected, "pole-flag census drift")
    require(all_positive == 8241, "total pole-flag count drift")
    return supports, tuple(stats), all_positive, tuple(cancellation_records)


def controls():
    examples = []
    for invariant, expected_poly, expected_flag in (
        (0, [36, -108, -72], (36, 216, 288)),
        (1, [32, 24, 40], (32, 312, 960)),
    ):
        numerator, denominator, remaining, _, _, counters = (
            reduced_row_fraction(invariant, 1, 2))
        require(numerator[:5] == [0] * 5, "(1,2) valuation drift")
        require(numerator[5:] == expected_poly, "(1,2) numerator drift")
        roots = sorted(remaining.elements(), reverse=True)
        flag, _ = pole_flag(expected_poly, roots)
        require(flag == expected_flag, "(1,2) flag drift")
        require(denominator == product_from_counter(Counter({
            1: 4, 2: 3, 3: 2, 4: 1, 5: 1,
        })), "(1,2) denominator drift")
        absent = sum(counter[5] == 0 for counter in counters)
        require(absent == (21, 21)[invariant],
                "top virtual-pole absence count drift")
        examples.append((invariant + 1, expected_poly, expected_flag, absent))

    require(cycle_123_coefficient(0, 1, 2) == -84,
            "I1 independent-cycle hostile drift")
    require(cycle_123_coefficient(1, 1, 2) == -84,
            "I2 independent-cycle hostile drift")

    numerator, denominator, _, _, _, _ = reduced_row_fraction(0, 1, 2)
    twice_stripped = multiply(
        numerator, multiply_linear(multiply_linear([1], 5), 4))
    hostile_series = series_coefficients(twice_stripped, denominator, 17)
    require(hostile_series[17] == -5901696,
            "sequential pole-stripping hostile drift")
    return tuple(examples), hostile_series[17]


def main():
    supports, stats, positive, finite_cancellations = exact_census()
    family = cancellation_audit()
    examples, strip_hostile = controls()

    print(f"support_universe={len(supports)}:"
          "1<=a<=10:a<b<=min(3a+4,21):banks=2:cases=230")
    for invariant, record in enumerate(stats, 1):
        print(
            f"I{invariant}_flag=cases:{record['cases']}:"
            f"coefficients:{record['flags']}:positive:{record['flags']}:"
            "zero:0:negative:0:"
            f"degree:{record['degree_min']}..{record['degree_max']}:"
            f"recurrence_order:{record['order_min']}..{record['order_max']}:"
            f"order_minus_degree_min:{record['gap_min']}:"
            f"minimum:{record['minimum_record']}"
        )
    print(f"all_degree_row_positivity=n>=5:certificates={positive}:PASS")
    print("newton_chain=ordinary_and_confluent_divided_differences:PASS")
    print("positive_rational_decomposition=denominator_suffixes:PASS")
    print("minimal_recurrence=reduced_split_denominator:PASS")
    print("finite_cancellations=" + ";".join(
        f"I{i}:{a}-{b}:root{root}:multiplicity{multiplicity}"
        for i, a, b, root, multiplicity in finite_cancellations))
    print("cancellation_family=" + ";".join(
        f"I{i}:{a}-{b}:root{root}" for i, a, b, root in family))
    print("example_1_2=" + ";".join(
        f"I{i}:P={poly}:flag={flag}:top_pole_absent_atoms={absent}"
        for i, poly, flag, absent in examples))
    print("virtual_pole_commutator_data=top5_absent:I1=21/24:I2=21/25")
    print("hostile_independent_cycle_z1z2z3=I1:-84:I2:-84")
    print(f"hostile_second_prefix_strip=I1:1-2:n17:{strip_hostile}")
    print("status=PASS")


if __name__ == "__main__":
    main()
