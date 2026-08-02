#!/usr/bin/env python3
"""Exact N>=8 scout for the THM-3115 refinement-flow mechanism.

This is deliberately a theorem-ID-free continuation of the N=5,6,7
certificate.  It uses power sums and the monomial-basis recurrence instead of
enumerating residual roots one variable at a time.  For both product-Gamma
banks and both integer chambers it (at N=8 by default):

* reconstructs every chamber-Newton coefficient of

      G_mu=(Phi(h_8)m_mu(Q)-Phi(m_mu)h_8(Q))/base;

* checks the predicted total-degree bound and three off-grid values;
* transports every negative coefficient to positive coarsening types; and
* tests whether hooks, two-row types, or their union already suffice.

Only exact integers and fractions are used.
"""

import ast
import hashlib
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from math import comb
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
UPSTREAM = HERE / "gmc_product_gamma_arbitrary_anchored_schur_thm3110.py"


def load_bank_prefix():
    """Load only the small exact bank constructor, not the upstream audit."""

    tree = ast.parse(UPSTREAM.read_text(encoding="utf-8"))
    prefix = []
    for node in tree.body:
        prefix.append(node)
        if (isinstance(node, ast.Assign)
                and any(isinstance(target, ast.Name) and target.id == "BANKS"
                        for target in node.targets)):
            break
    namespace = {}
    exec(compile(ast.Module(body=prefix, type_ignores=[]),
                 str(UPSTREAM), "exec"), namespace)
    return namespace["BANKS"]


BANKS = load_bank_prefix()
N = int(sys.argv[1]) if len(sys.argv) > 1 else 8
require(N >= 5, "the low-histogram refinement starts at N=5")
THEORETICAL_DEGREE = 4 * N - 9
SHELL_DEGREE = THEORETICAL_DEGREE + 1


@lru_cache(maxsize=None)
def partitions(total, largest=None):
    if total == 0:
        return ((),)
    if largest is None or largest > total:
        largest = total
    return tuple(
        (first,) + tail
        for first in range(largest, 0, -1)
        for tail in partitions(total - first, first)
    )


SHAPES = partitions(N)
ALL_SHAPES = tuple(sorted(
    (shape for degree in range(1, N + 1) for shape in partitions(degree)),
    key=lambda shape: (sum(shape), len(shape), shape),
))
ZERO_TYPE = (1,) * N


@lru_cache(maxsize=None)
def prefix_power(length, exponent):
    return sum(value**exponent for value in range(length))


def residual_power_sums(invariant, row, a_value, b_value):
    middle = max(2 * a_value, b_value)
    return tuple(
        sum(prefix_power(alpha * a_value + beta * b_value, exponent)
            for alpha, beta in row)
        - 3 * prefix_power(a_value, exponent)
        - (invariant + 1) * prefix_power(middle, exponent)
        for exponent in range(1, N + 1)
    )


def all_monomial_values(power_sums):
    """Recover m_lambda from p_1,...,p_N by a triangular recurrence."""

    values = {(): 1}
    for shape in ALL_SHAPES:
        exponent = shape[-1]
        remainder = shape[:-1]
        value = power_sums[exponent - 1] * values[remainder]
        for old_exponent in set(remainder):
            merged = list(remainder)
            merged.remove(old_exponent)
            merged.append(old_exponent + exponent)
            merged = tuple(sorted(merged, reverse=True))
            value -= Counter(merged)[old_exponent + exponent] * values[merged]
        multiplicity = Counter(shape)[exponent]
        require(value % multiplicity == 0,
                "monomial power-sum recurrence lost integrality")
        values[shape] = value // multiplicity
    return values


def chamber_parameters(chamber, first, second):
    if chamber == "tight":
        a_value = first + second + 2
        b_value = first + 2 * second + 3
        return a_value, b_value
    require(chamber == "wide", "unknown chamber")
    a_value = first + 1
    b_value = 2 * a_value + second
    return a_value, b_value


def dominant_row(invariant):
    if invariant == 0:
        return ((0, 3), (2, 0), (2, 0), (1, 0))
    return ((0, 3), (1, 1), (2, 0), (1, 0), (1, 0))


def forced_divisor(invariant, a_value, b_value):
    if invariant == 0:
        return a_value**4 * b_value**2 * (b_value - a_value)**3
    return a_value**3 * b_value**2 * (b_value - a_value)**4


def coefficient_vector(invariant, chamber, first, second):
    a_value, b_value = chamber_parameters(chamber, first, second)
    quotient = all_monomial_values(residual_power_sums(
        invariant, dominant_row(invariant), a_value, b_value))
    quotient_h = sum(quotient[shape] for shape in SHAPES)
    phi_h = 0
    phi_monomial = {shape: 0 for shape in SHAPES}
    for coefficient, row in BANKS[invariant]:
        values = all_monomial_values(residual_power_sums(
            invariant, row, a_value, b_value))
        phi_h += coefficient * sum(values[shape] for shape in SHAPES)
        for shape in SHAPES:
            phi_monomial[shape] += coefficient * values[shape]
    base = forced_divisor(invariant, a_value, b_value)
    answer = {
        shape: Fraction(
            phi_h * quotient[shape]
            - phi_monomial[shape] * quotient_h,
            base,
        )
        for shape in SHAPES
    }
    require(sum(answer.values()) == 0, "coefficient mass did not cancel")
    return answer


def forward_heads(values):
    row = list(values)
    heads = []
    while row:
        heads.append(row[0])
        row = [row[index + 1] - row[index]
               for index in range(len(row) - 1)]
    return heads


def newton_triangle(values, maximum_degree):
    order = maximum_degree + 1
    x_heads = {}
    for second in range(order + 1):
        sequence = [values[first, second]
                    for first in range(order + 1 - second)]
        for first_order, value in enumerate(forward_heads(sequence)):
            x_heads[first_order, second] = value
    coefficients = {}
    for first_order in range(order + 1):
        sequence = [x_heads[first_order, second]
                    for second in range(order + 1 - first_order)]
        for second_order, value in enumerate(forward_heads(sequence)):
            coefficients[first_order, second_order] = value
    return coefficients


def evaluate_newton(coefficients, first, second, degree):
    return sum(
        coefficient * comb(first, i) * comb(second, j)
        for (i, j), coefficient in coefficients.items()
        if i + j <= degree and i <= first and j <= second
    )


@lru_cache(maxsize=None)
def coarsens(fine, coarse):
    """Whether coarse is obtained by merging parts of fine."""

    if sum(fine) != sum(coarse) or len(coarse) > len(fine):
        return False
    if not coarse:
        return not fine
    target = coarse[0]
    remainder = coarse[1:]
    choices = set()

    def choose(start, subtotal, selected):
        if subtotal == target:
            choices.add(tuple(sorted(
                (fine[index] for index in range(len(fine))
                 if index not in selected),
                reverse=True,
            )))
            return
        if subtotal > target:
            return
        previous = None
        for index in range(start, len(fine)):
            if fine[index] == previous:
                continue
            previous = fine[index]
            choose(index + 1, subtotal + fine[index], selected + (index,))

    choose(0, 0, ())
    return any(coarsens(choice, remainder) for choice in choices)


def transport(coefficients, allowed=lambda shape: True):
    """Return an exact max flow and the used coarsening edges."""

    positive = tuple(
        (shape, value) for shape, value in sorted(coefficients.items())
        if value > 0 and shape != ZERO_TYPE and allowed(shape)
    )
    negative = tuple(
        (shape, -value) for shape, value in sorted(coefficients.items())
        if value < 0 and shape != ZERO_TYPE
    )
    demand = sum(value for _, value in negative)
    source = ("source",)
    sink = ("sink",)
    adjacency = defaultdict(list)
    capacity = {}

    def add_edge(left, right, value):
        adjacency[left].append(right)
        adjacency[right].append(left)
        capacity[left, right] = value
        capacity[right, left] = Fraction(0)

    for shape, value in positive:
        add_edge(source, ("positive", shape), value)
    for shape, value in negative:
        add_edge(("negative", shape), sink, value)
    for high, _ in positive:
        for low, _ in negative:
            if coarsens(low, high):
                add_edge(("positive", high), ("negative", low), demand)

    original = dict(capacity)
    value = Fraction(0)
    while True:
        parent = {source: None}
        queue = [source]
        for vertex in queue:
            for neighbour in adjacency[vertex]:
                if neighbour not in parent and capacity[vertex, neighbour] > 0:
                    parent[neighbour] = vertex
                    queue.append(neighbour)
            if sink in parent:
                break
        if sink not in parent:
            break
        increment = demand
        vertex = sink
        while parent[vertex] is not None:
            increment = min(increment, capacity[parent[vertex], vertex])
            vertex = parent[vertex]
        vertex = sink
        while parent[vertex] is not None:
            previous = parent[vertex]
            capacity[previous, vertex] -= increment
            capacity[vertex, previous] += increment
            vertex = previous
        value += increment

    used = []
    for high, _ in positive:
        for low, _ in negative:
            edge = (("positive", high), ("negative", low))
            if edge in original:
                amount = original[edge] - capacity[edge]
                if amount:
                    used.append((high, low, amount))
    return value == demand, tuple(used), len(positive), len(negative)


def is_hook(shape):
    return all(part == 1 for part in shape[1:])


def run_case(invariant, chamber):
    grid = {
        (first, second): coefficient_vector(
            invariant, chamber, first, second)
        for second in range(SHELL_DEGREE + 1)
        for first in range(SHELL_DEGREE + 1 - second)
    }
    polynomials = {
        shape: newton_triangle(
            {(first, second): vector[shape]
             for (first, second), vector in grid.items()},
            THEORETICAL_DEGREE,
        )
        for shape in SHAPES
    }
    require(all(
        polynomials[shape][first, second] == 0
        for shape in SHAPES
        for first in range(SHELL_DEGREE + 1)
        for second in range(SHELL_DEGREE + 1 - first)
        if first + second == SHELL_DEGREE
    ), "degree-bound shell failed")

    off_grid = (
        (THEORETICAL_DEGREE + 2, 2),
        (2, THEORETICAL_DEGREE + 2),
        (THEORETICAL_DEGREE + 3, THEORETICAL_DEGREE + 1),
    )
    for first, second in off_grid:
        exact = coefficient_vector(invariant, chamber, first, second)
        for shape in SHAPES:
            require(evaluate_newton(
                polynomials[shape], first, second, THEORETICAL_DEGREE,
            ) == exact[shape], "off-grid reconstruction failed")

    census = Counter()
    edge_union = set()
    sign_patterns = set()
    singleton = Counter()
    certificate = []
    families = {
        "row": lambda shape: len(shape) == 1,
        "hooks": is_hook,
        "two_row": lambda shape: len(shape) <= 2,
        "hook_or_two_row": lambda shape: is_hook(shape) or len(shape) <= 2,
    }
    for first in range(THEORETICAL_DEGREE + 1):
        for second in range(THEORETICAL_DEGREE + 1 - first):
            vector = {
                shape: Fraction(polynomials[shape][first, second])
                for shape in SHAPES
            }
            require(sum(vector.values()) == 0,
                    "Newton coefficient mass did not cancel")
            singleton["positive" if vector[ZERO_TYPE] > 0 else
                      "negative" if vector[ZERO_TYPE] < 0 else "zero"] += 1
            if not any(value for shape, value in vector.items()
                       if shape != ZERO_TYPE):
                census["zero"] += 1
                continue
            passed, used, positive_count, negative_count = transport(vector)
            require(passed, f"full N={N} refinement flow failed")
            census["flow"] += 1
            census["used_edges"] += len(used)
            census["max_edges"] = max(census["max_edges"], len(used))
            census["max_positive"] = max(census["max_positive"], positive_count)
            census["max_negative"] = max(census["max_negative"], negative_count)
            sign_patterns.add(tuple(
                (shape, 1 if value > 0 else -1)
                for shape, value in sorted(vector.items())
                if value and shape != ZERO_TYPE
            ))
            for high, low, amount in used:
                edge_union.add((high, low))
                certificate.append(
                    f"{invariant}|{chamber}|{first},{second}|"
                    f"{high}>{low}|{amount}"
                )
            for name, predicate in families.items():
                if transport(vector, predicate)[0]:
                    census[name] += 1

    digest = hashlib.sha256("\n".join(certificate).encode("ascii")).hexdigest()
    return (
        f"N{N}:I{invariant + 1}:{chamber}:degree={THEORETICAL_DEGREE}:"
        f"flow={census['flow']}:zero={census['zero']}:"
        f"sign_patterns={len(sign_patterns)}:used_edges={census['used_edges']}:"
        f"max_edges={census['max_edges']}:edge_union={len(edge_union)}:"
        f"max_posneg={census['max_positive']}/{census['max_negative']}:"
        f"restricted=row:{census['row']},hooks:{census['hooks']},"
        f"two_row:{census['two_row']},union:{census['hook_or_two_row']}:"
        f"singleton={dict(singleton)}:shell=PASS:offgrid=PASS:"
        f"digest={digest}"
    )


def main():
    print("degree_reason=deg(m_mu)<=N+length(mu);base_degree=9")
    print("relevant_degree=4N-10 because Lbar_(1^N)=0")
    for invariant in (0, 1):
        for chamber in ("tight", "wide"):
            print(run_case(invariant, chamber), flush=True)
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
