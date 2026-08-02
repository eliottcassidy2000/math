#!/usr/bin/env python3
"""Exact chamber-Newton refinement certificates for THM-3115.

For N=5,6,7 this script expands the row-normalized THM-3110 operator in
the monomial-fibre Young-subgroup basis from THM-3112.  Every chamber-Newton
coefficient vector is certified positive by an exact rational transport from
negative fibre types to positive coarsenings.
"""

import hashlib
import importlib.util
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
SPEC = importlib.util.spec_from_file_location("thm3110_schur", UPSTREAM)
THM3110 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(THM3110)

DEGREES = (5, 6, 7)
MAXIMUM_DEGREE = max(DEGREES)
INVARIANTS = (0, 1)
CHAMBERS = ("tight", "wide")


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


ALL_PARTITIONS = tuple(sorted(
    (
        partition
        for total in range(1, MAXIMUM_DEGREE + 1)
        for partition in partitions(total)
    ),
    key=lambda partition: (sum(partition), len(partition), partition),
))


def all_monomial_values(power_sums):
    """Recover every m_mu from p_1,...,p_MAXIMUM_DEGREE."""

    values = {(): 1}
    for partition in ALL_PARTITIONS:
        exponent = partition[-1]
        remainder = partition[:-1]
        value = power_sums[exponent - 1] * values[remainder]
        for old_exponent in set(remainder):
            merged = list(remainder)
            merged.remove(old_exponent)
            merged.append(old_exponent + exponent)
            merged = tuple(sorted(merged, reverse=True))
            value -= (Counter(merged)[old_exponent + exponent]
                      * values[merged])
        multiplicity = Counter(partition)[exponent]
        require(value % multiplicity == 0,
                "monomial power-sum recurrence lost integrality")
        values[partition] = value // multiplicity
    return values


@lru_cache(maxsize=None)
def prefix_power(length, exponent):
    return sum(value**exponent for value in range(length))


def residual_power_sums(invariant, row, a, b):
    """Power sums of E(row) minus the exact THM-3110 common roots."""

    middle = max(2 * a, b)
    return tuple(
        sum(prefix_power(alpha * a + beta * b, exponent)
            for alpha, beta in row)
        - 3 * prefix_power(a, exponent)
        - (invariant + 1) * prefix_power(middle, exponent)
        for exponent in range(1, MAXIMUM_DEGREE + 1)
    )


@lru_cache(maxsize=None)
def coarsens(fine, coarse):
    """Return whether the block type coarse is a coarsening of fine."""

    if sum(fine) != sum(coarse) or len(coarse) > len(fine):
        return False
    fine = tuple(sorted(fine, reverse=True))
    coarse = tuple(sorted(coarse, reverse=True))

    @lru_cache(maxsize=None)
    def recurse(targets, remaining):
        if not targets:
            return not remaining
        target = targets[0]
        count = len(remaining)
        remainders = set()

        def choose(start, subtotal, selected):
            if subtotal == target:
                remainders.add(tuple(sorted(
                    (remaining[index] for index in range(count)
                     if index not in selected),
                    reverse=True,
                )))
                return
            if subtotal > target:
                return
            previous = None
            for index in range(start, count):
                if remaining[index] == previous:
                    continue
                previous = remaining[index]
                choose(
                    index + 1,
                    subtotal + remaining[index],
                    selected + (index,),
                )

        choose(0, 0, ())
        return any(recurse(targets[1:], rest) for rest in remainders)

    return recurse(coarse, fine)


def rational_transport(coefficients):
    """Realize the coefficient vector as a positive refinement boundary."""

    total = sum(next(iter(coefficients)))
    zero_type = (1,) * total
    positive = tuple(sorted(
        (partition, value)
        for partition, value in coefficients.items()
        if value > 0
    ))
    negative = tuple(sorted(
        (partition, -value)
        for partition, value in coefficients.items()
        if value < 0 and partition != zero_type
    ))
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

    for partition, value in positive:
        add_edge(source, ("positive", partition), value)
    for partition, value in negative:
        add_edge(("negative", partition), sink, value)
    for high, _ in positive:
        for low, _ in negative:
            if coarsens(low, high):
                add_edge(("positive", high), ("negative", low), demand)

    original = dict(capacity)
    flow = Fraction(0)
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
            increment = min(
                increment, capacity[parent[vertex], vertex])
            vertex = parent[vertex]
        vertex = sink
        while parent[vertex] is not None:
            previous = parent[vertex]
            capacity[previous, vertex] -= increment
            capacity[vertex, previous] += increment
            vertex = previous
        flow += increment

    require(flow == demand, "refinement transport did not cover all debt")
    require(coefficients[zero_type] <= 0,
            "the trivial-gap root carried positive unmatched mass")
    used = []
    for high, _ in positive:
        for low, _ in negative:
            edge = (("positive", high), ("negative", low))
            if edge not in original:
                continue
            amount = original[edge] - capacity[edge]
            if amount:
                require(coarsens(low, high), "transport used a non-coarsening")
                used.append((high, low, amount))

    # The total coefficient mass is zero.  Once every nontrivial negative
    # type is paid, the unused positive supply is exactly the remaining debt
    # at the singleton (trivial-subgroup) root.  Route it there to obtain an
    # exact positive boundary, not a decomposition with a positive leftover.
    singleton_debt = -coefficients[zero_type]
    singleton_flow = Fraction(0)
    for high, _ in positive:
        residual = capacity[source, ("positive", high)]
        if residual:
            require(high != zero_type,
                    "positive singleton cannot refine a nontrivial gap")
            require(coarsens(zero_type, high),
                    "positive residual did not coarsen the singleton root")
            used.append((high, zero_type, residual))
            singleton_flow += residual
    require(singleton_flow == singleton_debt,
            "positive residual did not exactly pay singleton debt")
    return tuple(used)


def dominant_row(invariant):
    if invariant == 0:
        return ((0, 3), (2, 0), (2, 0), (1, 0))
    return ((0, 3), (1, 1), (2, 0), (1, 0), (1, 0))


def evaluation_vectors(invariant, chamber, first, second):
    """Return G_mu for all three degrees at one chamber-grid point."""

    a, b = THM3110.chamber_parameters(chamber, first, second)
    bank = THM3110.BANKS[invariant]
    quotient = all_monomial_values(residual_power_sums(
        invariant, dominant_row(invariant), a, b))
    atoms = tuple(
        (coefficient, all_monomial_values(residual_power_sums(
            invariant, row, a, b)))
        for coefficient, row in bank
    )
    base = THM3110.forced_divisor(invariant, a, b)
    answer = {}
    for degree in DEGREES:
        shapes = partitions(degree)
        phi_h = sum(
            sum(coefficient * values[shape]
                for coefficient, values in atoms)
            for shape in shapes
        )
        quotient_h = sum(quotient[shape] for shape in shapes)
        vector = {}
        for shape in shapes:
            phi_monomial = sum(
                coefficient * values[shape]
                for coefficient, values in atoms
            )
            vector[shape] = Fraction(
                phi_h * quotient[shape] - phi_monomial * quotient_h,
                base,
            )
        require(sum(vector.values()) == 0,
                "monomial coefficient mass did not cancel")
        answer[degree] = vector
    return answer


def evaluate_newton(polynomial, first, second, degree):
    return sum(
        coefficient * comb(first, i) * comb(second, j)
        for (i, j), coefficient in polynomial.items()
        if i + j <= degree and i <= first and j <= second
    )


def fraction_text(value):
    return (str(value.numerator) if value.denominator == 1
            else f"{value.numerator}/{value.denominator}")


records = []
all_certificate_lines = []
for invariant in INVARIANTS:
    for chamber in CHAMBERS:
        grid_points = set()
        off_grid = {}
        for degree in DEGREES:
            # A type mu has endpoint degree at most N+ell(mu).  After the
            # degree-nine collision divisor, every nontrivial Young gap has
            # degree at most 4N-10 because ell(mu)<=N-1.  The sole length-N
            # type is the trivial-subgroup zero gap; coefficient-mass
            # cancellation removes its possible leading term as well.
            polynomial_degree = 4 * degree - 10
            grid_points.update(
                (first, second)
                for second in range(polynomial_degree + 2)
                for first in range(polynomial_degree + 2 - second)
            )
            off_grid[degree] = (
                (polynomial_degree + 2, 2),
                (2, polynomial_degree + 2),
                (polynomial_degree + 3, polynomial_degree + 1),
            )
            grid_points.update(off_grid[degree])

        evaluations = {
            point: evaluation_vectors(invariant, chamber, *point)
            for point in sorted(grid_points)
        }

        for degree in DEGREES:
            polynomial_degree = 4 * degree - 10
            shapes = partitions(degree)
            polynomials = {}
            for shape in shapes:
                values = {
                    (first, second):
                        evaluations[first, second][degree][shape]
                    for second in range(polynomial_degree + 2)
                    for first in range(polynomial_degree + 2 - second)
                }
                polynomials[shape] = THM3110.newton_triangle(
                    values, polynomial_degree)

            require(all(
                polynomials[shape][first, second] == 0
                for shape in shapes
                for first in range(polynomial_degree + 2)
                for second in range(polynomial_degree + 2 - first)
                if first + second == polynomial_degree + 1
            ), "excess-degree Newton shell did not vanish")

            for first, second in off_grid[degree]:
                exact = evaluations[first, second][degree]
                for shape in shapes:
                    require(
                        evaluate_newton(
                            polynomials[shape], first, second,
                            polynomial_degree,
                        ) == exact[shape],
                        "off-grid Newton reconstruction failed",
                    )

            slots = (polynomial_degree + 1) * (polynomial_degree + 2) // 2
            nonzero = zero = transport_slots = direct_slots = 0
            transport_edges = maximum_edges = 0
            edge_union = set()
            for first in range(polynomial_degree + 1):
                for second in range(polynomial_degree + 1 - first):
                    coefficient_vector = {
                        shape: Fraction(polynomials[shape][first, second])
                        for shape in shapes
                    }
                    require(sum(coefficient_vector.values()) == 0,
                            "Newton-slot coefficient mass did not cancel")
                    for shape in shapes:
                        all_certificate_lines.append(
                            f"C|{invariant + 1}|{chamber}|{degree}|"
                            f"{first},{second}|{shape}|"
                            f"{fraction_text(coefficient_vector[shape])}"
                        )
                    if not any(coefficient_vector.values()):
                        zero += 1
                        continue
                    flow = rational_transport(coefficient_vector)
                    nonzero += 1
                    if flow:
                        transport_slots += 1
                    else:
                        direct_slots += 1
                    transport_edges += len(flow)
                    maximum_edges = max(maximum_edges, len(flow))
                    for high, low, amount in flow:
                        edge_union.add((high, low))
                        all_certificate_lines.append(
                            f"{invariant + 1}|{chamber}|{degree}|"
                            f"{first},{second}|{high}>{low}|"
                            f"{fraction_text(amount)}"
                        )

            require(nonzero + zero == slots, "Newton slot census drift")
            require(zero == 0,
                    "unexpected zero inside the sharp Newton triangle")
            records.append(
                f"N{degree}:I{invariant + 1}:{chamber}:"
                f"degree={polynomial_degree}:slots={slots}:"
                f"vectors={nonzero}:transport={transport_slots}:"
                f"direct={direct_slots}:zero={zero}:"
                f"used_edges={transport_edges}:max_edges={maximum_edges}:"
                f"edge_union={len(edge_union)}:shell=PASS:offgrid=PASS"
            )

certificate_digest = hashlib.sha256(
    "\n".join(all_certificate_lines).encode("ascii")
).hexdigest()

print("fibre_identity=beta_N(S)=sum_mu m_mu(S)*Lbar_mu")
print("refinement_order=nu_coarsens_mu_implies_Lbar_mu<=Lbar_nu")
print("coefficient=G_mu=(Phi(hN)*m_mu(Q)-Phi(m_mu)*hN(Q))/base")
print("degree_bound=4N-10;excess_shell_zero")
for record in records:
    print(record)
print(f"certificate_sha256={certificate_digest}")
print("global_scope=N=5,6,7;both_banks;all_integer_0<a<b")
print("all_exact_checks=PASS")
