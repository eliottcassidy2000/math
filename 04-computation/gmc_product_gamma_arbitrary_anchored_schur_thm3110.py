#!/usr/bin/env python3
"""Exact Schur/Newton companion for THM-3110.

For anchored support {0,a,b}, 0<a<b, this script derives the two cleared
24/25-atom product-Gamma response banks.  After their chamber-dependent
common rising-factor deletion, it evaluates the signed residual-alphabet
functional

    Phi_i(f) = sum_R c_R f(S_R).

For every partition lambda of size 5 through 12 it proves, by exact
Gregory--Newton interpolation, that Phi_i(s_lambda) divided by its forced
degree-nine divisor has nonnegative coefficients in both integer chambers.
It also runs a separate six-support hostile scan through size 14.

Only Python integers and fractions are used.
"""

from collections import Counter, defaultdict, deque
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import comb, gcd, prod


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ZERO_INDEX = (0, 0)
A_INDEX = (1, 0)
B_INDEX = (0, 1)


def clean(poly):
    return {exponent: coefficient for exponent, coefficient in poly.items()
            if coefficient}


def add(*polys):
    out = defaultdict(int)
    for poly in polys:
        for exponent, coefficient in poly.items():
            out[exponent] += coefficient
    return clean(out)


def scale(poly, scalar):
    return clean({exponent: scalar * coefficient
                  for exponent, coefficient in poly.items()})


def multiply(left, right):
    out = defaultdict(int)
    for exponent, coefficient in left.items():
        for other, scalar in right.items():
            merged = defaultdict(int, dict(exponent))
            for index, power in other:
                merged[index] += power
            key = tuple(sorted((index, power) for index, power in merged.items()
                               if power))
            out[key] += coefficient * scalar
    return clean(out)


def response_moment(indices):
    exponent = defaultdict(int)
    exponent[tuple(map(sum, zip(*indices)))] += 1
    for index in indices:
        exponent[index] -= 1
    # w_0=1, so its Laurent exponent is immaterial.
    exponent[ZERO_INDEX] = 0
    key = tuple(sorted((index, power) for index, power in exponent.items()
                       if power))
    return {key: 1}


def contract(vectors):
    out = {}
    for choices in product(*vectors):
        out = add(
            out,
            scale(
                response_moment([item[0] for item in choices]),
                prod(item[1] for item in choices),
            ),
        )
    return out


def cleared_banks():
    """Return the exact oriented 24/25 moment-row banks."""

    u = ((A_INDEX, 1), (ZERO_INDEX, -1))
    v = ((B_INDEX, 1), (A_INDEX, -1))
    g11 = contract((u, u))
    g12 = contract((u, v))
    g22 = contract((v, v))
    t111 = contract((u, u, u))
    t112 = contract((u, u, v))
    t122 = contract((u, v, v))
    t222 = contract((v, v, v))

    i1 = add(
        scale(multiply(multiply(t112, g11), g22), 3),
        scale(multiply(t222, multiply(g11, g11)), -1),
        scale(multiply(multiply(t111, g12), g22), -2),
    )
    i2 = add(
        scale(multiply(multiply(t122, g11), g22), 3),
        scale(multiply(multiply(t222, g12), g11), -2),
        scale(multiply(t111, multiply(g22, g22)), -1),
    )

    banks = []
    for invariant, denominator in (
        (i1, {A_INDEX: 5, B_INDEX: 3}),
        (i2, {A_INDEX: 5, B_INDEX: 4}),
    ):
        collected = defaultdict(int)
        for exponent, coefficient in invariant.items():
            cleared = defaultdict(int, dict(exponent))
            for index, power in denominator.items():
                cleared[index] += power
            require(all(power >= 0 for power in cleared.values()),
                    "moment denominator did not clear")
            rows = []
            for index, power in cleared.items():
                rows.extend([index] * power)
            # The orientation is F_i=-w_a^5 w_b^(i+2) I_i.
            collected[tuple(sorted(rows, reverse=True))] -= coefficient
        bank = tuple(sorted(
            ((coefficient, rows) for rows, coefficient in collected.items()
             if coefficient),
            key=lambda item: (item[1], item[0]),
        ))
        banks.append(bank)

    require(tuple(map(len, banks)) == (24, 25), "generic atom count drift")
    require(tuple(sum(c for c, _ in bank) for bank in banks) == (0, 0),
            "signed atom mass drift")
    require(tuple(sum(c for c, _ in bank if c > 0) for bank in banks)
            == (37, 39), "positive atom mass drift")
    return tuple(banks)


BANKS = cleared_banks()


def chamber_parameters(chamber, first, second):
    """Parametrize every integer point in one of the two b/2a chambers."""

    if chamber == "tight":
        # d=b-a=second+1 and a-d=first+1.
        a = first + second + 2
        b = first + 2 * second + 3
        require(a < b < 2 * a, "tight chamber parametrization drift")
        return a, b
    require(chamber == "wide", "unknown chamber")
    a = first + 1
    b = 2 * a + second
    require(b >= 2 * a, "wide chamber parametrization drift")
    return a, b


def forced_divisor(invariant, a, b):
    if invariant == 0:
        return a ** 4 * b ** 2 * (b - a) ** 3
    return a ** 3 * b ** 2 * (b - a) ** 4


def expected_residual_degree(invariant, a, b):
    middle = max(2 * a, b)
    if invariant == 0:
        return 2 * a + 3 * b - middle
    return 2 * a + 4 * b - 2 * middle


def residual_roots(invariant, rows, a, b):
    """Delete the exact common rising-factor root multiset."""

    middle = max(2 * a, b)
    common = 3 * list(range(a)) + (invariant + 1) * list(range(middle))
    counts = Counter()
    for alpha, beta in rows:
        counts.update(range(alpha * a + beta * b))
    counts.subtract(common)
    require(all(multiplicity >= 0 for multiplicity in counts.values()),
            "claimed common root divisor is not present")
    roots = sorted(counts.elements())
    require(len(roots) == expected_residual_degree(invariant, a, b),
            "residual degree drift")
    # Some non-dominant atoms retain additional zero roots.  They count in the
    # monic theta-degree but disappear from the reciprocal response degree.
    return roots


def complete_and_elementary(roots, degree):
    complete = [0] * (degree + 1)
    elementary = [0] * (degree + 1)
    complete[0] = elementary[0] = 1
    for root in roots:
        # Multiplication by (1-root*t)^(-1), using updated lower entries.
        for index in range(1, degree + 1):
            complete[index] += root * complete[index - 1]
        # Multiplication by (1+root*t), using old lower entries.
        for index in range(degree, 0, -1):
            elementary[index] += root * elementary[index - 1]
    return tuple(complete), tuple(elementary)


def bareiss_determinant(matrix):
    size = len(matrix)
    if size == 0:
        return 1
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for pivot_index in range(size - 1):
        if work[pivot_index][pivot_index] == 0:
            swap = next((row for row in range(pivot_index + 1, size)
                         if work[row][pivot_index]), None)
            if swap is None:
                return 0
            work[pivot_index], work[swap] = work[swap], work[pivot_index]
            sign = -sign
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = (
                    work[row][column] * pivot
                    - work[row][pivot_index] * work[pivot_index][column]
                )
                require(numerator % previous == 0,
                        "Bareiss exact division failed")
                work[row][column] = numerator // previous
        previous = pivot
    return sign * work[-1][-1]


def conjugate_partition(partition):
    return tuple(sum(part >= column for part in partition)
                 for column in range(1, partition[0] + 1))


def schur_value(symmetric_data, partition):
    complete, elementary = symmetric_data
    # Use the shorter of the two Jacobi--Trudi determinants.
    if len(partition) <= partition[0]:
        shape = partition
        values = complete
    else:
        shape = conjugate_partition(partition)
        values = elementary
    size = len(shape)
    matrix = []
    for row in range(size):
        entries = []
        for column in range(size):
            index = shape[row] - row + column
            entries.append(values[index] if index >= 0 else 0)
        matrix.append(entries)
    return bareiss_determinant(matrix)


@lru_cache(maxsize=None)
def partitions(total, largest=None):
    if total == 0:
        return ((),)
    if largest is None or largest > total:
        largest = total
    out = []
    for first in range(largest, 0, -1):
        for rest in partitions(total - first, first):
            out.append((first,) + rest)
    return tuple(out)


def atom_symmetric_data(invariant, a, b, degree):
    return tuple(
        (coefficient,
         complete_and_elementary(
             residual_roots(invariant, rows, a, b), degree))
        for coefficient, rows in BANKS[invariant]
    )


def phi_from_data(data, partition):
    return sum(coefficient * schur_value(symmetric_data, partition)
               for coefficient, symmetric_data in data)


def forward_heads(values):
    """Return Delta^j values[0] for every possible order j."""

    row = list(values)
    heads = []
    while row:
        heads.append(row[0])
        row = [row[index + 1] - row[index]
               for index in range(len(row) - 1)]
    return heads


def newton_triangle(values, maximum_degree):
    """Recover c_ij in sum c_ij binom(A,i)binom(D,j)."""

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


GENERIC_EXPECTED = {
    5: (7, 1, 84, 84, 0, Fraction(2)),
    6: (11, 3, 440, 440, 0, Fraction(5)),
    7: (15, 5, 1260, 1258, 2, Fraction(156)),
    8: (22, 7, 3168, 3164, 4, Fraction(628, 3)),
    9: (30, 9, 6600, 6589, 11, Fraction(420)),
    10: (42, 11, 13104, 13083, 21, Fraction(9984)),
    11: (56, 13, 23520, 23479, 41, Fraction(8400)),
    12: (77, 15, 41888, 41819, 69, Fraction(2672)),
}


def generic_schur_audit():
    maximum_partition_degree = 12
    maximum_newton_degree = 2 * maximum_partition_degree - 9
    audit_points = {
        (first, second)
        for first in range(maximum_newton_degree + 2)
        for second in range(maximum_newton_degree + 2 - first)
    }
    for degree in range(5, maximum_partition_degree + 1):
        newton_degree = 2 * degree - 9
        audit_points.update({
            (newton_degree + 2, newton_degree + 1),
            (newton_degree + 1, newton_degree + 2),
            (newton_degree + 2, 2),
        })

    by_degree = {
        degree: {"slots": 0, "positive": 0, "zero": 0,
                 "minimum": None}
        for degree in range(5, maximum_partition_degree + 1)
    }
    degree_five_newton = {}

    for invariant in (0, 1):
        for chamber in ("tight", "wide"):
            cached = {}
            for first, second in sorted(audit_points):
                a, b = chamber_parameters(chamber, first, second)
                cached[first, second] = (
                    forced_divisor(invariant, a, b),
                    atom_symmetric_data(
                        invariant, a, b, maximum_partition_degree),
                )

            for degree in range(5, maximum_partition_degree + 1):
                newton_degree = 2 * degree - 9
                for partition in partitions(degree):
                    values = {}
                    for first in range(newton_degree + 2):
                        for second in range(newton_degree + 2 - first):
                            divisor, data = cached[first, second]
                            values[first, second] = Fraction(
                                phi_from_data(data, partition), divisor)

                    coefficients = newton_triangle(values, newton_degree)
                    for (first_order, second_order), coefficient in coefficients.items():
                        if first_order + second_order <= newton_degree:
                            by_degree[degree]["slots"] += 1
                            require(coefficient >= 0,
                                    "negative Schur/Newton coefficient")
                            if coefficient == 0:
                                by_degree[degree]["zero"] += 1
                            else:
                                by_degree[degree]["positive"] += 1
                                current = by_degree[degree]["minimum"]
                                if current is None or coefficient < current:
                                    by_degree[degree]["minimum"] = coefficient
                        elif first_order + second_order == newton_degree + 1:
                            require(coefficient == 0,
                                    "chamber degree guard failed")

                    for first, second in (
                        (newton_degree + 2, newton_degree + 1),
                        (newton_degree + 1, newton_degree + 2),
                        (newton_degree + 2, 2),
                    ):
                        divisor, data = cached[first, second]
                        exact = Fraction(phi_from_data(data, partition), divisor)
                        reconstructed = evaluate_newton(
                            coefficients, first, second, newton_degree)
                        require(exact == reconstructed,
                                "off-grid Newton reconstruction failed")

                    if degree == 5 and partition == (1, 1, 1, 1, 1):
                        key = (invariant, chamber)
                        degree_five_newton[key] = tuple(
                            coefficients[index]
                            for index in ((0, 0), (1, 0), (0, 1)))

    expected_degree_five = {
        (0, "tight"): (Fraction(8), Fraction(4), Fraction(13, 2)),
        (0, "wide"): (Fraction(4), Fraction(13, 2), Fraction(5, 2)),
        (1, "tight"): (Fraction(13, 2), Fraction(7, 2), Fraction(11, 2)),
        (1, "wide"): (Fraction(3), Fraction(11, 2), Fraction(2)),
    }
    require(degree_five_newton == expected_degree_five,
            "degree-five kappa table drift")

    total_slots = total_positive = total_zero = 0
    lines = []
    for degree, expected in GENERIC_EXPECTED.items():
        partition_count, newton_degree, slots, positive, zero, minimum = expected
        result = by_degree[degree]
        require(len(partitions(degree)) == partition_count,
                "partition census drift")
        require(newton_degree == 2 * degree - 9, "degree formula drift")
        require((result["slots"], result["positive"], result["zero"],
                 result["minimum"]) == (slots, positive, zero, minimum),
                f"generic degree-{degree} census drift")
        total_slots += slots
        total_positive += positive
        total_zero += zero
        lines.append(
            f"schur_n={degree}:partitions={partition_count}:degree={newton_degree}:"
            f"slots={slots}:positive={positive}:zero={zero}:min={minimum}"
        )

    require((total_slots, total_positive, total_zero)
            == (90064, 89916, 148), "generic total census drift")
    return lines, (total_slots, total_positive, total_zero)


FINITE_SUPPORTS = ((1, 2), (1, 3), (2, 3), (1, 4), (2, 5), (4, 5))


def finite_hostile_scan():
    checked = positive = zero = negative = 0
    for invariant in (0, 1):
        for a, b in FINITE_SUPPORTS:
            data = atom_symmetric_data(invariant, a, b, 14)
            alphabet_degree = expected_residual_degree(invariant, a, b)
            for degree in range(5, 15):
                for partition in partitions(degree):
                    value = phi_from_data(data, partition)
                    checked += 1
                    if value < 0:
                        negative += 1
                    elif value == 0:
                        zero += 1
                    else:
                        positive += 1
                    require((value == 0) == (len(partition) > alphabet_degree),
                            "finite support equality boundary drift")
    require((checked, positive, zero, negative) == (5952, 5520, 432, 0),
            "finite hostile census drift")
    return checked, positive, zero, negative


def row_value(row, ratio):
    return row[0] + row[1] * ratio


def row_walls(positive_rows, negative_rows):
    forms = set(positive_rows) | set(negative_rows) | {ZERO_INDEX}
    walls = set()
    for left in forms:
        for right in forms:
            delta_a = left[0] - right[0]
            delta_b = left[1] - right[1]
            if delta_b:
                wall = Fraction(-delta_a, delta_b)
                if wall > 0:
                    walls.add(wall)
    return walls


def rows_majorize_in_chamber(positive_rows, negative_rows, chamber):
    """Exact row-majorization test on the whole ratio chamber.

    Row majorization implies coordinatewise dominance of the corresponding
    interval-root multisets: every tail count is a sum of the convex functions
    (row_length-t)_+.  The common deleted roots cancel from both sides.
    """

    length = max(len(positive_rows), len(negative_rows))
    positive_rows = tuple(positive_rows) + (ZERO_INDEX,) * (
        length - len(positive_rows))
    negative_rows = tuple(negative_rows) + (ZERO_INDEX,) * (
        length - len(negative_rows))
    if tuple(map(sum, zip(*positive_rows))) != tuple(map(sum, zip(*negative_rows))):
        return False

    walls = row_walls(positive_rows, negative_rows) | {Fraction(1), Fraction(2)}
    if chamber == "tight":
        cuts = sorted(wall for wall in walls if 1 <= wall <= 2)
        cells = [("point", wall, wall) for wall in cuts]
        cells += [("interval", left, right)
                  for left, right in zip(cuts, cuts[1:]) if left < right]
    else:
        require(chamber == "wide", "unknown transport chamber")
        cuts = sorted(wall for wall in walls if wall >= 2)
        if not cuts or cuts[0] != 2:
            cuts.insert(0, Fraction(2))
        cells = [("point", wall, wall) for wall in cuts]
        cells += [("interval", left, right)
                  for left, right in zip(cuts, cuts[1:]) if left < right]
        cells.append(("ray", cuts[-1], None))

    for kind, left, right in cells:
        if kind == "point":
            sorting_ratio = left
            endpoints = (left,)
        elif kind == "interval":
            sorting_ratio = (left + right) / 2
            endpoints = (left, right)
        else:
            sorting_ratio = left + 1
            endpoints = ()

        positive = sorted(
            positive_rows, key=lambda row: row_value(row, sorting_ratio),
            reverse=True)
        negative = sorted(
            negative_rows, key=lambda row: row_value(row, sorting_ratio),
            reverse=True)
        delta_a = delta_b = 0
        for positive_row, negative_row in zip(positive, negative):
            delta_a += positive_row[0] - negative_row[0]
            delta_b += positive_row[1] - negative_row[1]
            if kind == "ray":
                if delta_b < 0 or delta_a + delta_b * left < 0:
                    return False
            elif any(delta_a + delta_b * endpoint < 0
                     for endpoint in endpoints):
                return False
    return True


class Dinic:
    def __init__(self, size):
        self.graph = [[] for _ in range(size)]

    def add_edge(self, source, target, capacity):
        forward_index = len(self.graph[source])
        reverse_index = len(self.graph[target])
        self.graph[source].append([target, capacity, reverse_index])
        self.graph[target].append([source, 0, forward_index])

    def maximum_flow(self, source, target):
        total = 0
        while True:
            level = [-1] * len(self.graph)
            level[source] = 0
            queue = deque([source])
            while queue:
                vertex = queue.popleft()
                for other, capacity, _ in self.graph[vertex]:
                    if capacity and level[other] < 0:
                        level[other] = level[vertex] + 1
                        queue.append(other)
            if level[target] < 0:
                return total

            cursor = [0] * len(self.graph)

            def augment(vertex, amount):
                if vertex == target:
                    return amount
                while cursor[vertex] < len(self.graph[vertex]):
                    edge = self.graph[vertex][cursor[vertex]]
                    other, capacity, reverse = edge
                    if capacity and level[other] == level[vertex] + 1:
                        sent = augment(other, min(amount, capacity))
                        if sent:
                            edge[1] -= sent
                            self.graph[other][reverse][1] += sent
                            return sent
                    cursor[vertex] += 1
                return 0

            while True:
                sent = augment(source, 10 ** 9)
                if not sent:
                    break
                total += sent


def first_order_transport_audit():
    results = []
    expected = {
        (0, "tight"): (58, 28),
        (0, "wide"): (61, 28),
        (1, "tight"): (66, 28),
        (1, "wide"): (66, 28),
    }
    for invariant, bank in enumerate(BANKS):
        positive = tuple((coefficient, rows) for coefficient, rows in bank
                         if coefficient > 0)
        negative = tuple((-coefficient, rows) for coefficient, rows in bank
                         if coefficient < 0)
        positive_mass = sum(coefficient for coefficient, _ in positive)
        negative_mass = sum(coefficient for coefficient, _ in negative)
        require(positive_mass == negative_mass == (37, 39)[invariant],
                "transport mass balance drift")

        for chamber in ("tight", "wide"):
            edges = tuple(
                (positive_index, negative_index)
                for positive_index, (_, positive_rows) in enumerate(positive)
                for negative_index, (_, negative_rows) in enumerate(negative)
                if rows_majorize_in_chamber(
                    positive_rows, negative_rows, chamber)
            )
            size = 2 + len(positive) + len(negative)
            source = size - 2
            target = size - 1
            flow = Dinic(size)
            for index, (capacity, _) in enumerate(positive):
                flow.add_edge(source, index, capacity)
            for index, (capacity, _) in enumerate(negative):
                flow.add_edge(len(positive) + index, target, capacity)
            for positive_index, negative_index in edges:
                flow.add_edge(positive_index, len(positive) + negative_index,
                              positive_mass)
            capacity = flow.maximum_flow(source, target)
            require((len(edges), capacity) == expected[invariant, chamber],
                    "first-order transport capacity drift")
            results.append((invariant, chamber, len(edges), capacity,
                            positive_mass - capacity))
    return tuple(results)


def gcd_of_integers(values):
    return prod(()) if not values else __import__("functools").reduce(gcd, values)


def main():
    # A small exact primitive-content guard on both signed coefficient banks.
    contents = tuple(gcd_of_integers([abs(c) for c, _ in bank]) for bank in BANKS)
    require(contents == (1, 1), "atom bank content drift")

    schur_lines, totals = generic_schur_audit()
    finite = finite_hostile_scan()
    transport = first_order_transport_audit()

    print("atom_banks=I1:24:mass37;I2:25:mass39;primitive=1,1")
    print("common_roots=I1:3*[0,a)+1*[0,max(2a,b));"
          "I2:3*[0,a)+2*[0,max(2a,b))")
    print("divisors=I1:a^4*b^2*(b-a)^3;I2:a^3*b^2*(b-a)^4")
    print("degree5_kappa=I1:(3a+5b-5)/2;I2:(3a+4b-5)/2")
    print("dual_cauchy=F_i(X)=sum_lambda s_lambda(X)*Phi_i(s_lambda_conjugate)")
    for line in schur_lines:
        print(line)
    print(f"generic_total=slots:{totals[0]}:positive:{totals[1]}:zero:{totals[2]}:negative:0")
    print("finite_supports=" + ",".join(f"{a}-{b}" for a, b in FINITE_SUPPORTS))
    print(f"finite_schur_n5_14=checked:{finite[0]}:positive:{finite[1]}:"
          f"zero:{finite[2]}:negative:{finite[3]}")
    print("first_order_transport=" + ";".join(
        f"I{invariant + 1}-{chamber}:edges{edges}:capacity{capacity}:residual{residual}"
        for invariant, chamber, edges, capacity, residual in transport))
    print("degree_and_offgrid_guards=PASS")
    print("status=PASS")


if __name__ == "__main__":
    main()
