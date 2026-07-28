#!/usr/bin/env python3
"""Exact checks for THM-2783.

No third-party packages, floating point, or truth-bearing ``assert`` nodes.
"""

from collections import Counter, defaultdict
from itertools import combinations, product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def determinant(matrix):
    """Fraction-free Bareiss determinant."""
    n = len(matrix)
    if n == 0:
        return 1
    a = [list(row) for row in matrix]
    sign = 1
    previous = 1
    for col in range(n - 1):
        pivot = next((r for r in range(col, n) if a[r][col]), None)
        if pivot is None:
            return 0
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            sign = -sign
        p = a[col][col]
        for i in range(col + 1, n):
            for j in range(col + 1, n):
                a[i][j] = (a[i][j] * p - a[i][col] * a[col][j]) // previous
            a[i][col] = 0
        previous = p
    return sign * a[-1][-1]


def root_bank(k):
    roots = []
    for i in range(k):
        vector = [0] * k
        vector[i] = 1
        roots.append((tuple(vector), ("half", i)))
    for i in range(k):
        for j in range(i + 1, k):
            for s in (1, -1):
                vector = [0] * k
                vector[i] = 1
                vector[j] = -s
                roots.append((tuple(vector), ("edge", i, j, s)))
    return roots


def component_data(k, selected):
    parent = list(range(k))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        x, y = find(x), find(y)
        if x != y:
            parent[y] = x

    for _, meta in selected:
        if meta[0] == "edge":
            union(meta[1], meta[2])

    vertices = defaultdict(list)
    for i in range(k):
        vertices[find(i)].append(i)
    rows = defaultdict(list)
    for vector, meta in selected:
        rows[find(meta[1])].append((vector, meta))

    deficient = []
    cycles = 0
    for key, verts in vertices.items():
        defect = len(verts) - len(rows[key])
        if defect == 1:
            deficient.append((key, verts, rows[key]))
        elif defect == 0:
            edge_count = sum(meta[0] == "edge" for _, meta in rows[key])
            half_count = len(rows[key]) - edge_count
            if edge_count == len(verts) and half_count == 0:
                cycles += 1
        else:
            return None
    if len(deficient) != 1:
        return None

    _, verts, tree_rows = deficient[0]
    adjacency = defaultdict(list)
    for _, meta in tree_rows:
        require(meta[0] == "edge", "deficient component contains a half-edge")
        _, i, j, s = meta
        adjacency[i].append((j, s))
        adjacency[j].append((i, s))
    epsilon = {verts[0]: 1}
    stack = [verts[0]]
    while stack:
        i = stack.pop()
        for j, s in adjacency[i]:
            value = s * epsilon[i]
            if j in epsilon:
                require(epsilon[j] == value, "inconsistent tree switching")
            else:
                epsilon[j] = value
                stack.append(j)
    require(set(epsilon) == set(verts), "tree switching did not reach every vertex")
    return cycles, tuple(sorted(epsilon.items()))


def signed_value(delta, h):
    return sum(d * x for d, x in zip(delta, h))


def subset_sums(h):
    return {sum(bit * x for bit, x in zip(bits, h)) for bits in product((0, 1), repeat=len(h))}


def balanced_sums(h):
    return {signed_value(delta, h) for delta in product((-1, 0, 1), repeat=len(h))}


def digit_sums(h, alphabet):
    return {
        sum(digit * x for digit, x in zip(digits, h))
        for digits in product(range(alphabet), repeat=len(h))
    }


def partitions(total, length, minimum=1):
    if length == 0:
        if total == 0:
            yield ()
        return
    maximum = total // length
    for first in range(minimum, maximum + 1):
        for tail in partitions(total - first, length - 1, first):
            yield (first,) + tail


def frame_from_delta(delta):
    k = len(delta)
    support = [i for i, d in enumerate(delta) if d]
    require(support, "zero delta has no witness frame")
    rows = []
    for i, j in zip(support, support[1:]):
        s = delta[i] * delta[j]
        row = [0] * k
        row[i] = 1
        row[j] = -s
        rows.append(tuple(row))
    for j, d in enumerate(delta):
        if d == 0:
            row = [0] * k
            row[j] = 1
            rows.append(tuple(row))
    require(len(rows) == k - 1, "witness frame has wrong row count")
    return rows


def polynomial_product(exponents, alphabet):
    coefficients = [1]
    for h in exponents:
        factor = [0] * ((alphabet - 1) * h + 1)
        for digit in range(alphabet):
            factor[digit * h] = 1
        result = [0] * (len(coefficients) + len(factor) - 1)
        for i, a in enumerate(coefficients):
            for j, b in enumerate(factor):
                result[i + j] += a * b
        coefficients = result
    return coefficients


def main():
    frame_counts = {}
    formula_checks = 0
    zero_linear = 0
    tested_walls = {}

    for k in range(2, 6):
        bank = root_bank(k)
        full_rank = 0
        wall_list = (
            tuple(2**i for i in range(k)),
            tuple(3**i for i in range(k)),
            tuple(range(1, k + 1)),
            tuple(2 * i + 1 for i in range(k)),
        )
        wall_zero = Counter()
        for indices in combinations(range(len(bank)), k - 1):
            selected = [bank[i] for i in indices]
            rows = [vector for vector, _ in selected]
            gram = [[sum(rows[a][t] * rows[b][t] for t in range(k))
                     for b in range(k - 1)] for a in range(k - 1)]
            if determinant(gram) == 0:
                continue
            full_rank += 1
            data = component_data(k, selected)
            require(data is not None, f"component classification failed at k={k}")
            cycles, epsilon_pairs = data
            for h in wall_list:
                actual = abs(determinant(rows + [h]))
                imbalance = sum(h[i] * e for i, e in epsilon_pairs)
                predicted = (2**cycles) * abs(imbalance)
                require(actual == predicted, f"weighted determinant mismatch at k={k}")
                formula_checks += 1
                if actual == 0:
                    wall_zero[h] += 1
        frame_counts[k] = full_rank
        binary = wall_list[0]
        ternary = wall_list[1]
        require(wall_zero[binary] == 0, f"binary wall vanished at k={k}")
        require(wall_zero[ternary] == 0, f"ternary wall vanished at k={k}")
        zero_linear += wall_zero[wall_list[2]]
        tested_walls[k] = tuple(wall_zero[h] for h in wall_list)

    witness_checks = 0
    for k in range(1, 7):
        for delta in product((-1, 0, 1), repeat=k):
            if not any(delta):
                continue
            rows = frame_from_delta(delta)
            gram = [[sum(rows[a][t] * rows[b][t] for t in range(k))
                     for b in range(k - 1)] for a in range(k - 1)]
            require(determinant(gram) != 0, "constructed witness frame lost rank")
            for h in (tuple(2**i for i in range(k)), tuple(range(1, k + 1))):
                actual = abs(determinant(rows + [h]))
                require(actual == abs(signed_value(delta, h)),
                        "constructed witness frame has wrong determinant")
            witness_checks += 1

    binary_minima = {}
    for k in range(1, 7):
        total = 2**k - 1
        good = [h for h in partitions(total, k) if len(subset_sums(h)) == 2**k]
        require(good == [tuple(2**i for i in range(k))],
                f"binary equality classification failed at k={k}: {good}")
        require(polynomial_product(good[0], 2) == [1] * (total + 1),
                "binary radix polynomial identity failed")
        binary_minima[k] = len(good)

    ternary_minima = {}
    for k in range(1, 5):
        total = (3**k - 1) // 2
        good = [h for h in partitions(total, k) if len(balanced_sums(h)) == 3**k]
        require(good == [tuple(3**i for i in range(k))],
                f"ternary equality classification failed at k={k}: {good}")
        require(polynomial_product(good[0], 3) == [1] * (2 * total + 1),
                "ternary radix polynomial identity failed")
        ternary_minima[k] = len(good)

    radix_minima = {}
    for alphabet in range(2, 7):
        for k in range(1, 4):
            total = (alphabet**k - 1) // (alphabet - 1)
            good = [
                h for h in partitions(total, k)
                if len(digit_sums(h, alphabet)) == alphabet**k
            ]
            expected = tuple(alphabet**i for i in range(k))
            require(good == [expected],
                    f"radix equality classification failed at q={alphabet}, k={k}: {good}")
            require(polynomial_product(expected, alphabet) ==
                    [1] * ((alphabet - 1) * total + 1),
                    "general radix polynomial identity failed")
            radix_minima[(alphabet, k)] = len(good)

    require(len(subset_sums((1, 2, 4))) == 8, "binary positive control failed")
    require(len(balanced_sums((1, 2, 4))) < 27,
            "binary state-reconstruction hostile did not collide")
    require(signed_value((1, 1, -1), (1, 2, 3)) == 0,
            "linear-wall null hostile failed")
    hostile_rows = frame_from_delta((1, 1, -1))
    require(determinant(hostile_rows + [(1, 2, 3)]) == 0,
            "linear-wall determinant hostile failed")

    print("THM-2783 WEIGHTED LONG-WALL CODING AUDIT")
    print(f"full_rank_frames_k2_to_k5={frame_counts}")
    print(f"weighted_determinant_formula_checks={formula_checks}")
    print(f"wall_zero_counts_binary_ternary_linear_odd={tested_walls}")
    print(f"constructed_signed_state_witnesses_k1_to_k6={witness_checks}")
    print(f"binary_minimum_sum_uniqueness_k1_to_k6={binary_minima}")
    print(f"ternary_minimum_sum_uniqueness_k1_to_k4={ternary_minima}")
    print(f"general_radix_uniqueness_q2_to_q6_k1_to_k3={radix_minima}")
    print("binary_wall=(1,2,4): regular_but_signed_states_alias")
    print("linear_wall=(1,2,3): delta=(1,1,-1), determinant=0")
    print(f"linear_zero_frames_total_k2_to_k5={zero_linear}")
    print("normal_form=binary_null_avoidance/ternary_full_state_reconstruction")
    print("PASS")


if __name__ == "__main__":
    main()
