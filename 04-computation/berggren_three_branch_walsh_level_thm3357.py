#!/usr/bin/env python3
"""Exact companion for THM-3357.

The program checks the simultaneous three-branch consequence objects rather
than fitting scalar prefixes.  It uses independent explicit tree enumeration
and matrix powering, audits the local determinant circuit, the weighted
transfer polynomial, the four Walsh specializations, branch-parity class
sums, the Pell/Markov mass compiler, the prime-five sibling clock, the level
Lorentz second moment and even/odd parameter-square split, the factorial
moments of the triple-transfer determinant, and the hostiles showing what the
aggregate state cannot recover.
"""

from __future__ import annotations

from itertools import combinations, product
from math import factorial, gcd


Vec2 = tuple[int, int]
Vec3 = tuple[int, int, int]
Mat2 = tuple[tuple[int, int], tuple[int, int]]
Mat3 = tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


I2: Mat2 = ((1, 0), (0, 1))
L: Mat2 = ((0, 1), (-1, 2))
M: Mat2 = ((0, 1), (1, 2))
R: Mat2 = ((1, 0), (2, 1))
H: Mat2 = ((-1, 1), (1, 1))
BRANCHES: tuple[tuple[str, Mat2], ...] = (("L", L), ("M", M), ("R", R))
ROOT: Vec2 = (1, 2)

TL: Mat3 = ((1, -2, 2), (2, -1, 2), (2, -2, 3))
TM: Mat3 = ((1, 2, 2), (2, 1, 2), (2, 2, 3))
TR: Mat3 = ((-1, 2, 2), (-2, 1, 2), (-2, 2, 3))
TRIPLE_BRANCHES: tuple[Mat3, ...] = (TL, TM, TR)


def v2_add(*vectors: Vec2) -> Vec2:
    return sum(v[0] for v in vectors), sum(v[1] for v in vectors)


def v2_scale(k: int, v: Vec2) -> Vec2:
    return k * v[0], k * v[1]


def m2_add(*matrices: Mat2) -> Mat2:
    return tuple(
        tuple(sum(matrix[i][j] for matrix in matrices) for j in range(2))
        for i in range(2)
    )  # type: ignore[return-value]


def m2_scale(k: int, matrix: Mat2) -> Mat2:
    return tuple(
        tuple(k * matrix[i][j] for j in range(2)) for i in range(2)
    )  # type: ignore[return-value]


def m2_mul(a: Mat2, b: Mat2) -> Mat2:
    return tuple(
        tuple(sum(a[i][k] * b[k][j] for k in range(2)) for j in range(2))
        for i in range(2)
    )  # type: ignore[return-value]


def m2_pow(matrix: Mat2, exponent: int) -> Mat2:
    result = I2
    base = matrix
    while exponent:
        if exponent & 1:
            result = m2_mul(result, base)
        base = m2_mul(base, base)
        exponent //= 2
    return result


def m2_apply(matrix: Mat2, vector: Vec2) -> Vec2:
    return tuple(
        sum(matrix[i][j] * vector[j] for j in range(2)) for i in range(2)
    )  # type: ignore[return-value]


def m2_det(matrix: Mat2) -> int:
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def m2_inv_unimodular(matrix: Mat2) -> Mat2:
    determinant = m2_det(matrix)
    require(abs(determinant) == 1, "conjugating matrix must be unimodular")
    return (
        (matrix[1][1] // determinant, -matrix[0][1] // determinant),
        (-matrix[1][0] // determinant, matrix[0][0] // determinant),
    )


def m3_add(*matrices: Mat3) -> Mat3:
    return tuple(
        tuple(sum(matrix[i][j] for matrix in matrices) for j in range(3))
        for i in range(3)
    )  # type: ignore[return-value]


def m3_scale(k: int, matrix: Mat3) -> Mat3:
    return tuple(
        tuple(k * matrix[i][j] for j in range(3)) for i in range(3)
    )  # type: ignore[return-value]


def m3_mul(a: Mat3, b: Mat3) -> Mat3:
    return tuple(
        tuple(sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3))
        for i in range(3)
    )  # type: ignore[return-value]


def m3_pow(matrix: Mat3, exponent: int) -> Mat3:
    result: Mat3 = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    base = matrix
    while exponent:
        if exponent & 1:
            result = m3_mul(result, base)
        base = m3_mul(base, base)
        exponent //= 2
    return result


def m3_apply(matrix: Mat3, vector: Vec3) -> Vec3:
    return tuple(
        sum(matrix[i][j] * vector[j] for j in range(3)) for i in range(3)
    )  # type: ignore[return-value]


def m3_det(matrix: Mat3) -> int:
    a, b, c = matrix[0]
    d, e, f = matrix[1]
    g, h, i = matrix[2]
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def det2(u: Vec2, v: Vec2) -> int:
    return u[0] * v[1] - u[1] * v[0]


def det3_columns(u: Vec3, v: Vec3, w: Vec3) -> int:
    return m3_det(
        ((u[0], v[0], w[0]), (u[1], v[1], w[1]), (u[2], v[2], w[2]))
    )


def psi(u: Vec2) -> Vec3:
    m, n = u
    return n * n - m * m, 2 * m * n, n * n + m * m


def lorentz(x: Vec3, y: Vec3) -> int:
    return x[2] * y[2] - x[0] * y[0] - x[1] * y[1]


def v3_sum(vectors: list[Vec3]) -> Vec3:
    return tuple(sum(v[i] for v in vectors) for i in range(3))  # type: ignore[return-value]


def level(root: Vec2, depth: int) -> list[tuple[str, Vec2]]:
    current = [("", root)]
    for _ in range(depth):
        current = [
            (word + name, m2_apply(matrix, vector))
            for word, vector in current
            for name, matrix in BRANCHES
        ]
    return current


def level_sum(root: Vec2, depth: int) -> Vec2:
    nodes = level(root, depth)
    return v2_add(*(vector for _, vector in nodes))


def weighted_level_sum(root: Vec2, depth: int, weights: dict[str, int]) -> Vec2:
    total = (0, 0)
    for word, vector in level(root, depth):
        weight = 1
        for letter in word:
            weight *= weights[letter]
        total = v2_add(total, v2_scale(weight, vector))
    return total


def pell(index: int) -> int:
    a, b = 0, 1
    for _ in range(index):
        a, b = b, 2 * b + a
    return a


def triangular(index: int) -> int:
    return index * (index + 1) // 2


def determinant_energy(vectors: tuple[Vec2, ...] | list[Vec2]) -> int:
    return sum(det2(u, v) ** 2 for u, v in combinations(vectors, 2))


def determinant_pattern(vectors: tuple[Vec2, ...]) -> tuple[int, ...]:
    return tuple(sorted(abs(det2(u, v)) for u, v in combinations(vectors, 2)))


def quadratic_state(vectors: tuple[Vec2, ...]) -> tuple[int, int, int]:
    return (
        sum(m * m for m, _ in vectors),
        sum(m * n for m, n in vectors),
        sum(n * n for _, n in vectors),
    )


def audit_local_circuit() -> int:
    parents = 0
    for m in range(1, 61):
        for n in range(m + 1, 101):
            if gcd(m, n) != 1 or (m + n) % 2 != 1:
                continue
            parent = (m, n)
            a, b, c = psi(parent)
            left, middle, right = (m2_apply(matrix, parent) for _, matrix in BRANCHES)
            require(det2(left, middle) == b, "left/middle determinant")
            require(det2(left, right) == c, "left/right determinant")
            require(det2(middle, right) == a, "middle/right determinant")
            require(v2_add(v2_scale(a, left), v2_scale(b, right)) == v2_scale(c, middle),
                    "parent circuit")
            require(left[1] * middle[0] < middle[1] * left[0], "left slope order")
            require(middle[1] * right[0] < right[1] * middle[0], "right slope order")
            triples = tuple(psi(child) for child in (left, middle, right))
            norms = tuple(x * x + y * y for x, y in (left, middle, right))
            require(norms[1] - norms[0] == 4 * b, "middle/left norm gap")
            require(norms[1] - norms[2] == 4 * a, "middle/right norm gap")
            require(norms[0] + norms[2] == 6 * c, "outer norm sum")
            require(norms[2] - norms[0] == 4 * (b - a), "outer norm arc")
            require(len(set(norms)) == 3 and norms[1] == max(norms),
                    "intrinsic transitive norm tournament")
            child_leg_signs = tuple(
                (triple[1] > triple[0]) - (triple[1] < triple[0])
                for triple in triples
            )
            parent_leg_sign = (b > a) - (b < a)
            require(child_leg_signs == (1, -parent_leg_sign, -1),
                    "reset/flip/reset leg-sign automaton")
            require(lorentz(triples[0], triples[1]) == 2 * b * b, "LM shell")
            require(lorentz(triples[0], triples[2]) == 2 * c * c, "LR shell")
            require(lorentz(triples[1], triples[2]) == 2 * a * a, "MR shell")
            require(det3_columns(*triples) == -4 * a * b * c, "Lorentz volume")
            require(determinant_energy((left, middle, right)) == 2 * c * c,
                    "local energy")
            parents += 1
    return parents


def determinant_gauge(vector: Vec2, deck: tuple[Vec2, ...]) -> int:
    return max(abs(det2(vector, column)) for column in deck)


def audit_convex_gate_horn() -> int:
    """Check the owner-free three-sibling determinant-gate implication."""
    deck: tuple[Vec2, ...] = ((1, 0), (0, 1), (1, 1), (2, -1), (3, 2))
    checks = 0
    for m in range(1, 31):
        for n in range(m + 1, 51):
            if gcd(m, n) != 1 or (m + n) % 2 != 1:
                continue
            parent = (m, n)
            a, b, c = psi(parent)
            left, middle, right = (m2_apply(matrix, parent) for _, matrix in BRANCHES)
            norms = tuple(x * x + y * y for x, y in (left, middle, right))
            correction = c * norms[1] - a * norms[0] - b * norms[2]
            require(
                correction == 2 * m * (n - m) * (3 * n * n + 4 * m * n - m * m),
                "strict Horn correction",
            )
            require(correction > 0, "strict Horn positivity")
            for kappa in (0, 1, 7, 91, 137):
                gauge = tuple(determinant_gauge(v, deck) for v in (left, middle, right))
                score = tuple(norms[i] - kappa * gauge[i] for i in range(3))
                require(
                    c * score[1] >= a * score[0] + b * score[2] + correction,
                    "convex determinant-gate Horn inequality",
                )
                checks += 1
    return checks


def audit_prime_five_sibling_clock() -> list[int]:
    infinity = None

    def inverse_mod_five(state: int | None) -> int | None:
        if state is None:
            return 0
        if state == 0:
            return infinity
        return pow(state, -1, 5)

    def step(name: str, state: int | None) -> int | None:
        inverse = inverse_mod_five(state)
        if name == "R":
            return None if state is None else (state + 2) % 5
        if inverse is None:
            return infinity
        return (2 + (-inverse if name == "L" else inverse)) % 5

    classes: tuple[tuple[int | None, ...], ...] = ((None, 1), (0, 4), (2, 3))
    transition_counts: list[list[int]] = []
    for source_class in classes:
        source_patterns: list[list[int]] = []
        for source in source_class:
            totals = [0, 0, 0]
            for name, _ in BRANCHES:
                target = step(name, source)
                totals[next(i for i, block in enumerate(classes) if target in block)] += 1
            source_patterns.append(totals)
        require(source_patterns[0] == source_patterns[1], "equitable class total")
        transition_counts.append(source_patterns[0])
    # Columns are source classes and rows are target classes.
    quotient = [[transition_counts[column][row] for column in range(3)]
                for row in range(3)]
    require(quotient == [[1, 2, 0], [0, 0, 3], [2, 1, 0]],
            "prime-five equitable quotient")
    for state in (None, 0, 1, 2, 3, 4):
        if state is not None:
            require(((state * state + 1) % 5 == 0) == (state in classes[2]),
                    "parent norm-five class")

    counts: list[int] = []
    current = [ROOT]
    for depth in range(11):
        collision_count = 0
        for parent in current:
            m, n = parent
            left, middle, right = (m2_apply(matrix, parent) for _, matrix in BRANCHES)
            norms = tuple(x * x + y * y for x, y in (left, middle, right))
            require(gcd(norms[0], norms[1]) == (5 if m % 5 == 0 else 1),
                    "left/middle prime-five gcd")
            require(gcd(norms[1], norms[2]) == (5 if (n - m) % 5 == 0 else 1),
                    "middle/right prime-five gcd")
            require(gcd(norms[0], norms[2]) == 1, "outer sibling coprimality")
            if gcd(norms[0], norms[1]) == 5 or gcd(norms[1], norms[2]) == 5:
                collision_count += 1
        counts.append(collision_count)
        current = [m2_apply(matrix, parent) for parent in current for _, matrix in BRANCHES]
    require(counts[:8] == [0, 0, 6, 6, 24, 96, 222, 726],
            "prime-five collision prefix")
    for depth in range(8):
        require(counts[depth + 3] == counts[depth + 2] + 3 * counts[depth + 1]
                + 9 * counts[depth], "prime-five collision recurrence")
    return counts


def audit_outer_hadamard() -> int:
    require(m2_mul(H, H) == m2_scale(2, I2), "Hadamard square")
    require(m2_mul(H, L) == m2_mul(R, H), "Hadamard L/R conjugacy")
    require(m2_mul(H, M) == m2_mul(M, H), "Hadamard middle fixed")
    require(m2_mul(H, R) == m2_mul(L, H), "Hadamard R/L conjugacy")
    checked = 0
    for m in range(1, 61):
        for n in range(m + 1, 101):
            if gcd(m, n) != 1 or (m + n) % 2 != 1:
                continue
            a, b, c = psi((m, n))
            transformed = m2_apply(H, (m, n))
            require(gcd(*transformed) == 1 and sum(transformed) % 2 == 0,
                    "Hadamard primitive odd/odd chamber")
            require(psi(transformed) == (2 * b, 2 * a, 2 * c),
                    "Hadamard leg swap/content two")
            checked += 1
    return checked


def audit_weighted_transfer() -> int:
    weights_checked = 0
    for x, y, z in product(range(-2, 3), repeat=3):
        transfer = m2_add(m2_scale(x, L), m2_scale(y, M), m2_scale(z, R))
        require(transfer[0] == (z, x + y), "weighted first row")
        require(transfer[1] == (-x + y + 2 * z, 2 * x + 2 * y + z),
                "weighted second row")
        sigma = x + y + z
        delta = x * x - y * y + z * z
        require(transfer[0][0] + transfer[1][1] == 2 * sigma, "weighted trace")
        require(m2_det(transfer) == delta, "weighted determinant")
        lhs = m2_mul(transfer, transfer)
        rhs = m2_add(m2_scale(2 * sigma, transfer), m2_scale(-delta, I2))
        require(lhs == rhs, "weighted Cayley-Hamilton")
        if (x, y, z) in ((1, 1, 1), (1, 1, -1), (-1, 1, 1), (1, -1, 1),
                         (2, -1, 3), (-2, 1, 0)):
            for depth in range(6):
                explicit = weighted_level_sum(ROOT, depth, {"L": x, "M": y, "R": z})
                powered = m2_apply(m2_pow(transfer, depth), ROOT)
                require(explicit == powered, "weighted enumeration/power mismatch")
        weights_checked += 1
    return weights_checked


def audit_walsh_and_parity() -> list[tuple[int, Vec2, tuple[int, ...]]]:
    require(m2_add(L, M, R) == m2_pow(M, 2), "unsigned collapse")
    require(m2_add(L, M, m2_scale(-1, R)) == m2_pow(L, 2), "R-sign collapse")
    require(m2_add(m2_scale(-1, L), M, R) == m2_pow(R, 2), "L-sign collapse")
    require(m2_add(L, m2_scale(-1, M), R) == I2, "M-sign conservation")
    require(m2_pow(M, 2) == m2_add(m2_scale(2, M), I2), "middle silver law")
    require(m2_pow(L, 2) == m2_add(m2_scale(2, L), m2_scale(-1, I2)),
            "left parabolic law")
    require(m2_pow(R, 2) == m2_add(m2_scale(2, R), m2_scale(-1, I2)),
            "right parabolic law")
    require(m2_add(L, R) == m2_add(M, I2), "quartet additive certificate")

    rows: list[tuple[int, Vec2, tuple[int, ...]]] = []
    for depth in range(9):
        nodes = level(ROOT, depth)
        require(len(nodes) == 3 ** depth, "level node count")
        require(len({vector for _, vector in nodes}) == len(nodes), "tree collision")
        transforms = {
            "all": weighted_level_sum(ROOT, depth, {"L": 1, "M": 1, "R": 1}),
            "R": weighted_level_sum(ROOT, depth, {"L": 1, "M": 1, "R": -1}),
            "L": weighted_level_sum(ROOT, depth, {"L": -1, "M": 1, "R": 1}),
            "M": weighted_level_sum(ROOT, depth, {"L": 1, "M": -1, "R": 1}),
        }
        require(transforms["all"] == m2_apply(m2_pow(M, 2 * depth), ROOT),
                "all-level unary collapse")
        require(transforms["R"] == m2_apply(m2_pow(L, 2 * depth), ROOT),
                "R-parity unary collapse")
        require(transforms["L"] == m2_apply(m2_pow(R, 2 * depth), ROOT),
                "L-parity unary collapse")
        require(transforms["M"] == ROOT, "M-parity root conservation")

        grouped: dict[tuple[int, int, int], tuple[int, Vec2]] = {}
        for word, vector in nodes:
            parity = tuple(word.count(letter) % 2 for letter in "LMR")
            count, total = grouped.get(parity, (0, (0, 0)))
            grouped[parity] = count + 1, v2_add(total, vector)
        allowed = [p for p in product((0, 1), repeat=3) if sum(p) % 2 == depth % 2]
        require(set(grouped).issubset(set(allowed)), "parity support")
        for parity in allowed:
            grouped.setdefault(parity, (0, (0, 0)))
        for parity in allowed:
            p_l, p_m, p_r = parity
            formula = v2_scale(
                1,
                v2_add(
                    transforms["all"],
                    v2_scale((-1) ** p_r, transforms["R"]),
                    v2_scale((-1) ** p_l, transforms["L"]),
                    v2_scale((-1) ** p_m, transforms["M"]),
                ),
            )
            require(formula[0] % 4 == formula[1] % 4 == 0, "Walsh integrality")
            formula = formula[0] // 4, formula[1] // 4
            require(grouped[parity][1] == formula, "parity-class sum")
            expected_count = (
                3 ** depth + (-1) ** p_l + (-1) ** p_m + (-1) ** p_r
            ) // 4
            require(grouped[parity][0] == expected_count, "parity-class count")
        rows.append((depth, transforms["all"], tuple(sorted(count for count, _ in grouped.values()))))
    return rows


def audit_conjugation_and_boundary() -> tuple[Mat2, Mat2]:
    for change in (((1, 1), (0, 1)), ((0, 1), (-1, 0)), ((2, 1), (1, 1))):
        inverse = m2_inv_unimodular(change)
        branches = tuple(m2_mul(m2_mul(change, matrix), inverse) for matrix in (L, M, R))
        left, middle, right = branches
        require(m2_add(left, middle, right) == m2_pow(middle, 2), "conjugate total")
        require(m2_add(left, middle, m2_scale(-1, right)) == m2_pow(left, 2),
                "conjugate R-sign")
        require(m2_add(m2_scale(-1, left), middle, right) == m2_pow(right, 2),
                "conjugate L-sign")
        require(m2_add(left, m2_scale(-1, middle), right) == I2,
                "conjugate M-sign")

    bad_right: Mat2 = ((1, 1), (2, 3))
    bad_sum = m2_add(L, M, bad_right)
    require(m2_det(bad_right) == 1, "bad control is unimodular")
    require(bad_sum != m2_pow(M, 2), "arbitrary ternary collapse hostile")
    require(m2_add(L, m2_scale(-1, M), bad_right) != I2,
            "arbitrary ternary conservation hostile")
    return bad_right, bad_sum


def audit_pell_markov(rows: list[tuple[int, Vec2, tuple[int, ...]]]) -> list[tuple[int, Vec2, int, int]]:
    compiled: list[tuple[int, Vec2, int, int]] = []
    for depth in range(11):
        total = m2_apply(m2_pow(M, 2 * depth), ROOT)
        next_total = m2_apply(m2_pow(M, 2 * (depth + 1)), ROOT)
        lower, carry = total
        upper = next_total[0]
        require(total == (pell(2 * depth + 1), pell(2 * depth + 2)), "Pell total")
        require(upper == pell(2 * depth + 3), "Pell next lower")
        require(lower * upper == carry * carry + 1, "Cassini/Markov carry")
        require(lower * lower + upper * upper + 4 == 6 * lower * upper,
                "fixed-two Markov equation")
        require(carry % 2 == 0 and (lower + upper - 2) % 4 == 0,
                "selector parity")
        index = (lower + upper - 2) // 4
        root = carry // 2
        require(triangular(index) == root * root, "square-triangular compiler")
        if depth < 6:
            compiled.append((depth, total, upper, index))
    require(rows[2][1] == (29, 70), "depth-two cannonball total")
    require(sum(r * r for r in range(1, 25)) == rows[2][1][1] ** 2,
            "cannonball identity")
    require(
        sum((2 * r) ** 2 for r in range(1, 13))
        + sum((2 * r - 1) ** 2 for r in range(1, 13))
        == rows[2][1][1] ** 2,
        "even/odd cannonball split",
    )
    return compiled


def audit_triple_levels() -> tuple[list[tuple[int, Vec3, int]], list[int]]:
    transfer = m3_add(TL, TM, TR)
    require(transfer == ((1, 2, 6), (2, 1, 6), (2, 2, 9)), "triple transfer")
    require(m3_det(transfer) == -3, "triple transfer determinant")

    for x, y, z in product(range(-2, 3), repeat=3):
        weighted = m3_add(m3_scale(x, TL), m3_scale(y, TM), m3_scale(z, TR))
        expected = (x - y - z) * (x + y - z) * (x + y + z)
        require(m3_det(weighted) == expected, "weighted triple determinant")

    level_rows: list[tuple[int, Vec3, int]] = []
    energies: list[int] = []
    for depth in range(13):
        aggregate = m3_apply(m3_pow(transfer, depth), psi(ROOT))
        if depth <= 7:
            nodes = [vector for _, vector in level(ROOT, depth)]
            explicit = v3_sum([psi(vector) for vector in nodes])
            require(explicit == aggregate, "triple enumeration/power mismatch")
            even_square_total = sum(
                (m if m % 2 == 0 else n) ** 2 for m, n in nodes
            )
            odd_square_total = sum(
                (m if m % 2 == 1 else n) ** 2 for m, n in nodes
            )
            signed_difference = (-1) ** depth * (4 * 3 ** depth - 1)
            require(even_square_total + odd_square_total == aggregate[2],
                    "Euclid-parameter square split total")
            require(even_square_total - odd_square_total == signed_difference,
                    "Euclid-parameter even/odd square split")
        a_sum, b_sum, c_sum = aggregate
        require(a_sum - b_sum == (-1) ** (depth + 1), "alternating leg total")
        energy = (c_sum * c_sum - a_sum * a_sum - b_sum * b_sum) // 4
        require(4 * energy == c_sum * c_sum - a_sum * a_sum - b_sum * b_sum,
                "energy divisibility")
        if depth <= 6:
            nodes = [vector for _, vector in level(ROOT, depth)]
            require(determinant_energy(nodes) == energy, "Cauchy-Binet energy")
        level_rows.append((depth, aggregate, energy))
        energies.append(energy)

    for depth in range(10):
        for coordinate in (0, 1, 2):
            values = [row[1][coordinate] for row in level_rows]
            require(
                values[depth + 3]
                == 11 * values[depth + 2] + 9 * values[depth + 1] - 3 * values[depth],
                "triple order-three recurrence",
            )
    leg_sums = [a + b for _, (a, b, _), _ in level_rows]
    hyp_sums = [c for _, (_, _, c), _ in level_rows]
    for values in (leg_sums, hyp_sums):
        for depth in range(11):
            require(values[depth + 2] == 12 * values[depth + 1] - 3 * values[depth],
                    "order-two aggregate recurrence")
    for depth in range(9):
        require(
            energies[depth + 4]
            == 142 * energies[depth + 3]
            - 564 * energies[depth + 2]
            + 450 * energies[depth + 1]
            - 27 * energies[depth],
            "energy recurrence",
        )
    require(energies[:4] == [0, 50, 7012, 967510], "energy seeds")
    return level_rows[:6], energies[:6]


Poly3 = dict[tuple[int, int, int], int]


def poly3_mul(first: Poly3, second: Poly3) -> Poly3:
    result: Poly3 = {}
    for (i, j, k), coefficient in first.items():
        for (a, b, c), other in second.items():
            exponent = i + a, j + b, k + c
            result[exponent] = result.get(exponent, 0) + coefficient * other
    return {exponent: coefficient for exponent, coefficient in result.items() if coefficient}


def factorial_functional(poly: Poly3) -> int:
    return sum(coefficient * factorial(i) * factorial(j) * factorial(k)
               for (i, j, k), coefficient in poly.items())


def audit_factorial_transfer_determinant() -> list[int]:
    # (x-y-z)(x+y-z)(x+y+z)
    linear_factors: tuple[Poly3, ...] = (
        {(1, 0, 0): 1, (0, 1, 0): -1, (0, 0, 1): -1},
        {(1, 0, 0): 1, (0, 1, 0): 1, (0, 0, 1): -1},
        {(1, 0, 0): 1, (0, 1, 0): 1, (0, 0, 1): 1},
    )
    determinant_poly: Poly3 = {(0, 0, 0): 1}
    for factor in linear_factors:
        determinant_poly = poly3_mul(determinant_poly, factor)
    power: Poly3 = {(0, 0, 0): 1}
    values: list[int] = []
    for exponent in range(9):
        value = factorial_functional(power)
        expected = 0 if exponent % 2 else factorial(3 * exponent + 2) // (
            2 * (exponent + 1) ** 2
        )
        require(value == expected, "factorial determinant moment")
        values.append(value)
        power = poly3_mul(power, determinant_poly)
    require(values[:5] == [1, 0, 2240, 0, 1743565824],
            "factorial determinant controls")
    return values


def audit_hostiles() -> tuple[Vec2, Vec3, Vec3, tuple[object, ...]]:
    children = tuple(m2_apply(matrix, ROOT) for _, matrix in BRANCHES)
    parameter_sum = v2_add(*children)
    triple_sum = v3_sum([psi(child) for child in children])
    nonlinear_image = psi(parameter_sum)
    require(parameter_sum == (5, 12), "nonlinear hostile parameter sum")
    require(triple_sum == (41, 40, 59), "nonlinear hostile triple sum")
    require(nonlinear_image == (119, 120, 169), "nonlinear hostile image")
    require(triple_sum != nonlinear_image, "Gaussian square must not linearize")

    deck_first: tuple[Vec2, ...] = ((1, 4), (2, 3), (2, 9))
    deck_second: tuple[Vec2, ...] = ((1, 6), (2, 3), (2, 7))
    for packet in (deck_first, deck_second):
        require(all(gcd(*u) == 1 and sum(u) % 2 == 1 for u in packet),
                "hostile primitive parameters")
    first_sum, second_sum = v2_add(*deck_first), v2_add(*deck_second)
    first_energy = determinant_energy(deck_first)
    second_energy = determinant_energy(deck_second)
    first_pattern = determinant_pattern(deck_first)
    second_pattern = determinant_pattern(deck_second)
    require(first_sum == second_sum == (5, 16), "same first moment hostile")
    require(first_energy == second_energy == 170, "same energy hostile")
    require(first_pattern == (1, 5, 12) and second_pattern == (5, 8, 9),
            "determinant patterns")
    require(max(first_pattern) != max(second_pattern), "maximum determinant loss")

    depth_three = dict(level(ROOT, 3))
    moment_first = tuple(depth_three[word] for word in ("LLM", "MLR", "MRL"))
    moment_second = tuple(depth_three[word] for word in ("LMM", "LMR", "LRL"))
    require(v2_add(*moment_first) == v2_add(*moment_second) == (18, 45),
            "full-moment hostile first moment")
    require(quadratic_state(moment_first) == quadratic_state(moment_second)
            == (122, 278, 701), "full-moment hostile quadratic state")
    require(determinant_energy(moment_first) == determinant_energy(moment_second) == 8238,
            "full-moment hostile energy")
    moment_patterns = determinant_pattern(moment_first), determinant_pattern(moment_second)
    require(moment_patterns == ((17, 35, 82), (37, 55, 62)),
            "full-moment hostile determinant patterns")
    lm = m2_apply(M, m2_apply(L, ROOT))
    ml = m2_apply(L, m2_apply(M, ROOT))
    require(lm == (3, 8) and ml == (5, 8) and lm != ml,
            "branch-count aggregate loses ancestry order")
    return parameter_sum, triple_sum, nonlinear_image, (
        first_sum, first_energy, first_pattern, second_pattern,
        quadratic_state(moment_first), moment_patterns, lm, ml
    )


def main() -> None:
    parents = audit_local_circuit()
    horn_checks = audit_convex_gate_horn()
    hadamard = audit_outer_hadamard()
    weights = audit_weighted_transfer()
    walsh_rows = audit_walsh_and_parity()
    bad_right, bad_sum = audit_conjugation_and_boundary()
    compiled = audit_pell_markov(walsh_rows)
    prime_five_counts = audit_prime_five_sibling_clock()
    triple_rows, energy_rows = audit_triple_levels()
    factorial_values = audit_factorial_transfer_determinant()
    nonlinear = audit_hostiles()

    print("THM-3357 BERGGREN THREE-BRANCH WALSH / LEVEL AUDIT")
    print(f"LOCAL_CIRCUIT parents={parents} PASS")
    print(f"CONVEX_GATE_HORN checks={horn_checks} kappa_includes_91 PASS")
    print(f"OUTER_HADAMARD parents={hadamard} swaps=L/R fixes=M raw_content=2 PASS")
    print(f"WEIGHTED_TRANSFER integer_weight_triples={weights} PASS")
    print("WALSH_IDENTITIES all=M^2 Rsign=L^2 Lsign=R^2 Msign=I PASS")
    for depth, total, counts in walsh_rows[:7]:
        print(f"LEVEL depth={depth} nodes={3 ** depth} spinor_sum={total} parity_counts={counts}")
    print("PELL_MARKOV_MASS")
    for depth, total, upper, index in compiled:
        print(f"  depth={depth} total={total} next_first={upper} square_triangular_index={index}")
    print("CANNONBALL depth=2 second_mass=70 even_square_sum=2600 odd_square_sum=2300")
    print(f"PRIME_FIVE_SIBLING_COLLISIONS counts={prime_five_counts[:8]} PASS")
    print("TRIPLE_LEVELS")
    for depth, aggregate, energy in triple_rows:
        print(f"  depth={depth} triple_sum={aggregate} determinant_energy={energy}")
    print(f"ENERGY_PREFIX {energy_rows}")
    print(f"FACTORIAL_TRIPLE_DETERMINANT prefix={factorial_values[:5]} PASS")
    print(f"NONLINEAR_HOSTILE parameter_sum={nonlinear[0]} triple_sum={nonlinear[1]} psi_sum={nonlinear[2]}")
    print(f"MOMENT_MAX_ORDER_HOSTILES state={nonlinear[3]}")
    print(f"CONJUGATION_COVARIANCE PASS arbitrary_unimodular_hostile_R={bad_right} total={bad_sum}")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
