"""Independent exact audit for the rational-edge/p-adic/C6 candidate.

This scratch verifier deliberately uses only the Python standard library.  It
does not import SymPy or any function from the candidate script.  Its finite
checks support, but do not replace, the symbolic proofs in REPORT.md.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, permutations
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def identity(n: int) -> list[list[Fraction]]:
    return [
        [Fraction(int(i == j)) for j in range(n)]
        for i in range(n)
    ]


def matmul(
    left: list[list[Fraction]], right: list[list[Fraction]]
) -> list[list[Fraction]]:
    rows = len(left)
    middle = len(right)
    columns = len(right[0])
    return [
        [
            sum((left[i][k] * right[k][j] for k in range(middle)), Fraction(0))
            for j in range(columns)
        ]
        for i in range(rows)
    ]


def determinant_fraction(matrix: list[list[Fraction]]) -> Fraction:
    work = [[Fraction(value) for value in row] for row in matrix]
    n = len(work)
    answer = Fraction(1)
    for column in range(n):
        pivot_row = next(
            (row for row in range(column, n) if work[row][column] != 0),
            None,
        )
        if pivot_row is None:
            return Fraction(0)
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            answer = -answer
        pivot = work[column][column]
        answer *= pivot
        for row in range(column + 1, n):
            if work[row][column] == 0:
                continue
            factor = work[row][column] / pivot
            for j in range(column + 1, n):
                work[row][j] -= factor * work[column][j]
            work[row][column] = Fraction(0)
    return answer


def determinant_bareiss(matrix: list[list[int]]) -> int:
    """Fraction-free exact determinant, with honest row pivoting."""

    work = [row[:] for row in matrix]
    n = len(work)
    if n == 1:
        return work[0][0]
    sign = 1
    previous = 1
    for column in range(n - 1):
        pivot_row = next(
            (row for row in range(column, n) if work[row][column] != 0),
            None,
        )
        if pivot_row is None:
            return 0
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            sign = -sign
        pivot = work[column][column]
        for row in range(column + 1, n):
            for j in range(column + 1, n):
                numerator = work[row][j] * pivot - work[row][column] * work[column][j]
                require(numerator % previous == 0, "Bareiss division lost exactness")
                work[row][j] = numerator // previous
            work[row][column] = 0
        previous = pivot
    return sign * work[-1][-1]


def gauge_matrix(
    matrix: list[list[int | Fraction]], weights: tuple[Fraction, ...]
) -> list[list[Fraction]]:
    if any(weight == 0 for weight in weights):
        raise ValueError("vertex weights must be units/nonzero in the field")
    n = len(matrix)
    require(len(weights) == n, "weight count")
    return [
        [Fraction(matrix[i][j]) * weights[i] / weights[j] for j in range(n)]
        for i in range(n)
    ]


def i_minus_scalar_matrix(
    matrix: list[list[int | Fraction]], scalar: Fraction
) -> list[list[Fraction]]:
    n = len(matrix)
    return [
        [Fraction(int(i == j)) - scalar * Fraction(matrix[i][j]) for j in range(n)]
        for i in range(n)
    ]


def audit_gauge() -> tuple[int, int]:
    """All 3-vertex digraphs, including loops and digons, over Q."""

    weight_rows = (
        (Fraction(1), Fraction(2), Fraction(3)),
        (Fraction(-2), Fraction(3), Fraction(5)),
        (Fraction(1, 2), Fraction(-3, 4), Fraction(5, 7)),
    )
    scalar_rows = (Fraction(-2), Fraction(1, 3), Fraction(2), Fraction(5, 2))
    determinant_gates = 0
    power_gates = 0
    for code in range(1 << 9):
        adjacency = [
            [(code >> (3 * i + j)) & 1 for j in range(3)]
            for i in range(3)
        ]
        adjacency_q = [[Fraction(value) for value in row] for row in adjacency]
        for weights in weight_rows:
            weighted = gauge_matrix(adjacency, weights)
            for scalar in scalar_rows:
                ordinary_det = determinant_fraction(
                    i_minus_scalar_matrix(adjacency, scalar)
                )
                weighted_det = determinant_fraction(
                    i_minus_scalar_matrix(weighted, scalar)
                )
                require(ordinary_det == weighted_det, "digraph determinant gauge")
                determinant_gates += 1

            ordinary_power = identity(3)
            weighted_power = identity(3)
            for _power in range(1, 6):
                ordinary_power = matmul(ordinary_power, adjacency_q)
                weighted_power = matmul(weighted_power, weighted)
                expected = gauge_matrix(ordinary_power, weights)
                require(weighted_power == expected, "power similarity gauge")
                require(
                    sum(weighted_power[i][i] for i in range(3))
                    == sum(ordinary_power[i][i] for i in range(3)),
                    "trace gauge",
                )
                power_gates += 1

    multigraph = [[2, 1, 0], [3, 0, 2], [1, 4, 1]]
    multigraph_weights = (Fraction(-7, 3), Fraction(5, 2), Fraction(11, 13))
    multigraph_gauge = gauge_matrix(multigraph, multigraph_weights)
    for scalar in scalar_rows:
        require(
            determinant_fraction(i_minus_scalar_matrix(multigraph, scalar))
            == determinant_fraction(i_minus_scalar_matrix(multigraph_gauge, scalar)),
            "multiple-edge/loop determinant gauge",
        )
        determinant_gates += 1

    try:
        gauge_matrix(multigraph, (Fraction(1), Fraction(0), Fraction(2)))
    except ValueError:
        pass
    else:
        raise RuntimeError("zero-weight hostile was not rejected")
    return determinant_gates, power_gates


def tournament(n: int, code: int) -> list[list[int]]:
    adjacency = [[0 for _ in range(n)] for _ in range(n)]
    bit = 0
    for i, j in combinations(range(n), 2):
        if (code >> bit) & 1:
            adjacency[i][j] = 1
        else:
            adjacency[j][i] = 1
        bit += 1
    return adjacency


def directed_triangles(adjacency: list[list[int]]) -> int:
    answer = 0
    for i, j, k in combinations(range(len(adjacency)), 3):
        answer += (
            adjacency[i][j] * adjacency[j][k] * adjacency[k][i]
            + adjacency[i][k] * adjacency[k][j] * adjacency[j][i]
        )
    return answer


def permutation_sign(permutation: tuple[int, ...]) -> int:
    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(len(permutation))
        for j in range(i + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def polynomial_multiply(left: list[int], right: list[int]) -> list[int]:
    answer = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            answer[i + j] += a * b
    return answer


def tournament_determinant_polynomial(adjacency: list[list[int]]) -> list[int]:
    """Leibniz expansion of det(I-uA), independent of Bareiss evaluation."""

    n = len(adjacency)
    answer = [0] * (n + 1)
    for permutation in permutations(range(n)):
        term = [1]
        for row, column in enumerate(permutation):
            if row == column:
                entry = [1]
            elif adjacency[row][column]:
                entry = [0, -1]
            else:
                entry = [0]
            term = polynomial_multiply(term, entry)
            if not any(term):
                break
        sign = permutation_sign(permutation)
        for degree, coefficient in enumerate(term):
            answer[degree] += sign * coefficient
    return answer


def tournament_polynomial_at(adjacency: list[list[int]], x: int) -> int:
    n = len(adjacency)
    return determinant_bareiss(
        [
            [int(i == j) - x * adjacency[i][j] for j in range(n)]
            for i in range(n)
        ]
    )


def vp_integer(value: int, prime: int) -> int | None:
    if value == 0:
        return None
    value = abs(value)
    answer = 0
    while value % prime == 0:
        value //= prime
        answer += 1
    return answer


def audit_padic_tangent() -> tuple[int, int, tuple[int, int, int, int, int]]:
    coefficient_gates = 0
    valuation_gates = 0
    strict_hostile: tuple[int, int, int, int, int] | None = None
    for n in range(1, 7):
        edge_count = n * (n - 1) // 2
        for code in range(1 << edge_count):
            adjacency = tournament(n, code)
            c3 = directed_triangles(adjacency)

            if n <= 5:
                polynomial = tournament_determinant_polynomial(adjacency)
                require(polynomial[0] == 1, "constant determinant coefficient")
                coefficient_u3 = polynomial[3] if len(polynomial) > 3 else 0
                require(coefficient_u3 == -c3, "u^3 coefficient is not -c3")
                coefficient_gates += 1

            for prime in (2, 3, 5):
                x = prime
                polynomial_value = tournament_polynomial_at(adjacency, x)
                require(polynomial_value % prime == 1 % prime, "P(x) is not a p-unit")
                valuation = vp_integer(1 - polynomial_value, prime)
                require(
                    valuation is None or valuation >= 3,
                    "unconditional 3m lower bound failed",
                )
                if c3 % prime != 0:
                    require(valuation == 3, "unit-c3 sharp tangent failed")
                elif c3 > 0 and valuation != 3 and strict_hostile is None:
                    strict_hostile = (n, code, c3, prime, -1 if valuation is None else valuation)
                valuation_gates += 1

            if n <= 5:
                # A nontrivial unit factor and m=2 exercise the actual valuation
                # statement rather than only the special inputs x=p^m.
                prime = 7
                m = 2
                x = (prime**m) * (prime + 1)
                polynomial_value = tournament_polynomial_at(adjacency, x)
                valuation = vp_integer(1 - polynomial_value, prime)
                require(valuation is None or valuation >= 3 * m, "m=2 lower bound")
                if c3 % prime != 0:
                    require(valuation == 3 * m, "m=2 unit-c3 sharp tangent")
                valuation_gates += 1

    require(strict_hostile is not None, "no p|c3 strict hostile found")
    return coefficient_gates, valuation_gates, strict_hostile


def quotient_representative(a: int, prime: int = 13) -> int:
    a %= prime
    require(a != 0, "unit class required")
    return min(a, prime - a)


def odd_pm_representative(a: int, prime: int = 13) -> int:
    a %= prime
    require(a != 0, "unit class required")
    return a if a % 2 else prime - a


def legendre(a: int, prime: int) -> int:
    residue = pow(a % prime, (prime - 1) // 2, prime)
    require(residue in (1, prime - 1), "Legendre unit")
    return 1 if residue == 1 else -1


def ordinal_children(r: int, s: int) -> tuple[tuple[str, int, int], ...]:
    return (
        ("L", r + 2 * s - 1, s),
        ("M", 2 * r + s - 1, r),
        ("R", 2 * r - s, r),
    )


def shell_words(max_r: int) -> dict[tuple[int, int], str]:
    queue = [(2, 1, "")]
    words: dict[tuple[int, int], str] = {}
    while queue:
        r, s, word = queue.pop(0)
        require((r, s) not in words, "Berggren ordinal duplicate")
        words[(r, s)] = word
        for letter, child_r, child_s in ordinal_children(r, s):
            if child_r <= max_r:
                queue.append((child_r, child_s, word + letter))
    return words


def audit_c6_shell() -> tuple[list[int], list[tuple[int, tuple[int, int, int], str, int]]]:
    prime = 13
    classes = tuple(range(1, 7))
    require(
        {quotient_representative(a, prime) for a in range(1, prime)} == set(classes),
        "U13/+-1 classes",
    )
    require(
        {odd_pm_representative(a, prime) for a in classes}
        == {1, 3, 5, 7, 9, 11},
        "odd representative section",
    )

    generator_cycle: list[int] = []
    residue = 1
    quotient_powers: list[int] = []
    for _ in range(6):
        quotient_powers.append(quotient_representative(residue, prime))
        generator_cycle.append(odd_pm_representative(residue, prime))
        residue = (2 * residue) % prime
    require(len(set(quotient_powers)) == 6, "class of two does not generate C6")
    require(residue in (1, prime - 1), "class of two order is not six")
    require(generator_cycle == [1, 11, 9, 5, 3, 7], "generator cycle")

    for a in classes:
        require(legendre(a, prime) == legendre(-a, prime), "character descent")
        for b in classes:
            product_class = quotient_representative(a * b, prime)
            require(
                legendre(product_class, prime) == legendre(a, prime) * legendre(b, prime),
                "quotient character law",
            )

    words = shell_words(7)
    shell = {(r, s): word for (r, s), word in words.items() if r == 7}
    require(len(shell) == 6, "r=7 shell size")
    expected_words = {
        1: "LLLLL",
        3: "ML",
        5: "RM",
        7: "LLR",
        9: "LRR",
        11: "RRRRR",
    }
    rows: list[tuple[int, tuple[int, int, int], str, int]] = []
    for d in (1, 3, 5, 7, 9, 11):
        s = (d + 1) // 2
        word = shell[(7, s)]
        require(word == expected_words[d], "Berggren word table")
        triple = (13 * d, (169 - d * d) // 2, (169 + d * d) // 2)
        a, b, c = triple
        require(a * a + b * b == c * c, "Pythagorean identity")
        require(b % 2 == 0 and b + c == 169 and c - b == d * d, "shell equations")
        require(gcd(a, gcd(b, c)) == 1, "primitive shell triple")
        rows.append((d, triple, word, legendre(d, prime)))

    cycle_signs = [legendre(d, prime) for d in generator_cycle]
    require(cycle_signs == [1, -1, 1, -1, 1, -1], "unique order-two C6 row")
    involution = quotient_representative(5, prime)
    require(
        involution != 1 and quotient_representative(involution * involution, prime) == 1,
        "nonidentity order-two class",
    )
    return generator_cycle, rows


def main() -> None:
    determinant_gates, power_gates = audit_gauge()
    coefficient_gates, valuation_gates, hostile = audit_padic_tangent()
    cycle, rows = audit_c6_shell()

    print("INDEPENDENT STANDARD-LIBRARY EXACT AUDIT")
    print(
        "gauge: all 512 looped/digoned 3-vertex digraphs x 3 nonzero Q-weight rows; "
        f"determinant_gates={determinant_gates} power_trace_gates={power_gates}"
    )
    print(
        "p-adic tangent: all labelled tournaments n<=6; "
        f"coefficient_gates={coefficient_gates} valuation_gates={valuation_gates}"
    )
    print(
        "necessity hostile p|c3: "
        f"n={hostile[0]} code={hostile[1]} c3={hostile[2]} p={hostile[3]} "
        f"v_p(zeta(p)-1)={hostile[4]} (not 3)"
    )
    print("U13/{+-1} generator-two odd cycle:", cycle)
    for d, triple, word, sign in rows:
        print(f"d={d:2d} triple={triple!s:16s} word={word:5s} chi2={sign:+d}")
    print("all independent exact gates PASS")


if __name__ == "__main__":
    main()
