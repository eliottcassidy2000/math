"""Independent referee: unique degree-two boards, direct grid cells, Q algebra.

No producer imports or producer code reads. Normal and -O retain every gate.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import comb
import hashlib
import json
import sys

sys.stdout.reconfigure(newline="\n")
gates = 0


def need(test, label):
    global gates
    gates += 1
    if not test:
        raise RuntimeError(label)


def boards(n):
    """Every labeled simple 2-regular board once, by its unordered row pairs."""
    degree = [0] * n
    rows = []

    def visit(r):
        if r == n:
            yield tuple(rows)
            return
        remaining = n - r - 1
        for a, b in combinations([c for c in range(n) if degree[c] < 2], 2):
            degree[a] += 1
            degree[b] += 1
            if all(2 - value <= remaining for value in degree):
                rows.append((a, b))
                yield from visit(r + 1)
                rows.pop()
            degree[a] -= 1
            degree[b] -= 1

    yield from visit(0)


def cycle_type(rows):
    n = len(rows)
    adjacency = [set() for _ in range(2 * n)]
    for r, cols in enumerate(rows):
        for c in cols:
            adjacency[r].add(n + c)
            adjacency[n + c].add(r)
    remaining, result = set(range(2 * n)), []
    while remaining:
        stack, seen = [min(remaining)], set()
        while stack:
            vertex = stack.pop()
            if vertex not in seen:
                seen.add(vertex)
                stack.extend(adjacency[vertex] - seen)
        remaining -= seen
        result.append(len(seen) // 2)
    return tuple(sorted(result))


def grid_triples(points):
    return sum((b[0] - a[0]) * (c[1] - a[1]) == (c[0] - a[0]) * (b[1] - a[1])
               for a, b, c in combinations(points, 3))


def literal_census(n):
    bank = defaultdict(lambda: [0] * 13)
    omitted_witness = None
    for row_pairs in boards(n):
        need(all(sum(c in pair for pair in row_pairs) == 2 for c in range(n)),
             "exact column degree, independent board carrier")
        points = [(r, c) for r, cols in enumerate(row_pairs) for c in cols]
        plus, minus = [0] * (2 * n - 1), [0] * (2 * n - 1)
        for r, c in points:
            plus[c - r + n - 1] += 1
            minus[r + c] += 1
        A = sum(comb(z, 3) for z in plus if z >= 3)
        B = sum(comb(z, 3) for z in minus if z >= 3)
        row = bank[cycle_type(row_pairs)]
        values = (1, A, B, A * A, B * B, A * B, A == 0, A + B == 0,
                  plus[n - 1] * minus[n - 1], plus[n - 1] * minus[n - 2])
        for i, value in enumerate(values):
            row[i] += value
        if n <= 5:
            X = grid_triples(points)
            need(A + B <= X, "two-direction triples are actual distinct collinear events")
            row[10] += X == 0
            if n <= 4:
                need(A + B == X, "all nonaxis grid triples have slopes +-1 at n<=4")
            if A + B == 0 and X > 0 and omitted_witness is None:
                omitted_witness = (points, X)
        row[11] += comb(plus[n - 1], 3) * comb(minus[n - 1], 3)
        row[12] += comb(plus[n - 1], 3) * comb(minus[n - 2], 3)
    need(sum(row[0] for row in bank.values()) == {3: 6, 4: 90, 5: 2040, 6: 67950}[n],
         "complete unique-board count")
    output = []
    for shape, row in sorted(bank.items()):
        N, a, b, aa, bb, ab = row[:6]
        mean = F(2 * n - 5, 3)
        need(F(a, N) == F(b, N) == mean, "both exact direction means")
        varA, varB, covariance = F(aa, N) - mean ** 2, F(bb, N) - mean ** 2, F(ab, N) - mean ** 2
        need(varA == varB, "column reflection preserves the cycle type")
        p1 = F(2, n)
        p2 = F(4 * n - 6, n * (n - 1) ** 2)
        q2 = F(2, n * (n - 1))
        for column, M, Z in ((8, n, n % 2), (9, n - 1, 1 - n % 2)):
            L, R, C = n, M, M
            want = Z * p1 + (R + C - 2 * Z) * q2 + (L * M - R - C + Z) * p2
            need(F(row[column], N) == want, "exact common-cell/row/column two-edge formula")
        output.append({"n": n, "shape": shape, "boards": N, "mean_each": str(mean),
                       "variance_each": str(varA), "cross_covariance": str(covariance),
                       "variance_total": str(2 * varA + 2 * covariance),
                       "P_both_zero": str(F(row[7], N)),
                       "P_full_zero": str(F(row[10], N)) if n <= 5 else None,
                       "selected_triple_pair_moments": (str(F(row[11], N)), str(F(row[12], N)))})
    return output, omitted_witness


def multiply(p, q):
    result = defaultdict(F)
    for (a, b), v in p.items():
        for (c, d), w in q.items():
            result[(a + c, b + d)] += v * w
    return dict(result)


def triangle_integral(p):
    """0<=x<=1/2, x<=y<=1-x; exact monomial antiderivatives."""
    answer = F(0)
    for (i, j), coefficient in p.items():
        upper = sum(F((-1) ** k * comb(j + 1, k), (i + k + 1) * 2 ** (i + k + 1))
                    for k in range(j + 2))
        lower = F(1, (i + j + 2) * 2 ** (i + j + 2))
        answer += coefficient * (upper - lower) / (j + 1)
    return answer


def lagrange(points):
    result = [F(0)] * len(points)
    for i, (x, y) in enumerate(points):
        basis, denom = [F(1)], F(1)
        for j, (xx, _) in enumerate(points):
            if i == j:
                continue
            updated = [F(0)] * (len(basis) + 1)
            for k, v in enumerate(basis):
                updated[k] -= xx * v
                updated[k + 1] += v
            basis = updated
            denom *= x - xx
        for k, v in enumerate(basis):
            result[k] += y * v / denom
    return result


def geometry(n):
    lengths = {d: n - abs(d) for d in range(1 - n, n)}
    anti = {e: n - abs(e - (n - 1)) for e in range(2 * n - 1)}
    grid = sum((n - abs(r - c)) ** 2 * (n - abs(r + c - (n - 1))) ** 2
               for r in range(n) for c in range(n))
    by_lines = 0
    overlap = 0
    for d, L in lengths.items():
        for e, M in anti.items():
            Z = (e - d) % 2 == 0 and 0 <= e - d <= 2 * (n - 1) and 0 <= e + d <= 2 * (n - 1)
            by_lines += L * L * M * M * Z
            overlap += L * L * M * M * (min(L, M) + max(L + M - n, 0))
    need(grid == by_lines, "every common cell counted once; crossing parity exact")
    need(grid == F(19 * n ** 6 + 25 * n ** 4 - 44 * n ** 2, 90) + (n % 2) * n ** 2,
         "exact crossing polynomial, parity retained")
    need(overlap == F(46 * n ** 7 + 85 * n ** 5 + 49 * n ** 3, 90),
         "exact row-plus-column overlap polynomial")
    cube_sum = sum(L ** 3 for L in lengths.values())
    need(cube_sum == F(n ** 2 * (n ** 2 + 1), 2), "exact cubic length sum")
    kernel = 8 * (F(cube_sum ** 2, n ** 8) - F(overlap, n ** 7) + F(grid, n ** 6))
    need(kernel == -F(2, 5) - F(4, 3 * n ** 2) + F(8 * (n % 2), n ** 4) - F(94, 15 * n ** 4),
         "exact finite leading covariance kernel")
    return grid, overlap


def main():
    # Direct unit-square integration, no checkerboard-density assumption.
    u = {(0, 0): F(1), (1, 0): F(1), (0, 1): F(-1)}
    v = {(1, 0): F(1), (0, 1): F(1)}
    common = 4 * triangle_integral(multiply(multiply(u, u), multiply(v, v)))
    need(common == F(19, 90), "direct grid-cell integral")
    same = F(1, 14)
    # Integral x^2 y^2 (x+y-1)+: integrate y from 1-x to1 by binomial expansion.
    # Its upper and lower antiderivatives are integrated in x here literally.
    poly = {(3, 2): F(1), (2, 3): F(1), (2, 2): F(-1)}
    opposite = F(0)
    for (i, j), coefficient in poly.items():
        lower = sum(F((-1) ** k * comb(j + 1, k), i + k + 1) for k in range(j + 2))
        opposite += coefficient * (F(1, i + 1) - lower) / (j + 1)
    need(opposite == F(71, 1260), "independent opposite-interval integral")
    row_col = 4 * (same + opposite)
    need(row_col == F(23, 45), "shared row/column limit")
    cross = 8 * (F(1, 4) - row_col + common)
    total = 2 * F(40, 9) + 2 * cross
    need(cross == -F(2, 5), "cross-direction covariance coefficient")
    need(total == F(364, 45), "two-direction variance coefficient")
    need(total / F(4, 3) ** 2 == F(91, 20), "one-sided Chebyshev coefficient")

    geometries = {n: geometry(n) for n in range(1, 51)}
    # These interpolation identities are exact because chamber/Faulhaber degree
    # bounds give degree6 period2 for crossings and degree7 for row/column overlap.
    for parity in (0, 1):
        ns = list(range(2 if parity == 0 else 1, 15 if parity == 0 else 14, 2))
        coeff = lagrange([(n, geometries[n][0]) for n in ns])
        want = [F(0)] * 7
        want[2], want[4], want[6] = F(-44, 90) + parity, F(25, 90), F(19, 90)
        need(coeff == want, "independent bounded-degree parity reconstruction")
    coeff = lagrange([(n, geometries[n][1]) for n in range(1, 9)])
    want = [F(0)] * 8
    want[3], want[5], want[7] = F(49, 90), F(85, 90), F(46, 90)
    need(coeff == want, "independent bounded-degree overlap reconstruction")

    census = []
    witnesses = []
    for n in range(3, 7):
        rows, witness = literal_census(n)
        census.extend(rows)
        if witness:
            witnesses.append((n, witness))
    bank = {(row["n"], tuple(row["shape"])): row for row in census}
    need(bank[(4, (4,))]["P_both_zero"] == "1/36", "frozen C8 full-zero hostile")
    need(bank[(4, (2, 2))]["P_both_zero"] == "1/2", "frozen 2C4 full-zero hostile")
    need(bank[(5, (5,))]["P_full_zero"] == "1/45", "n5 C10 literal determinant control")
    need(bank[(5, (2, 3))]["P_full_zero"] == "0", "n5 C4+C6 literal determinant control")
    need(len(witnesses) > 0, "two directions still do not characterize full safety")
    payload = {"census": census, "omitted_direction_witness": witnesses[0]}
    print("STATUS: PASS; independent proof, unique-board census, and exact crossing geometry")
    print("BOARD UNIVERSE: n3..6, 70086 unique simple degree-two boards, all cycle types")
    for row in census:
        print(json.dumps(row, sort_keys=True))
    print("OTHER-DIRECTION HOSTILE:", witnesses[0])
    print("CROSSING NUMERATOR: (19n^6+25n^4-44n^2)/90 + (n mod2)n^2")
    print("ROW/COLUMN NUMERATOR: (46n^7+85n^5+49n^3)/90")
    print("LIMITS: common=19/90; row+col=23/45; cross=-2/5; total variance=364/45; probability coefficient=91/20")
    print("SCOPE: covariance theorem and uniform random-density bound, not nonexistence or optimality")
    print("SEMANTIC SHA256:", hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest())
    print("ACTIVE GATES:", gates)


if __name__ == "__main__":
    main()
