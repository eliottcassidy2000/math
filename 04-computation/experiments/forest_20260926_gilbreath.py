"""Exact probes for forest_20260926_gilbreath.md; standard library only.

Run from repository root with python -X utf8, also with -O.  All validation
uses explicit exceptions, so optimized execution retains every check.
"""

from itertools import combinations, permutations, product
from math import comb, prod


def check(condition, label):
    if not condition:
        raise RuntimeError(label)


def triangle(row):
    rows = [list(row)]
    while len(rows[-1]) > 1:
        old = rows[-1]
        rows.append([abs(b - a) for a, b in zip(old, old[1:])])
    return rows


def edge(row):
    return [r[0] for r in triangle(row)]


def pascal_bits(bits):
    return [sum(comb(i, j) * bits[j] for j in range(i + 1)) % 2
            for i in range(len(bits))]


def newton_seed(values):
    return [sum(comb(i, j) * values[j] for j in range(i + 1))
            for i in range(len(values))]


def lift_matrix(signs, n):
    rows = [[int(i == j) for j in range(n)] for i in range(n)]
    result = [rows[0]]
    for layer in signs:
        rows = [[s * (b - a) for a, b in zip(rows[i], rows[i + 1])]
                for i, s in enumerate(layer)]
        result.append(rows[0])
    return result


def inverse_lower(matrix, values):
    result = []
    for i, value in enumerate(values):
        numerator = value - sum(matrix[i][j] * result[j] for j in range(i))
        check(abs(matrix[i][i]) == 1, "unimodular diagonal")
        result.append(numerator // matrix[i][i])
    return result


def signs_of(rows):
    return [[1 if b >= a else -1 for a, b in zip(row, row[1:])]
            for row in rows[:-1]]


def full_order(rows):
    return tuple(tuple((row[j] > row[i]) - (row[j] < row[i])
                       for i, j in combinations(range(len(row)), 2))
                 for row in rows)


def shortcut(n):
    return (3 * n + 1) // 2 if n % 2 else n // 2


def orbit_prefix(n, length):
    values = []
    for _ in range(length):
        values.append(n)
        n = shortcut(n)
    return values


def hamilton_count(arcs):
    n = len(arcs)
    return sum(all(arcs[p[i]][p[i + 1]] for i in range(n - 1))
               for p in permutations(range(n)))


def sieve(limit):
    flags = bytearray(b"\x01") * (limit + 1)
    flags[:2] = b"\x00\x00"
    for p in range(2, int(limit ** 0.5) + 1):
        if flags[p]:
            flags[p * p::p] = b"\x00" * (((limit - p * p) // p) + 1)
    return [i for i in range(limit + 1) if flags[i]]


def main():
    binary_count = 0
    for n in range(1, 13):
        for bits in product(range(2), repeat=n):
            seed = pascal_bits(bits)
            check(pascal_bits(seed) == list(bits), "Pascal involution")
            check(edge(seed) == list(bits), "binary CA edge universality")
            binary_count += 1
    print(f"Pascal involution and binary edge: {binary_count} traces, lengths 1..12")

    newton_count = 0
    for n in range(1, 8):
        for values in product(range(3), repeat=n):
            seed = newton_seed(values)
            rows = triangle(seed)
            check([row[0] for row in rows] == list(values), "integer edge universality")
            for r, row in enumerate(rows):
                expected = [sum(comb(i, j) * values[r + j] for j in range(i + 1))
                            for i in range(n - r)]
                check(row == expected, "Newton shifted-row formula")
                check(all(b >= a for a, b in zip(row, row[1:])), "monotone rows")
            newton_count += 1
    print(f"Nonnegative integer edge universality: {newton_count} traces over 0..2, lengths 1..7")
    check(newton_seed([2] + [1] * 7) == [1 + 2 ** i for i in range(8)], "unit-edge seed")
    print("Unit-edge nonprime seed:", newton_seed([2] + [1] * 7))

    dyadic_count = 0
    for m in range(5):
        jump = 2 ** m
        for mask in range(256):
            bits = [(mask >> (i % 8)) & 1 for i in range(2 * jump + 8)]
            work = bits[:]
            for _ in range(jump):
                work = [a ^ b for a, b in zip(work, work[1:])]
            check(work == [bits[i] ^ bits[i + jump] for i in range(len(work))], "dyadic recursion")
            dyadic_count += 1
    print(f"Dyadic D^(2^m)=1+S^(2^m): {dyadic_count} exact controls")

    finite_count = 0
    for n in range(1, 9):
        period = 1 << (n - 1).bit_length()
        for seed in product(range(2), repeat=n):
            trace = [sum(comb(r, j) * seed[j] for j in range(min(n, r + 1))) % 2
                     for r in range(3 * period)]
            check(trace[:2 * period] == trace[period:], "finite support periodicity")
            finite_count += 1
    print(f"Finite-support binary seeds: {finite_count} controls, power-of-two global periods")

    actual_count = 0
    for n in range(1, 7):
        for initial in product(range(5), repeat=n):
            rows = triangle(initial)
            signs = signs_of(rows)
            matrix = lift_matrix(signs, n)
            trace = [row[0] for row in rows]
            check(inverse_lower(matrix, trace) == list(initial), "signed edge inversion")
            check(prod(matrix[i][i] for i in range(n)) == prod(s for layer in signs for s in layer),
                  "signed determinant formula")
            check(all(matrix[i][j] % 2 == (comb(i, j) % 2 if j <= i else 0)
                      for i in range(n) for j in range(n)), "signed matrix Pascal reduction")
            actual_count += 1
    schedule_count = 0
    for n in range(1, 6):
        for flat in product((-1, 1), repeat=n * (n - 1) // 2):
            cursor = 0
            signs = []
            for size in range(n - 1, 0, -1):
                signs.append(flat[cursor:cursor + size])
                cursor += size
            matrix = lift_matrix(signs, n)
            check(prod(matrix[i][i] for i in range(n)) == prod(flat), "arbitrary sign determinant")
            check(all(matrix[i][j] % 2 == (comb(i, j) % 2 if j <= i else 0)
                      for i in range(n) for j in range(n)), "arbitrary signs modulo two")
            schedule_count += 1
    print(f"Signed unimodular lift: {actual_count} actual rows over 0..4, lengths 1..6")
    print(f"Arbitrary sign schedules: {schedule_count} at lengths 1..5; feasibility not inferred")
    check(edge([2, 5, 7]) == edge([2, 5, 9]) == [2, 3, 1], "unsigned edge noninjectivity")
    print("Unsigned edge collision: [2,5,7] and [2,5,9] both map to [2,3,1]")

    for modulus in range(3, 33):
        check(1 % modulus != abs(1 - modulus) % modulus, "residue quotient hostile")
    print("Absolute-difference quotient modulo m: canonical obstruction for every m=3..32")

    parity_count = 0
    for n in range(2, 8):
        for tail in product((1, 3, 5), repeat=n - 1):
            rows = triangle((2,) + tail)
            for row in rows[1:]:
                check(row[0] % 2 == 1 and all(x % 2 == 0 for x in row[1:]), "Gilbreath parity")
            parity_count += 1
    primes = sieve(8000)[:1000]
    check(len(primes) == 1000, "prime test size")
    check(all(x == 1 for x in edge(primes)[1:]), "finite Gilbreath control")
    check(edge([2, 3, 7]) == [2, 1, 3], "nonconsecutive prime hostile")
    print(f"Gilbreath parity: {parity_count} synthetic rows; unit edge checked for first 1000 primes only")
    print("Prime but nonconsecutive hostile: [2,3,7] has edge [2,1,3]")

    reference = full_order(triangle(orbit_prefix(23, 5)))
    for s in range(1001):
        values = orbit_prefix(23 + 32 * s, 5)
        expected = [23 + 32 * s, 35 + 48 * s, 53 + 72 * s, 80 + 108 * s, 40 + 54 * s]
        check(values == expected, "Collatz affine family")
        rows = triangle(values)
        check(full_order(rows) == reference, "complete orientations fixed")
        check([x % 2 for x in values] == [1, 1, 1, 0, 0], "parity cylinder")
        check(rows[-1] == [1 + 2 * s], "unbounded magnitude despite orientations")
    print("Actual Collatz family: 1001 starts n=23+32s, same 11100 parity and every row orientation")
    print("Final absolute difference is 1+2s; n=23 gives 1, n=55 gives 3")

    sources = [1, 2, 3, 7, 23, 27, 55, 97, 871, 6171, 63728127]
    for source in sources:
        values = orbit_prefix(source, 128)
        bits = [n % 2 for n in values]
        check(edge(pascal_bits(bits)) == bits, "causal Collatz parity code")
        check(edge(newton_seed(values)) == values, "causal Collatz integer code")
    print(f"Causal source-derived codes: {len(sources)} actual Collatz sources, 128 values each")

    tournament_count = 0
    for n in range(1, 6):
        pairs = list(combinations(range(n), 2))
        for bits in product(range(2), repeat=len(pairs)):
            arcs = [[False] * n for _ in range(n)]
            for (i, j), bit in zip(pairs, bits):
                arcs[i][j] = bool(bit)
                arcs[j][i] = not bit
            check(hamilton_count(arcs) % 2 == 1, "Redei finite control")
            tournament_count += 1
    scalar_count = 0
    for n in range(1, 7):
        for values in permutations(range(n)):
            arcs = [[a < b for b in values] for a in values]
            check(hamilton_count(arcs) == 1, "scalar tournament unique Hamilton path")
            for i, j, k in product(range(n), repeat=3):
                check(values[k] - values[i] == (values[j] - values[i]) + (values[k] - values[j]),
                      "metric cocycle")
            scalar_count += 1
    print(f"Redei parity: all {tournament_count} labelled tournaments with 1..5 vertices")
    print(f"Scalar orientations and metric cocycles: {scalar_count} distinct-value rows with 1..6 vertices")
    print("PASS: all checks use explicit exceptions and remain active under python -O")


if __name__ == "__main__":
    main()
