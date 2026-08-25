#!/usr/bin/env python3
"""Exact order-seven audit of the open THM-1950 strong base.

Primary path: Hamiltonian subset DP plus Bareiss determinant and exact
Fraction Gaussian elimination.  The existing 456 canonical representatives
are a proved complete isomorphism-class universe (THM-1370).
"""

from fractions import Fraction
from hashlib import sha256


N = 7
DEN = 1 << (N - 1)
REPS_PATH = "05-knowledge/results/n7_class_reps_boxeph_S152.txt"
PAIRS = tuple((i, j) for i in range(N) for j in range(i + 1, N))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def out_masks(bits):
    out = [0] * N
    for k, (i, j) in enumerate(PAIRS):
        if (bits >> k) & 1:
            out[i] |= 1 << j
        else:
            out[j] |= 1 << i
    return tuple(out)


def reverse_masks(out):
    full = (1 << N) - 1
    return tuple((full ^ (1 << i)) & ~out[i] for i in range(N))


def reachable(out, start):
    seen = 1 << start
    frontier = seen
    while frontier:
        bit = frontier & -frontier
        frontier ^= bit
        v = bit.bit_length() - 1
        new = out[v] & ~seen
        seen |= new
        frontier |= new
    return seen


def is_strong(out):
    full = (1 << N) - 1
    return reachable(out, 0) == full and reachable(reverse_masks(out), 0) == full


def hamiltonian_paths(out):
    size = 1 << N
    dp = [[0] * N for _ in range(size)]
    for v in range(N):
        dp[1 << v][v] = 1
    for mask in range(size):
        for v in range(N):
            count = dp[mask][v]
            if not count:
                continue
            choices = out[v] & ~mask
            while choices:
                bit = choices & -choices
                choices ^= bit
                w = bit.bit_length() - 1
                dp[mask | bit][w] += count
    return sum(dp[-1])


def matrix_i_plus_k(out):
    matrix = [[0] * N for _ in range(N)]
    for i in range(N):
        matrix[i][i] = 1
        for j in range(i + 1, N):
            sign = 1 if (out[i] >> j) & 1 else -1
            matrix[i][j] = sign
            matrix[j][i] = -sign
    return matrix


def bareiss_det(matrix):
    a = [row[:] for row in matrix]
    sign = 1
    previous = 1
    n = len(a)
    for k in range(n - 1):
        pivot = next((r for r in range(k, n) if a[r][k]), None)
        require(pivot is not None, "singular I+K")
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        p = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * p - a[i][k] * a[k][j]
                require(numerator % previous == 0, "nonexact Bareiss division")
                a[i][j] = numerator // previous
        for i in range(k + 1, n):
            a[i][k] = 0
        previous = p
    return sign * a[-1][-1]


def solve_ones(matrix):
    n = len(matrix)
    a = [[Fraction(x) for x in row] + [Fraction(1)] for row in matrix]
    for col in range(n):
        pivot = next((r for r in range(col, n) if a[r][col]), None)
        require(pivot is not None, "singular solve")
        a[col], a[pivot] = a[pivot], a[col]
        p = a[col][col]
        a[col] = [x / p for x in a[col]]
        for r in range(n):
            if r == col or not a[r][col]:
                continue
            q = a[r][col]
            a[r] = [a[r][j] - q * a[col][j] for j in range(n + 1)]
    return tuple(a[i][-1] for i in range(n))


def invariants(bits):
    out = out_masks(bits)
    h = hamiltonian_paths(out)
    matrix = matrix_i_plus_k(out)
    determinant = bareiss_det(matrix)
    require(determinant > 0, "det(I+K) must be positive")
    even = Fraction(determinant, DEN)
    response = sum(solve_ones(matrix), Fraction(0))
    odd = even * response
    require(even.denominator == 1 and odd.denominator == 1, "nonintegral energy")
    return is_strong(out), h, int(even), int(odd), response


def paley_bits():
    residues = {1, 2, 4}
    bits = 0
    for k, (i, j) in enumerate(PAIRS):
        if (j - i) % 7 in residues:
            bits |= 1 << k
    return bits


def main():
    with open(REPS_PATH, encoding="utf-8") as handle:
        reps = tuple(int(line) for line in handle if line.strip())
    require(len(reps) == 456 and len(set(reps)) == 456, "order-seven universe changed")
    require(all(0 <= bits < (1 << len(PAIRS)) for bits in reps), "invalid representative")

    rows = []
    complement_checks = 0
    for bits in reps:
        strong, h, even, odd, response = invariants(bits)
        target = max(even, odd)
        require(h >= target, f"strong-base inequality failed at {bits}")
        comp = ((1 << len(PAIRS)) - 1) ^ bits
        cstrong, ch, ceven, codd, cresponse = invariants(comp)
        require((strong, h, even, odd, response) ==
                (cstrong, ch, ceven, codd, cresponse), "converse invariance failed")
        complement_checks += 1
        rows.append((bits, strong, h, even, odd, response))

    strong_rows = [row for row in rows if row[1]]
    require(len(strong_rows) == 353, "strong-class count changed")
    equality = tuple(row[0] for row in strong_rows if row[2] == max(row[3], row[4]))
    slack_row = min(strong_rows, key=lambda row: row[2] - max(row[3], row[4]))
    ratio_row = min(
        strong_rows,
        key=lambda row: Fraction(row[2], max(row[3], row[4])),
    )
    min_response = min(strong_rows, key=lambda row: row[5])
    max_response = max(strong_rows, key=lambda row: row[5])

    pstrong, ph, peven, podd, presponse = invariants(paley_bits())
    require((pstrong, ph, peven, podd, presponse) == (True, 189, 8, 56, Fraction(7)),
            "Paley-seven positive control changed")

    serialization = "\n".join(
        f"{bits}:{int(strong)}:{h}:{even}:{odd}:{response}"
        for bits, strong, h, even, odd, response in rows
    ) + "\n"
    digest = sha256(serialization.encode()).hexdigest()

    print("TOURNAMENT STRONG-BASE ORDER-SEVEN EXACT FRONTIER")
    print("status=FINITE-EXACT; all-order H>=disc remains OPEN")
    print(f"universe=THM1370_canonical_order7_reps;classes={len(rows)};strong={len(strong_rows)}")
    print(f"method=Hamiltonian_subset_DP+Bareiss+Fraction_inverse;converse_checks={complement_checks}")
    print(f"strong_equalities={equality}")
    print("minimum_strong_slack="
          f"{slack_row[2]-max(slack_row[3],slack_row[4])};"
          f"rep={slack_row[0]};H={slack_row[2]};Eeven={slack_row[3]};Eodd={slack_row[4]}")
    print("minimum_strong_ratio="
          f"{Fraction(ratio_row[2],max(ratio_row[3],ratio_row[4]))};"
          f"rep={ratio_row[0]};H={ratio_row[2]};Eeven={ratio_row[3]};Eodd={ratio_row[4]}")
    print(f"strong_response_range={min_response[5]}..{max_response[5]};"
          f"min_rep={min_response[0]};max_rep={max_response[0]}")
    print("paley7_control=strong:True,H:189,Eeven:8,Eodd:56,response:7")
    print(f"semantic_digest={digest}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
