#!/usr/bin/env python3
"""Independent order-seven tournament strong-base audit.

No production code is imported.  This route counts Hamiltonian paths by
literal permutations, detects strongness by dominating cuts, and computes
both energies from principal Pfaffian-square sums.
"""

from hashlib import sha256
from itertools import combinations, permutations
from fractions import Fraction


N = 7
DEN = 1 << (N - 1)
REPS_PATH = "05-knowledge/results/n7_class_reps_boxeph_S152.txt"
PAIRS = tuple((i, j) for i in range(N) for j in range(i + 1, N))
ALL_PERMS = tuple(permutations(range(N)))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def skew_matrix(bits):
    k_matrix = [[0] * N for _ in range(N)]
    for index, (i, j) in enumerate(PAIRS):
        sign = 1 if (bits >> index) & 1 else -1
        k_matrix[i][j] = sign
        k_matrix[j][i] = -sign
    return k_matrix


def is_strong_by_cuts(k_matrix):
    full = (1 << N) - 1
    for mask in range(1, full):
        complement = full ^ mask
        dominates = True
        for i in range(N):
            if not (mask >> i) & 1:
                continue
            for j in range(N):
                if (complement >> j) & 1 and k_matrix[i][j] != 1:
                    dominates = False
                    break
            if not dominates:
                break
        if dominates:
            return False
    return True


def hamiltonian_paths_by_permutations(k_matrix):
    return sum(
        all(k_matrix[order[i]][order[i + 1]] == 1 for i in range(N - 1))
        for order in ALL_PERMS
    )


def pfaffian(matrix):
    n = len(matrix)
    if n == 0:
        return 1
    require(n % 2 == 0, "Pfaffian needs even order")
    total = 0
    for j in range(1, n):
        minor = [
            [matrix[r][c] for c in range(n) if c not in (0, j)]
            for r in range(n) if r not in (0, j)
        ]
        total += (1 if j % 2 else -1) * matrix[0][j] * pfaffian(minor)
    return total


def principal(matrix, subset):
    return [[matrix[i][j] for j in subset] for i in subset]


def rooted_matrix(k_matrix, subset):
    size = len(subset)
    result = [[0] * (size + 1) for _ in range(size + 1)]
    for a, i in enumerate(subset):
        for b, j in enumerate(subset):
            result[a][b] = k_matrix[i][j]
        result[a][-1] = 1
        result[-1][a] = -1
    return result


def energies(k_matrix):
    even_numerator = 0
    odd_numerator = 0
    vertices = tuple(range(N))
    for size in range(N + 1):
        for subset in combinations(vertices, size):
            if size % 2 == 0:
                value = pfaffian(principal(k_matrix, subset))
                even_numerator += value * value
            else:
                value = pfaffian(rooted_matrix(k_matrix, subset))
                odd_numerator += value * value
    require(even_numerator % DEN == 0 and odd_numerator % DEN == 0,
            "nonintegral Pfaffian energy")
    return even_numerator // DEN, odd_numerator // DEN


def main():
    with open(REPS_PATH, encoding="utf-8") as handle:
        reps = tuple(int(line) for line in handle if line.strip())
    require(len(reps) == 456 and len(set(reps)) == 456, "order-seven universe changed")

    rows = []
    for bits in reps:
        matrix = skew_matrix(bits)
        strong = is_strong_by_cuts(matrix)
        h = hamiltonian_paths_by_permutations(matrix)
        even, odd = energies(matrix)
        require(h >= max(even, odd), f"inequality failed at {bits}")
        rows.append((bits, strong, h, even, odd))

    strong_rows = [row for row in rows if row[1]]
    require(len(strong_rows) == 353, "strong-class count changed")
    equality = tuple(row[0] for row in strong_rows if row[2] == max(row[3], row[4]))
    slack_row = min(strong_rows, key=lambda row: row[2] - max(row[3], row[4]))
    ratio_row = min(
        strong_rows,
        key=lambda row: (Fraction(row[2], max(row[3], row[4])), row[0]),
    )

    # The primary serialization has one additional response field.  Since
    # response=Eodd/Eeven, reconstruct its reduced rational spelling here.
    serialization = "\n".join(
        f"{bits}:{int(strong)}:{h}:{even}:{odd}:{Fraction(odd, even)}"
        for bits, strong, h, even, odd in rows
    ) + "\n"
    digest = sha256(serialization.encode()).hexdigest()

    print("TOURNAMENT STRONG-BASE ORDER-SEVEN INDEPENDENT AUDIT")
    print("status=FINITE-EXACT; no all-order consequence")
    print(f"universe=456_canonical_reps;strong={len(strong_rows)}")
    print("method=literal_permutations+dominating_cuts+Pfaffian_square_sums")
    print(f"strong_equalities={equality}")
    print("minimum_strong_slack="
          f"{slack_row[2]-max(slack_row[3],slack_row[4])};"
          f"rep={slack_row[0]};H={slack_row[2]};Eeven={slack_row[3]};Eodd={slack_row[4]}")
    print("minimum_strong_ratio="
          f"{Fraction(ratio_row[2],max(ratio_row[3],ratio_row[4]))};"
          f"rep={ratio_row[0]};H={ratio_row[2]};Eeven={ratio_row[3]};Eodd={ratio_row[4]}")
    print(f"semantic_digest={digest}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
