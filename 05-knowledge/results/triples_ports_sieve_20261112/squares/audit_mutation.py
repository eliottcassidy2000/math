"""Independent small symbolic-case audit of the sibling mutation proof.

This checks its claimed complete rank-three return rows and rank-four
obstruction table. The all-size step remains the proof by restriction.
"""

from itertools import combinations, product
from pathlib import Path
import json


def mutate(matrix, k):
    n = len(matrix)
    result = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == k or j == k:
                result[i][j] = -matrix[i][j]
            else:
                x, y = matrix[i][k], matrix[k][j]
                result[i][j] = (matrix[i][j] + max(x, 0) * max(y, 0)
                                - max(-x, 0) * max(-y, 0))
    return result


def matrix(n, entries):
    result = [[0] * n for _ in range(n)]
    for (i, j), value in entries.items():
        result[i][j], result[j][i] = value, -value
    return result


def apply(initial, word):
    for k in word:
        initial = mutate(initial, k)
    return initial


def tournament(matrix):
    return all(abs(matrix[i][j]) == 1 for i in range(len(matrix))
               for j in range(i + 1, len(matrix)))


def require(value, message):
    if not value:
        raise RuntimeError(message)


def four_step_audit():
    words = ((0, 1, 0, 1), (0, 1, 0, 2), (0, 1, 2, 0),
             (0, 1, 2, 1), (0, 1, 2, 3))
    counts = []
    for n in range(3, 6):
        pairs = list(combinations(range(n), 2))
        returns = {str(word): 0 for word in words if max(word) < n}
        for signs in product((-1, 1), repeat=len(pairs)):
            initial = matrix(n, dict(zip(pairs, signs)))
            swapped = [row[:] for row in initial]
            swap = [1, 0] + list(range(2, n))
            swapped = [[initial[swap[i]][swap[j]] for j in range(n)] for i in range(n)]
            require(apply(initial, (0, 1, 0, 1, 0)) == swapped,
                    "full-matrix rank-two pentagon fails")
            for word in words:
                if max(word) >= n:
                    continue
                stages = [initial]
                for k in word:
                    stages.append(mutate(stages[-1], k))
                final = stages[-1]
                if word == (0, 1, 2, 0) and n >= 4:
                    declared = (initial[0][1] == initial[1][2] == initial[2][0]
                                and len({initial[i][j] for i in range(3) for j in range(3, n)}) == 1)
                    require(tournament(final) == declared, "0120 exact return criterion fails")
                if word == (0, 1, 2, 3):
                    sign = initial[0][1]
                    declared = (all(initial[i][j] == sign for i in range(4) for j in range(i + 1, 4))
                                and all(initial[i][j] == sign for i in range(4) for j in range(4, n)))
                    require(tournament(final) == declared, "0123 exact return criterion fails")
                if not tournament(final):
                    continue
                returns[str(word)] += 1
                if word == (0, 1, 0, 1):
                    require(len({initial[1][j] for j in range(n) if j != 1}) == 1,
                            "0101 did not force extreme vertex1")
                elif word == (0, 1, 0, 2):
                    require(tournament(stages[3]), "0102 intermediate implication fails")
                elif word == (0, 1, 2, 1):
                    require(tournament(stages[1]), "0121 first-native implication fails")
                elif word == (0, 1, 2, 3):
                    sign = initial[0][1]
                    require(all(initial[i][j] == sign for i in range(4) for j in range(i + 1, 4)),
                            "0123 named block not prescribed TT4")
                    require(all(initial[i][j] == sign for i in range(4) for j in range(4, n)),
                            "0123 exterior not on prescribed side")
                elif n >= 4:
                    require(initial[0][1] == initial[1][2] == initial[2][0],
                            "0120 named block not C3")
                    cross = {initial[i][j] for i in range(3) for j in range(3, n)}
                    require(len(cross) == 1, "0120 exterior is not on one uniform side")
                    require(all(final[i][j] == (-initial[i][j] if i < 3 or j < 3 else initial[i][j])
                                for i in range(n) for j in range(i + 1, n)),
                            "0120 output is not declared triangle transfer")
        counts.append({"n": n, "tournaments": 1 << len(pairs), "return_counts": returns})
    return counts


def main():
    rows = {}
    for word in ((0, 1, 0), (0, 1, 2)):
        valid = []
        for a, b, c in product((-1, 1), repeat=3):
            initial = matrix(3, {(0, 1): a, (0, 2): b, (1, 2): c})
            if tournament(apply(initial, word)):
                valid.append((a, b, c))
        rows[str(word)] = valid
    require(rows[str((0, 1, 0))] == [(-1, 1, 1), (1, -1, -1)],
            "rank-three repeated-endpoint return rows differ")
    require(rows[str((0, 1, 2))] == [(-1, -1, -1), (1, 1, 1)],
            "rank-three distinct-name return rows differ")
    obstruction = []
    for e, f in product((-1, 1), repeat=2):
        initial = matrix(4, {(0, 1): 1, (0, 2): 1, (1, 2): 1,
                             (0, 3): -1, (1, 3): e, (2, 3): f})
        final = apply(initial, (0, 1, 2))
        obstruction.append({"e_f": [e, f], "remaining_column": [final[i][3] for i in range(3)]})
        require(not tournament(final), "rank-four hostile returns a tournament")
    require([r["remaining_column"] for r in obstruction]
            == [[-5, -2, 4], [-3, 0, 2], [-1, -2, 2], [1, 0, 0]],
            "rank-four obstruction table differs")
    receipt = {"status": "FINITE-EXACT independent case-table audit",
               "source": "../ports/PROOF.md, propositions4-5",
               "rank_three_complete_return_rows": rows,
               "rank_four_obstruction_rows": obstruction,
               "four_step_local_implications": four_step_audit()}
    destination = Path(__file__).with_name("mutation_audit.json")
    destination.write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    print("PASS: rank-three/ rank-four tables and all five four-step local implications")


if __name__ == "__main__":
    main()
