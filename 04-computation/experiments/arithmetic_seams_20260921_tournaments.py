"""Exact replay for pointed triangular identities and fixed-path tournaments."""
from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
import json
from math import comb
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def triangular(n: int) -> int:
    return n * (n + 1) // 2


def frame(n: int, mask: int = 0) -> list[list[bool]]:
    out = [[i < j for j in range(n)] for i in range(n)]
    bit = 0
    for i in range(n):
        for j in range(i + 2, n):
            if mask >> bit & 1:
                out[i][j], out[j][i] = False, True
            bit += 1
    return out


def families(n: int) -> tuple[list[list[bool]], list[list[bool]]]:
    intervals = frame(n, (1 << comb(n - 1, 2)) - 1)
    endpoint = frame(n)
    endpoint[n - 1][0], endpoint[0][n - 1] = True, False
    return intervals, endpoint


def triangle_count(out: list[list[bool]]) -> int:
    return sum(out[i][j] == out[j][k] == out[k][i]
               for i, j, k in combinations(range(len(out)), 3))


def chord_polynomial(out: list[list[bool]]) -> int:
    total = 0
    for i, j, k in combinations(range(len(out)), 3):
        a, b, c = int(out[j][i]), int(out[k][j]), int(out[k][i])
        total += c - c * a - c * b + a * b
    return total


def hamilton_dp(out: list[list[bool]]) -> int:
    n = len(out)
    if n == 0:
        return 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v, count in enumerate(dp[mask]):
            if count:
                for w in range(n):
                    if not (mask >> w & 1) and out[v][w]:
                        dp[mask | 1 << w][w] += count
    return sum(dp[-1])


def cycles(out: list[list[bool]]) -> list[tuple[int, ...]]:
    """Simple directed cycles, with their minimum vertex first: no rotations."""
    n = len(out)
    result = []
    for start in range(n):
        def visit(path: tuple[int, ...], used: int) -> None:
            if len(path) >= 3 and out[path[-1]][start]:
                result.append(path)
            for target in range(start + 1, n):
                if not (used >> target & 1) and out[path[-1]][target]:
                    visit(path + (target,), used | 1 << target)
        visit((start,), 1 << start)
    return result


def odd_collection_weight(n: int, all_cycles: list[tuple[int, ...]]) -> int:
    masks = [sum(1 << v for v in cycle) for cycle in all_cycles if len(cycle) % 2]
    containing = [[mask for mask in masks if mask >> v & 1] for v in range(n)]
    @lru_cache(None)
    def evaluate(available: int) -> int:
        if not available:
            return 1
        bit = available & -available
        v = bit.bit_length() - 1
        return evaluate(available ^ bit) + 2 * sum(
            evaluate(available ^ mask) for mask in containing[v]
            if mask & available == mask)
    return evaluate((1 << n) - 1)


def identity_controls() -> dict:
    checks = 0
    for a in range(21):
        for b in range(21):
            n = (a + 1) * (b + 1)
            proposed = (triangular(a) + triangular(b) + triangular(a * b)
                        + a * (b * b - 1) + b * (a * a - 1))
            corrected = (triangular(a) + triangular(b) + triangular(a * b)
                         + a * b * (a + b + 1))
            require(proposed == triangular(n - 2), "root-deleted polynomial failed")
            require(corrected == triangular(n - 1), "root-retained polynomial failed")
            require(corrected - proposed == n - 1, "root star difference failed")
            require(triangular(a + b) == triangular(a) + triangular(b) + a * b,
                    "additive identity failed")
            require(triangular(a + b + 1) == triangular(a) + triangular(b) + n,
                    "shifted additive identity failed")
            require(comb(n, 2) == (a + 1) * triangular(b) + (b + 1)**2 * triangular(a),
                    "product fibre identity failed")
            checks += 1
    return {"A_B_universe": [0, 20], "pairs_checked": checks,
            "smallest_positive_hostile": {"A": 1, "B": 1, "claimed_left": 6,
                                          "proposed_right": 3}}


def cube_controls() -> dict:
    rows = []
    for n in range(3, 8):
        free = comb(n - 1, 2)
        histogram = Counter()
        for mask in range(1 << free):
            out = frame(n, mask)
            count = triangle_count(out)
            require(count == chord_polynomial(out), "chord interaction identity failed")
            require(count == comb(n, 3) - sum(comb(sum(row), 2) for row in out),
                    "independent score identity failed")
            histogram[count] += 1
            if n <= 5:
                direct = sum(all(out[p[i]][p[i + 1]] for i in range(n - 1))
                             for p in permutations(range(n)))
                require(direct == hamilton_dp(out), "Hamilton DP versus permutation failed")
        rows.append({"n": n, "free_chords": free, "completions": 1 << free,
                     "triangle_histogram": dict(sorted(histogram.items()))})
    return {"fixed_path_cube_universe": "n=3..7, every free-chord bit pattern",
            "total_completions": sum(row["completions"] for row in rows),
            "permutation_HP_cross_check": "all fixed-path cubes n=3..5", "rows": rows}


def family_controls() -> list[dict]:
    tribonacci = [1, 1, 1]
    for n in range(3, 11):
        tribonacci.append(sum(tribonacci[-3:]))
    result = []
    for n in range(3, 11):
        left, right = families(n)
        require(sorted(map(sum, left)) == sorted(map(sum, right)), "score hostile failed")
        record = {"n": n, "shared_scores": sorted(map(sum, left))}
        for name, out, expected_h in [("all_backward_chords", left, tribonacci[n]),
                                      ("single_backward_endpoint", right, 1 + 2**(n - 2))]:
            all_cycles = cycles(out)
            counts = Counter(map(len, all_cycles))
            expected = {ell: n - ell + 1 if name == "all_backward_chords"
                        else comb(n - 2, ell - 2) for ell in range(3, n + 1)}
            require(dict(counts) == expected, "complete simple-cycle law failed")
            require(triangle_count(out) == n - 2, "family triangle law failed")
            observed_h = hamilton_dp(out)
            require(observed_h == expected_h, "family Hamilton formula failed")
            require(observed_h == odd_collection_weight(n, all_cycles),
                    "independent odd-cycle collection count failed")
            record[name] = {"cycle_counts": dict(sorted(counts.items())),
                            "total_cycles": len(all_cycles), "hamilton_paths": observed_h}
        result.append(record)
    return result


def substitution(outer: list[list[bool]], inner: list[list[bool]]) -> list[list[bool]]:
    a, b = len(outer), len(inner)
    return [[inner[i % b][j % b] if i // b == j // b else outer[i // b][j // b]
             for j in range(a * b)] for i in range(a * b)]


def substitution_controls() -> dict:
    c3 = [[False, True, False], [False, False, True], [True, False, False]]
    examples = [frame(2), frame(3), c3, families(4)[0]]
    cases = 0
    for outer in examples:
        for inner in examples:
            out = substitution(outer, inner)
            require(triangle_count(out) == len(outer) * triangle_count(inner)
                    + len(inner)**3 * triangle_count(outer), "substitution triangle law failed")
            cases += 1
    squared = substitution(c3, c3)
    observed_h = hamilton_dp(squared)
    require(observed_h == 3159, "C3[C3] path hostile changed")
    return {"triangle_cases_checked": cases,
            "C3_C3": {"vertices": 9, "directed_triangles": triangle_count(squared),
                      "hamilton_paths": observed_h, "false_naive_product": 81}}


def balanced_ten_control() -> dict:
    out = [[i != j and (j - i) % 11 in range(1, 6) for j in range(10)] for i in range(10)]
    require(sorted(map(sum, out)) == [4] * 5 + [5] * 5, "balanced10 scores failed")
    require(triangle_count(out) == 40, "balanced10 triangle maximum failed")
    return {"construction": "cyclic regular tournament11 restricted to vertices0..9",
            "scores": sorted(map(sum, out)), "directed_triangles": 40}


def run() -> dict:
    fermat = [{"r": r, "tournament_order": 2**r + 2,
               "hamilton_paths_from_proved_formula": 2**(2**r) + 1} for r in range(6)]
    require(fermat[5]["hamilton_paths_from_proved_formula"] == 641 * 6700417,
            "Fermat5 factorization failed")
    return {"status": "FINITE-EXACT; universal proofs and cited dependencies in note",
            "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
            "identities": identity_controls(), "path_cubes": cube_controls(),
            "strong_families": family_controls(), "substitution": substitution_controls(),
            "balanced_ten": balanced_ten_control(), "Fermat_subsequence": fermat,
            "n10_fundamental_cycle_length_counts": {ell: 11 - ell for ell in range(3, 11)},
            "checks": "PASS; explicit checks survive python -O",
            "limitations": ["No arithmetic pairwise orientation inferred from scalar labels",
                            "No Collatz convergence claim", "No primality sufficiency claim",
                            "No enumeration of all 2^36 order10 path completions"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_suffix(".json"))
    args = parser.parse_args()
    result = run()
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8", newline="\n")
    print(json.dumps({"checks": result["checks"], "output": str(args.output),
                      "path_cube_count": result["path_cubes"]["total_completions"]}, sort_keys=True))
