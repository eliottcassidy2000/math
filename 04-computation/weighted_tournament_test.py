"""
Test: Does the signed position identity hold for weighted tournaments?
i.e., T[a][b] + T[b][a] = 1 but T[a][b] can be any real in [0,1]?

If yes: the identity is a POLYNOMIAL identity in arc variables subject to
T[a][b]+T[b][a]=1, and can be proved by algebraic methods.

If no: it requires the Boolean constraint T[a][b] in {0,1}.

Instance: opus-2026-03-05-S4
"""

from itertools import permutations, combinations
import random


def h_end_weighted(T, verts, target):
    """Sum over paths of product of weights, ending at target."""
    total = 0.0
    for perm in permutations(verts):
        if perm[-1] != target:
            continue
        prod = 1.0
        for k in range(len(perm) - 1):
            prod *= T[perm[k]][perm[k + 1]]
        total += prod
    return total


def h_start_weighted(T, verts, source):
    """Sum over paths of product of weights, starting at source."""
    total = 0.0
    for perm in permutations(verts):
        if perm[0] != source:
            continue
        prod = 1.0
        for k in range(len(perm) - 1):
            prod *= T[perm[k]][perm[k + 1]]
        total += prod
    return total


def test_weighted(n, num_tests=500):
    print(f"=== Weighted tournament test at n={n} ===")
    max_err = 0.0
    passes = 0

    for trial in range(num_tests):
        T = [[0.0]*n for _ in range(n)]
        for a in range(n):
            for b in range(a+1, n):
                w = random.random()  # T[a][b] in [0,1]
                T[a][b] = w
                T[b][a] = 1.0 - w  # tournament constraint

        i, j = 0, 1
        others = [x for x in range(n) if x != i and x != j]
        m = len(others)

        alt_sum = 0.0
        for size in range(m + 1):
            for S in combinations(others, size):
                S_set = set(S)
                R = [x for x in others if x not in S_set]
                S_list = list(S)
                sign = (-1) ** size

                ei = h_end_weighted(T, S_list + [i], i)
                bj = h_start_weighted(T, [j] + R, j)
                ej = h_end_weighted(T, S_list + [j], j)
                bi = h_start_weighted(T, [i] + R, i)

                alt_sum += sign * (ei * bj - ej * bi)

        err = abs(alt_sum)
        max_err = max(max_err, err)
        if err < 1e-10:
            passes += 1

    print(f"  Result: {passes}/{num_tests} pass (max error = {max_err:.2e})")
    return max_err < 1e-8


def test_antisymmetric(n, num_tests=200):
    """Test with T[a][b] = 1/2 + epsilon * A[a][b] where A is antisymmetric."""
    print(f"=== Antisymmetric perturbation test at n={n} ===")
    max_err = 0.0

    for trial in range(num_tests):
        T = [[0.0]*n for _ in range(n)]
        for a in range(n):
            for b in range(a+1, n):
                T[a][b] = 0.5 + random.uniform(-0.5, 0.5)
                T[b][a] = 1.0 - T[a][b]

        i, j = 0, 1
        others = list(range(2, n))
        m = len(others)

        alt_sum = 0.0
        for size in range(m + 1):
            for S in combinations(others, size):
                S_set = set(S)
                R = [x for x in others if x not in S_set]
                S_list = list(S)
                sign = (-1) ** size

                ei = h_end_weighted(T, S_list + [i], i)
                bj = h_start_weighted(T, [j] + R, j)
                ej = h_end_weighted(T, S_list + [j], j)
                bi = h_start_weighted(T, [i] + R, i)

                alt_sum += sign * (ei * bj - ej * bi)

        max_err = max(max_err, abs(alt_sum))

    print(f"  Max error: {max_err:.2e} (should be ~0 if identity holds)")


if __name__ == "__main__":
    random.seed(42)
    for n in [4, 5, 6]:
        test_weighted(n, num_tests=500 if n <= 5 else 100)

    print()
    for n in [4, 5, 6]:
        test_antisymmetric(n, num_tests=200 if n <= 5 else 50)
