#!/usr/bin/env python3
"""Exact checks for THM-2786.

Dependency-free integer verification; no truth-bearing ``assert`` nodes.
"""

from itertools import combinations


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def is_prime(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def next_odd_prime(n):
    p = max(3, n)
    while not (p % 2 and is_prime(p)):
        p += 1
    return p


def difference_map(marks):
    result = {}
    for i, j in combinations(range(len(marks)), 2):
        d = marks[j] - marks[i]
        if d in result:
            return None
        result[d] = (i, j)
    return result


def sum_map(marks):
    result = {}
    for i in range(len(marks)):
        for j in range(i, len(marks)):
            s = marks[i] + marks[j]
            if s in result and result[s] != (i, j):
                return None
            result[s] = (i, j)
    return result


def v2(n):
    value = 0
    while n % 2 == 0:
        n //= 2
        value += 1
    return value, n


def binary_marks(n):
    return tuple(2**i - 1 for i in range(n))


def erdos_turan_marks(n):
    p = next_odd_prime(n)
    marks = tuple(2 * p * i + (i * i % p) for i in range(n))
    return p, marks


def is_golomb(marks):
    return difference_map(marks) is not None


def optimal_small_ruler(n):
    lower = n * (n - 1) // 2
    length = lower
    while True:
        witnesses = []
        for inside in combinations(range(1, length), n - 2):
            marks = (0,) + inside + (length,)
            if is_golomb(marks):
                witnesses.append(marks)
        if witnesses:
            return length, tuple(witnesses)
        length += 1


def main():
    binary_checks = 0
    binary_spans = {}
    for n in range(2, 14):
        marks = binary_marks(n)
        differences = difference_map(marks)
        require(differences is not None, f"binary difference collision at n={n}")
        for d, (i, j) in differences.items():
            valuation, odd = v2(d)
            require(valuation == i, "binary valuation failed to recover lower endpoint")
            require(odd + 1 == 2 ** (j - i),
                    "binary odd quotient failed to recover endpoint gap")
            binary_checks += 1
        binary_spans[n] = marks[-1]

    et_checks = 0
    et_rows = {}
    for n in range(2, 51):
        p, marks = erdos_turan_marks(n)
        require(n <= p < 2 * n, f"Bertrand range failed at n={n}, p={p}")
        require(all(marks[i] < marks[i + 1] for i in range(n - 1)),
                "Erdos-Turan marks are not increasing")
        require(sum_map(marks) is not None, "Erdos-Turan Sidon sum collision")
        differences = difference_map(marks)
        require(differences is not None, "Erdos-Turan difference collision")
        require(marks[-1] < 2 * p * p <= 8 * n * n,
                "Erdos-Turan span invoice failed")
        et_checks += len(differences)
        if n <= 12:
            et_rows[n] = (p, marks[-1])

    small_optima = {}
    small_witnesses = {}
    for n in range(2, 7):
        length, witnesses = optimal_small_ruler(n)
        small_optima[n] = length
        small_witnesses[n] = witnesses[0]
        require(length >= n * (n - 1) // 2, "counting lower bound failed")
        require(is_golomb(witnesses[0]), "stored small witness is not Golomb")
    require(small_optima == {2: 1, 3: 3, 4: 6, 5: 11, 6: 17},
            f"unexpected small Golomb optima: {small_optima}")

    tree_debts = {}
    for n in range(2, 14):
        graceful_span = n - 1
        binary_span = binary_spans[n]
        p, marks = erdos_turan_marks(n)
        tree_debts[n] = (graceful_span, binary_span, marks[-1])

    print("THM-2786 BINARY GOLOMB / GRACEFUL COMPRESSION AUDIT")
    print(f"binary_pair_decoder_checks_n2_to_n13={binary_checks}")
    print(f"binary_spans_n2_to_n13={binary_spans}")
    print(f"erdos_turan_pair_checks_n2_to_n50={et_checks}")
    print(f"erdos_turan_prime_span_rows_n2_to_n12={et_rows}")
    print(f"exact_optimal_golomb_lengths_n2_to_n6={small_optima}")
    print(f"first_optimal_witnesses_n2_to_n6={small_witnesses}")
    print(f"tree_compression_debt_graceful_binary_ET_n2_to_n13={tree_debts}")
    print("K_n_lower_bound=binom(n,2); equality_holds_n<=4_and_fails_n=5,6")
    print("scope=global_pair_separation_not_tree_graceful_compression")
    print("PASS")


if __name__ == "__main__":
    main()
