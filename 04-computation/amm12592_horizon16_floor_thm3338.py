#!/usr/bin/env python3
"""Exact audit for THM-3338: horizon-16 cross-shell AMM surgery.

The rule stops at the first disagreement for every critical value 1..15
except n=5.  On n=5 it reads nine continuation bits, through flip 15, and
uses a lexicographic subset of each continuation Hamming class.  The exact
counts are hard-coded below.  For n>=16 the rule retains THM-2225 unchanged.

This program uses only standard-library integer arithmetic and checks:

  * causality of every finite-prefix verdict;
  * exact bisection of every nonconstant length-16 Hamming layer;
  * the equivalent Bernoulli polynomial identity;
  * equality with the aggregate THM-2225 checksum prefix;
  * the nonzero dyadic-shell deficits and their global cancellation.

Reproduce:
    python 04-computation/amm12592_horizon16_floor_thm3338.py
    python -O 04-computation/amm12592_horizon16_floor_thm3338.py
"""

from collections import defaultdict
from hashlib import sha256
from itertools import product
from math import comb
from pathlib import Path


M = 16
DEADLINE = {
    1: 2, 2: 3, 3: 4, 4: 5, 5: 15, 6: 7, 7: 8, 8: 9,
    9: 10, 10: 11, 11: 12, 12: 13, 13: 14, 14: 15, 15: 16,
}

# Away from the donor n=5, the verdict depends only on the initial bit.
# Empty means tails on both branches; {0,1} means heads on both branches.
HEAD_INITIAL_BITS = {
    1: {1}, 2: {1}, 3: {0}, 4: {0},
    6: {1}, 7: {1}, 8: {0}, 9: {1}, 10: {1},
    11: {0, 1}, 12: set(), 13: set(), 14: {0, 1}, 15: {0},
}

# At n=5, after 000001 or 111110, read positions 7..15.  In the class with
# k ones, declare heads on the first COUNT[initial][k] words in lexicographic
# order.  Every entry lies in [0, binom(9,k)].
COUNT = {
    0: (1, 9, 19, 25, 10, 0, 0, 9, 9, 0),
    1: (1, 0, 27, 84, 126, 116, 58, 15, 0, 0),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def critical_value(word: tuple[int, ...]) -> int | None:
    initial = word[0]
    for index, bit in enumerate(word[1:], start=1):
        if bit != initial:
            return index
    return None


def lexicographic_weight_rank(bits: tuple[int, ...]) -> int:
    """Zero-based lexicographic rank among same-weight binary words."""
    rank = 0
    remaining_ones = sum(bits)
    for position, bit in enumerate(bits):
        remaining_places = len(bits) - position - 1
        if bit == 1:
            # All words obtained by putting 0 here precede this word.
            if 0 <= remaining_ones <= remaining_places:
                rank += comb(remaining_places, remaining_ones)
            remaining_ones -= 1
    return rank


def horizon16_head(word: tuple[int, ...]) -> bool:
    n = critical_value(word)
    require(n is not None and 1 <= n < M, "finite rule outside its prefix")
    initial = word[0]
    if n != 5:
        return initial in HEAD_INITIAL_BITS[n]
    continuation = word[6:15]
    k = sum(continuation)
    return lexicographic_weight_rank(continuation) < COUNT[initial][k]


def shell_length(n: int) -> int:
    length = 2
    while length < n + 1:
        length *= 2
    return length


def checksum_head(word: tuple[int, ...]) -> bool:
    """THM-2225 cyclic-checksum verdict on n<16."""
    n = critical_value(word)
    require(n is not None and 1 <= n < M, "checksum outside prefix")
    length = shell_length(n)
    if length == 2:
        return word[:2] == (0, 1)
    half = length // 2
    tail = word[half:length]
    residue = sum((i + 1) * bit for i, bit in enumerate(tail)) % half
    return residue < half // 2


def words16():
    return product((0, 1), repeat=M)


def audit_boxes_and_ranks() -> None:
    for initial in (0, 1):
        for k, count in enumerate(COUNT[initial]):
            require(0 <= count <= comb(9, k), f"invalid count side={initial}, k={k}")
            ranks = [lexicographic_weight_rank(bits)
                     for bits in product((0, 1), repeat=9) if sum(bits) == k]
            require(sorted(ranks) == list(range(comb(9, k))),
                    f"rank map failed side={initial}, k={k}")


def audit_causality() -> None:
    seen: dict[tuple[int, tuple[int, ...]], bool] = {}
    for word in words16():
        n = critical_value(word)
        if n is None:
            continue
        prefix = word[:DEADLINE[n]]
        key = (n, prefix)
        verdict = horizon16_head(word)
        if key in seen:
            require(seen[key] == verdict, f"causality failure at n={n}")
        seen[key] = verdict
    require(set(n for n, _ in seen) == set(range(1, M)), "missing branch")


def layer_table(rule, predicate=lambda n: True):
    table = defaultdict(lambda: [0, 0])
    for word in words16():
        n = critical_value(word)
        if n is None or not predicate(n):
            continue
        weight = sum(word)
        table[weight][1] += 1
        if rule(word):
            table[weight][0] += 1
    return dict(sorted(table.items()))


def doubled_deficits(table):
    return {weight: 2 * heads - total for weight, (heads, total) in table.items()}


def power_coefficients(table):
    """Coefficients after q=1-p in sum heads_w p^(16-w)q^w."""
    coefficients = [0] * (M + 1)
    for ones, (heads, _total) in table.items():
        zeros = M - ones
        for j in range(ones + 1):
            coefficients[zeros + j] += heads * comb(ones, j) * (-1) ** j
    return coefficients


def target_coefficients():
    doubled = [0] * (M + 1)
    doubled[0] = 1
    doubled[M] -= 1
    for j in range(M + 1):
        doubled[j] -= comb(M, j) * (-1) ** j
    require(all(value % 2 == 0 for value in doubled), "target parity failure")
    return [value // 2 for value in doubled]


def main() -> None:
    audit_boxes_and_ranks()
    audit_causality()

    table = layer_table(horizon16_head)
    checksum = layer_table(checksum_head)
    expected = {k: [comb(M, k) // 2, comb(M, k)] for k in range(1, M)}
    require(table == expected, "new prefix does not bisect all Hamming layers")
    require(checksum == expected, "THM-2225 prefix control drifted")
    require(power_coefficients(table) == target_coefficients(),
            "Bernoulli polynomial identity failed")

    shells = ((1,), (2, 3), tuple(range(4, 8)), tuple(range(8, 16)))
    deficits = []
    for shell in shells:
        current = doubled_deficits(layer_table(horizon16_head, lambda n, s=shell: n in s))
        deficits.append(tuple(current.get(k, 0) for k in range(1, M)))
    require(any(any(value for value in vector) for vector in deficits),
            "all shells unexpectedly balanced")
    require(all(sum(vector[k] for vector in deficits) == 0 for k in range(M - 1)),
            "shell deficits do not cancel")

    floor_hits = tuple(n for n in range(1, M) if DEADLINE[n] == n + 1)
    require(floor_hits == tuple(n for n in range(1, M) if n != 5),
            "deadline profile drifted")

    source_hash = sha256(Path(__file__).read_bytes()).hexdigest()
    print("status=THM-3338_VERIFIED_EXACT")
    print("deadline_vector_T1_to_T15=" + ",".join(str(DEADLINE[n]) for n in range(1, M)))
    print("floor_hits=" + repr(floor_hits))
    print("donor_n5_counts_side0=" + repr(COUNT[0]))
    print("donor_n5_counts_side1=" + repr(COUNT[1]))
    print("layer_heads=" + repr({k: table[k][0] for k in table}))
    for shell, vector in zip(shells, deficits):
        print("shell_" + "_".join(map(str, shell)) + "_doubled_deficit=" + repr(vector))
    print("target=(1-p^16-(1-p)^16)/2")
    print("checksum_prefix_control=PASS")
    print("causality=PASS")
    print("source_sha256=" + source_hash)


if __name__ == "__main__":
    main()
