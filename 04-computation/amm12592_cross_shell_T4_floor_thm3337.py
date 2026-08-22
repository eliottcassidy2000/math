#!/usr/bin/env python3
"""Exact audit for THM-3337: a cross-shell AMM 12592 extractor with T(4)=5.

For every nonconstant binary word of length eight, let n be the initial
constant-run length.  The candidate rule on n <= 7 is

    n=1: H iff X1=0,             stop at 2;
    n=2: H iff X1=1,             stop at 3;
    n=3: H iff X5=0,             stop at 5;
    n=4: H iff X1=0,             stop at 5;
    n=5: H iff X7=1,             stop at 7;
    n=6: H iff X1=1,             stop at 7;
    n=7: H iff X1=0,             stop at 8.

This program checks, using only integer arithmetic:

  1. the stated decisions are causal at their deadlines;
  2. every Hamming layer of nonconstant eight-bit words is bisected;
  3. the old dyadic shells n in {2,3} and n in {4,5,6,7} are not separately
     balanced, but their doubled deficits are exact opposites;
  4. the resulting head polynomial is (1-p^8-(1-p)^8)/2;
  5. the THM-2225 cyclic-checksum prefix has the same aggregate polynomial.

The theorem then splices this prefix into THM-2225 and leaves all n >= 8
branches unchanged.  The program audits the finite replacement; THM-2225
supplies the already-proved infinite tail and its deadlines.

Reproduce:
    python 04-computation/amm12592_cross_shell_T4_floor_thm3337.py
    python -O 04-computation/amm12592_cross_shell_T4_floor_thm3337.py
"""

from collections import defaultdict
from hashlib import sha256
from math import comb
from pathlib import Path


DEADLINE = {1: 2, 2: 3, 3: 5, 4: 5, 5: 7, 6: 7, 7: 8}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def critical_value(word: tuple[int, ...]) -> int | None:
    """Initial constant-run length; None only for a constant word."""
    initial = word[0]
    for index, bit in enumerate(word[1:], start=1):
        if bit != initial:
            return index
    return None


def cross_shell_head(word: tuple[int, ...]) -> bool:
    n = critical_value(word)
    require(n is not None and 1 <= n <= 7, "rule called outside its prefix")
    if n == 1:
        return word[0] == 0
    if n == 2:
        return word[0] == 1
    if n == 3:
        return word[4] == 0
    if n == 4:
        return word[0] == 0
    if n == 5:
        return word[6] == 1
    if n == 6:
        return word[0] == 1
    return word[0] == 0  # n == 7


def shell_length(n: int) -> int:
    length = 2
    while length < n + 1:
        length *= 2
    return length


def checksum_head(word: tuple[int, ...]) -> bool:
    """THM-2225 cyclic-checksum verdict on the applicable dyadic prefix."""
    n = critical_value(word)
    require(n is not None and 1 <= n <= 7, "checksum called outside prefix")
    length = shell_length(n)
    if length == 2:
        return word[:2] == (0, 1)
    half = length // 2
    tail = word[half:length]
    residue = sum((i + 1) * bit for i, bit in enumerate(tail)) % half
    return residue < half // 2


def words8():
    for value in range(1 << 8):
        yield tuple((value >> (7 - i)) & 1 for i in range(8))


def audit_causality() -> None:
    # A deterministic stopping rule is causal exactly when completions of a
    # stopped prefix receive the same verdict whenever they retain that n.
    for n, deadline in DEADLINE.items():
        seen: dict[tuple[int, ...], bool] = {}
        for word in words8():
            if critical_value(word) != n:
                continue
            prefix = word[:deadline]
            verdict = cross_shell_head(word)
            if prefix in seen:
                require(seen[prefix] == verdict, f"causality failure n={n}")
            seen[prefix] = verdict
        require(seen, f"empty critical branch n={n}")


def layer_table(rule, predicate=lambda n: True):
    table = defaultdict(lambda: [0, 0])  # weight -> [heads,total]
    for word in words8():
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


def bernoulli_coefficients(table):
    """Power-basis coefficients of sum_k heads_k p^(8-k)(1-p)^k."""
    coefficients = [0] * 9
    for ones, (heads, _total) in table.items():
        zeros = 8 - ones
        for j in range(ones + 1):
            coefficients[zeros + j] += heads * comb(ones, j) * (-1) ** j
    return coefficients


def target_coefficients():
    # (1-p^8-(1-p)^8)/2 is integral coefficientwise because the two constant
    # endpoints cancel the odd endpoint coefficients of (1-p)^8.
    target = [0] * 9
    target[0] = 1
    target[8] -= 1
    for j in range(9):
        target[j] -= comb(8, j) * (-1) ** j
    require(all(value % 2 == 0 for value in target), "target not even")
    return [value // 2 for value in target]


def main() -> None:
    audit_causality()

    new_all = layer_table(cross_shell_head)
    old_all = layer_table(checksum_head)
    expected = {k: [comb(8, k) // 2, comb(8, k)] for k in range(1, 8)}
    require(new_all == expected, "cross-shell prefix does not bisect every layer")
    require(old_all == expected, "THM-2225 prefix control drifted")

    shell_1 = layer_table(cross_shell_head, lambda n: n == 1)
    shell_2 = layer_table(cross_shell_head, lambda n: 2 <= n <= 3)
    shell_4 = layer_table(cross_shell_head, lambda n: 4 <= n <= 7)
    d1 = doubled_deficits(shell_1)
    d2 = doubled_deficits(shell_2)
    d4 = doubled_deficits(shell_4)
    require(all(value == 0 for value in d1.values()), "first shell lost balance")
    require(any(value != 0 for value in d2.values()), "hostile imbalance vanished")
    require(any(value != 0 for value in d4.values()), "hostile imbalance vanished")
    require(all(d2[k] + d4[k] == 0 for k in range(1, 8)),
            "cross-shell deficits do not cancel")

    coefficients = bernoulli_coefficients(new_all)
    target = target_coefficients()
    require(coefficients == target, "Bernoulli polynomial identity failed")

    require(DEADLINE[4] == 5, "T(4) is not the claimed floor")
    require(all(DEADLINE[n] >= n + 1 for n in DEADLINE), "deadline precedes disagreement")

    source_hash = sha256(Path(__file__).read_bytes()).hexdigest()
    print("status=THM-3337_VERIFIED_EXACT")
    print("deadline_vector_T1_to_T7=" + ",".join(str(DEADLINE[n]) for n in range(1, 8)))
    print("all_layer_heads=" + repr({k: new_all[k][0] for k in new_all}))
    print("all_layer_totals=" + repr({k: new_all[k][1] for k in new_all}))
    print("shell_2_doubled_deficit=" + repr(d2))
    print("shell_4_doubled_deficit=" + repr(d4))
    print("power_basis_head_polynomial=" + repr(coefficients))
    print("target=(1-p^8-(1-p)^8)/2")
    print("checksum_prefix_control=PASS")
    print("causality=PASS")
    print("source_sha256=" + source_hash)


if __name__ == "__main__":
    main()
