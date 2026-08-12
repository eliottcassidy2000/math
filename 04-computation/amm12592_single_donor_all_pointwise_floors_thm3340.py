#!/usr/bin/env python3
"""Exact audit for THM-3340: one donor gives every AMM pointwise floor.

For a dyadic horizon M, all critical values 2<=n<M stop at their first
disagreement: even n are heads on both initial-bit branches and odd n are
tails.  The n=1 branch reads through M and repairs the remaining layer
defect using a lexicographic subset of the 01-prefix words.  Constant
length-M prefixes retain the THM-2225 continuation.

The proof uses cyclic rotation.  This audit checks the exact coefficient
identities through M=512 and exhaustively checks the rotation maps, causal
rule, and all words for M=2,4,8,16.

Reproduce:
    python 04-computation/amm12592_single_donor_all_pointwise_floors_thm3340.py
    python -O 04-computation/amm12592_single_donor_all_pointwise_floors_thm3340.py
"""

from collections import defaultdict
from hashlib import sha256
from itertools import product
from math import comb
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def critical_value(word: tuple[int, ...]) -> int | None:
    initial = word[0]
    for index, bit in enumerate(word[1:], start=1):
        if bit != initial:
            return index
    return None


def branch_layer_count(M: int, n: int, weight: int) -> int:
    """Number of length-M, weight-w words with critical value n."""
    free = M - n - 1
    total = 0
    # Initial branch 0^n 1.
    if 0 <= weight - 1 <= free:
        total += comb(free, weight - 1)
    # Initial branch 1^n 0.
    if 0 <= weight - n <= free:
        total += comb(free, weight - n)
    return total


def layer_data(M: int, weight: int):
    require(M >= 2 and M & (M - 1) == 0, "M must be a power of two")
    require(1 <= weight < M, "nonconstant layer required")
    donor_capacity = comb(M - 2, weight - 1)  # one of the two n=1 branches
    even = sum(branch_layer_count(M, n, weight) for n in range(2, M, 2))
    odd = sum(branch_layer_count(M, n, weight) for n in range(3, M, 2))
    require(comb(M, weight) % 2 == 0, "Lucas parity failure")
    repair = comb(M, weight) // 2 - even
    defect = even - odd
    return donor_capacity, even, odd, defect, repair


def audit_coefficients(M: int) -> None:
    m = M - 2
    a = []
    for k in range(m + 1):
        value = sum(comb(2 * i, k - 1) for i in range(m // 2)
                    if 0 <= k - 1 <= 2 * i)
        a.append(value)
        if m >= 2 and k >= 1:
            require(2 * value + (a[k - 1] if k else 0) == comb(m, k),
                    f"A_m coefficient recurrence failed M={M}, k={k}")
    for weight in range(1, M):
        capacity, even, odd, defect, repair = layer_data(M, weight)
        k = weight - 1
        require(defect == a[k] + a[m - k],
                f"closed discrepancy factorization failed M={M}, k={k}")
        require(a[k] <= comb(m - 1, k) if k <= m - 1 else a[k] == 0,
                f"first hockey-stick bound failed M={M}, k={k}")
        require(a[m - k] <= comb(m - 1, k - 1) if k >= 1 else a[m] == 0,
                f"second hockey-stick bound failed M={M}, k={k}")
        require(even + odd + 2 * capacity == comb(M, weight),
                f"partition failure M={M}, w={weight}")
        require(defect % 2 == 0, f"defect parity failure M={M}, w={weight}")
        require(repair == capacity - defect // 2,
                f"repair identity failure M={M}, w={weight}")
        require(0 <= defect <= capacity,
                f"strong one-sided capacity bound failed M={M}, w={weight}")
        require(0 <= repair <= capacity,
                f"repair outside one donor branch M={M}, w={weight}")
        require(2 * repair >= capacity,
                f"repair used less than half the donor class M={M}, w={weight}")
        require(even + repair == comb(M, weight) // 2,
                f"layer bisection failed M={M}, w={weight}")


def left_rotate(word: tuple[int, ...], places: int = 1) -> tuple[int, ...]:
    places %= len(word)
    return word[places:] + word[:places]


def lexicographic_weight_rank(bits: tuple[int, ...]) -> int:
    rank = 0
    remaining_ones = sum(bits)
    for position, bit in enumerate(bits):
        remaining_places = len(bits) - position - 1
        if bit == 1:
            if 0 <= remaining_ones <= remaining_places:
                rank += comb(remaining_places, remaining_ones)
            remaining_ones -= 1
    return rank


def donor_head(word: tuple[int, ...], M: int) -> bool:
    """Verdict on a nonconstant length-M word under the finite rule."""
    n = critical_value(word)
    require(n is not None and n < M, "finite rule outside its event")
    if n >= 2:
        return n % 2 == 0
    # Only the 01 orientation is used for repair.  Within its continuation
    # weight class choose the first repair words lexicographically.
    if word[:2] != (0, 1):
        return False
    weight = sum(word)
    repair = layer_data(M, weight)[4]
    return lexicographic_weight_rank(word[2:]) < repair


def audit_rotation_and_rule(M: int) -> None:
    by_weight = defaultdict(list)
    for word in product((0, 1), repeat=M):
        if critical_value(word) is not None:
            by_weight[sum(word)].append(word)

    for weight, words in by_weight.items():
        odd = {word for word in words
               if (critical_value(word) or 0) >= 3 and critical_value(word) % 2 == 1}
        even_same = {word for word in words
                     if (critical_value(word) or 0) >= 2
                     and critical_value(word) % 2 == 0 and word[-1] == word[0]}
        even_opposite = {word for word in words
                         if (critical_value(word) or 0) >= 2
                         and critical_value(word) % 2 == 0 and word[-1] != word[0]}

        # Left rotation is a literal bijection O -> E_same.
        images = {left_rotate(word) for word in odd}
        require(images == even_same, f"odd/even rotation bijection failed M={M}, w={weight}")

        # Rotate an unmatched even word left by n-1.  The image begins 01 or
        # 10, and its terminal run recovers n-1, so this is injective.
        injected = set()
        for word in even_opposite:
            n = critical_value(word)
            image = left_rotate(word, n - 1)
            require(critical_value(image) == 1, "unmatched image missed donor")
            require(sum(image) == weight, "rotation changed weight")
            require(image not in injected, "unmatched injection collided")
            injected.add(image)

        capacity, even_count, odd_count, defect, repair = layer_data(M, weight)
        require(len(even_opposite) == defect == even_count - odd_count,
                "unmatched defect count failed")
        require(len(injected) <= 2 * capacity, "donor injection exceeded capacity")

        heads = sum(donor_head(word, M) for word in words)
        require(heads == comb(M, weight) // 2,
                f"exhaustive layer bisection failed M={M}, w={weight}")
        require(repair <= capacity, "one-sided donor overflow")


def main() -> None:
    for exponent in range(1, 10):
        audit_coefficients(1 << exponent)
    for M in (2, 4, 8, 16):
        audit_rotation_and_rule(M)

    samples = {}
    for M in (2, 4, 8, 16, 32):
        samples[M] = tuple(layer_data(M, weight)[4] for weight in range(1, M))

    source_hash = sha256(Path(__file__).read_bytes()).hexdigest()
    print("status=THM-3340_VERIFIED_EXACT")
    print("coefficient_horizons=2,4,8,16,32,64,128,256,512")
    print("exhaustive_rotation_horizons=2,4,8,16")
    for M, repair in samples.items():
        print(f"M{M}_one_sided_donor_counts=" + repr(repair))
    print("finite_profile=T_M(1)=M; T_M(n)=n+1 for 2<=n<M")
    print("rotation_bijection=PASS")
    print("unmatched_injection=PASS")
    print("strong_one_sided_half_capacity=PASS")
    print("closed_discrepancy_factorization=PASS")
    print("layer_bisection=PASS")
    print("source_sha256=" + source_hash)


if __name__ == "__main__":
    main()
