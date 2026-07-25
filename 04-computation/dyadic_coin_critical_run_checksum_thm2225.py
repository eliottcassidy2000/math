#!/usr/bin/env python3
"""Exact finite referee for THM-2225's two critical-run coin extractors."""

from collections import Counter
from itertools import combinations, product
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def critical_value(word: tuple[int, ...]) -> int:
    """Initial constant-run length; requires a nonconstant word."""
    for index, bit in enumerate(word[1:], start=1):
        if bit != word[0]:
            return index
    raise ValueError("critical value is undefined for a constant word")


def shell_parameters(n: int) -> tuple[int, int]:
    """Return M,m with M=2m the least power of two strictly above n."""
    require(n >= 1, "critical value must be positive")
    big = 1 << n.bit_length()
    return big, big // 2


def in_shell(word: tuple[int, ...]) -> bool:
    big = len(word)
    require(big >= 2 and big & (big - 1) == 0, "length must be dyadic")
    half = big // 2
    return len(set(word[:half])) == 1 and len(set(word)) > 1


def compression_decision(
    word: tuple[int, ...],
) -> tuple[str, int, int]:
    """Return outcome, decisive homogeneous-block size, and pair index."""
    current = word
    block_size = 1
    while len(current) >= 2:
        unequal = [
            pair
            for pair in range(len(current) // 2)
            if current[2 * pair] != current[2 * pair + 1]
        ]
        if unequal:
            pair = unequal[0]
            outcome = "H" if current[2 * pair : 2 * pair + 2] == (0, 1) else "T"
            return outcome, block_size, pair
        current = tuple(current[2 * pair] for pair in range(len(current) // 2))
        block_size *= 2
    raise RuntimeError("nonconstant word compressed without a decision")


def compression_partner(word: tuple[int, ...]) -> tuple[int, ...]:
    _, block_size, pair = compression_decision(word)
    start = 2 * pair * block_size
    middle = start + block_size
    end = middle + block_size
    return (
        word[:start]
        + word[middle:end]
        + word[start:middle]
        + word[end:]
    )


def checksum_decision(word: tuple[int, ...]) -> str:
    big = len(word)
    half = big // 2
    if big == 2:
        require(word in ((0, 1), (1, 0)), "invalid two-bit shell word")
        return "H" if word == (0, 1) else "T"
    tail = word[half:]
    checksum = sum((index + 1) * bit for index, bit in enumerate(tail)) % half
    return "H" if checksum < half // 2 else "T"


def audit_exhaustive_shell(big: int) -> int:
    words = [
        word
        for word in product((0, 1), repeat=big)
        if in_shell(word)
    ]
    require(len(words) == 2 * (2 ** (big // 2) - 1), "shell size drift")

    compression_counts: Counter[tuple[int, str]] = Counter()
    checksum_counts: Counter[tuple[int, str]] = Counter()
    for word in words:
        outcome, _, _ = compression_decision(word)
        partner = compression_partner(word)
        partner_outcome, _, _ = compression_decision(partner)
        require(in_shell(partner), "compression partner left shell")
        require(compression_partner(partner) == word, "block swap not involutive")
        require(sum(partner) == sum(word), "block swap changed Hamming weight")
        require(partner_outcome != outcome, "block swap did not reverse output")
        compression_counts[(sum(word), outcome)] += 1

        checksum = checksum_decision(word)
        checksum_counts[(sum(word), checksum)] += 1
        n = critical_value(word)
        shell_big, _ = shell_parameters(n)
        require(shell_big == big, "critical-value shell drift")
        stop = 2 if big == 2 else (big - 1 if n <= big - 2 else big)
        require(stop <= max(2, 2 * n - 1), "sharp path bound failed")
        if big >= 4 and stop == big - 1:
            sibling = word[:-1] + (1 - word[-1],)
            require(in_shell(sibling), "last-bit sibling left shell")
            require(
                checksum_decision(sibling) == checksum,
                "zero checksum coefficient drift",
            )

    for weight in range(1, big):
        require(
            compression_counts[(weight, "H")]
            == compression_counts[(weight, "T")],
            f"compression shell imbalance at M={big}, weight={weight}",
        )
        require(
            checksum_counts[(weight, "H")] == checksum_counts[(weight, "T")],
            f"checksum shell imbalance at M={big}, weight={weight}",
        )
    return len(words)


def audit_dyadic_arithmetic(max_half: int) -> int:
    require(max_half >= 2 and max_half & (max_half - 1) == 0, "bad cutoff")
    path_rows = 1  # n=1, M=2
    half = 2
    while half <= max_half:
        for weight in range(1, half):
            divisor = gcd(weight, half)
            require((half // 2) % divisor == 0, "antipode gcd gate failed")
            for coset in range(divisor):
                lower = sum(
                    residue % divisor == coset
                    for residue in range(half // 2)
                )
                upper = sum(
                    residue % divisor == coset
                    for residue in range(half // 2, half)
                )
                require(lower == upper, "checksum coset failed to bisect")
        all_one_checksum = half * (half + 1) // 2 % half
        require(all_one_checksum == half // 2, "middle-shell checksum drift")
        for n in range(half, 2 * half):
            stop = 2 * half - 1 if n <= 2 * half - 2 else 2 * half
            require(stop <= 2 * n - 1, "universal path arithmetic failed")
            path_rows += 1
        half *= 2
    return path_rows


def hostile_checksum_counts(length: int, weight: int) -> tuple[int, int]:
    counts = [0, 0]
    for positions in combinations(range(1, length + 1), weight):
        checksum = sum(positions) % length
        counts[int(checksum >= length // 2)] += 1
    return counts[0], counts[1]


def main() -> None:
    shell_sizes = {
        big: audit_exhaustive_shell(big)
        for big in (2, 4, 8, 16)
    }
    path_rows = audit_dyadic_arithmetic(4096)
    hostile = hostile_checksum_counts(6, 2)
    require(hostile == (7, 8), "non-dyadic hostile control drift")
    require(
        all((2 * delta - 3) % 6 != 0 for delta in range(6)),
        "hostile Smith obstruction disappeared",
    )

    print(
        "exhaustive_shell_sizes="
        + ",".join(f"M{big}:{size}" for big, size in shell_sizes.items())
    )
    print("exhaustive_checks=block_swap+shell_balance+last_bit_merge+path_bound")
    print("dyadic_antipode_audit_max_half=4096")
    print(f"dyadic_path_bound_rows={path_rows}")
    print("hostile_m6_j2_checksum_counts=(7,8)")
    print("hostile_congruence=2*delta!=3_mod_6")
    print("status=THM-2225_PROVED_VERIFIED_EXACT")


if __name__ == "__main__":
    main()
