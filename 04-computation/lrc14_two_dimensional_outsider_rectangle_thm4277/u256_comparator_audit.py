#!/usr/bin/env python3
"""Detached exact audit of the C++ 128x128 -> 256 gap comparator."""

from __future__ import annotations

import hashlib

MASK32 = (1 << 32) - 1
MASK64 = (1 << 64) - 1
MASK127 = (1 << 127) - 1
CASES = 100_000


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def multiply_limbs(left: int, right: int) -> tuple[int, ...]:
    require(0 <= left <= MASK127 and 0 <= right <= MASK127,
            "factor outside nonnegative signed-128 range")
    a = [(left >> (32 * i)) & MASK32 for i in range(4)]
    b = [(right >> (32 * i)) & MASK32 for i in range(4)]
    product = [0] * 8
    for i in range(4):
        carry = 0
        for j in range(4):
            k = i + j
            current = a[i] * b[j] + product[k] + carry
            require(current <= MASK64, "C++ uint64 accumulator bound failed")
            product[k] = current & MASK32
            carry = current >> 32
        k = i + 4
        while carry:
            require(k < 8, "limb carry overflow")
            current = product[k] + carry
            require(current <= MASK64, "carry accumulator bound failed")
            product[k] = current & MASK32
            carry = current >> 32
            k += 1
    return tuple(product)


def limbs_to_int(limbs: tuple[int, ...]) -> int:
    require(len(limbs) == 8, "wrong limb count")
    return sum(limb << (32 * i) for i, limb in enumerate(limbs))


def limbs_less(left: tuple[int, ...], right: tuple[int, ...]) -> bool:
    for i in range(7, -1, -1):
        if left[i] != right[i]:
            return left[i] < right[i]
    return False


class SplitMix64:
    def __init__(self, state: int) -> None:
        self.state = state & MASK64

    def next(self) -> int:
        self.state = (self.state + 0x9E3779B97F4A7C15) & MASK64
        value = self.state
        value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
        value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & MASK64
        return value ^ (value >> 31)

    def signed128_nonnegative(self) -> int:
        return self.next() | ((self.next() & ((1 << 63) - 1)) << 64)


def check_product(left: int, right: int) -> tuple[int, ...]:
    limbs = multiply_limbs(left, right)
    require(limbs_to_int(limbs) == left * right, "product reconstruction failed")
    return limbs


def main() -> None:
    digest = hashlib.sha256()
    boundaries = (0, 1, MASK32, 1 << 32, (1 << 64) - 1,
                  1 << 64, 1 << 110, 1 << 120, MASK127)
    boundary_checks = 0
    for left in boundaries:
        for right in boundaries:
            check_product(left, right)
            boundary_checks += 1

    generator = SplitMix64(0x4255494C44553235)
    for _ in range(CASES):
        a = generator.signed128_nonnegative()
        b = generator.signed128_nonnegative()
        c = generator.signed128_nonnegative()
        d = generator.signed128_nonnegative()
        ab = check_product(a, b)
        cd = check_product(c, d)
        require(limbs_less(ab, cd) == (a * b < c * d),
                "lexicographic product comparison failed")
        for value in (a, b, c, d):
            digest.update(value.to_bytes(16, "little"))
        digest.update(bytes([a * b < c * d]))

    print("LRC14_RECTANGLE_U256_COMPARATOR_AUDIT_V1")
    print(f"BOUNDARY_PRODUCTS {boundary_checks}")
    print(f"DETERMINISTIC_RANDOM_COMPARISONS {CASES}")
    print(f"STREAM_SHA256 {digest.hexdigest()}")
    print("ACCUMULATOR_MAX_LE_2^64_MINUS_1 PASS")
    print("VERDICT PASS EXACT_PYTHON_BIGINT_ORACLE")


if __name__ == "__main__":
    main()
