#!/usr/bin/env python3
"""Exact referee for THM-3428's rough maximal-order rank-seven exclusion."""

from __future__ import annotations

from hashlib import sha256
from math import gcd, isqrt


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(label)


def least_prime_factor(n: int) -> int:
    require(n >= 2, ("least prime factor domain", n))
    if n % 2 == 0:
        return 2
    for divisor in range(3, isqrt(n) + 1, 2):
        if n % divisor == 0:
            return divisor
    return n


def half_mask(modulus: int, residue: int) -> int:
    mask = 0
    for sheet in range(modulus):
        value = residue * (2 * sheet + 1) % (2 * modulus)
        distance = min(value, 2 * modulus - value)
        if 14 * distance < 2 * modulus:
            mask |= 1 << sheet
    return mask


def reflected_mask(mask: int, modulus: int) -> int:
    out = 0
    for sheet in range(modulus):
        if mask >> sheet & 1:
            out |= 1 << ((-1 - sheet) % modulus)
    return out


def unique_maximal_masks(modulus: int, parity: int) -> list[tuple[int, int]]:
    by_mask: dict[int, int] = {}
    for residue in range(1, modulus + 1):
        if residue % 2 == parity and gcd(residue, modulus) == 1:
            mask = half_mask(modulus, residue)
            by_mask.setdefault(mask, residue)
    return sorted((residue, mask) for mask, residue in by_mask.items())


def maximum_disjoint_with_base(items: list[tuple[int, int]], base: int) -> tuple[int, tuple[int, ...]]:
    candidates = [(residue, mask) for residue, mask in items if mask & base == 0]
    best: tuple[int, ...] = ()

    def search(index: int, used: int, chosen: tuple[int, ...]) -> None:
        nonlocal best
        if len(chosen) + len(candidates) - index <= len(best):
            return
        if index == len(candidates):
            if len(chosen) > len(best):
                best = chosen
            return
        residue, mask = candidates[index]
        if not used & mask:
            search(index + 1, used | mask, chosen + (residue,))
        search(index + 1, used, chosen)

    search(0, 0, ())
    return 1 + len(best), (1,) + best


def first_rough_composites(residue_class: int, count: int) -> list[int]:
    out: list[int] = []
    candidate = 513
    if candidate % 2 == 0:
        candidate += 1
    while len(out) < count:
        factor = least_prime_factor(candidate)
        if candidate % 14 == residue_class and factor > 7 and factor < candidate:
            out.append(candidate)
        candidate += 2
    return out


def check_control(modulus: int) -> tuple[int, int, int, int, tuple[int, ...], int, int]:
    require(modulus >= 512 and modulus % 2 == 1, ("control range", modulus))
    factor = least_prime_factor(modulus)
    require(factor > 7 and factor < modulus, ("rough composite", modulus, factor))
    k, residue_class = divmod(modulus, 14)
    require(residue_class in (1, 3, 5, 9, 11, 13), ("residue class", modulus))
    fixed_sheet = (modulus - 1) // 2
    fixed_bit = 1 << fixed_sheet
    even_items = unique_maximal_masks(modulus, 0)
    odd_items = unique_maximal_masks(modulus, 1)
    require(even_items and odd_items, ("parity banks", modulus))

    expected_odd = 2 * k if residue_class in (1, 3, 5) else 2 * k + 2
    for residue, mask in even_items:
        require(mask.bit_count() == 2 * k + 1, ("even size", modulus, residue))
        require(mask & fixed_bit, ("even fixed sheet", modulus, residue))
        require(reflected_mask(mask, modulus) == mask, ("even reflection", modulus, residue))
    for residue, mask in odd_items:
        require(mask.bit_count() == expected_odd, ("odd size", modulus, residue))
        require(not mask & fixed_bit, ("odd fixed sheet", modulus, residue))
        require(reflected_mask(mask, modulus) == mask, ("odd reflection", modulus, residue))

    odd_base = half_mask(modulus, 1)
    odd_packing, odd_witness = maximum_disjoint_with_base(odd_items, odd_base)
    require(odd_packing == 3, ("odd packing", modulus, odd_packing, odd_witness))

    even_base = half_mask(modulus, 2) & ~fixed_bit
    even_collisions = 0
    for residue, mask in even_items:
        if even_base & (mask & ~fixed_bit):
            even_collisions += 1
    if residue_class == 1:
        require(k >= 13, ("even gap range", modulus, k))
        require(even_collisions == len(even_items),
                ("even off-fixed collision", modulus, even_collisions, len(even_items)))

    if residue_class in (3, 5):
        forced_disjoint = 0
    elif residue_class == 1:
        forced_disjoint = 6
    elif residue_class == 9:
        forced_disjoint = 4
    elif residue_class == 11:
        forced_disjoint = 5
    else:
        forced_disjoint = 6
    if forced_disjoint:
        require(forced_disjoint > odd_packing,
                ("failed contradiction", modulus, forced_disjoint, odd_packing))
    return (modulus, factor, residue_class, k, odd_witness,
            len(even_items), len(odd_items))


def check_positive_and_scope_controls() -> tuple[tuple[int, tuple[int, ...]], tuple[int, tuple[int, ...], tuple[int, ...]]]:
    prime_modulus = 29
    prime_witness = (1, 5, 7, 8, 12, 13, 22)
    covered = 0
    for residue in prime_witness:
        covered |= half_mask(prime_modulus, residue)
    require(covered == (1 << prime_modulus) - 1, "order-29 positive boundary")
    require(all(gcd(residue, prime_modulus) == 1 for residue in prime_witness),
            "order-29 maximal orders")

    lower_modulus = 513
    lower_witness = (1, 57, 285, 342, 399)
    covered = 0
    for residue in lower_witness:
        covered |= half_mask(lower_modulus, residue)
    require(covered == (1 << lower_modulus) - 1, "lower-order scope hostile")
    orders = tuple(lower_modulus // gcd(lower_modulus, residue) for residue in lower_witness)
    require(orders == (513, 9, 9, 3, 9), ("lower-order profile", orders))
    return (prime_modulus, prime_witness), (lower_modulus, lower_witness, orders)


def check_even_gap_arithmetic() -> int:
    cells = 0
    for k in range(13, 2001):
        modulus = 14 * k + 1
        require(modulus // (k + 1) <= 13, ("even shortest gap", k, modulus))
        cells += 1
    return cells


def main() -> None:
    records = []
    for residue_class in (1, 3, 5, 9, 11, 13):
        for modulus in first_rough_composites(residue_class, 2):
            records.append(check_control(modulus))
    even_gap_cells = check_even_gap_arithmetic()
    prime_boundary, lower_scope = check_positive_and_scope_controls()
    semantic_payload = (tuple(records), even_gap_cells, prime_boundary, lower_scope)
    semantic = sha256(repr(semantic_payload).encode("ascii")).hexdigest()
    print(f"rough_composite_controls={len(records)}")
    print("control_records=" + ";".join(
        f"Q{q}:spf{factor}:r{residue}:k{k}:odd{odd}:E{even_count}:O{odd_count}"
        for q, factor, residue, k, odd, even_count, odd_count in records
    ))
    print(f"even_shortest_gap_cells={even_gap_cells}")
    print(f"small_positive_boundary={prime_boundary}")
    print(f"lower_quotient_order_scope_hostile={lower_scope}")
    print(f"semantic_sha256={semantic}")
    print("rough_maximal_order_rank7_exclusion=PASS")


if __name__ == "__main__":
    main()
