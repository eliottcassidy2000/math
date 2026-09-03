#!/usr/bin/env python3
"""Exact audit of a divisor-complete body escaping clocks 15..43.

All inequalities are integer inequalities.  Tail residues are the odd classes
modulo 2N, and pairs are unordered with repetition, matching the literal
two-lift semantics of THM-2066/THM-4347.
"""

from itertools import combinations_with_replacement

H = (1, 5, 11, 13, 17, 19, 23, 37, 41, 70, 72)


def abs_residue(x: int, modulus: int) -> int:
    r = x % modulus
    return min(r, modulus - r)


def packet(n: int) -> tuple[int, ...]:
    return tuple(
        r for r in range(n)
        if all(14 * abs_residue(h * r, n) >= n for h in H)
    )


def danger(z: int, r: int, lift: int, n: int) -> bool:
    return 14 * abs_residue(z * (r + lift * n), 2 * n) < 2 * n


def covering_pairs(n: int, a: tuple[int, ...]) -> list[tuple[int, int]]:
    odds = range(1, 2 * n, 2)
    return [
        (u, v)
        for u, v in combinations_with_replacement(odds, 2)
        if all(danger(u, r, j, n) or danger(v, r, j, n)
               for r in a for j in (0, 1))
    ]


def main() -> None:
    covered_divisors = tuple(
        d for d in range(2, 15) if any(h % d == 0 for h in H)
    )
    print(f"body={H}")
    print(f"distinct={len(set(H)) == 11}")
    print(f"covered_divisors={covered_divisors}")
    print("clock packet_size packet least_cover covering_pair_count")
    fixed_escape = True
    for n in range(15, 44):
        a = packet(n)
        covers = covering_pairs(n, a)
        fixed_escape &= bool(covers)
        print(n, len(a), a, covers[0] if covers else None, len(covers))
    print(f"fixed_bank_15_43_escape={fixed_escape}")

    first_repair = None
    for n in range(44, 501):
        a = packet(n)
        if a and not covering_pairs(n, a):
            first_repair = (n, a)
            break
    print(f"first_nonempty_repair_44_500={first_repair}")


if __name__ == "__main__":
    main()
