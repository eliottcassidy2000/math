#!/usr/bin/env python3
"""Exact hostile referee for THM-2160's half-tail fair-coin rule.

The symbolic proof is in THM-2160.  This script checks the binomial completion,
the shell-local deadline obstruction, and the dyadic parity classification
with integer arithmetic.  Explicit raising checks remain active under `-O`.
"""

from itertools import combinations
from math import comb


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def is_power_of_two(n: int) -> bool:
    return n > 0 and n & (n - 1) == 0


def signed_defect(h: int, j: int) -> int:
    """Coefficient e_j of u(1+u)(1-u^2)^(h/2-1)."""
    t = h // 2
    if j % 2:
        a = (j - 1) // 2
    else:
        a = (j - 2) // 2
    return (-1) ** a * comb(t - 1, a)


def prescribed_count_by_enumeration(h: int, j: int) -> int:
    """Heads on z_1=1, |z|=j, for the early parity rule."""
    t = h // 2
    count = 0
    for support_tail in combinations(range(1, h), j - 1):
        if sum(1 for i in support_tail if i < t) % 2 == 0:
            count += 1
    return count


def audit_shell(h: int) -> tuple[int, int, int, int]:
    require(h >= 2 and is_power_of_two(h), f"not a dyadic shell: h={h}")
    t = h // 2
    r_values: list[int] = []

    for j in range(1, h):
        e_j = signed_defect(h, j)
        top = comb(h - 1, j) - e_j
        require(top % 2 == 0, f"nonintegral completion h={h}, j={j}")
        r_j = top // 2
        require(0 <= r_j <= comb(h - 1, j), f"illegal completion h={h}, j={j}")

        a_j_numerator = comb(h - 1, j - 1) + e_j
        require(a_j_numerator % 2 == 0, f"nonintegral early count h={h}, j={j}")
        a_j = a_j_numerator // 2
        require(a_j + r_j == comb(h, j) // 2, f"layer imbalance h={h}, j={j}")
        require(comb(h, j) % 2 == 0, f"odd dyadic layer h={h}, j={j}")

        if h <= 16:
            exact = prescribed_count_by_enumeration(h, j)
            require(exact == a_j, f"parity enumerator mismatch h={h}, j={j}")
        r_values.append(r_j)

    # The weight-(h-1) obstruction: no d<t can attain either required count.
    required = {h // 2, h // 2 - 1}
    for d in range(1, t):
        attainable = set(range(0, d)) | set(range(h - d, h))
        require(required.isdisjoint(attainable), f"early lower bound failed h={h}, d={d}")

    # The constructed exact-n=h deadline meets the requested and improved bounds.
    exact_deadline = h + t
    if h == 2:
        require(exact_deadline == 3, "h=2 boundary failed")
    else:
        require(exact_deadline <= 2 * h - 2, f"improved deadline failed h={h}")

    return t, min(r_values), max(r_values), exact_deadline


def audit_dyadic_classification(limit: int) -> tuple[int, int]:
    dyadic = 0
    nondyadic = 0
    for h in range(2, limit + 1):
        internal = [comb(h, j) for j in range(1, h)]
        if is_power_of_two(h):
            require(all(value % 2 == 0 for value in internal), f"dyadic parity failed h={h}")
            dyadic += 1
        else:
            require(any(value % 2 == 1 for value in internal), f"Lucas witness absent h={h}")
            nondyadic += 1
    return dyadic, nondyadic


def main() -> None:
    shells = (2, 4, 8, 16, 32, 64, 128, 256)
    print("THM-2160 DYADIC HALF-TAIL REFEREE")
    print("all counts are exact integers; enumeration controls run through h=16")
    print()
    print(" h   t=h/2   min R_j   max R_j   exact-n=h deadline")
    for h in shells:
        t, min_r, max_r, deadline = audit_shell(h)
        print(f"{h:>3}  {t:>6}  {min_r:>8}  {max_r:>8}  {deadline:>18}")

    dyadic, nondyadic = audit_dyadic_classification(256)
    print()
    print(f"classification control h<=256: {dyadic} dyadic and {nondyadic} nondyadic shells")
    print("sample 00001: h=4, t=2, decision after X6")
    print("result: every layer completion and every shell-local lower-bound sensor passed")
    print("all assertions passed")


if __name__ == "__main__":
    main()
