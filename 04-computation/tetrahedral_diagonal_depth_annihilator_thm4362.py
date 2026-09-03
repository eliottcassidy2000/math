#!/usr/bin/env python3
"""Exact finite checks for THM-4362's universal depth annihilator.

The theorem itself is proved generator-by-generator by a binomial finite-
difference identity.  This dependency-free certificate enumerates the exact
THM-4308 monomial columns over a hostile range, verifies the sharp boundary
generator, and identifies the row-nine and row-ten stencils used later.
"""

from __future__ import annotations

from math import comb
import sys


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def columns(depth: int, max_row: int):
    """Canonical columns (a,b,c,e) of pi_max_row(P_depth)."""

    for b in range(depth + 1):
        for a in range(depth - b + 1):
            for e in range(max_row // 2 + 1):
                for c in range(max_row + 1):
                    if b + c + 2 * e <= max_row:
                        yield a, b, c, e


def stencil(max_row: int) -> tuple[int, ...]:
    return tuple(
        (-1) ** (n - 5) * comb(max_row + 3 - n, 3)
        for n in range(5, max_row + 1)
    )


def on_column(max_row: int, column: tuple[int, int, int, int]) -> int:
    """Apply L_m to one expanded source-normal monomial column."""

    a, b, c, e = column
    if a != 2 * c + 3 * e - 10:
        return 0
    total = 0
    n0 = b + c + 2 * e
    for k in range(c + e + 1):
        n = n0 + k
        if 5 <= n <= max_row:
            total += (
                comb(c + e, k)
                * (-1) ** (n - 5)
                * comb(max_row + 3 - n, 3)
            )
    return total


def hostile(depth: int) -> tuple[int, int, int, int]:
    """One-depth-beyond generator for m=d+6."""

    parity = depth % 2
    c = (10 + depth - 3 * parity) // 2
    return depth, 0, c, parity


def triangular(z: int) -> int:
    return z * (z + 1) // 2


def main() -> None:
    # The finite-difference core: N>=4 annihilates every cubic binomial
    # coefficient.  The proof uses [s^3] s^N(1+s)^(M-N)=0.
    difference_checks = 0
    for n_terms in range(4, 31):
        for top in range(n_terms, n_terms + 31):
            value = sum(
                (-1) ** k * comb(n_terms, k) * comb(top - k, 3)
                for k in range(n_terms + 1)
            )
            check(value == 0, f"finite difference N={n_terms} M={top}")
            difference_checks += 1

    # Enumerate every canonical generator through a substantial exact range.
    # The universal proof in the theorem handles arbitrary m,d.
    enumerated_columns = 0
    for max_row in range(7, 19):
        check(stencil(max_row)[-1] in (-1, 1), f"primitive m={max_row}")
        for depth in range(max_row - 6):
            for column in columns(depth, max_row):
                check(on_column(max_row, column) == 0, f"annihilation m={max_row} d={depth} col={column}")
                enumerated_columns += 1

    # Sharpness: at d=m-6 an explicit generator has value (-1)^(d+7).
    hostile_checks = 0
    for depth in range(31):
        max_row = depth + 6
        column = hostile(depth)
        a, b, c, e = column
        check(a + b == depth, f"hostile depth d={depth}")
        check(b + c + 2 * e <= max_row, f"hostile retained d={depth}")
        check(a == 2 * c + 3 * e - 10, f"hostile diagonal d={depth}")
        check(on_column(max_row, column) == (-1) ** (depth + 7), f"hostile value d={depth}")
        hostile_checks += 1

    # The two consecutive stencils used by THM-4358 and proposed THM-4361.
    check(stencil(9) == (35, -20, 10, -4, 1), "row-nine stencil")
    check(stencil(10) == (56, -35, 20, -10, 4, -1), "row-ten stencil")
    check(sum(1 for _ in columns(2, 9)) == 160, "pi9 P2 columns")
    check(sum(1 for _ in columns(3, 10)) == 304, "pi10 P3 columns")

    # Figurate-number bridge: tetrahedral numbers are cumulative triangular
    # numbers, and the user's integer extension of T obeys its centered law.
    for r in range(3, 51):
        check(comb(r, 3) == sum(triangular(q) for q in range(1, r - 1)), f"tetrahedral r={r}")
    for z in range(-50, 51):
        check(triangular(-z) == triangular(z - 1), f"integer triangle z={z}")
        check(triangular(z + 2) - triangular(z - 2) == 4 * z + 2, f"centered triangle z={z}")

    print("THM-4362 TETRAHEDRAL DIAGONAL DEPTH ANNIHILATOR: PASS")
    print("L_m=sum_(n=5)^m (-1)^(n-5) C(m+3-n,3) h_(n,2n-10)")
    print("UNIVERSAL_RANGE=m>=7 and d<=m-7")
    print("SHARP_HOSTILE=d=m-6; (a,b,c,e)=(d,0,(10+d-3(d mod 2))/2,d mod 2)")
    print("SHARP_HOSTILE_VALUE=(-1)^(d+7)")
    print("M9_P2_STENCIL=(35,-20,10,-4,1); columns=160")
    print("M10_P3_STENCIL=(56,-35,20,-10,4,-1); columns=304")
    print("FIGURATE=C(r,3)=sum_(q=1)^(r-2) T(q); fourth finite difference vanishes")
    print(f"FINITE_DIFFERENCE_CHECKS={difference_checks}")
    print(f"ENUMERATED_ZERO_COLUMNS={enumerated_columns}")
    print(f"SHARP_HOSTILES={hostile_checks}")
    print(f"CHECKS={CHECKS}")
    print("SCOPE finite projected depth modules only; no all-row membership, seam entry, Keller pair, JC2, or DC2")


if __name__ == "__main__":
    main()
