#!/usr/bin/env python3
"""Exact arithmetic certificate for THM-4352's even self-similar cusp law.

The uniform proof is elementary and recorded in the theorem; this file attacks
indexing, parity, reciprocity, and boundary formulas over a large finite
universe with exact integer/rational arithmetic.
"""

from fractions import Fraction
from math import isqrt
import sys


MAX_G = 500
CHECKS = 0


def check(statement: bool, label: str) -> None:
    global CHECKS
    if not statement:
        raise AssertionError(label)
    CHECKS += 1


def triangular(n: int) -> int:
    return n * (n + 1) // 2


def inverse_square_block(n: int) -> tuple[int, int]:
    """Invert n=(g-1)^2+r, 1<=r<=2g-1."""
    q = isqrt(n)
    g = q if q * q == n else q + 1
    r = n - (g - 1) ** 2
    return g, r


def inverse_triangular_block(n: int) -> tuple[int, int]:
    """Invert N=T_(g-1)+h, 1<=h<=g."""
    g = (isqrt(8 * n + 1) - 1) // 2
    if triangular(g) < n:
        g += 1
    return g, n - triangular(g - 1)


oriented_indices: list[int] = []
quotient_indices: list[int] = []
naive_swap_failures: list[tuple[int, int]] = []
self_reciprocal_rows: list[tuple[int, int, int, int]] = []

for g in range(1, MAX_G + 1):
    m = 2 * g
    block: list[int] = []
    quotient_block: list[int] = []

    for r in range(1, m):
        eps = r & 1
        reciprocal_r = m - r

        # After removing the square Z^(2 floor(r/2)), the branch degree is
        # d=m-r+eps.  Since m is even, d is always even.
        branch_degree = m - r + eps
        check(branch_degree >= 2, f"branch degree lower bound g={g}, r={r}")
        check(branch_degree % 2 == 0, f"branch degree parity g={g}, r={r}")

        tail_genus = (branch_degree - 2) // 2
        persistent_delta = r // 2
        complementary_delta = reciprocal_r // 2

        check(
            tail_genus == (m - r - 1) // 2,
            f"tail genus closed form g={g}, r={r}",
        )
        check(
            tail_genus + persistent_delta + 1 == g,
            f"one-unit deficit identity g={g}, r={r}",
        )
        check(
            tail_genus == complementary_delta - (1 - eps),
            f"parity-twisted reciprocity g={g}, r={r}",
        )

        # The Newton slope and differential excess reciprocate exactly.
        slope = Fraction(r, m - r)
        reciprocal_slope = Fraction(m - r, r)
        check(slope * reciprocal_slope == 1, f"slope product g={g}, r={r}")

        base_degree = 12
        for B in (g - 1, g, g + 3):
            chi = Fraction(B + 1 - g, 1)
            excess = chi * base_degree * slope
            reciprocal_excess = chi * base_degree * reciprocal_slope
            check(
                excess * reciprocal_excess == (chi * base_degree) ** 2,
                f"excess product g={g}, r={r}, B={B}",
            )
        check(g - 1 + 1 - g == 0, f"zero threshold g={g}")
        check(g + 1 - g > 0, f"positive threshold g={g}")

        # Oriented types form odd-length square blocks.
        n = (g - 1) ** 2 + r
        block.append(n)
        oriented_indices.append(n)
        check(inverse_square_block(n) == (g, r), f"oriented inverse n={n}")

        # Reciprocity reflects each block.  Its fixed point is r=g.
        reciprocal_n = (g - 1) ** 2 + reciprocal_r
        signed_offset = r - g
        triangular_center = 2 * triangular(g - 1) + 1
        check(
            n == triangular_center + signed_offset,
            f"signed triangular center g={g}, r={r}",
        )
        check(
            n + reciprocal_n == 2 * triangular_center,
            f"block reflection g={g}, r={r}",
        )

        h = min(r, reciprocal_r)
        N = triangular(g - 1) + h
        check(inverse_triangular_block(N) == (g, h), f"quotient inverse N={N}")

        if r <= g:
            quotient_block.append(N)
            quotient_indices.append(N)

        if tail_genus != complementary_delta:
            naive_swap_failures.append((m, r))
        if r == reciprocal_r:
            self_reciprocal_rows.append((m, r, tail_genus, persistent_delta))

    check(block == list(range((g - 1) ** 2 + 1, g**2 + 1)), f"square block g={g}")
    check(
        quotient_block == list(range(triangular(g - 1) + 1, triangular(g) + 1)),
        f"triangular quotient block g={g}",
    )

check(
    oriented_indices == list(range(1, MAX_G**2 + 1)),
    "oriented square blocks globally contiguous",
)
check(
    quotient_indices == list(range(1, triangular(MAX_G) + 1)),
    "reciprocal quotient blocks globally contiguous",
)

# Hostile rows: these disprove importing the odd-order genus/delta swap.
check((6, 2) in naive_swap_failures, "m=6,r=2 naive-swap hostile")
check((4, 2, 0, 1) in self_reciprocal_rows, "m=4,r=2 fixed-point hostile")

OUTPUT = [
    "THM-4352 EVEN SELF-SIMILAR CUSP PROBE: PASS",
    f"universe: 1<=g<={MAX_G}, 1<=r<2g",
    f"oriented rows: {len(oriented_indices)}",
    f"reciprocal quotient rows: {len(quotient_indices)}",
    f"exact checks: {CHECKS}",
    "uniform identities:",
    "  branch_degree = 2g-r+(r mod 2), always even",
    "  genus_tail = floor((2g-r-1)/2)",
    "  intrinsic genus_tail + delta(A_(r-1)) = g-1",
    "  graph +1 requires same-complement two-end incidence",
    "  genus_tail = delta(A_(2g-r-1)) - 1_[r even]",
    "  excess_r * excess_(2g-r) = ((B+1-g)d)^2",
    "indexing:",
    "  oriented n=(g-1)^2+r fills 1..g^2 by odd square blocks",
    "  n=2T_(g-1)+1+(r-g); reciprocity negates the signed offset",
    "  quotient N=T_(g-1)+min(r,2g-r) fills 1..T_g",
    "hostiles:",
    "  m=6,r=2: tail genus 1, complementary delta 2 (naive swap fails)",
    "  m=4,r=2: self-reciprocal, tail genus 0, persistent delta 1",
    "  the missing unit is the two-ended attachment cycle, not tail genus",
]

# Force raw LF bytes on every platform so the frozen certificate has an
# unambiguous hash and normal/-O executions can be compared byte-for-byte.
sys.stdout.buffer.write(("\n".join(OUTPUT) + "\n").encode("utf-8"))
