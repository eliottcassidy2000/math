#!/usr/bin/env python3
"""Exact scout for ordinal task ranks and the normal-strip top Kummer sieve.

This is a FINITE-EXACT companion to a research reflection, not a theorem
dependency.  The global statements used in the reflection have elementary
integer/UFD proofs; this program freezes their boundary cases and the first
composite transverse depth n=6.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from math import gcd
import sys


Polynomial = tuple[Fraction, ...]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def triangular(z: int) -> int:
    return z * (z + 1) // 2


def task_rank(n: int, j: int) -> int:
    """One-based triangular address for 1 <= j < n, n >= 2."""
    require(n >= 2 and 1 <= j < n, "normal-strip cell outside triangular cone")
    return triangular(n - 2) + j


def decode_task_rank(rank: int) -> tuple[int, int]:
    require(rank >= 1, "task rank must be positive")
    n = 2
    while triangular(n - 1) < rank:
        n += 1
    return n, rank - triangular(n - 2)


def trim(poly: Polynomial) -> Polynomial:
    data = list(poly)
    while len(data) > 1 and data[-1] == 0:
        data.pop()
    return tuple(data)


def poly_add(left: Polynomial, right: Polynomial) -> Polynomial:
    size = max(len(left), len(right))
    data = [Fraction(0) for _ in range(size)]
    for index, value in enumerate(left):
        data[index] += value
    for index, value in enumerate(right):
        data[index] += value
    return trim(tuple(data))


def poly_scale(scalar: Fraction, poly: Polynomial) -> Polynomial:
    return trim(tuple(scalar * value for value in poly))


def poly_mul(left: Polynomial, right: Polynomial) -> Polynomial:
    data = [Fraction(0) for _ in range(len(left) + len(right) - 1)]
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            data[i + j] += left_value * right_value
    return trim(tuple(data))


def poly_pow(poly: Polynomial, exponent: int) -> Polynomial:
    require(exponent >= 0, "negative polynomial exponent")
    result: Polynomial = (Fraction(1),)
    base = trim(poly)
    power = exponent
    while power:
        if power & 1:
            result = poly_mul(result, base)
        base = poly_mul(base, base)
        power //= 2
    return result


def poly_derivative(poly: Polynomial) -> Polynomial:
    if len(poly) == 1:
        return (Fraction(0),)
    return trim(tuple(Fraction(i) * poly[i] for i in range(1, len(poly))))


def top_bucket(n: int, j: int, w: Polynomial, c: Polynomial) -> Polynomial:
    """Coefficient of z^(n+j-1) in J(A,C) for tops wz^n, cz^j."""
    return poly_add(
        poly_scale(Fraction(n), poly_mul(w, poly_derivative(c))),
        poly_scale(Fraction(-j), poly_mul(poly_derivative(w), c)),
    )


def kummer_exponents(n: int, j: int) -> tuple[int, int]:
    divisor = gcd(n, j)
    return n // divisor, j // divisor


def divisor_count(n: int) -> int:
    return sum(1 for divisor in range(1, n + 1) if n % divisor == 0)


def main() -> None:
    sys.stdout.reconfigure(newline="\n")

    centered_checks = 0
    odd_square_checks = 0
    for z in range(-2000, 2001):
        for h in range(1, 10):
            require(
                triangular(z + h) - triangular(z - h) == h * (2 * z + 1),
                "signed triangular centered identity failed",
            )
            centered_checks += 1
        require(
            triangular(z + 2) - triangular(z - 2) == 4 * z + 2,
            "four-step signed triangular identity failed",
        )
    for rank in range(1, 2001):
        require(
            (2 * rank - 1) ** 2 == 8 * triangular(rank - 1) + 1,
            "odd-square ordinal peel failed",
        )
        odd_square_checks += 1

    # The cells (n,j), 1<=j<n, pack consecutively by increasing n.  This is
    # a scheduler; adjacent cells need not share a structural type.
    address_checks = 0
    addresses: list[int] = []
    conductor_edge_checks = 0
    for n in range(2, 101):
        for j in range(1, n):
            rank = task_rank(n, j)
            require(decode_task_rank(rank) == (n, j), "task-rank round trip failed")
            addresses.append(rank)
            address_checks += 1
    require(addresses == list(range(1, triangular(99) + 1)), "task ranks have a gap")
    # Through depth N the cells are literally the edges {j,n} of K_N.  This
    # identifies their count with THM-3745's ordinary N-fold-point delta,
    # but carries no Jacobian or conductor predicate by itself.
    for size in range(2, 101):
        cells = {(j, n) for n in range(2, size + 1) for j in range(1, n)}
        edges = {(left, right) for left in range(1, size) for right in range(left + 1, size + 1)}
        require(cells == edges, "normal-strip cells do not match complete-graph edges")
        require(len(cells) == triangular(size - 1), "triangular edge count failed")
        conductor_edge_checks += 1

    # THM-3756 uses the identical ambient address on Pythagorean ordinal
    # pairs, but its primitive filter is gcd(2n-1,2j-1)=1.  The JC top-row
    # shear filter is instead j|n.  Rank eight is the first Pythagorean hole
    # and simultaneously a genuine (not shear-removable) quintic task.
    cross_predicate_checks = 0
    first_pythagorean_hole: tuple[int, int, int] | None = None
    for rank in range(1, triangular(40) + 1):
        n, j = decode_task_rank(rank)
        pythagorean_primitive = gcd(2 * n - 1, 2 * j - 1) == 1
        jc_shear = n % j == 0
        if not pythagorean_primitive and first_pythagorean_hole is None:
            first_pythagorean_hole = (rank, n, j)
        # Merely execute both predicates at each shared address; equality is
        # explicitly neither expected nor asserted.
        require(isinstance(pythagorean_primitive, bool), "primitive predicate not Boolean")
        require(isinstance(jc_shear, bool), "JC shear predicate not Boolean")
        cross_predicate_checks += 1
    require(first_pythagorean_hole == (8, 5, 2), "shared-address hostile moved")
    require(5 % 2 != 0, "rank-eight JC cell stopped being genuine")

    # At a UFD prime p, the top ODE forces n*v_p(c)=j*v_p(w).
    # Coprimality of n/g and j/g gives the primitive Kummer exponent ray.
    valuation_checks = 0
    shear_checks = 0
    for n in range(2, 61):
        for j in range(1, n):
            w_exponent, c_exponent = kummer_exponents(n, j)
            require(gcd(w_exponent, c_exponent) == 1, "Kummer ray not primitive")
            for valuation_w in range(0, 121):
                for valuation_c in range(0, 121):
                    solves = n * valuation_c == j * valuation_w
                    parameterized = (
                        valuation_w % w_exponent == 0
                        and valuation_c % c_exponent == 0
                        and valuation_w // w_exponent
                        == valuation_c // c_exponent
                    )
                    require(solves == parameterized, "valuation ray classification failed")
                    valuation_checks += 1
            # A constant target shear removes the A-top precisely when the
            # primitive c exponent is one, equivalently j divides n.
            require((c_exponent == 1) == (n % j == 0), "shear criterion failed")
            shear_checks += 1

    # Freeze the polynomial identity on several non-squarefree scale factors.
    polynomial_checks = 0
    scales: list[Polynomial] = [
        (Fraction(1), Fraction(1)),
        (Fraction(2), Fraction(-3), Fraction(1)),
        (Fraction(-6), Fraction(1), Fraction(1)),
        (Fraction(1), Fraction(0), Fraction(2), Fraction(1)),
    ]
    for n in range(2, 21):
        for j in range(1, n):
            w_exponent, c_exponent = kummer_exponents(n, j)
            for scale in scales:
                w = poly_scale(Fraction(2), poly_pow(scale, w_exponent))
                c = poly_scale(Fraction(-3), poly_pow(scale, c_exponent))
                require(top_bucket(n, j, w, c) == (Fraction(0),), "top ODE failed")
                if n % j == 0:
                    quotient = n // j
                    # w/c^q is the constant 2/(-3)^q, so the leading term of
                    # A-rho*C^q vanishes exactly.
                    rho = Fraction(2, (-3) ** quotient)
                    require(
                        poly_add(w, poly_scale(-rho, poly_pow(c, quotient)))
                        == (Fraction(0),),
                        "target-shear top cancellation failed",
                    )
                polynomial_checks += 1

    sextic = []
    for j in range(1, 6):
        w_exponent, c_exponent = kummer_exponents(6, j)
        sextic.append(
            {
                "rank": task_rank(6, j),
                "j": j,
                "gcd": gcd(6, j),
                "ray": (w_exponent, c_exponent),
                "shear_power": 6 // j if 6 % j == 0 else None,
                "status": "shear_to_depth_le_5" if 6 % j == 0 else "genuine_sextic",
            }
        )
    require([row["rank"] for row in sextic] == [11, 12, 13, 14, 15], "sextic ranks moved")
    require(
        [row["j"] for row in sextic if row["status"] == "genuine_sextic"] == [4, 5],
        "sextic genuine-row sieve failed",
    )
    require(
        task_rank(6, 4) == task_rank(6, 3) + 1,
        "ordinal-successor hostile moved",
    )

    genuine_counts = {
        n: sum(1 for j in range(1, n) if n % j != 0) for n in range(2, 31)
    }
    require(
        all(genuine_counts[n] == n - divisor_count(n) for n in genuine_counts),
        "genuine-row count formula failed",
    )

    semantic = {
        "centered_checks": centered_checks,
        "odd_square_checks": odd_square_checks,
        "address_checks": address_checks,
        "conductor_edge_checks": conductor_edge_checks,
        "cross_predicate_checks": cross_predicate_checks,
        "first_pythagorean_hole": first_pythagorean_hole,
        "valuation_checks": valuation_checks,
        "shear_checks": shear_checks,
        "polynomial_checks": polynomial_checks,
        "sextic": sextic,
        "genuine_counts": genuine_counts,
    }
    digest = sha256(repr(semantic).encode("ascii")).hexdigest()

    print("JC2 NORMAL-STRIP ORDINAL/GCD SIDECAR PROBE")
    print(f"signed_triangular_centered_checks={centered_checks}")
    print(f"odd_square_ordinal_checks={odd_square_checks}")
    print(f"triangular_task_address_checks={address_checks}; max_rank={addresses[-1]}")
    print(
        f"complete_graph_conductor_alignment_checks={conductor_edge_checks}; "
        "cells_through_N=edges(K_N)=T(N-1)"
    )
    print(
        f"shared_pythagorean_address_checks={cross_predicate_checks}; "
        "rank8=(5,2):odd_root_gcd=3:Pythagorean_hole:JC_genuine_quintic"
    )
    print(f"ufd_valuation_ray_checks={valuation_checks}")
    print(f"target_shear_criterion_checks={shear_checks}")
    print(f"top_bucket_polynomial_checks={polynomial_checks}")
    for row in sextic:
        print(
            "sextic_cell="
            f"rank{row['rank']}:j{row['j']}:gcd{row['gcd']}:"
            f"ray{row['ray']}:shear{row['shear_power']}:status={row['status']}"
        )
    print("ordinal_successor_hostile=rank13:(6,3):shear -> rank14:(6,4):genuine")
    print(f"genuine_row_counts_n2_to_30={genuine_counts}")
    print(f"semantic_sha256={digest}")
    print("PASS")


if __name__ == "__main__":
    main()
