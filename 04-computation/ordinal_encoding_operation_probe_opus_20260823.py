#!/usr/bin/env python3
"""Cross-frontier exact probe for operation-compatible ordinal encodings.

This is a reflection companion, not a theorem dependency.  It checks the
algebraic consequences used in the 2026-08-23 ordinal-encoding synthesis and
prints the hostiles that prevent untyped transfers.
"""

from __future__ import annotations

from collections import defaultdict
from hashlib import sha256
import sys


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def T(z: int) -> int:
    return z * (z + 1) // 2


def binom2(z: int) -> int:
    return z * (z - 1) // 2


def quadratic(a: int, b: int, c: int, z: int) -> int:
    return a * binom2(z) + b * z + c


def lrc_profile_rank(b: int, c: int) -> int:
    require(5 <= c <= 19 and 1 <= b < c, "profile outside THM-2349 ledger")
    return T(c - 2) - 6 + b


def odd_rank_product(r: int, s: int) -> int:
    return 2 * r * s - r - s + 1


def walk_address(k: int, endpoint: int) -> int:
    require(-k <= endpoint <= k, "endpoint outside ternary-walk layer")
    return k * k + endpoint + k


def exponent_rank(a: int, b: int) -> int:
    require(a >= 0 and b >= 0, "negative monomial exponent")
    return T(a + b) + b


def monomial_bracket(
    left: tuple[int, int, int], right: tuple[int, int, int]
) -> tuple[tuple[int, int], int] | None:
    a, b, alpha = left
    c, d, beta = right
    coefficient = alpha * beta * (a * d - b * c)
    if coefficient == 0:
        return None
    return (a + c - 1, b + d - 1), coefficient


def main() -> None:
    # Keep the checked-in transcript byte-identical across operating systems.
    sys.stdout.reconfigure(newline="\n")

    # Every integer-valued quadratic has a linear centered difference in the
    # Newton basis.  The parity hostile has the same centered differences.
    quadratic_checks = 0
    periodic_kernel_checks = 0
    graph_laplacian_checks = 0
    for a in range(-4, 5):
        for b in range(-4, 5):
            for c0 in range(-2, 3):
                for z in range(-30, 31):
                    require(
                        quadratic(a, b, c0, z + 1)
                        - 2 * quadratic(a, b, c0, z)
                        + quadratic(a, b, c0, z - 1)
                        == a,
                        "constant forward second difference failed",
                    )
                    for step in range(1, 7):
                        difference = quadratic(a, b, c0, z + step) - quadratic(
                            a, b, c0, z - step
                        )
                        require(
                            difference == step * (a * (2 * z - 1) + 2 * b),
                            "quadratic centered difference failed",
                        )
                        require(
                            difference
                            == (
                                quadratic(a, b, c0, z + step) + (-1) ** (z + step)
                            )
                            - (
                                quadratic(a, b, c0, z - step) + (-1) ** (z - step)
                            ),
                            "parity-kernel hostile failed",
                        )
                        quadratic_checks += 1
                        periodic_kernel_checks += 1

                for height in range(-8, 9):
                    for child_count in range(0, 6):
                        for parent_count in range(0, 6):
                            centre = quadratic(a, b, c0, height)
                            direct = child_count * (
                                quadratic(a, b, c0, height + 1) - centre
                            ) + parent_count * (
                                quadratic(a, b, c0, height - 1) - centre
                            )
                            formula = (a * height + b) * (
                                child_count - parent_count
                            ) + a * parent_count
                            require(direct == formula, "graded graph Laplacian failed")
                            graph_laplacian_checks += 1

    # A quadratic prefix potential is exactly the prefix size of an affine-
    # width layered family, and local child moves have affine address change.
    layered_checks = 0
    for a in range(0, 6):
        for b in range(1, 8):
            for k in range(0, 60):
                qk = a * binom2(k) + b * k
                width = a * k + b
                require(
                    a * binom2(k + 1) + b * (k + 1) - qk == width,
                    "layered prefix/width law failed",
                )
                for slot in range(width):
                    for delta in range(-2, 3):
                        old_address = qk + slot
                        new_address = (
                            a * binom2(k + 1) + b * (k + 1) + slot + delta
                        )
                        require(
                            new_address - old_address == width + delta,
                            "layered child address law failed",
                        )
                        layered_checks += 1

    # THM-2349's 165 valuation profiles have a literal triangular shortlex
    # rank; its centered c-difference is blind to b and ordinal successor is
    # not a native row operation.
    profiles = [(b, c) for c in range(5, 20) for b in range(1, c)]
    profile_addresses = [lrc_profile_rank(*profile) for profile in profiles]
    require(len(profiles) == 165, "THM-2349 profile count moved")
    require(profile_addresses == list(range(1, 166)), "profile rank is not bijective")
    profile_difference_checks = 0
    for c in range(7, 18):
        for b in range(1, c - 2):
            require(
                lrc_profile_rank(b, c + 2) - lrc_profile_rank(b, c - 2)
                == 4 * c - 6,
                "profile centered difference failed",
            )
            profile_difference_checks += 1
    require(lrc_profile_rank(4, 5) == 4, "profile successor left control moved")
    require(lrc_profile_rank(1, 6) == 5, "profile successor right control moved")
    repeated_first_lane = [lrc_profile_rank(1, c) for c in range(5, 20)]
    require(
        repeated_first_lane
        == [(n * n + 5 * n - 4) // 2 for n in range(1, 16)],
        "repeated-first quadratic lane failed",
    )

    # The triangular centered-difference observer on THM-3713's balanced
    # deep offsets is exactly 4L+2A and has a rational anchored kernel.
    defect_checks = 0
    for seed in range(-4, 5):
        defect = {u: ((u + 2 * seed) % 7) - 3 for u in range(-6, 7)}
        augmentation = sum(defect.values())
        first_moment = sum(u * coefficient for u, coefficient in defect.items())
        triangular_flux = sum(
            (T(u + 2) - T(u - 2)) * coefficient
            for u, coefficient in defect.items()
        )
        require(
            triangular_flux == 4 * first_moment + 2 * augmentation,
            "triangular flux reduction failed",
        )
        defect_checks += 1
    flux_hostile = {1: 5, 2: -3}
    require(
        sum((T(u + 2) - T(u - 2)) * v for u, v in flux_hostile.items()) == 0,
        "triangular flux hostile moved",
    )
    positive_cone_control = {1: 2, 2: 3, -1: -5, -2: -7, -3: -11}
    require(
        sum(
            (T(u + 2) - T(u - 2)) * v
            for u, v in positive_cone_control.items()
        )
        > 0,
        "typed five-colour cone lost positivity",
    )

    # Odd-value selection is operation-compatible: multiplication pulls back
    # to a commutative associative law on odd ordinals.
    odd_product_checks = 0
    for r in range(1, 81):
        for s in range(1, 81):
            star = odd_rank_product(r, s)
            require(
                2 * star - 1 == (2 * r - 1) * (2 * s - 1),
                "odd-rank product law failed",
            )
            require(odd_rank_product(r, 1) == r, "odd-rank identity failed")
            for t in range(1, 9):
                require(
                    odd_rank_product(odd_rank_product(r, s), t)
                    == odd_rank_product(r, odd_rank_product(s, t)),
                    "odd-rank associativity failed",
                )
                odd_product_checks += 1

    # Linear odd widths cumulate to squares.  Endpoint rank preserves the
    # ternary-walk endpoint operation but loses word order and prefix barriers.
    walk_checks = 0
    walk_addresses: set[int] = set()
    for k in range(0, 101):
        for endpoint in range(-k, k + 1):
            address = walk_address(k, endpoint)
            walk_addresses.add(address)
            for epsilon in (-1, 0, 1):
                require(
                    walk_address(k + 1, endpoint + epsilon) - address
                    == 2 * k + 2 + epsilon,
                    "ternary-walk address update failed",
                )
                walk_checks += 1
    require(walk_addresses == set(range(10201)), "ternary-walk shell packing failed")
    walk_hostiles = ((1, -1), (-1, 1))
    require(sum(walk_hostiles[0]) == sum(walk_hostiles[1]) == 0, "walk endpoint hostile")
    require(
        min(0, walk_hostiles[0][0], sum(walk_hostiles[0]))
        != min(0, walk_hostiles[1][0], sum(walk_hostiles[1])),
        "walk prefix hostile failed",
    )

    # Diagonal triangular rank is a bijection on monomial exponents and
    # transports multiplication/Poisson landing.  Coefficient cancellation
    # remains a labelled convolution sidecar.
    exponent_checks = 0
    exponent_addresses: set[int] = set()
    for a in range(0, 31):
        for b in range(0, 31 - a):
            exponent_addresses.add(exponent_rank(a, b))
            for c in range(0, 16):
                for d in range(0, 16 - c):
                    h = a + b
                    hp = c + d
                    require(
                        exponent_rank(a + c, b + d)
                        == exponent_rank(a, b) + exponent_rank(c, d) + h * hp,
                        "monomial multiplication rank failed",
                    )
                    determinant = a * d - b * c
                    if determinant != 0:
                        require(
                            exponent_rank(a + c - 1, b + d - 1)
                            == exponent_rank(a, b)
                            + exponent_rank(c, d)
                            + h * hp
                            - 2 * (h + hp),
                            "Poisson landing rank failed",
                        )
                    exponent_checks += 1
    require(
        exponent_addresses == set(range(T(31))),
        "monomial diagonal rank failed to pack its finite triangle",
    )

    f_terms = [(1, 0, 1), (0, 1, 1)]
    g_terms = [(2, 0, 1), (1, 1, 2), (0, 2, 1)]
    raw_brackets: list[tuple[tuple[int, int], int]] = []
    aggregated: defaultdict[tuple[int, int], int] = defaultdict(int)
    for left in f_terms:
        for right in g_terms:
            term = monomial_bracket(left, right)
            if term is not None:
                raw_brackets.append(term)
                aggregated[term[0]] += term[1]
    require(raw_brackets, "Poisson cancellation hostile has no raw channels")
    require(all(value == 0 for value in aggregated.values()), "Poisson hostile did not cancel")

    semantic = {
        "quadratic_checks": quadratic_checks,
        "periodic_kernel_checks": periodic_kernel_checks,
        "graph_laplacian_checks": graph_laplacian_checks,
        "layered_checks": layered_checks,
        "profile_addresses": profile_addresses,
        "profile_difference_checks": profile_difference_checks,
        "repeated_first_lane": repeated_first_lane,
        "flux_hostile": flux_hostile,
        "odd_product_checks": odd_product_checks,
        "walk_checks": walk_checks,
        "walk_hostiles": walk_hostiles,
        "exponent_checks": exponent_checks,
        "poisson_raw": raw_brackets,
        "poisson_aggregated": sorted(aggregated.items()),
    }
    digest = sha256(repr(semantic).encode("ascii")).hexdigest()

    print("ORDINAL ENCODING OPERATION PROBE")
    print(f"integer_quadratic_centered_checks={quadratic_checks}")
    print(f"periodic_kernel_hostile_checks={periodic_kernel_checks}")
    print(f"graded_graph_laplacian_checks={graph_laplacian_checks}")
    print(f"affine_width_layered_address_checks={layered_checks}")
    print(f"lrc_profile_rank=1..{len(profile_addresses)}")
    print(f"lrc_centered_profile_checks={profile_difference_checks}")
    print("lrc_ordinal_successor_hostile=4:(4,5) -> 5:(1,6)")
    print(f"lrc_repeated_first_lane={repeated_first_lane}")
    print("triangular_flux=4L+2A; hostile=5*delta_1-3*delta_2 -> 0")
    print(f"odd_rank_product_checks={odd_product_checks}; star(r,s)=2rs-r-s+1")
    print(f"ternary_walk_address_checks={walk_checks}; shells_end_at_squares")
    print("ternary_walk_hostile=(+1,-1) versus (-1,+1): same endpoint, different prefix")
    print(f"monomial_rank_landing_checks={exponent_checks}")
    print(f"poisson_raw_channels={raw_brackets}")
    print(f"poisson_aggregated={sorted(aggregated.items())}")
    print(f"semantic_sha256={digest}")
    print("PASS")


if __name__ == "__main__":
    main()
