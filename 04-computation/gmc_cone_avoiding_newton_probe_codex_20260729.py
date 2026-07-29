#!/usr/bin/env python3
"""Exact scout for the remaining consecutive cone-avoiding chart.

This is not a theorem companion.  It asks whether the Gregory--Newton
positivity mechanism of THM-2890 applies directly to the cubic eliminant
on the remaining high-chart sector X<0,Y<0.

The answer is no for a structural reason: that sector contains genuine
cubic-divisible branches.  Every Newton coefficient polynomial of the
negative-Y eliminant has mixed signs and positive roots.  Fixed-depth
endpoint controls nevertheless show no shared cubic--quartic point at
the sampled depths, so a future proof should apply discrete positivity
to a branch-reduced exit invariant, not to the cubic eliminant itself.
"""

from __future__ import annotations

import importlib.util
from hashlib import sha256
from pathlib import Path

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_thm2879():
    dependency = Path(__file__).with_name(
        "gmc_all_shift_cubic_null_endpoint_holonomy_thm2879.py"
    )
    dependency_bytes = dependency.read_bytes().replace(b"\r\n", b"\n")
    require(
        sha256(dependency_bytes).hexdigest()
        == "44012d84c88a22f246ef350f7f9a364116ac1fc839347361dee64c0a9c4a6e27",
        "THM-2879 exact dependency hash changed",
    )
    spec = importlib.util.spec_from_file_location("thm2879_exact", dependency)
    require(spec is not None and spec.loader is not None, "cannot load THM-2879")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    source = load_thm2879()
    n, x, y = source.n, source.x, source.y
    a = sp.symbols("a")

    invariant_numerators = tuple(
        sp.together(invariant).as_numer_denom()[0]
        for invariant in (source.cubic_one, source.cubic_two)
    )
    subresultants = sp.subresultants(*invariant_numerators, x)
    require(
        [sp.degree(item, x) for item in subresultants] == [4, 3, 2, 1, 0],
        "cubic subresultant profile changed",
    )
    eliminant = next(
        factor
        for factor, _ in sp.factor_list(subresultants[-1])[1]
        if sp.degree(factor, y) == 15
    )
    if sp.Poly(eliminant, y).LC().subs(n, 1) < 0:
        eliminant = -eliminant

    linear = sp.Poly(subresultants[-2], x)
    content = sp.gcd(linear.nth(1), linear.nth(0))
    selector_a = sp.cancel(linear.nth(1) / content)
    selector_n = sp.cancel(-linear.nth(0) / content)
    if selector_a.subs({n: 1, y: 1}) < 0:
        selector_a = -selector_a
        selector_n = -selector_n

    # Gregory--Newton expansion of the negative-Y cubic eliminant.
    negative_y_eliminant = sp.expand(eliminant.subs(y, -a))
    depth_degree = sp.degree(negative_y_eliminant, n)
    current = negative_y_eliminant
    newton_profiles: list[tuple[int, int, int, int, int]] = []
    for order in range(depth_degree + 1):
        coefficient = sp.Poly(current.subs(n, 1), a)
        positive = sum(1 for value in coefficient.all_coeffs() if value > 0)
        negative = sum(1 for value in coefficient.all_coeffs() if value < 0)
        positive_roots = int(coefficient.count_roots(0, sp.oo))
        newton_profiles.append(
            (order, len(coefficient.terms()), positive, negative, positive_roots)
        )
        current = sp.expand(current.subs(n, n + 1) - current)
    require(
        depth_degree == 21
        and len(newton_profiles) == 22
        and all(
            terms == 16 and positive > 0 and negative > 0 and roots >= 2
            for _, terms, positive, negative, roots in newton_profiles
        ),
        "negative-Y Newton obstruction profile changed",
    )

    # Fixed-depth branch and endpoint controls.
    endpoint_numerator = sp.together(source.endpoint_holonomy).as_numer_denom()[0]
    endpoint_in_x = sp.Poly(endpoint_numerator, x)
    endpoint_degree = endpoint_in_x.degree()
    require(endpoint_degree == 5, "endpoint x-degree changed")

    depth_profiles: list[tuple[int, int, int, tuple[int, ...], int]] = []
    for depth in (0, 1, 2, 8, 64):
        fixed_eliminant = sp.Poly(eliminant.subs(n, depth), y)
        fixed_a = sp.Poly(selector_a.subs(n, depth), y)
        fixed_n = sp.Poly(selector_n.subs(n, depth), y)
        negative_roots = [
            sp.CRootOf(fixed_eliminant.as_expr(), index)
            for index in range(fixed_eliminant.degree())
        ]
        negative_roots = [
            root for root in negative_roots if root.is_real and root < 0
        ]

        cleared_endpoint = sum(
            (
                sp.Poly(endpoint_in_x.nth(index).subs(n, depth), y)
                * fixed_n**index
                * fixed_a ** (endpoint_degree - index)
                for index in range(endpoint_degree + 1)
            ),
            sp.Poly(0, y),
        )
        endpoint_remainder = cleared_endpoint.rem(fixed_eliminant)
        endpoint_gcd_degree = sp.gcd(
            fixed_eliminant, endpoint_remainder
        ).degree()
        require(
            endpoint_gcd_degree == 0,
            f"depth {depth}: sampled endpoint acquired a common factor",
        )

        cone_avoiding_signs: list[int] = []
        for root in negative_roots:
            sign_a = int(sp.sign(fixed_a.as_expr().subs(y, root)))
            sign_n = int(sp.sign(fixed_n.as_expr().subs(y, root)))
            sign_x = sign_a * sign_n
            if sign_x < 0:
                sign_cleared = int(
                    sp.sign(endpoint_remainder.as_expr().subs(y, root))
                )
                # cleared_endpoint=A^5*endpoint_numerator.
                cone_avoiding_signs.append(sign_cleared * sign_a)

        expected_count = 2
        require(
            len(cone_avoiding_signs) == expected_count
            and all(sign != 0 for sign in cone_avoiding_signs),
            f"depth {depth}: cone-avoiding branch profile changed",
        )
        depth_profiles.append(
            (
                depth,
                len(negative_roots),
                len(cone_avoiding_signs),
                tuple(cone_avoiding_signs),
                endpoint_gcd_degree,
            )
        )

    root_histogram: dict[int, int] = {}
    for _, _, _, _, roots in newton_profiles:
        root_histogram[roots] = root_histogram.get(roots, 0) + 1

    print("GMC CONE-AVOIDING GREGORY-NEWTON PROBE")
    print("status=FINITE-EXACT SCOUT / NO CLOSURE")
    print("chart=consecutive high coordinates X<0,Y<0")
    print(f"cubic_eliminant_degrees=n:{depth_degree},Y:15")
    print("newton_base=n=1;variable=a=-Y>0")
    print(f"newton_coefficients={len(newton_profiles)}")
    print("newton_all_mixed_sign=true")
    print(
        "newton_positive_root_histogram="
        + ",".join(f"{roots}:{count}" for roots, count in sorted(root_histogram.items()))
    )
    print(
        "newton_sign_profiles="
        + ",".join(
            f"k{order}:{positive}+/{negative}-"
            for order, _, positive, negative, _ in newton_profiles
        )
    )
    print(
        "fixed_depth_profiles="
        + ";".join(
            f"n{depth}:negativeY={negative_count}:Xnegative={cone_count}:"
            f"endpoint_signs={endpoint_signs}:gcd_degree={gcd_degree}"
            for (
                depth,
                negative_count,
                cone_count,
                endpoint_signs,
                gcd_degree,
            ) in depth_profiles
        )
    )
    print(
        "assessment=integer depth n is genuine on each fixed-gap translation "
        "slice, but direct Newton positivity of the cubic eliminant is impossible "
        "because two cone-avoiding cubic branches are real"
    )
    print(
        "next_target=reduce endpoint holonomy modulo the cubic branch first, "
        "then seek a Newton/Bernstein certificate with the branch-order sidecar"
    )
    print(
        "scope=no all-n endpoint exit; no general mixed four-slot, GMC2, DvdK, "
        "or JC consequence"
    )


if __name__ == "__main__":
    main()
