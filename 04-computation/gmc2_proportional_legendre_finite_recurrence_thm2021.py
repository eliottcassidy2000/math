#!/usr/bin/env python3
"""Exact audit for THM-2021 and HYP-8771.

The script verifies the primitive-return polynomial, direct Gaussian moments on
two realizable proportional slices, the Legendre recurrence, and exact pairwise
gcd data for several higher-charge return-polynomial families.
"""

from __future__ import annotations

from collections import Counter
from math import factorial, gcd

import sympy as sp


K = sp.symbols("kappa")
S = sp.symbols("s")
Z, W = sp.symbols("Z W")


def return_poly(m: int, p0: int, q0: int) -> sp.Poly:
    """A_m(kappa) for coprime primitive charges p0,-q0."""
    assert gcd(p0, q0) == 1
    r = p0 + q0
    value = 0
    for k in range(m // r + 1):
        value += (
            sp.Integer(factorial(m))
            / (
                factorial(q0 * k)
                * factorial(p0 * k)
                * factorial(m - r * k)
            )
            * K**k
        )
    return sp.Poly(value, K, domain=sp.QQ)


def radial_L(poly: sp.Expr) -> sp.Expr:
    expanded = sp.Poly(sp.expand(poly), S)
    return sp.expand(
        sum(coeff * factorial(power[0]) for power, coeff in expanded.terms())
    )


def gaussian_L(poly: sp.Expr) -> sp.Expr:
    expanded = sp.Poly(sp.expand(poly), Z, W)
    value = 0
    for (a, b), coeff in expanded.terms():
        if a == b:
            value += coeff * factorial(a)
    return sp.expand(value)


def channel_moment(
    m: int, p0: int, q0: int, b: sp.Expr, h: sp.Expr
) -> sp.Expr:
    r = p0 + q0
    value = 0
    for k in range(m // r + 1):
        value += (
            sp.Integer(factorial(m))
            / (
                factorial(q0 * k)
                * factorial(p0 * k)
                * factorial(m - r * k)
            )
            * radial_L(h**k * b ** (m - r * k))
        )
    return sp.expand(value)


def verify_direct_slices() -> list[str]:
    reports: list[str] = []

    # p=q=1, b=s*d, a=d, c=kappa*s*d, hence h=s*a*c=kappa*b^2.
    d = 1 - 2 * S + 2 * S**2
    b = S * d
    kappa = sp.Rational(-2, 3)
    a = d
    c = kappa * S * d
    P = Z * a.subs(S, Z * W) + b.subs(S, Z * W) + W * c.subs(S, Z * W)
    for m in range(1, 8):
        direct = gaussian_L(P**m)
        channel = channel_moment(m, 1, 1, b, kappa * b**2)
        factored = return_poly(m, 1, 1).as_expr().subs(K, kappa) * radial_L(b**m)
        assert sp.expand(direct - channel) == 0
        assert sp.expand(channel - factored) == 0
    reports.append("p=q=1 realizable polynomial slice: direct=channel=factored for m=1..7")

    # Primitive charges +1,-2: a=s^2, b=s^2, c=kappa gives
    # h=s^2*a^2*c=kappa*b^3.
    kappa = sp.Rational(3, 5)
    a = S**2
    b = S**2
    c = kappa
    P = Z * a.subs(S, Z * W) + b.subs(S, Z * W) + W**2 * c
    for m in range(1, 8):
        direct = gaussian_L(P**m)
        channel = channel_moment(m, 1, 2, b, kappa * b**3)
        factored = return_poly(m, 1, 2).as_expr().subs(K, kappa) * radial_L(b**m)
        assert sp.expand(direct - channel) == 0
        assert sp.expand(channel - factored) == 0
    reports.append("charges (+1,-2) realizable slice: direct=channel=factored for m=1..7")

    # MISTAKE-212 witness: all channels have the same radial word b^m and hence
    # the same factorial degree, but noncancellation does not require a source.
    b = S * (S - 2)
    a = S - 2
    c = S * (S - 2)
    P = Z * a.subs(S, Z * W) + b.subs(S, Z * W) + W * c.subs(S, Z * W)
    moment_1 = gaussian_L(P)
    moment_2 = gaussian_L(P**2)
    assert moment_1 == 0
    assert moment_2 == 24
    assert sp.expand(S * a * c - b**2) == 0
    reports.append("MISTAKE-212 tied-channel witness: h=b^2, M_1=0, M_2=24")
    return reports


def verify_legendre() -> str:
    polys = [return_poly(m, 1, 1).as_expr() for m in range(0, 31)]
    assert polys[0] == 1 and polys[1] == 1
    for m in range(1, 30):
        residual = (m + 1) * polys[m + 1] - (2 * m + 1) * polys[m]
        residual += m * (1 - 4 * K) * polys[m - 1]
        assert sp.expand(residual) == 0

    t = sp.symbols("t")
    truncated = sum(polys[m] * t**m for m in range(21))
    target = sp.series((1 - 2 * t + (1 - 4 * K) * t**2) ** sp.Rational(-1, 2), t, 0, 21).removeO()
    assert sp.expand(truncated - target) == 0
    return "Legendre recurrence and OGF verified exactly through m=30 / t^20"


def verify_computed_root_localization() -> str:
    checked = 0
    for m in range(2, 31):
        poly = return_poly(m, 1, 1)
        for root in sp.nroots(poly.as_expr(), n=40, maxsteps=200):
            value = complex(root)
            assert abs(value.imag) < 1e-28
            assert value.real < 0
            checked += 1
    return f"all {checked} computed roots of S_m(kappa), m=2..30, are negative real"


def strongly_connected_component_sizes(n: int, edges: set[tuple[int, int]]) -> list[int]:
    adjacency = [[] for _ in range(n)]
    reverse = [[] for _ in range(n)]
    for u, v in edges:
        adjacency[u].append(v)
        reverse[v].append(u)

    seen: set[int] = set()
    order: list[int] = []

    def visit(u: int) -> None:
        seen.add(u)
        for v in adjacency[u]:
            if v not in seen:
                visit(v)
        order.append(u)

    for u in range(n):
        if u not in seen:
            visit(u)

    seen.clear()
    sizes: list[int] = []

    def collect(u: int) -> int:
        seen.add(u)
        size = 1
        for v in reverse[u]:
            if v not in seen:
                size += collect(v)
        return size

    for u in reversed(order):
        if u not in seen:
            sizes.append(collect(u))
    return sorted(sizes, reverse=True)


def tournament_audit(p0: int, q0: int, max_m: int) -> str:
    r = p0 + q0
    levels = list(range(r, max_m + 1))
    polys = {m: return_poly(m, p0, q0) for m in levels}
    flips: list[tuple[int, int, sp.Expr]] = []
    edges: set[tuple[int, int]] = set()

    # Gauge: chronological edge i->j, flipped exactly when the pairwise gcd has
    # positive degree. Thus an edge flip is an exact shared-root event.
    for i, m in enumerate(levels):
        for j in range(i + 1, len(levels)):
            n = levels[j]
            common = sp.gcd(polys[m], polys[n])
            if common.degree() > 0:
                edges.add((j, i))
                flips.append((m, n, sp.factor(common.as_expr())))
            else:
                edges.add((i, j))

    n_vertices = len(levels)
    scores = [0] * n_vertices
    for u, _ in edges:
        scores[u] += 1
    score_hist = sorted(Counter(scores).items())
    expected_transitive_hist = [(j, 1) for j in range(n_vertices)]
    score_report = (
        f"0..{n_vertices-1}(each once)"
        if score_hist == expected_transitive_hist
        else str(score_hist)
    )
    triangles = 0
    for i in range(n_vertices):
        for j in range(i + 1, n_vertices):
            for k in range(j + 1, n_vertices):
                cycle = (
                    ((i, j) in edges and (j, k) in edges and (k, i) in edges)
                    or ((j, i) in edges and (k, j) in edges and (i, k) in edges)
                )
                triangles += int(cycle)
    scc = strongly_connected_component_sizes(n_vertices, edges)
    scc_report = f"1^{n_vertices}" if scc == [1] * n_vertices else str(scc)
    tie_path_valid = all((i, i + 1) in edges for i in range(n_vertices - 1))
    hp_count = 1 if not flips else "not-counted-after-flips"
    assert not flips, flips
    assert tie_path_valid

    return (
        f"(+{p0},-{q0}) levels={levels[0]}..{levels[-1]}: "
        f"pairs={n_vertices*(n_vertices-1)//2}, edge_flips={len(flips)}, "
        f"score_hist={score_report}, directed_3cycles={triangles}, "
        f"SCC_sizes={scc_report}, tie_path=chronological(valid), Hamiltonian_paths={hp_count}"
    )


def main() -> None:
    print("THM-2021 / HYP-8771 EXACT AUDIT")
    for report in verify_direct_slices():
        print("PASS", report)
    print("PASS", verify_legendre())
    print("PASS", verify_computed_root_localization(), "(Legendre localization used in proof)")
    for args in [(1, 1, 36), (1, 2, 30), (1, 3, 26), (2, 3, 22)]:
        print("TOURNAMENT", tournament_audit(*args))
    print("ALTERNATE VERTICES considered: return channels, individual roots, moment levels, proof obligations")
    print("CHOSEN quotient: moment levels; observable=degree gcd(A_m,A_n); chronological gauge flips on shared roots")
    print("PRESERVES: exact pairwise root recurrence; DESTROYS: root identity, multiplicity, and analytic location")
    print("CHALLENGED ASSUMPTION / MISTAKE-213: finite zero recurrence is not needed for NC2; no-consecutive-zero already gives unbounded toral nonzeros, which meet EMP's cofinite radial tail")
    print("NOTE external 2026 Legendre finiteness is not proved by this computation")


if __name__ == "__main__":
    main()
