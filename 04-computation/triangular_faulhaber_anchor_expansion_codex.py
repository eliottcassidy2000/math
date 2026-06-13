#!/usr/bin/env python3
"""Faulhaber anchor expansion for triangular power-balance towers.

codex-2026-06-13

For row n and anchor a, compare

    sum_{j=0..n} (a+j)^p    and    sum_{j=1..n} (a+n+j)^p.

Writing c=a+n puts the comparison around the midpoint:

    D_p(c,n) = sum_{i=0..n} (c-i)^p - sum_{i=1..n} (c+i)^p.

The binomial expansion kills all even powers and leaves only odd Faulhaber
moments S_r(n)=sum_{i=1..n} i^r:

    D_p(c,n) = c^p - 2 * sum_{r odd} binom(p,r) c^(p-r) S_r(n).

This script verifies that identity, derives/checks the first expansion terms
for the real root, records the square-pyramidal cuboid carrier at p=2, and
adds Tournament Analysis over the possible proof-transfer carriers.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import comb

import mpmath as mp
import sympy as sp


def T(n: int) -> int:
    return n * (n + 1) // 2


def U(n: int) -> int:
    return n * (n + 1)


def power_sum(n: int, r: int) -> int:
    return sum(i**r for i in range(1, n + 1))


def balance_defect(p: int, c: mp.mpf, n: int) -> mp.mpf:
    return sum((c - i) ** p for i in range(n + 1)) - sum(
        (c + i) ** p for i in range(1, n + 1)
    )


def exact_balance_defect(p: int, c: int, n: int) -> int:
    return sum((c - i) ** p for i in range(n + 1)) - sum(
        (c + i) ** p for i in range(1, n + 1)
    )


def balance_derivative(p: int, c: mp.mpf, n: int) -> mp.mpf:
    return sum(p * (c - i) ** (p - 1) for i in range(n + 1)) - sum(
        p * (c + i) ** (p - 1) for i in range(1, n + 1)
    )


def odd_moment_defect(p: int, c: int, n: int) -> int:
    return c**p - 2 * sum(
        comb(p, r) * c ** (p - r) * power_sum(n, r)
        for r in range(1, p + 1, 2)
    )


def alpha(p: int) -> Fraction:
    return Fraction((p - 1) * (p - 2), 12 * p)


def beta(p: int) -> Fraction:
    return -Fraction((p - 1) * (p - 2) * (2 * p * p - 4 * p - 1), 180 * p**3)


def gamma(p: int) -> Fraction:
    return Fraction(
        (p - 1)
        * (p - 2)
        * (p + 2)
        * (40 * p**3 - 128 * p**2 + 53 * p + 11),
        30240 * p**5,
    )


def frac_to_mpf(x: Fraction) -> mp.mpf:
    return mp.mpf(x.numerator) / mp.mpf(x.denominator)


def center_expansion(p: int, n: int, order: int = 3) -> mp.mpf:
    u = mp.mpf(U(n))
    c = mp.mpf(p) * u
    if order >= 1:
        c += frac_to_mpf(alpha(p))
    if order >= 2:
        c += frac_to_mpf(beta(p)) / u
    if order >= 3:
        c += frac_to_mpf(gamma(p)) / (u * u)
    return c


def newton_center_root(p: int, n: int) -> mp.mpf:
    c = center_expansion(p, n, order=3)
    for _ in range(40):
        step = balance_defect(p, c, n) / balance_derivative(p, c, n)
        c -= step
        if abs(step) < mp.mpf("1e-60") * max(1, abs(c)):
            break
    return c


def symbolic_odd_moment_identity(max_p: int = 10) -> list[tuple[int, bool]]:
    c, n, i = sp.symbols("c n i", integer=True, positive=True)
    rows = []
    for p in range(1, max_p + 1):
        direct = sp.summation((c - i) ** p, (i, 0, n)) - sp.summation(
            (c + i) ** p, (i, 1, n)
        )
        odd = c**p - 2 * sum(
            sp.binomial(p, r)
            * c ** (p - r)
            * sp.summation(i**r, (i, 1, n))
            for r in range(1, p + 1, 2)
        )
        rows.append((p, sp.simplify(direct - odd) == 0))
    return rows


def reduced_series_residuals(max_p: int = 12) -> list[tuple[int, list[sp.Expr]]]:
    """Check the u-expansion coefficients through gamma_p for fixed p.

    Put q=1/u and x=c/u.  The relevant reduced equation through q^3 is

        x^p - p x^(p-1)
        - C(p,3) x^(p-3) q / 2
        - C(p,5) x^(p-5) q^2 / 3
        + C(p,5) x^(p-5) q^3 / 6
        - C(p,7) x^(p-7) q^3 / 4

    using S1=u/2, S3=u^2/4, S5=u^3/6-u^2/12, and
    S7=u^4/8-u^3/6+u^2/12.  Coefficients q^0..q^3 vanish.
    """
    q = sp.symbols("q")
    rows = []
    for p in range(1, max_p + 1):
        a = sp.Rational(alpha(p).numerator, alpha(p).denominator)
        b = sp.Rational(beta(p).numerator, beta(p).denominator)
        g = sp.Rational(gamma(p).numerator, gamma(p).denominator)
        x = sp.Integer(p) + a * q + b * q**2 + g * q**3
        expr = x**p - p * x ** (p - 1)
        if p >= 3:
            expr -= sp.binomial(p, 3) * x ** (p - 3) * q / 2
        if p >= 5:
            expr -= sp.binomial(p, 5) * x ** (p - 5) * q**2 / 3
            expr += sp.binomial(p, 5) * x ** (p - 5) * q**3 / 6
        if p >= 7:
            expr -= sp.binomial(p, 7) * x ** (p - 7) * q**3 / 4
        series = sp.series(expr, q, 0, 4).removeO().expand()
        coeffs = [sp.simplify(series.coeff(q, k)) for k in range(4)]
        rows.append((p, coeffs))
    return rows


def strongly_connected_components(vertices: list[str], edges: dict[tuple[str, str], str]) -> list[list[str]]:
    adjacency = {v: [] for v in vertices}
    for (a, b), winner in edges.items():
        loser = b if winner == a else a
        adjacency[winner].append(loser)

    index = 0
    stack: list[str] = []
    on_stack: set[str] = set()
    indices: dict[str, int] = {}
    low: dict[str, int] = {}
    out: list[list[str]] = []

    def visit(v: str) -> None:
        nonlocal index
        indices[v] = index
        low[v] = index
        index += 1
        stack.append(v)
        on_stack.add(v)
        for w in adjacency[v]:
            if w not in indices:
                visit(w)
                low[v] = min(low[v], low[w])
            elif w in on_stack:
                low[v] = min(low[v], indices[w])
        if low[v] == indices[v]:
            comp = []
            while True:
                w = stack.pop()
                on_stack.remove(w)
                comp.append(w)
                if w == v:
                    break
            out.append(comp)

    for v in vertices:
        if v not in indices:
            visit(v)
    return out


def count_hamiltonian_paths(vertices: list[str], edges: dict[tuple[str, str], str]) -> int:
    idx = {v: i for i, v in enumerate(vertices)}
    n = len(vertices)
    beats = [[False] * n for _ in range(n)]
    for (a, b), winner in edges.items():
        loser = b if winner == a else a
        beats[idx[winner]][idx[loser]] = True
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            value = dp[mask][last]
            if not value:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if beats[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += value
    return sum(dp[(1 << n) - 1])


def tournament_analysis() -> dict[str, object]:
    tie_path = [
        "odd_moment_projection",
        "u_anchor_expansion",
        "p1_p2_exact_towers",
        "square_pyramid_cuboid",
        "simplex_cuboid_packing",
        "bernoulli_denominator_address",
        "lrc14_moment_resource",
        "support_lift_transfer",
    ]
    scores = {
        # exactness, proof_readiness, geometry, lrc_transfer, computability, novelty
        "odd_moment_projection": (5, 5, 3, 4, 5, 4),
        "u_anchor_expansion": (5, 4, 3, 4, 5, 5),
        "p1_p2_exact_towers": (5, 5, 4, 3, 5, 3),
        "square_pyramid_cuboid": (5, 4, 5, 3, 5, 4),
        "simplex_cuboid_packing": (3, 3, 5, 2, 4, 4),
        "bernoulli_denominator_address": (4, 4, 3, 5, 4, 5),
        "lrc14_moment_resource": (3, 3, 3, 5, 3, 5),
        "support_lift_transfer": (3, 4, 4, 5, 3, 5),
    }
    edges: dict[tuple[str, str], str] = {}
    vertices = list(tie_path)
    flips_vs_exactness = 0
    for a, b in combinations(vertices, 2):
        a_score = scores[a]
        b_score = scores[b]
        votes_a = sum(1 for x, y in zip(a_score, b_score) if x > y)
        votes_b = sum(1 for x, y in zip(a_score, b_score) if y > x)
        if votes_a > votes_b:
            winner = a
        elif votes_b > votes_a:
            winner = b
        else:
            winner = a if tie_path.index(a) < tie_path.index(b) else b
        edges[(a, b)] = winner
        exact_winner = a if a_score[0] >= b_score[0] else b
        if winner != exact_winner:
            flips_vs_exactness += 1

    outdegrees = Counter({v: 0 for v in vertices})
    for winner in edges.values():
        outdegrees[winner] += 1
    cycles = []
    for a, b, c in combinations(vertices, 3):
        ab = edges[tuple(sorted((a, b), key=vertices.index))]
        ac = edges[tuple(sorted((a, c), key=vertices.index))]
        bc = edges[tuple(sorted((b, c), key=vertices.index))]
        if (ab, bc, ac) in ((a, b, c), (b, c, a)):
            cycles.append((a, b, c))
    sccs = strongly_connected_components(vertices, edges)
    return {
        "criteria": (
            "exactness",
            "proof_readiness",
            "geometry",
            "lrc_transfer",
            "computability",
            "novelty",
        ),
        "tie_path": tie_path,
        "scores": scores,
        "ranking": sorted(vertices, key=lambda v: (-outdegrees[v], tie_path.index(v))),
        "outdegrees": dict(outdegrees),
        "score_hist": dict(sorted(Counter(outdegrees.values()).items())),
        "directed_3cycles": cycles,
        "scc_sizes": sorted((len(c) for c in sccs), reverse=True),
        "hamiltonian_paths": count_hamiltonian_paths(vertices, edges),
        "edge_flips_vs_exactness": flips_vs_exactness,
    }


def print_identity_section() -> None:
    print("[1] Odd Faulhaber moment identity")
    checks = symbolic_odd_moment_identity(10)
    print(f"    Symbolic identity checks p=1..10: {checks}")
    print("    Exact formula:")
    print("      D_p(c,n)=c^p - 2*sum_{r odd} binom(p,r)c^(p-r)S_r(n)")
    print("      S_r(n)=sum_{i=1..n} i^r")
    print("    Sample integer checks:")
    for p, n, c in [(3, 4, 60), (4, 5, 121), (7, 6, 300)]:
        direct = exact_balance_defect(p, c, n)
        odd = odd_moment_defect(p, c, n)
        print(f"      p={p}, n={n}, c={c}: direct={direct}, odd_moment={odd}")
    print()


def print_expansion_section() -> None:
    print("[2] Anchor expansion in u=n(n+1)")
    print("    c=a+n, u=n(n+1), c_p(n)=p*u + alpha_p + beta_p/u + gamma_p/u^2 + ...")
    print("    a_p(n)=c_p(n)-n.")
    print(f"    alpha_p = {sp.factor((sp.Symbol('p') - 1) * (sp.Symbol('p') - 2) / (12 * sp.Symbol('p')))}")
    p = sp.Symbol("p")
    beta_expr = -((p - 1) * (p - 2) * (2 * p * p - 4 * p - 1)) / (180 * p**3)
    gamma_expr = (
        (p - 1)
        * (p - 2)
        * (p + 2)
        * (40 * p**3 - 128 * p**2 + 53 * p + 11)
        / (30240 * p**5)
    )
    print(f"    beta_p  = {sp.factor(beta_expr)}")
    print(f"    gamma_p = {sp.factor(gamma_expr)}")
    residuals = reduced_series_residuals(12)
    ok = all(all(x == 0 for x in coeffs) for _, coeffs in residuals)
    print(f"    Reduced series coefficients q^0..q^3 vanish for p=1..12: {ok}")
    print("    First coefficients by p:")
    for p0 in range(1, 9):
        print(
            f"      p={p0}: alpha={alpha(p0)}, beta={beta(p0)}, gamma={gamma(p0)}"
        )
    print()
    print("    Exact p=1 and p=2 anchors recovered:")
    for n in range(1, 7):
        c1 = U(n)
        c2 = 2 * U(n)
        print(
            f"      n={n}: p=1 a=c-n={c1-n}=n^2; "
            f"p=2 a=c-n={c2-n}=2n^2+n"
        )
    print()


def print_numerical_section() -> None:
    print("[3] Numerical root check")
    mp.mp.dps = 80
    for p in (3, 4, 5, 8):
        print(f"    p={p}: target gamma={gamma(p)}")
        for n in (8, 32, 128):
            u = mp.mpf(U(n))
            root = newton_center_root(p, n)
            approx2 = center_expansion(p, n, order=2)
            approx3 = center_expansion(p, n, order=3)
            scaled2 = (root - approx2) * u**2
            scaled3 = (root - approx3) * u**3
            print(
                "      n={:3d}: (c-c2)u^2={} ; (c-c3)u^3={}".format(
                    n,
                    mp.nstr(scaled2, 18),
                    mp.nstr(scaled3, 18),
                )
            )
    print()


def print_square_pyramid_section() -> None:
    print("[4] Square-pyramidal cuboid carrier")
    print("    P2(n)=1^2+...+n^2=n(n+1)(2n+1)/6.")
    print("    Six square pyramids tile the cuboid n*(n+1)*(2n+1).")
    print("    HYP-2454's first tower common sum is S1=3*P2, so the cuboid is 2*S1.")
    for n in range(1, 8):
        p2 = power_sum(n, 2)
        cuboid = n * (n + 1) * (2 * n + 1)
        s1 = n * (n + 1) * (2 * n + 1) // 2
        print(
            f"      n={n}: P2={p2}, 6P2={6*p2}, cuboid={cuboid}, "
            f"S1={s1}, 2S1={2*s1}"
        )
    print()
    print("    Packing reading:")
    print("      p=1 is the square-shell carrier; p=2 is the square-pyramidal cuboid carrier.")
    print("      Higher p keeps the same u=n(n+1) carrier, but Bernoulli denominators")
    print("      move the real anchor off the integer tower by alpha_p, beta_p/u, ...")
    print()


def print_tournament_section() -> None:
    print("[5] Tournament Analysis")
    ta = tournament_analysis()
    print(f"    Criteria: {ta['criteria']}")
    print(f"    Tie Hamiltonian path: {ta['tie_path']}")
    print("    Vertex scores:")
    for name in ta["tie_path"]:
        print(f"      {name}: {ta['scores'][name]}")
    print(f"    Ranking: {ta['ranking']}")
    print(f"    Outdegrees: {ta['outdegrees']}")
    print(f"    score_hist={ta['score_hist']}")
    print(f"    directed_3cycles={len(ta['directed_3cycles'])}: {ta['directed_3cycles'][:8]}")
    print(f"    scc_sizes={ta['scc_sizes']}")
    print(f"    hamiltonian_paths={ta['hamiltonian_paths']}")
    print(f"    edge_flips_vs_exactness={ta['edge_flips_vs_exactness']}")
    print()
    print("    Assumption challenge:")
    print("      Tournament vertices need not be row entries or runners.  Considered vertices")
    print("      were odd moments, u-carriers, Bernoulli denominators, exact towers, cuboid")
    print("      packings, simplex-packing analogues, LRC proof obligations, and support lifts.")
    print("      The chosen quotient preserves the balance equation, odd-moment ledger,")
    print("      exact p=1/p=2 tower recovery, and transfer channels.  It destroys individual")
    print("      endpoint residues and any actual LRC/unit-distance/code witness.")
    print("      Challenged assumption: the scalar anchor is not the proof object; the proof")
    print("      object is the retained odd-moment/support address that explains when the")
    print("      anchor stays integral and when it drifts into Bernoulli fractions.")
    print()


def main() -> None:
    print("Triangular Faulhaber anchor expansion")
    print("HYP-2457 / T801 / OPEN-Q-079 candidate")
    print()
    print_identity_section()
    print_expansion_section()
    print_numerical_section()
    print_square_pyramid_section()
    print_tournament_section()
    print("[6] Next proof tasks")
    print("    1. Prove the expansion with a uniform remainder in u for fixed p.")
    print("    2. Prove the HYP-2454 integer bracket for p>=3 using the odd-moment formula.")
    print("    3. Lift the odd-moment ledger to LRC14: odd walls, owner/carry fibers, and support.")
    print("    4. Compare simplex-in-cuboid packing constants with higher Faulhaber carriers.")
    print("    5. Search for a finite-dimensional residue/tournament quotient that detects")
    print("       Bernoulli denominator drift before the full polynomial is expanded.")


if __name__ == "__main__":
    main()
