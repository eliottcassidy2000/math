#!/usr/bin/env python3
"""Irreducibility, prime values, and tournament witness prisms.

codex-2026-06-12

External seeds:
  * arXiv:2411.18366, Prime numbers and factorization of polynomials.
  * arXiv:2410.15880, Integer Polynomial Factorization by Recombination
    of Real Factors.

This script keeps four related directions typed:

1. Bunyakovsky/Buniakowski: irreducible + fixed divisor 1 should produce
   infinitely many prime values.
2. Murty/Girstmair/Singh: a sufficiently large value with few prime factors
   bounds the number of irreducible factors of the polynomial.
3. Cohn: a prime written in base b produces an irreducible digit polynomial.
4. Real-factor recombination: factor over a looser field, then recover integer
   factors by subset-sum recombination of root-trace atoms.

The repo transfer is that a tournament should not be built on raw integers
alone.  Vertices are proof obligations / factor atoms / witness slots, and
edges record which side channel preserves irreducibility information.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import reduce
from itertools import combinations, permutations
from math import gcd, isqrt

import sympy as sp


x = sp.Symbol("x")


def is_prime(n: int) -> bool:
    n = abs(int(n))
    if n < 2:
        return False
    if n in (2, 3):
        return True
    if n % 2 == 0:
        return False
    r = isqrt(n)
    d = 3
    while d <= r:
        if n % d == 0:
            return False
        d += 2
    return True


def factor_count_with_multiplicity(n: int) -> int:
    n = abs(int(n))
    if n <= 1:
        return 0
    return sum(sp.factorint(n).values())


def fixed_divisor(poly: sp.Poly) -> int:
    """Fixed divisor via deg+1 consecutive values."""
    deg = poly.degree()
    vals = [abs(int(poly.eval(k))) for k in range(deg + 1)]
    return reduce(gcd, vals)


def primitive_part(poly: sp.Poly) -> sp.Poly:
    _, prim = sp.primitive(poly.as_expr(), expand=True)
    return sp.Poly(prim, x, domain=sp.ZZ)


def factor_count_over_Z(poly: sp.Poly) -> int:
    coeff, factors = sp.factor_list(poly.as_expr(), x)
    return sum(exp for _, exp in factors)


def is_irreducible_over_Z(poly: sp.Poly) -> bool:
    return factor_count_over_Z(poly) == 1


def height_H(poly: sp.Poly) -> float:
    coeffs = [int(c) for c in poly.all_coeffs()]  # high to low
    lead = abs(coeffs[0])
    return max(abs(c) / lead for c in coeffs[1:]) if coeffs[1:] else 0.0


def prime_run(poly: sp.Poly, start: int = 0, limit: int = 500) -> int:
    t = start
    while t < limit and is_prime(poly.eval(t)):
        t += 1
    return t - start


@dataclass(frozen=True)
class PolyRow:
    name: str
    expr: sp.Expr
    reading: str


POLY_ROWS = [
    PolyRow("Euler41", x * x - x + 41, "irreducible/admissible; 40-value Euler run at x=0..39"),
    PolyRow("Rabinowitsch41", x * x + x + 41, "same horizon shifted; boundary p^2 at x=40"),
    PolyRow("BlockedEven", x * x + x + 2, "irreducible but fixed divisor 2; Bunyakovsky condition fails"),
    PolyRow("ReduciblePair", (x + 1) * (x + 2), "reducible but often has small prime-looking values"),
    PolyRow("Heegner17", x * x + x + 17, "Euler-Heegner prime horizon p=17"),
    PolyRow("NonHeegner7", x * x + x + 7, "irreducible/admissible but prime run fails immediately"),
]


def q_p(p: int) -> sp.Poly:
    return sp.Poly(x * x + x + p, x, domain=sp.ZZ)


def first_composite_x_for_qp(p: int) -> int:
    poly = q_p(p)
    t = 1
    while is_prime(poly.eval(t)):
        t += 1
    return t


def cohn_digit_poly(n: int, base: int = 10) -> sp.Poly:
    coeffs = []
    m = n
    if m == 0:
        coeffs.append(0)
    while m:
        coeffs.append(m % base)
        m //= base
    expr = sum(c * x**i for i, c in enumerate(coeffs))
    return sp.Poly(expr, x, domain=sp.ZZ)


@dataclass(frozen=True)
class SinghCase:
    name: str
    poly: sp.Poly
    m: int
    d: int
    reading: str


SINGH_CASES = [
    SinghCase(
        "paper_f2",
        sp.Poly(x**4 - x**2 - 2 * x - 1, x, domain=sp.ZZ),
        3,
        1,
        "Singh table: value 65=5*13 and exactly two irreducible factors",
    ),
    SinghCase(
        "paper_f3",
        sp.Poly(x**6 - x**2 - 2 * x - 1, x, domain=sp.ZZ),
        3,
        1,
        "Singh table: value 23*31 and exactly two irreducible factors",
    ),
    SinghCase(
        "quartic_a2",
        sp.Poly(x**4 + 10 * x**2 + 1, x, domain=sp.ZZ),
        12,
        1,
        "paper example: value has two prime factors; coefficient symmetry proves irreducible",
    ),
]


def unitary_divisor(d: int, n: int) -> bool:
    if d <= 0 or n % d != 0:
        return False
    return gcd(d, n // d) == 1


@dataclass(frozen=True)
class Atom:
    name: str
    degree: int
    trace: sp.Expr
    factor_label: str


RFR_ATOMS = [
    Atom("A_plus", 2, sp.sqrt(2), "x^4+1"),
    Atom("A_minus", 2, -sp.sqrt(2), "x^4+1"),
    Atom("B_unit", 2, 0, "x^2+1"),
]


def trace_integral(expr: sp.Expr) -> bool:
    return bool(sp.simplify(expr).is_integer)


def rfr_subset_table(atoms: list[Atom]) -> list[dict[str, object]]:
    rows = []
    n = len(atoms)
    group_sizes = Counter(a.factor_label for a in atoms)
    for mask in range(1, 1 << n):
        subset = [atoms[i] for i in range(n) if mask & (1 << i)]
        trace = sp.simplify(sum(a.trace for a in subset))
        labels = tuple(sorted({a.factor_label for a in subset}))
        counts = Counter(a.factor_label for a in subset)
        complete_groups = all(counts[label] == group_sizes[label] for label in counts)
        rows.append(
            {
                "subset": tuple(a.name for a in subset),
                "degree": sum(a.degree for a in subset),
                "trace": trace,
                "trace_integral": trace_integral(trace),
                "integer_factor": trace_integral(trace) and complete_groups,
                "labels": labels,
            }
        )
    return rows


@dataclass(frozen=True)
class Route:
    name: str
    domain: str
    scalar: str
    side_channel: str
    exactness: int
    bridge: int
    computable: int
    risk: int

    def score_tuple(self) -> tuple[int, int, int, int]:
        return (self.side_strength(), self.exactness, self.bridge, self.computable - self.risk)

    def side_strength(self) -> int:
        return {
            "prime_value_factorization": 5,
            "fixed_divisor_admissibility": 5,
            "digit_place_value": 4,
            "root_trace_subset_sum": 4,
            "heegner_class_group": 4,
            "raw_prime_run": 1,
        }[self.side_channel]


ROUTES = [
    Route(
        "bunyakovsky_admissible_irreducible",
        "prime values of polynomials",
        "irreducible polynomial",
        "fixed_divisor_admissibility",
        2,
        5,
        5,
        4,
    ),
    Route(
        "singh_value_factor_bound",
        "irreducibility criteria",
        "f(m) has few prime factors",
        "prime_value_factorization",
        5,
        5,
        5,
        1,
    ),
    Route(
        "cohn_digit_prime",
        "prime -> irreducible polynomial",
        "base-b digits of a prime",
        "digit_place_value",
        5,
        4,
        5,
        1,
    ),
    Route(
        "real_factor_recombination",
        "integer factorization algorithms",
        "factor over R into atoms",
        "root_trace_subset_sum",
        4,
        5,
        4,
        2,
    ),
    Route(
        "heegner_euler_horizon",
        "prime-generating quadratics",
        "long prime run",
        "heegner_class_group",
        4,
        5,
        5,
        2,
    ),
    Route(
        "raw_prime_run_length",
        "numerology",
        "how long values stay prime",
        "raw_prime_run",
        2,
        2,
        5,
        4,
    ),
]


def route_tournament(routes: list[Route]) -> dict[str, object]:
    wins: dict[int, set[int]] = {i: set() for i in range(len(routes))}
    flips_vs_scalar = 0
    for i, j in combinations(range(len(routes)), 2):
        a, b = routes[i], routes[j]
        if a.score_tuple() >= b.score_tuple():
            winner, loser = i, j
        else:
            winner, loser = j, i
        wins[winner].add(loser)
        scalar_pref = i if (a.exactness, a.computable) >= (b.exactness, b.computable) else j
        if scalar_pref != winner:
            flips_vs_scalar += 1

    cycles = 0
    for i, j, k in combinations(range(len(routes)), 3):
        if j in wins[i] and k in wins[j] and i in wins[k]:
            cycles += 1
        if k in wins[i] and j in wins[k] and i in wins[j]:
            cycles += 1

    ham = 0
    first_path = None
    for perm in permutations(range(len(routes))):
        if all(perm[t + 1] in wins[perm[t]] for t in range(len(perm) - 1)):
            ham += 1
            if first_path is None:
                first_path = tuple(routes[i].name for i in perm)

    return {
        "score_hist": dict(sorted(Counter(len(w) for w in wins.values()).items())),
        "directed_3cycles": cycles,
        "hamiltonian_paths": ham,
        "first_path": first_path,
        "flips_vs_scalar": flips_vs_scalar,
        "scores": sorted(((len(wins[i]), routes[i].name) for i in range(len(routes))), reverse=True),
    }


def main() -> None:
    print("=" * 78)
    print("Irreducibility-prime tournament prism")
    print("=" * 78)
    print("HYP-2447 / T791")
    print("External anchors: arXiv:2411.18366 and arXiv:2410.15880")
    print()

    print("[1] Bunyakovsky admissibility and prime-run rows")
    print("    fixed_divisor uses deg+1 consecutive values; run starts at x=0.")
    print("    name             deg  irreducible  fixed_divisor  prime_run  factors  reading")
    for row in POLY_ROWS:
        poly = primitive_part(sp.Poly(row.expr, x, domain=sp.ZZ))
        print(
            f"    {row.name:16s} {poly.degree():3d} {str(is_irreducible_over_Z(poly)):11s}"
            f" {fixed_divisor(poly):14d} {prime_run(poly):10d}"
            f" {factor_count_over_Z(poly):8d}  {row.reading}"
        )
    print()

    print("[2] Euler-Heegner boundary check")
    euler_primes = [2, 3, 5, 11, 17, 41]
    hits = [p for p in range(2, 80) if is_prime(p) and first_composite_x_for_qp(p) == p - 1]
    print(f"  primes p<80 whose Q_p(x)=x^2+x+p first positive composite is x=p-1: {hits}")
    print(f"  matches Euler set under 80? {hits == euler_primes}")
    print("  p  d=4p-1  first_bad  Q(first_bad)  positive_horizon")
    for p in euler_primes + [7, 13, 23, 29]:
        bad = first_composite_x_for_qp(p)
        print(f"  {p:2d} {4*p-1:7d} {bad:10d} {int(q_p(p).eval(bad)):12d} {bad-1:17d}")
    print()

    print("[3] Singh value-factor bound examples")
    print("    plus1 is the printed H_f+d+1 threshold; paper_ok follows the examples' H_f+d usage.")
    print("    name         H_f   m  d  unitary plus1 paper_ok |f(m)/d| factors actual_Z_factors")
    for case in SINGH_CASES:
        poly = primitive_part(case.poly)
        val = abs(int(poly.eval(case.m)))
        quotient = val // case.d
        nu = factor_count_with_multiplicity(quotient)
        actual = factor_count_over_Z(poly)
        H = height_H(poly)
        plus1 = case.m >= H + case.d + 1
        paper_threshold = case.m >= H + case.d
        unitary = unitary_divisor(case.d, val)
        paper_ok = actual <= nu and paper_threshold and unitary
        print(
            f"    {case.name:12s} {H:5.1f} {case.m:3d} {case.d:2d}"
            f" {str(unitary):7s} {str(plus1):5s} {str(paper_ok):8s}"
            f" {quotient:10d} {nu:7d} {actual:16d}"
        )
        print(f"      {case.reading}; factorization={sp.factor(poly.as_expr())}")
    print()

    print("[4] Cohn digit-prime criterion samples")
    print("    prime/base digits -> polynomial; irreducibility verified by sympy")
    for p in [2339, 1091, 6551, 7643, 10007]:
        poly = cohn_digit_poly(p, 10)
        print(
            f"    p={p:5d}  digit_poly={str(poly.as_expr()):24s}"
            f" irreducible={is_irreducible_over_Z(poly)}"
        )
    print()

    print("[5] Real-factor recombination as subset-sum over trace atoms")
    print("    Toy polynomial: (x^4+1)(x^2+1).  Over R, x^4+1 splits into")
    print("    two quadratic atoms with traces +/-sqrt(2); integer factors are")
    print("    recovered by recombining trace-compatible subsets.")
    for row in rfr_subset_table(RFR_ATOMS):
        print(
            f"    subset={row['subset']!s:28s} degree={row['degree']}"
            f" trace={str(row['trace']):8s} traceZ={row['trace_integral']}"
            f" integer_factor={row['integer_factor']} labels={row['labels']}"
        )
    print()

    print("[6] Tournament Analysis over proof-route vertices")
    analysis = route_tournament(ROUTES)
    print("  vertices:", ", ".join(r.name for r in ROUTES))
    print("  observable = (side_channel_strength, exactness, repo_bridge, computability-risk)")
    print("  score histogram:", analysis["score_hist"])
    print("  directed 3-cycles:", analysis["directed_3cycles"])
    print("  Hamiltonian paths:", analysis["hamiltonian_paths"])
    print("  first Hamiltonian path:", analysis["first_path"])
    print("  edge flips versus scalar-only ranking:", analysis["flips_vs_scalar"])
    print("  scores:")
    for score, name in analysis["scores"]:
        print(f"    {score}: {name}")
    print()

    print("[7] Repo transfer")
    print("  prime values -> irreducibility is a finite witness channel (Singh/Murty).")
    print("  irreducible + fixed divisor 1 -> prime values is the Bunyakovsky direction.")
    print("  prime -> irreducible polynomial is the Cohn digit-place channel.")
    print("  reducibility over Z -> recombination of looser atoms is the RFR channel.")
    print("  Tournament analogue: vertices should be witness slots or factor atoms;")
    print("  a Hamiltonian path is a recombination order, and directed cycles mark")
    print("  incompatible side-channel commitments.")


if __name__ == "__main__":
    main()
