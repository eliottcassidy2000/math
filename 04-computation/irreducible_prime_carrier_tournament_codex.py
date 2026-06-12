#!/usr/bin/env python3
"""Irreducible polynomials, prime values, and proof-carrier tournaments.

codex-2026-06-12

This script is a small computational atlas for the user's prompt:

    irreducible integer polynomials <-> prime values

It merges two papers:

* Jitender Singh, "Prime numbers and factorization of polynomials",
  arXiv:2411.18366.
* Shahriar Iravanian, "Integer Polynomial Factorization by Recombination of
  Real Factors", arXiv:2410.15880.

The purpose is not to prove Bunyakovsky's conjecture.  The purpose is to make
the two-way relationship testable inside the repo's carrier/tournament grammar:

* Bunyakovsky points forward: primitive irreducible polynomials should emit
  infinitely many integer primes.
* Murty/Singh and Cohn point backward: prime or low-prime-factor integer values
  certify irreducibility or bound the number of irreducible factors.
* Real-factor recombination turns factorization into a subset-selection problem,
  so irreducibility becomes absence of a nontrivial subset satisfying all
  integer trace/coefficient constraints.

Tournament Analysis is deliberately over proof carriers, not over runners,
arcs, or individual polynomials.  The assumption challenge is printed at the
end.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from functools import reduce
from itertools import combinations
from math import gcd

import sympy as sp


x = sp.Symbol("x")


def igcd(values: list[int]) -> int:
    return abs(reduce(gcd, [abs(v) for v in values], 0))


def bigomega(n: int) -> int:
    n = abs(int(n))
    if n <= 1:
        return 0
    return sum(sp.factorint(n).values())


def factor_degrees(poly: sp.Poly) -> list[int]:
    _unit, factors = sp.factor_list(poly, domain=sp.ZZ)
    degrees: list[int] = []
    for fac, exp in factors:
        degrees.extend([fac.degree()] * exp)
    return degrees


def fixed_divisor(poly: sp.Poly) -> int:
    # For an integer polynomial of degree d, gcd f(0),...,f(d) is the fixed
    # divisor gcd_{n in Z} f(n).
    d = poly.degree()
    return igcd([int(poly.eval(k)) for k in range(d + 1)])


def hf(poly: sp.Poly) -> sp.Rational:
    coeffs = poly.all_coeffs()
    leading = abs(int(coeffs[0]))
    lower = [abs(int(c)) for c in coeffs[1:]]
    return sp.Rational(max(lower) if lower else 0, leading)


def ceil_q(q: sp.Rational) -> int:
    return int(sp.ceiling(q))


def unitary_divisors_from_factorization(factors: dict[int, int]) -> list[int]:
    divs = [1]
    for p, e in factors.items():
        pe = p**e
        divs = divs + [d * pe for d in divs]
    return sorted(divs)


def singh_bound(poly: sp.Poly, max_m: int = 180) -> dict[str, object]:
    """Search a small Singh/Murty value-factor certificate.

    The d=1 case says that if m >= H_f+2, then the number of irreducible
    factors is at most Omega(|f(m)|).  The more flexible unitary-divisor version
    can strip a unitary divisor d provided m >= H_f+d+1.
    """
    h = hf(poly)
    factor_count = len(factor_degrees(poly))
    start = max(0, ceil_q(h + 2))
    best_plain: tuple[int, int, int] | None = None
    best_unitary: tuple[int, int, int, int] | None = None

    for m in range(start, max_m + 1):
        val = abs(int(poly.eval(m)))
        if val <= 1:
            continue
        omega = bigomega(val)
        if best_plain is None or (omega, m) < (best_plain[1], best_plain[0]):
            best_plain = (m, omega, val)

        factors = sp.factorint(val)
        for d in unitary_divisors_from_factorization(factors):
            if d == val:
                continue
            if m < h + d + 1:
                continue
            quotient = val // d
            if quotient <= 1:
                continue
            qomega = bigomega(quotient)
            if best_unitary is None or (qomega, m, d) < (
                best_unitary[2],
                best_unitary[0],
                best_unitary[1],
            ):
                best_unitary = (m, d, qomega, quotient)

    tight = None
    for m in range(start, max_m + 1):
        val = abs(int(poly.eval(m)))
        if val > 1 and bigomega(val) == factor_count:
            tight = (m, val)
            break

    return {
        "H": h,
        "threshold_d1": start,
        "plain_best": best_plain,
        "unitary_best": best_unitary,
        "tight_plain_m": tight,
    }


def prime_value_count(poly: sp.Poly, nmax: int = 160) -> tuple[int, list[tuple[int, int]]]:
    hits: list[tuple[int, int]] = []
    for n in range(nmax + 1):
        val = int(poly.eval(n))
        if val > 1 and sp.isprime(val):
            hits.append((n, val))
    return len(hits), hits[:8]


@dataclass(frozen=True)
class PolyRow:
    name: str
    expr: sp.Expr
    note: str


POLY_ROWS = [
    PolyRow("linear_x", x, "prime values are exactly the primes in the domain"),
    PolyRow("cohn_101_x2_plus_1", x**2 + 1, "Cohn from prime 101 in base 10"),
    PolyRow("euler_x2_minus_x_plus_41", x**2 - x + 41, "Euler/Singh prime window"),
    PolyRow("bunyakovsky_x2_plus_x_plus_1", x**2 + x + 1, "primitive irreducible candidate"),
    PolyRow("local_obstruction_x2_plus_x_plus_2", x**2 + x + 2, "irreducible but fixed divisor 2"),
    PolyRow("cohn_9841_repunit_base3", sum(x**i for i in range(9)), "Cohn bound tight: 9841=13*757"),
    PolyRow("singh_degree8_example", x**8 + 6 * x**7 + 5, "Singh example with x+1 factor"),
    PolyRow("irreducible_quartic_trace_rich", x**4 - 10 * x**2 + 1, "irreducible but many trace-zero subsets"),
    PolyRow("reducible_two_quadratics", (x**2 + x + 1) * (x**2 + 3 * x + 1), "recombination should find two factors"),
]


def digits_in_base(n: int, base: int) -> list[int]:
    if n == 0:
        return [0]
    digs: list[int] = []
    while n:
        digs.append(n % base)
        n //= base
    return digs


@dataclass(frozen=True)
class CohnRow:
    n: int
    base: int
    label: str


COHN_ROWS = [
    CohnRow(101, 10, "prime 101 -> x^2+1"),
    CohnRow(131, 10, "prime 131 -> x^2+3x+1"),
    CohnRow(9841, 3, "9841=13*757 -> 1+x+...+x^8"),
    CohnRow(2047, 2, "2047=23*89 -> binary repunit degree 10"),
]


def cohn_polynomial(row: CohnRow) -> sp.Poly:
    digits = digits_in_base(row.n, row.base)
    expr = sum(d * x**i for i, d in enumerate(digits))
    return sp.Poly(expr, x, domain=sp.ZZ)


def root_trace_blocks(poly: sp.Poly, prec: int = 50) -> list[float]:
    roots = [complex(r) for r in sp.nroots(poly.as_expr(), n=prec, maxsteps=200)]
    used = [False] * len(roots)
    blocks: list[float] = []
    tol = 10**-18
    for i, root in enumerate(roots):
        if used[i]:
            continue
        used[i] = True
        if abs(root.imag) < tol:
            blocks.append(root.real)
            continue
        target = complex(root.real, -root.imag)
        partner = None
        best_dist = float("inf")
        for j, other in enumerate(roots):
            if used[j]:
                continue
            dist = abs(other - target)
            if dist < best_dist:
                best_dist = dist
                partner = j
        if partner is None:
            blocks.append(root.real)
        else:
            used[partner] = True
            blocks.append((root + roots[partner]).real)
    return blocks


def trace_subset_candidates(poly: sp.Poly) -> dict[str, object]:
    blocks = root_trace_blocks(poly)
    n = len(blocks)
    candidate_sums: list[float] = []
    for mask in range(1, (1 << n) - 1):
        total = sum(blocks[i] for i in range(n) if mask & (1 << i))
        if abs(total - round(total)) < 10**-8:
            candidate_sums.append(total)
    true_count = len(factor_degrees(poly))
    true_nontrivial_subsets = 0 if true_count <= 1 else (2**true_count - 2)
    return {
        "blocks": [round(v, 10) for v in blocks],
        "integer_trace_candidates": len(candidate_sums),
        "sample_candidate_sums": [int(round(v)) for v in candidate_sums[:8]],
        "true_nontrivial_factor_subsets": true_nontrivial_subsets,
        "false_trace_excess": max(0, len(candidate_sums) - true_nontrivial_subsets),
    }


def hamiltonian_paths(adj: list[list[bool]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if not count:
                continue
            for nxt in range(n):
                bit = 1 << nxt
                if mask & bit:
                    continue
                if adj[last][nxt]:
                    dp[mask | bit][nxt] += count
    return sum(dp[-1])


def tournament_fingerprint(adj: list[list[bool]]) -> dict[str, object]:
    n = len(adj)
    scores = [sum(row) for row in adj]
    cycles = 0
    for a, b, c in combinations(range(n), 3):
        out = [0, 0, 0]
        triple = [a, b, c]
        for i, j in combinations(range(3), 2):
            u, v = triple[i], triple[j]
            if adj[u][v]:
                out[i] += 1
            else:
                out[j] += 1
        if sorted(out) == [1, 1, 1]:
            cycles += 1

    def reach(start: int, reverse: bool = False) -> set[int]:
        seen = {start}
        q: deque[int] = deque([start])
        while q:
            u = q.popleft()
            for v in range(n):
                edge = adj[v][u] if reverse else adj[u][v]
                if edge and v not in seen:
                    seen.add(v)
                    q.append(v)
        return seen

    remaining = set(range(n))
    scc_sizes: list[int] = []
    while remaining:
        start = min(remaining)
        comp = reach(start) & reach(start, reverse=True)
        scc_sizes.append(len(comp))
        remaining -= comp

    return {
        "score_hist": dict(sorted(Counter(scores).items())),
        "directed_3cycles": cycles,
        "scc_sizes": sorted(scc_sizes, reverse=True),
        "hamiltonian_paths": hamiltonian_paths(adj),
    }


@dataclass(frozen=True)
class Carrier:
    name: str
    forward_primes: int
    reverse_irreducibility: int
    recombination: int
    support_retention: int
    repo_bridge: int
    risk: int
    note: str

    def votes_against(self, other: "Carrier") -> int:
        wins = 0
        for field in (
            "forward_primes",
            "reverse_irreducibility",
            "recombination",
            "support_retention",
            "repo_bridge",
        ):
            if getattr(self, field) > getattr(other, field):
                wins += 1
            elif getattr(self, field) == getattr(other, field) and self.risk < other.risk:
                wins += 1
        return wins


CARRIERS = [
    Carrier(
        "bunyakovsky_forward_atom",
        5,
        1,
        1,
        2,
        4,
        5,
        "irreducible primitive f should hit infinitely many primes",
    ),
    Carrier(
        "singh_murty_value_factor_bound",
        2,
        5,
        2,
        4,
        5,
        1,
        "large value plus prime factorization bounds irreducible factors",
    ),
    Carrier(
        "cohn_digit_prime_certificate",
        2,
        5,
        1,
        3,
        4,
        1,
        "base-b digits of prime N give irreducible digit polynomial",
    ),
    Carrier(
        "iravanian_real_trace_recombination",
        1,
        3,
        5,
        4,
        4,
        2,
        "factor over R, then solve integer trace/subset constraints",
    ),
    Carrier(
        "finite_field_hensel_recombination",
        1,
        4,
        4,
        3,
        3,
        1,
        "standard modular factor-lift-recombine pipeline",
    ),
    Carrier(
        "newton_nonarch_bivariate_gate",
        1,
        4,
        2,
        5,
        4,
        2,
        "Singh bivariate/non-Archimedean extension keeps degree support",
    ),
    Carrier(
        "local_fixed_divisor_sieve",
        4,
        2,
        1,
        4,
        5,
        2,
        "fixed divisor / local obstruction is Bunyakovsky's admissibility gate",
    ),
    Carrier(
        "repo_support_gate_LRC72",
        2,
        3,
        2,
        5,
        5,
        2,
        "same scalar/support split seen in LRC14 and [72,36,16]",
    ),
]


def carrier_tournament(carriers: list[Carrier]) -> tuple[list[list[bool]], dict[str, object]]:
    n = len(carriers)
    adj = [[False] * n for _ in range(n)]
    edge_flips_vs_reverse_only = 0
    for i, j in combinations(range(n), 2):
        a, b = carriers[i], carriers[j]
        av = a.votes_against(b)
        bv = b.votes_against(a)
        if av > bv or (av == bv and (a.repo_bridge, -a.risk, a.name) > (b.repo_bridge, -b.risk, b.name)):
            winner, loser = i, j
        else:
            winner, loser = j, i
        adj[winner][loser] = True

        reverse_winner = i if (a.reverse_irreducibility, -a.risk) >= (b.reverse_irreducibility, -b.risk) else j
        if reverse_winner != winner:
            edge_flips_vs_reverse_only += 1

    fp = tournament_fingerprint(adj)
    fp["edge_flips_vs_reverse_irreducibility_only"] = edge_flips_vs_reverse_only
    scores = [(carriers[i].name, sum(adj[i])) for i in range(n)]
    fp["ranking"] = [name for name, _score in sorted(scores, key=lambda t: (-t[1], t[0]))]
    return adj, fp


def print_poly_atlas() -> None:
    print("POLYNOMIAL ATLAS")
    print("================")
    for row in POLY_ROWS:
        poly = sp.Poly(row.expr, x, domain=sp.ZZ)
        degrees = factor_degrees(poly)
        fdiv = fixed_divisor(poly)
        pcount, first_hits = prime_value_count(poly)
        bound = singh_bound(poly)
        print(f"\n[{row.name}] {sp.sstr(poly.as_expr())}")
        print(f"  note: {row.note}")
        print(f"  factor_degrees={degrees}; fixed_divisor={fdiv}; bunyakovsky_admissible={len(degrees)==1 and fdiv==1}")
        print(f"  prime_values_in_n_0_160={pcount}; first_hits={first_hits}")
        print(f"  H_f={bound['H']}; d=1 threshold m>={bound['threshold_d1']}")
        print(f"  Singh plain best (m,Omega,|f(m)|)={bound['plain_best']}")
        print(f"  Singh unitary best (m,d,Omega(|f(m)|/d),quotient)={bound['unitary_best']}")
        print(f"  first plain tight m with Omega(|f(m)|)=factor_count: {bound['tight_plain_m']}")


def print_cohn_atlas() -> None:
    print("\nCOHN DIGIT-POLYNOMIAL ATLAS")
    print("===========================")
    for row in COHN_ROWS:
        poly = cohn_polynomial(row)
        degrees = factor_degrees(poly)
        omega_n = bigomega(row.n)
        print(f"\n[{row.label}] N={row.n}, base={row.base}, digits={digits_in_base(row.n, row.base)}")
        print(f"  digit_poly={sp.sstr(poly.as_expr())}")
        print(f"  Omega(N)={omega_n}; factor_degrees={degrees}; bound_ok={len(degrees) <= omega_n}")


def print_recombination_atlas() -> None:
    print("\nREAL-FACTOR RECOMBINATION SCOUT")
    print("===============================")
    examples = [
        ("irreducible trace-rich quartic", x**4 - 10 * x**2 + 1),
        ("irreducible generic quartic", x**4 + x + 1),
        ("reducible two quadratics", (x**2 + x + 1) * (x**2 + 3 * x + 1)),
        ("Euler quadratic", x**2 - x + 41),
    ]
    for name, expr in examples:
        poly = sp.Poly(expr, x, domain=sp.ZZ)
        data = trace_subset_candidates(poly)
        print(f"\n[{name}] {sp.sstr(poly.as_expr())}")
        print(f"  factor_degrees={factor_degrees(poly)}")
        print(f"  trace_blocks={data['blocks']}")
        print(f"  integer_trace_candidates={data['integer_trace_candidates']}")
        print(f"  true_nontrivial_factor_subsets={data['true_nontrivial_factor_subsets']}")
        print(f"  false_trace_excess={data['false_trace_excess']}")
        print(f"  sample_candidate_sums={data['sample_candidate_sums']}")


def print_carrier_tournament() -> None:
    print("\nPROOF-CARRIER TOURNAMENT")
    print("========================")
    print("Pairwise observable: majority vote over")
    print("  forward_primes, reverse_irreducibility, recombination, support_retention, repo_bridge")
    print("Tie path: lower risk, then declaration/name.")
    adj, fp = carrier_tournament(CARRIERS)
    for i, carrier in enumerate(CARRIERS):
        print(
            f"  {carrier.name}: score={sum(adj[i])}, metrics="
            f"({carrier.forward_primes},{carrier.reverse_irreducibility},"
            f"{carrier.recombination},{carrier.support_retention},{carrier.repo_bridge}; risk={carrier.risk})"
        )
        print(f"    note: {carrier.note}")
    print(f"fingerprint={fp}")


def print_infinite_tournament_hypothesis() -> None:
    print("\nINFINITE-TOURNAMENT HYPOTHESIS")
    print("==============================")
    print("One useful infinite tournament is not on integers themselves but on certificate states.")
    print("A vertex can be (f, local data retained): fixed divisor, residue obstructions,")
    print("least Singh/Cohn certificate depth, trace-subset survivor set, Newton polygon edge data.")
    print("Orient A -> B when A has the smaller retained certificate depth after normalizing degree")
    print("and fixed divisor; orient ties by the richer side channel.  Finite truncations then")
    print("record edge flips as the prime window grows.  This makes Bunyakovsky a statement")
    print("about an infinite supply of forward prime-hit witnesses, while Singh/Cohn/Iravanian")
    print("supply reverse and algorithmic edges in the same tournament.")


def print_assumption_challenge() -> None:
    print("\nASSUMPTION CHALLENGE")
    print("====================")
    print("Alternate vertex sets considered:")
    print("  polynomials; integer values; primes; residues modulo p; fixed-divisor obstructions;")
    print("  roots; real/complex conjugate blocks; recombination subsets; digit positions;")
    print("  Newton polygon edges; bivariate specialization channels; LRC runners/gaps;")
    print("  codewords, matroid circuits, minimum-design blocks, and proof obligations.")
    print("Chosen vertices: proof carriers.")
    print("Preserved: the two-way prime/irreducible direction, support retention,")
    print("computability, and transfer potential to LRC14 and [72,36,16].")
    print("Destroyed: individual polynomial arithmetic, actual asymptotic prime production,")
    print("root geometry beyond first trace, and code-support realizability.")
    print("Challenged assumption: the tournament need not live on runners, arcs, or primes;")
    print("it can live on certificate channels and still preserve the obstruction predicate.")


def main() -> None:
    print("IRREDUCIBLE PRIME CARRIER TOURNAMENT ATLAS")
    print("==========================================")
    print("Sources:")
    print("  Singh arXiv:2411.18366 - prime factorization of large values bounds factors.")
    print("  Iravanian arXiv:2410.15880 - real-factor recombination becomes subset-sum.")
    print()
    print_poly_atlas()
    print_cohn_atlas()
    print_recombination_atlas()
    print_carrier_tournament()
    print_infinite_tournament_hypothesis()
    print_assumption_challenge()


if __name__ == "__main__":
    main()
