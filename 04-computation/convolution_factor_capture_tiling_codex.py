#!/usr/bin/env python3
"""Hidden convolution tilings for the prime <-> irreducible bridge.

codex-2026-06-12

This is an extension of HYP-2450.  There the user's triangular coefficient
picture was installed as a diagonal quotient of fixed-path tournament tilings.
Here the quotient is lifted one layer downward:

    target coefficient vector a_k
        asks for hidden factor vectors b_i,c_j with
    a_k = sum_{i+j=k} b_i c_j.

For primitive polynomials of degree at most 5, a reducible polynomial has a
linear or quadratic factor.  The script therefore implements an exact
convolution-lift search in that range, then wraps the lift with:

* factor-capture witness scores from integer values f(m),
* residue-class tournaments modulo small primes,
* sign-cube chamber tournaments at fixed coefficient magnitudes, and
* a proof-route tournament whose tie path is the arbitrary fixed Hamiltonian
  path, not mathematical evidence.

The output is a reproducible hypothesis generator.  It does not prove
Bunyakovsky or any new general irreducibility criterion.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import reduce
from itertools import combinations, product
from math import gcd, prod, sqrt
from typing import Iterable

import sympy as sp


x = sp.Symbol("x")
t = sp.Symbol("t")


def trim(coeffs: Iterable[int]) -> tuple[int, ...]:
    out = list(coeffs)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return tuple(int(c) for c in out)


def poly_from_coeffs(coeffs: tuple[int, ...]) -> sp.Poly:
    return sp.Poly(sum(c * x**i for i, c in enumerate(coeffs)), x, domain=sp.ZZ)


def coeffs_from_poly(poly: sp.Poly) -> tuple[int, ...]:
    return trim(reversed([int(c) for c in poly.all_coeffs()]))


def expr_str(coeffs: tuple[int, ...]) -> str:
    return sp.sstr(poly_from_coeffs(coeffs).as_expr())


def signed_divisors(n: int) -> list[int]:
    if n == 0:
        raise ValueError("zero has infinitely many divisors in this scout")
    ds: list[int] = []
    for d in sp.divisors(abs(int(n))):
        ds.extend([-int(d), int(d)])
    return sorted(ds)


def coeff_gcd(coeffs: tuple[int, ...]) -> int:
    return reduce(gcd, (abs(c) for c in coeffs), 0)


def primitive_coeffs(coeffs: tuple[int, ...]) -> tuple[int, ...]:
    g = coeff_gcd(coeffs)
    if g <= 1:
        return trim(coeffs)
    return trim(c // g for c in coeffs)


def convolve(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return trim(out)


def divide_if_exact(f: tuple[int, ...], g: tuple[int, ...]) -> tuple[int, ...] | None:
    """Return h with f=g*h if h has integral coefficients.

    The low-to-high recurrence is deliberately the coefficient-tiling equation:
    after the constant coefficient of g is fixed, each next coefficient of h is
    forced by a diagonal sum.
    """
    f = trim(f)
    g = trim(g)
    if len(g) < 2 or g[0] == 0 or len(g) > len(f):
        return None
    q_len = len(f) - len(g) + 1
    q: list[int] = []
    for k in range(q_len):
        rhs = f[k]
        for i in range(1, min(k + 1, len(g))):
            rhs -= g[i] * q[k - i]
        if rhs % g[0] != 0:
            return None
        q.append(rhs // g[0])
    h = trim(q)
    if convolve(g, h) == f:
        return h
    return None


def normalize_pair(g: tuple[int, ...], h: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    g = trim(g)
    h = trim(h)
    if (len(g), g) > (len(h), h):
        g, h = h, g
    if g[-1] < 0:
        g = tuple(-c for c in g)
        h = tuple(-c for c in h)
    if (len(g), g) > (len(h), h):
        g, h = h, g
    return g, h


def linear_lifts(coeffs: tuple[int, ...]) -> set[tuple[tuple[int, ...], tuple[int, ...]]]:
    coeffs = trim(coeffs)
    if len(coeffs) < 3 or coeffs[0] == 0:
        return set()
    lifts: set[tuple[tuple[int, ...], tuple[int, ...]]] = set()
    for g0 in signed_divisors(coeffs[0]):
        for g1 in signed_divisors(coeffs[-1]):
            g = (g0, g1)
            h = divide_if_exact(coeffs, g)
            if h is not None and len(h) >= 2:
                lifts.add(normalize_pair(g, h))
    return lifts


def integer_roots(poly: sp.Poly) -> set[int]:
    if poly.is_zero or poly.degree() <= 0:
        return set()
    roots: set[int] = set()
    for r, _mult in sp.roots(poly.as_expr(), t).items():
        if bool(r.is_integer):
            roots.add(int(r))

    cleared = sp.Poly(poly.as_expr(), t, domain=sp.QQ).clear_denoms()[1]
    const = int(cleared.nth(0))
    candidates = {0} if const == 0 else set(signed_divisors(const))
    for cand in candidates:
        if cleared.eval(cand) == 0:
            roots.add(int(cand))
    return roots


def mignotteish_bound(coeffs: tuple[int, ...], factor_degree: int) -> int:
    norm = sqrt(sum(c * c for c in coeffs))
    return int((2 ** factor_degree) * norm + 3)


def candidate_b1_degree4(coeffs: tuple[int, ...], b0: int, b2: int) -> set[int]:
    a0, a1, _a2, a3, a4 = coeffs
    c0 = a0 // b0
    c2 = a4 // b2
    denom = b0 * c2 - b2 * c0
    numer = b0 * a3 - b2 * a1
    if denom:
        return {numer // denom} if numer % denom == 0 else set()
    if numer:
        return set()
    bound = mignotteish_bound(coeffs, 2)
    return set(range(-bound, bound + 1))


def candidate_b1_degree5(coeffs: tuple[int, ...], b0: int, b2: int) -> set[int]:
    a0, a1, a2, _a3, a4, a5 = coeffs
    c0 = a0 // b0
    c3 = a5 // b2
    # From c1=(a1-b1*c0)/b0, c2=(a4-b1*c3)/b2 and
    # a2=b0*c2+b1*c1+b2*c0, after multiplying by b0*b2:
    # A*b1^2 + B*b1 + C = 0.
    A = -b2 * c0
    B = -b0 * b0 * c3 + b2 * a1
    C = b0 * b0 * a4 + b0 * b2 * b2 * c0 - a2 * b0 * b2
    if A == 0 and B == 0:
        if C:
            return set()
        bound = mignotteish_bound(coeffs, 2)
        return set(range(-bound, bound + 1))
    if A == 0:
        return {-C // B} if C % B == 0 else set()
    return integer_roots(sp.Poly(A * t**2 + B * t + C, t, domain=sp.ZZ))


def quadratic_lifts(coeffs: tuple[int, ...]) -> set[tuple[tuple[int, ...], tuple[int, ...]]]:
    """Find quadratic-factor convolution lifts for nonzero constant terms.

    For degree 4 and 5 this is exact: after fixing the boundary tiles b0,b2,
    the remaining unknown b1 is forced by inward diagonal equations.  For
    larger sample polynomials the bounded fallback is only a scout.
    """
    coeffs = trim(coeffs)
    n = len(coeffs) - 1
    if n < 4 or coeffs[0] == 0:
        return set()
    lifts: set[tuple[tuple[int, ...], tuple[int, ...]]] = set()
    for b0 in signed_divisors(coeffs[0]):
        for b2 in signed_divisors(coeffs[-1]):
            if n == 4:
                roots = candidate_b1_degree4(coeffs, b0, b2)
            elif n == 5:
                roots = candidate_b1_degree5(coeffs, b0, b2)
            else:
                # Larger rows are not used for exact certification in this
                # script; the cap keeps Cohn samples cheap and deterministic.
                bound = min(mignotteish_bound(coeffs, 2), 200)
                roots = set(range(-bound, bound + 1))
            for b1 in roots:
                g = (int(b0), int(b1), int(b2))
                h = divide_if_exact(coeffs, g)
                if h is not None and len(h) >= 2:
                    lifts.add(normalize_pair(g, h))
    return lifts


def convolution_lifts_degree_le5(coeffs: tuple[int, ...]) -> set[tuple[tuple[int, ...], tuple[int, ...]]]:
    """Exact reducibility witness search for primitive degree <= 5 polynomials."""
    coeffs = primitive_coeffs(coeffs)
    if len(coeffs) - 1 > 5:
        raise ValueError("linear/quadratic lift search is complete only through degree 5")
    return linear_lifts(coeffs) | quadratic_lifts(coeffs)


def factor_degrees(coeffs: tuple[int, ...]) -> list[int]:
    poly = poly_from_coeffs(primitive_coeffs(coeffs))
    if poly.degree() <= 0:
        return []
    _unit, factors = sp.factor_list(poly.as_expr(), x)
    degrees: list[int] = []
    for fac, exp in factors:
        degrees.extend([sp.Poly(fac, x).degree()] * exp)
    return sorted(degrees)


def is_irreducible(coeffs: tuple[int, ...]) -> bool:
    return factor_degrees(coeffs) == [len(trim(primitive_coeffs(coeffs))) - 1]


def fixed_divisor(coeffs: tuple[int, ...]) -> int:
    poly = poly_from_coeffs(coeffs)
    vals = [abs(int(poly.eval(k))) for k in range(poly.degree() + 1)]
    return reduce(gcd, vals, 0)


def prime_hits(coeffs: tuple[int, ...], limit: int = 40) -> list[tuple[int, int]]:
    poly = poly_from_coeffs(coeffs)
    hits: list[tuple[int, int]] = []
    for m in range(limit + 1):
        val = int(poly.eval(m))
        if val > 1 and sp.isprime(val):
            hits.append((m, val))
    return hits


def omega(n: int) -> int:
    if n == 0:
        return 10**9
    return sum(sp.factorint(abs(int(n))).values())


def allocation_count_two(n: int) -> int:
    if abs(n) <= 1:
        return 0
    exps = sp.factorint(abs(int(n))).values()
    return prod(e + 1 for e in exps) - 2


@dataclass(frozen=True)
class CaptureWitness:
    m: int
    value: int
    omega: int
    allocation_count: int
    factor_values: tuple[int, ...]
    factor_degrees: tuple[int, ...]


def best_capture_witness(coeffs: tuple[int, ...], limit: int = 50) -> CaptureWitness | None:
    poly = poly_from_coeffs(coeffs)
    _unit, facs = sp.factor_list(poly_from_coeffs(primitive_coeffs(coeffs)).as_expr(), x)
    best: CaptureWitness | None = None
    for m in range(1, limit + 1):
        val = int(poly.eval(m))
        if val <= 1:
            continue
        factor_values: list[int] = []
        factor_degrees_at_m: list[int] = []
        usable = True
        for fac, exp in facs:
            fpoly = sp.Poly(fac, x)
            fval = abs(int(fpoly.eval(m)))
            if fval <= 1:
                usable = False
                break
            factor_values.extend([fval] * exp)
            factor_degrees_at_m.extend([fpoly.degree()] * exp)
        if not usable:
            continue
        witness = CaptureWitness(
            m=m,
            value=val,
            omega=omega(val),
            allocation_count=allocation_count_two(val),
            factor_values=tuple(factor_values),
            factor_degrees=tuple(factor_degrees_at_m),
        )
        key = (witness.omega, witness.allocation_count, witness.m)
        if best is None or key < (best.omega, best.allocation_count, best.m):
            best = witness
    return best


def vp(n: int, p: int, cap: int = 8) -> int:
    n = abs(int(n))
    if n == 0:
        return cap + 1
    count = 0
    while n % p == 0 and count < cap:
        count += 1
        n //= p
    return count


def tournament_from_scores(scores: list[tuple[int, int]]) -> list[list[bool]]:
    """Lower score wins; the second coordinate is the fixed Hamiltonian path."""
    n = len(scores)
    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        if scores[i] < scores[j]:
            adj[i][j] = True
        else:
            adj[j][i] = True
    return adj


def directed_3cycles(adj: list[list[bool]]) -> int:
    count = 0
    for a, b, c in combinations(range(len(adj)), 3):
        edges = int(adj[a][b]) + int(adj[b][c]) + int(adj[c][a])
        if edges in {0, 3}:
            count += 1
    return count


def scc_sizes(adj: list[list[bool]]) -> list[int]:
    n = len(adj)

    def reach_from(start: int, graph: list[list[bool]]) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            v = stack.pop()
            for w, edge in enumerate(graph[v]):
                if edge and w not in seen:
                    seen.add(w)
                    stack.append(w)
        return seen

    radj = [[adj[j][i] for j in range(n)] for i in range(n)]
    remaining = set(range(n))
    sizes: list[int] = []
    while remaining:
        v = next(iter(remaining))
        comp = reach_from(v, adj) & reach_from(v, radj)
        sizes.append(len(comp))
        remaining -= comp
    return sorted(sizes, reverse=True)


def hamiltonian_paths(adj: list[list[bool]]) -> int:
    n = len(adj)
    if n > 12:
        return -1
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
    score_hist = Counter(sum(row) for row in adj)
    return {
        "score_hist": dict(sorted(score_hist.items())),
        "directed_3cycles": directed_3cycles(adj),
        "scc_sizes": scc_sizes(adj),
        "hamiltonian_paths": hamiltonian_paths(adj),
    }


def residue_tournament_summary(coeffs: tuple[int, ...], p: int) -> dict[str, object]:
    poly = poly_from_coeffs(coeffs)
    vals = [vp(int(poly.eval(r)), p) for r in range(p)]
    scores = [(vals[r], r) for r in range(p)]
    adj = tournament_from_scores(scores)
    return {
        "p": p,
        "valuations": vals,
        "all_bad_mod_p": all(v > 0 for v in vals),
        **tournament_fingerprint(adj),
    }


def scan_family(degree: int, max_abs: int) -> dict[str, object]:
    values = [v for v in range(-max_abs, max_abs + 1) if v]
    rows = 0
    reducible = 0
    irreducible = 0
    lift_yes = 0
    mismatches: list[tuple[tuple[int, ...], str, list[int], int]] = []
    fixed_divisor_blocked = 0
    examples: dict[str, object] = {}
    for leading in range(1, max_abs + 1):
        for rest in product(values, repeat=degree):
            coeffs = tuple(rest) + (leading,)
            if coeff_gcd(coeffs) != 1:
                continue
            rows += 1
            degs = factor_degrees(coeffs)
            truth_reducible = len(degs) > 1
            lifts = convolution_lifts_degree_le5(coeffs)
            has_lift = bool(lifts)
            reducible += int(truth_reducible)
            irreducible += int(not truth_reducible)
            lift_yes += int(has_lift)
            fixed_divisor_blocked += int(fixed_divisor(coeffs) > 1)
            if truth_reducible and "reducible_lift" not in examples:
                examples["reducible_lift"] = {
                    "poly": expr_str(coeffs),
                    "degrees": degs,
                    "lift": tuple(expr_str(f) for f in next(iter(lifts))),
                }
            if not truth_reducible and fixed_divisor(coeffs) > 1 and "irreducible_fixed_divisor" not in examples:
                examples["irreducible_fixed_divisor"] = {
                    "poly": expr_str(coeffs),
                    "fixed_divisor": fixed_divisor(coeffs),
                    "prime_hits_0_40": len(prime_hits(coeffs)),
                }
            if truth_reducible != has_lift and len(mismatches) < 5:
                mismatches.append((coeffs, expr_str(coeffs), degs, len(lifts)))
    return {
        "degree": degree,
        "max_abs": max_abs,
        "primitive_rows": rows,
        "sympy_reducible": reducible,
        "sympy_irreducible": irreducible,
        "convolution_lift_yes": lift_yes,
        "mismatches": mismatches,
        "fixed_divisor_blocked": fixed_divisor_blocked,
        "examples": examples,
    }


def sample_rows() -> dict[str, tuple[int, ...]]:
    return {
        "fixed_divisor_block_x2_x_2": (2, 1, 1),
        "quartic_irred_all_minus": (-1, -1, -1, -1, 1),
        "quartic_reducible_hyp2449": (-1, -1, -1, -2, 2),
        "trace_false_quartic": (1, 0, -10, 0, 1),
        "parity_min_slice_plus": (1, 1, 0, 1, 0, 1),
        "cohn_9841_repunit8": tuple([1] * 9),
        "cohn_2047_repunit10": tuple([1] * 11),
    }


def summarize_samples() -> dict[str, object]:
    out: dict[str, object] = {}
    for name, coeffs in sample_rows().items():
        deg = len(coeffs) - 1
        if deg <= 5:
            lifts = convolution_lifts_degree_le5(coeffs)
        else:
            lifts = linear_lifts(primitive_coeffs(coeffs)) | quadratic_lifts(primitive_coeffs(coeffs))
        witness = best_capture_witness(coeffs, limit=50)
        out[name] = {
            "poly": expr_str(coeffs),
            "degree": deg,
            "factor_degrees": factor_degrees(coeffs),
            "fixed_divisor": fixed_divisor(coeffs),
            "prime_hits_0_40": prime_hits(coeffs, 40)[:8],
            "lift_count_linear_or_quadratic": len(lifts),
            "first_lift": tuple(expr_str(f) for f in next(iter(lifts))) if lifts else None,
            "best_capture": witness.__dict__ if witness else None,
            "residue_tournaments": [residue_tournament_summary(coeffs, p) for p in (2, 3, 5, 7)],
        }
    return out


def sign_cube_summary(magnitudes: tuple[int, ...], base: int, hit_limit: int = 30) -> dict[str, object]:
    rows_by_coeffs: dict[tuple[int, ...], dict[str, object]] = {}
    for signs in product((-1, 1), repeat=len(magnitudes)):
        coeffs = tuple(s * m for s, m in zip(signs, magnitudes))
        if coeffs[-1] == 0:
            continue
        poly = poly_from_coeffs(coeffs)
        rows_by_coeffs.setdefault(
            coeffs,
            {
                "signs": signs,
                "coeffs": coeffs,
                "value": int(poly.eval(base)),
                "irreducible": is_irreducible(coeffs),
                "prime_hits": len(prime_hits(coeffs, hit_limit)),
                "fixed_divisor": fixed_divisor(coeffs),
            },
        )
    rows = list(rows_by_coeffs.values())
    scores = [(int(row["value"]), idx) for idx, row in enumerate(rows)]
    adj = tournament_from_scores(scores)
    value_fp = tournament_fingerprint(adj)
    irreducible_rows = [row for row in rows if row["irreducible"]]
    prime_rows = [row for row in rows if row["prime_hits"]]
    fixed_blocked = [row for row in rows if int(row["fixed_divisor"]) > 1]
    value_buckets: dict[int, Counter] = defaultdict(Counter)
    for row in rows:
        value_buckets[int(row["value"])]["irreducible" if row["irreducible"] else "reducible"] += 1
    mixed_value_buckets = sum(1 for c in value_buckets.values() if len(c) > 1)
    top_by_value = sorted(rows, key=lambda r: (-int(r["value"]), r["signs"]))[:3]
    top_by_prime_hits = sorted(rows, key=lambda r: (-int(r["prime_hits"]), r["signs"]))[:3]
    return {
        "magnitudes": magnitudes,
        "base": base,
        "raw_sign_patterns": 2 ** len(magnitudes),
        "distinct_chambers": len(rows),
        "distinct_base_values": len(value_buckets),
        "irreducible": len(irreducible_rows),
        "has_prime_hit": len(prime_rows),
        "fixed_divisor_blocked": len(fixed_blocked),
        "mixed_value_buckets": mixed_value_buckets,
        "value_tournament": value_fp,
        "top_by_value": [
            {"poly": expr_str(row["coeffs"]), "value": row["value"], "prime_hits": row["prime_hits"]}
            for row in top_by_value
        ],
        "top_by_prime_hits": [
            {"poly": expr_str(row["coeffs"]), "value": row["value"], "prime_hits": row["prime_hits"]}
            for row in top_by_prime_hits
        ],
    }


def proof_route_tournament() -> dict[str, object]:
    routes = {
        "convolution_lift_disprover": (5, 5, 5, 5, 4, 3),
        "factor_capture_hypertournament": (4, 5, 5, 4, 4, 4),
        "residue_winner_tournaments": (4, 5, 4, 4, 3, 4),
        "sign_cube_chamber_tournament": (4, 5, 4, 5, 3, 5),
        "newton_slope_boundary_lift": (5, 3, 5, 4, 4, 4),
        "real_factor_recombination_gpu": (4, 3, 5, 3, 5, 3),
        "cohn_diagonal_count_certificate": (5, 5, 4, 5, 3, 2),
        "code_support_weight_enumerator_fiber": (3, 3, 4, 3, 5, 4),
        "lrc14_denominator_resource_lift": (3, 4, 4, 3, 5, 4),
    }
    names = list(routes)
    adj = [[False] * len(names) for _ in names]
    for i, j in combinations(range(len(names)), 2):
        ri, rj = routes[names[i]], routes[names[j]]
        votes_i = sum(a > b for a, b in zip(ri, rj))
        votes_j = sum(b > a for a, b in zip(ri, rj))
        if votes_i > votes_j:
            adj[i][j] = True
        elif votes_j > votes_i:
            adj[j][i] = True
        else:
            # Arbitrary fixed Hamiltonian path tiebreak, per the user's note.
            adj[min(i, j)][max(i, j)] = True
    ranking = sorted(range(len(names)), key=lambda i: sum(adj[i]), reverse=True)
    return {
        "criteria": "exactness,testability,prime_bridge,tiling_bridge,LRC72_transfer,search_breadth",
        "routes": routes,
        **tournament_fingerprint(adj),
        "ranking": [names[i] for i in ranking],
    }


def main() -> None:
    print("CONVOLUTION FACTOR-CAPTURE TILING PILOT")
    print("=======================================")
    print("HYP-2452/T796/OPEN-Q-074 addendum to HYP-2450 and HYP-2451.")
    print()
    print("SOURCE ANCHORS")
    print("==============")
    print("Singh arXiv:2411.18366: value factorizations bound polynomial factor complexity.")
    print("Iravanian arXiv:2410.15880: real-factor recombination as subset/trace assembly.")
    print("Abu Salem-Gao-Lauder: Newton polytope decomposition + boundary-to-interior factor lift.")
    print()

    print("EXACT CONVOLUTION-LIFT SCANS")
    print("============================")
    for degree, max_abs in ((4, 3), (5, 2)):
        summary = scan_family(degree, max_abs)
        print(summary)
    print()

    print("SAMPLE FACTOR-CAPTURE AND RESIDUE TOURNAMENTS")
    print("=============================================")
    samples = summarize_samples()
    for name, data in samples.items():
        compact = {
            "poly": data["poly"],
            "degree": data["degree"],
            "factor_degrees": data["factor_degrees"],
            "fixed_divisor": data["fixed_divisor"],
            "prime_hits_0_40": data["prime_hits_0_40"],
            "lift_count_linear_or_quadratic": data["lift_count_linear_or_quadratic"],
            "first_lift": data["first_lift"],
            "best_capture": data["best_capture"],
            "residue_digest": [
                {
                    "p": r["p"],
                    "valuations": r["valuations"],
                    "all_bad_mod_p": r["all_bad_mod_p"],
                    "score_hist": r["score_hist"],
                    "cycles": r["directed_3cycles"],
                    "scc": r["scc_sizes"],
                    "Hpaths": r["hamiltonian_paths"],
                }
                for r in data["residue_tournaments"]
            ],
        }
        print(f"{name}: {compact}")
    print()

    print("SIGN-CUBE CHAMBER TOURNAMENTS")
    print("=============================")
    for mags, base in (((1, 1, 1, 1, 1, 1), 5), ((1, 1, 0, 1, 0, 1), 5), ((1, 2, 1, 2, 1, 1), 7)):
        print(sign_cube_summary(mags, base))
    print()

    print("PROOF-ROUTE TOURNAMENT")
    print("======================")
    print(proof_route_tournament())
    print()

    print("ASSUMPTION CHALLENGE")
    print("====================")
    print("Alternate vertex sets considered: factor coefficients, convolution grid cells,")
    print("diagonal sums, witness integers m, prime tokens of f(m), residue classes,")
    print("sign chambers, Newton boundary edges, real-factor recombination subsets,")
    print("LRC denominator resources, and [72,36,16] support/design obligations.")
    print("Chosen quotient for this pilot: hidden 2D coefficient convolution lifts.")
    print("Preserved: exact reducibility witnesses through degree 5, diagonal-sum")
    print("tiling geometry, local residue winners, and witness-value factor complexity.")
    print("Destroyed: root geometry, Galois action, large-degree factor search, and")
    print("asymptotic prime-value information.  Those are side channels to reattach.")
    print("Challenged assumption: the tournament vertices need not be polynomials,")
    print("arcs, or primes; the best vertices here are proof obligations and hidden")
    print("factor tiles, while tiebreaks are only the fixed Hamiltonian path.")
    print()

    print("NEXT DIRECTIONS")
    print("===============")
    print("1. Add bounded ILP/SAT for degree >=6 convolution lifts with coefficient")
    print("   bounds, using the exact degree<=5 solver as the regression oracle.")
    print("2. Attach Newton-slope boundary layers to sparse/multivariate examples:")
    print("   boundary factorisation first, then inward quadratic constraints.")
    print("3. Use factor-capture witnesses as search pruning: low Omega(f(m)) means")
    print("   few token allocations and few possible irreducible factor slots.")
    print("4. Transfer the scalar/fiber split to [72,36,16]: weight enumerator")
    print("   coefficients are boundary totals, while support/design incidence is")
    print("   the hidden lift that must exist.")


if __name__ == "__main__":
    main()
