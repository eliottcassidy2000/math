#!/usr/bin/env python3
"""Coefficient-layer tilings for the prime <-> irreducible-polynomial bridge.

codex-2026-06-12

This script follows the user's suggested picture:

  degree n, coefficients a_0..a_n, and a triangular tournament/tiling array.

For N=n+1 vertices, the gap-d diagonal has N-d arcs.  The user's layer counts
for n=5 are exactly 1,2,3,4,5 from the apex gap down to the base path.

The script tests three ways tournaments can hide inside the bridge:

1. Diagonal-count Cohn map:
   coefficient c_d = number of forward arcs at gap d.  In base b>=N, the
   profile integer has valid base-b digits.  If it is prime, Cohn gives an
   irreducible digit polynomial.

2. Centered magnetization map:
   coefficient A_d = #forward - #backward on gap d.  The sign of a coefficient
   is the majority direction in that layer; |A_d| is the layer imbalance.  A
   fixed coefficient-magnitude vector is therefore a slice through the tiling
   hypercube.

3. Fixed-base-path quotient:
   force every gap-1 edge to agree with the Hamiltonian path.  This is the
   repo tiling gauge.  The coefficient polynomial only sees diagonal counts,
   while H(T) and SCC/fiber data vary inside the same coefficient profile.

The output is intentionally exploratory.  It is a reproducible idea generator,
not a proof of Bunyakovsky or an irreducibility criterion.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations, product
from math import comb

import sympy as sp


x = sp.Symbol("x")


def is_irreducible(poly: sp.Poly) -> bool:
    if poly.degree() <= 0:
        return False
    _unit, factors = sp.factor_list(poly.as_expr(), x)
    return sum(exp for _fac, exp in factors) == 1


def factor_degrees(poly: sp.Poly) -> list[int]:
    if poly.degree() <= 0:
        return []
    _unit, factors = sp.factor_list(poly.as_expr(), x)
    out: list[int] = []
    for fac, exp in factors:
        out.extend([sp.Poly(fac, x).degree()] * exp)
    return out


def fixed_divisor(poly: sp.Poly) -> int:
    if poly.degree() < 0:
        return 0
    vals = [abs(int(poly.eval(k))) for k in range(poly.degree() + 1)]
    g = 0
    for val in vals:
        g = sp.igcd(g, val)
    return int(g)


def prime_hits(poly: sp.Poly, limit: int = 80) -> int:
    hits = 0
    for t in range(limit + 1):
        val = int(poly.eval(t))
        if val > 1 and sp.isprime(val):
            hits += 1
    return hits


def caps_for_N(N: int) -> dict[int, int]:
    return {d: N - d for d in range(1, N)}


def profile_weight(profile: dict[int, int], fixed_path: bool = False) -> int:
    weight = 1
    for d, c in profile.items():
        cap = max(profile.keys()) + 1 - d
        if fixed_path and d == 1:
            continue
        weight *= comb(cap, c)
    return weight


def all_profiles(N: int, fixed_path: bool = False) -> list[dict[int, int]]:
    caps = caps_for_N(N)
    ranges = []
    ds = list(range(1, N))
    for d in ds:
        if fixed_path and d == 1:
            ranges.append([caps[d]])
        else:
            ranges.append(range(caps[d] + 1))
    return [dict(zip(ds, counts)) for counts in product(*ranges)]


def count_poly(profile: dict[int, int], constant: int = 1) -> sp.Poly:
    expr = constant
    for d, c in profile.items():
        expr += c * x**d
    return sp.Poly(expr, x, domain=sp.ZZ)


def fixed_path_count_poly(profile: dict[int, int], N: int) -> sp.Poly:
    # Treat the Hamiltonian path layer as the constant term, matching the
    # user's "a_0 is the fixed path" suggestion.
    expr = N - 1
    for d in range(2, N):
        expr += profile[d] * x ** (d - 1)
    return sp.Poly(expr, x, domain=sp.ZZ)


def centered_poly(profile: dict[int, int], constant: int = 1) -> sp.Poly:
    N = max(profile.keys()) + 1
    expr = constant
    for d, c in profile.items():
        cap = N - d
        expr += (2 * c - cap) * x**d
    return sp.Poly(expr, x, domain=sp.ZZ)


def profile_from_centered_magnitudes(N: int, mags: dict[int, int], signs: dict[int, int]) -> dict[int, int] | None:
    profile: dict[int, int] = {}
    for d, cap in caps_for_N(N).items():
        mag = mags[d]
        val = signs[d] * mag
        if abs(val) > cap or (val + cap) % 2:
            return None
        profile[d] = (val + cap) // 2
    return profile


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


def fixed_path_tilings(N: int):
    variable_pairs = [(i, j) for i in range(N) for j in range(i + 2, N)]
    total = 1 << len(variable_pairs)
    for mask in range(total):
        adj = [[False] * N for _ in range(N)]
        for i in range(N - 1):
            adj[i][i + 1] = True
        profile = {d: 0 for d in range(1, N)}
        profile[1] = N - 1
        for bit, (i, j) in enumerate(variable_pairs):
            d = j - i
            if mask & (1 << bit):
                adj[i][j] = True
                profile[d] += 1
            else:
                adj[j][i] = True
        yield profile, adj


def profile_summary(N: int, fixed_path: bool = False) -> dict[str, object]:
    profiles = all_profiles(N, fixed_path=fixed_path)
    base = N
    total_weight = 0
    count_irred = count_prime_value = count_prime_mismatch = 0
    centered_irred = centered_admissible = centered_prime_hits = 0
    weighted_centered_irred = 0
    extreme_total = extreme_centered_irred = 0

    for profile in profiles:
        w = profile_weight(profile, fixed_path=fixed_path)
        total_weight += w
        poly = fixed_path_count_poly(profile, N) if fixed_path else count_poly(profile)
        val = abs(int(poly.eval(base)))
        irred = is_irreducible(poly)
        if irred:
            count_irred += 1
        if poly.degree() > 0 and val > 1 and sp.isprime(val):
            count_prime_value += 1
            if not irred:
                count_prime_mismatch += 1

        cpoly = centered_poly(profile)
        cirred = is_irreducible(cpoly)
        if cirred:
            centered_irred += 1
            weighted_centered_irred += w
            if fixed_divisor(cpoly) == 1:
                centered_admissible += 1
        if prime_hits(cpoly, limit=50):
            centered_prime_hits += 1

        if all(profile[d] in {0, N - d} for d in range(1, N)):
            extreme_total += 1
            if cirred:
                extreme_centered_irred += 1

    return {
        "N": N,
        "fixed_path": fixed_path,
        "profiles": len(profiles),
        "weighted_tilings": total_weight,
        "count_poly_irred_profiles": count_irred,
        "cohn_prime_value_profiles": count_prime_value,
        "cohn_prime_mismatches": count_prime_mismatch,
        "centered_irred_profiles": centered_irred,
        "centered_admissible_profiles": centered_admissible,
        "centered_has_prime_hit_profiles": centered_prime_hits,
        "weighted_centered_irred_tilings": weighted_centered_irred,
        "extreme_profiles": extreme_total,
        "extreme_centered_irred": extreme_centered_irred,
    }


def fixed_path_fiber_H_summary(N: int = 6) -> dict[str, object]:
    by_profile: dict[tuple[int, ...], list[int]] = defaultdict(list)
    buckets: dict[str, list[int]] = defaultdict(list)
    examples: dict[str, object] = {}

    for profile, adj in fixed_path_tilings(N):
        key = tuple(profile[d] for d in range(1, N))
        H = hamiltonian_paths(adj)
        by_profile[key].append(H)

        cpoly = centered_poly(profile)
        path_poly = fixed_path_count_poly(profile, N)
        cohn_prime = path_poly.degree() > 0 and sp.isprime(abs(int(path_poly.eval(N))))
        centered_irred = is_irreducible(cpoly)
        buckets[f"cohn_prime={cohn_prime}"].append(H)
        buckets[f"centered_irred={centered_irred}"].append(H)
        if cohn_prime and "cohn_prime" not in examples:
            examples["cohn_prime"] = {
                "profile": key,
                "path_poly": sp.sstr(path_poly.as_expr()),
                "value_base_N": int(path_poly.eval(N)),
                "H": H,
            }
        if centered_irred and "centered_irred" not in examples:
            examples["centered_irred"] = {
                "profile": key,
                "centered_poly": sp.sstr(cpoly.as_expr()),
                "factor_degrees": factor_degrees(cpoly),
                "H": H,
            }

    mixed = 0
    max_spread = 0
    max_spread_key: tuple[int, ...] | None = None
    for key, vals in by_profile.items():
        spread = max(vals) - min(vals)
        if spread:
            mixed += 1
            if spread > max_spread:
                max_spread = spread
                max_spread_key = key

    bucket_stats = {}
    for label, vals in buckets.items():
        bucket_stats[label] = {
            "count": len(vals),
            "avg_H": round(sum(vals) / len(vals), 3) if vals else None,
            "min_H": min(vals) if vals else None,
            "max_H": max(vals) if vals else None,
        }

    return {
        "N": N,
        "tilings": sum(len(v) for v in by_profile.values()),
        "profiles": len(by_profile),
        "profiles_with_H_variation": mixed,
        "max_H_spread_in_one_profile": max_spread,
        "max_spread_profile": max_spread_key,
        "bucket_stats": bucket_stats,
        "examples": examples,
    }


def magnitude_slice_summary(N: int = 6) -> list[dict[str, object]]:
    caps = caps_for_N(N)
    slices = {
        "max_layer_magnitudes": {d: caps[d] for d in caps},
        "minimum_parity_magnitudes": {d: caps[d] % 2 for d in caps},
        "alternating_rigid_soft": {d: caps[d] if d % 2 == 0 else caps[d] % 2 for d in caps},
        "apex_rigid_base_soft": {d: (caps[d] if d == N - 1 else caps[d] % 2) for d in caps},
    }
    out: list[dict[str, object]] = []
    for name, mags in slices.items():
        valid = irred = admissible = prime_hit = 0
        polys_seen: set[str] = set()
        examples: list[str] = []
        for sign_bits in product([-1, 1], repeat=N - 1):
            signs = {d: sign_bits[d - 1] for d in range(1, N)}
            profile = profile_from_centered_magnitudes(N, mags, signs)
            if profile is None:
                continue
            valid += 1
            poly = centered_poly(profile)
            poly_s = sp.sstr(poly.as_expr())
            polys_seen.add(poly_s)
            if is_irreducible(poly):
                irred += 1
                if fixed_divisor(poly) == 1:
                    admissible += 1
                if len(examples) < 2 and poly_s not in examples:
                    examples.append(poly_s)
            if prime_hits(poly, limit=50):
                prime_hit += 1
        out.append(
            {
                "slice": name,
                "magnitudes": tuple(mags[d] for d in range(1, N)),
                "valid_sign_assignments": valid,
                "distinct_polys": len(polys_seen),
                "irreducible": irred,
                "admissible": admissible,
                "has_prime_hit": prime_hit,
                "examples": examples,
            }
        )
    return out


@dataclass(frozen=True)
class Idea:
    name: str
    novelty: int
    testability: int
    prime_poly_bridge: int
    tiling_bridge: int
    support_transfer: int
    risk: int
    note: str

    def criteria(self) -> tuple[int, int, int, int, int, int]:
        return (
            self.novelty,
            self.testability,
            self.prime_poly_bridge,
            self.tiling_bridge,
            self.support_transfer,
            6 - self.risk,
        )


IDEAS = [
    Idea("diagonal_count_cohn_map", 4, 5, 5, 5, 3, 1, "layer counts are base-N digits; prime profile integer certifies irreducibility"),
    Idea("centered_magnetization_slice", 5, 5, 4, 5, 4, 2, "coefficient magnitudes are layer-imbalance slices through tiling space"),
    Idea("fiber_entropy_vs_irreducibility", 4, 4, 3, 5, 5, 2, "same polynomial profile can hide many H-values and SCC shapes"),
    Idea("root_argument_tournament", 5, 3, 4, 2, 3, 3, "orient roots by argument/real part; reducible factors become conjugation-stable modules"),
    Idea("newton_polygon_edge_tournament", 4, 4, 4, 3, 4, 2, "vertices are Newton edges; orientations by slope/support dominance"),
    Idea("finite_field_factor_race", 4, 5, 5, 2, 3, 2, "orient primes by factorization type of f mod p and edge flips across p"),
    Idea("evaluation_time_tournament", 3, 5, 4, 3, 3, 2, "vertices are n-values; edge by which value gives stronger Singh certificate"),
    Idea("derivative_trace_tournament", 4, 4, 4, 3, 4, 2, "Singh derivative conditions become a hierarchy of trace/support gates"),
    Idea("galois_orbit_recombination", 5, 2, 5, 2, 4, 4, "subsets of roots oriented by stabilizer size and trace integrality"),
    Idea("code_support_polynomial_lift", 4, 3, 3, 4, 5, 3, "length-72 support moves become coefficient-layer tilings with fixed low terms"),
]


def idea_tournament() -> dict[str, object]:
    n = len(IDEAS)
    wins = {i: set() for i in range(n)}
    cycles = 0
    edge_flips_vs_novelty = 0
    for i, j in combinations(range(n), 2):
        a, b = IDEAS[i], IDEAS[j]
        a_votes = sum(x > y for x, y in zip(a.criteria(), b.criteria()))
        b_votes = sum(y > x for x, y in zip(a.criteria(), b.criteria()))
        if a_votes > b_votes:
            wins[i].add(j)
            if a.novelty < b.novelty:
                edge_flips_vs_novelty += 1
        elif b_votes > a_votes:
            wins[j].add(i)
            if b.novelty < a.novelty:
                edge_flips_vs_novelty += 1
        else:
            # Tie Hamiltonian path: prefer the lower-risk idea, then declaration order.
            if (a.risk, i) <= (b.risk, j):
                wins[i].add(j)
                if a.novelty < b.novelty:
                    edge_flips_vs_novelty += 1
            else:
                wins[j].add(i)
                if b.novelty < a.novelty:
                    edge_flips_vs_novelty += 1
    for i, j, k in combinations(range(n), 3):
        if j in wins[i] and k in wins[j] and i in wins[k]:
            cycles += 1
        if k in wins[i] and j in wins[k] and i in wins[j]:
            cycles += 1
    scores = {IDEAS[i].name: len(wins[i]) for i in range(n)}
    return {
        "score_hist": dict(sorted(Counter(scores.values()).items())),
        "directed_3cycles": cycles,
        "edge_flips_vs_novelty_only": edge_flips_vs_novelty,
        "ranking": [name for name, _ in sorted(scores.items(), key=lambda kv: (-kv[1], kv[0]))],
    }


def print_profile_summaries() -> None:
    print("DIAGONAL PROFILE SUMMARIES")
    print("==========================")
    for N in range(3, 8):
        for fixed in (False, True):
            summary = profile_summary(N, fixed_path=fixed)
            print(summary)
        print()


def print_fiber_summary() -> None:
    print("FIXED-PATH FIBER H SUMMARY")
    print("==========================")
    summary = fixed_path_fiber_H_summary(6)
    print(summary)


def print_magnitude_slices() -> None:
    print("COEFFICIENT MAGNITUDE SLICES (N=6)")
    print("==================================")
    for row in magnitude_slice_summary(6):
        print(row)


def print_generated_ideas() -> None:
    print("PROCEDURALLY GENERATED HIDING PLACES")
    print("====================================")
    for idea in IDEAS:
        print(f"{idea.name}: criteria={idea.criteria()} :: {idea.note}")
    print(f"idea_tournament={idea_tournament()}")


def print_assumption_challenge() -> None:
    print("ASSUMPTION CHALLENGE")
    print("====================")
    print("Considered vertex sets: arcs, coefficient layers, layer profiles,")
    print("fixed coefficient-magnitude slices, roots, factor subsets, finite primes,")
    print("evaluation times, Newton polygon edges, derivative obligations, code")
    print("support moves, and LRC denominator resources.")
    print("Chosen primary quotient: diagonal profiles of fixed-path tilings.")
    print("Preserved: the user's triangular coefficient-layer geometry, Cohn prime")
    print("certificates, sign/magnitude slices, and the fact that many tournaments can")
    print("share one coefficient profile.")
    print("Destroyed: individual arc arrangement, most H/SCC data, and actual")
    print("asymptotic Bunyakovsky information.  The fiber-H audit measures that loss.")


def main() -> None:
    print("COEFFICIENT TILING PRIME BRIDGE ATLAS")
    print("=====================================")
    print("Layer d has N-d arcs.  Counts c_d and magnetizations A_d=2c_d-(N-d)")
    print("turn tournament tilings into integer polynomials.\n")
    print_profile_summaries()
    print_fiber_summary()
    print()
    print_magnitude_slices()
    print()
    print_generated_ideas()
    print()
    print_assumption_challenge()


if __name__ == "__main__":
    main()
