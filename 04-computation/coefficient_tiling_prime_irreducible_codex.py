#!/usr/bin/env python3
"""Coefficient tiling tournaments for prime/irreducible polynomial carriers.

codex-2026-06-12

This is a procedural addendum to HYP-2447/HYP-2448.  The user's suggested
picture is:

  degree 5:
    a5 on the apex tile,
    a4 on the two next tiles,
    a3 on the next three,
    a2 on the next four,
    a1 on the next five,
    a0 as the fixed Hamiltonian-path / constant spine.

That row count is exactly the skip-row decomposition of a fixed ordered
tournament.  This script tests two literal encodings:

  top_no_constant:
    degree d polynomial -> tournament on d+1 vertices.
    coefficient a_s orients every edge of skip s, for s=1..d.
    The constant a0 is an external scalar.

  constant_spine:
    degree d polynomial -> tournament on d+2 vertices.
    coefficient a_0 orients the adjacent Hamiltonian-path row, and a_d is
    the apex skip d+1 row.  This is the "turn the idea on its head" version:
    the constant term is not outside the tournament; it is the spine.

The point is not to prove Bunyakovsky.  The point is to discover what the
coefficient-sign tournament preserves and what side channels it destroys:
fixed divisors, value-prime certificates, Cohn digit addresses, and
recombination ambiguity.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache, reduce
from itertools import permutations, product
from math import gcd

import sympy as sp


x = sp.Symbol("x")


def sign(n: int) -> int:
    return 1 if n >= 0 else -1


def igcd(values: list[int]) -> int:
    return abs(reduce(gcd, [abs(v) for v in values], 0))


def fixed_divisor(poly: sp.Poly) -> int:
    d = poly.degree()
    return igcd([int(poly.eval(k)) for k in range(d + 1)])


def factor_degrees(poly: sp.Poly) -> list[int]:
    _unit, factors = sp.factor_list(poly.as_expr(), x)
    degrees: list[int] = []
    for fac, exp in factors:
        degrees.extend([sp.Poly(fac, x).degree()] * exp)
    return sorted(degrees)


def is_irreducible(poly: sp.Poly) -> bool:
    return len(factor_degrees(poly)) == 1


def prime_hits(poly: sp.Poly, limit: int = 50) -> list[tuple[int, int]]:
    hits: list[tuple[int, int]] = []
    for n in range(limit + 1):
        val = int(poly.eval(n))
        if val > 1 and sp.isprime(val):
            hits.append((n, val))
    return hits


def local_zero_primes(poly: sp.Poly, primes: tuple[int, ...] = (2, 3, 5, 7)) -> tuple[int, ...]:
    """Primes p for which f(r)=0 mod p for every residue r mod p."""
    bad: list[int] = []
    for p in primes:
        if all(int(poly.eval(r)) % p == 0 for r in range(p)):
            bad.append(p)
    return tuple(bad)


def poly_from_coeffs(coeffs: tuple[int, ...]) -> sp.Poly:
    expr = sum(c * x**i for i, c in enumerate(coeffs))
    return sp.Poly(expr, x, domain=sp.ZZ)


def adjacency_from_skip_signs(skip_signs: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    """Tournament adjacency for vertices 0..N-1.

    skip_signs[s-1] controls all pairs i<j with j-i=s.  Positive means
    i -> j, negative means j -> i.
    """
    n = len(skip_signs) + 1
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            s = j - i
            if skip_signs[s - 1] > 0:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
    return tuple(tuple(row) for row in adj)


def coefficient_tournament(coeffs: tuple[int, ...], mode: str) -> tuple[tuple[int, ...], ...]:
    d = len(coeffs) - 1
    if mode == "top_no_constant":
        # a_s controls skip s, for s=1..d.  a_d is the apex.
        skip_signs = tuple(sign(coeffs[s]) for s in range(1, d + 1))
    elif mode == "constant_spine":
        # a_{s-1} controls skip s, for s=1..d+1.  a_0 is the base path.
        skip_signs = tuple(sign(coeffs[s - 1]) for s in range(1, d + 2))
    else:
        raise ValueError(mode)
    return adjacency_from_skip_signs(skip_signs)


def bitstring(adj: tuple[tuple[int, ...], ...], perm: tuple[int, ...] | None = None) -> str:
    n = len(adj)
    if perm is None:
        perm = tuple(range(n))
    bits: list[str] = []
    for a in range(n):
        for b in range(a + 1, n):
            i, j = perm[a], perm[b]
            bits.append("1" if adj[i][j] else "0")
    return "".join(bits)


@lru_cache(maxsize=None)
def canonical_key_from_bits(bits: str, n: int) -> str:
    adj = [[0] * n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits[k] == "1":
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            k += 1
    tadj = tuple(tuple(row) for row in adj)
    return min(bitstring(tadj, p) for p in permutations(range(n)))


def canonical_key(adj: tuple[tuple[int, ...], ...]) -> str:
    return canonical_key_from_bits(bitstring(adj), len(adj))


def score_sequence(adj: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    return tuple(sorted(sum(row) for row in adj))


def directed_3cycles(adj: tuple[tuple[int, ...], ...]) -> int:
    n = len(adj)
    c = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                # A triple is cyclic iff every vertex has outdegree 1 inside it.
                out = []
                for v in (i, j, k):
                    out.append(sum(adj[v][u] for u in (i, j, k) if u != v))
                if out == [1, 1, 1]:
                    c += 1
    return c


def hamiltonian_paths(adj: tuple[tuple[int, ...], ...]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += val
    return sum(dp[(1 << n) - 1])


def row_geometry(degree: int, mode: str) -> list[tuple[str, int, int]]:
    """Return (coefficient, skip, row_size), printed high-to-low."""
    if mode == "top_no_constant":
        n_vertices = degree + 1
        rows = [(f"a{s}", s, n_vertices - s) for s in range(1, degree + 1)]
    elif mode == "constant_spine":
        n_vertices = degree + 2
        rows = [(f"a{s-1}", s, n_vertices - s) for s in range(1, degree + 2)]
    else:
        raise ValueError(mode)
    return list(reversed(rows))


@dataclass(frozen=True)
class PolyRecord:
    coeffs: tuple[int, ...]
    poly: sp.Poly
    signs: tuple[int, ...]
    top_key: str
    spine_key: str
    top_score: tuple[int, ...]
    spine_score: tuple[int, ...]
    irreducible: bool
    factor_degrees: tuple[int, ...]
    fixed_divisor: int
    local_zero_primes: tuple[int, ...]
    prime_hit_count: int


def build_degree4_sweep() -> list[PolyRecord]:
    """Degree-4 sweep over coefficient magnitudes 1..3 and leading sign +."""
    rows: list[PolyRecord] = []
    # coeffs are a0..a4.  Keep leading sign positive; multiply-by--1 is a
    # global gauge for irreducibility/fixed divisor, but not for positive
    # prime values, so we keep the ordinary forward orientation.
    for mags in product((1, 2, 3), repeat=5):
        for signs_low in product((-1, 1), repeat=4):
            coeffs = tuple(signs_low[i] * mags[i] for i in range(4)) + (mags[4],)
            poly = poly_from_coeffs(coeffs)
            top_adj = coefficient_tournament(coeffs, "top_no_constant")
            spine_adj = coefficient_tournament(coeffs, "constant_spine")
            rows.append(
                PolyRecord(
                    coeffs=coeffs,
                    poly=poly,
                    signs=tuple(sign(c) for c in coeffs),
                    top_key=canonical_key(top_adj),
                    spine_key=canonical_key(spine_adj),
                    top_score=score_sequence(top_adj),
                    spine_score=score_sequence(spine_adj),
                    irreducible=is_irreducible(poly),
                    factor_degrees=tuple(factor_degrees(poly)),
                    fixed_divisor=fixed_divisor(poly),
                    local_zero_primes=local_zero_primes(poly),
                    prime_hit_count=len(prime_hits(poly, 40)),
                )
            )
    return rows


def mixed_summary(rows: list[PolyRecord], key_name: str, label_name: str) -> dict[str, object]:
    buckets: dict[object, list[PolyRecord]] = defaultdict(list)
    for row in rows:
        key = getattr(row, key_name)
        buckets[key].append(row)
    mixed = []
    for key, vals in buckets.items():
        labels = {getattr(v, label_name) for v in vals}
        if len(labels) > 1:
            mixed.append((key, vals, labels))
    max_bucket = max(len(v) for v in buckets.values()) if buckets else 0
    example = None
    if mixed:
        key, vals, labels = mixed[0]
        seen: dict[object, PolyRecord] = {}
        for v in vals:
            lab = getattr(v, label_name)
            if lab not in seen:
                seen[lab] = v
        example = {
            "key": str(key)[:80],
            "labels": sorted(str(l) for l in labels),
            "examples": [
                {
                    "label": str(lab),
                    "coeffs": rec.coeffs,
                    "poly": str(rec.poly.as_expr()),
                    "factors": rec.factor_degrees,
                    "fixed_divisor": rec.fixed_divisor,
                    "local_zero_primes": rec.local_zero_primes,
                }
                for lab, rec in list(seen.items())[:3]
            ],
        }
    return {
        "groups": len(buckets),
        "mixed": len(mixed),
        "max_bucket": max_bucket,
        "example": example,
    }


def mixed_by_custom_key(rows: list[PolyRecord], key_func, label_func) -> dict[str, object]:
    buckets: dict[object, list[PolyRecord]] = defaultdict(list)
    for row in rows:
        buckets[key_func(row)].append(row)
    mixed = []
    for key, vals in buckets.items():
        labels = {label_func(v) for v in vals}
        if len(labels) > 1:
            mixed.append((key, vals, labels))
    return {
        "groups": len(buckets),
        "mixed": len(mixed),
        "max_bucket": max(len(v) for v in buckets.values()) if buckets else 0,
    }


def find_same_sign_examples(rows: list[PolyRecord]) -> dict[str, object]:
    buckets: dict[tuple[int, ...], list[PolyRecord]] = defaultdict(list)
    for row in rows:
        buckets[row.signs].append(row)
    out: dict[str, object] = {}
    for signs, vals in buckets.items():
        irreducible_vals = [v for v in vals if v.irreducible]
        reducible_vals = [v for v in vals if not v.irreducible]
        if irreducible_vals and reducible_vals and "irreducibility" not in out:
            out["irreducibility"] = (signs, irreducible_vals[0], reducible_vals[0])
        admissible_vals = [v for v in vals if v.fixed_divisor == 1]
        blocked_vals = [v for v in vals if v.fixed_divisor > 1]
        if admissible_vals and blocked_vals and "fixed_divisor" not in out:
            out["fixed_divisor"] = (signs, admissible_vals[0], blocked_vals[0])
        prime_rich = [v for v in vals if v.prime_hit_count >= 5]
        prime_poor = [v for v in vals if v.prime_hit_count == 0]
        if prime_rich and prime_poor and "prime_hits" not in out:
            out["prime_hits"] = (signs, prime_rich[0], prime_poor[0])
    return out


def centered_digit_signs(n: int, base: int) -> tuple[int, ...]:
    digits: list[int] = []
    m = n
    while m:
        digits.append(m % base)
        m //= base
    if not digits:
        digits = [0]
    center2 = base - 1
    return tuple(1 if 2 * d >= center2 else -1 for d in digits)


def cohn_rows() -> list[dict[str, object]]:
    samples = [
        (101, 10),
        (131, 10),
        (2339, 10),
        (9841, 3),
        (2047, 2),
        (6551, 10),
        (7643, 10),
    ]
    rows = []
    for n, base in samples:
        digits: list[int] = []
        m = n
        while m:
            digits.append(m % base)
            m //= base
        expr = sum(d * x**i for i, d in enumerate(digits))
        poly = sp.Poly(expr, x, domain=sp.ZZ)
        signs = centered_digit_signs(n, base)
        adj = adjacency_from_skip_signs(signs)
        rows.append(
            {
                "n": n,
                "base": base,
                "digits_low_to_high": digits,
                "centered_signs": signs,
                "poly": str(poly.as_expr()),
                "omega_n": sum(sp.factorint(n).values()),
                "factor_degrees": factor_degrees(poly),
                "irreducible": is_irreducible(poly),
                "c3": directed_3cycles(adj),
                "H": hamiltonian_paths(adj),
                "score": score_sequence(adj),
            }
        )
    return rows


@dataclass(frozen=True)
class Carrier:
    name: str
    preserves_row_address: int
    detects_fixed_divisor: int
    sees_reverse_certificate: int
    recombination_power: int
    algorithmic_cost: int  # lower is better
    wildness: int
    note: str


def carrier_tournament() -> tuple[list[Carrier], dict[str, object]]:
    carriers = [
        Carrier("raw_unmarked_skip_tournament", 1, 1, 1, 1, 1, 2, "beautiful but loses coefficient labels"),
        Carrier("marked_skip_row_tiling", 4, 2, 2, 2, 1, 4, "user tiling with coefficient row labels retained"),
        Carrier("constant_spine_tiling", 5, 3, 2, 2, 1, 5, "a0 becomes Hamiltonian spine rather than external scalar"),
        Carrier("row_sign_plus_local_residues", 5, 5, 3, 2, 2, 4, "adds fixed-divisor/local obstruction address"),
        Carrier("cohn_digit_carry_address", 4, 2, 5, 2, 2, 4, "base-b evaluation weights rows exponentially"),
        Carrier("singh_value_depth_state", 3, 3, 5, 3, 3, 3, "evaluation depth gives reverse irreducibility certificate"),
        Carrier("trace_subset_recombination", 3, 2, 3, 5, 4, 5, "loose atoms first, integer factor later"),
        Carrier("newton_polygon_valuation_tiling", 4, 5, 3, 4, 2, 5, "edge slopes from p-adic coefficient valuations"),
        Carrier("coefficient_vertex_tournament", 2, 1, 2, 1, 1, 3, "vertices are coefficients, usually too coarse"),
    ]
    n = len(carriers)
    adj = [[0] * n for _ in range(n)]
    # Majority over five criteria.  For algorithmic_cost, lower wins.
    for i, a in enumerate(carriers):
        for j, b in enumerate(carriers):
            if i == j:
                continue
            votes = 0
            votes += a.preserves_row_address > b.preserves_row_address
            votes += a.detects_fixed_divisor > b.detects_fixed_divisor
            votes += a.sees_reverse_certificate > b.sees_reverse_certificate
            votes += a.recombination_power > b.recombination_power
            votes += a.algorithmic_cost < b.algorithmic_cost
            if votes > 2:
                adj[i][j] = 1
            elif votes < 2:
                adj[j][i] = 1
            else:
                # Tie path: prefer wildness, then name for determinism.
                if (a.wildness, a.name) >= (b.wildness, b.name):
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1
    tadj = tuple(tuple(row) for row in adj)
    ranking = sorted(range(n), key=lambda i: sum(adj[i]), reverse=True)
    c3 = directed_3cycles(tadj)
    h = hamiltonian_paths(tadj)
    return carriers, {
        "score_hist": dict(sorted(Counter(sum(row) for row in adj).items())),
        "directed_3cycles": c3,
        "hamiltonian_paths": h,
        "ranking": [carriers[i].name for i in ranking],
    }


def fmt_record(rec: PolyRecord) -> str:
    return (
        f"coeffs={rec.coeffs}; f={rec.poly.as_expr()}; "
        f"factors={rec.factor_degrees}; fd={rec.fixed_divisor}; "
        f"local0={rec.local_zero_primes}; prime_hits={rec.prime_hit_count}"
    )


def main() -> None:
    print("=" * 78)
    print("Coefficient tiling tournaments for prime <-> irreducible polynomials")
    print("=" * 78)
    print("HYP-2449 / T793 / OPEN-Q-071 candidate")
    print()

    print("[1] Degree-5 row geometry")
    for mode in ("top_no_constant", "constant_spine"):
        print(f"  mode={mode}")
        for coeff, skip, size in row_geometry(5, mode):
            label = "apex" if size == 1 else ("spine/HP row" if skip == 1 else "row")
            print(f"    {coeff:>2} -> skip {skip}, row_size {size:2d} ({label})")
    print()

    print("[2] Degree-5 sign-row fingerprints for a concrete polynomial")
    coeffs5 = (1, -2, 3, -4, 5, 6)  # a0..a5
    for mode in ("top_no_constant", "constant_spine"):
        adj = coefficient_tournament(coeffs5, mode)
        print(
            f"  {mode}: vertices={len(adj)}, score={score_sequence(adj)}, "
            f"c3={directed_3cycles(adj)}, H={hamiltonian_paths(adj)}, "
            f"canonical={canonical_key(adj)[:24]}..."
        )
    print()

    print("[3] Exhaustive degree-4 sign+magnitude sweep")
    rows = build_degree4_sweep()
    print(f"  rows={len(rows)} (degree 4, |a_i| in {{1,2,3}}, leading sign positive)")
    print(f"  irreducible={sum(r.irreducible for r in rows)}, reducible={sum(not r.irreducible for r in rows)}")
    print(f"  fixed_divisor>1={sum(r.fixed_divisor > 1 for r in rows)}")
    print(f"  no_prime_hit_0..40={sum(r.prime_hit_count == 0 for r in rows)}")
    print()

    print("[4] What leaks under different coefficient-tournament quotients?")
    summaries = [
        ("top_key", "irreducible"),
        ("spine_key", "irreducible"),
        ("top_score", "irreducible"),
        ("spine_score", "irreducible"),
        ("top_key", "fixed_divisor"),
        ("spine_key", "fixed_divisor"),
    ]
    for key_name, label_name in summaries:
        s = mixed_summary(rows, key_name, label_name)
        print(
            f"  key={key_name:10s} label={label_name:13s} "
            f"groups={s['groups']:3d} mixed={s['mixed']:3d} max_bucket={s['max_bucket']:4d}"
        )
        if s["example"] and label_name == "irreducible" and key_name in {"top_key", "spine_key"}:
            print(f"    example mixed key={s['example']['key']}")
            for ex in s["example"]["examples"]:
                print(f"      {ex['label']}: coeffs={ex['coeffs']} f={ex['poly']} factors={ex['factors']}")
    custom = [
        (
            "marked_signs",
            lambda r: r.signs,
            lambda r: r.irreducible,
        ),
        (
            "marked_signs+local_zero_primes",
            lambda r: (r.signs, r.local_zero_primes),
            lambda r: r.fixed_divisor,
        ),
        (
            "marked_signs+factor_degrees",
            lambda r: (r.signs, r.factor_degrees),
            lambda r: r.irreducible,
        ),
    ]
    for name, key_func, label_func in custom:
        s = mixed_by_custom_key(rows, key_func, label_func)
        print(
            f"  key={name:31s} groups={s['groups']:3d} mixed={s['mixed']:3d} max_bucket={s['max_bucket']:4d}"
        )
    print()

    print("[5] Same sign tournament, different arithmetic")
    examples = find_same_sign_examples(rows)
    for label, payload in examples.items():
        signs, good, bad = payload
        print(f"  {label}: signs={signs}")
        print(f"    A: {fmt_record(good)}")
        print(f"    B: {fmt_record(bad)}")
    print()

    print("[6] Cohn digit rows show why sign alone is too weak")
    for row in cohn_rows():
        print(
            f"  N={row['n']:5d} base={row['base']:2d} omega={row['omega_n']} "
            f"digits={row['digits_low_to_high']} signs={row['centered_signs']}"
        )
        print(
            f"    poly={row['poly']}; factor_degrees={row['factor_degrees']}; "
            f"irreducible={row['irreducible']}; sign-tournament score={row['score']}, "
            f"c3={row['c3']}, H={row['H']}"
        )
    print()

    print("[7] Procedurally generated carrier tournament")
    carriers, fp = carrier_tournament()
    print(f"  score_hist={fp['score_hist']}")
    print(f"  directed_3cycles={fp['directed_3cycles']}")
    print(f"  Hamiltonian_paths={fp['hamiltonian_paths']}")
    print("  ranking:")
    for name in fp["ranking"]:
        c = next(c for c in carriers if c.name == name)
        print(
            f"    {name:34s} "
            f"(row={c.preserves_row_address}, fixed={c.detects_fixed_divisor}, "
            f"rev={c.sees_reverse_certificate}, recomb={c.recombination_power}, "
            f"cost={c.algorithmic_cost})"
        )
        print(f"      {c.note}")
    print()

    print("[8] New directions extracted")
    directions = [
        "Treat coefficient rows as modules: polynomial multiplication becomes convolution of skip-row modules, not just scalar coefficient multiplication.",
        "Read Gauss lemma as switching: content is a global gauge; primitive part is the quotient where row gcd debt has been removed.",
        "Read Eisenstein as a valuation tournament on coefficient rows: all lower rows point into the p-adic sink except the leading escape and p^2 constant guard.",
        "Read Newton polygons as coefficient-tile tournaments whose edges are slope comparisons between p-adic row valuations.",
        "Read Cohn as exponential row weighting: base b gives the skip rows a place-value address, and primality says no nontrivial row subset recombines at that address.",
        "For LRC14, copy fixed-divisor detection: a denominator family is dangerous only if some runner/resource row vanishes on every residue in a finite local gate.",
        "For A000568/switching classes, coefficient signs are tiling-cube coordinates; irreducibility should live in marked switching classes, not unmarked tournaments.",
    ]
    for i, text in enumerate(directions, 1):
        print(f"  {i}. {text}")
    print()

    print("Conclusion:")
    print("  The coefficient-sign tournament is real, but unmarked it is too lossy.")
    print("  The useful carrier is a marked skip-row tiling plus local residue/valuation")
    print("  and evaluation-depth addresses.  That is exactly the prime<->irreducible")
    print("  bridge in tournament clothing: signs give the tiling cube, while primes")
    print("  and irreducibility live in the addresses that survive recombination.")


if __name__ == "__main__":
    main()
