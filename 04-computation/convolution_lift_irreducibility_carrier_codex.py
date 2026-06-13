#!/usr/bin/env python3
"""Convolution-lift carriers for integer polynomial irreducibility.

codex-2026-06-12

This continues HYP-2449.  The coefficient-sign tiling is a real fixed-path
tournament carrier, but reducibility itself is not a sign property.  It is the
existence of a hidden two-dimensional lift:

    a_k = sum_{i+j=k} b_i c_j.

So a one-dimensional coefficient row is the diagonal shadow of a coefficient
grid.  Irreducibility means no nontrivial integer grid lift exists.

The high-leverage move is to test local shadows of that lift:

* residue/convolution blockers: if f mod p has no nontrivial diagonal lift,
  then a primitive f is irreducible over Z;
* valuation/Newton blockers: if the p-adic lower hull has one primitive edge,
  the residue picture may look reducible while the valuation lift is rigid.

Tournament Analysis here is not over primes-as-numbers.  Vertices are local
certificate channels (small primes and the unresolved channel).  The pairwise
observable is certificate strength on the degree-4 coefficient-row scout.
Ties are resolved by the fixed prime order, the analogue of an arbitrary fixed
Hamiltonian path.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import reduce
from itertools import product
from math import gcd
import warnings

import sympy as sp
from sympy.utilities.exceptions import SymPyDeprecationWarning


warnings.filterwarnings("ignore", category=SymPyDeprecationWarning)

x = sp.Symbol("x")
SMALL_PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31)


def poly_from_coeffs(coeffs: tuple[int, ...]) -> sp.Poly:
    return sp.Poly(sum(c * x**i for i, c in enumerate(coeffs)), x, domain=sp.ZZ)


def primitive_part(poly: sp.Poly) -> sp.Poly:
    _content, prim = sp.primitive(poly.as_expr(), expand=True)
    return sp.Poly(prim, x, domain=sp.ZZ)


def factor_degrees_zz(poly: sp.Poly) -> tuple[int, ...]:
    _unit, factors = sp.factor_list(poly.as_expr(), x)
    out: list[int] = []
    for fac, exp in factors:
        out.extend([sp.Poly(fac, x).degree()] * exp)
    return tuple(sorted(out))


def is_irreducible_zz(poly: sp.Poly) -> bool:
    return len(factor_degrees_zz(poly)) == 1


def sign_tuple(coeffs: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(1 if c >= 0 else -1 for c in coeffs)


def mod_factor_degrees(poly: sp.Poly, p: int) -> tuple[int, ...] | None:
    """Factor degrees of primitive f mod p.

    None means the leading coefficient drops degree mod p, so this prime is not
    a clean residue certificate for the degree of f.
    """
    pp = primitive_part(poly)
    if int(pp.LC()) % p == 0:
        return None
    _unit, factors = sp.factor_list(pp.as_expr(), x, modulus=p)
    out: list[int] = []
    for fac, exp in factors:
        out.extend([sp.Poly(fac, x, modulus=p).degree()] * exp)
    return tuple(sorted(out))


def split_survivors_from_degrees(degrees: tuple[int, ...] | None, degree: int) -> tuple[int, ...] | None:
    """Nontrivial degree splits surviving after a local factorization.

    The value () means no split survives, i.e. the polynomial is irreducible
    over this residue field.  A reducible factorization records possible
    smaller-side split sizes, since r x (d-r) and (d-r) x r are the same grid
    up to swapping the hidden factors.
    """
    if degrees is None:
        return None
    if len(degrees) <= 1:
        return ()
    survivors: set[int] = set()
    n = len(degrees)
    for mask in range(1, (1 << n) - 1):
        s = sum(degrees[i] for i in range(n) if mask & (1 << i))
        if 0 < s < degree:
            survivors.add(min(s, degree - s))
    return tuple(sorted(survivors))


def modp_split_survivors(poly: sp.Poly, p: int) -> tuple[int, ...] | None:
    pp = primitive_part(poly)
    return split_survivors_from_degrees(mod_factor_degrees(pp, p), pp.degree())


def least_modp_convolution_blocker(poly: sp.Poly, primes: tuple[int, ...] = SMALL_PRIMES) -> int | None:
    pp = primitive_part(poly)
    for p in primes:
        if modp_split_survivors(pp, p) == ():
            return p
    return None


def brute_force_split_survivors_modp(poly: sp.Poly, p: int) -> tuple[int, ...] | None:
    """Tiny exact diagonal-lift enumerator over F_p.

    This is used only for small examples.  It directly asks whether coefficient
    arrays b,c exist with convolution equal to f mod p for each nontrivial
    degree split.
    """
    pp = primitive_part(poly)
    d = pp.degree()
    if int(pp.LC()) % p == 0:
        return None
    coeffs = [int(pp.nth(i)) % p for i in range(d + 1)]
    survivors: set[int] = set()
    for r in range(1, d):
        s = d - r
        found = False
        for b in product(range(p), repeat=r + 1):
            if b[-1] == 0:
                continue
            for c in product(range(p), repeat=s + 1):
                if c[-1] == 0:
                    continue
                conv = [0] * (d + 1)
                for i, bi in enumerate(b):
                    for j, cj in enumerate(c):
                        conv[i + j] = (conv[i + j] + bi * cj) % p
                if conv == coeffs:
                    survivors.add(min(r, s))
                    found = True
                    break
            if found:
                break
    return tuple(sorted(survivors))


def v_p(n: int, p: int) -> int | None:
    n = abs(int(n))
    if n == 0:
        return None
    count = 0
    while n % p == 0:
        count += 1
        n //= p
    return count


def newton_points(poly: sp.Poly, p: int) -> list[tuple[int, int]]:
    pp = primitive_part(poly)
    points: list[tuple[int, int]] = []
    for i in range(pp.degree() + 1):
        coeff = int(pp.nth(i))
        if coeff:
            points.append((i, v_p(coeff, p) or 0))
    return points


def cross(o: tuple[int, int], a: tuple[int, int], b: tuple[int, int]) -> int:
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def lower_hull(points: list[tuple[int, int]]) -> list[tuple[int, int]]:
    hull: list[tuple[int, int]] = []
    for pt in points:
        while len(hull) >= 2 and cross(hull[-2], hull[-1], pt) <= 0:
            hull.pop()
        hull.append(pt)
    return hull


def newton_one_edge_certificate(poly: sp.Poly, p: int) -> dict[str, object]:
    """A small Dumas/Eisenstein-style one-edge certificate scout."""
    pts = newton_points(poly, p)
    hull = lower_hull(pts)
    cert = False
    edge_gcd = None
    if len(hull) == 2:
        (x0, y0), (x1, y1) = hull
        edge_gcd = gcd(x1 - x0, abs(y1 - y0))
        cert = x0 == 0 and x1 == primitive_part(poly).degree() and edge_gcd == 1
    return {"points": pts, "hull": hull, "edge_gcd": edge_gcd, "certificate": cert}


@dataclass(frozen=True)
class SweepRow:
    coeffs: tuple[int, ...]
    signs: tuple[int, ...]
    poly: sp.Poly
    irreducible: bool
    factor_degrees: tuple[int, ...]
    blocker: int | None


def build_degree4_sweep() -> list[SweepRow]:
    rows: list[SweepRow] = []
    for mags in product((1, 2, 3), repeat=5):
        for signs_low in product((-1, 1), repeat=4):
            coeffs = tuple(signs_low[i] * mags[i] for i in range(4)) + (mags[4],)
            poly = poly_from_coeffs(coeffs)
            rows.append(
                SweepRow(
                    coeffs=coeffs,
                    signs=sign_tuple(coeffs),
                    poly=poly,
                    irreducible=is_irreducible_zz(poly),
                    factor_degrees=factor_degrees_zz(poly),
                    blocker=least_modp_convolution_blocker(poly),
                )
            )
    return rows


def mixed_bucket_summary(rows: list[SweepRow], key_func) -> dict[str, int]:
    buckets: dict[object, list[SweepRow]] = defaultdict(list)
    for row in rows:
        buckets[key_func(row)].append(row)
    mixed = sum(1 for vals in buckets.values() if len({v.irreducible for v in vals}) > 1)
    return {
        "groups": len(buckets),
        "mixed": mixed,
        "max_bucket": max(len(vals) for vals in buckets.values()),
    }


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
    return sum(dp[-1])


def directed_3cycles(adj: tuple[tuple[int, ...], ...]) -> int:
    n = len(adj)
    total = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                outs = [sum(adj[v][u] for u in (i, j, k) if u != v) for v in (i, j, k)]
                if outs == [1, 1, 1]:
                    total += 1
    return total


def prime_certificate_tournament(rows: list[SweepRow]) -> dict[str, object]:
    counts = Counter(row.blocker for row in rows if row.irreducible)
    primes = SMALL_PRIMES
    n = len(primes)
    adj = [[0] * n for _ in range(n)]
    for i, p in enumerate(primes):
        for j, q in enumerate(primes):
            if i == j:
                continue
            # Pairwise observable: which prime is the least decisive local
            # certificate for more irreducible rows.  Tie path: smaller prime.
            if (counts[p], -p) >= (counts[q], -q):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
    tadj = tuple(tuple(row) for row in adj)
    ranking = sorted(primes, key=lambda p: (counts[p], -p), reverse=True)
    return {
        "least_blocker_counts": dict((p, counts[p]) for p in primes),
        "unresolved_irreducibles": counts[None],
        "score_hist": dict(sorted(Counter(sum(row) for row in tadj).items())),
        "directed_3cycles": directed_3cycles(tadj),
        "hamiltonian_paths": hamiltonian_paths(tadj),
        "ranking": ranking,
    }


def fmt_survivors(value: tuple[int, ...] | None) -> str:
    if value is None:
        return "degree-drop"
    if value == ():
        return "none"
    return str(value)


def main() -> None:
    print("=" * 78)
    print("Convolution-lift carriers for prime <-> irreducible polynomials")
    print("=" * 78)
    print("HYP-2451 / T795 / OPEN-Q-073 candidate")
    print()

    print("[1] Reducibility as a hidden diagonal-sum lift")
    print("  f = g*h means every coefficient a_k is a diagonal total:")
    print("    a_k = sum_{i+j=k} b_i*c_j")
    print("  Degree splits are proof obligations; the split-check order is only")
    print("  a fixed Hamiltonian path used for ties, not the invariant itself.")
    print()

    print("[2] Exact small residue lifts: symbolic factorization vs brute force")
    examples = [
        (sp.Poly(x**4 - x**3 - x**2 - x - 1, x, domain=sp.ZZ), 2, "irreducible quartic, mod-2 blocker"),
        (sp.Poly(2 * x**4 - 2 * x**3 - x**2 - x - 1, x, domain=sp.ZZ), 3, "reducible quartic, split survives"),
    ]
    for poly, p, note in examples:
        degrees = mod_factor_degrees(poly, p)
        symbolic = modp_split_survivors(poly, p)
        brute = brute_force_split_survivors_modp(poly, p)
        print(f"  {poly.as_expr()}  mod {p}  ({note})")
        print(f"    factor_degrees={degrees}, symbolic_survivors={fmt_survivors(symbolic)}, brute_survivors={fmt_survivors(brute)}")
    print()

    print("[3] Residue blocker vs Newton valuation blocker")
    newton_cases = [
        (sp.Poly(x**5 + 6 * x + 3, x, domain=sp.ZZ), 3, "Eisenstein-like: residue split survives, valuation blocks"),
        (sp.Poly(x**4 + 10 * x**2 + 5, x, domain=sp.ZZ), 5, "one-edge Newton hull blocks at p=5"),
        (sp.Poly(x**4 + 10 * x**2 + 1, x, domain=sp.ZZ), 2, "irreducible but neither this residue nor this hull blocks"),
    ]
    for poly, p, note in newton_cases:
        residue = modp_split_survivors(poly, p)
        cert = newton_one_edge_certificate(poly, p)
        print(f"  {poly.as_expr()}  p={p}  ({note})")
        print(f"    mod-p split survivors={fmt_survivors(residue)}")
        print(
            f"    newton_points={cert['points']}, lower_hull={cert['hull']}, "
            f"edge_gcd={cert['edge_gcd']}, certificate={cert['certificate']}"
        )
    print()

    print("[4] Degree-4 coefficient-row scout, |a_i| in {1,2,3}")
    rows = build_degree4_sweep()
    irr = [r for r in rows if r.irreducible]
    red = [r for r in rows if not r.irreducible]
    certified = [r for r in irr if r.blocker is not None]
    false_certified = [r for r in red if r.blocker is not None]
    print(f"  rows={len(rows)}, irreducible={len(irr)}, reducible={len(red)}")
    print(
        f"  least mod-p convolution blocker in {SMALL_PRIMES}: "
        f"{len(certified)}/{len(irr)} irreducibles certified "
        f"({100 * len(certified) / len(irr):.2f}%), false positives={len(false_certified)}"
    )
    print(f"  blocker histogram={dict(sorted(Counter(r.blocker for r in irr).items(), key=lambda kv: (999 if kv[0] is None else kv[0])))}")
    for name, key_func in [
        ("signs", lambda r: r.signs),
        ("signs+least_blocker", lambda r: (r.signs, r.blocker)),
    ]:
        summary = mixed_bucket_summary(rows, key_func)
        print(
            f"  key={name:20s} groups={summary['groups']:3d} "
            f"mixed={summary['mixed']:3d} max_bucket={summary['max_bucket']:3d}"
        )
    leftovers = [r for r in irr if r.blocker is None][:4]
    print("  sample irreducibles with no blocker up to 31:")
    for row in leftovers:
        sig = tuple((p, modp_split_survivors(row.poly, p)) for p in (2, 3, 5, 7))
        print(f"    coeffs={row.coeffs}; f={row.poly.as_expr()}; small_signature={sig}")
    print()

    print("[5] Prime-certificate tournament")
    tour = prime_certificate_tournament(rows)
    print("  pairwise observable: more least-blocker wins on the degree-4 scout")
    print("  switch/gauge: primitive part by Gauss lemma; degree-drop primes ignored")
    print("  tie Hamiltonian path: 2 -> 3 -> 5 -> ... -> 31")
    print(f"  least_blocker_counts={tour['least_blocker_counts']}")
    print(f"  unresolved_irreducibles={tour['unresolved_irreducibles']}")
    print(f"  score_hist={tour['score_hist']}")
    print(f"  directed_3cycles={tour['directed_3cycles']}")
    print(f"  Hamiltonian_paths={tour['hamiltonian_paths']}")
    print(f"  ranking={tour['ranking']}")
    print()

    print("[6] Extracted proof/program directions")
    directions = [
        "Use mod-p diagonal-lift infeasibility as a cheap exact irreducibility prefilter before expensive Z-factorization.",
        "Use Newton/valuation hulls for the cases where residue reduction collapses to a reducible shadow, especially Eisenstein/Dumas rows.",
        "Treat the leftover no-small-blocker rows as the proper target for Singh-depth and recombination certificates.",
        "For HYP-2449, replace sign-only chambers by split-survivor signatures: which hidden degree rectangles remain possible?",
        "For LRC14, copy the same split-survivor language: a dangerous residue fiber is one whose local lift survivors never empty.",
    ]
    for i, direction in enumerate(directions, 1):
        print(f"  {i}. {direction}")
    print()
    print("Conclusion:")
    print("  The coefficient row is a diagonal projection of a hidden factor grid.")
    print("  Residue blockers empty the possible split set; Newton blockers see")
    print("  valuation rigidity when the residue shadow is misleading.  The next")
    print("  carrier should track split survivors, not just coefficient signs.")


if __name__ == "__main__":
    main()
