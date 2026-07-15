#!/usr/bin/env python3
"""Exact closure of the THM-841 ladder and toothpick-recursion refutation.

For consecutive Farey denominators i,j of order k, i+j>k.  Hence

    min(floor(k/i), floor(k/j)) = 1.

Only one end of a Farey gap can therefore carry a nontrivial multiple
staircase.  Moreover, a collision between the m-th left copy and the n-th
right copy occurs at

    lambda = mn / (n*j + m*i).

If one of m,n is 1 and the other is at least 2, this lies strictly above
1/(k+1).  Thus, in the THM-841 window, only the primitive collision
lambda=1/(i+j) occurs.  The apparent divisor-refined breakpoint tree
collapses to the original Farey breakpoint path.

The script verifies the resulting formulas for every order statistic m_r
and every factorial moment S_rho against an independent exact interval
sweep.  It also tests the literal Farey-end version of A139250 and records
its first recurrence and dyadic failures.

No theorem or hypothesis identifier is claimed by this computation.
"""

from __future__ import annotations

import hashlib
from fractions import Fraction as F
from functools import lru_cache
from math import comb, gcd
from pathlib import Path


THM841 = Path("01-canon/theorems/THM-841-violation-ladder-order-statistic-profiles.md")
THM826 = Path("01-canon/theorems/THM-826-farey-profile-theorem-interval-core-measure.md")
REFEREE = Path("04-computation/thm841_ladder_ouH_diffusion_referee_kps_S128c12.py")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def phi(n: int) -> int:
    if n == 1:
        return 1
    return sum(gcd(a, n) == 1 for a in range(1, n + 1))


def farey(k: int) -> list[F]:
    """Circular Farey representatives in [0,1), with 0/1 appearing once."""
    return sorted({F(a, d) for d in range(1, k + 1) for a in range(d)})


def farey_gaps(k: int) -> list[tuple[int, int, F]]:
    fractions = farey(k)
    answer = []
    for left, right in zip(fractions, fractions[1:] + [F(1)]):
        i = left.denominator
        j = 1 if right == 1 else right.denominator
        answer.append((i, j, right - left))
    return answer


def overlap_sum(k: int, lam: F) -> F:
    """Total overlap of the two primitive endpoint arcs across Farey gaps."""
    return sum(
        (max(F(0), F(lam * (i + j) - 1, i * j)) for i, j, _ in farey_gaps(k)),
        F(0),
    )


def formula_m(k: int, lam: F, r: int) -> F:
    """Closed form for m_r(lambda)=measure{W<r}."""
    if r == 1:
        return sum(
            (max(F(0), F(1 - lam * (i + j), i * j)) for i, j, _ in farey_gaps(k)),
            F(0),
        )
    if r == 2:
        nested = lam * sum((F(phi(d), d) for d in range(1, k // 2 + 1)), F(0))
        return 1 - nested - overlap_sum(k, lam)
    tail = F(2) * lam / r * sum(
        (F(phi(d), d) for d in range(1, k // r + 1)), F(0)
    )
    return 1 - tail


def formula_s(k: int, lam: F, rho: int) -> F:
    """Closed form for S_rho(lambda)=integral binom(W,rho)."""
    if rho == 1:
        return 2 * k * lam
    answer = overlap_sum(k, lam) if rho == 2 else F(0)
    for d in range(1, k // rho + 1):
        p = k // d
        chain = sum(
            (F(comb(m - 1, rho - 1), m) for m in range(rho, p + 1)),
            F(0),
        )
        answer += 2 * lam * F(phi(d), d) * chain
    return answer


def occupancy_histogram(k: int, lam: F) -> list[F]:
    """Independent exact circle sweep of W(t)=#{j<=k: ||jt||<lambda}."""
    events: dict[F, int] = {}

    def event(x: F, delta: int) -> None:
        events[x] = events.get(x, 0) + delta

    for speed in range(1, k + 1):
        radius = lam / speed
        for numerator in range(speed):
            center = F(numerator, speed)
            lo = center - radius
            hi = center + radius
            if lo < 0:
                event(F(0), 1)
                event(hi, -1)
                event(lo + 1, 1)
                event(F(1), -1)
            elif hi > 1:
                event(lo, 1)
                event(F(1), -1)
                event(F(0), 1)
                event(hi - 1, -1)
            else:
                event(lo, 1)
                event(hi, -1)

    histogram = [F(0) for _ in range(k + 1)]
    depth = 0
    previous = F(0)
    for position in sorted(events):
        if position > previous:
            histogram[depth] += position - previous
            previous = position
        depth += events[position]
    if previous < 1:
        histogram[depth] += 1 - previous
    if sum(histogram, F(0)) != 1 or depth != 0:
        raise RuntimeError("circle sweep failed to close")
    return histogram


def audit_structure(max_k: int = 64) -> tuple[int, int]:
    gaps_checked = 0
    nonprimitive_collisions_checked = 0
    for k in range(2, max_k + 1):
        denominator_sums = set()
        for i, j, gap_length in farey_gaps(k):
            gaps_checked += 1
            denominator_sums.add(i + j)
            if gap_length != F(1, i * j) or i + j <= k:
                raise RuntimeError("Farey-neighbour identity failed")
            p, q = k // i, k // j
            if min(p, q) != 1:
                raise RuntimeError("both Farey endpoints have nontrivial chains")
            for m in range(1, p + 1):
                for n in range(1, q + 1):
                    collision = F(m * n, n * j + m * i)
                    if (m, n) != (1, 1):
                        nonprimitive_collisions_checked += 1
                        if collision <= F(1, k + 1):
                            raise RuntimeError("nonprimitive collision entered the THM-841 window")

        # THM-826's claimed complete range has one systematic omission.
        if k % 2 == 0 and 2 * k - 2 in denominator_sums:
            raise RuntimeError("even-k forbidden denominator sum appeared")
        if k % 2 == 1 and 2 * k - 2 not in denominator_sums:
            raise RuntimeError("odd-k denominator sum unexpectedly absent")
    return gaps_checked, nonprimitive_collisions_checked


def audit_formulas(max_k: int = 20) -> tuple[int, int]:
    rows = 0
    equalities = 0
    for k in range(2, max_k + 1):
        lambdas = {F(a, 12 * (k + 1)) for a in range(1, 12)}
        lambdas.update(F(1, s) for s in range(k + 2, 2 * k))
        for lam in sorted(lambdas):
            if not 0 < lam < F(1, k + 1):
                continue
            histogram = occupancy_histogram(k, lam)
            for r in range(1, k + 1):
                direct_m = sum(histogram[:r], F(0))
                if direct_m != formula_m(k, lam, r):
                    raise RuntimeError(f"m_r mismatch at k={k}, lambda={lam}, r={r}")
                equalities += 1
            for rho in range(1, k + 1):
                direct_s = sum(
                    (comb(w, rho) * measure for w, measure in enumerate(histogram)),
                    F(0),
                )
                if direct_s != formula_s(k, lam, rho):
                    raise RuntimeError(f"S_rho mismatch at k={k}, lambda={lam}, rho={rho}")
                equalities += 1
            rows += 1
    return rows, equalities


@lru_cache(maxsize=None)
def toothpick_added(n: int) -> int:
    """A139251 recurrence from Applegate--Pol--Sloane."""
    if n == 0:
        return 0
    if n & (n - 1) == 0:
        return n
    power = 1 << (n.bit_length() - 1)
    remainder = n - power
    return 2 * toothpick_added(remainder) + toothpick_added(remainder + 1)


def toothpick_total(n: int) -> int:
    return sum(toothpick_added(stage) for stage in range(1, n + 1))


def farey_end_added(n: int) -> int:
    return 1 if n == 1 else 2 * phi(n)


def farey_end_total(n: int) -> int:
    return 2 * sum(phi(d) for d in range(1, n + 1)) - 1


def audit_dyadic_refutation() -> None:
    for p in (1, 2, 4, 8, 16, 32):
        doubled = {F(1, m) for m in range(1, 2 * p + 1)}
        inherited = {F(1, m) for m in range(1, p + 1)} | {
            F(1, 2 * m) for m in range(1, p + 1)
        }
        defect = doubled - inherited
        expected = {F(1, m) for m in range(p + 1, 2 * p + 1) if m % 2}
        if defect != expected or len(defect) != p // 2:
            raise RuntimeError("reciprocal-switch dyadic defect formula failed")

    if farey_end_added(5) != 8:
        raise RuntimeError("Farey-end count changed")
    if 2 * farey_end_added(1) + farey_end_added(2) != 4:
        raise RuntimeError("toothpick recurrence comparison changed")
    dyadic = (1, 2, 4, 8, 16)
    agreement = [farey_end_total(n) == toothpick_total(n) for n in dyadic]
    if agreement != [True, True, True, True, False]:
        raise RuntimeError("dyadic coincidence/refutation profile changed")


def main() -> None:
    gaps, collisions = audit_structure()
    rows, equalities = audit_formulas()
    audit_dyadic_refutation()

    print("=" * 88)
    print("THM-841: NO DYADIC BREAKPOINT TREE; EXACT ALL-ORDER CLOSURE")
    print("=" * 88)
    print("local lemma: i+j>k => min(floor(k/i),floor(k/j))=1")
    print("collision: lambda_(m,n)=mn/(n*j+m*i)")
    print("nonprimitive collision: lambda_(m,n)>1/(k+1)")
    print("actual lambda-breakpoints: primitive lambda=1/(i+j) only")
    print("parameter graph per Farey gap: a path with at most one kink")
    print()
    print("m_1 = sum_g [1-lambda(i+j)]_+/(ij)")
    print("m_2 = 1-lambda*sum_(d<=floor(k/2)) phi(d)/d-O_k(lambda)")
    print("m_r = 1-(2lambda/r)*sum_(d<=floor(k/r)) phi(d)/d, r>=3")
    print("O_k(lambda)=sum_g [lambda(i+j)-1]_+/(ij)")
    print("S_1 = 2k lambda")
    print("S_rho = 2lambda*sum_d phi(d)/d*sum_(m=rho..floor(k/d)) C(m-1,rho-1)/m")
    print("        + 1_(rho=2) O_k(lambda), rho>=2")
    print()
    print(f"structure audit: Farey gaps={gaps}, nonprimitive collisions={collisions}, k<=64")
    print(f"formula audit: (k,lambda) rows={rows}, exact m_r/S_rho equalities={equalities}, k<=20")
    print("formula verdict: ALL EXACT against independent circle interval sweep")
    print()
    print("DYADIC / TOOTHPICK REFUTATION")
    print("D_(2p)=D_p union (1/2)D_p union {1/m: p<m<=2p, m odd}")
    print("first local defect: D_4 \\ (D_2 union (1/2)D_2) = {1/3}")
    print("actual witness gap: (0/1,1/4)")
    print("literal Farey-end additions b(1)=1, b(n)=2phi(n) for n>=2")
    print("first recurrence failure: b(5)=8, but 2b(1)+b(2)=4=A139251(5)")
    print("dyadic totals (n, Farey-end total, A139250 total):")
    for n in (1, 2, 4, 8, 16):
        print(f"  {n:2d}: {farey_end_total(n):3d}  {toothpick_total(n):3d}")
    print("the 1,2,4,8 agreement is a finite mirage; 159 != 171 at n=16")
    print()
    print("THM-826 CROSS-CORRECTION")
    print("for even k, denominator sum 2k-2 is absent: candidates (k-2,k),(k-1,k-1),(k,k-2)")
    print("are non-coprime/equal; it is not an actual breakpoint")
    print()
    print("ASSUMPTION CHALLENGE / TOURNAMENT ANALYSIS")
    print("vertices tested: endpoint multiple-events, spatial cells, lambda-collisions, Farey gaps")
    print("predicate preserved by the endpoint carrier: the full violation count W(t)")
    print("information destroyed by dyadic quotient: odd multipliers p<m<=2p")
    print("no tournament is formed: the exact switch relation is laminar/path-like, with no cycles or SCC choice")
    print("challenged assumption: two-sided chains can branch; Farey adjacency forces one side to be primitive-only")
    print()
    print("PROVENANCE")
    print(f"THM841_sha256={sha256(THM841)}")
    print(f"THM826_sha256={sha256(THM826)}")
    print(f"legacy_referee_sha256={sha256(REFEREE)}")
    print(f"script_sha256={sha256(Path(__file__))}")


if __name__ == "__main__":
    main()
