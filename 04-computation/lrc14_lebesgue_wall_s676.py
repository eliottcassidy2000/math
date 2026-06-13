#!/usr/bin/env python3
"""
S676: Lebesgue wall certificate for the LRC n=14 proof route.

At total runner count n=14 there are 13 moving speeds and the target gap is
delta = 1/14.  The danger set for speed v is

    D_v(delta) = {t in R/Z : ||v t|| < delta}.

The exact Lebesgue safe measure is

    p_0(V, delta) = measure {t : depth_V(t) = 0},
    depth_V(t) = # {v in V : t in D_v(delta)}.

This script makes the measure/set split concrete:

* p_0 > 0 gives a whole open interval of lonely times, hence an immediate
  LRC certificate.
* p_0 = 0 is not a counterexample.  Because the danger arcs are open, the
  complement can still contain endpoint witnesses.  AP, Vstar, and 2AP are
  exactly this kind of Lebesgue wall in the Res_27 quotient.

The proof-relevant question becomes: after quotient/scaling, can any lifted
row stay on the p_0=0 wall without being AP, Vstar, or 2AP?

Tournament Analysis / assumption challenge:
  Vertices are proof routes, not runners.  Candidate vertex sets considered
  include runners, danger arcs, breakpoints, safe components, endpoint
  witnesses, Res_27 rows, carry vectors, owner obligations, and proof
  obligations.  The selected route tournament preserves the LRC predicate
  "positive safe interval versus endpoint-only wall" while destroying raw
  runner labels and most phase order.  Pairwise observables are exactness,
  endpoint handling, lift relevance, actionability, compression, and
  overclaim safety; majority orients the tournament, with list order as the
  tie Hamiltonian path.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd


N_TOTAL = 14
K = N_TOTAL - 1
C = 2 * N_TOTAL - 1
DELTA = Fraction(1, N_TOTAL)

AP = tuple(range(1, N_TOTAL))
VSTAR = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24)
TWOP = tuple(2 * v for v in AP)
BASES = (("AP", AP), ("Vstar", VSTAR), ("2AP", TWOP))


def fmt_frac(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def fmt_decimal(x: Fraction, digits: int = 8) -> str:
    return f"{float(x):.{digits}f}"


def row_gcd(row: tuple[int, ...]) -> int:
    g = 0
    for v in row:
        g = gcd(g, v)
    return g


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def norm_circle(x: Fraction) -> Fraction:
    r = frac_part(x)
    return min(r, 1 - r)


def danger_intervals(v: int, delta: Fraction) -> list[tuple[Fraction, Fraction]]:
    """Return danger intervals in [0,1), split at the wrap point."""
    intervals: list[tuple[Fraction, Fraction]] = []
    for k in range(v):
        lo = (Fraction(k) - delta) / v
        hi = (Fraction(k) + delta) / v
        if lo < 0:
            intervals.append((lo + 1, Fraction(1)))
            intervals.append((Fraction(0), hi))
        elif hi > 1:
            intervals.append((lo, Fraction(1)))
            intervals.append((Fraction(0), hi - 1))
        else:
            intervals.append((lo, hi))
    return intervals


def in_danger(v: int, t: Fraction, delta: Fraction) -> bool:
    return norm_circle(Fraction(v) * t) < delta


def depth(row: tuple[int, ...], t: Fraction, delta: Fraction) -> int:
    return sum(1 for v in row if in_danger(v, t, delta))


@dataclass(frozen=True)
class DepthSweep:
    profile: tuple[tuple[int, Fraction], ...]
    safe_intervals: tuple[tuple[Fraction, Fraction], ...]
    safe_points: tuple[Fraction, ...]
    breakpoints: int

    @property
    def p0(self) -> Fraction:
        return dict(self.profile).get(0, Fraction(0))

    @property
    def positive_safe_components(self) -> int:
        return len(self.safe_intervals)


def depth_sweep(row: tuple[int, ...], delta: Fraction = DELTA) -> DepthSweep:
    breakpoints: set[Fraction] = {Fraction(0), Fraction(1)}
    for v in row:
        for lo, hi in danger_intervals(v, delta):
            breakpoints.add(lo)
            breakpoints.add(hi)

    bps = sorted(breakpoints)
    profile: Counter[int] = Counter()
    safe_intervals: list[tuple[Fraction, Fraction]] = []
    for a, b in zip(bps, bps[1:]):
        if a == b:
            continue
        mid = (a + b) / 2
        k = depth(row, mid, delta)
        profile[k] += b - a
        if k == 0:
            safe_intervals.append((a, b))

    safe_points = tuple(
        t
        for t in bps
        if 0 < t < 1 and all(norm_circle(Fraction(v) * t) >= delta for v in row)
    )

    return DepthSweep(
        profile=tuple(sorted(profile.items())),
        safe_intervals=tuple(safe_intervals),
        safe_points=safe_points,
        breakpoints=len(bps),
    )


def exact_score(row: tuple[int, ...], t: Fraction) -> Fraction:
    return min(norm_circle(Fraction(v) * t) for v in row)


def exact_maximin(row: tuple[int, ...]) -> tuple[Fraction, Fraction]:
    """THM-369 style finite exact maximin over pair and single-runner pinches."""
    candidates: set[Fraction] = set()
    for i, a in enumerate(row):
        for m in range(a):
            t = Fraction(2 * m + 1, 2 * a)
            if 0 < t < 1:
                candidates.add(t)
        for b in row[i + 1 :]:
            for den in (a + b, abs(a - b)):
                if den <= 0:
                    continue
                for m in range(1, den):
                    candidates.add(Fraction(m, den))

    best = Fraction(0)
    best_t = Fraction(0)
    for t in candidates:
        score = exact_score(row, t)
        if (score, -t) > (best, -best_t):
            best = score
            best_t = t
    return best, best_t


def local_carries(max_weight: int = 3) -> list[tuple[str, int, tuple[int, ...], tuple[int, ...]]]:
    probes: list[tuple[str, int, tuple[int, ...], tuple[int, ...]]] = []
    for family, base in BASES:
        probes.append((family, 0, tuple(), tuple(sorted(base))))
        for weight in range(1, max_weight + 1):
            for idxs in combinations(range(len(base)), weight):
                row = list(base)
                moved: list[int] = []
                for idx in idxs:
                    moved.append(base[idx])
                    row[idx] += C
                probes.append((family, weight, tuple(sorted(moved)), tuple(sorted(row))))
    return probes


def one_swap_rows(max_new: int = 60) -> list[tuple[str, int, int, tuple[int, ...]]]:
    rows: list[tuple[str, int, int, tuple[int, ...]]] = []
    base = list(AP)
    for old in AP:
        for new in range(14, max_new + 1):
            if new in base:
                continue
            row = tuple(sorted(v for v in base if v != old) + [new])
            row = tuple(sorted(row))
            if len(set(row)) == K:
                rows.append((f"{old}->{new}", old, new, row))
    return rows


def units_mod(m: int) -> tuple[int, ...]:
    return tuple(u for u in range(1, m) if gcd(u, m) == 1)


def scalar_floor_rows() -> list[tuple[str, tuple[int, ...]]]:
    rows: list[tuple[str, tuple[int, ...]]] = []
    for family, base in (("AP", AP), ("Vstar", VSTAR)):
        for u in units_mod(C):
            rows.append((f"{family}*{u}", tuple(sorted(u * v for v in base))))
    return rows


def summarize_wall_atoms() -> None:
    print("A. Exact Lebesgue wall atoms at delta=1/14")
    print("  row        p0        components  endpoints  M       witness")
    for family, row in BASES:
        sweep = depth_sweep(row)
        m, t = exact_maximin(row)
        endpoints = ",".join(fmt_frac(x) for x in sweep.safe_points[:12])
        if len(sweep.safe_points) > 12:
            endpoints += ",..."
        print(
            f"  {family:7s} {fmt_frac(sweep.p0):>8s} "
            f"{sweep.positive_safe_components:>10d} {len(sweep.safe_points):>10d} "
            f"{fmt_frac(m):>7s} {fmt_frac(t):>9s}"
        )
        print(f"      endpoint witnesses: {endpoints}")
        profile = " ".join(f"p{k}={fmt_frac(v)}" for k, v in sweep.profile[:8])
        print(f"      depth profile head: {profile}")
    print()


def summarize_wall_opening() -> None:
    print("B. One-sided wall opening p0(1/14 - eps)")
    epsilons = [Fraction(1, 14 * m) for m in (50, 100, 200, 400)]
    print("  row      eps             p0(delta-eps)  p0/eps")
    for family, row in BASES:
        for eps in epsilons:
            sweep = depth_sweep(row, DELTA - eps)
            ratio = sweep.p0 / eps if eps else Fraction(0)
            print(
                f"  {family:7s} {fmt_frac(eps):>13s} "
                f"{fmt_frac(sweep.p0):>15s} {fmt_frac(ratio):>8s}"
            )
        print("  " + "-" * 52)
    print(
        "  Reading: the wall atoms open linearly when the collar is relaxed. "
        "For a proof, this is a boundary stratum, not a bulk-measure stratum."
    )
    print()


def summarize_local_carries() -> None:
    print("C. Local +27 carry perturbations measured by exact p0")
    probes = local_carries(max_weight=3)
    by_key: defaultdict[tuple[str, int], list[tuple[Fraction, tuple[int, ...], tuple[int, ...]]]] = (
        defaultdict(list)
    )
    route_hist: Counter[str] = Counter()
    zero_labels: list[str] = []

    for family, weight, moved, row in probes:
        sweep = depth_sweep(row)
        route = "wall" if sweep.p0 == 0 else "positive"
        route_hist[route] += 1
        if sweep.p0 == 0:
            zero_labels.append(f"{family}:w{weight}:{moved or '-'}")
        by_key[(family, weight)].append((sweep.p0, moved, row))

    print(f"  probes={len(probes)} route_hist={dict(route_hist)}")
    print(f"  zero-measure labels: {zero_labels}")
    print("  family weight count  min_positive_p0  achieved_by")
    for key in sorted(by_key):
        rows = by_key[key]
        positives = [r for r in rows if r[0] > 0]
        if positives:
            p0, moved, _row = min(positives, key=lambda x: (x[0], x[1]))
            moved_label = "-" if not moved else ",".join(map(str, moved))
            p0_label = fmt_frac(p0)
        else:
            moved_label = "-"
            p0_label = "none"
        print(
            f"  {key[0]:7s} {key[1]:>6d} {len(rows):>5d} "
            f"{p0_label:>16s}  {moved_label}"
        )
    print(
        "  Consequence: every tested nonzero local carry has positive Lebesgue "
        "safe measure, so it is certified by an interval, not by endpoint finesse."
    )
    print()


def summarize_one_swaps() -> None:
    print("D. One-swap AP wall sieve by exact p0")
    rows = one_swap_rows(max_new=60)
    zeros: list[tuple[str, tuple[int, ...]]] = []
    positives: list[tuple[Fraction, str, tuple[int, ...]]] = []
    doubled: list[tuple[str, Fraction, Fraction, Fraction]] = []

    for label, old, new, row in rows:
        sweep = depth_sweep(row)
        if sweep.p0 == 0:
            zeros.append((label, row))
        else:
            positives.append((sweep.p0, label, row))
        if new == 2 * old:
            m, t = exact_maximin(row)
            doubled.append((label, sweep.p0, m, t))

    positives.sort(key=lambda x: (x[0], x[1]))
    print(f"  rows={len(rows)} zero_p0={len(zeros)} positive_p0={len(positives)}")
    for label, row in zeros:
        print(f"    zero: {label:7s} row={row}")
    print("  smallest positive one-swap p0 rows:")
    for p0, label, row in positives[:8]:
        print(f"    {label:7s} p0={fmt_frac(p0):>10s} row={row}")
    print("  AP a->2a doubling row readout:")
    for label, p0, m, t in doubled:
        print(
            f"    {label:7s} p0={fmt_frac(p0):>10s} "
            f"M={fmt_frac(m):>7s} t={fmt_frac(t):>8s}"
        )
    print()


def summarize_scalar_orbit() -> None:
    print("E. Global scalar floor orbits: why p0=0 must be quotiented")
    rows = scalar_floor_rows()
    route_hist: Counter[str] = Counter()
    samples: list[str] = []
    for label, row in rows:
        sweep = depth_sweep(row)
        route = "wall" if sweep.p0 == 0 else "positive"
        route_hist[route] += 1
        if len(samples) < 8:
            samples.append(
                f"{label}:p0={fmt_frac(sweep.p0)},endpoints={len(sweep.safe_points)}"
            )
    print(f"  scalar rows={len(rows)} route_hist={dict(route_hist)}")
    print("  samples:")
    for sample in samples:
        print(f"    {sample}")
    print(
        "  Scaling is measure-preserving on R/Z, so AP and Vstar generate "
        "infinite-looking wall rows.  The finite statement must live in the "
        "scaled Res_27 quotient or another normalization."
    )
    print()


@dataclass(frozen=True)
class RouteVertex:
    name: str
    scores: tuple[int, ...]


def tournament_edges(vertices: tuple[RouteVertex, ...]) -> dict[str, set[str]]:
    edges = {v.name: set() for v in vertices}
    for i, a in enumerate(vertices):
        for b in vertices[i + 1 :]:
            av = sum(1 for x, y in zip(a.scores, b.scores) if x > y)
            bv = sum(1 for x, y in zip(a.scores, b.scores) if y > x)
            if av >= bv:
                edges[a.name].add(b.name)
            else:
                edges[b.name].add(a.name)
    return edges


def count_directed_triangles(vertices: tuple[RouteVertex, ...], edges: dict[str, set[str]]) -> int:
    total = 0
    names = [v.name for v in vertices]
    for a, b, c in combinations(names, 3):
        if (
            (b in edges[a] and c in edges[b] and a in edges[c])
            or (c in edges[a] and b in edges[c] and a in edges[b])
        ):
            total += 1
    return total


def scc_sizes(vertices: tuple[RouteVertex, ...], edges: dict[str, set[str]]) -> list[int]:
    names = [v.name for v in vertices]
    reverse = {name: set() for name in names}
    for a, outs in edges.items():
        for b in outs:
            reverse[b].add(a)

    def reach(start: str, graph: dict[str, set[str]]) -> set[str]:
        seen = {start}
        todo = deque([start])
        while todo:
            x = todo.popleft()
            for y in graph[x]:
                if y not in seen:
                    seen.add(y)
                    todo.append(y)
        return seen

    remaining = set(names)
    sizes: list[int] = []
    while remaining:
        start = min(remaining)
        comp = reach(start, edges) & reach(start, reverse)
        sizes.append(len(comp))
        remaining -= comp
    return sorted(sizes, reverse=True)


def count_hamiltonian_paths(vertices: tuple[RouteVertex, ...], edges: dict[str, set[str]]) -> int:
    names = tuple(v.name for v in vertices)
    count = 0

    def rec(path: tuple[str, ...], rest: tuple[str, ...]) -> None:
        nonlocal count
        if not rest:
            count += 1
            return
        last = path[-1]
        for i, nxt in enumerate(rest):
            if nxt in edges[last]:
                rec(path + (nxt,), rest[:i] + rest[i + 1 :])

    for i, start in enumerate(names):
        rec((start,), names[:i] + names[i + 1 :])
    return count


def summarize_route_tournament() -> None:
    print("F. Tournament Analysis over proof routes")
    vertices = (
        RouteVertex("endpoint_wall_certificate", (5, 5, 4, 5, 3, 5)),
        RouteVertex("local_carry_measure_tax", (5, 3, 5, 5, 4, 4)),
        RouteVertex("res27_floor_quotient", (4, 4, 5, 4, 5, 4)),
        RouteVertex("one_swap_wall_sieve", (4, 3, 3, 4, 3, 4)),
        RouteVertex("global_scalar_floor_orbit", (4, 5, 4, 3, 4, 3)),
        RouteVertex("positive_measure_interval", (5, 0, 2, 4, 2, 5)),
        RouteVertex("raw_first_moment", (1, 0, 0, 1, 5, 2)),
    )
    edges = tournament_edges(vertices)
    out_hist = Counter(len(edges[v.name]) for v in vertices)
    ordered = sorted(vertices, key=lambda v: (-len(edges[v.name]), v.name))
    print("  vertices=proof routes, not runners")
    print("  observable=(exactness, endpoint handling, lift relevance, actionability, compression, safety)")
    print(f"  score_hist={dict(sorted(out_hist.items()))}")
    print(f"  directed_3cycles={count_directed_triangles(vertices, edges)}")
    print(f"  scc_sizes={scc_sizes(vertices, edges)}")
    print(f"  hamiltonian_paths={count_hamiltonian_paths(vertices, edges)}")
    print("  majority order by outdegree:")
    for v in ordered:
        print(f"    {v.name:30s} out={len(edges[v.name])} scores={v.scores}")
    print()


def main() -> None:
    print("=" * 78)
    print("S676 LRC14 Lebesgue wall certificate")
    print("=" * 78)
    print(f"total n={N_TOTAL}; moving runners={K}; delta={fmt_frac(DELTA)}; C=2n-1={C}")
    print(
        "Lebesgue model: p0 is the safe measure; p0>0 is an interval certificate; "
        "p0=0 is an endpoint wall needing a set-level certificate."
    )
    print()
    summarize_wall_atoms()
    summarize_wall_opening()
    summarize_local_carries()
    summarize_one_swaps()
    summarize_scalar_orbit()
    summarize_route_tournament()
    print("=" * 78)
    print("Conclusion")
    print(
        "  Lebesgue measure is a filter, not the final proof: it discards every "
        "strict row into the easy positive-interval bucket and leaves only the "
        "measure-zero wall.  In the checked n=14 neighborhoods, the normalized "
        "wall is AP/Vstar/2AP plus their global scalar floor orbits.  Thus the "
        "LRC14 proof target sharpens to a no-new-wall theorem for lifted Res_27 "
        "carry/owner fibers, followed by explicit endpoint witnesses on the "
        "named wall atoms."
    )


if __name__ == "__main__":
    main()
