#!/usr/bin/env python3
"""
S677: Improve the HYP-2252 no-new-wall target by splitting off apex debt.

S676 says the LRC14 Lebesgue proof target is:

    no new normalized p_0=0 wall in the Res_27 carry/owner fiber.

This script sharpens that target.  AP and Vstar are p_0=0 walls, but they are
no-multiple rows, so the old n-clock t=1/14 already witnesses them.  The
dangerous branch for a full LRC14 proof is the primitive multiple-of-14 branch:
once a speed is divisible by 14, all six unit n-clock witnesses are killed.

In the carry model v = r + 27k, 27 == -1 (mod 14), so

    14 | v  <=>  k == r (mod 14).

That congruence is the apex debt.  S677 audits exact p_0 on coherent families
that deliberately create or move apex debt:

  * single_apex: one coordinate is carried to its first multiple of 14;
  * apex_plus_one: a single apex debt plus one ordinary +27 side carry;
  * interval_block: contiguous carry blocks of height 1 or 2;
  * affine_mod14: small affine carry laws k = a*r+b (mod 14).

The result is not a proof of all carry fibers.  It improves the theorem target
by separating harmless endpoint walls from primitive-multiple hidden walls.

Tournament Analysis / assumption challenge:
  Vertices are proof filters, not runners.  Candidate LRC vertex sets
  considered include runners, n-clock endpoints, gcd-scaled endpoints, apex
  congruence sites, carry vectors, safe intervals, owner obligations, and proof
  obligations.  The selected quotient preserves "which proof filter discharges
  a row" and destroys raw phase order.  Pairwise observables are exactness,
  multiple-branch relevance, endpoint handling, side-channel content,
  compression, and proof actionability.
"""

from __future__ import annotations

import sys
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_lebesgue_wall_s676 as s676  # noqa: E402


N = s676.N_TOTAL
C = s676.C
DELTA = s676.DELTA
AP = s676.AP
VSTAR = s676.VSTAR
BASES = (("AP", AP), ("Vstar", VSTAR))
UNIT_CLOCKS = (1, 3, 5, 9, 11, 13)


def fmt_frac(x: Fraction) -> str:
    return s676.fmt_frac(x)


def row_gcd(row: tuple[int, ...]) -> int:
    g = 0
    for v in row:
        g = gcd(g, v)
    return g


def carry_to_first_multiple(v: int) -> int:
    """Smallest nonnegative k with v + 27k divisible by 14."""
    return v % N


def lift(base: tuple[int, ...], carry: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(v + C * k for v, k in zip(base, carry)))


def has_multiple(row: tuple[int, ...]) -> bool:
    return any(v % N == 0 for v in row)


def n_clock_witness_count(row: tuple[int, ...]) -> int:
    return sum(
        1
        for j in UNIT_CLOCKS
        if all(s676.norm_circle(Fraction(v * j, N)) >= DELTA for v in row)
    )


def gcd_scaled_clock_witness_count(row: tuple[int, ...]) -> int:
    g = row_gcd(row)
    den = N * g
    return sum(
        1
        for j in range(1, den)
        if gcd(j, den) == 1
        and all(s676.norm_circle(Fraction(v * j, den)) >= DELTA for v in row)
    )


@dataclass(frozen=True)
class Probe:
    family: str
    base_name: str
    label: str
    carry: tuple[int, ...]
    row: tuple[int, ...]
    p0: Fraction
    components: int
    safe_points: int
    gcd_value: int
    multiple: bool
    n_clock: int
    gcd_clock: int

    @property
    def route(self) -> str:
        if self.p0 == 0:
            return "wall"
        return "positive"


def make_probe(family: str, base_name: str, label: str, base: tuple[int, ...], carry: tuple[int, ...]) -> Probe:
    row = lift(base, carry)
    sweep = s676.depth_sweep(row)
    return Probe(
        family=family,
        base_name=base_name,
        label=label,
        carry=carry,
        row=row,
        p0=sweep.p0,
        components=sweep.positive_safe_components,
        safe_points=len(sweep.safe_points),
        gcd_value=row_gcd(row),
        multiple=has_multiple(row),
        n_clock=n_clock_witness_count(row),
        gcd_clock=gcd_scaled_clock_witness_count(row),
    )


def single_apex_vectors(base: tuple[int, ...]):
    for i, v in enumerate(base):
        k = carry_to_first_multiple(v)
        if k == 0:
            continue
        carry = [0] * len(base)
        carry[i] = k
        yield f"i{i}:v{v}:+{k}", tuple(carry)


def apex_plus_one_vectors(base: tuple[int, ...]):
    for i, v in enumerate(base):
        k = carry_to_first_multiple(v)
        if k == 0:
            continue
        for j in range(len(base)):
            if j == i:
                continue
            carry = [0] * len(base)
            carry[i] = k
            carry[j] += 1
            yield f"apex{i}:v{v}:+{k};side{j}:+1", tuple(carry)


def interval_block_vectors(base: tuple[int, ...]):
    n = len(base)
    for height in (1, 2):
        for lo in range(n):
            for hi in range(lo, n):
                carry = [0] * n
                for i in range(lo, hi + 1):
                    carry[i] = height
                yield f"[{lo},{hi}]:h{height}", tuple(carry)


def affine_mod14_vectors(base: tuple[int, ...]):
    # Small coherent laws, including identity, shifts, doubling, and negation.
    laws = (
        (1, 0),
        (1, 1),
        (1, 7),
        (2, 0),
        (2, 1),
        (3, 0),
        (5, 0),
        (13, 0),
        (13, 1),
    )
    seen: set[tuple[int, ...]] = set()
    for a, b in laws:
        carry = tuple((a * (v % N) + b) % N for v in base)
        if carry in seen or not any(carry):
            continue
        seen.add(carry)
        yield f"k={a}r+{b} mod14", carry


def gather_probes() -> list[Probe]:
    probes: list[Probe] = []
    families = (
        ("single_apex", single_apex_vectors),
        ("apex_plus_one", apex_plus_one_vectors),
        ("interval_block", interval_block_vectors),
        ("affine_mod14", affine_mod14_vectors),
    )
    for base_name, base in BASES:
        for family, builder in families:
            for label, carry in builder(base):
                probes.append(make_probe(family, base_name, label, base, carry))
    return probes


def summarize_family(probes: list[Probe]) -> None:
    print("A. Exact p0 audit over apex-debt and coherent carry families")
    by_family: defaultdict[tuple[str, str], list[Probe]] = defaultdict(list)
    for probe in probes:
        by_family[(probe.base_name, probe.family)].append(probe)

    print("  base   family          rows walls positive multiples min_p0       min_label")
    for key in sorted(by_family):
        rows = by_family[key]
        walls = [p for p in rows if p.p0 == 0]
        positives = [p for p in rows if p.p0 > 0]
        multiples = [p for p in rows if p.multiple]
        min_probe = min(positives, key=lambda p: (p.p0, p.label)) if positives else None
        print(
            f"  {key[0]:6s} {key[1]:14s} {len(rows):4d} {len(walls):5d} "
            f"{len(positives):8d} {len(multiples):9d} "
            f"{fmt_frac(min_probe.p0) if min_probe else '-':>12s} "
            f"{min_probe.label if min_probe else '-'}"
        )
    print()


def summarize_multiple_branch(probes: list[Probe]) -> None:
    print("B. Primitive multiple branch readout")
    primitive_multiple = [p for p in probes if p.multiple and p.gcd_value == 1]
    walls = [p for p in primitive_multiple if p.p0 == 0]
    positives = [p for p in primitive_multiple if p.p0 > 0]
    print(f"  primitive multiple probes={len(primitive_multiple)}")
    print(f"  primitive multiple walls={len(walls)}")
    print(f"  primitive multiple positive={len(positives)}")
    if positives:
        min_probe = min(positives, key=lambda p: (p.p0, p.label))
        print(
            "  minimum primitive-multiple p0: "
            f"{fmt_frac(min_probe.p0)} via {min_probe.base_name}/{min_probe.family}/{min_probe.label}"
        )
        print(f"    row={min_probe.row}")
        print(
            f"    components={min_probe.components}; n_clock={min_probe.n_clock}; "
            f"gcd_scaled_clock={min_probe.gcd_clock}; gcd={min_probe.gcd_value}"
        )

    print("  smallest primitive-multiple p0 rows:")
    for probe in sorted(positives, key=lambda p: (p.p0, p.base_name, p.family, p.label))[:10]:
        print(
            f"    {probe.base_name:5s} {probe.family:14s} "
            f"p0={fmt_frac(probe.p0):>14s} nclk={probe.n_clock} "
            f"gclk={probe.gcd_clock} label={probe.label}"
        )
    print()


def summarize_endpoint_split(probes: list[Probe]) -> None:
    print("C. Endpoint split: no-multiple versus apex debt")
    hist = Counter((p.multiple, p.gcd_value, p.n_clock, p.gcd_clock, p.route) for p in probes)
    print("  top endpoint buckets (multiple,gcd,n_clock,gcd_clock,route):")
    for key, count in hist.most_common(12):
        print(f"    {key}: {count}")
    print()
    print(
        "  Interpretation: no-multiple rows retain the six n-clock witnesses; "
        "primitive multiple rows lose them and must be discharged by positive "
        "measure or a new hidden endpoint.  No hidden p0=0 primitive-multiple "
        "endpoint appears in this audit."
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


def summarize_tournament() -> None:
    print("D. Tournament Analysis over improved proof filters")
    vertices = (
        RouteVertex("primitive_multiple_positive_measure", (5, 5, 2, 4, 3, 5)),
        RouteVertex("apex_congruence_debt", (5, 5, 3, 5, 4, 4)),
        RouteVertex("owner_private_derivative", (4, 4, 3, 5, 4, 4)),
        RouteVertex("gcd_scaled_endpoint_wall", (5, 3, 5, 3, 3, 4)),
        RouteVertex("no_multiple_n_clock", (5, 1, 5, 2, 5, 3)),
        RouteVertex("raw_res27_shadow", (3, 2, 2, 1, 5, 2)),
        RouteVertex("raw_first_moment", (1, 1, 1, 1, 5, 1)),
    )
    edges = tournament_edges(vertices)
    hist = Counter(len(edges[v.name]) for v in vertices)
    ordered = sorted(vertices, key=lambda v: (-len(edges[v.name]), v.name))
    print("  vertices=proof filters, not runners")
    print("  observable=(exactness,multiple relevance,endpoint handling,side-channel,compression,action)")
    print(f"  score_hist={dict(sorted(hist.items()))}")
    print(f"  directed_3cycles={count_directed_triangles(vertices, edges)}")
    print(f"  scc_sizes={scc_sizes(vertices, edges)}")
    print(f"  hamiltonian_paths={count_hamiltonian_paths(vertices, edges)}")
    print("  outdegree order:")
    for v in ordered:
        print(f"    {v.name:36s} out={len(edges[v.name])} scores={v.scores}")
    print()


def main() -> None:
    print("=" * 78)
    print("S677 LRC14 apex-debt Lebesgue sieve")
    print("=" * 78)
    print(f"n={N}; C={C}; delta={fmt_frac(DELTA)}")
    print("Apex debt identity: for v=r+27k, 14|v iff k == r (mod 14).")
    print()

    probes = gather_probes()
    summarize_family(probes)
    summarize_multiple_branch(probes)
    summarize_endpoint_split(probes)
    summarize_tournament()

    print("=" * 78)
    print("Conclusion")
    print(
        "  The no-new-wall target should be split.  No-multiple p0=0 rows are "
        "harmless n-clock endpoint walls; nonprimitive gcd-scaled rows such as "
        "2AP have scaled endpoint walls.  The real residual is primitive apex "
        "debt: rows with a multiple of 14 and gcd 1.  Across the exact S677 "
        "apex/coherent carry families, every primitive apex-debt row has "
        "positive Lebesgue safe measure.  A future proof should therefore aim "
        "directly at primitive-multiple positive measure, using apex congruence "
        "and owner-private deletion as the retained side channels."
    )


if __name__ == "__main__":
    main()
