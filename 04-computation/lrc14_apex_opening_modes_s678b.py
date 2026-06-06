#!/usr/bin/env python3
"""
S678b: local opening modes for the LRC14 primitive apex-debt branch.

S677 sharpened the remaining LRC14 target to:

    primitive apex debt => p_0(V, 1/14) > 0.

Here "apex debt" means a carried speed v=r+27k has become divisible by 14,
equivalently k == r mod 14.  This script asks for a more proof-like witness
than the scalar fact p_0>0.  For the exact S677 coherent carry probes, it
classifies why positive measure opens.

The four certificate modes are:

  1. clock_shutter:
       delete the apex speed; at a unit clock t=j/14 the non-apex row has a
       one-sided safe cone wider than the apex danger shutter radius 1/(14a).
  2. apex_free_side_door:
       the full row has a positive safe interval whose two endpoint owners are
       non-apex speeds.
  3. one_apex_hinge:
       the full row has a positive safe interval with exactly one apex-owned
       endpoint and one non-apex endpoint.
  4. apex_period_aperture:
       even the endpoint owners are both apex speeds; the witness is a literal
       gap between neighbouring apex danger slits while the non-apex row stays
       safe.

Tournament Analysis / assumption challenge:
  Vertices are proof-certificate modes, not runners.  Candidate vertex sets
  considered include runners, unit clocks, apex-deletion rows, safe intervals,
  endpoint owners, carry congruence sites, and proof obligations.  The selected
  quotient preserves "which local opening lemma certifies p_0>0" and destroys
  raw runner order.  Pairwise observables are locality, exactness, coverage,
  deletion/private-owner content, quotient stability, and proof actionability.
"""

from __future__ import annotations

import sys
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_apex_debt_lebesgue_s677 as s677  # noqa: E402
import lrc14_lebesgue_wall_s676 as s676  # noqa: E402


N = s676.N_TOTAL
DELTA = s676.DELTA
UNIT_CLOCKS = s677.UNIT_CLOCKS


def fmt_frac(x: Fraction) -> str:
    return s676.fmt_frac(x)


def remove_index(row: tuple[int, ...], index: int) -> tuple[int, ...]:
    return tuple(v for i, v in enumerate(row) if i != index)


def endpoint_owners(row: tuple[int, ...], t: Fraction) -> tuple[int, ...]:
    if t in (0, 1):
        return tuple()
    return tuple(v for v in row if s676.norm_circle(Fraction(v) * t) == DELTA)


def local_safe_cone(row: tuple[int, ...], t0: Fraction) -> tuple[bool, Fraction, Fraction]:
    """Return whether row is safe at t0 and exact one-sided local cone widths.

    Widths are measured left and right from t0 before some non-apex runner
    enters the open danger collar.  The computation is local but exact: it is
    just the linear derivative of ||v t|| on the circle at t0.
    """
    left = Fraction(10**18)
    right = Fraction(10**18)
    half = Fraction(1, 2)
    one = Fraction(1)
    for v in row:
        r = s676.frac_part(Fraction(v) * t0)
        dist = min(r, 1 - r)
        if dist < DELTA:
            return False, Fraction(0), Fraction(0)
        if r == 0:
            return False, Fraction(0), Fraction(0)
        if r < half:
            left = min(left, (r - DELTA) / v)
            right = min(right, (one - DELTA - r) / v)
        elif r == half:
            width = (half - DELTA) / v
            left = min(left, width)
            right = min(right, width)
        else:
            left = min(left, (r - DELTA) / v)
            right = min(right, (one - r - DELTA) / v)
    return True, left, right


@dataclass(frozen=True)
class ClockCert:
    apex_speed: int
    unit_clock: int
    shutter: Fraction
    left_width: Fraction
    right_width: Fraction
    surplus: Fraction

    @property
    def side(self) -> str:
        if self.left_width >= self.right_width:
            return "left"
        return "right"


@dataclass(frozen=True)
class SafeIntervalCert:
    left: Fraction
    right: Fraction
    left_owners: tuple[int, ...]
    right_owners: tuple[int, ...]

    @property
    def width(self) -> Fraction:
        return self.right - self.left

    @property
    def left_apex(self) -> bool:
        return any(v % N == 0 for v in self.left_owners)

    @property
    def right_apex(self) -> bool:
        return any(v % N == 0 for v in self.right_owners)

    @property
    def endpoint_class(self) -> str:
        if not self.left_apex and not self.right_apex:
            return "apex_free"
        if self.left_apex ^ self.right_apex:
            return "one_apex"
        return "apex_both"


@dataclass(frozen=True)
class RowAudit:
    probe: s677.Probe
    route: str
    clock_cert: ClockCert | None
    endpoint_class: str
    endpoint_cert: SafeIntervalCert
    safe_intervals: tuple[SafeIntervalCert, ...]


def clock_certificates(probe: s677.Probe) -> tuple[ClockCert, ...]:
    certs: list[ClockCert] = []
    for index, apex_speed in enumerate(probe.row):
        if apex_speed % N != 0:
            continue
        deleted = remove_index(probe.row, index)
        shutter = DELTA / apex_speed
        for j in UNIT_CLOCKS:
            t0 = Fraction(j, N)
            safe, left, right = local_safe_cone(deleted, t0)
            if not safe:
                continue
            surplus = max(left, right) - shutter
            if surplus > 0:
                certs.append(ClockCert(apex_speed, j, shutter, left, right, surplus))
    return tuple(sorted(certs, key=lambda c: (c.surplus, c.apex_speed, c.unit_clock), reverse=True))


def safe_interval_certs(row: tuple[int, ...]) -> tuple[SafeIntervalCert, ...]:
    sweep = s676.depth_sweep(row)
    certs: list[SafeIntervalCert] = []
    for left, right in sweep.safe_intervals:
        certs.append(
            SafeIntervalCert(
                left=left,
                right=right,
                left_owners=endpoint_owners(row, left),
                right_owners=endpoint_owners(row, right),
            )
        )
    return tuple(certs)


def endpoint_class(certs: tuple[SafeIntervalCert, ...]) -> str:
    if any(c.endpoint_class == "apex_free" for c in certs):
        return "apex_free"
    if any(c.endpoint_class == "one_apex" for c in certs):
        return "one_apex"
    return "apex_both_only"


def choose_endpoint_cert(endpoint_mode: str, certs: tuple[SafeIntervalCert, ...]) -> SafeIntervalCert:
    if endpoint_mode == "apex_free":
        pool = [c for c in certs if c.endpoint_class == "apex_free"]
    elif endpoint_mode == "one_apex":
        pool = [c for c in certs if c.endpoint_class == "one_apex"]
    else:
        pool = list(certs)
    return max(pool, key=lambda c: (c.width, c.left))


def audit_rows() -> tuple[RowAudit, ...]:
    probes = [p for p in s677.gather_probes() if p.multiple and p.gcd_value == 1]
    audits: list[RowAudit] = []
    for probe in probes:
        clocks = clock_certificates(probe)
        intervals = safe_interval_certs(probe.row)
        mode = endpoint_class(intervals)
        interval_cert = choose_endpoint_cert(mode, intervals)
        if clocks:
            route = "clock_shutter"
        elif mode == "apex_free":
            route = "apex_free_side_door"
        elif mode == "one_apex":
            route = "one_apex_hinge"
        else:
            route = "apex_period_aperture"
        audits.append(
            RowAudit(
                probe=probe,
                route=route,
                clock_cert=clocks[0] if clocks else None,
                endpoint_class=mode,
                endpoint_cert=interval_cert,
                safe_intervals=intervals,
            )
        )
    return tuple(audits)


def summarize_opening_modes(audits: tuple[RowAudit, ...]) -> None:
    print("A. Primitive apex-debt opening modes")
    print(f"  primitive apex-debt rows audited: {len(audits)}")
    print(f"  p0-positive rows: {sum(1 for a in audits if a.probe.p0 > 0)}")
    print(f"  p0-wall rows: {sum(1 for a in audits if a.probe.p0 == 0)}")
    print()
    print("  priority route histogram:")
    for route, count in Counter(a.route for a in audits).most_common():
        print(f"    {route:24s} {count:4d}")
    print()
    print("  endpoint owner classes, independent of clock priority:")
    for mode, count in Counter(a.endpoint_class for a in audits).most_common():
        print(f"    {mode:24s} {count:4d}")
    print()
    clocked = [a for a in audits if a.clock_cert is not None]
    print(f"  clock-shutter rows: {len(clocked)}")
    if clocked:
        min_clock = min(clocked, key=lambda a: (a.clock_cert.surplus, a.probe.p0))
        cert = min_clock.clock_cert
        assert cert is not None
        print(
            "  smallest clock surplus: "
            f"{fmt_frac(cert.surplus)} via {min_clock.probe.base_name}/"
            f"{min_clock.probe.family}/{min_clock.probe.label}"
        )
        print(
            f"    apex={cert.apex_speed}; clock={cert.unit_clock}/14; "
            f"shutter={fmt_frac(cert.shutter)}; left={fmt_frac(cert.left_width)}; "
            f"right={fmt_frac(cert.right_width)}; side={cert.side}"
        )
    print()


def summarize_endpoint_geometry(audits: tuple[RowAudit, ...]) -> None:
    print("B. Endpoint-owner geometry")
    interval_hist = Counter(c.endpoint_class for audit in audits for c in audit.safe_intervals)
    print("  safe interval endpoint-owner histogram:")
    for key, count in interval_hist.most_common():
        print(f"    {key:12s} {count:5d}")
    print()
    print("  widest representative by priority route:")
    for route in ("clock_shutter", "apex_free_side_door", "one_apex_hinge", "apex_period_aperture"):
        rows = [a for a in audits if a.route == route]
        if not rows:
            continue
        row = max(rows, key=lambda a: (a.endpoint_cert.width, a.probe.p0))
        cert = row.endpoint_cert
        clock = row.clock_cert
        print(
            f"    {route:24s} p0={fmt_frac(row.probe.p0):>16s} "
            f"width={fmt_frac(cert.width):>12s} row={row.probe.base_name}/"
            f"{row.probe.family}/{row.probe.label}"
        )
        print(
            f"      interval=({fmt_frac(cert.left)}, {fmt_frac(cert.right)}) "
            f"left_owners={cert.left_owners or '-'} right_owners={cert.right_owners or '-'}"
        )
        if clock is not None:
            print(
                f"      clock={clock.unit_clock}/14 apex={clock.apex_speed} "
                f"surplus={fmt_frac(clock.surplus)}"
            )
    print()
    print("  narrowest safe component by endpoint class:")
    by_class: defaultdict[str, list[tuple[SafeIntervalCert, RowAudit]]] = defaultdict(list)
    for audit in audits:
        for cert in audit.safe_intervals:
            by_class[cert.endpoint_class].append((cert, audit))
    for key in ("apex_free", "one_apex", "apex_both"):
        if key not in by_class:
            continue
        cert, audit = min(by_class[key], key=lambda item: (item[0].width, item[1].probe.p0))
        print(
            f"    {key:12s} min_width={fmt_frac(cert.width):>16s} "
            f"p0={fmt_frac(audit.probe.p0):>16s} "
            f"{audit.probe.base_name}/{audit.probe.family}/{audit.probe.label}"
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
    print("C. Tournament Analysis over certificate modes")
    vertices = (
        RouteVertex("clock_shutter", (5, 5, 4, 5, 4, 5)),
        RouteVertex("apex_free_side_door", (4, 5, 5, 3, 5, 4)),
        RouteVertex("one_apex_hinge", (4, 5, 4, 4, 4, 4)),
        RouteVertex("apex_period_aperture", (5, 4, 2, 3, 3, 3)),
        RouteVertex("raw_p0_scalar", (5, 1, 5, 1, 2, 2)),
        RouteVertex("raw_res27_shadow", (2, 1, 1, 1, 5, 1)),
    )
    edges = tournament_edges(vertices)
    hist = Counter(len(edges[v.name]) for v in vertices)
    ordered = sorted(vertices, key=lambda v: (-len(edges[v.name]), v.name))
    print("  vertices=certificate modes, not runners")
    print("  observable=(locality,exactness,coverage,deletion-owner,quotient-stability,action)")
    print(f"  score_hist={dict(sorted(hist.items()))}")
    print(f"  directed_3cycles={count_directed_triangles(vertices, edges)}")
    print(f"  scc_sizes={scc_sizes(vertices, edges)}")
    print(f"  hamiltonian_paths={count_hamiltonian_paths(vertices, edges)}")
    print("  outdegree order:")
    for v in ordered:
        print(f"    {v.name:24s} out={len(edges[v.name])} scores={v.scores}")
    print()


def summarize_proof_route(audits: tuple[RowAudit, ...]) -> None:
    print("D. Proof-route synthesis")
    uncovered = [a for a in audits if a.probe.p0 <= 0]
    assert not uncovered
    print(
        "  The S677 primitive apex-debt branch is not one phenomenon.  It "
        "splits into local opening modes.  Clock shutter covers rows whose "
        "deleted-apex unit-clock cone outlives the apex collar.  Rows not "
        "caught there still expose endpoint-owner intervals: an apex-free "
        "side door, a one-apex hinge, or a pure apex-period aperture."
    )
    print()
    print("  Lemma package suggested by the audit:")
    print("    L1 clock shutter: deleted-apex cone width > 1/(14a) opens p0.")
    print("    L2 side door: a safe interval with non-apex owners on both ends survives apex debt.")
    print("    L3 hinge: one apex owner plus one non-apex owner gives a transversal opening.")
    print("    L4 aperture: two apex endpoints still leave an apex-period gap if non-apex depth is zero.")
    print()
    print(
        "  This is the owner-private deletion idea in a usable form: do not "
        "ask whether the apex is globally harmless.  Ask which boundary of a "
        "specific safe interval it owns."
    )


def main() -> None:
    print("=" * 78)
    print("S678b LRC14 primitive apex-debt opening modes")
    print("=" * 78)
    print("HYP-2256 / T752")
    print(f"n={N}; C={s676.C}; delta={fmt_frac(DELTA)}")
    print("Parent route: HYP-2253 primitive apex debt positive measure.")
    print()

    audits = audit_rows()
    summarize_opening_modes(audits)
    summarize_endpoint_geometry(audits)
    summarize_tournament()
    summarize_proof_route(audits)


if __name__ == "__main__":
    main()
