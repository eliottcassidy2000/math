#!/usr/bin/env python3
"""Safe-component stalk descent scout for LRC14.

HYP-3026 fuses many sidecars into a large proof-interface packet.  This pass
asks whether a smaller local object can carry the same information on the
first hard automatic fiber:

    the stalk of one strict safe component.

The stalk is attached to the largest strict safe interval.  It records endpoint
owners, peak bottleneck owners, exact or coarse length/height data, and the
open/boundary status.  The tested vertices are proof carriers, not runners.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
from pathlib import Path
import argparse
import os
import sys


REPO = Path(__file__).resolve().parents[1]
TARGET_WORD = "MFCMMCCFFFCCC"
THRESHOLD = Fraction(1, 14)


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


lp = load_module(
    "lrc14_labelled_packet_classifier_s193",
    REPO / "04-computation" / "lrc14_labelled_packet_counterexample_classifier_codex_20260624.py",
)
s189 = load_module(
    "lrc14_carrier_fusion_s193",
    REPO / "04-computation" / "lrc14_carrier_fusion_switchboard_codex_s189.py",
)


def fmt(fr: Fraction) -> str:
    return str(fr.numerator) if fr.denominator == 1 else f"{fr.numerator}/{fr.denominator}"


def status(packet) -> str:
    return "boundary" if packet.strict_safe_mu == 0 else "open"


def residue_mfc_pairs(packet) -> tuple[tuple[int, str], ...]:
    return tuple(sorted((speed % 14, lp.automatic_letter(speed)) for speed in packet.speeds))


def terminal_word(packet) -> tuple[tuple[str, str], ...]:
    return tuple(
        (
            "Fok" if lp.is_fibbinary(speed) else "Fdead",
            "Mok" if lp.is_moser(speed) else "Mdead",
        )
        for speed in packet.speeds
    )


def residue_terminal(packet):
    return packet.automatic_word, residue_mfc_pairs(packet), terminal_word(packet)


def magnitude_key(packet) -> tuple[Fraction, int, int, Fraction]:
    return packet.M, packet.q_threshold, packet.farey_excess, packet.lacunary_tail_ratio


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def signed_residue(speed: int, t: Fraction, height: Fraction) -> str:
    y = frac_part(Fraction(speed) * t)
    if y <= Fraction(1, 2):
        side = "+"
    else:
        side = "-"
    return f"{speed % 14}{side}"


def threshold_owner_word(speeds: tuple[int, ...], t: Fraction) -> tuple[str, ...]:
    return tuple(
        sorted(
            signed_residue(speed, t, THRESHOLD)
            for speed in speeds
            if s189.circle_distance(speed * t) == THRESHOLD
        )
    )


def active_owner_word(speeds: tuple[int, ...], t: Fraction, height: Fraction) -> tuple[str, ...]:
    return tuple(
        sorted(
            signed_residue(speed, t, height)
            for speed in speeds
            if s189.circle_distance(speed * t) == height
        )
    )


def boundary_owner_word(speeds: tuple[int, ...]) -> tuple[tuple[str, ...], ...]:
    words = []
    for point, active in s189.boundary_safe_points(speeds):
        words.append(
            tuple(
                sorted(
                    signed_residue(speed, point, THRESHOLD)
                    for speed in active
                )
            )
        )
    return tuple(sorted(words))


def zero_sum_boundary_pairs(speeds: tuple[int, ...]) -> int:
    count = 0
    for _point, active in s189.boundary_safe_points(speeds):
        for a, b in combinations(active, 2):
            if (a + b) % 14 == 0:
                count += 1
    return count


def peak_on_component(
    speeds: tuple[int, ...],
    component: tuple[Fraction, Fraction],
) -> tuple[Fraction, Fraction, tuple[str, ...]]:
    lo, hi = component
    pts = {lo, hi}
    for point in s189.breakpoints(speeds):
        if lo <= point <= hi:
            pts.add(point)
    ordered = sorted(pts)
    best_height = Fraction(-1)
    best_t = lo
    best_owners: tuple[str, ...] = ()
    for left, right in zip(ordered, ordered[1:]):
        if left == right:
            continue
        affines = s189.cell_affines(speeds, left, right)
        candidates = {left, right}
        for a, b in combinations(affines, 2):
            if a.slope == b.slope:
                continue
            t = Fraction(b.intercept - a.intercept, a.slope - b.slope)
            if left <= t <= right:
                candidates.add(t)
        for t in candidates:
            height = min(s189.circle_distance(speed * t) for speed in speeds)
            owners = active_owner_word(speeds, t, height)
            if height > best_height or (height == best_height and t < best_t):
                best_height = height
                best_t = t
                best_owners = owners
    return best_height, best_t, best_owners


def length_bucket(length: Fraction) -> str:
    if length < Fraction(1, 2000):
        return "tiny<1/2000"
    if length < Fraction(1, 500):
        return "small<1/500"
    if length < Fraction(1, 100):
        return "medium<1/100"
    return "wide>=1/100"


def height_bucket(height: Fraction) -> str:
    if height == Fraction(3, 41):
        return "h=3/41"
    if height == Fraction(2, 27):
        return "h=2/27"
    if height >= Fraction(1, 12):
        return "h>=1/12"
    if height >= Fraction(7, 89):
        return "h>=7/89"
    return "h-mid"


@dataclass(frozen=True)
class Stalk:
    status: str
    component_count: int
    length: Fraction
    left: Fraction | None
    right: Fraction | None
    peak_time: Fraction | None
    peak_height: Fraction
    left_owners: tuple[str, ...]
    peak_owners: tuple[str, ...]
    right_owners: tuple[str, ...]
    boundary_owners: tuple[tuple[str, ...], ...]
    boundary_zero_pairs: int

    @property
    def exact_key(self) -> tuple[object, ...]:
        if self.status == "boundary":
            return (
                "boundary",
                self.component_count,
                self.boundary_owners,
                self.boundary_zero_pairs,
            )
        return (
            "open",
            self.component_count,
            self.length,
            self.left_owners,
            self.peak_owners,
            self.right_owners,
            self.peak_height,
        )

    @property
    def coarse_key(self) -> tuple[object, ...]:
        if self.status == "boundary":
            return (
                "boundary",
                self.component_count,
                self.boundary_owners,
                self.boundary_zero_pairs,
            )
        return (
            "open",
            self.component_count,
            length_bucket(self.length),
            self.left_owners,
            self.peak_owners,
            self.right_owners,
            height_bucket(self.peak_height),
        )

    @property
    def owner_only_key(self) -> tuple[object, ...]:
        if self.status == "boundary":
            return ("boundary", self.boundary_owners, self.boundary_zero_pairs)
        return ("open", self.left_owners, self.peak_owners, self.right_owners)

    def short(self) -> str:
        if self.status == "boundary":
            return f"boundary|B={len(self.boundary_owners)}|Z={self.boundary_zero_pairs}"
        return (
            f"open|bars={self.component_count}|len={fmt(self.length)}|"
            f"owners={self.left_owners}->{self.peak_owners}->{self.right_owners}|"
            f"h={fmt(self.peak_height)}"
        )


def largest_component_stalk(packet) -> Stalk:
    components = s189.safe_components(packet.speeds)
    if not components:
        return Stalk(
            status="boundary",
            component_count=0,
            length=Fraction(0),
            left=None,
            right=None,
            peak_time=None,
            peak_height=THRESHOLD,
            left_owners=(),
            peak_owners=(),
            right_owners=(),
            boundary_owners=boundary_owner_word(packet.speeds),
            boundary_zero_pairs=zero_sum_boundary_pairs(packet.speeds),
        )

    component = max(components, key=lambda gap: (gap[1] - gap[0], -gap[0]))
    peak_height, peak_time, peak_owners = peak_on_component(packet.speeds, component)
    return Stalk(
        status="open",
        component_count=len(components),
        length=component[1] - component[0],
        left=component[0],
        right=component[1],
        peak_time=peak_time,
        peak_height=peak_height,
        left_owners=threshold_owner_word(packet.speeds, component[0]),
        peak_owners=peak_owners,
        right_owners=threshold_owner_word(packet.speeds, component[1]),
        boundary_owners=(),
        boundary_zero_pairs=0,
    )


@dataclass(frozen=True)
class Row:
    packet: object
    stalk: Stalk


@dataclass(frozen=True)
class Carrier:
    name: str
    description: str
    key_func: object
    topology: int
    owner_data: int
    exact_local: int
    exact_magnitude: int
    fusion_size: int
    proof_cost: int


def carriers() -> tuple[Carrier, ...]:
    return (
        Carrier(
            "automatic_word",
            "raw Moser/fibbinary/carry word",
            lambda row: (row.packet.automatic_word,),
            0,
            0,
            0,
            0,
            1,
            1,
        ),
        Carrier(
            "residue_terminal_fiber",
            "automatic word plus mod-14 terminal-state fiber",
            lambda row: residue_terminal(row.packet),
            0,
            1,
            0,
            0,
            2,
            2,
        ),
        Carrier(
            "owner_only_stalk",
            "largest-component endpoint and peak owner residues only",
            lambda row: (residue_terminal(row.packet), row.stalk.owner_only_key),
            3,
            4,
            0,
            0,
            3,
            2,
        ),
        Carrier(
            "coarse_component_stalk",
            "largest-component owner residues plus coarse length/height buckets",
            lambda row: (residue_terminal(row.packet), row.stalk.coarse_key),
            4,
            4,
            2,
            0,
            4,
            3,
        ),
        Carrier(
            "exact_component_stalk",
            "one exact safe-component stalk, without full barcode or magnitude cocycle",
            lambda row: (residue_terminal(row.packet), row.stalk.exact_key),
            5,
            5,
            5,
            0,
            4,
            4,
        ),
        Carrier(
            "magnitude_cocycle",
            "exact M/q/Farey/lacunary magnitude cocycle",
            lambda row: (row.packet.automatic_word, magnitude_key(row.packet)),
            1,
            1,
            0,
            5,
            3,
            3,
        ),
        Carrier(
            "stalk_plus_magnitude",
            "exact component stalk with magnitude fallback",
            lambda row: (row.stalk.exact_key, magnitude_key(row.packet)),
            5,
            5,
            5,
            5,
            5,
            5,
        ),
    )


@dataclass(frozen=True)
class Report:
    carrier: Carrier
    fibers: int
    mixed_status: int
    mixed_route: int
    mixed_family: int
    max_fiber: int
    max_mixed: int
    mixed_rows: int
    largest_key: object
    largest_rows: tuple[Row, ...]

    @property
    def route_purity(self) -> int:
        if self.fibers == 0:
            return 0
        return (1000 * (self.fibers - self.mixed_route)) // self.fibers

    @property
    def status_purity(self) -> int:
        if self.fibers == 0:
            return 0
        return (1000 * (self.fibers - self.mixed_status)) // self.fibers

    @property
    def vector(self) -> tuple[int, ...]:
        return (
            self.route_purity,
            self.status_purity,
            -self.max_mixed,
            self.carrier.topology,
            self.carrier.owner_data,
            self.carrier.exact_local,
            5 - self.carrier.exact_magnitude,
            6 - self.carrier.fusion_size,
            6 - self.carrier.proof_cost,
        )


def carrier_report(rows: list[Row], carrier: Carrier) -> Report:
    groups: dict[object, list[Row]] = defaultdict(list)
    for row in rows:
        groups[carrier.key_func(row)].append(row)

    mixed_status = []
    mixed_route = []
    mixed_family = []
    for key, group in groups.items():
        if len({status(row.packet) for row in group}) > 1:
            mixed_status.append((key, group))
        if len({row.packet.route for row in group}) > 1:
            mixed_route.append((key, group))
        if len({row.packet.family for row in group}) > 1:
            mixed_family.append((key, group))

    if mixed_route:
        largest_key, largest_rows = max(mixed_route, key=lambda item: (len(item[1]), str(item[0])))
    else:
        largest_key, largest_rows = None, []
    return Report(
        carrier=carrier,
        fibers=len(groups),
        mixed_status=len(mixed_status),
        mixed_route=len(mixed_route),
        mixed_family=len(mixed_family),
        max_fiber=max((len(group) for group in groups.values()), default=0),
        max_mixed=max((len(group) for _, group in mixed_route), default=0),
        mixed_rows=sum(len(group) for _, group in mixed_route),
        largest_key=largest_key,
        largest_rows=tuple(largest_rows),
    )


def route_counts(rows: list[Row]) -> str:
    return ", ".join(f"{key}:{value}" for key, value in Counter(row.packet.route for row in rows).most_common())


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge")
    print("  considered vertices:")
    print("    runners, gaps, all bars, single bars, endpoints, peak bottlenecks,")
    print("    wall-crossing events, residues, cover arcs, Cech nerves, Fourier modes,")
    print("    matroid cocircuits, and proof obligations.")
    print("  chosen vertices:")
    print("    proof carriers built from the largest safe-component stalk.")
    print("  preserved LRC predicate:")
    print("    strict-open status and route purity on the target automatic fiber when")
    print("    the exact local stalk is retained.")
    print("  destroyed information:")
    print("    all non-largest bars, global barcode multiplicity beyond the count, and")
    print("    exact magnitude unless the stalk theorem reconstructs it locally.")
    print("  challenged assumption:")
    print("    the fusion packet may not need every sidecar as a primitive field; some")
    print("    fields can descend from a local safe-component stalk.")
    print()


def print_reports(reports: list[Report]) -> None:
    print("[1] Stalk descent ladder")
    print("  carrier                       fibers mixS mixR mixF maxF maxMix mixRows purity")
    for report in reports:
        purity = 100.0
        if report.fibers:
            purity = 100.0 * (report.fibers - report.mixed_route) / report.fibers
        print(
            "  {name:29s} {fibers:6d} {mixs:4d} {mixr:4d} {mixf:4d} "
            "{maxf:4d} {maxm:6d} {mixrows:7d} {purity:5.1f}%".format(
                name=report.carrier.name,
                fibers=report.fibers,
                mixs=report.mixed_status,
                mixr=report.mixed_route,
                mixf=report.mixed_family,
                maxf=report.max_fiber,
                maxm=report.max_mixed,
                mixrows=report.mixed_rows,
                purity=purity,
            )
        )
    print()


def print_coarse_residual(report: Report) -> None:
    print("[2] Coarse-stalk residuals")
    if report.mixed_route == 0:
        print("  coarse stalk is route-pure on this target fiber.")
        print()
        return

    groups: dict[object, list[Row]] = defaultdict(list)
    for row in report.largest_rows:
        groups[report.carrier.key_func(row)].append(row)

    print(
        f"  coarse_component_stalk leaves {report.mixed_route} mixed route fibers, "
        f"largest mixed size {report.max_mixed}, mixed rows {report.mixed_rows}."
    )

    all_groups: dict[object, list[Row]] = defaultdict(list)
    for row in ALL_ROWS:
        all_groups[report.carrier.key_func(row)].append(row)
    examples = [
        (key, group)
        for key, group in all_groups.items()
        if len({row.packet.route for row in group}) > 1
    ]
    examples.sort(key=lambda item: (-len(item[1]), str(item[0])))
    for idx, (key, group) in enumerate(examples[:5], 1):
        labels = ", ".join(row.packet.name for row in group[:6])
        print(f"  residual {idx}: size={len(group)} routes={route_counts(group)}")
        print(f"    rows={labels}")
        print(f"    key={key}")
    print("  Readout: the residuals are open-route collisions, not boundary/open leaks.")
    print()


def print_exact_stalk_examples(rows: list[Row]) -> None:
    print("[3] Exact stalk examples")
    picks = []
    wanted = (
        "AP",
        "GW",
        "12->36",
        "12->84",
        "13->104",
        "13->118",
    )
    for row in rows:
        if any(token in row.packet.name for token in wanted):
            picks.append(row)
    seen = set()
    for row in picks[:10]:
        if row.packet.name in seen:
            continue
        seen.add(row.packet.name)
        print(f"  {row.packet.name:30s} route={row.packet.route:24s} {row.stalk.short()}")
    print()


def orient(a: Report, b: Report) -> str:
    av = a.vector
    bv = b.vector
    awins = sum(x > y for x, y in zip(av, bv))
    bwins = sum(x < y for x, y in zip(av, bv))
    if awins != bwins:
        return a.carrier.name if awins > bwins else b.carrier.name
    names = [carrier.name for carrier in carriers()]
    return a.carrier.name if names.index(a.carrier.name) < names.index(b.carrier.name) else b.carrier.name


def tournament_fingerprint(reports: list[Report]) -> tuple[dict[int, int], int, list[int], int, list[str]]:
    names = [report.carrier.name for report in reports]
    edges: dict[tuple[str, str], bool] = {}
    scores = Counter()
    for a, b in combinations(reports, 2):
        winner = orient(a, b)
        loser = b.carrier.name if winner == a.carrier.name else a.carrier.name
        edges[(winner, loser)] = True
        scores[winner] += 1
        scores.setdefault(loser, 0)

    def beats(x: str, y: str) -> bool:
        return (x, y) in edges

    c3 = 0
    for a, b, c in combinations(names, 3):
        if (
            beats(a, b)
            and beats(b, c)
            and beats(c, a)
            or beats(a, c)
            and beats(c, b)
            and beats(b, a)
        ):
            c3 += 1

    graph: dict[str, list[str]] = defaultdict(list)
    reverse: dict[str, list[str]] = defaultdict(list)
    for winner, loser in edges:
        graph[winner].append(loser)
        reverse[loser].append(winner)
    seen: set[str] = set()
    order: list[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in graph[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for name in names:
        if name not in seen:
            dfs(name)

    seen.clear()
    scc: list[int] = []

    def rdfs(v: str) -> int:
        seen.add(v)
        size = 1
        for w in reverse[v]:
            if w not in seen:
                size += rdfs(w)
        return size

    for name in reversed(order):
        if name not in seen:
            scc.append(rdfs(name))

    hp = 0
    first_path: tuple[str, ...] | None = None
    for perm in permutations(names):
        if all(beats(perm[i], perm[i + 1]) for i in range(len(perm) - 1)):
            hp += 1
            if first_path is None:
                first_path = perm

    ranking = sorted(names, key=lambda name: scores[name], reverse=True)
    if first_path is not None:
        ranking = list(first_path)
    return dict(sorted(Counter(scores.values()).items())), c3, sorted(scc, reverse=True), hp, ranking


def print_tournament(reports: list[Report]) -> None:
    hist, c3, scc, hp, ranking = tournament_fingerprint(reports)
    print("[4] Tournament Analysis")
    print("  vertices_are=proof carriers and local stalks, not runners")
    print("  observable=route_purity,status_purity,max_mixed,topology,owner_data,")
    print("             exact_local_geometry,avoid_exact_magnitude,small_fusion_size,proof_cost")
    print("  switch=majority comparison of observable vectors")
    print("  tie_hamiltonian_path=" + " > ".join(carrier.name for carrier in carriers()))
    print(f"  score_hist={hist}")
    print(f"  directed_3cycles={c3}")
    print(f"  scc_sizes={scc}")
    print(f"  hamiltonian_path_count={hp}")
    print("  first_hamiltonian_path=" + " > ".join(ranking))
    print()


def print_interpretation(reports: list[Report]) -> None:
    lookup = {report.carrier.name: report for report in reports}
    residue = lookup["residue_terminal_fiber"]
    coarse = lookup["coarse_component_stalk"]
    exact = lookup["exact_component_stalk"]
    magnitude = lookup["magnitude_cocycle"]
    print("[5] Interpretation")
    print(
        f"  residue-terminal fibers have {residue.mixed_route} mixed route fibers "
        f"with max mixed size {residue.max_mixed}."
    )
    print(
        f"  A coarse largest-component stalk contracts this to {coarse.mixed_route} "
        f"mixed route fibers with max mixed size {coarse.max_mixed}."
    )
    print(
        f"  The exact largest-component stalk has {exact.mixed_route} mixed route "
        f"fibers; exact magnitude has {magnitude.mixed_route}."
    )
    print("  This suggests a descent lemma: derive barcode/normal-fan/Cech sidecars")
    print("  from the stalk of one strict safe component whenever the automatic fiber")
    print("  is already in the AP/GW target word.")


ALL_ROWS: list[Row] = []


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--single-limit", type=int, default=180)
    parser.add_argument("--two-swap-limit", type=int, default=36)
    parser.add_argument("--alias-depth", type=int, default=4)
    parser.add_argument("--lcm-tail-max", type=int, default=5)
    parser.add_argument("--target-word", default=TARGET_WORD)
    parser.add_argument("--workers", type=int, default=max(1, min(os.cpu_count() or 1, 8)))
    args = parser.parse_args()

    bank = lp.build_bank(args.single_limit, args.two_swap_limit, args.alias_depth, args.lcm_tail_max)
    target = [row for row in bank if lp.automatic_word(row[2]) == args.target_word]
    packets = lp.compute_packets(target, args.workers)

    global ALL_ROWS
    ALL_ROWS = [Row(packet=packet, stalk=largest_component_stalk(packet)) for packet in packets]
    reports = [carrier_report(ALL_ROWS, carrier) for carrier in carriers()]

    print("=== LRC14 safe-component stalk descent S193 ===")
    print(
        f"bank=HYP-2963 target automatic word {args.target_word} "
        f"single_limit={args.single_limit} two_swap_limit={args.two_swap_limit} "
        f"alias_depth={args.alias_depth} lcm_tail_max={args.lcm_tail_max}"
    )
    print(f"candidate_rows={len(target)} of {len(bank)}")
    print(f"packets={len(packets)}")
    print()
    print_assumption_challenge()
    print_reports(reports)
    print_coarse_residual(next(report for report in reports if report.carrier.name == "coarse_component_stalk"))
    print_exact_stalk_examples(ALL_ROWS)
    print_tournament(reports)
    print_interpretation(reports)


if __name__ == "__main__":
    main()
