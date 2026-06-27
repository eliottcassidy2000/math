#!/usr/bin/env python3
"""Creative LRC14 angle atlas.

This is a synthesis computation, not a proof.  It takes several speculative
LRC angles seriously enough to attach exact packet data to them: Cech nerves,
tropical slack, CRT/solenoid charts, endpoint chip-firing currents,
danger-count distributions, tope/cocircuit walls, and automaton/divisor
sidecars.

Tournament Analysis is over proof lenses and packet side channels, not runners.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, permutations
from math import gcd


THRESHOLD = Fraction(1, 14)
MAX_CRT_DENOM = 140


@dataclass(frozen=True)
class NamedRow:
    name: str
    speeds: tuple[int, ...]
    route: str


ROWS = (
    NamedRow("AP13", tuple(range(1, 14)), "AP/GW boundary control"),
    NamedRow("GW_12_to_24", tuple(list(range(1, 12)) + [13, 24]), "AP/GW boundary control"),
    NamedRow("K33_12_to_36", tuple(list(range(1, 12)) + [13, 36]), "K33 state-lift seed"),
    NamedRow("petal_10_to_20", tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 20]), "C27/unit-petal seed"),
    NamedRow("petal_13_to_26", tuple(list(range(1, 13)) + [26]), "C27/unit-petal seed"),
    NamedRow("P10_plus_GW", tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 20, 24]), "petal plus GW splice"),
    NamedRow("cover_12_to_84", tuple(list(range(1, 12)) + [13, 84]), "covering/lcm-tail seed"),
)


def fmt(frac: Fraction) -> str:
    return str(frac.numerator) if frac.denominator == 1 else f"{frac.numerator}/{frac.denominator}"


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def circle_distance(x: Fraction) -> Fraction:
    y = frac_part(x)
    return min(y, 1 - y)


def danger_intervals(speed: int) -> list[tuple[Fraction, Fraction]]:
    radius = THRESHOLD / speed
    out: list[tuple[Fraction, Fraction]] = []
    for m in range(speed):
        center = Fraction(m, speed)
        lo = center - radius
        hi = center + radius
        if lo < 0:
            out.append((lo + 1, Fraction(1, 1)))
            out.append((Fraction(0, 1), hi))
        elif hi > 1:
            out.append((lo, Fraction(1, 1)))
            out.append((Fraction(0, 1), hi - 1))
        else:
            out.append((lo, hi))
    return out


def merge_intervals(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged: list[tuple[Fraction, Fraction]] = []
    cur_lo, cur_hi = intervals[0]
    for lo, hi in intervals[1:]:
        if lo <= cur_hi:
            cur_hi = max(cur_hi, hi)
        else:
            merged.append((cur_lo, cur_hi))
            cur_lo, cur_hi = lo, hi
    merged.append((cur_lo, cur_hi))
    return merged


def safe_components(speeds: tuple[int, ...]) -> list[tuple[Fraction, Fraction]]:
    intervals: list[tuple[Fraction, Fraction]] = []
    for speed in speeds:
        intervals.extend(danger_intervals(speed))
    danger = merge_intervals(intervals)
    gaps: list[tuple[Fraction, Fraction]] = []
    cur = Fraction(0, 1)
    for lo, hi in danger:
        if lo > cur:
            gaps.append((cur, lo))
        cur = max(cur, hi)
    if cur < 1:
        gaps.append((cur, Fraction(1, 1)))
    return gaps


def safe_measure(components: list[tuple[Fraction, Fraction]]) -> Fraction:
    return sum((hi - lo for lo, hi in components), Fraction(0, 1))


def all_endpoints(speeds: tuple[int, ...]) -> list[Fraction]:
    pts = {Fraction(0, 1), Fraction(1, 1)}
    for speed in speeds:
        for lo, hi in danger_intervals(speed):
            pts.add(lo)
            pts.add(hi)
    return sorted(pts)


def boundary_safe_points(speeds: tuple[int, ...]) -> list[tuple[Fraction, tuple[int, ...]]]:
    out: list[tuple[Fraction, tuple[int, ...]]] = []
    for point in all_endpoints(speeds):
        if point == 1:
            continue
        distances = {speed: circle_distance(speed * point) for speed in speeds}
        if all(dist >= THRESHOLD for dist in distances.values()):
            active = tuple(speed for speed, dist in distances.items() if dist == THRESHOLD)
            out.append((point, active))
    return out


def danger_count_distribution(speeds: tuple[int, ...]) -> dict[int, Fraction]:
    pts = all_endpoints(speeds)
    hist: dict[int, Fraction] = defaultdict(Fraction)
    for a, b in zip(pts, pts[1:]):
        if a == b:
            continue
        mid = (a + b) / 2
        count = sum(1 for speed in speeds if circle_distance(speed * mid) < THRESHOLD)
        hist[count] += b - a
    return dict(sorted(hist.items()))


def first_unit_witness(speeds: tuple[int, ...], max_denom: int = MAX_CRT_DENOM) -> tuple[int | None, int | None, int]:
    first_d = None
    first_a = None
    first_count = 0
    for denom in range(2, max_denom + 1):
        hits = []
        for residue in range(1, denom):
            if gcd(residue, denom) != 1:
                continue
            t = Fraction(residue, denom)
            if all(circle_distance(speed * t) >= THRESHOLD for speed in speeds):
                hits.append(residue)
        if hits:
            first_d = denom
            first_a = hits[0]
            first_count = len(hits)
            break
    return first_d, first_a, first_count


def largest_component(components: list[tuple[Fraction, Fraction]]) -> tuple[Fraction, Fraction] | None:
    if not components:
        return None
    return max(components, key=lambda gap: (gap[1] - gap[0], -gap[0]))


def midpoint_slack(speeds: tuple[int, ...], component: tuple[Fraction, Fraction] | None) -> Fraction:
    if component is None:
        return Fraction(0, 1)
    mid = (component[0] + component[1]) / 2
    return min(circle_distance(speed * mid) - THRESHOLD for speed in speeds)


def zero_sum_active_pair_count(boundary: list[tuple[Fraction, tuple[int, ...]]]) -> int:
    count = 0
    for _point, active in boundary:
        for a, b in combinations(active, 2):
            if (a + b) % 14 == 0:
                count += 1
    return count


def distribution_word(hist: dict[int, Fraction]) -> str:
    return ",".join(f"N{k}:{fmt(v)}" for k, v in hist.items())


def row_readout() -> None:
    print("Exact named-row creative-angle readout")
    print(
        "name route safe_mu components largest_component midpoint_slack "
        "boundary_safe_points zero_sum_active_pairs first_unit_witness danger_count_distribution"
    )
    for row in ROWS:
        components = safe_components(row.speeds)
        boundary = boundary_safe_points(row.speeds)
        largest = largest_component(components)
        first_d, first_a, first_count = first_unit_witness(row.speeds)
        witness = "NA" if first_d is None else f"{first_a}/{first_d}({first_count})"
        print(
            row.name,
            row.route,
            fmt(safe_measure(components)),
            len(components),
            "NA" if largest is None else f"{fmt(largest[0])}-{fmt(largest[1])}",
            fmt(midpoint_slack(row.speeds, largest)),
            len(boundary),
            zero_sum_active_pair_count(boundary),
            witness,
            distribution_word(danger_count_distribution(row.speeds)),
            sep=" | ",
        )
    print()


@dataclass(frozen=True)
class Lens:
    name: str
    lrc_predicate: int
    exactness: int
    topology: int
    arithmetic: int
    endpoint_owner: int
    dual_certificate: int
    computability: int
    quotient_guard: int
    destroys: str

    def vector(self) -> dict[str, int]:
        return {
            "lrc_predicate": self.lrc_predicate,
            "exactness": self.exactness,
            "topology": self.topology,
            "arithmetic": self.arithmetic,
            "endpoint_owner": self.endpoint_owner,
            "dual_certificate": self.dual_certificate,
            "computability": self.computability,
            "quotient_guard": self.quotient_guard,
        }


FIELDS = (
    "lrc_predicate",
    "exactness",
    "topology",
    "arithmetic",
    "endpoint_owner",
    "dual_certificate",
    "computability",
    "quotient_guard",
)

TIE_PATH = (
    "labelled_packet_sheaf",
    "cech_nerve_cover_homology",
    "endpoint_chip_firing_current",
    "tropical_slack_potential",
    "crt_solenoid_chart",
    "matroid_tope_cocircuit_wall",
    "danger_count_entropy_dual",
    "automaton_divisor_sidecar",
    "raw_creative_metaphor",
)


def lens_bank() -> list[Lens]:
    return [
        Lens("labelled_packet_sheaf", 4, 4, 4, 4, 4, 4, 4, 4, "nothing if packet labels are retained"),
        Lens("cech_nerve_cover_homology", 3, 3, 4, 1, 3, 2, 3, 3, "arithmetic period and qdiv labels"),
        Lens("endpoint_chip_firing_current", 3, 3, 3, 2, 4, 2, 3, 3, "global period charts unless coupled to CRT"),
        Lens("tropical_slack_potential", 3, 4, 2, 1, 2, 4, 4, 3, "boundary-owner identity at equal slack"),
        Lens("crt_solenoid_chart", 3, 4, 1, 4, 2, 2, 4, 3, "real interval topology between denominator charts"),
        Lens("matroid_tope_cocircuit_wall", 3, 2, 4, 1, 4, 3, 2, 3, "exact denominator arithmetic"),
        Lens("danger_count_entropy_dual", 2, 3, 2, 1, 1, 4, 4, 2, "endpoint owner and route labels"),
        Lens("automaton_divisor_sidecar", 2, 3, 1, 4, 1, 2, 4, 3, "circle topology and endpoint currents"),
        Lens("raw_creative_metaphor", 0, 0, 0, 0, 0, 0, 0, 0, "all LRC predicate data"),
    ]


def compare(a: Lens, b: Lens) -> str:
    va = a.vector()
    vb = b.vector()
    aw = sum(1 for field in FIELDS if va[field] > vb[field])
    bw = sum(1 for field in FIELDS if vb[field] > va[field])
    if aw > bw:
        return a.name
    if bw > aw:
        return b.name
    order = {name: i for i, name in enumerate(TIE_PATH)}
    return a.name if order[a.name] < order[b.name] else b.name


def edge_map(lenses: list[Lens]) -> dict[tuple[str, str], str]:
    edges: dict[tuple[str, str], str] = {}
    for a, b in combinations(lenses, 2):
        winner = compare(a, b)
        loser = b.name if winner == a.name else a.name
        edges[(winner, loser)] = winner
    return edges


def beats(edges: dict[tuple[str, str], str], a: str, b: str) -> bool:
    return (a, b) in edges


def score_histogram(names: list[str], edges: dict[tuple[str, str], str]) -> dict[int, int]:
    scores = {name: 0 for name in names}
    for winner, _loser in edges:
        scores[winner] += 1
    return dict(sorted(Counter(scores.values()).items()))


def directed_3cycles(names: list[str], edges: dict[tuple[str, str], str]) -> int:
    count = 0
    for a, b, c in combinations(names, 3):
        if (
            beats(edges, a, b)
            and beats(edges, b, c)
            and beats(edges, c, a)
            or beats(edges, a, c)
            and beats(edges, c, b)
            and beats(edges, b, a)
        ):
            count += 1
    return count


def scc_sizes(names: list[str], edges: dict[tuple[str, str], str]) -> list[int]:
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
    sizes: list[int] = []

    def rdfs(v: str) -> int:
        seen.add(v)
        size = 1
        for w in reverse[v]:
            if w not in seen:
                size += rdfs(w)
        return size

    for name in reversed(order):
        if name not in seen:
            sizes.append(rdfs(name))
    return sorted(sizes, reverse=True)


def hamiltonian_paths(names: list[str], edges: dict[tuple[str, str], str]) -> tuple[int, tuple[str, ...] | None]:
    count = 0
    first = None
    for perm in permutations(names):
        if all(beats(edges, perm[i], perm[i + 1]) for i in range(len(perm) - 1)):
            count += 1
            if first is None:
                first = perm
    return count, first


def tournament_readout() -> None:
    lenses = lens_bank()
    names = [lens.name for lens in lenses]
    edges = edge_map(lenses)
    hp_count, hp_first = hamiltonian_paths(names, edges)
    print("Creative Lens Tournament Analysis")
    print("vertices=proof lenses and side channels, not runners/arcs/raw sequence entries")
    print(f"pairwise_observable={','.join(FIELDS)}")
    print(f"tie_hamiltonian_path={'>'.join(TIE_PATH)}")
    print(f"score_hist={score_histogram(names, edges)}")
    print(f"directed_3cycles={directed_3cycles(names, edges)}")
    print(f"scc_sizes={scc_sizes(names, edges)}")
    print(f"hamiltonian_path_count={hp_count}")
    print(f"first_hamiltonian_path={'>'.join(hp_first) if hp_first else 'NONE'}")
    print()


def angle_manifest() -> None:
    print("Creative angle manifest")
    print("cech_nerve_cover_homology: use positive complement components and zero-open boundary cocircuits as topological packet data")
    print("tropical_slack_potential: maximize minimum distance slack, but keep ties and endpoint owners")
    print("crt_solenoid_chart: treat rational witnesses as charts on a profinite/solenoid atlas, not as one fixed denominator ladder")
    print("endpoint_chip_firing_current: boundary active-speed pairs are currents; AP/GW zero-sum pairs are calibration, not generic scalar data")
    print("danger_count_entropy_dual: danger multiplicity distribution feeds moment/LP duals, but loses route labels unless packetized")
    print("matroid_tope_cocircuit_wall: endpoint arrangements give topes and boundary cocircuits; no-tope/no-cocircuit is the residual shape")
    print("automaton_divisor_sidecar: fibbinary/Moser/divisor labels should be sidecars on exact packets, not substitute witnesses")
    print()
    print("Assumption challenge")
    print(
        "Alternate vertices considered: runners, gaps, fixed circle sections, "
        "section boundaries, wall-crossing events, residues, cover arcs, Fourier "
        "modes, matroid circuits, divisor atoms, CRT charts, Cech nerve classes, "
        "and proof obligations."
    )
    print(
        "Chosen vertices are proof lenses because the preserved predicate is the "
        "labelled LRC14 packet verdict.  The quotient intentionally destroys raw "
        "runner identity and rejects metaphor-only carriers."
    )


def main() -> None:
    print("LRC14 creative angle atlas (codex-2026-06-26-S175)")
    print("source_threads=HYP-2963,HYP-2970,HYP-2973,HYP-2974,HYP-3008,HYP-3012,HYP-3013")
    print()
    row_readout()
    tournament_readout()
    angle_manifest()


if __name__ == "__main__":
    main()
