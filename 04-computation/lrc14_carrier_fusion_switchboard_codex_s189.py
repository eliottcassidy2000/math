#!/usr/bin/env python3
"""LRC14 carrier-fusion switchboard scout.

This is a proof-interface computation, not a proof of LRC14.  It takes the
recent barcode, automatic-language, perfect-number/divisor, and creative
packet-lens work and asks a narrow question:

    which carrier bundles stop mixing AP/Goddyn-Wong boundary atoms with
    strictly open rows in the finite stress bank?

The row bank is deliberately small and named.  Rows are test packets; the
Tournament Analysis vertices are proof carriers, not runners.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, permutations
from math import gcd


N = 14
THRESHOLD = Fraction(1, N)
MAX_CRT_DENOM = 180


@dataclass(frozen=True)
class Row:
    name: str
    speeds: tuple[int, ...]
    route: str


def first_fibbinary(count: int) -> tuple[int, ...]:
    out: list[int] = []
    n = 1
    while len(out) < count:
        if is_fibbinary(n):
            out.append(n)
        n += 1
    return tuple(out)


def first_moser(count: int) -> tuple[int, ...]:
    out: list[int] = []
    n = 1
    while len(out) < count:
        if is_moser(n):
            out.append(n)
        n += 1
    return tuple(out)


def named_rows() -> tuple[Row, ...]:
    return (
        Row("AP13", tuple(range(1, 14)), "BOUNDARY-AP-GW"),
        Row("GW_12_to_24", tuple(list(range(1, 12)) + [13, 24]), "BOUNDARY-AP-GW"),
        Row("K33_12_to_36", tuple(list(range(1, 12)) + [13, 36]), "K33-STATE-LIFT"),
        Row("loose_12_to_26", tuple(list(range(1, 12)) + [13, 26]), "Q-WITNESS"),
        Row("loose_12_to_96", tuple(list(range(1, 12)) + [13, 96]), "Q-WITNESS"),
        Row("GW_fiber_tail_12_to_38", tuple(list(range(1, 12)) + [13, 38]), "Q-WITNESS"),
        Row("GW_fiber_tail_12_to_52", tuple(list(range(1, 12)) + [13, 52]), "Q-WITNESS"),
        Row("GW_fiber_tail_12_to_150", tuple(list(range(1, 12)) + [13, 150]), "Q-WITNESS"),
        Row("petal_10_to_20", tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 20]), "BOUNDARY-PETAL-SPORADIC"),
        Row("petal_13_to_26", tuple(list(range(1, 13)) + [26]), "BOUNDARY-PETAL-SPORADIC"),
        Row("P10_plus_GW", tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 20, 24]), "BOUNDARY-PETAL-SPORADIC"),
        Row("P10_plus_K33", tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 20, 36]), "K33-STATE-LIFT"),
        Row("cover_12_to_84", tuple(list(range(1, 12)) + [13, 84]), "COVERING-MOMENT"),
        Row("fibbinary_first13", first_fibbinary(13), "AUTOMATIC-OPEN-CONTROL"),
        Row("moser_de_bruijn_first13", first_moser(13), "AUTOMATIC-OPEN-CONTROL"),
    )


def fmt(fr: Fraction) -> str:
    return str(fr.numerator) if fr.denominator == 1 else f"{fr.numerator}/{fr.denominator}"


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def circle_distance(x: Fraction) -> Fraction:
    y = frac_part(x)
    return min(y, Fraction(1) - y)


def danger_intervals(speed: int) -> list[tuple[Fraction, Fraction]]:
    radius = THRESHOLD / speed
    out: list[tuple[Fraction, Fraction]] = []
    for m in range(speed):
        center = Fraction(m, speed)
        lo = center - radius
        hi = center + radius
        if lo < 0:
            out.append((lo + 1, Fraction(1)))
            out.append((Fraction(0), hi))
        elif hi > 1:
            out.append((lo, Fraction(1)))
            out.append((Fraction(0), hi - 1))
        else:
            out.append((lo, hi))
    return out


def merge_intervals(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    out: list[tuple[Fraction, Fraction]] = []
    cur_lo, cur_hi = intervals[0]
    for lo, hi in intervals[1:]:
        if lo <= cur_hi:
            cur_hi = max(cur_hi, hi)
        else:
            out.append((cur_lo, cur_hi))
            cur_lo, cur_hi = lo, hi
    out.append((cur_lo, cur_hi))
    return out


def safe_components(speeds: tuple[int, ...]) -> list[tuple[Fraction, Fraction]]:
    intervals: list[tuple[Fraction, Fraction]] = []
    for speed in speeds:
        intervals.extend(danger_intervals(speed))
    danger = merge_intervals(intervals)
    gaps: list[tuple[Fraction, Fraction]] = []
    cur = Fraction(0)
    for lo, hi in danger:
        if cur < lo:
            gaps.append((cur, lo))
        cur = max(cur, hi)
    if cur < 1:
        gaps.append((cur, Fraction(1)))
    return gaps


def safe_measure(components: list[tuple[Fraction, Fraction]]) -> Fraction:
    return sum((hi - lo for lo, hi in components), Fraction(0))


def all_endpoints(speeds: tuple[int, ...]) -> list[Fraction]:
    pts = {Fraction(0), Fraction(1)}
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


def zero_sum_active_pair_count(boundary: list[tuple[Fraction, tuple[int, ...]]]) -> int:
    count = 0
    for _point, active in boundary:
        for a, b in combinations(active, 2):
            if (a + b) % N == 0:
                count += 1
    return count


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


def distribution_word(hist: dict[int, Fraction]) -> str:
    return ",".join(f"N{k}:{fmt(v)}" for k, v in hist.items())


def first_unit_witness(speeds: tuple[int, ...]) -> tuple[int | None, int | None, int]:
    for denom in range(2, MAX_CRT_DENOM + 1):
        hits: list[int] = []
        for residue in range(1, denom):
            if gcd(residue, denom) != 1:
                continue
            t = Fraction(residue, denom)
            if all(circle_distance(speed * t) >= THRESHOLD for speed in speeds):
                hits.append(residue)
        if hits:
            return denom, hits[0], len(hits)
    return None, None, 0


def midpoint_slack(speeds: tuple[int, ...], component: tuple[Fraction, Fraction] | None) -> Fraction:
    if component is None:
        return Fraction(0)
    mid = (component[0] + component[1]) / 2
    return min(circle_distance(speed * mid) - THRESHOLD for speed in speeds)


def largest_component(components: list[tuple[Fraction, Fraction]]) -> tuple[Fraction, Fraction] | None:
    if not components:
        return None
    return max(components, key=lambda gap: (gap[1] - gap[0], -gap[0]))


def breakpoints(row: tuple[int, ...]) -> list[Fraction]:
    pts = {Fraction(0), Fraction(1)}
    for speed in row:
        for k in range(2 * speed + 1):
            pts.add(Fraction(k, 2 * speed))
    return sorted(pts)


@dataclass(frozen=True)
class Affine:
    slope: int
    intercept: int

    def value(self, t: Fraction) -> Fraction:
        return Fraction(self.slope) * t + self.intercept


def cell_affines(row: tuple[int, ...], left: Fraction, right: Fraction) -> list[Affine]:
    mid = (left + right) / 2
    out: list[Affine] = []
    for speed in row:
        x = Fraction(speed) * mid
        base = x.numerator // x.denominator
        residue = x - base
        if residue < Fraction(1, 2):
            out.append(Affine(speed, -base))
        else:
            out.append(Affine(-speed, base + 1))
    return out


def min_affine_value(affines: list[Affine], t: Fraction) -> Fraction:
    return min(affine.value(t) for affine in affines)


def exact_maximin(row: tuple[int, ...]) -> tuple[Fraction, Fraction]:
    best_h = Fraction(-1)
    best_t = Fraction(0)
    pts = breakpoints(tuple(sorted(row)))
    for left, right in zip(pts, pts[1:]):
        if left == right:
            continue
        affines = cell_affines(row, left, right)
        candidates = {left, right}
        for a, b in combinations(affines, 2):
            if a.slope == b.slope:
                continue
            t = Fraction(b.intercept - a.intercept, a.slope - b.slope)
            if left <= t <= right:
                candidates.add(t)
        for t in candidates:
            h = min_affine_value(affines, t)
            if h > best_h or (h == best_h and t < best_t):
                best_h, best_t = h, t
    return best_h, best_t


def bits_lsb(n: int) -> tuple[int, ...]:
    if n == 0:
        return (0,)
    out: list[int] = []
    while n:
        out.append(n & 1)
        n >>= 1
    return tuple(out)


def is_fibbinary(n: int) -> bool:
    return (n & (n >> 1)) == 0


def is_moser(n: int) -> bool:
    return all(bit == 0 for pos, bit in enumerate(bits_lsb(n)) if pos % 2 == 1)


def language_letter(n: int) -> str:
    if is_moser(n):
        return "M"
    if is_fibbinary(n):
        return "F"
    return "C"


def automatic_word(speeds: tuple[int, ...]) -> str:
    return "".join(language_letter(speed) for speed in sorted(speeds))


def doubling_transition_word(speeds: tuple[int, ...]) -> str:
    counts = Counter(f"{language_letter(v)}>{language_letter(2 * v)}" for v in speeds)
    return ",".join(f"{key}:{counts[key]}" for key in sorted(counts))


def magnitude_cocycle(speeds: tuple[int, ...]) -> tuple[str, int, int]:
    base = set(range(1, 14))
    current = set(speeds)
    dropped = tuple(sorted(base - current))
    added = tuple(sorted(current - base))
    if not dropped and not added:
        return "id", 0, 0
    word = f"drop{'.'.join(map(str, dropped)) or 'none'}_add{'.'.join(map(str, added)) or 'none'}"
    height = max(added or (0,))
    delta_sum = sum(added) - sum(dropped)
    return word, height, delta_sum


@dataclass(frozen=True)
class RowFeatures:
    row: Row
    status: str
    exact_m: Fraction
    exact_t: Fraction
    safe_mu: Fraction
    bar_count: int
    longest_bar: Fraction
    midpoint_slack: Fraction
    first_chart: str
    first_chart_den: int | None
    boundary_count: int
    zero_sum_pairs: int
    endpoint_current_word: str
    danger_word: str
    automatic_word: str
    doubling_word: str
    magnitude_word: str
    magnitude_height: int
    magnitude_delta: int
    fusion_signature: tuple[object, ...]


def row_features(row: Row) -> RowFeatures:
    components = safe_components(row.speeds)
    mu = safe_measure(components)
    status = "boundary" if mu == 0 else "open"
    m_val, m_time = exact_maximin(row.speeds)
    largest = largest_component(components)
    longest = Fraction(0) if largest is None else largest[1] - largest[0]
    slack = midpoint_slack(row.speeds, largest)
    boundary = boundary_safe_points(row.speeds)
    zero_pairs = zero_sum_active_pair_count(boundary)
    chart_den, chart_res, chart_count = first_unit_witness(row.speeds)
    chart = "NA" if chart_den is None else f"{chart_res}/{chart_den}({chart_count})"
    mag_word, mag_height, mag_delta = magnitude_cocycle(row.speeds)
    endpoint_word = f"B{len(boundary)}Z{zero_pairs}"
    danger_word = distribution_word(danger_count_distribution(row.speeds))
    auto_word = automatic_word(row.speeds)
    double_word = doubling_transition_word(row.speeds)
    fusion = (
        status,
        len(components),
        fmt(mu),
        fmt(longest),
        chart,
        endpoint_word,
        mag_word,
        auto_word,
        danger_word,
    )
    return RowFeatures(
        row=row,
        status=status,
        exact_m=m_val,
        exact_t=m_time,
        safe_mu=mu,
        bar_count=len(components),
        longest_bar=longest,
        midpoint_slack=slack,
        first_chart=chart,
        first_chart_den=chart_den,
        boundary_count=len(boundary),
        zero_sum_pairs=zero_pairs,
        endpoint_current_word=endpoint_word,
        danger_word=danger_word,
        automatic_word=auto_word,
        doubling_word=double_word,
        magnitude_word=mag_word,
        magnitude_height=mag_height,
        magnitude_delta=mag_delta,
        fusion_signature=fusion,
    )


def key_bank(features: RowFeatures) -> dict[str, object]:
    return {
        "raw_automatic_word": features.automatic_word,
        "automatic_plus_chart_den": (features.automatic_word, features.first_chart_den),
        "automatic_plus_barcode_shape": (
            features.automatic_word,
            features.status,
            features.bar_count,
            fmt(features.longest_bar),
        ),
        "automatic_plus_magnitude": (
            features.automatic_word,
            features.magnitude_word,
            features.magnitude_height,
        ),
        "endpoint_current_plus_chart": (
            features.endpoint_current_word,
            features.first_chart,
            features.status,
        ),
        "fusion_signature": features.fusion_signature,
    }


def group_mixing(rows: list[RowFeatures]) -> None:
    print("Carrier quotient stress test")
    print("key rows fibers mixed_status_fibers mixed_route_fibers largest_fiber")
    for key_name in key_bank(rows[0]):
        groups: dict[object, list[RowFeatures]] = defaultdict(list)
        for features in rows:
            groups[key_bank(features)[key_name]].append(features)
        mixed_status = 0
        mixed_route = 0
        largest = max(len(group) for group in groups.values())
        examples: list[str] = []
        for key, group in groups.items():
            statuses = {item.status for item in group}
            routes = {item.row.route for item in group}
            if len(statuses) > 1:
                mixed_status += 1
                labels = ",".join(item.row.name for item in group[:5])
                examples.append(f"{repr(key)}=>{labels}")
            if len(routes) > 1:
                mixed_route += 1
        print(key_name, len(rows), len(groups), mixed_status, mixed_route, largest, sep=" | ")
        if examples:
            print(f"  mixed_status_examples: {' ; '.join(examples[:3])}")
    print()


def row_readout(rows: list[RowFeatures]) -> None:
    print("Exact row readout")
    print(
        "name route status M@t safe_mu bars longest_stick midpoint_slack "
        "chart endpoint_current magnitude automatic_word doubling_word"
    )
    for item in rows:
        print(
            item.row.name,
            item.row.route,
            item.status,
            f"{fmt(item.exact_m)}@{fmt(item.exact_t)}",
            fmt(item.safe_mu),
            item.bar_count,
            fmt(item.longest_bar),
            fmt(item.midpoint_slack),
            item.first_chart,
            item.endpoint_current_word,
            f"{item.magnitude_word}|h={item.magnitude_height}|d={item.magnitude_delta}",
            item.automatic_word,
            item.doubling_word,
            sep=" | ",
        )
    print()


@dataclass(frozen=True)
class Carrier:
    name: str
    predicate: int
    exactness: int
    topology: int
    arithmetic: int
    endpoint: int
    route_split: int
    family_transfer: int
    computable: int
    anti_scalar: int
    destroys: str

    def vector(self) -> dict[str, int]:
        return {
            "predicate": self.predicate,
            "exactness": self.exactness,
            "topology": self.topology,
            "arithmetic": self.arithmetic,
            "endpoint": self.endpoint,
            "route_split": self.route_split,
            "family_transfer": self.family_transfer,
            "computable": self.computable,
            "anti_scalar": self.anti_scalar,
        }


FIELDS = (
    "predicate",
    "exactness",
    "topology",
    "arithmetic",
    "endpoint",
    "route_split",
    "family_transfer",
    "computable",
    "anti_scalar",
)

TIE_PATH = (
    "labelled_packet_fusion_signature",
    "large_stick_potato_safe_body",
    "lonely_profile_barcode",
    "magnitude_cocycle_fiber",
    "endpoint_current_cocycle",
    "crt_solenoid_lift_chart",
    "danger_count_dual_distribution",
    "hurwitz_doubling_state",
    "automatic_word_sidecar",
    "raw_row_name",
)


def carrier_bank() -> list[Carrier]:
    return [
        Carrier("labelled_packet_fusion_signature", 4, 4, 4, 4, 4, 4, 4, 4, 4, "nothing if all sidecars are retained"),
        Carrier("large_stick_potato_safe_body", 4, 4, 4, 1, 2, 3, 3, 4, 4, "arithmetic period and automaton fiber identity"),
        Carrier("lonely_profile_barcode", 4, 4, 4, 1, 2, 3, 3, 4, 4, "route labels unless attached to packets"),
        Carrier("magnitude_cocycle_fiber", 3, 4, 2, 4, 2, 4, 4, 4, 4, "safe topology inside the fiber"),
        Carrier("endpoint_current_cocycle", 3, 3, 3, 2, 4, 3, 3, 4, 4, "off-boundary open bars"),
        Carrier("crt_solenoid_lift_chart", 3, 4, 2, 4, 2, 3, 3, 4, 3, "real interval topology between charts"),
        Carrier("danger_count_dual_distribution", 2, 4, 2, 1, 1, 2, 2, 4, 3, "endpoint owners and packet route"),
        Carrier("hurwitz_doubling_state", 2, 2, 1, 4, 1, 2, 3, 4, 3, "archimedean safe interval geometry"),
        Carrier("automatic_word_sidecar", 1, 2, 1, 4, 1, 1, 2, 4, 3, "magnitude and endpoint ownership"),
        Carrier("raw_row_name", 0, 0, 0, 0, 0, 0, 0, 1, 0, "all theorem-facing structure"),
    ]


def compare(a: Carrier, b: Carrier) -> str:
    av = a.vector()
    bv = b.vector()
    aw = sum(1 for field in FIELDS if av[field] > bv[field])
    bw = sum(1 for field in FIELDS if bv[field] > av[field])
    if aw > bw:
        return a.name
    if bw > aw:
        return b.name
    order = {name: idx for idx, name in enumerate(TIE_PATH)}
    return a.name if order[a.name] < order[b.name] else b.name


def edge_map(carriers: list[Carrier]) -> dict[tuple[str, str], str]:
    edges: dict[tuple[str, str], str] = {}
    for a, b in combinations(carriers, 2):
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
    carriers = carrier_bank()
    names = [carrier.name for carrier in carriers]
    edges = edge_map(carriers)
    hp_count, hp_first = hamiltonian_paths(names, edges)
    print("Carrier-Fusion Tournament Analysis")
    print("vertices=proof carriers and sidecar bundles, not runners/arcs")
    print(f"pairwise_observable={','.join(FIELDS)}")
    print(f"tie_hamiltonian_path={'>'.join(TIE_PATH)}")
    print(f"score_hist={score_histogram(names, edges)}")
    print(f"directed_3cycles={directed_3cycles(names, edges)}")
    print(f"scc_sizes={scc_sizes(names, edges)}")
    print(f"hamiltonian_path_count={hp_count}")
    print(f"first_hamiltonian_path={'>'.join(hp_first) if hp_first else 'NONE'}")
    print()


def switchboard_readout(rows: list[RowFeatures]) -> None:
    print("Switchboard readout")
    print("best_creative_carriers:")
    print("  large_stick_potato_safe_body = largest strict safe interval plus total safe body; geometric analogue of largest inscribed safe subset")
    print("  labelled_packet_fusion_signature = status, barcode, CRT chart, endpoint current, magnitude cocycle, automatic word, danger-count dual")
    print("  magnitude_cocycle_fiber = missing coordinate inside the AP/GW automatic-language fibers")
    print("  hurwitz_doubling_state = useful 2-adic sidecar, but not an archimedean certificate")
    print("leakage_rule:")
    print("  automatic words and doubling states remain telemetry until fusion_signature or an equivalent packet theorem proves fiber purity")
    open_rows = [item for item in rows if item.status == "open"]
    smallest = min(open_rows, key=lambda item: (item.safe_mu, item.row.name))
    print(f"smallest_positive_open_row={smallest.row.name}|mu={fmt(smallest.safe_mu)}|stick={fmt(smallest.longest_bar)}|slack={fmt(smallest.midpoint_slack)}")
    cover = next(item for item in rows if item.row.name == "cover_12_to_84")
    print(f"covering_chart_shift={cover.row.name}|chart={cover.first_chart}|endpoint_current={cover.endpoint_current_word}|magnitude={cover.magnitude_word}")
    print()


def assumption_challenge() -> None:
    print("Assumption challenge")
    print(
        "Alternate vertices considered: runners, gaps, fixed circle sections, "
        "section boundaries, wall-crossing events, residues, cover arcs, Fourier "
        "modes, matroid circuits, persistence bars, CRT charts, doubling states, "
        "large safe sticks/potatoes, quotient fibers, and proof obligations."
    )
    print(
        "Chosen vertices are proof carriers.  The quotient preserves the LRC14 "
        "predicate only after strict-open status, exact safe topology, chart, "
        "endpoint current, magnitude cocycle, and sidecar words stay attached.  "
        "It destroys raw runner identity and rejects the assumption that Moser/"
        "fibbinary or 2-adic state is theorem-safe by itself."
    )


def main() -> None:
    print("LRC14 carrier-fusion switchboard (codex-2026-06-26-S189)")
    print("source_threads=HYP-3025,HYP-3024,HYP-3023,HYP-3022,HYP-3021,HYP-3020,HYP-3018,HYP-3017,HYP-3016,HYP-3015,HYP-3014,HYP-3013,HYP-3009,HYP-2963")
    print("external_prompts=2-adic_Littlewood_Hurwitz_doubling,Ostrowski-Hadamard_lacunarity,large_sticks_and_potatoes")
    print()
    rows = [row_features(row) for row in named_rows()]
    row_readout(rows)
    group_mixing(rows)
    tournament_readout()
    switchboard_readout(rows)
    assumption_challenge()


if __name__ == "__main__":
    main()
