#!/usr/bin/env python3
"""S674: signed LRC clocks and unit-distance threshold tournaments.

User reframes:

1. Let some runners move in the opposite direction.  As an observer predicate,
   LRC is invariant under independent sign flips because dist(x)=dist(-x) on
   the circle.  But the pair clocks change from differences to sums when two
   runners have opposite signs.  The useful object is therefore a sign lift of
   the same observer orbit, not a new LRC instance.

2. For unit-distance point sets, treat distances <1, =1, and >1 as a
   three-state tournament switch rather than a binary unit/nonunit flip.  Unit
   edges become a trienerment/tie layer whose sign chooses one of two directed
   resolutions.

This is a finite scout, not a proof.  It records exact finite invariances and
the side channels they suggest for HYP-2249.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, product
from math import gcd


N = 14
C = 2 * N - 1
GRID_Q = 756  # lcm(14, 27, 28), enough to see AP, Vstar, and 2AP wall ticks.
DELTA = Fraction(1, N)

AP = tuple(range(1, N))
VSTAR = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24)
TWOP = tuple(2 * v for v in AP)


def mod_dist_num(value: int, modulus: int) -> int:
    r = value % modulus
    return min(r, modulus - r)


def grid_safe_mask(row: tuple[int, ...], signs: tuple[int, ...]) -> tuple[int, ...]:
    """Return q-grid ticks where all signed runners are at distance >= 1/N."""

    threshold_num = GRID_Q // N
    mask = []
    for tick in range(1, GRID_Q):
        ok = True
        for sign, speed in zip(signs, row):
            if mod_dist_num(sign * speed * tick, GRID_Q) < threshold_num:
                ok = False
                break
        if ok:
            mask.append(tick)
    return tuple(mask)


def shell(residue: int, modulus: int = C) -> int:
    residue %= modulus
    if residue == 0:
        return 0
    return min(residue, modulus - residue)


def sign_patterns(m: int) -> dict[str, tuple[int, ...]]:
    return {
        "all_plus": tuple([1] * m),
        "alternating": tuple(1 if i % 2 == 0 else -1 for i in range(m)),
        "split": tuple(1 if i < (m + 1) // 2 else -1 for i in range(m)),
        "mod3_plus": tuple(1 if (i + 1) % 3 else -1 for i in range(m)),
    }


@dataclass(frozen=True)
class PairClockStats:
    same_pairs: int
    opposite_pairs: int
    zero_clocks: int
    gcd_hist: tuple[tuple[int, int], ...]
    shell_hist: tuple[tuple[int, int], ...]
    exposed_sum_pairs: tuple[tuple[int, int, int], ...]


def pair_clock_stats(row: tuple[int, ...], signs: tuple[int, ...]) -> PairClockStats:
    same_pairs = 0
    opposite_pairs = 0
    zero_clocks = 0
    gcd_hist: Counter[int] = Counter()
    shell_hist: Counter[int] = Counter()
    exposed_sum_pairs = []
    residues = tuple(v % C for v in row)
    for i, j in combinations(range(len(row)), 2):
        a, b = residues[i], residues[j]
        if signs[i] == signs[j]:
            same_pairs += 1
        else:
            opposite_pairs += 1
        clock = (signs[i] * a - signs[j] * b) % C
        sh = shell(clock)
        if sh == 0:
            zero_clocks += 1
            gcd_hist[C] += 1
        else:
            gcd_hist[gcd(sh, C)] += 1
        shell_hist[sh] += 1
        if signs[i] != signs[j]:
            exposed_sum_pairs.append((row[i], row[j], sh))
    return PairClockStats(
        same_pairs=same_pairs,
        opposite_pairs=opposite_pairs,
        zero_clocks=zero_clocks,
        gcd_hist=tuple(sorted(gcd_hist.items())),
        shell_hist=tuple(sorted(shell_hist.items())),
        exposed_sum_pairs=tuple(exposed_sum_pairs[:8]),
    )


def best_sign_lifts(row: tuple[int, ...], limit: int = 5) -> list[tuple[tuple[int, int, int, int], str, PairClockStats]]:
    """Search sign gauges, fixing the first sign to +, for pair-clock pressure."""

    rows = []
    for tail in product((-1, 1), repeat=len(row) - 1):
        signs = (1,) + tail
        stats = pair_clock_stats(row, signs)
        hist = dict(stats.gcd_hist)
        objective = (
            stats.zero_clocks,
            hist.get(9, 0),
            hist.get(3, 0),
            stats.opposite_pairs,
        )
        word = "".join("+" if s > 0 else "-" for s in signs)
        rows.append((objective, word, stats))
    return sorted(rows, reverse=True)[:limit]


def lrc_signed_report() -> list[str]:
    lines = []
    lines.append("== Signed LRC gauge audit ==")
    lines.append(
        "Fact used: circle distance is even, so each independent sign flip "
        "preserves the observer loneliness predicate for every real time."
    )
    for name, row in (("AP", AP), ("Vstar", VSTAR), ("2AP", TWOP)):
        patterns = sign_patterns(len(row))
        base_mask = grid_safe_mask(row, patterns["all_plus"])
        lines.append("")
        lines.append(f"{name}: speeds={row}")
        lines.append(f"  grid q={GRID_Q} observer-safe ticks for all_plus: {len(base_mask)}")
        for pattern_name, signs in patterns.items():
            mask = grid_safe_mask(row, signs)
            same_observer = mask == base_mask
            stats = pair_clock_stats(row, signs)
            lines.append(
                "  "
                + f"{pattern_name:12s} observer_same={same_observer} "
                + f"same/opposite={stats.same_pairs}/{stats.opposite_pairs} "
                + f"zero_pair_clocks={stats.zero_clocks} gcd_hist={stats.gcd_hist}"
            )
        lines.append("  best sign lifts for zero/gcd9/gcd3/opposite pressure:")
        for objective, word, stats in best_sign_lifts(row):
            lines.append(
                "    "
                + f"obj={objective} signs={word} "
                + f"gcd_hist={stats.gcd_hist} first_sum_pairs={stats.exposed_sum_pairs}"
            )
    return lines


# Exact Eisenstein/triangular-lattice support.
AxialPoint = tuple[int, int]
RationalPoint = tuple[Fraction, Fraction]


def axial_norm2(delta: AxialPoint) -> int:
    a, b = delta
    return a * a + a * b + b * b


def axial_distance_class(p: AxialPoint, q: AxialPoint) -> str:
    n2 = axial_norm2((p[0] - q[0], p[1] - q[1]))
    if n2 < 1:
        return "close"
    if n2 == 1:
        return "unit"
    return "far"


def rational_distance_class(p: RationalPoint, q: RationalPoint) -> str:
    dx = p[0] - q[0]
    dy = p[1] - q[1]
    d2 = dx * dx + dy * dy
    if d2 < 1:
        return "close"
    if d2 == 1:
        return "unit"
    return "far"


def triangular_cluster(n: int) -> tuple[AxialPoint, ...]:
    pts = []
    radius = 4
    for a in range(-radius, radius + 1):
        for b in range(-radius, radius + 1):
            pts.append((a, b))
    pts.sort(key=lambda p: (axial_norm2(p), abs(p[0]) + abs(p[1]), p[0], p[1]))
    return tuple(pts[:n])


def undirected_mask(points: tuple, classifier) -> list[int]:
    n = len(points)
    adj = [0] * n
    for i, j in combinations(range(n), 2):
        if classifier(points[i], points[j]) == "unit":
            adj[i] |= 1 << j
            adj[j] |= 1 << i
    return adj


def find_hamiltonian_path(adj: list[int]) -> tuple[int, ...] | None:
    n = len(adj)
    parent: dict[tuple[int, int], tuple[int, int] | None] = {}
    states = []
    for v in range(n):
        state = (1 << v, v)
        parent[state] = None
        states.append(state)
    for mask in range(1 << n):
        for last in range(n):
            state = (mask, last)
            if state not in parent:
                continue
            avail = adj[last] & ~mask
            while avail:
                bit = avail & -avail
                avail -= bit
                nxt = bit.bit_length() - 1
                new_state = (mask | bit, nxt)
                if new_state not in parent:
                    parent[new_state] = state
    full = (1 << n) - 1
    for last in range(n):
        state = (full, last)
        if state not in parent:
            continue
        path = []
        while state is not None:
            path.append(state[1])
            state = parent[state]
        return tuple(reversed(path))
    return None


def reorder_by_path(points: tuple, path: tuple[int, ...] | None) -> tuple:
    if path is None:
        return points
    return tuple(points[i] for i in path)


def threshold_tournament(points: tuple, classifier, unit_sign: int) -> list[int]:
    """Orient i<j by close/unit/far threshold states.

    close: i -> j
    far:   j -> i
    unit:  i -> j if unit_sign=+1, else j -> i
    """

    n = len(points)
    adj = [0] * n
    for i, j in combinations(range(n), 2):
        cls = classifier(points[i], points[j])
        forward = cls == "close" or (cls == "unit" and unit_sign > 0)
        if forward:
            adj[i] |= 1 << j
        else:
            adj[j] |= 1 << i
    return adj


def category_counts(points: tuple, classifier) -> tuple[tuple[str, int], ...]:
    counts: Counter[str] = Counter()
    for i, j in combinations(range(len(points)), 2):
        counts[classifier(points[i], points[j])] += 1
    return tuple(sorted(counts.items()))


def score_hist(adj: list[int]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(a.bit_count() for a in adj).items()))


def directed_3cycles(adj: list[int]) -> int:
    n = len(adj)
    total = 0
    for a, b, c in combinations(range(n), 3):
        ab = (adj[a] >> b) & 1
        bc = (adj[b] >> c) & 1
        ca = (adj[c] >> a) & 1
        ac = (adj[a] >> c) & 1
        cb = (adj[c] >> b) & 1
        ba = (adj[b] >> a) & 1
        if (ab and bc and ca) or (ac and cb and ba):
            total += 1
    return total


def hamiltonian_path_count(adj: list[int]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if not count:
                continue
            avail = adj[last] & ~mask
            while avail:
                bit = avail & -avail
                avail -= bit
                nxt = bit.bit_length() - 1
                dp[mask | bit][nxt] += count
    return sum(dp[-1])


def directed_all_unit_hp(adj: list[int], unit_adj: list[int]) -> bool:
    n = len(adj)
    dp = [0] * (1 << n)
    for v in range(n):
        dp[1 << v] |= 1 << v
    for mask in range(1 << n):
        ends = dp[mask]
        while ends:
            bit = ends & -ends
            ends -= bit
            last = bit.bit_length() - 1
            avail = adj[last] & unit_adj[last] & ~mask
            while avail:
                nb = avail & -avail
                avail -= nb
                nxt = nb.bit_length() - 1
                dp[mask | nb] |= nb
    return dp[-1] != 0


@dataclass(frozen=True)
class ThresholdRow:
    name: str
    n: int
    categories: tuple[tuple[str, int], ...]
    unit_graph_traceable: bool
    unit_sign: str
    score_hist: tuple[tuple[int, int], ...]
    c3: int
    hp_count: int
    directed_all_unit_hp: bool


def threshold_rows() -> list[ThresholdRow]:
    rows: list[ThresholdRow] = []
    for n in range(3, 11):
        raw = triangular_cluster(n)
        unit_adj_raw = undirected_mask(raw, axial_distance_class)
        path = find_hamiltonian_path(unit_adj_raw)
        pts = reorder_by_path(raw, path)
        unit_adj = undirected_mask(pts, axial_distance_class)
        for unit_sign, sign_name in ((1, "unit_positive"), (-1, "unit_negative")):
            adj = threshold_tournament(pts, axial_distance_class, unit_sign)
            rows.append(
                ThresholdRow(
                    name="triangular_spine",
                    n=n,
                    categories=category_counts(pts, axial_distance_class),
                    unit_graph_traceable=path is not None,
                    unit_sign=sign_name,
                    score_hist=score_hist(adj),
                    c3=directed_3cycles(adj),
                    hp_count=hamiltonian_path_count(adj),
                    directed_all_unit_hp=directed_all_unit_hp(adj, unit_adj),
                )
            )

    toys: tuple[tuple[str, tuple[RationalPoint, ...]], ...] = (
        (
            "compressed_line",
            (
                (Fraction(0), Fraction(0)),
                (Fraction(1, 2), Fraction(0)),
                (Fraction(1), Fraction(0)),
                (Fraction(2), Fraction(0)),
                (Fraction(3), Fraction(0)),
            ),
        ),
        (
            "square_with_center",
            (
                (Fraction(0), Fraction(0)),
                (Fraction(1), Fraction(0)),
                (Fraction(1), Fraction(1)),
                (Fraction(0), Fraction(1)),
                (Fraction(1, 2), Fraction(1, 2)),
            ),
        ),
    )
    for name, pts in toys:
        unit_adj = undirected_mask(pts, rational_distance_class)
        path = find_hamiltonian_path(unit_adj)
        for unit_sign, sign_name in ((1, "unit_positive"), (-1, "unit_negative")):
            adj = threshold_tournament(pts, rational_distance_class, unit_sign)
            rows.append(
                ThresholdRow(
                    name=name,
                    n=len(pts),
                    categories=category_counts(pts, rational_distance_class),
                    unit_graph_traceable=path is not None,
                    unit_sign=sign_name,
                    score_hist=score_hist(adj),
                    c3=directed_3cycles(adj),
                    hp_count=hamiltonian_path_count(adj),
                    directed_all_unit_hp=directed_all_unit_hp(adj, unit_adj),
                )
            )
    return rows


def threshold_report() -> list[str]:
    lines = []
    lines.append("")
    lines.append("== Unit-distance threshold tournament audit ==")
    lines.append(
        "Rule: for i<j, d<1 gives i->j, d>1 gives j->i, and d=1 is the "
        "trienerment layer resolved by a unit sign."
    )
    for row in threshold_rows():
        lines.append(
            f"  {row.name:18s} n={row.n:2d} {row.unit_sign:13s} "
            f"cats={row.categories} unit_traceable={row.unit_graph_traceable} "
            f"score={row.score_hist} c3={row.c3} H={row.hp_count} "
            f"all_unit_directed_hp={row.directed_all_unit_hp}"
        )
    return lines


def tarjan_scc(adj: list[int]) -> list[list[int]]:
    n = len(adj)
    index = 0
    stack: list[int] = []
    on_stack = [False] * n
    indices = [-1] * n
    low = [0] * n
    comps: list[list[int]] = []

    def strongconnect(v: int) -> None:
        nonlocal index
        indices[v] = index
        low[v] = index
        index += 1
        stack.append(v)
        on_stack[v] = True
        out = adj[v]
        for w in range(n):
            if not ((out >> w) & 1):
                continue
            if indices[w] == -1:
                strongconnect(w)
                low[v] = min(low[v], low[w])
            elif on_stack[w]:
                low[v] = min(low[v], indices[w])
        if low[v] == indices[v]:
            comp = []
            while True:
                w = stack.pop()
                on_stack[w] = False
                comp.append(w)
                if w == v:
                    break
            comps.append(comp)

    for v in range(n):
        if indices[v] == -1:
            strongconnect(v)
    return comps


def route_tournament() -> tuple[list[str], list[int]]:
    names = [
        "signed_observer_gauge",
        "opposite_pair_sum_clock",
        "alternating_sign_search",
        "owner_carry_rank",
        "unit_threshold_tristate",
        "unit_spine_tie_order",
        "raw_negative_speed",
        "raw_unit_flip",
    ]
    metrics = {
        "signed_observer_gauge": (5, 1, 2, 0, 3, 3),
        "opposite_pair_sum_clock": (3, 5, 4, 0, 4, 5),
        "alternating_sign_search": (2, 4, 3, 0, 3, 4),
        "owner_carry_rank": (3, 3, 5, 0, 5, 2),
        "unit_threshold_tristate": (1, 2, 1, 5, 3, 4),
        "unit_spine_tie_order": (1, 1, 2, 5, 4, 3),
        "raw_negative_speed": (2, 0, 0, 0, 1, 2),
        "raw_unit_flip": (0, 1, 0, 3, 1, 1),
    }
    adj = [0] * len(names)
    for i, j in combinations(range(len(names)), 2):
        a = metrics[names[i]]
        b = metrics[names[j]]
        votes = sum(1 if x > y else -1 if x < y else 0 for x, y in zip(a, b))
        if votes >= 0:
            adj[i] |= 1 << j
        else:
            adj[j] |= 1 << i
    return names, adj


def tournament_analysis_report() -> list[str]:
    names, adj = route_tournament()
    scores = score_hist(adj)
    comps = sorted((len(c) for c in tarjan_scc(adj)), reverse=True)
    order = sorted(range(len(names)), key=lambda i: adj[i].bit_count(), reverse=True)
    lines = []
    lines.append("")
    lines.append("== Route Tournament Analysis ==")
    lines.append("Pairwise observable: retained side-channel power across six proof-use metrics.")
    lines.append("Switch/gauge: majority comparison of metric vectors.")
    lines.append("Tie Hamiltonian path: score-descending route order.")
    lines.append(
        f"score_hist={scores} directed_3cycles={directed_3cycles(adj)} "
        f"scc_sizes={comps} H={hamiltonian_path_count(adj)}"
    )
    lines.append("route_order=" + " > ".join(names[i] for i in order))
    return lines


def synthesis_report() -> list[str]:
    return [
        "",
        "== Synthesis ==",
        "1. Negative speed is not a new observer predicate; it is a sign gauge.",
        "2. Opposite signs expose pair sums v_i+v_j in the same orbit where same signs expose differences.",
        "3. Vstar's exceptional nonunit seam is visible as a sign pattern with a zero pair clock on 3+24=27.",
        "4. Alternating signs are therefore a cheap way to turn the C=27 pair-sum ledger into pair-clock data.",
        "5. Unit distance wants a three-state threshold tournament, not just a unit/nonunit flip.",
        "6. In normalized lattice/unit-spine rows, there are no d<1 pairs, so the three-state gauge collapses back to the old binary spine gauge.",
        "7. The next proof target is a signed owner/carry rank: every sign-lifted sum-clock obstruction must either be AP/Vstar/2AP scalar-floor data or pay a strict owner-deletion tax.",
        "Challenged assumption: tournament vertices need not be runners or points; useful vertices here are sign gauges, pair clocks, threshold pair states, spine steps, and proof obligations.",
    ]


def main() -> None:
    for line in (
        lrc_signed_report()
        + threshold_report()
        + tournament_analysis_report()
        + synthesis_report()
    ):
        print(line)


if __name__ == "__main__":
    main()
