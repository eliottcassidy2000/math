#!/usr/bin/env python3
"""
lrc_n14_n18_tournament_pingpong_s481.py

codex-2026-06-01 S481

Alternate between the n=14 and n=18 Lonely Runner frontiers, using
Tournament Analysis whenever a pairwise object is available.

The working protocol is:

1. start from an n=14 obstruction;
2. translate the same question to n=18;
3. when the scalar route stalls, inject a deterministic "noise card";
4. turn the resulting idea into a small exact hypothesis check.

The two main comparisons are:

* local n-gate cover families;
* row-parent / gate / double-gate ladder debt cascades.

The Tournament Analysis layer records three gauges at selected times:

* semicircle orientation of the runner positions;
* close-pair threshold switch over the base Hamiltonian path;
* two-neighbor deletion-pressure tournament.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import gcd
from pathlib import Path
from random import Random


ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()
S360 = SourceFileLoader(
    "lonely_runner_endpoint_protection_s360",
    str(ROOT / "04-computation" / "lonely_runner_endpoint_protection_s360.py"),
).load_module()
S411 = SourceFileLoader(
    "lrc_column_row_modes_s411",
    str(ROOT / "04-computation" / "lrc_column_row_modes_s411.py"),
).load_module()
S420 = SourceFileLoader(
    "lrc_integer_programming_modes_s420",
    str(ROOT / "04-computation" / "lrc_integer_programming_modes_s420.py"),
).load_module()


ONE = Fraction(1, 1)
HALF = Fraction(1, 2)


@dataclass(frozen=True)
class SpeedRow:
    label: str
    n: int
    speeds: tuple[int, ...]
    mode: str
    scale: int | None = None
    skip: int | None = None


@dataclass(frozen=True)
class EndpointLedger:
    label: str
    n: int
    mode: str
    classification: str
    gap_ratio: Fraction
    unprotected: int
    gap_debt_product: Fraction
    first_unprotected: Fraction | None
    depth_hist: tuple[tuple[tuple[int, ...], int], ...]
    frontier_mass: Fraction
    denominator_pressure: int


@dataclass(frozen=True)
class CoverFamily:
    n: int
    exact_size: int | None
    forced: tuple[int, ...]
    cover_count: int
    free_parts: tuple[tuple[int, ...], ...]
    example: tuple[int, ...]
    search_calls: int


@dataclass(frozen=True)
class TournamentLedger:
    row_label: str
    n: int
    time_kind: str
    t: Fraction
    gauge: str
    info: str
    score_hist: str
    cyclic_triples: int
    scc_sizes: tuple[int, ...]
    hamiltonian_paths: int | None


def fmt(value: Fraction | None) -> str:
    return S356.fmt_frac(value)


def fmt_ratio(value: Fraction) -> str:
    return fmt(value)


def fmt_float(value: Fraction) -> str:
    return f"{float(value):.6f}"


def hist_text(counter: Counter[int]) -> str:
    if not counter:
        return "-"
    return " ".join(
        f"{key}:{counter[key]}" for key in sorted(counter)
    )


def circle(value: Fraction) -> Fraction:
    return value % ONE


def clockwise_delta(a: Fraction, b: Fraction) -> Fraction:
    return circle(b - a)


def circular_distance(a: Fraction, b: Fraction) -> Fraction:
    delta = clockwise_delta(a, b)
    return min(delta, ONE - delta)


def prime_factors(value: int) -> tuple[int, ...]:
    out: list[int] = []
    p = 2
    while p * p <= value:
        if value % p == 0:
            out.append(p)
            while value % p == 0:
                value //= p
        p += 1 if p == 2 else 2
    if value > 1:
        out.append(value)
    return tuple(out)


def vp(value: int, prime: int) -> int:
    value = abs(value)
    out = 0
    while value and value % prime == 0:
        out += 1
        value //= prime
    return out


def extra_depth(n: int, point: Fraction) -> tuple[int, ...]:
    return tuple(
        max(0, vp(point.denominator, prime) - vp(n, prime))
        for prime in prime_factors(n)
    )


def depth_scale(n: int, depth: tuple[int, ...]) -> int:
    scale = 1
    for prime, height in zip(prime_factors(n), depth):
        scale *= prime**height
    return scale


def fmt_depth(n: int, depth: tuple[int, ...]) -> str:
    primes = prime_factors(n)
    return "{" + ",".join(
        f"{prime}:+{height}" for prime, height in zip(primes, depth)
    ) + "}"


def fmt_depth_hist(n: int, hist: tuple[tuple[tuple[int, ...], int], ...]) -> str:
    if not hist:
        return "-"
    return " ".join(f"{fmt_depth(n, depth)}:{count}" for depth, count in hist)


def initial(n: int) -> tuple[int, ...]:
    return tuple(range(1, n))


def ladder(n: int, scale: int, skip: int) -> tuple[int, ...]:
    speeds = tuple(sorted({1} | {scale * q for q in range(1, n) if q != skip}))
    if len(speeds) != n - 1:
        raise ValueError(f"bad ladder length n={n} scale={scale} skip={skip}")
    if gcd_all(speeds) != 1:
        raise ValueError(f"nonprimitive ladder n={n} scale={scale} skip={skip}")
    return speeds


def gcd_all(values: tuple[int, ...]) -> int:
    g = 0
    for value in values:
        g = gcd(g, value)
    return g


@lru_cache(maxsize=None)
def summary(speeds: tuple[int, ...]):
    return S360.summarize(list(speeds))


@lru_cache(maxsize=None)
def gap_report(speeds: tuple[int, ...]):
    return S356.report("cached", list(speeds))


def unprotected_values(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    values = sorted({endpoint.value for endpoint in S360.endpoints(speeds)})
    return tuple(
        value
        for value in values
        if not any(S360.direct_protects(speeds, speed, value) for speed in speeds)
    )


def endpoint_ledger(row: SpeedRow) -> EndpointLedger:
    item = summary(row.speeds)
    gap_ratio = item.max_gap / item.threshold
    bad = unprotected_values(row.speeds)
    depth_counter = Counter(extra_depth(row.n, point) for point in bad)
    frontier_mass = sum(
        Fraction(1, depth_scale(row.n, extra_depth(row.n, point))) for point in bad
    )
    denominator_pressure = sum(
        depth_scale(row.n, extra_depth(row.n, point)) for point in bad
    )
    return EndpointLedger(
        label=row.label,
        n=row.n,
        mode=row.mode,
        classification=item.classification,
        gap_ratio=gap_ratio,
        unprotected=item.unprotected_count,
        gap_debt_product=gap_ratio * item.unprotected_count,
        first_unprotected=item.first_unprotected,
        depth_hist=tuple(sorted(depth_counter.items())),
        frontier_mass=frontier_mass,
        denominator_pressure=denominator_pressure,
    )


def speed_rows() -> tuple[SpeedRow, ...]:
    rows: list[SpeedRow] = []
    for n in (14, 18):
        rows.append(SpeedRow(f"n{n} initial", n, initial(n), "initial"))

    for level, factor in (("row-parent", Fraction(1, 2)), ("gate", Fraction(1, 1)), ("double-gate", Fraction(2, 1))):
        for n in (14, 18):
            best = S411.best_lpd_ladder(n)
            if best.skip is None:
                continue
            scale = int(n * factor)
            rows.append(
                SpeedRow(
                    f"n{n} {level}",
                    n,
                    ladder(n, scale, best.skip),
                    level,
                    scale=scale,
                    skip=best.skip,
                )
            )
    return tuple(rows)


def cover_family(n: int) -> CoverFamily:
    targets = S420.endpoint_values_for_owner(n, n)
    candidates = tuple(range(1, n))
    result = S420.set_cover_result(
        f"n{n} owner {n} endpoints, lower 1..{n-1}",
        n,
        targets,
        candidates,
    )
    columns = S420.build_cover_columns(n, targets, candidates)
    speed_to_mask = {column.speed: column.mask for column in columns}
    full_mask = (1 << len(targets)) - 1

    covers: list[tuple[int, ...]] = []
    if result.exact_size is not None:
        for combo in combinations(candidates, result.exact_size):
            mask = 0
            for speed in combo:
                mask |= speed_to_mask.get(speed, 0)
            if mask == full_mask:
                covers.append(combo)

    forced = set(result.forced_columns)
    free_parts = tuple(
        sorted({tuple(speed for speed in cover if speed not in forced) for cover in covers})
    )
    return CoverFamily(
        n=n,
        exact_size=result.exact_size,
        forced=result.forced_columns,
        cover_count=len(covers),
        free_parts=free_parts,
        example=result.exact_columns,
        search_calls=result.search_calls,
    )


def positions(runner_speeds: tuple[int, ...], t: Fraction) -> tuple[Fraction, ...]:
    return tuple(circle(speed * t) for speed in runner_speeds)


def score_hist(adj: list[list[bool]]) -> str:
    scores = [sum(row) for row in adj]
    counter = Counter(scores)
    return " ".join(
        f"{score}^{counter[score]}" if counter[score] > 1 else str(score)
        for score in sorted(counter)
    )


def directed_triangles(adj: list[list[bool]]) -> int:
    total = 0
    for a, b, c in combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            total += 1
    return total


def scc_sizes(adj: list[list[bool]]) -> tuple[int, ...]:
    n = len(adj)
    graph = [[] for _ in range(n)]
    reverse = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if adj[i][j]:
                graph[i].append(j)
                reverse[j].append(i)

    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for nxt in graph[v]:
            if nxt not in seen:
                dfs(nxt)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)

    sizes: list[int] = []
    seen.clear()
    for start in reversed(order):
        if start in seen:
            continue
        todo = deque([start])
        seen.add(start)
        size = 0
        while todo:
            v = todo.pop()
            size += 1
            for nxt in reverse[v]:
                if nxt not in seen:
                    seen.add(nxt)
                    todo.append(nxt)
        sizes.append(size)
    return tuple(sorted(sizes, reverse=True))


def tournament_key(adj: list[list[bool]]) -> tuple[int, tuple[int, ...]]:
    rows = []
    for row in adj:
        mask = 0
        for j, value in enumerate(row):
            if value:
                mask |= 1 << j
        rows.append(mask)
    return len(adj), tuple(rows)


@lru_cache(maxsize=None)
def hamiltonian_paths_from_key(n: int, row_masks: tuple[int, ...]) -> int:
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            total = dp[mask][last]
            if not total:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if row_masks[last] & (1 << nxt):
                    dp[mask | (1 << nxt)][nxt] += total
    return sum(dp[-1])


def hamiltonian_paths(adj: list[list[bool]]) -> int:
    n, row_masks = tournament_key(adj)
    return hamiltonian_paths_from_key(n, row_masks)


def summarize_tournament(
    row_label: str,
    n: int,
    time_kind: str,
    t: Fraction,
    gauge: str,
    adj: list[list[bool]],
    info: str,
) -> TournamentLedger:
    return TournamentLedger(
        row_label=row_label,
        n=n,
        time_kind=time_kind,
        t=t,
        gauge=gauge,
        info=info,
        score_hist=score_hist(adj),
        cyclic_triples=directed_triangles(adj),
        scc_sizes=scc_sizes(adj),
        hamiltonian_paths=hamiltonian_paths(adj) if len(adj) <= 14 else None,
    )


def semicircle_tournament(pos: tuple[Fraction, ...]) -> tuple[list[list[bool]], str]:
    n = len(pos)
    adj = [[False] * n for _ in range(n)]
    collisions = 0
    antipodal = 0
    for i, j in combinations(range(n), 2):
        delta = clockwise_delta(pos[i], pos[j])
        if delta == 0:
            collisions += 1
            winner = i
        elif delta == HALF:
            antipodal += 1
            winner = i
        elif delta < HALF:
            winner = i
        else:
            winner = j
        loser = j if winner == i else i
        adj[winner][loser] = True
    return adj, f"collision={collisions} antipodal={antipodal}"


def close_threshold_tournament(
    pos: tuple[Fraction, ...], threshold: Fraction
) -> tuple[list[list[bool]], str]:
    n = len(pos)
    adj = [[False] * n for _ in range(n)]
    close_flips = 0
    exact_threshold = 0
    for i, j in combinations(range(n), 2):
        d = circular_distance(pos[i], pos[j])
        if d < threshold:
            close_flips += 1
            winner = j
        else:
            if d == threshold:
                exact_threshold += 1
            winner = i
        loser = j if winner == i else i
        adj[winner][loser] = True
    return adj, f"close_flips={close_flips} exact_th={exact_threshold}"


def distance_matrix(pos: tuple[Fraction, ...]) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(
        tuple(Fraction(0) if i == j else circular_distance(pos[i], pos[j]) for j in range(len(pos)))
        for i in range(len(pos))
    )


def nearest_distances(dist: tuple[tuple[Fraction, ...], ...]) -> tuple[Fraction, ...]:
    return tuple(
        min(dist[i][j] for j in range(len(dist)) if i != j)
        for i in range(len(dist))
    )


def nearest_without(
    dist: tuple[tuple[Fraction, ...], ...], i: int, removed: int
) -> Fraction:
    return min(dist[i][j] for j in range(len(dist)) if j not in (i, removed))


def pressure_tournament(pos: tuple[Fraction, ...]) -> tuple[list[list[bool]], str]:
    n = len(pos)
    dist = distance_matrix(pos)
    d1 = nearest_distances(dist)
    adj = [[False] * n for _ in range(n)]
    relief_ties = 0
    for i, j in combinations(range(n), 2):
        relief_i_from_j = nearest_without(dist, i, j) - d1[i]
        relief_j_from_i = nearest_without(dist, j, i) - d1[j]
        if relief_i_from_j > relief_j_from_i:
            winner = j
        elif relief_j_from_i > relief_i_from_j:
            winner = i
        else:
            relief_ties += 1
            winner = i
        loser = j if winner == i else i
        adj[winner][loser] = True
    return adj, f"relief_ties={relief_ties}"


def mixed_two_neighbor_tournament(
    pos: tuple[Fraction, ...], threshold: Fraction
) -> tuple[list[list[bool]], str]:
    """Unequal-threshold noise gauge from first/second-nearest slack.

    This is intentionally a rank-like control gauge.  It asks whether a
    mixed-threshold/two-neighbor idea keeps cyclic information or collapses to
    an order.  The edge-local pressure gauge above is the non-collapsing
    alternative.
    """

    n = len(pos)
    dist = distance_matrix(pos)
    slacks: list[tuple[Fraction, Fraction]] = []
    for i in range(n):
        ordered = sorted(dist[i][j] for j in range(n) if i != j)
        d1, d2 = ordered[0], ordered[1]
        slacks.append((d1 - threshold, d2 - threshold))

    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        # Higher two-neighbor slack is safer; ties use lower label.
        if slacks[i] > slacks[j]:
            winner = i
        elif slacks[j] > slacks[i]:
            winner = j
        else:
            winner = i
        loser = j if winner == i else i
        adj[winner][loser] = True
    distinct_profiles = len(set(slacks))
    return adj, f"profiles={distinct_profiles}"


def selected_time(row: SpeedRow) -> tuple[str, Fraction]:
    gap = gap_report(row.speeds)
    item = summary(row.speeds)
    if gap.witness is not None:
        return "gap-mid", gap.witness
    if item.first_unprotected is not None:
        return "boundary", item.first_unprotected
    return "unit", Fraction(1, row.n)


def tournament_ledgers(rows: tuple[SpeedRow, ...]) -> tuple[TournamentLedger, ...]:
    out: list[TournamentLedger] = []
    for row in rows:
        time_kind, t = selected_time(row)
        pos = positions((0,) + row.speeds, t)
        threshold = Fraction(1, row.n)
        for gauge, builder in (
            ("semicircle", lambda p, th: semicircle_tournament(p)),
            ("close-threshold", close_threshold_tournament),
            ("pressure", lambda p, th: pressure_tournament(p)),
            ("mixed-two-neighbor", mixed_two_neighbor_tournament),
        ):
            adj, info = builder(pos, threshold)
            out.append(
                summarize_tournament(
                    row.label,
                    row.n,
                    time_kind,
                    t,
                    gauge,
                    adj,
                    info,
                )
            )
    return tuple(out)


def inspiration_cards() -> tuple[str, ...]:
    cards = (
        "product-divisibility sieve: bridge choices should be read as product rows, not decorations",
        "mixed thresholds: let different runners carry different local radii before collapsing to a tournament",
        "Zeckendorf carry: a legal repair is a non-adjacent support move, not a raw scalar improvement",
        "SC blowup: first-even rows create a twin sheet that hides inherited unit data",
        "basketball tie path: arbitrary labels are still a Hamiltonian path and should be declared",
        "Bruhat-Tits frontier: shrinking the real gap should push endpoint mass deeper, not erase it",
        "two-neighbor braid: the nearest runner alone is a projection of a left/right handoff pair",
    )
    rng = Random(481)
    return tuple(rng.sample(cards, 4))


def print_header(title: str) -> None:
    print("=" * 104)
    print(title)
    print("=" * 104)


def print_noise_cards() -> None:
    print_header("DETERMINISTIC FORCED-RANDOMNESS CARDS")
    print("Seed 481.  These are not evidence; they decide what to test next.")
    for index, card in enumerate(inspiration_cards(), start=1):
        print(f"{index}. {card}")
    print()


def print_cover_comparison() -> None:
    print_header("PING: LOCAL n-GATE COVER FAMILIES")
    print("Cover endpoints owned by the n-gate using only lower columns 1..n-1.")
    print("n   forced columns                 exact  covers  free bridge parts       example")
    print("-" * 104)
    for n in (14, 18):
        family = cover_family(n)
        free = ", ".join(str(part) for part in family.free_parts)
        print(
            f"{n:<3} {str(family.forced):<30} "
            f"{str(family.exact_size):>5} {family.cover_count:>7} "
            f"{free:<23} {family.example}"
        )
    print()
    print("Reading:")
    print("- n=14: forced odd fan plus one arbitrary even bridge.")
    print("- n=18: forced unit/half fan plus one square-torsion bridge, only 6 or 12.")
    print("- Both local invoices have size 8, but n=18 has a much thinner bridge fiber.")
    print()


def print_endpoint_ledgers(rows: tuple[SpeedRow, ...]) -> None:
    print_header("PONG: ROW-PARENT/GATE DEBT CASCADE")
    print("Rows are interleaved n=14 then n=18 at each repair depth.")
    print(
        "label              mode          class          gap/th      debt  product  "
        "frontier pressure depth histogram"
    )
    print("-" * 104)
    for row in rows:
        ledger = endpoint_ledger(row)
        print(
            f"{ledger.label:<18} {ledger.mode:<13} {ledger.classification:<14} "
            f"{fmt_ratio(ledger.gap_ratio):>9} {ledger.unprotected:>7} "
            f"{fmt_ratio(ledger.gap_debt_product):>8} "
            f"{fmt_ratio(ledger.frontier_mass):>8} {ledger.denominator_pressure:>8} "
            f"{fmt_depth_hist(ledger.n, ledger.depth_hist)}"
        )
    print()
    print("Check:")
    print("- n=14 ladder levels preserve gap*debt = 5/11 while depth moves down the 2-sheet.")
    print("- n=18 ladder levels preserve gap*debt = 1 with the same doubling pattern.")
    print("- The scalar gap improves only by exporting twice as many endpoint witnesses.")
    print()


def print_tournament_table(rows: tuple[SpeedRow, ...]) -> None:
    print_header("TOURNAMENT ANALYSIS AT SELECTED TIMES")
    print(
        "For positive-gap rows use the middle of the largest gap; for boundary rows "
        "use the first exposed endpoint. H is omitted for n=18 to keep the audit fast."
    )
    print(
        "row                t-kind   t             gauge              info"
        "                         c3   scc       H          score"
    )
    print("-" * 126)
    for item in tournament_ledgers(rows):
        h_text = "-" if item.hamiltonian_paths is None else str(item.hamiltonian_paths)
        scc_text = ",".join(map(str, item.scc_sizes[:4]))
        if len(item.scc_sizes) > 4:
            scc_text += ",..."
        print(
            f"{item.row_label:<18} {item.time_kind:<8} {fmt(item.t):<13} "
            f"{item.gauge:<18} {item.info:<29} "
            f"{item.cyclic_triples:>5} {scc_text:<9} {h_text:>10}  {item.score_hist}"
        )
    print()


def print_hypothesis_checks(rows: tuple[SpeedRow, ...]) -> None:
    print_header("HYPOTHESIS CHECKS GENERATED BY THE PING-PONG")
    ledgers = {row.label: endpoint_ledger(row) for row in rows}
    cover14 = cover_family(14)
    cover18 = cover_family(18)

    print("1. Same local tax, different bridge fiber.")
    print(
        f"   n=14 exact={cover14.exact_size}, free={cover14.free_parts}; "
        f"n=18 exact={cover18.exact_size}, free={cover18.free_parts}."
    )
    print(
        "   So n=18 is not locally larger at the n-gate invoice; it is locally "
        "more rigid.  Its difficulty is the square 3-torsion debt channel."
    )
    print()

    print("2. Gap-debt product is conserved on both first-even ladders.")
    for n in (14, 18):
        labels = [f"n{n} row-parent", f"n{n} gate", f"n{n} double-gate"]
        products = tuple(ledgers[label].gap_debt_product for label in labels)
        debts = tuple(ledgers[label].unprotected for label in labels)
        gaps = tuple(ledgers[label].gap_ratio for label in labels)
        print(
            f"   n={n}: gaps={tuple(fmt(g) for g in gaps)} "
            f"debts={debts} products={tuple(fmt(p) for p in products)}"
        )
    print(
        "   The product is not a proof yet, but it is a stable certificate target: "
        "visible improvement is paid for by exact endpoint multiplication."
    )
    print()

    print("3. Tournament gauges separate rank collapse from edge-local structure.")
    print(
        "   The mixed-two-neighbor gauge is deliberately rank-like; pressure is "
        "edge-local.  If mixed slack collapses but pressure has cycles/SCCs, the "
        "proof object should be deletion relief, not scalar two-neighbor score."
    )
    print()

    print("4. Working next lemma.")
    print(
        "   First-even LRC rows 2m should be attacked by a bridge-fiber lemma: "
        "unit/half fan rows force one bridge fiber, and every bridge fiber either "
        "reopens a small denominator row or exports a positive endpoint frontier "
        "mass under the pressure tournament."
    )


def main() -> None:
    rows = speed_rows()
    print("n=14 / n=18 LRC Tournament Analysis ping-pong (codex-2026-06-01 S481)")
    print("All LRC interval and endpoint computations are exact over Fraction.\n")
    print_noise_cards()
    print_cover_comparison()
    print_endpoint_ledgers(rows)
    print_tournament_table(rows)
    print_hypothesis_checks(rows)


if __name__ == "__main__":
    main()
