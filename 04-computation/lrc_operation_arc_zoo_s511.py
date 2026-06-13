#!/usr/bin/env python3
"""
lrc_operation_arc_zoo_s511.py

codex-2026-06-01 S511

Runner-level Tournament Analysis for LRC with operation-grid arc criteria.

This script is a companion to S506/S509/S510.  S506 tested LRC arc rules as a
loneliness metric vector; S509 lifted the vertices to runner-pairs; S510 placed
the resulting movies over the A000568 tournament-isomorphism quotient.

Here the vertices are runners again, including the stationary observer, but the
arc switchboard is deliberately richer:

* phase and endpoint-distance gauges;
* two-nearest-neighbor guard gauges;
* x*2 dyadic-height and x+2 odd-core arithmetic rows;
* product-sum gates using (a-1)(b-1)=N-1;
* Goldbach / Lemoine gates: even N via p+q=N, odd N via p+2q=N;
* a majority operation bundle that lets these channels vote.

The goal is not to prove Goldbach/Lemoine or LRC.  It is to test whether these
operation channels create tournament fingerprints whose H values, score
sequences, tie counts, and strict SCC/triangle data are useful loneliness
coordinates.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import gcd, sqrt
from statistics import mean


ONE = Fraction(1, 1)
HALF = Fraction(1, 2)
H_LIMIT = 14


@dataclass(frozen=True)
class ArcCriterion:
    name: str
    channel: str
    switch: str
    tie_path: str = "increasing runner label"


@dataclass(frozen=True)
class Snapshot:
    criterion: ArcCriterion
    out: tuple[int, ...]
    strict_out: tuple[int, ...]
    ties: int
    strict_arcs: int


@dataclass(frozen=True)
class Fingerprint:
    vertex_count: int
    hamiltonian_paths: int | None
    score_hist: tuple[tuple[int, int], ...]
    score_width: int
    observer_out: int
    directed_triangles: int
    strict_triangles: int
    largest_strict_scc: int
    ties: int
    strict_arcs: int


@dataclass(frozen=True)
class Geometry:
    origin_margin: Fraction
    origin_ratio: Fraction
    max_gap: Fraction
    min_gap: Fraction
    safe_gap_count: int
    local_lonely_count: int
    close_pair_count: int


CRITERIA: tuple[ArcCriterion, ...] = (
    ArcCriterion(
        "phase_half",
        "phase",
        "i -> j iff j is clockwise from i by less than 1/2",
    ),
    ArcCriterion(
        "origin_clearance_rank",
        "endpoint",
        "farther from stationary observer points to closer",
    ),
    ArcCriterion(
        "endpoint_deficit_rank",
        "endpoint",
        "smaller deficit below 1/N points to larger deficit",
    ),
    ArcCriterion(
        "two_neighbor_guard",
        "two-neighbor",
        "larger sum of two nearest circular distances points to smaller",
    ),
    ArcCriterion(
        "dyadic_height",
        "x*2 branch",
        "higher v2(speed) branch points to lower branch",
    ),
    ArcCriterion(
        "odd_core_row",
        "x+2 row",
        "smaller odd core points to larger odd core",
    ),
    ArcCriterion(
        "same_chain_dyadic_phase",
        "x*2 over phase",
        "same odd-core chain uses dyadic branch; other pairs use phase",
    ),
    ArcCriterion(
        "adjacent_oddrow_phase",
        "x+2 over phase",
        "adjacent odd rows differ by 2; otherwise use phase",
    ),
    ArcCriterion(
        "product_sum_gate",
        "addition*multiplication",
        "near (a-1)(b-1)=N-1 product-sum pairs override phase",
    ),
    ArcCriterion(
        "goldbach_lemoine_gate",
        "prime additive gate",
        "even N uses p+q; odd N uses p+2q; otherwise phase",
    ),
    ArcCriterion(
        "danger_then_grid",
        "LRC threshold + operation grid",
        "pairs closer than 1/N use endpoint danger then grid order",
    ),
    ArcCriterion(
        "operation_bundle_vote",
        "bundle",
        "phase, endpoint, two-neighbor, dyadic, odd, product-sum, and prime gates vote",
    ),
)


SMALL_FAMILIES: tuple[tuple[str, int, tuple[int, ...]], ...] = (
    ("N5 initial", 5, (1, 2, 3, 4)),
    ("N5 primes", 5, (2, 3, 5, 7)),
    ("N6 initial", 6, (1, 2, 3, 4, 5)),
    ("N6 prime-lac", 6, (1, 5, 7, 11, 13)),
    ("N7 initial", 7, (1, 2, 3, 4, 5, 6)),
    ("N7 dyadic", 7, (1, 2, 4, 8, 16, 32)),
)


def fmt_frac(value: Fraction | int | None) -> str:
    if value is None:
        return "-"
    if isinstance(value, int):
        return str(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def fmt_float(value: float | None) -> str:
    if value is None:
        return f"{'-':>10s}"
    return f"{value:+10.3f}"


def mod1(value: Fraction) -> Fraction:
    return value % ONE


def clockwise_delta(a: Fraction, b: Fraction) -> Fraction:
    return mod1(b - a)


def circular_distance(a: Fraction, b: Fraction) -> Fraction:
    delta = clockwise_delta(a, b)
    return min(delta, ONE - delta)


def positions(speeds_with_observer: tuple[int, ...], t: Fraction) -> tuple[Fraction, ...]:
    return tuple(mod1(Fraction(v) * t) for v in speeds_with_observer)


def v2(value: int) -> int:
    value = abs(value)
    if value == 0:
        return -1
    out = 0
    while value % 2 == 0:
        out += 1
        value //= 2
    return out


def odd_core(value: int) -> int:
    value = abs(value)
    if value == 0:
        return 0
    return value >> max(v2(value), 0)


def grid_coord(speed: int) -> tuple[int, int]:
    if speed == 0:
        return (-1, 0)
    return (v2(speed), odd_core(speed))


def residue(speed: int, modulus: int) -> int:
    out = abs(speed) % modulus
    return modulus if out == 0 else out


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value == 2:
        return True
    if value % 2 == 0:
        return False
    p = 3
    while p * p <= value:
        if value % p == 0:
            return False
        p += 2
    return True


def phase_winner(pos: tuple[Fraction, ...], i: int, j: int) -> tuple[int | None, bool]:
    delta = clockwise_delta(pos[i], pos[j])
    if delta == 0 or delta == HALF:
        return (None, False)
    return (i if delta < HALF else j, True)


def rank_winner(values: tuple, i: int, j: int, larger_wins: bool = True) -> tuple[int | None, bool]:
    if values[i] == values[j]:
        return (None, False)
    if larger_wins:
        return (i if values[i] > values[j] else j, True)
    return (i if values[i] < values[j] else j, True)


def dyadic_winner(speeds: tuple[int, ...], i: int, j: int) -> tuple[int | None, bool]:
    hi, _ = grid_coord(speeds[i])
    hj, _ = grid_coord(speeds[j])
    if hi == hj:
        return (None, False)
    return (i if hi > hj else j, True)


def odd_core_winner(speeds: tuple[int, ...], i: int, j: int) -> tuple[int | None, bool]:
    _, ci = grid_coord(speeds[i])
    _, cj = grid_coord(speeds[j])
    if ci == cj:
        return (None, False)
    return (i if ci < cj else j, True)


def same_chain_winner(speeds: tuple[int, ...], i: int, j: int) -> tuple[int | None, bool]:
    hi, ci = grid_coord(speeds[i])
    hj, cj = grid_coord(speeds[j])
    if ci == 0 or ci != cj or hi == hj:
        return (None, False)
    return (i if hi > hj else j, True)


def adjacent_oddrow_winner(speeds: tuple[int, ...], i: int, j: int) -> tuple[int | None, bool]:
    _, ci = grid_coord(speeds[i])
    _, cj = grid_coord(speeds[j])
    if ci == 0 or cj == 0 or abs(ci - cj) != 2:
        return (None, False)
    return (i if ci < cj else j, True)


def product_sum_winner(
    speeds: tuple[int, ...],
    total_n: int,
    i: int,
    j: int,
) -> tuple[int | None, bool]:
    if speeds[i] == 0 or speeds[j] == 0:
        return (None, False)
    a = residue(speeds[i], total_n)
    b = residue(speeds[j], total_n)
    left = (a - 1) * (b - 1)
    defect = abs(left - (total_n - 1))
    if defect > max(1, total_n // 6):
        return (None, False)
    if a == b:
        return (None, False)
    # The smaller shifted factor is the additive side; the larger factor is
    # the multiplicative expansion side.
    return (i if a < b else j, True)


def goldbach_lemoine_winner(
    speeds: tuple[int, ...],
    total_n: int,
    i: int,
    j: int,
) -> tuple[int | None, bool]:
    if speeds[i] == 0 or speeds[j] == 0:
        return (None, False)
    a = residue(speeds[i], total_n)
    b = residue(speeds[j], total_n)
    if total_n % 2 == 0:
        if is_prime(a) and is_prime(b) and a + b == total_n and a != b:
            return (i if a < b else j, True)
        return (None, False)
    if is_prime(a) and is_prime(b) and a + 2 * b == total_n:
        return (i, True)
    if is_prime(a) and is_prime(b) and b + 2 * a == total_n:
        return (j, True)
    return (None, False)


def danger_grid_winner(
    speeds: tuple[int, ...],
    pos: tuple[Fraction, ...],
    total_n: int,
    i: int,
    j: int,
) -> tuple[int | None, bool]:
    threshold = Fraction(1, total_n)
    if circular_distance(pos[i], pos[j]) >= threshold:
        return (None, False)
    deficits = tuple(max(Fraction(0), threshold - circular_distance(pos[0], p)) for p in pos)
    if deficits[i] != deficits[j]:
        return (i if deficits[i] > deficits[j] else j, True)
    winner, strict = dyadic_winner(speeds, i, j)
    if strict:
        return (winner, strict)
    return odd_core_winner(speeds, i, j)


def two_neighbor_scores(pos: tuple[Fraction, ...]) -> tuple[tuple[Fraction, Fraction], ...]:
    out: list[tuple[Fraction, Fraction]] = []
    for i, p in enumerate(pos):
        ds = sorted(circular_distance(p, q) for j, q in enumerate(pos) if j != i)
        out.append((ds[0] + ds[1], ds[0]))
    return tuple(out)


def decide_pair(
    criterion: ArcCriterion,
    speeds: tuple[int, ...],
    pos: tuple[Fraction, ...],
    i: int,
    j: int,
) -> tuple[int | None, bool]:
    total_n = len(speeds)
    threshold = Fraction(1, total_n)

    if criterion.name == "phase_half":
        return phase_winner(pos, i, j)

    if criterion.name == "origin_clearance_rank":
        values = tuple(circular_distance(pos[0], p) for p in pos)
        return rank_winner(values, i, j, larger_wins=True)

    if criterion.name == "endpoint_deficit_rank":
        values = tuple(max(Fraction(0), threshold - circular_distance(pos[0], p)) for p in pos)
        return rank_winner(values, i, j, larger_wins=False)

    if criterion.name == "two_neighbor_guard":
        return rank_winner(two_neighbor_scores(pos), i, j, larger_wins=True)

    if criterion.name == "dyadic_height":
        return dyadic_winner(speeds, i, j)

    if criterion.name == "odd_core_row":
        return odd_core_winner(speeds, i, j)

    if criterion.name == "same_chain_dyadic_phase":
        winner, strict = same_chain_winner(speeds, i, j)
        return (winner, strict) if strict else phase_winner(pos, i, j)

    if criterion.name == "adjacent_oddrow_phase":
        winner, strict = adjacent_oddrow_winner(speeds, i, j)
        return (winner, strict) if strict else phase_winner(pos, i, j)

    if criterion.name == "product_sum_gate":
        winner, strict = product_sum_winner(speeds, total_n, i, j)
        return (winner, strict) if strict else phase_winner(pos, i, j)

    if criterion.name == "goldbach_lemoine_gate":
        winner, strict = goldbach_lemoine_winner(speeds, total_n, i, j)
        return (winner, strict) if strict else phase_winner(pos, i, j)

    if criterion.name == "danger_then_grid":
        winner, strict = danger_grid_winner(speeds, pos, total_n, i, j)
        return (winner, strict) if strict else phase_winner(pos, i, j)

    if criterion.name == "operation_bundle_vote":
        votes = Counter()
        for helper in (
            lambda: phase_winner(pos, i, j),
            lambda: dyadic_winner(speeds, i, j),
            lambda: odd_core_winner(speeds, i, j),
            lambda: same_chain_winner(speeds, i, j),
            lambda: product_sum_winner(speeds, total_n, i, j),
            lambda: goldbach_lemoine_winner(speeds, total_n, i, j),
            lambda: danger_grid_winner(speeds, pos, total_n, i, j),
            lambda: rank_winner(two_neighbor_scores(pos), i, j, larger_wins=True),
        ):
            winner, strict = helper()
            if strict and winner is not None:
                votes[winner] += 1
        if votes[i] == votes[j]:
            return (None, False)
        return (i if votes[i] > votes[j] else j, True)

    raise ValueError(f"unknown criterion {criterion.name}")


def snapshot(criterion: ArcCriterion, speeds: tuple[int, ...], t: Fraction) -> Snapshot:
    pos = positions(speeds, t)
    n = len(speeds)
    out = [0] * n
    strict_out = [0] * n
    ties = 0
    strict_arcs = 0
    for i, j in combinations(range(n), 2):
        winner, strict = decide_pair(criterion, speeds, pos, i, j)
        if winner is None:
            winner = i
            ties += 1
        else:
            strict_arcs += 1
        loser = j if winner == i else i
        out[winner] |= 1 << loser
        if strict:
            strict_out[winner] |= 1 << loser
    return Snapshot(criterion, tuple(out), tuple(strict_out), ties, strict_arcs)


def hamiltonian_paths(out: tuple[int, ...]) -> int | None:
    n = len(out)
    if n > H_LIMIT:
        return None
    dp = [defaultdict(int) for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last, count in list(dp[mask].items()):
            available = out[last] & ~mask
            while available:
                bit = available & -available
                nxt = bit.bit_length() - 1
                dp[mask | bit][nxt] += count
                available ^= bit
    return sum(dp[(1 << n) - 1].values())


def directed_triangles(out: tuple[int, ...]) -> int:
    total = 0
    for a, b, c in combinations(range(len(out)), 3):
        degs = []
        for u in (a, b, c):
            degs.append(sum(1 for v in (a, b, c) if v != u and ((out[u] >> v) & 1)))
        if sorted(degs) == [1, 1, 1]:
            total += 1
    return total


def largest_scc(out: tuple[int, ...]) -> int:
    n = len(out)
    rev = [0] * n
    for u in range(n):
        bits = out[u]
        while bits:
            bit = bits & -bits
            v = bit.bit_length() - 1
            rev[v] |= 1 << u
            bits ^= bit

    seen = [False] * n
    order: list[int] = []

    def dfs1(u: int) -> None:
        seen[u] = True
        bits = out[u]
        while bits:
            bit = bits & -bits
            v = bit.bit_length() - 1
            if not seen[v]:
                dfs1(v)
            bits ^= bit
        order.append(u)

    def dfs2(u: int, comp: list[int]) -> None:
        seen[u] = True
        comp.append(u)
        bits = rev[u]
        while bits:
            bit = bits & -bits
            v = bit.bit_length() - 1
            if not seen[v]:
                dfs2(v, comp)
            bits ^= bit

    for v in range(n):
        if not seen[v]:
            dfs1(v)
    seen = [False] * n
    best = 0
    for v in reversed(order):
        if not seen[v]:
            comp: list[int] = []
            dfs2(v, comp)
            best = max(best, len(comp))
    return best


def fingerprint(snap: Snapshot) -> Fingerprint:
    scores = tuple(mask.bit_count() for mask in snap.out)
    hist = tuple(sorted(Counter(scores).items()))
    return Fingerprint(
        vertex_count=len(snap.out),
        hamiltonian_paths=hamiltonian_paths(snap.out),
        score_hist=hist,
        score_width=max(scores) - min(scores),
        observer_out=scores[0],
        directed_triangles=directed_triangles(snap.out),
        strict_triangles=directed_triangles(snap.strict_out),
        largest_strict_scc=largest_scc(snap.strict_out),
        ties=snap.ties,
        strict_arcs=snap.strict_arcs,
    )


def circular_gap_data(pos: tuple[Fraction, ...]) -> tuple[tuple[Fraction, ...], dict[int, tuple[Fraction, Fraction]]]:
    ordered = sorted((p, i) for i, p in enumerate(pos))
    n = len(ordered)
    gaps: list[Fraction] = []
    around: dict[int, tuple[Fraction, Fraction]] = {}
    for k, (p, i) in enumerate(ordered):
        q, _ = ordered[(k + 1) % n]
        gap = mod1(q - p)
        gaps.append(gap)
    for k, (_, i) in enumerate(ordered):
        left = gaps[(k - 1) % n]
        right = gaps[k]
        around[i] = (left, right)
    return tuple(gaps), around


def geometry(speeds_with_observer: tuple[int, ...], t: Fraction) -> Geometry:
    total_n = len(speeds_with_observer)
    threshold = Fraction(1, total_n)
    pos = positions(speeds_with_observer, t)
    moving_dists = [circular_distance(pos[0], pos[i]) for i in range(1, total_n)]
    origin_margin = min(moving_dists) if moving_dists else Fraction(0)
    gaps, around = circular_gap_data(pos)
    close_pair_count = sum(
        1
        for i, j in combinations(range(total_n), 2)
        if circular_distance(pos[i], pos[j]) < threshold
    )
    return Geometry(
        origin_margin=origin_margin,
        origin_ratio=origin_margin / threshold,
        max_gap=max(gaps),
        min_gap=min(gaps),
        safe_gap_count=sum(1 for g in gaps if g >= threshold),
        local_lonely_count=sum(1 for left, right in around.values() if left >= threshold and right >= threshold),
        close_pair_count=close_pair_count,
    )


def wall_midpoints(speeds: tuple[int, ...], total_n: int) -> tuple[Fraction, ...]:
    walls = {Fraction(0), Fraction(1)}
    speeds_with_observer = (0,) + speeds

    for s in speeds:
        s_abs = abs(s)
        for m in range(s_abs):
            walls.add(Fraction(m, s_abs))
            walls.add(mod1((Fraction(m) + Fraction(1, total_n)) / s_abs))
            walls.add(mod1((Fraction(m) - Fraction(1, total_n)) / s_abs))

    for i, j in combinations(range(len(speeds_with_observer)), 2):
        gap = abs(speeds_with_observer[i] - speeds_with_observer[j])
        if gap == 0:
            continue
        for m in range(gap):
            walls.add(Fraction(m, gap))
            walls.add(Fraction(2 * m + 1, 2 * gap))

    ordered = sorted(w for w in walls if 0 <= w <= 1)
    mids: list[Fraction] = []
    for a, b in zip(ordered, ordered[1:]):
        if a != b:
            mids.append((a + b) / 2)
    if ordered[0] == 0 and ordered[-1] == 1:
        return tuple(m for m in mids if 0 <= m < 1)
    return tuple(mids)


def best_origin_time(speeds: tuple[int, ...], total_n: int) -> tuple[Fraction, Geometry]:
    candidates = set(wall_midpoints(speeds, total_n))
    candidates.add(Fraction(1, 2 * total_n))
    candidates.add(Fraction(1, total_n))
    speeds_with_observer = (0,) + speeds
    best_t = Fraction(0)
    best_g = geometry(speeds_with_observer, best_t)
    for t in candidates:
        g = geometry(speeds_with_observer, t)
        key = (g.origin_margin, g.local_lonely_count, g.max_gap, -t)
        best_key = (best_g.origin_margin, best_g.local_lonely_count, best_g.max_gap, -best_t)
        if key > best_key:
            best_t, best_g = t, g
    return best_t, best_g


def ranks(values: list[Fraction | int | float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    out = [0.0] * len(values)
    k = 0
    while k < len(order):
        ell = k + 1
        while ell < len(order) and values[order[ell]] == values[order[k]]:
            ell += 1
        rank = (k + 1 + ell) / 2
        for idx in order[k:ell]:
            out[idx] = rank
        k = ell
    return out


def spearman(xs: list[Fraction | int | float], ys: list[Fraction | int | float]) -> float | None:
    if len(xs) < 3 or len(set(xs)) < 2 or len(set(ys)) < 2:
        return None
    rx = ranks(xs)
    ry = ranks(ys)
    mx = mean(rx)
    my = mean(ry)
    num = sum((x - mx) * (y - my) for x, y in zip(rx, ry))
    denx = sqrt(sum((x - mx) ** 2 for x in rx))
    deny = sqrt(sum((y - my) ** 2 for y in ry))
    if denx == 0 or deny == 0:
        return None
    return num / (denx * deny)


def mean_defined(values: list[float | None]) -> float | None:
    clean = [v for v in values if v is not None]
    if not clean:
        return None
    return mean(clean)


def score_small_clocks() -> dict[str, dict[str, float | None]]:
    per_criterion: dict[str, dict[str, list[float | None]]] = {
        c.name: defaultdict(list) for c in CRITERIA
    }
    print("== small exact clock scorecard ==")
    for label, total_n, speeds in SMALL_FAMILIES:
        samples = wall_midpoints(speeds, total_n)
        speeds_with_observer = (0,) + speeds
        print(f"{label}: {len(samples)} chamber midpoints")
        geoms = [geometry(speeds_with_observer, t) for t in samples]
        targets = {
            "origin": [g.origin_ratio for g in geoms],
            "max_gap": [g.max_gap for g in geoms],
            "close_pairs": [g.close_pair_count for g in geoms],
            "local_lonely": [g.local_lonely_count for g in geoms],
        }
        for criterion in CRITERIA:
            fps = [fingerprint(snapshot(criterion, speeds_with_observer, t)) for t in samples]
            metrics = {
                "H_to_max_gap": [fp.hamiltonian_paths or 0 for fp in fps],
                "c3_to_max_gap": [fp.directed_triangles for fp in fps],
                "observer_to_origin": [fp.observer_out for fp in fps],
                "ties_to_close_pairs": [fp.ties for fp in fps],
                "width_to_local_lonely": [fp.score_width for fp in fps],
            }
            per_criterion[criterion.name]["H_to_max_gap"].append(
                spearman(metrics["H_to_max_gap"], targets["max_gap"])
            )
            per_criterion[criterion.name]["c3_to_max_gap"].append(
                spearman(metrics["c3_to_max_gap"], targets["max_gap"])
            )
            per_criterion[criterion.name]["observer_to_origin"].append(
                spearman(metrics["observer_to_origin"], targets["origin"])
            )
            per_criterion[criterion.name]["ties_to_close_pairs"].append(
                spearman(metrics["ties_to_close_pairs"], targets["close_pairs"])
            )
            per_criterion[criterion.name]["width_to_local_lonely"].append(
                spearman(metrics["width_to_local_lonely"], targets["local_lonely"])
            )

    reduced: dict[str, dict[str, float | None]] = {}
    for name, metric_map in per_criterion.items():
        reduced[name] = {metric: mean_defined(vals) for metric, vals in metric_map.items()}

    print()
    print(
        "criterion                         H~maxgap  c3~maxgap  obs~origin  ties~close  width~local"
    )
    for criterion in CRITERIA:
        row = reduced[criterion.name]
        print(
            f"{criterion.name:33s}"
            f"{fmt_float(row.get('H_to_max_gap'))}"
            f"{fmt_float(row.get('c3_to_max_gap'))}"
            f"{fmt_float(row.get('observer_to_origin'))}"
            f"{fmt_float(row.get('ties_to_close_pairs'))}"
            f"{fmt_float(row.get('width_to_local_lonely'))}"
        )
    print()
    return reduced


def ladder(total_n: int, scale: int, skip: int) -> tuple[int, ...]:
    speeds = {1}
    for q in range(1, total_n):
        if q != skip:
            speeds.add(scale * q)
    return tuple(sorted(speeds))


def fmt_hist(hist: tuple[tuple[int, int], ...]) -> str:
    return ",".join(f"{score}:{count}" for score, count in hist)


def hard_row_report() -> None:
    rows = (
        ("n14 initial", 14, tuple(range(1, 14))),
        ("n14 row-parent 7*", 14, ladder(14, 7, 6)),
        ("n14 gate 14*", 14, ladder(14, 14, 6)),
        ("n18 initial", 18, tuple(range(1, 18))),
        ("n18 row-parent 9*", 18, ladder(18, 9, 8)),
        ("n18 gate 18*", 18, ladder(18, 18, 8)),
    )
    selected = {
        "phase_half",
        "origin_clearance_rank",
        "endpoint_deficit_rank",
        "two_neighbor_guard",
        "same_chain_dyadic_phase",
        "product_sum_gate",
        "goldbach_lemoine_gate",
        "danger_then_grid",
        "operation_bundle_vote",
    }
    print("== selected n=14 / n=18 hard-row fingerprints ==")
    for label, total_n, speeds in rows:
        speeds_with_observer = (0,) + speeds
        best_t, best_g = best_origin_time(speeds, total_n)
        half_g = geometry(speeds_with_observer, Fraction(1, 2 * total_n))
        unit_g = geometry(speeds_with_observer, Fraction(1, total_n))
        print()
        print(
            f"{label}: total N={total_n}, moving={len(speeds)}, "
            f"best_t={fmt_frac(best_t)}, best_origin={fmt_frac(best_g.origin_margin)} "
            f"(ratio {fmt_frac(best_g.origin_ratio)})"
        )
        print(
            "  geometry half/unit/best ratios="
            f"{fmt_frac(half_g.origin_ratio)}/"
            f"{fmt_frac(unit_g.origin_ratio)}/"
            f"{fmt_frac(best_g.origin_ratio)}; "
            f"best safe_gaps={best_g.safe_gap_count}, local_lonely={best_g.local_lonely_count}, "
            f"close_pairs={best_g.close_pair_count}, max_gap={fmt_frac(best_g.max_gap)}"
        )
        print("  criterion                         H          c3    width obs ties strictSCC hist")
        for criterion in CRITERIA:
            if criterion.name not in selected:
                continue
            fp = fingerprint(snapshot(criterion, speeds_with_observer, best_t))
            h = "-" if fp.hamiltonian_paths is None else str(fp.hamiltonian_paths)
            print(
                f"  {criterion.name:33s}"
                f"{h:>10s} {fp.directed_triangles:5d} {fp.score_width:6d}"
                f" {fp.observer_out:3d} {fp.ties:4d} {fp.largest_strict_scc:9d} "
                f"{fmt_hist(fp.score_hist)}"
            )
    print()


def goldbach_pairs(total: int) -> tuple[tuple[int, int], ...]:
    out = []
    for p in range(2, total // 2 + 1):
        q = total - p
        if p <= q and is_prime(p) and is_prime(q):
            out.append((p, q))
    return tuple(out)


def lemoine_pairs(total: int) -> tuple[tuple[int, int], ...]:
    out = []
    for p in range(2, total):
        q2 = total - p
        if q2 % 2 == 0:
            q = q2 // 2
            if is_prime(p) and is_prime(q):
                out.append((p, q))
    return tuple(out)


def product_sum_two_core(total_n: int) -> tuple[tuple[int, int], ...]:
    target = total_n - 1
    out = []
    for d in range(1, target + 1):
        if target % d == 0:
            a = d + 1
            b = target // d + 1
            if a <= b:
                out.append((a, b))
    return tuple(out)


def operation_grid_report() -> None:
    print("== addition / multiplication operation grid ==")
    print(
        "N   2^h*odd   additive prime gate       two-core product-sum gates (a-1)(b-1)=N-1"
    )
    for total_n in range(2, 31):
        h = v2(total_n)
        core = odd_core(total_n)
        if total_n % 2 == 0:
            gate = "Goldbach " + str(goldbach_pairs(total_n)[:4])
        else:
            gate = "Lemoine " + str(lemoine_pairs(total_n)[:4])
        ps = product_sum_two_core(total_n)
        print(f"{total_n:2d}  2^{h}*{core:<3d}  {gate:29s} {ps[:5]}")
    print()
    print("selected operation readings:")
    for total_n in (14, 15, 16, 18, 21, 24):
        print(
            f"  N={total_n}: odd_core={odd_core(total_n)}, v2={v2(total_n)}, "
            f"Goldbach={goldbach_pairs(total_n) if total_n % 2 == 0 else '-'}, "
            f"Lemoine={lemoine_pairs(total_n) if total_n % 2 else '-'}, "
            f"product_sum={product_sum_two_core(total_n)}"
        )
    print()


def criterion_manifest() -> None:
    print("== arc criteria manifest ==")
    for criterion in CRITERIA:
        print(f"{criterion.name:25s} | {criterion.channel:27s} | {criterion.switch}")
    print()


def interpretation() -> None:
    print("== interpretation ==")
    print(
        "A000568 is the unmarked quotient of complete binary relations.  Its Burnside "
        "formula only sees odd permutation cycles; the LRC operation-grid lift keeps "
        "a parallel odd/even memory by splitting each speed into an odd x+2 row and "
        "a vertical x*2 branch."
    )
    print(
        "Addition supplies fixed-sum and adjacent-row gates: Goldbach p+q=N for even "
        "N and Lemoine p+2q=N for odd N.  Multiplication supplies branch/factor gates: "
        "dyadic height, odd cores, divisibility, and the shifted product-sum equation "
        "(a-1)(b-1)=N-1."
    )
    print(
        "The useful loneliness object is therefore not a scalar H.  It is a marked "
        "fiber over A000568 carrying phase H, observer score, close-pair tie rate, "
        "two-neighbor width, strict pressure SCCs, and operation-grid labels."
    )
    print()


def main() -> None:
    criterion_manifest()
    score_small_clocks()
    hard_row_report()
    operation_grid_report()
    interpretation()


if __name__ == "__main__":
    main()
