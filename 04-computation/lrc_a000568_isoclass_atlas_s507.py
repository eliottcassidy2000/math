#!/usr/bin/env python3
"""
lrc_a000568_isoclass_atlas_s507.py

codex-2026-06-01 S507

Explore the user's hypothesis that the Lonely Runner Conjecture should be
read as an A000568-like tournament isomorphism-class problem.

The main refinement tested here is that LRC is not just an unrooted
tournament quotient.  The stationary runner, endpoint masks, and pressure
labels make it a rooted/marked quotient stack over the A000568 quotient.

Tournament Analysis declaration:

* pairwise observable: half-turn phase, endpoint safety, and deletion pressure;
* switch/gauge: phase_half, safe_phase_gate, and pressure_k2 from S506;
* tie Hamiltonian path: increasing speed/order label.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from importlib.machinery import SourceFileLoader
from itertools import combinations, permutations
from math import comb, factorial, gcd, log2
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S506 = SourceFileLoader(
    "lrc_tournament_gauge_lab_s506",
    str(ROOT / "04-computation" / "lrc_tournament_gauge_lab_s506.py"),
).load_module()

ONE = Fraction(1, 1)
HALF = Fraction(1, 2)


@dataclass(frozen=True)
class ClassCensus:
    n: int
    a000568: int
    rooted: int
    h_buckets: int
    score_buckets: int
    structural_buckets: int
    min_h: int
    max_h: int


@dataclass(frozen=True)
class ClockQuotient:
    n: int
    cells: int
    labeled_states: int
    unrooted_states: int
    rooted_states: int
    h_buckets: int
    a000568: int
    rooted_total: int | None


def fmt_frac(value: Fraction | int | None) -> str:
    if value is None:
        return "-"
    if isinstance(value, int):
        return str(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def fmt_float(value: float, places: int = 4) -> str:
    return f"{value:.{places}f}"


def partitions(n: int, max_part: int | None = None):
    if max_part is None:
        max_part = n
    if n == 0:
        yield ()
        return
    for first in range(min(n, max_part), 0, -1):
        for rest in partitions(n - first, first):
            yield (first,) + rest


def cycle_type_count(n: int, ct: tuple[int, ...]) -> int:
    counts = Counter(ct)
    out = factorial(n)
    for length, mult in counts.items():
        out //= (length**mult) * factorial(mult)
    return out


def fixed_tournaments_for_cycle_type(ct: tuple[int, ...]) -> int:
    if any(part % 2 == 0 for part in ct):
        return 0
    exp = sum((part - 1) // 2 for part in ct)
    for i, a in enumerate(ct):
        for b in ct[i + 1 :]:
            exp += gcd(a, b)
    return 2**exp


@lru_cache(maxsize=None)
def a000568(n: int) -> int:
    if n <= 2:
        return 1
    total = 0
    for ct in partitions(n):
        total += cycle_type_count(n, ct) * fixed_tournaments_for_cycle_type(ct)
    return total // factorial(n)


def odd_partition_fraction(n: int) -> Fraction:
    all_parts = 0
    odd_parts = 0
    for ct in partitions(n):
        all_parts += 1
        if all(part % 2 == 1 for part in ct):
            odd_parts += 1
    return Fraction(odd_parts, all_parts)


@lru_cache(maxsize=None)
def pair_list(n: int) -> tuple[tuple[int, int], ...]:
    return tuple(combinations(range(n), 2))


@lru_cache(maxsize=None)
def pair_index(n: int) -> dict[tuple[int, int], int]:
    return {pair: idx for idx, pair in enumerate(pair_list(n))}


def edge(bits: int, n: int, i: int, j: int) -> bool:
    if i < j:
        return bool((bits >> pair_index(n)[(i, j)]) & 1)
    return not bool((bits >> pair_index(n)[(j, i)]) & 1)


@lru_cache(maxsize=None)
def perm_group(n: int, rooted: bool) -> tuple[tuple[int, ...], ...]:
    if not rooted:
        return tuple(permutations(range(n)))
    return tuple((0,) + tail for tail in permutations(range(1, n)))


def relabel_bits(bits: int, n: int, perm: tuple[int, ...]) -> int:
    out = 0
    for idx, (i, j) in enumerate(pair_list(n)):
        if edge(bits, n, perm[i], perm[j]):
            out |= 1 << idx
    return out


def canonical_bits(bits: int, n: int, rooted: bool = False) -> int:
    return min(relabel_bits(bits, n, perm) for perm in perm_group(n, rooted))


def enumerate_orbit_reps(n: int, rooted: bool = False) -> tuple[int, ...]:
    seen: set[int] = set()
    reps: list[int] = []
    for bits in range(1 << comb(n, 2)):
        if bits in seen:
            continue
        orbit = {relabel_bits(bits, n, perm) for perm in perm_group(n, rooted)}
        seen.update(orbit)
        reps.append(min(orbit))
    return tuple(sorted(reps))


def masks_from_bits(bits: int, n: int) -> tuple[int, ...]:
    masks = [0] * n
    for idx, (i, j) in enumerate(pair_list(n)):
        if (bits >> idx) & 1:
            masks[i] |= 1 << j
        else:
            masks[j] |= 1 << i
    return tuple(masks)


@lru_cache(maxsize=None)
def hamiltonian_paths_from_masks(n: int, masks: tuple[int, ...]) -> int:
    full = (1 << n) - 1
    dp = [0] * ((1 << n) * n)
    for v in range(n):
        dp[((1 << v) * n) + v] = 1
    for mask in range(1 << n):
        base = mask * n
        for last in range(n):
            value = dp[base + last]
            if not value:
                continue
            allowed = masks[last] & ~mask
            while allowed:
                bit = allowed & -allowed
                nxt = bit.bit_length() - 1
                dp[((mask | bit) * n) + nxt] += value
                allowed -= bit
    return sum(dp[full * n : (full + 1) * n])


def score_sequence_from_masks(masks: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(mask.bit_count() for mask in masks))


def cyclic_triples_from_scores(n: int, scores: tuple[int, ...]) -> int:
    return comb(n, 3) - sum(comb(score, 2) for score in scores)


def scc_sizes_from_masks(n: int, masks: tuple[int, ...]) -> tuple[int, ...]:
    graph = [[j for j in range(n) if (masks[i] >> j) & 1] for i in range(n)]
    rev = [[] for _ in range(n)]
    for i, outs in enumerate(graph):
        for j in outs:
            rev[j].append(i)

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

    seen.clear()
    sizes: list[int] = []

    def rdfs(v: int) -> int:
        seen.add(v)
        size = 1
        for nxt in rev[v]:
            if nxt not in seen:
                size += rdfs(nxt)
        return size

    for v in reversed(order):
        if v not in seen:
            sizes.append(rdfs(v))
    return tuple(sorted(sizes, reverse=True))


def class_census(n: int) -> ClassCensus:
    reps = enumerate_orbit_reps(n, rooted=False)
    rooted_reps = enumerate_orbit_reps(n, rooted=True)
    h_values: Counter[int] = Counter()
    score_values: Counter[tuple[int, ...]] = Counter()
    structural_values: Counter[tuple[tuple[int, ...], int, tuple[int, ...]]] = Counter()

    for rep in reps:
        masks = masks_from_bits(rep, n)
        h = hamiltonian_paths_from_masks(n, masks)
        scores = score_sequence_from_masks(masks)
        c3 = cyclic_triples_from_scores(n, scores)
        scc = scc_sizes_from_masks(n, masks)
        h_values[h] += 1
        score_values[scores] += 1
        structural_values[(scores, c3, scc)] += 1

    return ClassCensus(
        n=n,
        a000568=len(reps),
        rooted=len(rooted_reps),
        h_buckets=len(h_values),
        score_buckets=len(score_values),
        structural_buckets=len(structural_values),
        min_h=min(h_values),
        max_h=max(h_values),
    )


def circle(value: Fraction) -> Fraction:
    return value % ONE


def phase_walls(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    walls: set[Fraction] = set()
    for a, b in combinations(speeds, 2):
        d = abs(a - b)
        if d == 0:
            continue
        for m in range(2 * d):
            walls.add(Fraction(m, 2 * d))
    return tuple(sorted(walls))


def clock_cell_midpoints(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    walls = phase_walls(speeds)
    out: list[Fraction] = []
    for idx, lo in enumerate(walls):
        hi = walls[(idx + 1) % len(walls)]
        if hi <= lo:
            hi += ONE
        out.append(circle(lo + (hi - lo) / 2))
    return tuple(out)


def phase_bits(speeds: tuple[int, ...], t: Fraction) -> int:
    pos = tuple(circle(speed * t) for speed in speeds)
    n = len(speeds)
    bits = 0
    for idx, (i, j) in enumerate(pair_list(n)):
        delta = circle(pos[j] - pos[i])
        if delta == 0 or delta == HALF:
            winner = i
        elif delta < HALF:
            winner = i
        else:
            winner = j
        if winner == i:
            bits |= 1 << idx
    return bits


def clock_quotient(n: int, rooted_total: int | None) -> ClockQuotient:
    speeds = tuple(range(n))
    cells = clock_cell_midpoints(speeds)
    labeled: set[int] = set()
    unrooted: set[int] = set()
    rooted: set[int] = set()
    h_values: set[int] = set()

    for t in cells:
        bits = phase_bits(speeds, t)
        labeled.add(bits)
        unrooted.add(canonical_bits(bits, n, rooted=False))
        rooted.add(canonical_bits(bits, n, rooted=True))
        h_values.add(hamiltonian_paths_from_masks(n, masks_from_bits(bits, n)))

    return ClockQuotient(
        n=n,
        cells=len(cells),
        labeled_states=len(labeled),
        unrooted_states=len(unrooted),
        rooted_states=len(rooted),
        h_buckets=len(h_values),
        a000568=a000568(n),
        rooted_total=rooted_total,
    )


def h_ratio(h: int, n: int) -> Fraction:
    return Fraction(h * (2 ** (n - 1)), factorial(n))


def print_burnside_table() -> None:
    print("A000568/Burnside quotient scale")
    print("=" * 112)
    print(
        f"{'n':>2} {'2^C(n,2)':>16} {'A000568':>14} {'avg orbit':>13} "
        f"{'label bits':>10} {'odd part frac':>14}"
    )
    print("-" * 112)
    for n in range(3, 11):
        labeled = 2 ** comb(n, 2)
        a = a000568(n)
        avg_orbit = Fraction(labeled, a)
        label_bits = comb(n, 2) - log2(a)
        print(
            f"{n:>2} {labeled:>16} {a:>14} {fmt_float(float(avg_orbit), 2):>13} "
            f"{fmt_float(label_bits, 3):>10} {fmt_frac(odd_partition_fraction(n)):>14}"
        )
    print()
    print(
        "Reading: A000568 is not just a count.  It is the quotient cost of "
        "forgetting vertex labels from binary pair decisions."
    )
    print(
        "The odd-part column is the Burnside fixed-point filter: even cycles "
        "kill orientation fixed points."
    )


def print_small_census(censuses: tuple[ClassCensus, ...]) -> None:
    print()
    print("Small exact tournament quotient census")
    print("=" * 112)
    print(
        f"{'n':>2} {'A000568':>8} {'rooted':>8} {'root/A':>8} "
        f"{'H buckets':>9} {'score':>7} {'struct':>7} {'H range':>16}"
    )
    print("-" * 112)
    for row in censuses:
        print(
            f"{row.n:>2} {row.a000568:>8} {row.rooted:>8} "
            f"{fmt_float(row.rooted / row.a000568, 2):>8} "
            f"{row.h_buckets:>9} {row.score_buckets:>7} "
            f"{row.structural_buckets:>7} {row.min_h:>7}..{row.max_h}"
        )
    print()
    print(
        "Reading: H and score sequence are useful quotient coordinates, but "
        "they are far coarser than isomorphism.  Rooting already expands the "
        "state space by a factor between 1 and n."
    )


def print_clock_table(rows: tuple[ClockQuotient, ...]) -> None:
    print()
    print("Initial LRC half-turn clock as a path in the quotient")
    print("=" * 112)
    print(
        f"{'n':>2} {'cells':>6} {'labeled':>8} {'unrooted':>9} {'rooted':>8} "
        f"{'H vals':>6} {'A000568':>8} {'rooted all':>10} {'unroot %':>9} {'root %':>8}"
    )
    print("-" * 112)
    for row in rows:
        rooted_all = "-" if row.rooted_total is None else str(row.rooted_total)
        unroot_pct = row.unrooted_states / row.a000568
        root_pct = None if row.rooted_total is None else row.rooted_states / row.rooted_total
        root_pct_s = "-" if root_pct is None else f"{fmt_float(100 * root_pct, 2)}%"
        print(
            f"{row.n:>2} {row.cells:>6} {row.labeled_states:>8} "
            f"{row.unrooted_states:>9} {row.rooted_states:>8} {row.h_buckets:>6} "
            f"{row.a000568:>8} {rooted_all:>10} {fmt_float(100 * unroot_pct, 2):>8}% "
            f"{root_pct_s:>8}"
        )
    print()
    print(
        "Reading: the LRC clock is a thin deterministic walk in the tournament "
        "quotient, not a census of all A000568 classes.  Rooting by the "
        "stationary runner recovers information that the unrooted quotient "
        "deliberately forgets."
    )


def shadow_key(adj: list[list[bool]], h: int, ties: int) -> tuple:
    scores = S506.score_sequence(adj)
    return (
        scores,
        h,
        S506.cyclic_triples(adj),
        S506.largest_scc(adj),
        sum(adj[0]),
        ties,
    )


def print_hard_row_shadows() -> None:
    print()
    print("n=14/n=18 hard-row rooted isoclass shadows")
    print("=" * 112)
    print(
        f"{'gauge':<16} {'row':<12} {'sid':>4} {'Hratio':>8} {'root':>5} "
        f"{'width':>5} {'c3':>5} {'SCC':>4} {'debt':>6} {'gap*debt':>10}"
    )
    print("-" * 112)

    rows = S506.build_rows()
    gauges = ("phase_half", "safe_phase_gate", "pressure_k2")
    for gauge in gauges:
        ids: dict[tuple, str] = {}
        next_id = 1
        for row in rows:
            adj, ties = S506.tournament_for_gauge(gauge, row)
            n = len(adj)
            h = S506.hamiltonian_paths(n, S506.row_masks(adj))
            key = shadow_key(adj, h, ties)
            if key not in ids:
                ids[key] = f"S{next_id}"
                next_id += 1
            scores = S506.score_sequence(adj)
            print(
                f"{gauge:<16} {row.label:<12} {ids[key]:>4} "
                f"{fmt_float(float(h_ratio(h, n)), 4):>8} {sum(adj[0]):>5} "
                f"{(max(scores) - min(scores)):>5} {S506.cyclic_triples(adj):>5} "
                f"{S506.largest_scc(adj):>4} {row.endpoint_debt:>6} "
                f"{fmt_frac(row.endpoint_product):>10}"
            )
        print("-" * 112)
    print()
    print(
        "Reading: after the first gate, the hard ladders often keep the same "
        "coarse rooted shadow while endpoint debt continues to move.  This is "
        "exactly the shape of a fiber problem over an A000568-like base."
    )


BRIDGE_ROUTES = (
    (
        "unrooted base",
        "ordinary tournament quotient, counted by A000568",
        "phase_half H lives here as a coarse class height",
        "cannot see the stationary observer or endpoint debt",
    ),
    (
        "rooted lift",
        "tournament quotient with one marked vertex",
        "observer 0 is fixed; root score becomes meaningful",
        "probably the first faithful LRC quotient",
    ),
    (
        "marked lift",
        "rooted class plus safe mask / endpoint labels",
        "safe_phase_gate and endpoint debt live in this fiber",
        "needs a Burnside formula with colored vertices/arcs",
    ),
    (
        "pressure lift",
        "rooted class plus deletion-relief orientation",
        "pressure_k2 SCCs model all-protected endpoint cores",
        "H is intentionally small; use SCC/source/sink data",
    ),
    (
        "fiber recursion",
        "A000568 roughly multiplies by 2^(n-1)/n",
        "LRC gates add binary pair data then quotient by time/speed symmetry",
        "first-gate H drop plus debt plateau is the diagnostic",
    ),
    (
        "fixed-point killing",
        "even permutation cycles contribute zero to Burnside",
        "unit endpoint protection requires divisibility/gate support",
        "turn endpoint impossibility into a fixed-point-free lemma",
    ),
    (
        "bad-class emptiness",
        "prove no isoclass has a forbidden invariant vector",
        "counterexample class would have phase alarm plus pressure SCC plus no witness",
        "requires a finite marked-species theorem, not a scalar inequality",
    ),
)


def print_bridge_routes() -> None:
    print()
    print("Creative bridge routes")
    print("=" * 112)
    print(f"{'route':<20} {'A000568 side':<42} {'LRC side':<48} risk")
    print("-" * 112)
    for route, a_side, lrc_side, risk in BRIDGE_ROUTES:
        print(f"{route:<20} {a_side:<42} {lrc_side:<48} {risk}")


def print_synthesis() -> None:
    print()
    print("S507 synthesis")
    print("=" * 112)
    print(
        "The strong version of the user's hypothesis should be phrased as a "
        "marked quotient theorem, not as plain equality with A000568."
    )
    print()
    print("Candidate theorem shape:")
    print()
    print(
        "  LRC counterexamples correspond to nonempty bad orbits in a rooted, "
        "endpoint-marked tournament species over the A000568 base."
    )
    print()
    print("The bad orbit must simultaneously have:")
    print()
    print("  1. no positive lonely gap and no unprotected endpoint;")
    print("  2. high or unstable phase_half H rather than a harmless plateau;")
    print("  3. a safe_phase endpoint alarm;")
    print("  4. a pressure_k2 SCC carrying the surviving endpoint core.")
    print()
    print(
        "The computations here support the opposite pattern on the known hard "
        "ladders: phase shadows stabilize after the first gate while endpoint "
        "debt is exported inside the rooted/marked fiber."
    )


def main() -> None:
    print("LRC/A000568 isoclass analogy atlas (codex-2026-06-01 S507)")
    print("=" * 112)
    print(
        "Question: can the LRC proof be organized as an isomorphism-class "
        "problem analogous to A000568?"
    )
    print()
    print_burnside_table()
    censuses = tuple(class_census(n) for n in range(3, 7))
    print_small_census(censuses)
    rooted_lookup = {row.n: row.rooted for row in censuses}
    clock_rows = tuple(clock_quotient(n, rooted_lookup.get(n)) for n in range(3, 8))
    print_clock_table(clock_rows)
    print_hard_row_shadows()
    print_bridge_routes()
    print_synthesis()


if __name__ == "__main__":
    main()
