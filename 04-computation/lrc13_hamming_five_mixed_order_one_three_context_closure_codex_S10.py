#!/usr/bin/env python3
"""Exact mixed order-one/order-three Hamming-five closure for THM-847.

THM-823 leaves a bounded common-sheet mixed branch consisting of one
effective-order-one replacement and an order-three quartet on a multiplicative
coset of {1,5,8,12}.  After dividing by the common gcd, the retained core is
3P.  The quartet replacements lie in labelled CRT progressions modulo 39,
while the proper order-one replacement at label b is

    3(b+13h),  h >= 1,

so its first allowed speed is 3b+39 and its step is also 39.

This verifier enumerates all 3*8*4=96 coset/flag/parity contexts and applies
the THM-815 longest-component cap at every labelled prefix.  Endpoints and
all decisions use fractions.Fraction.  There is no floating point, sampled
time grid, optimizer, lift-height cutoff, or heuristic pruning.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
from math import floor
from pathlib import Path
from typing import Iterable, Sequence


P = 13
HALF = F(1, 2)
LABELS = tuple(range(1, P))
QUARTET_BASE = (1, 5, 8, 12)
EXPECTED_CERTIFICATE = "1cab41e8d32b09d93e9548d1baa486c0e33fb7979695d1b10725829c3a4aeb75"

RunWord = tuple[F, ...]
Context = tuple[tuple[int, ...], int, tuple[int, int]]


def intersect_runs(left: RunWord, right: RunWord) -> RunWord:
    """Intersect sorted unions of open intervals in flat endpoint form."""

    out: list[F] = []
    i = j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i], right[j])
        hi = min(left[i + 1], right[j + 1])
        if lo < hi:
            out.extend((lo, hi))
        if left[i + 1] <= right[j + 1]:
            i += 2
        else:
            j += 2
    return tuple(out)


@lru_cache(maxsize=None)
def half_safe_runs(speed: int) -> RunWord:
    """Return {0<t<1/2: ||speed*t||>1/13} exactly."""

    bands: list[F] = []
    for k in range(speed):
        lo = max(F(P * k + 1, P * speed), F(0))
        hi = min(F(P * (k + 1) - 1, P * speed), HALF)
        if lo < hi:
            bands.extend((lo, hi))
    return tuple(bands)


def strict_half_core(speeds: Iterable[int]) -> RunWord:
    current: RunWord = (F(0), HALF)
    for speed in sorted(speeds):
        current = intersect_runs(current, half_safe_runs(speed))
        if not current:
            break
    return current


def component_count(runs: RunWord) -> int:
    assert len(runs) % 2 == 0
    return len(runs) // 2


def measure(runs: RunWord) -> F:
    return sum((runs[i + 1] - runs[i] for i in range(0, len(runs), 2)), F(0))


def longest_component(runs: RunWord) -> F:
    assert runs
    return max(runs[i + 1] - runs[i] for i in range(0, len(runs), 2))


def global_kl_bound(runs: RunWord, remaining_count: int) -> int:
    """The older total-residual K/L cap, retained as a comparison."""

    return floor(
        F(
            22 * remaining_count * component_count(runs),
            P * (P - 2 * remaining_count),
        )
        / measure(runs)
    )


def longest_component_bound(runs: RunWord, remaining_count: int) -> int:
    """THM-815's one-component next-speed cap."""

    return floor(
        F(22 * remaining_count, P * (P - 2 * remaining_count))
        / longest_component(runs)
    )


def opposite_pairs(labels: Sequence[int]) -> tuple[tuple[int, int], ...]:
    unused = set(labels)
    answer: list[tuple[int, int]] = []
    while unused:
        r = min(unused)
        pair = tuple(sorted((r, -r % P)))
        assert len(pair) == 2 and set(pair) <= unused
        answer.append(pair)
        unused -= set(pair)
    return tuple(answer)


def quartet_crt_base(label: int, unit: int) -> int:
    """Least positive u with u=3*label (mod 13), u=unit (mod 3)."""

    values = [
        speed
        for speed in range(1, 3 * P)
        if speed % P == 3 * label % P and speed % 3 == unit
    ]
    assert len(values) == 1
    return values[0]


def all_contexts() -> tuple[Context, ...]:
    cosets = sorted(
        {
            tuple(sorted((a * r) % P for r in QUARTET_BASE))
            for a in range(1, P)
        }
    )
    assert len(cosets) == 3
    contexts: list[Context] = []
    for quartet in cosets:
        flags = sorted(set(LABELS) - set(quartet))
        assert len(flags) == 8
        for flag in flags:
            for pair_bits in product((1, 2), repeat=2):
                contexts.append((quartet, flag, pair_bits))
    assert len(contexts) == 3 * 8 * 4 == 96
    return tuple(contexts)


def context_data(
    context: Context,
) -> tuple[tuple[int, ...], dict[int, int], tuple[int, ...]]:
    quartet, flag, pair_bits = context
    pairs = opposite_pairs(quartet)
    assert len(pairs) == 2
    units = {
        label: pair_bits[index]
        for index, pair in enumerate(pairs)
        for label in pair
    }
    missing = tuple(sorted((*quartet, flag)))
    bases = {label: quartet_crt_base(label, units[label]) for label in quartet}
    # The h=0 value 3*flag is the unchanged Hamming-four face.  Proper
    # Hamming-five replacement starts at h=1.
    bases[flag] = 3 * flag + 39
    core_labels = tuple(label for label in LABELS if label not in missing)
    assert len(missing) == 5 and len(core_labels) == 7
    assert bases[flag] % 39 == 3 * flag % 39 and bases[flag] > 3 * flag
    return missing, bases, core_labels


class Audit:
    def __init__(self) -> None:
        self.nodes: Counter[int] = Counter()
        self.dead_no_candidate: Counter[int] = Counter()
        self.covers: list[tuple[tuple[int, int], ...]] = []
        self.terminals: list[tuple[tuple[int, int], ...]] = []
        self.max_long_bound = 0
        self.max_global_bound = 0
        self.operator_speeds: set[int] = set()
        self.dominance_rows = 0
        self.strict_dominance_rows = 0
        self.digest = sha256()

    def record(self, row: object) -> None:
        self.digest.update((repr(row) + "\n").encode())


def recurse(
    *,
    runs: RunWord,
    remaining: tuple[int, ...],
    bases: dict[int, int],
    last_speed: int,
    chosen: tuple[tuple[int, int], ...],
    audit: Audit,
) -> None:
    depth = len(chosen)
    audit.nodes[depth] += 1

    if not runs:
        audit.covers.append(chosen)
        audit.record((depth, chosen, remaining, last_speed, (), F(0), "cover"))
        return

    if not remaining:
        audit.terminals.append(chosen)
        audit.record(
            (
                depth,
                chosen,
                remaining,
                last_speed,
                runs,
                measure(runs),
                "terminal-nonempty",
            )
        )
        return

    count = len(remaining)
    long_bound = longest_component_bound(runs, count)
    global_bound = global_kl_bound(runs, count)
    assert long_bound <= global_bound
    audit.dominance_rows += 1
    audit.strict_dominance_rows += long_bound < global_bound
    audit.max_long_bound = max(audit.max_long_bound, long_bound)
    audit.max_global_bound = max(audit.max_global_bound, global_bound)

    candidates: list[tuple[int, int]] = []
    for label in remaining:
        speed = bases[label]
        if speed <= last_speed:
            speed += 39 * ((last_speed - speed) // 39 + 1)
        while speed <= long_bound:
            candidates.append((speed, label))
            speed += 39
    candidates.sort()
    if not candidates:
        audit.dead_no_candidate[depth] += 1

    audit.record(
        (
            depth,
            chosen,
            remaining,
            last_speed,
            runs,
            measure(runs),
            longest_component(runs),
            long_bound,
            global_bound,
            tuple(candidates),
        )
    )

    for speed, label in candidates:
        audit.operator_speeds.add(speed)
        recurse(
            runs=intersect_runs(runs, half_safe_runs(speed)),
            remaining=tuple(r for r in remaining if r != label),
            bases=bases,
            last_speed=speed,
            chosen=(*chosen, (label, speed)),
            audit=audit,
        )


def comparison_edge(
    values: Sequence[F], tie_keys: Sequence[tuple[int, int]]
) -> tuple[tuple[bool, ...], ...]:
    """Orient larger values first, ties along the declared Hamiltonian path."""

    n = len(values)
    edge = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        if values[i] > values[j] or (
            values[i] == values[j] and tie_keys[i] < tie_keys[j]
        ):
            edge[i][j] = True
        else:
            edge[j][i] = True
    return tuple(tuple(row) for row in edge)


def tournament_fingerprint(edge: Sequence[Sequence[bool]]) -> tuple[object, ...]:
    n = len(edge)
    scores = tuple(sum(row) for row in edge)
    score_histogram = tuple(sorted(Counter(scores).items()))
    triangles = sum(
        (edge[i][j] and edge[j][k] and edge[k][i])
        or (edge[j][i] and edge[k][j] and edge[i][k])
        for i, j, k in combinations(range(n), 3)
    )
    reach = [[i == j or edge[i][j] for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                reach[i][j] = reach[i][j] or (reach[i][k] and reach[k][j])
    unseen = set(range(n))
    scc: list[int] = []
    while unseen:
        i = min(unseen)
        block = {j for j in unseen if reach[i][j] and reach[j][i]}
        scc.append(len(block))
        unseen -= block
    paths = 1 if triangles == 0 and sorted(scores) == list(range(n)) else 0
    return score_histogram, triangles, tuple(sorted(scc, reverse=True)), paths


def context_tournament(
    root: RunWord, missing: tuple[int, ...], bases: dict[int, int]
) -> tuple[object, ...]:
    root_measure = measure(root)
    root_longest = longest_component(root)
    total_values: list[F] = []
    longest_values: list[F] = []
    tie_keys: list[tuple[int, int]] = []
    for label in missing:
        speed = bases[label]
        child = intersect_runs(root, half_safe_runs(speed))
        assert child
        total_values.append(root_measure - measure(child))
        longest_values.append(root_longest - longest_component(child))
        tie_keys.append((speed, label))
    raw = comparison_edge(total_values, tie_keys)
    switched = comparison_edge(longest_values, tie_keys)
    flips = sum(
        raw[i][j] != switched[i][j]
        for i, j in combinations(range(len(missing)), 2)
    )
    return tournament_fingerprint(raw), tournament_fingerprint(switched), flips


def main() -> None:
    contexts = all_contexts()
    context_rows: list[object] = []
    aggregate_nodes: Counter[int] = Counter()
    aggregate_dead: Counter[int] = Counter()
    aggregate_states = 0
    aggregate_dominance = 0
    aggregate_strict_dominance = 0
    all_operator_speeds: set[int] = set()
    tournament_rows: list[tuple[object, ...]] = []

    for index, context in enumerate(contexts):
        quartet, flag, pair_bits = context
        missing, bases, core_labels = context_data(context)
        core_speeds = tuple(3 * label for label in core_labels)
        assert any(speed % 2 == 0 for speed in core_speeds)
        root = strict_half_core(core_speeds)
        assert root and component_count(root) > 0

        audit = Audit()
        recurse(
            runs=root,
            remaining=missing,
            bases=bases,
            last_speed=0,
            chosen=(),
            audit=audit,
        )
        assert not audit.covers
        assert not audit.terminals

        states = sum(audit.nodes.values())
        aggregate_nodes.update(audit.nodes)
        aggregate_dead.update(audit.dead_no_candidate)
        aggregate_states += states
        aggregate_dominance += audit.dominance_rows
        aggregate_strict_dominance += audit.strict_dominance_rows
        all_operator_speeds |= audit.operator_speeds
        tournament = context_tournament(root, missing, bases)
        tournament_rows.append(tournament)

        row = (
            index,
            quartet,
            flag,
            pair_bits,
            missing,
            tuple(sorted(bases.items())),
            core_labels,
            component_count(root),
            measure(root),
            longest_component(root),
            global_kl_bound(root, 5),
            longest_component_bound(root, 5),
            tuple(sorted(audit.nodes.items())),
            tuple(sorted(audit.dead_no_candidate.items())),
            states,
            audit.max_global_bound,
            audit.max_long_bound,
            len(audit.operator_speeds),
            audit.dominance_rows,
            audit.strict_dominance_rows,
            audit.digest.hexdigest(),
            tournament[2],
        )
        context_rows.append(row)

    assert len(context_rows) == 96
    assert aggregate_states == sum(row[14] for row in context_rows)
    assert all(row[14] > 0 for row in context_rows)

    certificate = sha256(
        "".join(repr(row) + "\n" for row in context_rows).encode()
    ).hexdigest()
    if EXPECTED_CERTIFICATE != "TO_BE_FILLED":
        assert certificate == EXPECTED_CERTIFICATE

    standard_fingerprint = (
        ((0, 1), (1, 1), (2, 1), (3, 1), (4, 1)),
        0,
        (1, 1, 1, 1, 1),
        1,
    )
    assert all(row[0] == standard_fingerprint for row in tournament_rows)
    assert all(row[1] == standard_fingerprint for row in tournament_rows)
    flip_hist = Counter(row[2] for row in tournament_rows)

    lines = [
        "THM-847 MIXED ORDER-ONE/ORDER-THREE HAMMING-FIVE CONTEXT CLOSURE",
        "arithmetic=integer+fractions.Fraction floating_point=none sampled_circle=none",
        "threshold=1/13 progression_modulus=39 reflection_quotient=(0,1/2)",
        "normalization=core_3P+quartet_D3_CRT+proper_D1_speed_3(b+13h),h>=1",
        "zero_height_D1_face=speed_3b_is_THM816_Hamming4_and_is_not_enumerated",
        "contexts=3_quartet_cosets*8_flags_outside_coset*4_quartet_pair_words=96",
        "carrier=active_open_endpoint_runs+remaining_labelled_progressions+last_speed",
        "bound.longest=22m/[13(13-2m)Lmax]",
        "bound.global=22mK/[13(13-2m)L]",
        "dominance=Lmax>=L/K_implies_longest_bound<=global_bound",
        f"aggregate.nodes={dict(sorted(aggregate_nodes.items()))}",
        f"aggregate.dead_no_candidate={dict(sorted(aggregate_dead.items()))}",
        f"aggregate.states={aggregate_states}",
        "aggregate.covering_prefixes=0 aggregate.depth5_terminals=0",
        f"aggregate.dominance_rows={aggregate_dominance}",
        f"aggregate.strict_dominance_rows={aggregate_strict_dominance}",
        f"aggregate.distinct_operator_speeds={len(all_operator_speeds)}",
        f"context.min_states={min(row[14] for row in context_rows)}",
        f"context.max_states={max(row[14] for row in context_rows)}",
        f"context.max_long_bound={max(row[16] for row in context_rows)}",
        f"context.max_global_bound={max(row[15] for row in context_rows)}",
        "CONTEXT_ROWS",
    ]

    for row in context_rows:
        lines.append(
            "context[{:02d}] C={} b={} pair_bits={} bases={} "
            "root=(K={},L={},Lmax={},global={},long={}) nodes={} dead={} "
            "states={} max_bounds=({}, {}) operators={} dominance=({},{}) "
            "state_sha256={} tournament_flips={}".format(
                row[0],
                row[1],
                row[2],
                row[3],
                row[5],
                row[7],
                row[8],
                row[9],
                row[10],
                row[11],
                row[12],
                row[13],
                row[14],
                row[15],
                row[16],
                row[17],
                row[18],
                row[19],
                row[20],
                row[21],
            )
        )

    lines.extend(
        (
            "TOURNAMENT_ANALYSIS",
            "vertices=five_labelled_least_allowed_comb_obligations",
            "pair_observable=root_total_measure_marginal_erosion",
            "switch=root_longest_component_marginal_erosion",
            "tie_hamiltonian_path=increasing_(least_allowed_speed,label)",
            f"raw_fingerprint={standard_fingerprint}",
            f"switched_fingerprint={standard_fingerprint}",
            f"edge_flip_histogram={dict(sorted(flip_hist.items()))}",
            f"edge_flips_total={sum(flips * count for flips, count in flip_hist.items())}",
            "preserves=planning_order_only",
            "destroys=absolute_residual_geometry,component_incidence,future_progressions,cover_truth",
            "assumption_challenge=comb_or_runner_vertices_are_not_the_state;component_obligation_incidence_is_faithful",
            "scope=all_THM823_mixed_D1_D3_contexts_not_unbounded_other_orders_not_global_H5",
            f"certificate_sha256={certificate}",
            f"source_sha256={sha256(Path(__file__).read_bytes()).hexdigest()}",
            "PASS: all 96 arbitrary-height mixed order-one/order-three contexts are strictly loose",
        )
    )
    payload = "\n".join(lines) + "\n"
    print(payload, end="")
    print(f"sha256={sha256(payload.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
