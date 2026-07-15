#!/usr/bin/env python3
"""Exact one-context five-comb closure benchmark for the THM-823 frontier.

This replay studies exactly one of the 96 directed-flag/parity contexts left
by THM-823:

    C = {1,5,8,12},  b = 10,  R = C union {b},
    parity bits = (1,1,1),
    CRT bases = {1:16, 5:28, 8:37, 10:4, 12:10} modulo 39.

It is deliberately not a uniform five-comb theorem.  The other 95 contexts
remain open.

It also does not recompute THM-820's scale-one collar/two-box caps.  Its input
is the distinct all-order-three directed-flag context classified by THM-823;
the exact root geometry and the five-comb discrepancy recursion are rebuilt
here from first principles.

The state is the current strict-safe interval union, the remaining labelled
CRT progressions, and the last chosen speed.  Every set involved is invariant
under t -> 1-t.  Since the seven-speed core contains an even speed, t=1/2 is
not safe, so components occur in disjoint mirror pairs.  We therefore retain
only the active endpoint runs in (0,1/2).  Both component count and measure are
halved, leaving K/L and hence every discrepancy bound unchanged.

This active-run representation is important.  A bare mask on a global
endpoint arrangement is not enough to recover K: endpoints belonging only to
unused combs would create artificial component cuts.  A cell implementation
must retain boundary-activation bits, or equivalently the active endpoint runs
used here.

All endpoints, measures, bounds, and tournament observables are exact
Fractions.  There is no floating point and no sampled time grid.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
from math import floor
from pathlib import Path
from typing import Iterable, Sequence


P = 13
HALF = F(1, 2)
BASE = tuple(range(1, P))
QUARTET = (1, 5, 8, 12)
FLAG = 10
MISSING = (1, 5, 8, 10, 12)
PAIR_BITS = (1, 1, 1)
EXPECTED_CERTIFICATE = "9d1024dc22acdc85ccd57afe74ac23d002765ec6ce0aee0d648da64dc23d041a"

RunWord = tuple[F, ...]


def intersect_runs(left: RunWord, right: RunWord) -> RunWord:
    """Intersect two sorted unions of open intervals in flat endpoint form."""

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
    """Return ``{0<t<1/2: ||speed*t||>1/13}`` as active runs."""

    bands: list[F] = []
    for k in range(speed):
        lo = F(P * k + 1, P * speed)
        hi = F(P * (k + 1) - 1, P * speed)
        lo = max(lo, F(0))
        hi = min(hi, HALF)
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


def recursive_bound(runs: RunWord, length: F, remaining_count: int) -> int:
    """Least-speed bound for ``remaining_count`` danger combs.

    The full-circle formula is

        floor(22*m*K / (13*(13-2*m)*L)).

    Mirror restriction halves K and L, so the same expression can be applied
    directly to the half-circle values.
    """

    assert runs and length > 0 and 1 <= remaining_count <= 5
    value = F(
        22 * remaining_count * component_count(runs),
        13 * (13 - 2 * remaining_count),
    ) / length
    return floor(value)


def negative_pairs(labels: Sequence[int]) -> tuple[tuple[int, int], ...]:
    unused = set(labels)
    answer: list[tuple[int, int]] = []
    while unused:
        r = min(unused)
        pair = tuple(sorted((r, -r % P)))
        assert len(pair) == 2 and set(pair) <= unused
        answer.append(pair)
        unused -= set(pair)
    return tuple(answer)


def crt_base(label: int, unit: int) -> int:
    values = [
        speed
        for speed in range(1, 3 * P)
        if speed % P == 3 * label % P and speed % 3 == unit
    ]
    assert len(values) == 1
    return values[0]


def context_bases() -> dict[int, int]:
    pairs = negative_pairs(QUARTET)
    parity = {
        label: PAIR_BITS[index]
        for index, pair in enumerate(pairs)
        for label in pair
    }
    parity[FLAG] = PAIR_BITS[2]
    answer = {label: crt_base(label, parity[label]) for label in MISSING}
    assert answer == {1: 16, 5: 28, 8: 37, 10: 4, 12: 10}
    return answer


def comparison_edge(values: Sequence[F]) -> tuple[tuple[bool, ...], ...]:
    """Orient larger erosion first; index order is the tie Hamiltonian path."""

    n = len(values)
    edge = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        if values[i] >= values[j]:
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
    paths = sum(
        all(edge[path[i]][path[i + 1]] for i in range(n - 1))
        for path in permutations(range(n))
    )
    return score_histogram, triangles, tuple(sorted(scc, reverse=True)), paths


def terminal_tournament(root: RunWord, chosen: tuple[tuple[int, int], ...]) -> tuple[object, ...]:
    """Raw and conditionally switched marginal-erosion tournaments."""

    speeds = tuple(speed for _, speed in chosen)
    assert speeds == tuple(sorted(speeds)) and len(speeds) == 5
    root_length = measure(root)
    raw_values = tuple(
        root_length - measure(intersect_runs(root, half_safe_runs(speed)))
        for speed in speeds
    )
    raw = comparison_edge(raw_values)

    switched = [[False] * 5 for _ in range(5)]
    for i, j in combinations(range(5), 2):
        conditioned = root
        for k, speed in enumerate(speeds):
            if k not in (i, j):
                conditioned = intersect_runs(conditioned, half_safe_runs(speed))
        conditioned_length = measure(conditioned)
        value_i = conditioned_length - measure(
            intersect_runs(conditioned, half_safe_runs(speeds[i]))
        )
        value_j = conditioned_length - measure(
            intersect_runs(conditioned, half_safe_runs(speeds[j]))
        )
        if value_i >= value_j:
            switched[i][j] = True
        else:
            switched[j][i] = True
    switched_edge = tuple(tuple(row) for row in switched)
    flips = sum(
        raw[i][j] != switched_edge[i][j]
        for i, j in combinations(range(5), 2)
    )
    return tournament_fingerprint(raw), tournament_fingerprint(switched_edge), flips


class Audit:
    def __init__(self) -> None:
        self.nodes = Counter()
        self.dead_no_candidate = Counter()
        self.covers: list[tuple[tuple[int, int], ...]] = []
        self.terminals: list[tuple[tuple[tuple[int, int], ...], RunWord, F]] = []
        self.max_bound = 0
        self.operator_speeds: set[int] = set()
        self.digest = sha256()

    def record(
        self,
        *,
        depth: int,
        chosen: tuple[tuple[int, int], ...],
        remaining: tuple[int, ...],
        last_speed: int,
        runs: RunWord,
        length: F,
        status: str,
        bound: int | None,
        candidate_count: int,
    ) -> None:
        row = (
            depth,
            chosen,
            remaining,
            last_speed,
            component_count(runs),
            length,
            status,
            bound,
            candidate_count,
        )
        self.digest.update((repr(row) + "\n").encode())


def recurse(
    *,
    runs: RunWord,
    length: F,
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
        audit.record(
            depth=depth,
            chosen=chosen,
            remaining=remaining,
            last_speed=last_speed,
            runs=runs,
            length=length,
            status="cover",
            bound=None,
            candidate_count=0,
        )
        return
    if not remaining:
        audit.terminals.append((chosen, runs, length))
        audit.record(
            depth=depth,
            chosen=chosen,
            remaining=remaining,
            last_speed=last_speed,
            runs=runs,
            length=length,
            status="terminal-nonempty",
            bound=None,
            candidate_count=0,
        )
        return

    bound = recursive_bound(runs, length, len(remaining))
    audit.max_bound = max(audit.max_bound, bound)
    candidates: list[tuple[int, int]] = []
    for label in remaining:
        speed = bases[label]
        if speed <= last_speed:
            speed += 39 * ((last_speed - speed) // 39 + 1)
        while speed <= bound:
            candidates.append((speed, label))
            speed += 39
    candidates.sort()
    if not candidates:
        audit.dead_no_candidate[depth] += 1
    audit.record(
        depth=depth,
        chosen=chosen,
        remaining=remaining,
        last_speed=last_speed,
        runs=runs,
        length=length,
        status="branch" if candidates else "dead-no-candidate",
        bound=bound,
        candidate_count=len(candidates),
    )

    for speed, label in candidates:
        audit.operator_speeds.add(speed)
        child = intersect_runs(runs, half_safe_runs(speed))
        recurse(
            runs=child,
            length=measure(child),
            remaining=tuple(r for r in remaining if r != label),
            bases=bases,
            last_speed=speed,
            chosen=(*chosen, (label, speed)),
            audit=audit,
        )


def main() -> None:
    bases = context_bases()
    core_labels = tuple(label for label in BASE if label not in MISSING)
    core_speeds = tuple(3 * label for label in core_labels)
    root = strict_half_core(core_speeds)
    root_length = measure(root)
    assert component_count(root) == 18
    assert root_length == F(2615, 18018)
    assert recursive_bound(root, root_length, 5) == 349

    audit = Audit()
    recurse(
        runs=root,
        length=root_length,
        remaining=MISSING,
        bases=bases,
        last_speed=0,
        chosen=(),
        audit=audit,
    )

    assert audit.nodes == Counter({0: 1, 1: 45, 2: 1262, 3: 20703, 4: 53303, 5: 57})
    assert sum(audit.nodes.values()) == 75_371
    assert audit.dead_no_candidate == Counter({3: 6_253, 4: 53_251})
    assert not audit.covers
    assert len(audit.terminals) == 57
    assert audit.max_bound == 1_025
    # The full safe-band cache also holds the seven core operators.  The
    # recursive replacement search itself reaches exactly 130 distinct speeds.
    assert len(audit.operator_speeds) == 130

    terminal_min = min(audit.terminals, key=lambda row: row[2])
    terminal_min_full_measure = 2 * terminal_min[2]
    terminal_min_full_components = 2 * component_count(terminal_min[1])
    assert terminal_min_full_measure == F(413_437, 5_714_280)
    assert terminal_min_full_components == 30
    assert terminal_min[0] == (
        (10, 4),
        (12, 10),
        (1, 16),
        (5, 28),
        (8, 37),
    )

    certificate = audit.digest.hexdigest()
    if EXPECTED_CERTIFICATE != "PENDING":
        assert certificate == EXPECTED_CERTIFICATE

    tournament_rows = [
        terminal_tournament(root, chosen)
        for chosen, _, _ in audit.terminals
    ]
    raw_hist = Counter(raw for raw, _, _ in tournament_rows)
    conditional_hist = Counter(switched for _, switched, _ in tournament_rows)
    flip_hist = Counter(flips for _, _, flips in tournament_rows)

    source_sha = sha256(Path(__file__).read_bytes()).hexdigest()
    search_cache = half_safe_runs.cache_info()
    lines = [
        "hamming-five single-context exact comb closure",
        "scope=one of 96 directed-flag/parity contexts; other_contexts_open=95",
        f"quartet={QUARTET} flag={FLAG} missing={MISSING} pair_bits={PAIR_BITS}",
        f"bases={dict(sorted(bases.items()))} progression_modulus=39",
        f"core_labels={core_labels} core_speeds={core_speeds}",
        "threshold=1/13 duty=2/13 discrepancy=22/169",
        "carrier=active open endpoint runs on (0,1/2) plus remaining labelled CRT progressions and last speed",
        "mirror=t->1-t; full_K=2*half_K; full_L=2*half_L; K/L and every recursive bound unchanged",
        "mask_guardrail=global cell bits without boundary activation do not determine K; unused endpoints create false cuts",
        f"root.half_components={component_count(root)} root.half_measure={root_length}",
        f"root.full_components={2 * component_count(root)} root.full_measure={2 * root_length}",
        f"root.first_speed_bound={recursive_bound(root, root_length, 5)}",
        f"nodes_by_depth={dict(sorted(audit.nodes.items()))}",
        f"total_states={sum(audit.nodes.values())}",
        f"dead_no_candidate_by_depth={dict(sorted(audit.dead_no_candidate.items()))}",
        f"covering_prefixes={len(audit.covers)}",
        f"terminal_nonempty_rows={len(audit.terminals)}",
        f"terminal_min_full_measure={terminal_min_full_measure}",
        f"terminal_min_full_components={terminal_min_full_components}",
        f"terminal_min_chosen_label_speed={terminal_min[0]}",
        f"max_recursive_bound={audit.max_bound}",
        f"distinct_cached_replacement_operators={len(audit.operator_speeds)}",
        f"cached_half_safe_runs={search_cache}",
        f"state_certificate_sha256={certificate}",
        "tournament.vertices=five chosen exception combs ordered by increasing speed",
        "tournament.raw_observable=marginal erosion of the root half-carrier",
        "tournament.switch=remove the other three combs before each pair comparison",
        "tournament.tie_hamiltonian_path=increasing speed",
        f"tournament.raw_fingerprint_histogram={dict(sorted(raw_hist.items(), key=repr))}",
        f"tournament.conditional_fingerprint_histogram={dict(sorted(conditional_hist.items(), key=repr))}",
        f"tournament.edge_flip_histogram={dict(sorted(flip_hist.items()))}",
        "assumption_challenge=comb vertices are overlap telemetry only; they forget absolute residual measure, active endpoint topology, and future CRT progression action",
        "closure_claim=this context only; no all-context or arbitrary-Hamming-five closure is claimed",
        "thm820_interface=does not recompute THM-820 scale-one bounds; independently rebuilds the THM-823 all-order-three context root and its five-comb bounds",
        f"source_sha256={source_sha}",
    ]
    output_hash = sha256(("\n".join(lines) + "\n").encode()).hexdigest()
    lines.append(f"sha256={output_hash}")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
