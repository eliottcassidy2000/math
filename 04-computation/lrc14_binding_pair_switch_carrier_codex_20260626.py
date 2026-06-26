#!/usr/bin/env python3
"""Binding-pair switch carrier for LRC14.

This is a proof-angle scout, not a proof.  THM-524 says an LRC gap maximum is
attained at a pair crossing.  The dangerous shortcut is to keep only the
pair-crossing value and forget the other-runner clearance scan.  This script
keeps the switch data as a labelled packet coordinate:

  binding pair, sum/difference lane, active blockers, other-runner clearance,
  grid/off-grid status, decoy pair switches, and strict safe measure.

Tournament Analysis uses proof/switch carriers as vertices, not runners.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from math import gcd


N = 14
THRESHOLD = F(1, N)
HALF = F(1, 2)


def fmt(x: F) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def primitive(row: tuple[int, ...]) -> bool:
    g = 0
    for v in row:
        g = gcd(g, abs(v))
    return g == 1


def dist(v: int, t: F) -> F:
    r = (v * t) % 1
    return min(r, 1 - r)


def is_int(x: F) -> bool:
    return x.denominator == 1


def candidate_times(row: tuple[int, ...]) -> tuple[F, ...]:
    """THM-524 pair-crossing candidates in [0, 1/2], plus the half-turn."""
    times: set[F] = {HALF}
    for a, b in combinations(row, 2):
        for d in {a + b, abs(a - b)}:
            if d <= 0:
                continue
            for k in range(1, d):
                t = F(k, d)
                if t <= HALF:
                    times.add(t)
    return tuple(sorted(times))


def strict_safe_measure(row: tuple[int, ...], threshold: F = THRESHOLD) -> F:
    """Measure of times with every runner strictly outside the threshold tube."""
    intervals: list[tuple[F, F]] = []
    for v in row:
        for m in range(v):
            center = F(m, v)
            lo = (center - threshold / v) % 1
            hi = (center + threshold / v) % 1
            if lo < hi:
                intervals.append((lo, hi))
            else:
                intervals.append((F(0), hi))
                intervals.append((lo, F(1)))
    if not intervals:
        return F(1)
    intervals.sort()
    total = F(0)
    a0, b0 = intervals[0]
    for a, b in intervals[1:]:
        if a <= b0:
            b0 = max(b0, b)
        else:
            total += b0 - a0
            a0, b0 = a, b
    total += b0 - a0
    return F(1) - total


@dataclass(frozen=True)
class TimeDigest:
    t: F
    gap: F
    active: tuple[int, ...]
    tied_pair_gap: F
    tied_pair_count: int
    active_pair_labels: tuple[str, ...]


@dataclass(frozen=True)
class RowDigest:
    name: str
    row: tuple[int, ...]
    family: str
    candidate_count: int
    M: F
    strict_safe_mu: F
    grid_witnesses: tuple[F, ...]
    witness_count: int
    decoy_count: int
    best_decoy_pair_gap: F
    best_decoy_actual_gap: F
    optimums: tuple[TimeDigest, ...]


def active_pair_label(a: int, b: int, t: F) -> str:
    modes: list[str] = []
    if is_int((a + b) * t):
        modes.append("sum")
    if is_int(abs(a - b) * t):
        modes.append("diff")
    if not modes:
        modes.append("same-distance")
    mod14 = (a + b) % N
    flags: list[str] = []
    if mod14 == 0:
        flags.append("sum0mod14")
    if a + b == N:
        flags.append("literal14")
    if gcd(a + b, N) > 1:
        flags.append(f"nonunit-sum-gcd{gcd(a + b, N)}")
    if abs(a - b) and gcd(abs(a - b), N) > 1:
        flags.append(f"nonunit-diff-gcd{gcd(abs(a - b), N)}")
    suffix = ",".join(flags) if flags else f"sum_mod14={mod14}"
    return f"{a},{b}:{'+'.join(modes)}:{suffix}"


def digest_time(row: tuple[int, ...], t: F) -> TimeDigest:
    distances = {v: dist(v, t) for v in row}
    gap = min(distances.values())
    active = tuple(v for v in row if distances[v] == gap)
    tied_pair_gap = F(0)
    tied_pair_count = 0
    for a, b in combinations(row, 2):
        if distances[a] == distances[b]:
            tied_pair_count += 1
            tied_pair_gap = max(tied_pair_gap, distances[a])
    labels = tuple(
        active_pair_label(a, b, t)
        for a, b in combinations(active, 2)
        if distances[a] == distances[b]
    )
    return TimeDigest(
        t=t,
        gap=gap,
        active=active,
        tied_pair_gap=tied_pair_gap,
        tied_pair_count=tied_pair_count,
        active_pair_labels=labels,
    )


def digest_row(name: str, row: tuple[int, ...], family: str) -> RowDigest:
    assert len(row) == 13 and primitive(row), row
    digests = tuple(digest_time(row, t) for t in candidate_times(row))
    M = max(d.gap for d in digests)
    optimums = tuple(d for d in digests if d.gap == M)
    grid_witnesses = tuple(
        t for t in (F(a, N) for a in range(1, N)) if min(dist(v, t) for v in row) >= THRESHOLD
    )
    witness_count = sum(1 for d in digests if d.gap >= THRESHOLD)
    decoys = tuple(d for d in digests if d.tied_pair_gap >= THRESHOLD and d.gap < THRESHOLD)
    if decoys:
        best_decoy = max(decoys, key=lambda d: (d.tied_pair_gap, d.gap, -d.t.denominator))
        best_pair_gap = best_decoy.tied_pair_gap
        best_actual = best_decoy.gap
    else:
        best_pair_gap = F(0)
        best_actual = F(0)
    return RowDigest(
        name=name,
        row=row,
        family=family,
        candidate_count=len(digests),
        M=M,
        strict_safe_mu=strict_safe_measure(row),
        grid_witnesses=grid_witnesses,
        witness_count=witness_count,
        decoy_count=len(decoys),
        best_decoy_pair_gap=best_pair_gap,
        best_decoy_actual_gap=best_actual,
        optimums=optimums,
    )


def named_rows() -> tuple[tuple[str, tuple[int, ...], str], ...]:
    ap = tuple(range(1, 14))
    return (
        ("AP13", ap, "AP/GW boundary"),
        ("GW 12->24", tuple(list(range(1, 12)) + [13, 24]), "AP/GW boundary"),
        ("petal 10->20", tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 20]), "unit petal"),
        ("petal 13->26", tuple(list(range(1, 13)) + [26]), "unit petal"),
        ("K33 12->36", tuple(list(range(1, 12)) + [13, 36]), "K33/state-lift"),
        ("K33 family drop12,13 add26,36", tuple(list(range(1, 12)) + [26, 36]), "K33/state-lift"),
        ("covering 12->84", tuple(list(range(1, 12)) + [13, 84]), "covering"),
        ("single far tail 12->200", tuple(list(range(1, 13)) + [200]), "covering"),
        ("drop6 core add180", tuple([1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 180]), "drop-core"),
        ("two-far doublet 20,21", tuple(list(range(1, 12)) + [20, 21]), "genuine-wide doublet"),
    )


@dataclass(frozen=True)
class Carrier:
    name: str
    score: tuple[int, ...]
    preserves: str
    destroys: str


def switch_carriers() -> tuple[Carrier, ...]:
    return (
        Carrier(
            "others_clear_binding_switch",
            (6, 6, 6, 6, 5, 6, 6),
            "exact M witness, pair lane, and clearance",
            "nothing material if endpoint owners are retained",
        ),
        Carrier(
            "active_blocker_deck",
            (5, 6, 5, 4, 5, 6, 5),
            "which runners actually attain the lower envelope",
            "non-active near misses",
        ),
        Carrier(
            "binding_denominator_lane",
            (5, 6, 6, 5, 4, 5, 6),
            "sum/difference denominator and Farey scale",
            "owner labels unless attached",
        ),
        Carrier(
            "grid_section_unit_witness",
            (6, 5, 5, 2, 5, 6, 6),
            "q=14 section witness and unit action",
            "off-grid covering witnesses",
        ),
        Carrier(
            "complement_sum14_boundary",
            (5, 6, 5, 2, 6, 5, 5),
            "AP/GW literal complement boundary pair",
            "positive open off-grid routes",
        ),
        Carrier(
            "covering_big_flank_switch",
            (6, 6, 6, 6, 3, 5, 5),
            "off-grid big-runner crossing and clearance",
            "small endpoint-owner cycle unless recorded",
        ),
        Carrier(
            "state_lift_nonunit_pair",
            (5, 5, 5, 5, 6, 4, 5),
            "K33/nonunit incidence side channel",
            "direct safe interval geometry",
        ),
        Carrier(
            "raw_pair_gap_shadow",
            (2, 1, 3, 2, 1, 0, 3),
            "pair gap before clearance",
            "the LRC predicate; produces many decoys",
        ),
        Carrier(
            "snapshot_position_tournament",
            (1, 1, 1, 1, 1, 0, 2),
            "total order of one time slice",
            "binding-pair identity and clearance",
        ),
    )


def tournament_fingerprint(carriers: tuple[Carrier, ...]) -> dict[str, object]:
    n = len(carriers)
    out = [0] * n
    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        si, sj = carriers[i].score, carriers[j].score
        wi = sum(a > b for a, b in zip(si, sj))
        wj = sum(b > a for a, b in zip(si, sj))
        if wi > wj or (wi == wj and i < j):
            a, b = i, j
        else:
            a, b = j, i
        adj[a][b] = True
        out[a] += 1

    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if (adj[i][j] and adj[j][k] and adj[k][i]) or (adj[i][k] and adj[k][j] and adj[j][i]):
            c3 += 1

    # Kosaraju SCC sizes.
    def dfs(v: int, graph: list[list[bool]], seen: set[int], order: list[int]) -> None:
        seen.add(v)
        for w, ok in enumerate(graph[v]):
            if ok and w not in seen:
                dfs(w, graph, seen, order)
        order.append(v)

    seen: set[int] = set()
    order: list[int] = []
    for v in range(n):
        if v not in seen:
            dfs(v, adj, seen, order)
    radj = [[adj[j][i] for j in range(n)] for i in range(n)]
    seen.clear()
    scc_sizes: list[int] = []
    for v in reversed(order):
        if v in seen:
            continue
        comp: list[int] = []
        dfs(v, radj, seen, comp)
        scc_sizes.append(len(comp))

    # Hamiltonian path count by subset DP.
    dp: dict[tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(1 << n):
        for last in range(n):
            ways = dp.get((mask, last), 0)
            if not ways:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + ways
    full = (1 << n) - 1
    hp = sum(dp.get((full, last), 0) for last in range(n))
    path = sorted(range(n), key=lambda i: carriers[i].score, reverse=True)
    return {
        "score_hist": sorted(Counter(out).items()),
        "directed_3cycles": c3,
        "scc_sizes": sorted(scc_sizes, reverse=True),
        "hamiltonian_path_count": hp,
        "tie_path": " > ".join(carriers[i].name for i in path),
    }


def print_assumption_challenge() -> None:
    print("Assumption challenge")
    print("  Candidate vertex sets considered:")
    print("    runners, gaps, fixed sections, section boundaries, wall crossings,")
    print("    residues, cover arcs, Fourier modes, matroid circuits, proof obligations,")
    print("    and binding-pair switches.")
    print("  Chosen vertex set:")
    print("    binding-pair/proof switches, because THM-524 says the max is there but")
    print("    the LRC predicate is only preserved after the other-runner clearance scan.")
    print("  Preserved predicate:")
    print("    existence of a time with min_v ||v t|| >= 1/14.")
    print("  Destroyed by raw pair quotients:")
    print("    active blocker deck, endpoint owners, off-grid witness intervals, and")
    print("    whether the promising pair is blocked by another runner.")
    print()


def main() -> None:
    print("LRC14 binding-pair switch carrier scout")
    print("=" * 72)
    print_assumption_challenge()

    digests = tuple(digest_row(name, row, family) for name, row, family in named_rows())
    print("[1] Named-row binding switch audit")
    for d in digests:
        grid = ",".join(fmt(t) for t in d.grid_witnesses[:6]) or "none"
        more = "" if len(d.grid_witnesses) <= 6 else f" (+{len(d.grid_witnesses)-6})"
        print(f"  {d.name}")
        print(f"    family={d.family}; row={d.row}")
        print(
            f"    candidates={d.candidate_count}; M={fmt(d.M)}; "
            f"strict_safe_mu={fmt(d.strict_safe_mu)}; grid_witnesses={grid}{more}"
        )
        print(
            f"    witness_switch_times={d.witness_count}; decoy_pair_times={d.decoy_count}; "
            f"best_decoy_pair_gap={fmt(d.best_decoy_pair_gap)} blocked_to={fmt(d.best_decoy_actual_gap)}"
        )
        for opt in d.optimums[:4]:
            labels = "; ".join(opt.active_pair_labels) if opt.active_pair_labels else "no active pair label"
            print(
                f"      optimum t={fmt(opt.t)} gap={fmt(opt.gap)} active={opt.active} "
                f"labels={labels}"
            )
        if len(d.optimums) > 4:
            print(f"      ... {len(d.optimums) - 4} additional optimum times")
    print()

    print("[2] Cross-row readout")
    offgrid = [d.name for d in digests if d.M >= THRESHOLD and not d.grid_witnesses]
    boundary = [d.name for d in digests if d.M == THRESHOLD and d.strict_safe_mu == 0]
    positive = [d.name for d in digests if d.strict_safe_mu > 0]
    decoy_heavy = sorted(digests, key=lambda d: d.decoy_count, reverse=True)[:4]
    print(f"  boundary atoms (M=1/14 and no strict open mass): {boundary}")
    print(f"  positive-open rows in this audit: {positive}")
    print(f"  off-grid-only witness rows: {offgrid}")
    print("  largest raw-pair decoy counts:")
    for d in decoy_heavy:
        print(f"    {d.name}: decoys={d.decoy_count}, best_pair_gap={fmt(d.best_decoy_pair_gap)} -> actual {fmt(d.best_decoy_actual_gap)}")
    print()

    print("[3] Tournament Analysis")
    carriers = switch_carriers()
    fp = tournament_fingerprint(carriers)
    print("  vertices are switch/proof carriers, not runners.")
    print("  pairwise observable:")
    print("    predicate retention, clearance retention, exact denominator scale,")
    print("    off-grid coverage, family separation, decoy resistance, proof maturity.")
    print("  binary gauge:")
    print("    carrier A -> B when A wins a majority of observable coordinates; ties")
    print("    use the declared Hamiltonian path below.")
    print(f"  score_hist={fp['score_hist']}")
    print(f"  directed_3cycles={fp['directed_3cycles']}")
    print(f"  scc_sizes={fp['scc_sizes']}")
    print(f"  hamiltonian_path_count={fp['hamiltonian_path_count']}")
    print(f"  tie_hamiltonian_path={fp['tie_path']}")
    print()

    print("[4] Proof-angle target")
    print("  Add binding-switch fields to HYP-2963-style packet records:")
    print("    binding_switch_type, active_pair_shell, denominator_lane,")
    print("    active_blocker_deck, other_clearance_margin, decoy_pair_gap,")
    print("    grid_or_offgrid_status, and switch_witness_count.")
    print("  The theorem-facing lemma is not 'large pair gap implies LRC witness'.")
    print("  It is the stronger switch statement:")
    print("    a quotient may keep a binding pair only when the active blocker deck and")
    print("    all other-runner clearance inequalities are retained, reconstructed,")
    print("    dual-annihilated, or routed to AP/GW, covering, petal, K33, or F7 debt.")


if __name__ == "__main__":
    main()
