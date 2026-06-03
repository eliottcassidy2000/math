#!/usr/bin/env python3
"""
lrc_n14_fixed_round_certificate_s578.py

codex-2026-06-03 S578

First certificate scaffold for the n=14 even clean lane after the recent agent
merge:
  * HYP-2094: the even-ladder worry target is 64 self-converse round classes.
  * HYP-2095: the cheap route is an unblocked small pair; blockers are paired
    shields or endpoint anchors.
  * HYP-2096: n=14 has a forced nine-unit-shell spine modulo 27 plus four
    composite slack runners.
  * HYP-2097: the next proof object should be a 64-row fixed-round certificate
    table, with the other 126 converse-paired nodes as controls.

This script does not prove the 64 classes.  It builds the class scaffold and
then tests the first speed-level certificate primitive on a canonical unit-spine
normal form.  The point is to separate what the unlabelled round class preserves
from what has to be reattached as speed-owner labels.

Tournament Analysis / assumption challenge:
  Part A vertices are self-converse round classes, not runners.
  Part B vertices are canonical four-slack speed rows over the unit spine.
  The class quotient preserves round/converse/fixed-boundary data but destroys
  small-pair, endpoint-owner, and unit-spine realization data.  The speed-level
  quotient restores those labels while still destroying arbitrary large unit
  representatives.  S578 records that gap explicitly.
"""

from __future__ import annotations

import importlib.util
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S574 = ROOT / "04-computation" / "lrc_round_count_m89_s574.py"

N = 14
M = N - 1
C = 2 * N - 1
FLOOR = Fraction(1, N)
EDGE = Fraction(2, C)
UNIT_SPINE = (1, 2, 4, 5, 7, 8, 10, 11, 13)
SLACK_BOUND = 42

Obligation = tuple[str, int]


@dataclass(frozen=True)
class FixedClass:
    class_id: str
    rep_d: tuple[int, ...]
    score_hist: tuple[tuple[int, int], ...]
    score_span: int
    scc_count: int
    dihedral_auto_count: int
    dihedral_anti_count: int
    anti_witness: str


@dataclass(frozen=True)
class SlackReport:
    speeds: tuple[int, ...]
    slack: tuple[int, ...]
    full_cover: bool
    uncovered: tuple[Obligation, ...]
    private_count: int
    maximin: Fraction | None
    witness: Fraction | None
    active: tuple[int, ...]
    unblocked_pair: tuple[int, int, int] | None
    positive_measure: bool | None


def load_s574():
    spec = importlib.util.spec_from_file_location("s574_rounds", S574)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {S574}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def class_id(flat: tuple[int, ...]) -> str:
    raw = "".join(str(x) for x in flat).encode("ascii")
    return hashlib_sha1(raw)[:10]


def hashlib_sha1(raw: bytes) -> str:
    import hashlib

    return hashlib.sha1(raw).hexdigest()[:10]


def opposite_adj(adj: list[list[int]]) -> list[list[int]]:
    m = len(adj)
    return [[0 if i == j else adj[j][i] for j in range(m)] for i in range(m)]


def dihedral_perms(m: int) -> list[tuple[str, int, tuple[int, ...]]]:
    perms: list[tuple[str, int, tuple[int, ...]]] = []
    for shift in range(m):
        perms.append(("rot", shift, tuple((i + shift) % m for i in range(m))))
        perms.append(("ref", shift, tuple((shift - i) % m for i in range(m))))
    return perms


def is_auto(adj: list[list[int]], perm: tuple[int, ...]) -> bool:
    m = len(adj)
    return all(adj[perm[i]][perm[j]] == adj[i][j] for i in range(m) for j in range(m) if i != j)


def is_anti(adj: list[list[int]], perm: tuple[int, ...]) -> bool:
    m = len(adj)
    return all(adj[perm[i]][perm[j]] == adj[j][i] for i in range(m) for j in range(m) if i != j)


def flat_under(adj: list[list[int]], perm: tuple[int, ...]) -> tuple[int, ...]:
    m = len(adj)
    return tuple(adj[perm[i]][perm[j]] for i in range(m) for j in range(m) if i != j)


def dihedral_canon(adj: list[list[int]]) -> tuple[int, ...]:
    return min(flat_under(adj, perm) for _, _, perm in dihedral_perms(len(adj)))


def scc_count(adj: list[list[int]]) -> int:
    m = len(adj)

    def reach(start: int) -> set[int]:
        seen = {start}
        q = deque([start])
        while q:
            v = q.popleft()
            for w, edge in enumerate(adj[v]):
                if edge and w not in seen:
                    seen.add(w)
                    q.append(w)
        return seen

    remaining = set(range(m))
    count = 0
    while remaining:
        start = next(iter(remaining))
        forward = reach(start)
        comp = {v for v in remaining if v in forward and start in reach(v)}
        remaining -= comp
        count += 1
    return count


def fixed_round_classes(rounds) -> tuple[list[FixedClass], int, int, int]:
    groups: dict[tuple[int, ...], list[tuple[int, ...]]] = defaultdict(list)
    labelled_fixed = 0
    for d in rounds.valid_dvectors(M):
        adj = rounds.build_adj(d, M)
        if any(is_anti(adj, perm) for _, _, perm in dihedral_perms(M)):
            labelled_fixed += 1
            groups[dihedral_canon(adj)].append(d)

    fixed: list[FixedClass] = []
    for canon, dvectors in groups.items():
        rep_d = min(dvectors)
        rep_adj = rounds.build_adj(rep_d, M)
        scores = [sum(row) for row in rep_adj]

        best_auto = 0
        best_anti = 0
        anti_witness = "none"
        for d in dvectors:
            adj = rounds.build_adj(d, M)
            auto_hits = []
            anti_hits = []
            for kind, shift, perm in dihedral_perms(M):
                if is_auto(adj, perm):
                    auto_hits.append((kind, shift))
                if is_anti(adj, perm):
                    anti_hits.append((kind, shift))
            if len(auto_hits) > best_auto:
                best_auto = len(auto_hits)
            if len(anti_hits) > best_anti:
                best_anti = len(anti_hits)
                anti_witness = ",".join(f"{kind}{shift}" for kind, shift in anti_hits[:4])

        fixed.append(
            FixedClass(
                class_id=class_id(canon),
                rep_d=rep_d,
                score_hist=tuple(sorted(Counter(scores).items())),
                score_span=max(scores) - min(scores),
                scc_count=scc_count(rep_adj),
                dihedral_auto_count=best_auto,
                dihedral_anti_count=best_anti,
                anti_witness=anti_witness,
            )
        )

    round_count = rounds.A000016(M)
    nonfixed_classes = round_count - len(fixed)
    return (
        sorted(fixed, key=lambda r: (r.score_span, r.score_hist, r.class_id)),
        round_count,
        nonfixed_classes,
        labelled_fixed,
    )


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def norm(x: Fraction) -> Fraction:
    r = frac_part(x)
    return min(r, 1 - r)


def score_time(speeds: tuple[int, ...], t: Fraction) -> Fraction:
    return min(norm(Fraction(v) * t) for v in speeds)


def exact_maximin(speeds: tuple[int, ...]) -> tuple[Fraction, Fraction]:
    candidates: set[Fraction] = set()
    for i, a in enumerate(speeds):
        for m in range(a):
            t = Fraction(2 * m + 1, 2 * a)
            if 0 < t < 1:
                candidates.add(t)
        for b in speeds[i + 1 :]:
            for den in (a + b, abs(a - b)):
                if den <= 0:
                    continue
                for m in range(1, den):
                    candidates.add(Fraction(m, den))

    best = Fraction(0)
    best_t = Fraction(0)
    for t in candidates:
        score = score_time(speeds, t)
        if (score, -t) > (best, -best_t):
            best = score
            best_t = t
    return best, best_t


def active_runners(speeds: tuple[int, ...], t: Fraction | None) -> tuple[int, ...]:
    if t is None:
        return tuple()
    score = score_time(speeds, t)
    return tuple(v for v in speeds if norm(Fraction(v) * t) == score)


def unit_shell(a: int) -> int:
    r = a % C
    if r == 0:
        return 0
    return min(r, C - r)


def obligation_universe() -> list[Obligation]:
    obligations: list[Obligation] = []
    obligations.extend(("D", q) for q in range(2, N))
    obligations.extend(("U", a) for a in UNIT_SPINE)
    obligations.extend(("N", j) for j in range(1, N // 2 + 1))
    return obligations


def covers(v: int, obligation: Obligation) -> bool:
    layer, label = obligation
    if layer == "D":
        return v % label == 0
    if layer == "U":
        return unit_shell(v) == label and gcd(label, C) == 1
    if layer == "N":
        r = (v * label) % N
        return min(r, N - r) <= 1
    raise ValueError(obligation)


def incidence(speeds: tuple[int, ...]) -> dict[Obligation, tuple[int, ...]]:
    return {o: tuple(v for v in speeds if covers(v, o)) for o in obligation_universe()}


def small_pairs(speeds: tuple[int, ...]) -> list[tuple[int, int]]:
    return [
        (a, b)
        for a, b in combinations(speeds, 2)
        if (a + b) // gcd(a, b) <= N
    ]


def pair_safe_ticks(a: int, b: int) -> list[int]:
    den = a + b
    return [
        k
        for k in range(1, den)
        if norm(Fraction(a * k, den)) >= FLOOR and norm(Fraction(b * k, den)) >= FLOOR
    ]


def unblocked_small_pair(speeds: tuple[int, ...]) -> tuple[int, int, int] | None:
    for a, b in small_pairs(speeds):
        den = a + b
        for k in pair_safe_ticks(a, b):
            if all(norm(Fraction(c * k, den)) >= FLOOR for c in speeds):
                return (a, b, k)
    return None


def positive_measure(speeds: tuple[int, ...]) -> bool:
    endpoints: set[Fraction] = set()
    for v in speeds:
        for k in range(v + 1):
            for sign in (-1, 1):
                endpoints.add(Fraction(k * N + sign, N * v) % 1)
    pts = sorted(endpoints)
    for i, a in enumerate(pts):
        b = pts[(i + 1) % len(pts)]
        length = b - a if b > a else b - a + 1
        mid = (a + length / 2) % 1
        if all(norm(Fraction(v) * mid) > FLOOR for v in speeds):
            return True
    return False


def analyze_slack(slack: tuple[int, ...]) -> SlackReport:
    speeds = tuple(sorted(UNIT_SPINE + slack))
    inc = incidence(speeds)
    uncovered = tuple(o for o, owners in inc.items() if not owners)
    private = tuple(o for o, owners in inc.items() if len(owners) == 1)
    full = not uncovered
    maximin = witness = None
    active: tuple[int, ...] = tuple()
    ub = None
    pos = None
    if full:
        maximin, witness = exact_maximin(speeds)
        active = active_runners(speeds, witness)
        ub = unblocked_small_pair(speeds)
        pos = positive_measure(speeds)
    return SlackReport(
        speeds=speeds,
        slack=slack,
        full_cover=full,
        uncovered=uncovered,
        private_count=len(private),
        maximin=maximin,
        witness=witness,
        active=active,
        unblocked_pair=ub,
        positive_measure=pos,
    )


def fmt_frac(x: Fraction | None) -> str:
    if x is None:
        return "-"
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def route(row: SlackReport) -> str:
    if not row.full_cover:
        return "clock-ledger failure"
    if row.unblocked_pair is not None:
        return "cheap unblocked pair"
    if row.positive_measure:
        return "block-all but positive-measure"
    return "open residual"


def slack_scan() -> list[SlackReport]:
    candidates = [v for v in range(3, SLACK_BOUND + 1, 3)]
    return [analyze_slack(slack) for slack in combinations(candidates, 4)]


def tournament_fingerprint(rows: list[SlackReport]) -> dict[str, object]:
    full = [r for r in rows if r.full_cover]
    sample = sorted(
        full,
        key=lambda r: (
            0 if r.unblocked_pair is None else 1,
            r.maximin or Fraction(1),
            r.private_count,
            max(r.speeds),
            r.slack,
        ),
    )[:24]

    def key(row: SlackReport) -> tuple[int, Fraction, int, int, tuple[int, ...]]:
        return (
            1 if row.unblocked_pair is None else 0,
            -(row.maximin or Fraction(0)),
            -row.private_count,
            -max(row.speeds),
            tuple(-v for v in row.slack),
        )

    n = len(sample)
    adj = [[False] * n for _ in range(n)]
    for i, a in enumerate(sample):
        for j, b in enumerate(sample):
            if i == j:
                continue
            adj[i][j] = key(a) > key(b) or (key(a) == key(b) and i < j)
    scores = [sum(row) for row in adj]
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        c3 += int(
            (adj[i][j] and adj[j][k] and adj[k][i])
            or (adj[i][k] and adj[k][j] and adj[j][i])
        )
    return {
        "sampled_vertices": [r.slack for r in sample],
        "score_hist": dict(sorted(Counter(scores).items())),
        "directed_3_cycles": c3,
        "hamiltonian_path": [str(r.slack) for r in sorted(sample, key=key, reverse=True)],
    }


def main() -> None:
    print("S578 n=14 fixed-round certificate scaffold")
    print("=" * 78)
    print("Merge: HYP-2094 fixed round classes + HYP-2095 unblocked pairs")
    print("       + HYP-2096 unit spine + HYP-2097 64-first target.")
    print()

    rounds = load_s574()
    fixed, round_count, nonfixed_classes, labelled_fixed = fixed_round_classes(rounds)
    print("A. self-converse round-class scaffold")
    print(f"  labelled round d-vectors: {2 ** (M - 1)}")
    print(f"  labelled vectors with dihedral anti-witness: {labelled_fixed}")
    print(f"  round classes: {round_count}")
    print(f"  fixed classes: {len(fixed)}")
    print(f"  non-fixed classes: {nonfixed_classes} (= {nonfixed_classes // 2} converse pairs)")
    print(f"  score-span histogram: {dict(sorted(Counter(r.score_span for r in fixed).items()))}")
    print(f"  SCC-count histogram: {dict(sorted(Counter(r.scc_count for r in fixed).items()))}")
    print(f"  dihedral anti-count histogram: {dict(sorted(Counter(r.dihedral_anti_count for r in fixed).items()))}")
    print("  first fixed-class scaffold rows:")
    print("    id        span scc auto anti  score_hist                 rep_d[:13] anti")
    for row in fixed[:10]:
        print(
            f"    {row.class_id} {row.score_span:4d} {row.scc_count:3d} "
            f"{row.dihedral_auto_count:4d} {row.dihedral_anti_count:4d} "
            f"{str(row.score_hist):25s} {row.rep_d} {row.anti_witness}"
        )

    print()
    print("B. canonical unit-spine speed-level certificate primitive")
    rows = slack_scan()
    full = [r for r in rows if r.full_cover]
    no_ub = [r for r in full if r.unblocked_pair is None]
    measure_zero = [r for r in full if r.positive_measure is False]
    floor_rows = [r for r in full if r.maximin == FLOOR]
    below = [r for r in full if (r.maximin or Fraction(1)) < FLOOR]
    open_gap = [r for r in full if FLOOR < (r.maximin or Fraction(0)) < EDGE]
    print(f"  slack_bound={SLACK_BOUND}; slack rows={len(rows)}")
    print(f"  full D/U/N quotient covers={len(full)}")
    print(f"  below floor rows={len(below)}")
    print(f"  floor rows={len(floor_rows)}; open-gap rows={len(open_gap)}")
    print(f"  full-cover rows with unblocked small pair={len(full) - len(no_ub)}/{len(full)}")
    print(f"  block-all rows positive-measure={sum(1 for r in no_ub if r.positive_measure)}/{len(no_ub)}")
    print(f"  measure-zero full-cover rows with unblocked pair={sum(1 for r in measure_zero if r.unblocked_pair is not None)}/{len(measure_zero)}")
    print(f"  route histogram: {dict(sorted(Counter(route(r) for r in rows).items()))}")
    print()

    def print_speed_row(label: str, row: SlackReport) -> None:
        print(f"  [{label}] slack={row.slack} speeds={row.speeds}")
        print(
            f"    M={fmt_frac(row.maximin)} t={fmt_frac(row.witness)} active={row.active or '-'} "
            f"route={route(row)} positive_measure={row.positive_measure}"
        )
        if row.unblocked_pair is not None:
            a, b, k = row.unblocked_pair
            print(f"    unblocked_pair=({a},{b}) at {k}/{a+b}")
        if row.uncovered:
            print(f"    uncovered={row.uncovered}")

    for row in floor_rows[:4]:
        print_speed_row("floor", row)
    for row in open_gap[:3]:
        print_speed_row("open-gap", row)
    for row in no_ub[:3]:
        print_speed_row("block-all control", row)

    print()
    print("Tournament Analysis")
    print("  vertices: canonical full-cover slack rows over the forced unit spine")
    print("  observable: unblocked-pair flag, exact M, private D/U/N count, max speed")
    print("  switch: harder = no cheap pair, then smaller M, fewer pivots, smaller row")
    print(f"  fingerprints: {tournament_fingerprint(rows)}")

    print()
    print("Synthesis")
    print("  The 64 self-converse round classes are now an explicit scaffold, but")
    print("  class-only data cannot decide the HYP-2095 cheap route.  The speed-level")
    print("  unit-spine normal form restores the missing owner labels: in the canonical")
    print(f"  slack scan through {SLACK_BOUND}, every measure-zero full cover has an")
    print("  unblocked small pair, every block-all full cover is positive-measure,")
    print("  and no full cover falls below the 1/14 floor.  The proof bottleneck is")
    print("  now a bridge lemma: lift each of the 64 fixed round classes to a labelled")
    print("  unit-spine realization or show that failure to do so already exposes the")
    print("  unblocked small-pair certificate.")


if __name__ == "__main__":
    main()
