#!/usr/bin/env python3
"""Exact information-preservation audit for LRC14 representations.

This is a small theorem-interface experiment, not a proof of LRC(14).  It asks
whether two natural binary carriers preserve three exact predicates:

* divisor covering (every q=2,...,14 divides a speed),
* the exact lonely-runner value M(S), and
* the THM-755 top-peel cap status v > r_P/(pi |G'_P|).

Carrier R is a tournament on runners using their mod-14 phase differences.
Carrier E is a tournament on the signed endpoints of the leave-top-out good
set, viewed at the peeled frequency.  Both are intentionally stripped of
metric labels.  The audit then attaches formula-facing sidecars:

* the divisor mask for covering;
* the projective cap ratio v|G'_P|/r_P for cap status; and
* an exact maximizing time/value pair for M(S).

Everything geometric and arithmetic is Fraction-exact.  The only appearance
of pi is decided rigorously between 333/106 < pi < 355/113; the generated bank
contains no unresolved comparison.  Endpoint Bernoulli reconstruction is
checked against the canonical lrc14_certificates.disc_exact implementation.

Tournament Analysis
-------------------
R pair observable: (v_j-v_i) mod 14.  Gauge: residues 1,...,6 orient i->j,
8,...,13 orient j->i.  Residues 0 and 7 are ties, oriented by the increasing-
speed Hamiltonian path.

E pair observable: {v(x_j-x_i)} for cyclically ordered signed endpoints x_i.
Gauge: the open half-circle (0,1/2) orients i->j; (1/2,1) reverses it.  Phase
0 and 1/2 ties follow the cyclic-endpoint Hamiltonian path.  Scores, directed
3-cycles, SCCs, Hamiltonian paths (R), and edge flips along controlled lifts
are reported.  Endpoint signs and exact phases are retained only in the B2
sidecar; the unweighted tournament does not retain discrepancy.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import comb, gcd
from typing import Callable, Hashable, Iterable

from lrc14_certificates import (
    LAM,
    M_exact,
    disc_exact,
    good_intervals,
    is_covering,
)


PI_LO = F(333, 106)
PI_HI = F(355, 113)


def fmt(x: F | None) -> str:
    if x is None:
        return "none"
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def gcd_all(values: Iterable[int]) -> int:
    g = 0
    for value in values:
        g = gcd(g, value)
    return g


def divisor_mask(speeds: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(int(any(v % q == 0 for v in speeds)) for q in range(2, 15))


def peak_witness(speeds: tuple[int, ...]) -> tuple[F, F]:
    """Return (first maximizing t, exact M) using the canonical peak set."""
    values = tuple(sorted(set(speeds)))
    denoms = {2 * v for v in values}
    for i, v in enumerate(values):
        for w in values[i + 1 :]:
            denoms.add(v + w)
            denoms.add(w - v)

    best_m = F(-1)
    best_t = F(0)
    for q in sorted(denoms):
        for m in range(1, q):
            t = F(m, q)
            clearance = min(
                min((v * m) % q, q - ((v * m) % q)) for v in values
            )
            candidate = F(clearance, q)
            if candidate > best_m or (candidate == best_m and t < best_t):
                best_m = candidate
                best_t = t
    assert best_m == M_exact(values)
    return best_t, best_m


@dataclass(frozen=True)
class Tournament:
    n: int
    edge_word: tuple[int, ...]
    score_hist: tuple[tuple[int, int], ...]
    cycles3: int
    scc_sizes: tuple[int, ...]
    tie_count: int
    hamiltonian_paths: int | None

    @property
    def fingerprint(self) -> tuple[Hashable, ...]:
        return (self.n, self.score_hist, self.cycles3, self.scc_sizes, self.hamiltonian_paths)


def adjacency_from_word(n: int, word: tuple[int, ...]) -> list[set[int]]:
    adj = [set() for _ in range(n)]
    cursor = 0
    for i in range(n):
        for j in range(i + 1, n):
            if word[cursor]:
                adj[i].add(j)
            else:
                adj[j].add(i)
            cursor += 1
    assert cursor == len(word)
    return adj


def scc_sizes(adj: list[set[int]]) -> tuple[int, ...]:
    n = len(adj)
    rev = [set() for _ in range(n)]
    for u, outs in enumerate(adj):
        for v in outs:
            rev[v].add(u)

    seen: set[int] = set()
    order: list[int] = []

    def visit(u: int) -> None:
        seen.add(u)
        for v in adj[u]:
            if v not in seen:
                visit(v)
        order.append(u)

    for u in range(n):
        if u not in seen:
            visit(u)

    seen.clear()
    sizes: list[int] = []

    def collect(u: int) -> int:
        seen.add(u)
        total = 1
        for v in rev[u]:
            if v not in seen:
                total += collect(v)
        return total

    for u in reversed(order):
        if u not in seen:
            sizes.append(collect(u))
    return tuple(sorted(sizes, reverse=True))


def hamiltonian_path_count(adj: list[set[int]]) -> int:
    """Exact number of directed Hamiltonian paths; used only at n=13."""
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if not count:
                continue
            for nxt in adj[last]:
                if not (mask >> nxt) & 1:
                    dp[mask | (1 << nxt)][nxt] += count
    return sum(dp[-1])


def tournament_from_word(n: int, word: tuple[int, ...], ties: int, hp: bool) -> Tournament:
    adj = adjacency_from_word(n, word)
    scores = [len(outs) for outs in adj]
    transitive_triples = sum(comb(score, 2) for score in scores)
    cycles = comb(n, 3) - transitive_triples
    return Tournament(
        n=n,
        edge_word=word,
        score_hist=tuple(sorted(Counter(scores).items())),
        cycles3=cycles,
        scc_sizes=scc_sizes(adj),
        tie_count=ties,
        hamiltonian_paths=hamiltonian_path_count(adj) if hp else None,
    )


@lru_cache(maxsize=None)
def runner_tournament(speeds: tuple[int, ...]) -> Tournament:
    word: list[int] = []
    ties = 0
    for i, j in combinations(range(len(speeds)), 2):
        delta = (speeds[j] - speeds[i]) % 14
        if delta in (0, 7):
            ties += 1
            word.append(1)  # increasing-speed tie Hamiltonian path
        else:
            word.append(int(1 <= delta <= 6))
    return tournament_from_word(len(speeds), tuple(word), ties, hp=True)


Endpoint = tuple[F, int, int]


@lru_cache(maxsize=None)
def endpoint_vertices(core: tuple[int, ...]) -> tuple[Endpoint, ...]:
    vertices: list[Endpoint] = []
    for interval_index, (left, right) in enumerate(good_intervals(core)):
        vertices.append((left, 1, interval_index))
        vertices.append((right, -1, interval_index))
    return tuple(sorted(vertices, key=lambda item: (item[0], -item[1], item[2])))


def endpoint_tournament(core: tuple[int, ...], peel: int) -> Tournament:
    endpoints = endpoint_vertices(core)
    word: list[int] = []
    ties = 0
    for i, j in combinations(range(len(endpoints)), 2):
        phase = (peel * (endpoints[j][0] - endpoints[i][0])) % 1
        if phase in (0, F(1, 2)):
            ties += 1
            word.append(1)  # cyclic endpoint tie Hamiltonian path
        else:
            word.append(int(phase < F(1, 2)))
    return tournament_from_word(len(endpoints), tuple(word), ties, hp=False)


def endpoint_b2_disc(core: tuple[int, ...], peel: int) -> F:
    endpoints = endpoint_vertices(core)

    def b2(x: F) -> F:
        return x * x - x + F(1, 6)

    total = F(0)
    for xp, sp, _ in endpoints:
        for xq, sq, _ in endpoints:
            total += sp * sq * b2((peel * (xp - xq)) % 1)
    return total / (2 * peel * peel)


def cap_status(cap_ratio: F) -> bool:
    """Decide pi*cap_ratio > 1 rigorously; fail if the narrow bounds straddle."""
    if PI_LO * cap_ratio > 1:
        return True
    if PI_HI * cap_ratio <= 1:
        return False
    raise ArithmeticError(f"pi bounds do not decide cap ratio {cap_ratio}")


@dataclass(frozen=True)
class Row:
    name: str
    family: str
    speeds: tuple[int, ...]
    covering: bool
    div_mask: tuple[int, ...]
    peak: tuple[F, F]
    cap: bool
    cap_ratio: F
    components: int
    good_measure: F
    discrepancy: F
    runner: Tournament
    endpoint: Tournament

    @property
    def M(self) -> F:
        return self.peak[1]


def build_bank(swap_max: int) -> list[tuple[str, str, tuple[int, ...]]]:
    ap = tuple(range(1, 14))
    raw: list[tuple[str, str, tuple[int, ...]]] = [
        ("AP", "boundary", ap),
        ("Goddyn-Wong", "boundary", tuple(list(range(1, 12)) + [13, 24])),
        ("deep-well", "single-killer", tuple(list(range(1, 13)) + [182])),
        ("residue-killer", "single-killer", tuple(list(range(1, 12)) + [13, 84])),
        ("translated-cluster", "cluster", tuple([1] + list(range(90, 102)))),
        (
            "multi-scale-comb",
            "multiscale",
            (1, 10, 21, 24, 56, 65, 77, 135, 219, 265, 335, 367, 390),
        ),
        (
            "compressed-ray-c26",
            "coherent-pack",
            tuple(list(range(26, 313, 26)) + [339]),
        ),
    ]
    for drop in range(1, 14):
        core = [v for v in ap if v != drop]
        for lift in range(14, swap_max + 1):
            speeds = tuple(sorted(core + [lift]))
            if len(set(speeds)) == 13 and gcd_all(speeds) == 1:
                raw.append((f"swap-d{drop}-a{lift}", f"single-swap-d{drop}", speeds))

    seen: set[tuple[int, ...]] = set()
    bank = []
    for name, family, speeds in raw:
        if speeds not in seen:
            seen.add(speeds)
            bank.append((name, family, speeds))
    return bank


def analyze(name: str, family: str, speeds: tuple[int, ...]) -> Row:
    peel = max(speeds)
    core = tuple(v for v in speeds if v != peel)
    intervals = endpoint_vertices(core)
    components = len(intervals) // 2
    good_measure = sum(
        intervals[i + 1][0] - intervals[i][0] for i in range(0, len(intervals), 2)
    )
    assert good_measure > 0 and components > 0
    ratio = F(peel) * good_measure / components
    peak = peak_witness(speeds)
    disc = endpoint_b2_disc(core, peel)
    assert disc == disc_exact(core, peel)
    mask = divisor_mask(speeds)
    return Row(
        name=name,
        family=family,
        speeds=speeds,
        covering=is_covering(speeds),
        div_mask=mask,
        peak=peak,
        cap=cap_status(ratio),
        cap_ratio=ratio,
        components=components,
        good_measure=good_measure,
        discrepancy=disc,
        runner=runner_tournament(speeds),
        endpoint=endpoint_tournament(core, peel),
    )


TARGETS: dict[str, Callable[[Row], Hashable]] = {
    "cover": lambda row: row.covering,
    "M": lambda row: row.M,
    "cap": lambda row: row.cap,
    "disc": lambda row: row.discrepancy,
}


@dataclass(frozen=True)
class Audit:
    name: str
    fibers: int
    collision_fibers: int
    mixed: tuple[int, int, int, int]


def audit(rows: list[Row], name: str, key: Callable[[Row], Hashable]) -> Audit:
    fibers: dict[Hashable, list[Row]] = defaultdict(list)
    for row in rows:
        fibers[key(row)].append(row)
    mixed = []
    for target in TARGETS.values():
        mixed.append(sum(len({target(row) for row in fiber}) > 1 for fiber in fibers.values()))
    return Audit(
        name=name,
        fibers=len(fibers),
        collision_fibers=sum(len(fiber) > 1 for fiber in fibers.values()),
        mixed=tuple(mixed),  # type: ignore[arg-type]
    )


def best_mixed_example(
    rows: list[Row], key: Callable[[Row], Hashable], targets: tuple[str, ...]
) -> tuple[Row, Row] | None:
    fibers: dict[Hashable, list[Row]] = defaultdict(list)
    for row in rows:
        fibers[key(row)].append(row)
    best: tuple[int, Row, Row] | None = None
    for fiber in fibers.values():
        for a, b in combinations(fiber, 2):
            differences = sum(TARGETS[target](a) != TARGETS[target](b) for target in targets)
            if differences and (best is None or differences > best[0]):
                best = (differences, a, b)
    return None if best is None else (best[1], best[2])


def flips(a: tuple[int, ...], b: tuple[int, ...]) -> int:
    assert len(a) == len(b)
    return sum(x != y for x, y in zip(a, b))


def flip_summary(rows: list[Row], step: int) -> tuple[int, int, int, F, tuple[Row, Row] | None]:
    by_family: dict[str, dict[int, Row]] = defaultdict(dict)
    for row in rows:
        if row.family.startswith("single-swap-d"):
            by_family[row.family][max(row.speeds)] = row
    values: list[int] = []
    cap_crossing: tuple[Row, Row] | None = None
    for family_rows in by_family.values():
        for lift, row in family_rows.items():
            other = family_rows.get(lift + step)
            if other is None:
                continue
            values.append(flips(row.endpoint.edge_word, other.endpoint.edge_word))
            if row.cap != other.cap and cap_crossing is None:
                cap_crossing = (row, other)
    return len(values), min(values), max(values), F(sum(values), len(values)), cap_crossing


def row_text(row: Row) -> str:
    return (
        f"{row.name}: S={row.speeds}, cover={row.covering}, M={fmt(row.M)} at "
        f"t={fmt(row.peak[0])}, cap={row.cap}, cap_ratio={fmt(row.cap_ratio)}, "
        f"r={row.components}, |G|={fmt(row.good_measure)}"
    )


def report(rows: list[Row], swap_max: int) -> None:
    runner_fp = lambda r: r.runner.fingerprint
    runner_full = lambda r: (r.runner.n, r.runner.edge_word)
    endpoint_fp = lambda r: r.endpoint.fingerprint
    endpoint_full = lambda r: (r.endpoint.n, r.endpoint.edge_word)
    div = lambda r: r.div_mask
    cap_sign = lambda r: r.cap
    cap_ratio = lambda r: r.cap_ratio
    peak = lambda r: r.peak

    audits = [
        audit(rows, "R fingerprint", runner_fp),
        audit(rows, "R full tournament", runner_full),
        audit(rows, "R + divisor mask", lambda r: (runner_full(r), div(r))),
        audit(rows, "R + cap sign", lambda r: (runner_full(r), cap_sign(r))),
        audit(rows, "R + divisor + cap sign", lambda r: (runner_full(r), div(r), cap_sign(r))),
        audit(rows, "R + div + cap sign + peak", lambda r: (runner_full(r), div(r), cap_sign(r), peak(r))),
        audit(rows, "R + exact cap ratio", lambda r: (runner_full(r), cap_ratio(r))),
        audit(rows, "E fingerprint", endpoint_fp),
        audit(rows, "E full tournament", endpoint_full),
        audit(rows, "E + divisor mask", lambda r: (endpoint_full(r), div(r))),
        audit(rows, "E + cap sign", lambda r: (endpoint_full(r), cap_sign(r))),
        audit(rows, "E + divisor + cap sign", lambda r: (endpoint_full(r), div(r), cap_sign(r))),
        audit(rows, "E + div + cap sign + peak", lambda r: (endpoint_full(r), div(r), cap_sign(r), peak(r))),
        audit(rows, "E + exact cap ratio", lambda r: (endpoint_full(r), cap_ratio(r))),
    ]

    print("LRC14 EXACT REPRESENTATION / SIDECAR AUDIT")
    print(f"bank: {len(rows)} primitive 13-speed rows (single swaps through {swap_max} + 7 landmarks)")
    print(f"covering rows: {sum(row.covering for row in rows)}")
    print(f"THM-755 cap-pass rows: {sum(row.cap for row in rows)}")
    print(f"distinct exact M values: {len({row.M for row in rows})}")
    print(f"endpoint B2 reconstructions checked exactly: {len(rows)}/{len(rows)}")
    print()
    print("MIXED-FIBER AUDIT (entries count fibers, not rows)")
    print(f"{'carrier':31s} {'fibers':>7s} {'coll':>6s} {'cover':>6s} {'M':>6s} {'cap':>6s} {'disc':>6s}")
    for item in audits:
        print(
            f"{item.name:31s} {item.fibers:7d} {item.collision_fibers:6d} "
            + " ".join(f"{value:6d}" for value in item.mixed)
        )

    print()
    print("EXACT COLLISION WITNESSES")
    examples = [
        ("R full", runner_full, ("cover", "M", "cap")),
        ("R + divisor + cap sign", lambda r: (runner_full(r), div(r), cap_sign(r)), ("M",)),
        ("E full", endpoint_full, ("cover", "M", "cap", "disc")),
        ("E + divisor + cap sign", lambda r: (endpoint_full(r), div(r), cap_sign(r)), ("M",)),
    ]
    for label, key, targets in examples:
        pair = best_mixed_example(rows, key, targets)
        print(f"[{label}]")
        if pair is None:
            print("  no mixed target fiber in this bank")
        else:
            print("  " + row_text(pair[0]))
            print("  " + row_text(pair[1]))

    runner_cycles = [row.runner.cycles3 for row in rows]
    runner_hps = [row.runner.hamiltonian_paths for row in rows]
    endpoint_cycles = [row.endpoint.cycles3 for row in rows]
    endpoint_vertices_n = [row.endpoint.n for row in rows]
    print()
    print("TOURNAMENT FINGERPRINTS")
    print("R observable=(v_j-v_i) mod 14; half-circle gauge; ties follow increasing speed")
    print(
        f"  score fingerprints={len({row.runner.score_hist for row in rows})}, "
        f"3-cycles range={min(runner_cycles)}..{max(runner_cycles)}, "
        f"SCC patterns={Counter(row.runner.scc_sizes for row in rows)}, "
        f"Hamiltonian paths range={min(runner_hps)}..{max(runner_hps)}"
    )
    print("E observable={v(x_j-x_i)}; half-circle gauge; ties follow cyclic endpoint order")
    print(
        f"  vertices range={min(endpoint_vertices_n)}..{max(endpoint_vertices_n)}, "
        f"3-cycles range={min(endpoint_cycles)}..{max(endpoint_cycles)}, "
        f"SCC patterns={Counter(row.endpoint.scc_sizes for row in rows).most_common(6)}, "
        f"ties range={min(row.endpoint.tie_count for row in rows)}..{max(row.endpoint.tie_count for row in rows)}"
    )
    print("  Hamiltonian-path counts omitted for E: endpoint tournaments reach hundreds of vertices.")

    base_core = tuple(range(1, 13))
    base_intervals = good_intervals(base_core)
    base_g = sum(right - left for left, right in base_intervals)
    base_r = len(base_intervals)
    base_ratio = F(13) * base_g / base_r
    scale_checks = []
    for scale in (1, 2, 3, 5):
        scaled_intervals = good_intervals(tuple(scale * v for v in base_core))
        scaled_g = sum(right - left for left, right in scaled_intervals)
        scaled_r = len(scaled_intervals)
        scaled_ratio = F(13 * scale) * scaled_g / scaled_r
        scale_checks.append(
            scaled_g == base_g and scaled_r == scale * base_r and scaled_ratio == base_ratio
        )
    print()
    print("PROJECTIVE SCALE CHECK")
    print(
        f"  P={{1,...,12}}, peel=13; scales 1,2,3,5: "
        f"|G_cP|=|G_P|, r_cP=c*r_P, cap_ratio invariant: {all(scale_checks)}"
    )
    print(f"  invariant cap_ratio={fmt(base_ratio)}; exact checks={sum(scale_checks)}/{len(scale_checks)}")

    print()
    print("EDGE FLIPS ON FIXED-CORE SINGLE-SWAP LIFTS")
    for step in (1, 14):
        count, minimum, maximum, mean, crossing = flip_summary(rows, step)
        runner_values = []
        by_family: dict[str, dict[int, Row]] = defaultdict(dict)
        for row in rows:
            if row.family.startswith("single-swap-d"):
                by_family[row.family][max(row.speeds)] = row
        for family_rows in by_family.values():
            for lift, row in family_rows.items():
                other = family_rows.get(lift + step)
                if other is not None:
                    runner_values.append(flips(row.runner.edge_word, other.runner.edge_word))
        print(
            f"step {step:2d}: pairs={count}, R flips={min(runner_values)}..{max(runner_values)} "
            f"(mean {fmt(F(sum(runner_values), len(runner_values)))}), "
            f"E flips={minimum}..{maximum} (mean {fmt(mean)})"
        )
        if crossing:
            print(
                f"  first cap-status flip: {crossing[0].name} -> {crossing[1].name}; "
                f"E edge flips={flips(crossing[0].endpoint.edge_word, crossing[1].endpoint.edge_word)}"
            )

    print()
    print("PRESERVATION VERDICT")
    print("  divisor mask recovers covering exactly; neither tournament alone does.")
    print("  projective ratio v|G|/r proves a one-bit cap-sign sidecar exactly.")
    print("  signed endpoint phases reconstruct Bernoulli discrepancy exactly; orientation alone does not.")
    print("  fixed-threshold endpoint data and both tournaments still mix exact M.")
    print("  the exact peak witness (t*,M) is the tested sidecar that restores M.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--swap-max", type=int, default=55)
    args = parser.parse_args()
    if args.swap_max < 28:
        parser.error("--swap-max must be at least 28 to expose mod-14 lift collisions")
    bank = build_bank(args.swap_max)
    rows = [analyze(*entry) for entry in bank]
    report(rows, args.swap_max)


if __name__ == "__main__":
    main()
