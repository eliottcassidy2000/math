#!/usr/bin/env python3
"""S586: extremality order for the LRC additive/sleeve branch.

The recent LRC spine now has several objects that can look "extreme":

* AP sleeves measure-saturate exactly at the last runner.
* visible 3-term folds depress the exact maximin and safe measure.
* shifted AP keeps maximal 4-term energy while becoming very safe.
* hidden 4-term fibers become real folds only after adding a virtual sum node.
* dyadic AP layers are self-similar, but only the odd/unit layer binds delta.

This script compares those axes on common rows.  The intended output is not a
counterexample search; it is a certificate ordering for which extremal features
should be allowed to drive a proof.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations, permutations
from math import gcd, sqrt
import random
import statistics


def dist(x: F) -> F:
    x %= 1
    return min(x, 1 - x)


def primitive(row: tuple[int, ...]) -> tuple[int, ...]:
    g = 0
    for v in row:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in row))


def v2(x: int) -> int:
    out = 0
    while x and x % 2 == 0:
        x //= 2
        out += 1
    return out


def exact_maximin(row: tuple[int, ...]) -> tuple[F, tuple[F, ...]]:
    """Exact maximin over the standard piecewise-linear candidate clocks."""
    candidates: set[F] = {F(0)}
    for a in row:
        for m in range(a):
            candidates.add(F(2 * m + 1, 2 * a))
    for a, b in combinations(row, 2):
        for denom in (a + b, abs(a - b)):
            if denom == 0:
                continue
            for m in range(1, denom):
                candidates.add(F(m, denom))

    best = F(0)
    witnesses: list[F] = []
    for t in candidates:
        margin = min(dist(v * t) for v in row)
        if margin > best:
            best = margin
            witnesses = [t]
        elif margin == best:
            witnesses.append(t)
    return best, tuple(sorted(witnesses))


def safe_measure(row: tuple[int, ...], n: int) -> F:
    delta = F(1, n)
    endpoints: set[F] = {F(0)}
    for v in row:
        for j in range(v + 1):
            endpoints.add((F(j) + delta) / v % 1)
            endpoints.add((F(j) - delta) / v % 1)

    pts = sorted(endpoints)
    total = F(0)
    for i, a in enumerate(pts):
        b = pts[(i + 1) % len(pts)]
        length = b - a if b > a else b - a + 1
        mid = (a + length / 2) % 1
        if all(dist(v * mid) >= delta for v in row):
            total += length
    return total


def saturation_curve(row: tuple[int, ...], n: int) -> tuple[F, ...]:
    return tuple(safe_measure(row[:j], n) for j in range(1, len(row) + 1))


def three_terms(row: tuple[int, ...]) -> tuple[tuple[int, int, int], ...]:
    vals = set(row)
    rels = []
    for a, b in combinations(sorted(row), 2):
        c = a + b
        if c in vals:
            rels.append((a, b, c))
    return tuple(rels)


def ordered_additive_energy(row: tuple[int, ...]) -> int:
    sums = Counter(a + b for a in row for b in row)
    return sum(count * count for count in sums.values())


def unordered_pair_sum_fibers(row: tuple[int, ...]) -> dict[int, tuple[tuple[int, int], ...]]:
    fibers: dict[int, list[tuple[int, int]]] = {}
    for a, b in combinations(sorted(row), 2):
        fibers.setdefault(a + b, []).append((a, b))
    return {s: tuple(pairs) for s, pairs in fibers.items() if len(pairs) >= 2}


def hidden_fibers(row: tuple[int, ...]) -> dict[int, tuple[tuple[int, int], ...]]:
    vals = set(row)
    return {s: pairs for s, pairs in unordered_pair_sum_fibers(row).items() if s not in vals}


@dataclass(frozen=True)
class HiddenPressure:
    sum_node: int | None
    multiplicity: int
    pair_edges: tuple[tuple[int, int], ...]
    augmented_M: F | None
    augmented_margin: F | None
    drop: F | None
    min_virtual_dist_at_witness: F | None


def hidden_pressure(row: tuple[int, ...], M: F, witnesses: tuple[F, ...]) -> HiddenPressure:
    hidden = hidden_fibers(row)
    if not hidden:
        return HiddenPressure(None, 0, tuple(), None, None, None, None)
    sum_node, pairs = max(hidden.items(), key=lambda item: (len(item[1]), item[0]))
    aug = tuple(sorted(row + (sum_node,)))
    aug_M, _ = exact_maximin(aug)
    aug_n = len(aug) + 1
    min_dist = min(dist(sum_node * t) for t in witnesses) if witnesses else None
    return HiddenPressure(
        sum_node=sum_node,
        multiplicity=len(pairs),
        pair_edges=pairs,
        augmented_M=aug_M,
        augmented_margin=aug_M - F(1, aug_n),
        drop=M - aug_M,
        min_virtual_dist_at_witness=min_dist,
    )


@dataclass(frozen=True)
class RowCase:
    label: str
    row: tuple[int, ...]
    n: int
    family: str


@dataclass(frozen=True)
class RowStats:
    case: RowCase
    M: F
    margin: F
    witnesses: tuple[F, ...]
    safe_mu: F
    curve: tuple[F, ...]
    t3: int
    energy: int
    hidden: HiddenPressure
    critical_saturation: bool
    every_sleeve_needed: bool
    clock_binders: tuple[int, ...]
    binder_v2: tuple[int, ...]


def clock_binders(row: tuple[int, ...], n: int) -> tuple[int, ...]:
    delta = F(1, n)
    return tuple(v for v in row if dist(F(v, n)) == delta)


def row_stats(case: RowCase) -> RowStats:
    M, witnesses = exact_maximin(case.row)
    mu = safe_measure(case.row, case.n)
    curve = saturation_curve(case.row, case.n)
    rels = three_terms(case.row)
    hidden = hidden_pressure(case.row, M, witnesses)
    sat_at = next((i + 1 for i, value in enumerate(curve) if value == 0), None)
    critical = sat_at == len(case.row)
    needed = False
    if curve and curve[-1] == 0:
        needed = all(safe_measure(tuple(v for v in case.row if v != drop), case.n) > 0 for drop in case.row)
    binders = clock_binders(case.row, case.n)
    return RowStats(
        case=case,
        M=M,
        margin=M - F(1, case.n),
        witnesses=witnesses,
        safe_mu=mu,
        curve=curve,
        t3=len(rels),
        energy=ordered_additive_energy(case.row),
        hidden=hidden,
        critical_saturation=critical,
        every_sleeve_needed=needed,
        clock_binders=binders,
        binder_v2=tuple(sorted({v2(v) for v in binders})),
    )


def fmt(x: F | None, digits: int = 5) -> str:
    if x is None:
        return "-"
    return f"{float(x):.{digits}f}"


def fmt_margin(x: F | None) -> str:
    if x is None:
        return "-"
    return f"{float(x):+.5f}"


def target_cases() -> list[RowCase]:
    cases: list[RowCase] = []
    for n in (8, 10, 14):
        cases.append(RowCase(f"AP_n{n}", tuple(range(1, n)), n, "critical_AP"))
        cases.append(RowCase(f"unit_shift_AP_n{n}", tuple(range(2, n + 1)), n, "fold_rich_but_clock_blocked"))
        cases.append(RowCase(f"far_shift_AP_n{n}", tuple(range(n + 1, 2 * n)), n, "raw_energy_decoy"))
    cases.extend(
        [
            RowCase("Vstar_n14", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24), 14, "known_floor_variant"),
            RowCase("doubled_apex_n14", tuple(range(1, 13)) + (26,), 14, "apex_stress"),
            RowCase("hidden_4fold_S584_k9", (6, 11, 14, 15, 16, 18, 19, 23, 28), 10, "hidden_virtual_fold"),
            RowCase("fixed_sum_hidden_k6", (2, 5, 6, 14, 15, 18), 7, "hidden_virtual_fold"),
            RowCase("second_gap_lift_n8", (1, 6, 7, 10, 17, 18, 19), 8, "open_gap_lift"),
        ]
    )
    return cases


def sample_rows(n: int, trials: int, rng: random.Random) -> list[RowCase]:
    k = n - 1
    rows: set[tuple[int, ...]] = set()
    universe = list(range(1, 5 * n + 1))
    attempts = 0
    while len(rows) < trials and attempts < 200 * trials:
        attempts += 1
        row = primitive(tuple(rng.sample(universe, k)))
        if len(row) == k:
            rows.add(row)
    return [RowCase(f"sample_n{n}_{i:03d}", row, n, "sample") for i, row in enumerate(sorted(rows))]


def pearson(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 2 or len(xs) != len(ys):
        return None
    mx = statistics.mean(xs)
    my = statistics.mean(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0:
        return None
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sqrt(vx * vy)


def family_table(stats: list[RowStats]) -> list[str]:
    lines = []
    lines.append("TARGETED EXTREMAL FAMILY TABLE")
    lines.append(
        "  label                    n  family                    M       M-delta   mu_safe  "
        "3t  energy hidden(s,r) aug_margin critical binders(v2)"
    )
    for s in stats:
        hp = s.hidden
        hidden_text = "-" if hp.sum_node is None else f"{hp.sum_node},{hp.multiplicity}"
        bind_text = "-" if not s.clock_binders else f"{list(s.clock_binders)}{list(s.binder_v2)}"
        lines.append(
            f"  {s.case.label:24s} {s.case.n:2d} {s.case.family:25s} "
            f"{fmt(s.M):>7s} {fmt_margin(s.margin):>9s} {fmt(s.safe_mu):>8s} "
            f"{s.t3:2d} {s.energy:6d} {hidden_text:>11s} {fmt_margin(hp.augmented_margin):>10s} "
            f"{str(s.critical_saturation and s.every_sleeve_needed):>8s} {bind_text}"
        )
    return lines


def saturation_report(stats: list[RowStats]) -> list[str]:
    lines = []
    lines.append("SATURATION-CURVE EXTREMALITY")
    for s in stats:
        if s.case.family not in {
            "critical_AP",
            "fold_rich_but_clock_blocked",
            "raw_energy_decoy",
            "known_floor_variant",
            "apex_stress",
        }:
            continue
        sat_at = next((i + 1 for i, value in enumerate(s.curve) if value == 0), None)
        last_drop = s.curve[-2] - s.curve[-1] if len(s.curve) >= 2 else s.curve[-1]
        curve_tail = [fmt(x, 3) for x in s.curve[-4:]]
        lines.append(
            f"  {s.case.label:24s}: final_mu={fmt(s.curve[-1])}, sat_at={sat_at}, "
            f"last_drop={fmt(last_drop)}, every_sleeve_needed={s.every_sleeve_needed}, tail={curve_tail}"
        )
    lines.append(
        "  reading: AP rows have the critical signature (zero final slack, first saturation at the last sleeve)."
    )
    lines.append(
        "  shifted AP keeps the same clock binders and raw energy as AP but has positive safe measure, so clock facade/energy are decoys."
    )
    return lines


def sample_audit() -> list[str]:
    rng = random.Random(585)
    lines = []
    lines.append("DETERMINISTIC SAMPLE AUDIT")
    lines.append("  hardness coordinate is -margin; covering hardness is -mu_safe.")
    for n, trials in ((8, 160), (10, 140), (14, 70)):
        stats = [row_stats(case) for case in sample_rows(n, trials, rng)]
        no3 = [s for s in stats if s.t3 == 0]
        has3 = [s for s in stats if s.t3 > 0]
        top4_cut = sorted((s.energy for s in no3), reverse=True)[max(0, len(no3) // 4 - 1)] if no3 else None
        high4_no3 = [s for s in no3 if top4_cut is not None and s.energy >= top4_cut]
        rich_cut = sorted((s.t3 for s in has3), reverse=True)[max(0, len(has3) // 4 - 1)] if has3 else None
        fold_rich = [s for s in has3 if rich_cut is not None and s.t3 >= rich_cut]

        def min_line(bucket: list[RowStats]) -> str:
            if not bucket:
                return "count=0"
            hard = min(bucket, key=lambda s: (s.margin, s.safe_mu, max(s.case.row)))
            return (
                f"count={len(bucket):3d}; min_margin={fmt_margin(hard.margin)}; "
                f"min_mu={fmt(hard.safe_mu)}; t3={hard.t3}; energy={hard.energy}; row={hard.case.row}"
            )

        hard_x = [-float(s.margin) for s in stats]
        cover_x = [-float(s.safe_mu) for s in stats]
        t3_y = [float(s.t3) for s in stats]
        energy_y = [float(s.energy) for s in stats]
        hidden_y = [float(s.hidden.multiplicity) for s in stats]
        lines.append(f"  n={n:2d} rows={len(stats)}")
        lines.append(f"    no_visible_3fold       {min_line(no3)}")
        lines.append(f"    high_4energy_no3       {min_line(high4_no3)}")
        lines.append(f"    visible_fold_rich      {min_line(fold_rich)}")
        lines.append(
            "    corr(-margin, 3t/energy/hidden_mult)="
            f"{pearson(hard_x, t3_y):+.3f}/{pearson(hard_x, energy_y):+.3f}/{pearson(hard_x, hidden_y):+.3f}; "
            "corr(-mu, 3t/energy)="
            f"{pearson(cover_x, t3_y):+.3f}/{pearson(cover_x, energy_y):+.3f}"
        )
    return lines


@dataclass(frozen=True)
class Axis:
    name: str
    certificate_power: int
    hardness_correlation: int
    extremal_specificity: int
    decoy_resistance: int
    maturity: int

    def key(self) -> tuple[int, int, int, int, int]:
        return (
            self.certificate_power,
            self.hardness_correlation,
            self.extremal_specificity,
            self.decoy_resistance,
            self.maturity,
        )


AXES = [
    Axis("critical_sleeve_saturation", 5, 5, 5, 5, 3),
    Axis("Phi_endpoint_gap_terminal", 5, 4, 4, 5, 5),
    Axis("visible_low_3fold_delta_clock", 4, 5, 4, 4, 4),
    Axis("circuit_free_margin_floor", 3, 4, 3, 5, 3),
    Axis("hidden_virtual_fold_pressure", 2, 3, 4, 3, 2),
    Axis("dyadic_fractal_binder_profile", 2, 2, 3, 3, 2),
    Axis("raw_4term_energy", 1, 1, 1, 1, 4),
]


def tournament_report() -> list[str]:
    names = [axis.name for axis in AXES]
    edge: dict[tuple[str, str], bool] = {}
    for a, b in combinations(AXES, 2):
        winner, loser = (a, b) if a.key() > b.key() else (b, a)
        edge[(winner.name, loser.name)] = True
        edge[(loser.name, winner.name)] = False

    scores = Counter()
    for u in names:
        scores[sum(1 for v in names if u != v and edge[(u, v)])] += 1

    c3 = 0
    for triple in combinations(names, 3):
        for a, b, c in permutations(triple):
            if edge[(a, b)] and edge[(b, c)] and edge[(c, a)]:
                c3 += 1
                break

    hpaths = 0
    first_path = None
    for path in permutations(names):
        if all(edge[(path[i], path[i + 1])] for i in range(len(path) - 1)):
            hpaths += 1
            if first_path is None:
                first_path = path

    ranking = sorted(AXES, key=lambda axis: axis.key(), reverse=True)
    lines = []
    lines.append("TOURNAMENT ANALYSIS: EXTREMALITY AXES")
    lines.append("  vertices: extremality/certificate axes, not runners.")
    lines.append(
        "  pairwise observable: (certificate_power, hardness_correlation, extremal_specificity, decoy_resistance, maturity)."
    )
    lines.append("  switch/gauge: lexicographically larger observable gets the arc; declaration order is the tie Hamiltonian path.")
    lines.append(f"  ranking: {[axis.name for axis in ranking]}")
    lines.append(f"  score_histogram: {dict(sorted(scores.items()))}")
    lines.append(f"  directed_3_cycles: {c3}")
    lines.append(f"  Hamiltonian_path_count: {hpaths}")
    lines.append(f"  first_Hamiltonian_path: {first_path}")
    lines.append(
        "  transitivity reading: raw energy cannot shortcut hidden pressure or visible folds without losing certificate power."
    )
    return lines


def assumption_challenge() -> list[str]:
    return [
        "ASSUMPTION CHALLENGE",
        "  Candidate vertices considered: runners, gaps, fixed circle sections, section boundaries, wall crossings, residues, cover arcs, Fourier modes, matroid circuits, hidden sum nodes, sleeve-order positions, proof obligations, and extremality axes.",
        "  Chosen quotient: extremality axes plus calibrated witness rows.",
        "  Preserved LRC predicate: which coordinate can actually force M or safe measure down to the delta floor.",
        "  Destroyed information: exact endpoint-owner geometry, residue languages, and component-by-component Phi ramps; those must re-enter through HYP-2112/HYP-2108/HYP-2107 terminal gates.",
        "  Challenged assumption: an extremal scalar statistic is automatically a hard LRC statistic. Shifted AP refutes this for raw 4-term energy.",
    ]


def proof_program() -> list[str]:
    return [
        "EXTREMALITY PROOF PROGRAM",
        "  E1. Prove critical sleeve saturation implies a visible low-fold/delta-clock scaffold, not merely high additive energy.",
        "  E2. Prove hidden-only virtual fold pressure stays in Lemma-A territory until the virtual sum becomes a real runner or endpoint gate.",
        "  E3. For n=14, classify floor rows by critical saturation plus endpoint/Phi terminal gates; AP and V* are calibration rows, while unit-shift and far-shift APs are decoys.",
        "  E4. Treat the dyadic/tetrational-looking AP fractal as a binder filter: even layers repeat structure, but odd/unit binders carry the delta contact.",
    ]


def main() -> None:
    targeted = [row_stats(case) for case in target_cases()]
    print("S586 LRC extremality-order session")
    print("Lens: separate certificate-bearing extremality from decoy extremality.")
    print()
    for block in (
        family_table(targeted),
        [""],
        saturation_report(targeted),
        [""],
        sample_audit(),
        [""],
        tournament_report(),
        [""],
        assumption_challenge(),
        [""],
        proof_program(),
    ):
        for line in block:
            print(line)


if __name__ == "__main__":
    main()
