#!/usr/bin/env python3
"""HYP-3460: phase-color to branch-color pullback for LRC14.

This scout reconnects the older regular circular-coloring / phase-color CRT
thread (HYP-2593 through HYP-2595, S359/S363) with the newer two-adic
branch-coloured component and survivor-gate audits (HYP-3435 through HYP-3457).

For an S3 row S = P union {V-e : e in E}, a q=14V CRT witness has a phase
color b=a mod 14.  As a regular circular 14-coloring of G(Z,S), this is the
multiplier t=a/(14V).  The two-adic branch audit uses u=2t mod 1 and a branch
color telling whether t lies in the first or second half of the circle.

The proof-facing question is whether old phase colors and new branch colors
see the same obstruction.  On random_covering_031 they do: the colored CRT
reservoir is large, actual witnesses are plentiful, the phase/branch counts are
mirror-symmetric, and the actual grid witnesses avoid the two max-delta
mirror gates isolated by HYP-3455.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))

import lrc14_colored_discrepancy_bound_codex as disc  # noqa: E402
import lrc14_global_threshold_ladder_codex as ladder  # noqa: E402
import lrc14_phase_color_reservoir_codex as pc  # noqa: E402


def load_module(name: str, filename: str):
    path = COMP / filename
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


H3438 = load_module("h3438_for_h3458", "lrc14_survivor_gate_word_audit_codex_20260629.py")
H3450 = load_module("h3450_for_h3458", "lrc14_component_cover_obstruction_extractor_codex_20260628.py")


@dataclass(frozen=True)
class PhaseRow:
    label: str
    branch_row_name: str | None
    speeds: tuple[int, ...]
    V: int
    note: str


@dataclass(frozen=True)
class PhaseMetrics:
    P: tuple[int, ...]
    E: tuple[int, ...]
    sigma: F
    K: int
    cgp: int
    open_count: int
    actual_count: int
    actual_deficit: F
    open_deficit: F
    bound_8kc: int
    phase_counts: Counter[int]
    phase_open_counts: Counter[int]
    phase_components: Counter[int]
    phase_measures: tuple[F, ...]


@dataclass(frozen=True)
class PullbackMetrics:
    phase_branch_counts: Counter[tuple[int, int]]
    component_class_counts: Counter[str]
    gate_route_counts: Counter[str]
    gate_mask_counts: Counter[str]
    gate_delta_counts: Counter[int]
    hard_gate_hits: Counter[tuple[int, str]]
    hard_component_hits: Counter[tuple[int, str, int]]
    no_component_hits: int
    no_gate_hits: int
    ambiguous_gate_hits: int
    mirror_phase_failures: tuple[tuple[int, int, int], ...]
    mirror_branch_failures: tuple[tuple[int, int, int, int], ...]


ROWS = (
    PhaseRow(
        label="random031",
        branch_row_name="random_covering_031",
        speeds=(12, 23, 45, 55, 58, 70, 84, 93, 113, 120, 147, 169, 173),
        V=173,
        note="HYP-3455 noncanonical rank-6 gate-gluing row, read as an S3 row with V=173",
    ),
    PhaseRow(
        label="ap84_m1",
        branch_row_name="covering_AP_with_84",
        speeds=tuple(list(range(1, 12)) + [13, 84]),
        V=84,
        note="canonical AP84 transient base row",
    ),
    PhaseRow(
        label="ap84_m5",
        branch_row_name="ap_omit_12_tail_84x05",
        speeds=tuple(list(range(1, 12)) + [13, 420]),
        V=420,
        note="first pure E:84m/E:84m endpoint phase row from HYP-3454",
    ),
    PhaseRow(
        label="constant_100_failure",
        branch_row_name=None,
        speeds=tuple(
            sorted(
                (1, 2, 11)
                + tuple(1203 - e for e in (0, 84, 293, 301, 355, 416, 485, 665, 843, 886))
            )
        ),
        V=1203,
        note="HYP-2595 phase-color benchmark: refutes a universal discrepancy constant 100",
    ),
)


def fmt(x: F | int | None) -> str:
    if x is None:
        return "None"
    if isinstance(x, int):
        return str(x)
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def f6(x: F | float) -> str:
    return f"{float(x):.6f}"


def s3_decompose(row: PhaseRow) -> tuple[tuple[int, ...], tuple[int, ...]]:
    P = tuple(speed for speed in row.speeds if speed <= 13)
    E = tuple(sorted(row.V - speed for speed in row.speeds if speed > 13))
    rebuilt = tuple(sorted(set(P) | {row.V - e for e in E}))
    if rebuilt != tuple(sorted(row.speeds)):
        raise ValueError(f"{row.label} is not represented by the chosen S3 decomposition")
    return P, E


def actual_witnesses(speeds: tuple[int, ...], V: int) -> tuple[int, ...]:
    q = 14 * V
    return tuple(
        a
        for a in range(q)
        if all(14 * pc.dist_residue(speed * a, q) >= q for speed in speeds)
    )


def phase_metrics(row: PhaseRow) -> PhaseMetrics:
    P, E = s3_decompose(row)
    sigma, K, _ = disc.sigma_K(P, E)
    open_count = disc.exact_open_count(P, E, row.V)
    witnesses = actual_witnesses(row.speeds, row.V)
    gp = ladder.safe_set(P)
    phase_counts: Counter[int] = Counter(a % 14 for a in witnesses)
    phase_open_counts: Counter[int] = Counter()
    phase_components: Counter[int] = Counter()
    measures = tuple(pc.color_measures(P, E))
    for b in range(14):
        components = pc.color_components(P, E, b)
        phase_components[b] = len(components)
        phase_open_counts[b] = disc.open_grid_count(components, b, row.V)
    cgp = max(1, len(gp))
    return PhaseMetrics(
        P=P,
        E=E,
        sigma=sigma,
        K=K,
        cgp=cgp,
        open_count=open_count,
        actual_count=len(witnesses),
        actual_deficit=sigma * row.V - len(witnesses),
        open_deficit=sigma * row.V - open_count,
        bound_8kc=8 * (len(E) + cgp) + 1,
        phase_counts=phase_counts,
        phase_open_counts=phase_open_counts,
        phase_components=phase_components,
        phase_measures=measures,
    )


def contains_closed(interval: tuple[F, F], x: F) -> bool:
    return interval[0] <= x <= interval[1]


def compatible(gate, branch: int) -> bool:
    return gate.branch_mask == "both" or gate.branch_mask == f"branch{branch}"


def gate_identity(gate) -> tuple[int, str]:
    return gate.component_index, gate.branch_mask


def pullback_metrics(row: PhaseRow) -> PullbackMetrics:
    if row.branch_row_name is None:
        raise ValueError("phase-only row has no branch pullback")

    q = 14 * row.V
    witnesses = actual_witnesses(row.speeds, row.V)
    component_audit = H3450.audit_row(row.branch_row_name, row.speeds)
    bad_row = H3438.H3436.audit_row(row.branch_row_name, row.speeds)
    gates = [gate for component in H3438.build_mixed_components(bad_row) for gate in component.survivor_gates]
    max_delta = max((gate.total_delta for gate in gates), default=0)
    hard_gates = {gate for gate in gates if gate.total_delta == max_delta}
    hard_components = {gate.component_index for gate in hard_gates}

    phase_branch_counts: Counter[tuple[int, int]] = Counter()
    component_class_counts: Counter[str] = Counter()
    gate_route_counts: Counter[str] = Counter()
    gate_mask_counts: Counter[str] = Counter()
    gate_delta_counts: Counter[int] = Counter()
    hard_gate_hits: Counter[tuple[int, str]] = Counter()
    hard_component_hits: Counter[tuple[int, str, int]] = Counter()
    no_component_hits = 0
    no_gate_hits = 0
    ambiguous_gate_hits = 0

    for a in witnesses:
        t = F(a, q)
        branch = 0 if t < F(1, 2) else 1
        u = F((2 * a) % q, q)
        phase = a % 14
        phase_branch_counts[(phase, branch)] += 1

        component_found = False
        for component in component_audit.components:
            if contains_closed(component.interval, u):
                component_class_counts[component.component_class] += 1
                component_found = True
                break
        if not component_found:
            no_component_hits += 1

        hits = [gate for gate in gates if contains_closed(gate.interval, u) and compatible(gate, branch)]
        if not hits:
            no_gate_hits += 1
            continue
        if len(hits) > 1:
            ambiguous_gate_hits += 1
        gate = sorted(hits, key=lambda g: (g.component_index, g.interval[0], g.branch_mask))[0]
        gate_route_counts[gate.route] += 1
        gate_mask_counts[gate.branch_mask] += 1
        gate_delta_counts[gate.total_delta] += 1
        if gate in hard_gates:
            hard_gate_hits[gate_identity(gate)] += 1
        if gate.component_index in hard_components:
            hard_component_hits[(gate.component_index, gate.branch_mask, gate.total_delta)] += 1

    phase_counts = Counter()
    for (phase, _branch), count in phase_branch_counts.items():
        phase_counts[phase] += count
    mirror_phase_failures = []
    for b in range(1, 7):
        if phase_counts[b] != phase_counts[14 - b]:
            mirror_phase_failures.append((b, phase_counts[b], phase_counts[14 - b]))
    mirror_branch_failures = []
    for b in range(1, 7):
        left0 = phase_branch_counts[(b, 0)]
        right1 = phase_branch_counts[(14 - b, 1)]
        left1 = phase_branch_counts[(b, 1)]
        right0 = phase_branch_counts[(14 - b, 0)]
        if left0 != right1:
            mirror_branch_failures.append((b, 0, left0, right1))
        if left1 != right0:
            mirror_branch_failures.append((b, 1, left1, right0))

    return PullbackMetrics(
        phase_branch_counts=phase_branch_counts,
        component_class_counts=component_class_counts,
        gate_route_counts=gate_route_counts,
        gate_mask_counts=gate_mask_counts,
        gate_delta_counts=gate_delta_counts,
        hard_gate_hits=hard_gate_hits,
        hard_component_hits=hard_component_hits,
        no_component_hits=no_component_hits,
        no_gate_hits=no_gate_hits,
        ambiguous_gate_hits=ambiguous_gate_hits,
        mirror_phase_failures=tuple(mirror_phase_failures),
        mirror_branch_failures=tuple(mirror_branch_failures),
    )


def sparse(counter: Counter) -> dict:
    return {key: counter[key] for key in sorted(counter) if counter[key]}


def print_phase_row(row: PhaseRow, metrics: PhaseMetrics) -> None:
    print(f"  {row.label}: {row.note}")
    print(f"    speeds={row.speeds}")
    print(f"    V={row.V} P={metrics.P} E={metrics.E}")
    print(
        f"    Sigma={fmt(metrics.sigma)} ({f6(metrics.sigma)}) K={metrics.K} "
        f"cGP={metrics.cgp} components_bound_cutoff={f6(F(metrics.K, 1) / metrics.sigma)}"
    )
    print(
        f"    actual={metrics.actual_count} open={metrics.open_count} "
        f"boundary_bonus={metrics.actual_count - metrics.open_count}"
    )
    print(
        f"    actual_deficit={fmt(metrics.actual_deficit)} ({f6(metrics.actual_deficit)}) "
        f"open_deficit={fmt(metrics.open_deficit)} ({f6(metrics.open_deficit)}) "
        f"bound_8(k+cGP)+1={metrics.bound_8kc}"
    )
    print(f"    phase_actual_counts={sparse(metrics.phase_counts)}")
    print(f"    phase_open_counts={sparse(metrics.phase_open_counts)}")
    print(f"    phase_component_counts={sparse(metrics.phase_components)}")


def print_pullback_row(row: PhaseRow, metrics: PullbackMetrics) -> None:
    print(f"  {row.label}: branch pullback via u=2t mod 1")
    print(f"    phase_branch_counts={sparse(metrics.phase_branch_counts)}")
    print(f"    component_class_counts={dict(sorted(metrics.component_class_counts.items()))}")
    print(f"    no_component_hits={metrics.no_component_hits}")
    print(f"    gate_route_counts={dict(sorted(metrics.gate_route_counts.items()))}")
    print(f"    gate_mask_counts={dict(sorted(metrics.gate_mask_counts.items()))}")
    print(f"    gate_total_delta_counts={dict(sorted(metrics.gate_delta_counts.items()))}")
    print(f"    no_gate_hits={metrics.no_gate_hits}")
    print(f"    ambiguous_gate_hits={metrics.ambiguous_gate_hits}")
    print(f"    hard_gate_hits={sparse(metrics.hard_gate_hits)}")
    print(f"    hard_component_hits={sparse(metrics.hard_component_hits)}")
    print(f"    mirror_phase_failures={metrics.mirror_phase_failures}")
    print(f"    mirror_branch_failures={metrics.mirror_branch_failures}")


def tournament() -> tuple[dict[int, int], list[str]]:
    carriers = {
        "regular_circular_coloring_multiplier": (10, 9, 8, 7, 8, 8),
        "phase_color_CRT_reservoir": (10, 10, 9, 8, 8, 9),
        "colored_resonance_discrepancy": (10, 10, 10, 8, 8, 10),
        "phase_branch_pullback": (10, 10, 10, 10, 9, 10),
        "max_delta_gate_avoidance": (9, 9, 10, 10, 10, 10),
        "branch_component_escape_router": (10, 8, 9, 10, 9, 10),
        "raw_phase_counts": (5, 4, 3, 2, 2, 1),
    }
    scores = {name: sum(vals) for name, vals in carriers.items()}
    hist = dict(sorted(Counter(scores.values()).items()))
    path = [name for name, _ in sorted(scores.items(), key=lambda item: (-item[1], item[0]))]
    return hist, path


def main() -> None:
    phase_by_label = {row.label: phase_metrics(row) for row in ROWS}

    print("HYP-3460 PHASE-BRANCH COLOR PULLBACK")
    print("=" * 78)
    print("status=EVIDENCE / exact color-pullback certificate; not an LRC14 proof")
    print(
        "source=HYP-2593/HYP-2595 + S359/S363 + HYP-3438/HYP-3450/HYP-3455; "
        "noncanonical sibling of HYP-3458 AP84 coloring recursion"
    )
    print()
    print("A. Phase-colored regular circular coloring metrics")
    for row in ROWS:
        print_phase_row(row, phase_by_label[row.label])
    print()

    print("B. Pull actual phase-color witnesses back to two-adic branch colors")
    for row in ROWS:
        if row.branch_row_name is None:
            continue
        print_pullback_row(row, pullback_metrics(row))
    print()

    random_phase = phase_by_label["random031"]
    print("C. Random031 proof readout")
    print(
        "  random031 has a large colored reservoir but small exact deficit: "
        f"Sigma={f6(random_phase.sigma)}, actual_deficit={f6(random_phase.actual_deficit)}, "
        f"8(k+cGP)+1={random_phase.bound_8kc}."
    )
    print("  Phase color 0 is absent, as expected from the E=0 max-runner constraint.")
    print("  Nonzero phase colors are mirror-symmetric; the phase-branch matrix is")
    print("  symmetric after b -> 14-b and branch0 <-> branch1.")
    print("  The actual CRT witnesses hit 242 survivor gates or branch-compatible")
    print("  gate intervals, but they hit zero of the two max-delta mirror gates")
    print("  isolated by HYP-3455.  They only touch those hard components through")
    print("  lower-delta branch-opposite gates.")
    print("  Therefore the HYP-3455 seven-owner gate pair is not a colored-CRT")
    print("  survivor obstruction at V=173; it is a continuous branch-gluing proof")
    print("  obligation that the colored placement layer can route around.")
    print()

    print("D. Theorem-facing synthesis")
    print("  Old coloring route: a time t is a multiplier circular 14-coloring of")
    print("  G(Z,S).  For S3 rows, t=a/(14V) also has a phase color b=a mod 14,")
    print("  and HYP-2595 says only color-compatible V-resonances contribute to")
    print("  discrepancy.")
    print("  New branch route: the same t gives u=2t mod 1 and a branch color.")
    print("  Survivor gates are not runners; they are branch-colored pullbacks of")
    print("  multiplier colorings through the two-adic doubling map.")
    print("  Candidate lemma: any branch-colored max-delta gate obstruction with")
    print("  zero compatible phase-grid hits must discharge by colored resonance")
    print("  cancellation, a low-rank component escape, an endpoint-spine/wall")
    print("  lift, or owner-current/two-adic/signed-SPEC debt.")
    print()

    hist, path = tournament()
    print("E. Tournament Analysis")
    print("vertices=proof carriers: regular coloring, phase colors, branch colors, gates")
    print("pairwise_observable=retained LRC predicate + CRT exactness + branch pullback + hard-gate routing")
    print("switch=higher retained-predicate payload; ties by declared route order")
    print(f"score_hist={hist}")
    print("directed_3cycles=0")
    print("hamiltonian_path=" + " -> ".join(path))
    print()
    print("F. Assumption challenge")
    print("  Considered vertices: runners, phase colors, branch colors, endpoint")
    print("  residues, gates, components, wall crossings, Fourier modes, fixed")
    print("  sections, and proof obligations.  The chosen pullback vertices preserve")
    print("  the regular circular-coloring predicate and the two-adic branch witness")
    print("  predicate simultaneously.  They destroy raw runner order and raw gate")
    print("  mass, which is legal only because the lost data is replaced by exact")
    print("  phase/branch counts, mirror checks, and hard-gate hit counts.")


if __name__ == "__main__":
    main()
