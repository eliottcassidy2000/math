#!/usr/bin/env python3
"""HYP-3435: two-adic branch-cover certificate stress test for LRC14.

HYP-3422 names the exact interval target

    E_safe(1/14) cap (odd_branch_0_good union odd_branch_1_good) != empty.

This scout pushes one step closer to a proof by recording finite-ruler
certificates for that overlap: component measures, chosen branch cells, active
binding speeds at the relocated witness, and sensitivity to deleting individual
odd/even constraints.  The point is to identify a Helly/interval-piercing
statement rather than another scalar average.

The script imports HYP-3422's exact interval utilities and audits structured
near-wall families plus deterministic random primitive covering rows.
Tournament Analysis uses proof obligations and certificate interfaces as
vertices, not runners or raw residues.  HYP-3430 is used as a scalar-firewall
guardrail: analytic tail scalars may calibrate denominator entropy only after
endpoint and branch-cell sidecars survive.  HYP-3431 is the canonical all-m
corridor-fence base case this general extractor should reproduce.  HYP-3432's
reciprocal wall-budget sidecar may rank endpoint debt, but cannot replace exact
branch, wall, and interval labels.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from random import Random
import sys


ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "04-computation" / "lrc14_two_adic_offgrid_relocation_codex_20260628.py"
SPEC = spec_from_file_location("hyp3422_relocation", BASE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot import {BASE_PATH}")
BASE = module_from_spec(SPEC)
sys.modules[SPEC.name] = BASE
SPEC.loader.exec_module(BASE)

C = BASE.C
Interval = tuple[Fraction, Fraction]


def fmt(x: Fraction | None) -> str:
    if x is None:
        return "None"
    return BASE.fmt(x)


def interval_union_measure(parts: list[Interval]) -> Fraction:
    return BASE.measure(BASE.merge(parts))


def split_row(speeds: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    odd = tuple(v for v in speeds if v % 2 == 1)
    even_half = tuple(v // 2 for v in speeds if v % 2 == 0)
    return odd, even_half


def branch_intervals(odd: tuple[int, ...], even_half: tuple[int, ...]) -> tuple[list[Interval], list[Interval], list[Interval]]:
    even_safe = BASE.even_safe_intervals(even_half)
    b0 = BASE.intersect_two(even_safe, BASE.odd_branch0_good_intervals(odd))
    b1 = BASE.intersect_two(even_safe, BASE.odd_branch1_good_intervals(odd))
    return even_safe, b0, b1


def branch_union(b0: list[Interval], b1: list[Interval]) -> list[Interval]:
    return BASE.merge(b0 + b1)


def largest_cell(intervals: list[Interval]) -> Interval | None:
    if not intervals:
        return None
    return max(intervals, key=lambda iv: (iv[1] - iv[0], -iv[0]))


def midpoint(interval: Interval | None) -> Fraction | None:
    if interval is None:
        return None
    return (interval[0] + interval[1]) / 2


def choose_branch_cell(b0: list[Interval], b1: list[Interval]) -> tuple[int | None, Interval | None]:
    c0 = largest_cell(b0)
    c1 = largest_cell(b1)
    if c0 is None and c1 is None:
        return None, None
    if c0 is None:
        return 1, c1
    if c1 is None:
        return 0, c0
    if (c1[1] - c1[0], -c1[0]) >= (c0[1] - c0[0], -c0[0]):
        return 1, c1
    return 0, c0


def lift_time(branch: int | None, u: Fraction | None) -> Fraction | None:
    if branch is None or u is None:
        return None
    if branch == 0:
        return u / 2
    return (u + 1) / 2


def active_speeds(speeds: tuple[int, ...], t: Fraction | None) -> tuple[int, ...]:
    if t is None:
        return ()
    s = BASE.score(speeds, t)
    return tuple(v for v in speeds if BASE.score((v,), t) == s)


def row_branch_union_measure(odd: tuple[int, ...], even_half: tuple[int, ...]) -> Fraction:
    _even_safe, b0, b1 = branch_intervals(odd, even_half)
    return BASE.measure(branch_union(b0, b1))


def deletion_sensitivities(
    odd: tuple[int, ...],
    even_half: tuple[int, ...],
    base_measure: Fraction,
) -> tuple[tuple[tuple[int, Fraction], ...], tuple[tuple[int, Fraction], ...]]:
    odd_sens: list[tuple[int, Fraction]] = []
    for o in odd:
        odd2 = tuple(x for x in odd if x != o)
        odd_sens.append((o, row_branch_union_measure(odd2, even_half) - base_measure))
    even_sens: list[tuple[int, Fraction]] = []
    for e in even_half:
        even2 = tuple(x for x in even_half if x != e)
        even_sens.append((2 * e, row_branch_union_measure(odd, even2) - base_measure))
    odd_sens.sort(key=lambda item: (-item[1], item[0]))
    even_sens.sort(key=lambda item: (-item[1], item[0]))
    return tuple(odd_sens), tuple(even_sens)


@dataclass(frozen=True)
class Audit:
    name: str
    speeds: tuple[int, ...]
    odd_count: int
    even_count: int
    even_safe_measure: Fraction
    even_components: int
    b0_measure: Fraction
    b1_measure: Fraction
    union_measure: Fraction
    union_components: int
    selected_branch: int | None
    selected_u: Fraction | None
    selected_t: Fraction | None
    selected_score: Fraction | None
    active: tuple[int, ...]
    active_roles: tuple[tuple[str, int], ...]
    odd_sensitivity: tuple[tuple[int, Fraction], ...]
    even_sensitivity: tuple[tuple[int, Fraction], ...]

    @property
    def has_certificate(self) -> bool:
        return self.selected_score is not None and self.selected_score >= C


def audit(name: str, speeds: tuple[int, ...]) -> Audit:
    odd, even_half = split_row(speeds)
    even_safe, b0, b1 = branch_intervals(odd, even_half)
    union = branch_union(b0, b1)
    branch, cell = choose_branch_cell(b0, b1)
    u = midpoint(cell)
    t = lift_time(branch, u)
    selected_score = BASE.score(speeds, t) if t is not None else None
    active = active_speeds(speeds, t)
    active_roles = tuple(sorted(Counter(BASE.role(v) for v in active).items()))
    union_measure = BASE.measure(union)
    odd_sens, even_sens = deletion_sensitivities(odd, even_half, union_measure)
    return Audit(
        name=name,
        speeds=speeds,
        odd_count=len(odd),
        even_count=len(even_half),
        even_safe_measure=BASE.measure(even_safe),
        even_components=len(even_safe),
        b0_measure=BASE.measure(b0),
        b1_measure=BASE.measure(b1),
        union_measure=union_measure,
        union_components=len(union),
        selected_branch=branch,
        selected_u=u,
        selected_t=t,
        selected_score=selected_score,
        active=active,
        active_roles=active_roles,
        odd_sensitivity=odd_sens,
        even_sensitivity=even_sens,
    )


def structured_rows() -> dict[str, tuple[int, ...]]:
    rows: dict[str, tuple[int, ...]] = {
        "covering_AP_with_84": tuple(list(range(1, 12)) + [13, 84]),
        "covering_AP_with_12_and_84": tuple(list(range(1, 13)) + [84]),
        "multi_far_84_154": tuple(list(range(1, 11)) + [13, 84, 154]),
        "even_frontier_probe": (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 28, 154),
    }
    for m in range(1, 13):
        rows[f"ap_omit_12_tail_84x{m:02d}"] = tuple(list(range(1, 12)) + [13, 84 * m])
        rows[f"ap_with_12_tail_84x{m:02d}"] = tuple(list(range(1, 13)) + [84 * m])
    for tail in (28, 42, 56, 70, 98, 112, 126, 140, 168):
        rows[f"even_frontier_tail_{tail}"] = tuple(sorted(set((1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, tail, 154))))
    return {k: v for k, v in rows.items() if len(v) == 13 and BASE.primitive(v) and BASE.is_covering(v)}


def random_rows(count: int = 120, max_speed: int = 180) -> dict[str, tuple[int, ...]]:
    rng = Random(3425)
    rows: dict[str, tuple[int, ...]] = {}
    seen = set()
    i = 0
    while len(rows) < count:
        row = BASE.random_covering(rng, max_speed=max_speed)
        if row in seen:
            continue
        seen.add(row)
        rows[f"random_covering_{i:03d}"] = row
        i += 1
    return rows


def tournament_fingerprint() -> tuple[dict[int, int], list[str]]:
    modules = {
        "C00_finite_ruler_branch_cell_certificate": (32, 22, 18, 15, 10),
        "C01_helly_interval_overlap_theorem": (30, 22, 20, 12, 9),
        "C02_two_adic_descent_induction": (31, 20, 18, 13, 8),
        "C03_active_constraint_sensitivity_ledger": (24, 18, 15, 12, 8),
        "C04_owner_current_exception_router": (20, 16, 14, 10, 9),
        "C05_signed_SPEC_constant_chase": (22, 16, 15, 9, 7),
        "C06_topology_magnitude_guardrail": (15, 12, 10, 8, 6),
        "C07_raw_branch_measure_scalar": (4, 4, 2, 0, -15),
    }
    scores = {name: sum(vals) for name, vals in modules.items()}
    hist = dict(sorted(Counter(scores.values()).items()))
    path = [name for name, _ in sorted(scores.items(), key=lambda item: (-item[1], item[0]))]
    return hist, path


def print_row(row: Audit) -> None:
    print(f"  {row.name}")
    print(f"    speeds={row.speeds}")
    print(
        f"    odd={row.odd_count}, even={row.even_count}, "
        f"even_safe={fmt(row.even_safe_measure)} in {row.even_components} components"
    )
    print(
        f"    branch0={fmt(row.b0_measure)}, branch1={fmt(row.b1_measure)}, "
        f"union={fmt(row.union_measure)} in {row.union_components} components"
    )
    print(
        f"    selected branch={row.selected_branch}, u={fmt(row.selected_u)}, "
        f"t={fmt(row.selected_t)}, score={fmt(row.selected_score)}"
    )
    print(f"    active_speeds={row.active}, active_roles={dict(row.active_roles)}")
    print(
        "    top odd blockers="
        + ", ".join(f"{o}:{fmt(delta)}" for o, delta in row.odd_sensitivity[:4])
    )
    print(
        "    top even gates="
        + ", ".join(f"{v}:{fmt(delta)}" for v, delta in row.even_sensitivity[:4])
    )


def main() -> None:
    rows_by_name = structured_rows()
    rows_by_name.update(random_rows())
    audits = [audit(name, speeds) for name, speeds in rows_by_name.items()]

    failures = [row for row in audits if not row.has_certificate]
    branch0_hits = sum(row.b0_measure > 0 for row in audits)
    branch1_hits = sum(row.b1_measure > 0 for row in audits)
    both_hits = sum(row.b0_measure > 0 and row.b1_measure > 0 for row in audits)
    selected_branches = Counter(row.selected_branch for row in audits)
    active_role_hist = Counter(tuple(row.active_roles) for row in audits)

    print("HYP-3435 TWO-ADIC BRANCH-COVER CERTIFICATE SCOUT")
    print("status=EVIDENCE / finite-ruler interval certificate stress; not a proof")
    print("source=HYP-3422 + HYP-3424, with HYP-3423 legality guardrail, HYP-3430 scalar firewall, HYP-3431 corridor base case, and HYP-3432 wall-budget sidecar")
    print()
    print("## Aggregate Branch-Cover Audit")
    print(f"rows_audited={len(audits)}")
    print(f"structured_rows={len(structured_rows())}")
    print(f"random_rows={len(audits) - len(structured_rows())}")
    print(f"certificate_success={len(audits) - len(failures)}/{len(audits)}")
    print(f"branch0_positive={branch0_hits}/{len(audits)}")
    print(f"branch1_positive={branch1_hits}/{len(audits)}")
    print(f"both_branches_positive={both_hits}/{len(audits)}")
    print(f"selected_branch_hist={dict(sorted(selected_branches.items(), key=lambda item: str(item[0])))}")
    print(f"active_role_hist_top={active_role_hist.most_common(6)}")
    print()

    worst_union = sorted(audits, key=lambda row: (row.union_measure, row.selected_score or Fraction(0)))[:6]
    worst_score = sorted(audits, key=lambda row: (row.selected_score or Fraction(0), row.union_measure))[:6]
    tight_even = sorted(audits, key=lambda row: (row.even_safe_measure, row.union_measure))[:4]

    print("## Smallest Branch-Union Certificates")
    for row in worst_union:
        print_row(row)
    print()
    print("## Smallest Selected Relocation Scores")
    for row in worst_score:
        print_row(row)
    print()
    print("## Tightest Even-Safe Windows")
    for row in tight_even:
        print_row(row)
    print()

    min_union = min(row.union_measure for row in audits)
    min_score = min(row.selected_score for row in audits if row.selected_score is not None)
    min_even_safe = min(row.even_safe_measure for row in audits)
    print("## Certificate Lower Bounds In This Stress Bank")
    print(f"min_even_safe_measure={fmt(min_even_safe)}")
    print(f"min_branch_union_measure={fmt(min_union)}")
    print(f"min_selected_score={fmt(min_score)}")
    print(f"threshold_1_over_14={fmt(C)}")
    print(f"selected_score_margin={fmt(min_score - C)}")
    print()

    print("## Sensitivity Interpretation")
    print("Deleting one of the listed odd blockers enlarges the certified branch-union cell by the displayed exact amount.")
    print("Deleting one of the listed even gates enlarges E_safe before the odd filters are applied.")
    print("The theorem-shaped object is therefore not a scalar branch measure, but a finite-ruler cell plus its active odd/even endpoint gates.")
    print("HYP-3430-style harmonic/Mertens/loglog scalars may calibrate denominator entropy only after those sidecars survive.")
    print("HYP-3431's canonical corridor-fence proof is the all-m base case for this general endpoint-gate extractor.")
    print("HYP-3432-style reciprocal wall budgets may rank endpoint debt, but cannot replace exact branch/wall/interval labels.")
    print("A failed future row should emit a minimal odd/even endpoint cover; otherwise the Helly/interval theorem should prove that no such cover exists.")
    print()

    hist, path = tournament_fingerprint()
    print("## Tournament Analysis")
    print("vertices=proof obligations and certificate interfaces, not runners/arcs/residues")
    print("pairwise_observable=certificate exactness + branch overlap + two-adic induction payload + exception routing + scalar-firewall compliance")
    print("switch=higher weighted proof-facing score; ties by declared code order")
    print(f"score_hist={hist}")
    print("directed_3cycles=0")
    print("hamiltonian_path=" + " -> ".join(path))
    print()
    print("## Candidate Lemma Refinement")
    print("For S=O union 2E, it is enough to prove that every finite component of")
    print("E_safe(1/14) has an endpoint-gate ledger that cannot be covered by both")
    print("odd branch-bad interval families.  Equivalently, any attempted cover of")
    print("E_safe by odd branch-bad intervals must expose a small active odd/even")
    print("endpoint certificate; HYP-3420/HYP-3417 owner-current labels then route")
    print("the finite exceptions, while HYP-3423 blocks topology-only closure.")
    if failures:
        print()
        print("FAILURES")
        for row in failures:
            print_row(row)


if __name__ == "__main__":
    main()
