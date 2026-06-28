#!/usr/bin/env python3
"""HYP-3432: Euler-Mascheroni harmonic wall-budget scout for LRC14.

Euler's constant enters here only as a warning: the expression
H_N - log(N) is a renormalized harmonic tail, not an exact LRC certificate.
For the finite LRC14 proof route we therefore keep the log/gamma analogy out
of the computation and audit exact rational reciprocal budgets attached to
HYP-3427 wall signatures and HYP-3429 endpoint spines.

The scout asks whether a Mertens-like harmonic owner budget can prioritize
endpoint debt while still proving that any scalar budget must remain a sidecar:
the actual finite certificate is still an endpoint window plus branch/wall
labels.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from importlib import util
from itertools import combinations
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
H3427_PATH = ROOT / "04-computation" / "lrc14_two_branch_wall_signature_atlas_codex_20260628.py"
H3429_PATH = ROOT / "04-computation" / "lrc14_component_spine_certificate_codex_20260628.py"


def load_module(path: Path, name: str):
    spec = util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


h3427 = load_module(H3427_PATH, "h3427_wall_signature_atlas")
h3429 = load_module(H3429_PATH, "h3429_component_spine")


def fmt(x: F | None) -> str:
    if x is None:
        return "n/a"
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def harmonic(n: int) -> F:
    return sum((F(1, k) for k in range(1, n + 1)), F(0))


def row_reciprocal(speeds: tuple[int, ...]) -> F:
    return sum((F(1, v) for v in speeds), F(0))


def label_text_structured(labels: tuple[tuple[str, int], ...]) -> str:
    if not labels:
        return "none"
    return ",".join(f"{kind}:{value}" for kind, value in labels)


def owner_values_structured(labels: tuple[tuple[str, int], ...]) -> tuple[int, ...]:
    return tuple(sorted({value for _kind, value in labels if value > 0}))


def owner_budget_structured(labels: tuple[tuple[str, int], ...]) -> F:
    return sum((F(1, value) for value in owner_values_structured(labels)), F(0))


def parse_wall_label(label: str) -> tuple[str, int] | None:
    if label in {"FREE", "END:0", "END:1"}:
        return None
    kind, value = label.split(":", 1)
    return kind, int(value)


def wall_label_tokens(window) -> tuple[tuple[str, int], ...]:
    labels: set[tuple[str, int]] = set()
    for label in window.left + window.right:
        parsed = parse_wall_label(label)
        if parsed is not None:
            labels.add(parsed)
    return tuple(sorted(labels))


def owner_values_wall(window) -> tuple[int, ...]:
    return tuple(sorted({value for _kind, value in wall_label_tokens(window)}))


def owner_budget_wall(window) -> F:
    return sum((F(1, value) for value in owner_values_wall(window)), F(0))


@dataclass(frozen=True)
class SpineBudget:
    name: str
    speeds: tuple[int, ...]
    branches: tuple[int, ...]
    labels: tuple[tuple[str, int], ...]
    budget: F
    row_sum: F
    owner_share: F
    length: F

    @property
    def kinds(self) -> tuple[str, ...]:
        return tuple(sorted({kind for kind, _value in self.labels}))

    @property
    def word(self) -> str:
        return label_text_structured(self.labels)


@dataclass(frozen=True)
class WallBudget:
    row: str
    mask: str
    left_types: str
    right_types: str
    labels: tuple[tuple[str, int], ...]
    budget: F
    width: F

    @property
    def signature(self) -> tuple[str, str, str]:
        return (self.mask, self.left_types, self.right_types)

    @property
    def word(self) -> str:
        return ",".join(f"{kind}:{value}" for kind, value in self.labels)


def spine_budget_rows() -> list[SpineBudget]:
    rows: list[SpineBudget] = []
    for name, speeds in h3429.audited_rows():
        audit = h3429.audit_row(name, speeds)
        budget = owner_budget_structured(audit.best.labels)
        reciprocal = row_reciprocal(audit.speeds)
        rows.append(
            SpineBudget(
                name=audit.name,
                speeds=audit.speeds,
                branches=audit.best.branches,
                labels=audit.best.labels,
                budget=budget,
                row_sum=reciprocal,
                owner_share=budget / reciprocal,
                length=audit.best.length,
            )
        )
    return rows


def wall_budget_rows() -> list[WallBudget]:
    rows: list[WallBudget] = []
    for name, speeds in h3427.audited_rows():
        audit = h3427.audit(name, speeds)
        for window in audit.windows:
            rows.append(
                WallBudget(
                    row=name,
                    mask=window.mask,
                    left_types=h3427.label_types(window.left),
                    right_types=h3427.label_types(window.right),
                    labels=wall_label_tokens(window),
                    budget=owner_budget_wall(window),
                    width=window.width,
                )
            )
    return rows


def budget_collision_summary(records: list[WallBudget]) -> tuple[int, list[tuple[F, int, int, list[tuple[str, str, str]], list[str]]]]:
    by_budget: dict[F, list[WallBudget]] = defaultdict(list)
    for record in records:
        by_budget[record.budget].append(record)

    collisions: list[tuple[F, int, int, list[tuple[str, str, str]], list[str]]] = []
    for budget, bucket in by_budget.items():
        signatures = sorted({record.signature for record in bucket})
        words = sorted({record.word for record in bucket})
        if len(signatures) > 1:
            sample_rows = sorted({record.row for record in bucket})[:4]
            collisions.append((budget, len(bucket), len(signatures), signatures[:5], sample_rows))
    collisions.sort(key=lambda item: (-item[2], -item[1], item[0]))
    return len(collisions), collisions


def spine_collision_summary(rows: list[SpineBudget]) -> tuple[int, list[tuple[F, int, list[tuple[tuple[int, ...], tuple[str, ...]]], list[str]]]]:
    by_budget: dict[F, list[SpineBudget]] = defaultdict(list)
    for row in rows:
        by_budget[row.budget].append(row)

    collisions: list[tuple[F, int, list[tuple[tuple[int, ...], tuple[str, ...]]], list[str]]] = []
    for budget, bucket in by_budget.items():
        shapes = sorted({(row.branches, row.kinds) for row in bucket})
        if len(shapes) > 1:
            names = sorted(row.name for row in bucket)[:5]
            collisions.append((budget, len(bucket), shapes[:5], names))
    collisions.sort(key=lambda item: (-len(item[2]), -item[1], item[0]))
    return len(collisions), collisions


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: tuple[int, ...]

    @property
    def total(self) -> int:
        return sum(self.scores)


AXES = (
    "predicate_retention",
    "wall_exactness",
    "tail_budget_signal",
    "two_adic_compatibility",
    "collision_safety",
    "proof_readiness",
)

CARRIERS = (
    Carrier("endpoint_wall_certificate", (10, 10, 7, 9, 10, 10)),
    Carrier("component_spine_rank2", (10, 9, 7, 9, 9, 9)),
    Carrier("two_adic_loss_ledger", (9, 8, 7, 10, 8, 9)),
    Carrier("harmonic_owner_budget_sidecar", (6, 5, 10, 7, 5, 6)),
    Carrier("mertens_reciprocal_tail_rank", (5, 4, 9, 6, 4, 5)),
    Carrier("gamma_shadow_normalizer", (3, 2, 8, 3, 2, 2)),
    Carrier("named_constant_scalar_shortcut", (1, 0, 6, 1, 0, 0)),
)


def tournament() -> tuple[dict[int, int], list[str], int]:
    hist = dict(sorted(Counter(carrier.total for carrier in CARRIERS).items()))
    ordered = sorted(CARRIERS, key=lambda carrier: (-carrier.total, CARRIERS.index(carrier)))
    cycles = 0
    for triple in combinations(CARRIERS, 3):
        ranked = sorted(triple, key=lambda carrier: (-carrier.total, CARRIERS.index(carrier)))
        if ranked[0].total < ranked[1].total < ranked[2].total:
            cycles += 1
    return hist, [carrier.name for carrier in ordered], cycles


def main() -> None:
    spine_rows = spine_budget_rows()
    wall_rows = wall_budget_rows()
    h14 = harmonic(14)
    resonant_14_budget = F(1, 2) + F(1, 7) + F(1, 14)
    unit_14_budget = sum((F(1, k) for k in range(1, 15) if h3427.gcd(k, 14) == 1), F(0))

    spine_budget_hist = Counter(row.budget for row in spine_rows)
    spine_kind_hist = Counter(row.kinds for row in spine_rows)
    spine_collision_count, spine_collisions = spine_collision_summary(spine_rows)
    wall_collision_count, wall_collisions = budget_collision_summary(wall_rows)
    wall_budget_hist = Counter(row.budget for row in wall_rows)

    min_spine = min(spine_rows, key=lambda row: (row.budget, row.length, row.name))
    max_spine = max(spine_rows, key=lambda row: (row.budget, -row.length, row.name))
    min_share = min(spine_rows, key=lambda row: (row.owner_share, row.name))
    max_share = max(spine_rows, key=lambda row: (row.owner_share, row.name))

    tight_wall_audit = h3427.audit(
        "covering_AP_with_84",
        tuple(list(range(1, 12)) + [13, 84]),
    )

    print("HYP-3432 EULER-MASCHERONI HARMONIC WALL-BUDGET SCOUT")
    print("=" * 78)
    print("Identity:")
    print("  Euler-Mascheroni enters only as the H_N - log(N) tail analogy.")
    print("  The executable audit uses exact rational reciprocal budgets over endpoint owners.")
    print("  Verdict target: harmonic budgets may rank wall debt, but cannot replace walls.")
    print()

    print("A. Finite harmonic normalizers")
    print(f"  H_14 exact:                         {fmt(h14)}")
    print(f"  reciprocal budget of 2,7,14:        {fmt(resonant_14_budget)}")
    print(f"  reciprocal budget of units mod 14:  {fmt(unit_14_budget)}")
    print("  No floating log or gamma value is used; this is a rational sidecar audit.")
    print()

    print("B. HYP-3429 endpoint-spine harmonic budgets")
    print(f"  endpoint-spine rows audited:        {len(spine_rows)}")
    print(f"  distinct best-spine budgets:        {len(spine_budget_hist)}")
    print(f"  best-spine kind histogram:          {dict(sorted(spine_kind_hist.items()))}")
    print(f"  scalar budget shape-collisions:     {spine_collision_count}")
    print(
        "  minimum best-spine budget:          "
        f"{fmt(min_spine.budget)} at {min_spine.name} labels={min_spine.word}"
    )
    print(
        "  maximum best-spine budget:          "
        f"{fmt(max_spine.budget)} at {max_spine.name} labels={max_spine.word}"
    )
    print(
        "  minimum endpoint reciprocal share:  "
        f"{fmt(min_share.owner_share)} at {min_share.name}"
    )
    print(
        "  maximum endpoint reciprocal share:  "
        f"{fmt(max_share.owner_share)} at {max_share.name}"
    )
    print("  top best-spine scalar collisions:")
    for budget, count, shapes, names in spine_collisions[:6]:
        print(f"    budget={fmt(budget)} rows={count} shapes={shapes} sample_rows={names}")
    if not spine_collisions:
        print("    none on best spines; collisions appear on the full wall atlas below.")
    print()

    print("C. HYP-3427 all-window wall-budget collision audit")
    print(f"  survivor windows audited:           {len(wall_rows)}")
    print(f"  distinct wall-owner budgets:        {len(wall_budget_hist)}")
    print(f"  budgets with >1 wall signature:     {wall_collision_count}")
    print("  largest scalar collisions:")
    for budget, count, sig_count, signatures, sample_rows in wall_collisions[:8]:
        print(
            f"    budget={fmt(budget)} windows={count} signature_types={sig_count} "
            f"sample_signatures={signatures} sample_rows={sample_rows}"
        )
    print()

    print("D. Tight canonical row as a harmonic sidecar")
    print("  row={1..11,13,84}; harmonic labels are transparent, walls are decisive.")
    for idx, window in enumerate(tight_wall_audit.windows, start=1):
        labels = wall_label_tokens(window)
        owners = owner_values_wall(window)
        print(
            f"  W{idx}: mask={window.mask} width={fmt(window.width)} "
            f"owners={owners} budget={fmt(owner_budget_wall(window))} labels={labels}"
        )
    print()

    print("E. Proof interpretation")
    print("  The harmonic/Mertens budget is useful as an owner-debt priority queue:")
    print("    small budget = far-tail endpoint owners, often E-only spine debt;")
    print("    large budget = low-owner mixed walls, often the true obstruction boundary.")
    print("  But the scalar budget collides across branch masks and wall signatures.")
    print("  Therefore an Euler-Mascheroni-style tail normalizer can only be a sidecar")
    print("  attached to HYP-3427/HYP-3429 endpoint windows, not a floor certificate.")
    print()

    hist, path, cycles = tournament()
    print("F. Tournament Analysis")
    print("  vertices are proof carriers and tail sidecars, not runners.")
    print(f"  axes={','.join(AXES)}")
    print(f"  score_hist={hist}")
    print(f"  directed_3cycles={cycles}")
    print("  hamiltonian_path=" + " -> ".join(path))
    print()

    print("G. Assumption challenge")
    print("  Considered vertices: runners, gaps, circle sections, endpoint walls,")
    print("  wall-crossing events, owner reciprocal budgets, Mertens tails, gamma")
    print("  shadows, component spines, and proof obligations.")
    print("  The chosen quotient preserves the LRC predicate only when paired with")
    print("  an actual survivor window and branch.  It destroys branch orientation,")
    print("  endpoint type, interval width, and the exact wall word if used alone.")


if __name__ == "__main__":
    main()
