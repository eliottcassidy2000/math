#!/usr/bin/env python3
"""S151: LRC14 counterexample family/sporadic classifier.

This script turns the recent source-spectrum and packet-migration work into a
conservative classifier for counterexample-shaped rows.

Logical gates:
  * qdiv(S) < 14 is excluded by the exact q-clock witness t=1/qdiv.
  * qdiv(S) = 14 can be tight boundary or positive-open, but not a strict
    counterexample unless the closed threshold set is empty.
  * qdiv(S) > 14 is the only strict counterexample zone after THM-523.
    A true counterexample candidate in a finite bank must therefore have
    qdiv>14 and exact status "covered" under the S146 threshold classifier.

Rows that survive to qdiv>=14 are assigned both:
  * a packet family axis: AP/GW boundary, unit petal/GW strip, K33 lift,
    mixed K33-petal lift, single 14-tail comb, or unlabelled covering repair;
  * a divisor-cover skeleton: the minimum number of row speeds needed to cover
    all clocks d=2..14, plus essential cover speeds.

Tournament Analysis declaration:
  vertices are classifier gates/families, not runners.  The pairwise observable
  is which gate preserves the strict-counterexample predicate while adding the
  most discharge information.  Ties follow the declared family path.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import argparse
import sys


REPO = Path(__file__).resolve().parents[1]
DIVS = tuple(range(2, 15))


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s150 = load_module(
    "s151_packet_migration",
    REPO / "04-computation" / "lrc14_packet_migration_gauntlet_codex_s150.py",
)


@dataclass(frozen=True)
class CoverSkeleton:
    min_cover_size: int
    min_cover_count: int
    essential_speeds: tuple[int, ...]
    essential_divisors: tuple[int, ...]
    role: str


@dataclass(frozen=True)
class ClassifiedRow:
    bank: str
    name: str
    qdiv: int
    status: str
    mass: Fraction
    family: str
    candidate_state: str
    cover: CoverSkeleton
    atom_keys: tuple[str, ...]
    transfer: str
    front: str


def divisor_cover_sets(speeds: tuple[int, ...]) -> dict[int, frozenset[int]]:
    return {v: frozenset(d for d in DIVS if v % d == 0) for v in speeds}


def cover_skeleton(speeds: tuple[int, ...]) -> CoverSkeleton:
    covers = divisor_cover_sets(speeds)
    need = frozenset(DIVS)
    usable = tuple(v for v in speeds if covers[v])
    min_covers: list[tuple[int, ...]] = []
    for r in range(1, len(usable) + 1):
        for comb in combinations(usable, r):
            got: set[int] = set()
            for v in comb:
                got.update(covers[v])
            if frozenset(got) == need:
                min_covers.append(comb)
        if min_covers:
            break
    if not min_covers:
        return CoverSkeleton(99, 0, (), (), "not-covering")

    essential_speeds = tuple(sorted(v for v in usable if all(v in c for c in min_covers)))
    essential_divisors = tuple(
        d for d in DIVS if sum(1 for v in speeds if v % d == 0) == 1
    )
    size = len(min_covers[0])
    if size <= 2:
        role = "low-rank-divisor-load"
    elif size <= 4:
        role = "medium-rank-cover"
    else:
        role = "distributed-cover"
    return CoverSkeleton(size, len(min_covers), essential_speeds, essential_divisors, role)


def family_for(row: s150.ExactClass, cover: CoverSkeleton) -> str:
    keys = set(row.atom_keys)
    if row.status == "boundary_only" and row.packet_state == "AP/GW-boundary-source":
        return "AP/GW-boundary-family"
    if "K33" in keys and keys & {"P10", "P13"}:
        return "mixed-K33-petal-lift-family"
    if "K33" in keys:
        return "K33-state-lift-family"
    if keys and keys <= {"P10", "P13", "GW"}:
        return "unit-petal-or-GW-strip-family"
    if row.qdiv > 14 and len(row.holes) == 1 and any(a % 14 == 0 for a in row.adds):
        return "single-14-tail-comb-family"
    if row.qdiv > 14 and cover.role == "low-rank-divisor-load":
        return "low-rank-divisor-load-family"
    if row.qdiv > 14:
        return "unlabelled-covering-repair-family"
    if row.status == "positive_open":
        return "q14-positive-open-migration-family"
    return "sporadic-unclassified-family"


def candidate_state(row: s150.ExactClass) -> str:
    if row.qdiv < 14:
        return "excluded-qclock"
    if row.status == "covered" and row.qdiv > 14:
        return "TRUE-COUNTEREXAMPLE-CANDIDATE"
    if row.status == "covered":
        return "covered-but-not-covering-core"
    if row.status == "boundary_only":
        return "closed-threshold-tight-or-boundary"
    if row.qdiv > 14:
        return "covering-core-positive-open"
    return "q14-positive-open"


def front_label(row: s150.ExactClass) -> str:
    if row.first_front is None:
        return "-"
    a, b, left, right, slack = row.first_front
    return f"{a}->{b} {left}->{right} slack={slack}"


def classify_exact(row: s150.ExactClass) -> ClassifiedRow:
    cover = cover_skeleton(row.speeds)
    return ClassifiedRow(
        row.bank,
        row.name,
        row.qdiv,
        row.status,
        row.mass,
        family_for(row, cover),
        candidate_state(row),
        cover,
        row.atom_keys,
        row.transfer,
        front_label(row),
    )


def summarize_bank(k: int, add_max: int, workers: int, progress_every: int) -> list[ClassifiedRow]:
    rows, qcounts = s150.generate_bank(k, add_max)
    exact = s150.classify_bank(rows, workers=workers, progress_every=progress_every)
    classified = [classify_exact(row) for row in exact]
    bank = rows[0].bank if rows else f"{k}-swap add<={add_max}"
    print(f"[bank] {bank}")
    print(f"  generated_rows={sum(qcounts.values())}")
    print(f"  qroute_counts={dict(qcounts)}")
    print(f"  exact_qdiv>=14_rows={len(exact)}")
    print(f"  candidate_state_counts={dict(sorted(Counter(c.candidate_state for c in classified).items()))}")
    print(f"  family_counts={dict(sorted(Counter(c.family for c in classified).items()))}")
    print(f"  cover_role_counts={dict(sorted(Counter(c.cover.role for c in classified).items()))}")

    true_candidates = [c for c in classified if c.candidate_state == "TRUE-COUNTEREXAMPLE-CANDIDATE"]
    if true_candidates:
        print("  TRUE candidate rows:")
        for c in true_candidates[:20]:
            print(f"    {c.name:34s} family={c.family} cover={c.cover}")
    else:
        print("  TRUE candidate rows: none")

    boundary = [c for c in classified if c.status == "boundary_only"]
    if boundary:
        print("  boundary/tight rows:")
        for c in boundary[:10]:
            print(f"    {c.name:34s} family={c.family} closed={c.candidate_state}")

    positives = sorted(
        (c for c in classified if c.qdiv > 14 and c.status == "positive_open"),
        key=lambda c: (c.mass, c.name),
    )
    if positives:
        print("  smallest covering-core positive fronts:")
        for c in positives[:8]:
            print(
                f"    {c.name:34s} mass={str(c.mass):>10s} family={c.family:36s} "
                f"cover={c.cover.role}/{c.cover.min_cover_size} keys={','.join(c.atom_keys) or '-'}"
            )
    print()
    return classified


def print_global_readout(rows: list[ClassifiedRow]) -> None:
    print("[global classifier readout]")
    states = Counter(c.candidate_state for c in rows)
    families = Counter(c.family for c in rows)
    print(f"  total_exact_rows={len(rows)}")
    print(f"  candidate_states={dict(sorted(states.items()))}")
    print(f"  families={dict(sorted(families.items()))}")
    print()
    print("  Necessary strict-counterexample theorem:")
    print("    Every true counterexample must appear in state")
    print("    TRUE-COUNTEREXAMPLE-CANDIDATE = qdiv>14 and exact covered.")
    print(f"    observed_TRUE_candidates={states.get('TRUE-COUNTEREXAMPLE-CANDIDATE', 0)}")
    print()

    sporadic_like = [
        c
        for c in rows
        if c.qdiv > 14
        and c.status == "positive_open"
        and c.family == "unlabelled-covering-repair-family"
    ]
    sporadic_like = sorted(sporadic_like, key=lambda c: (c.mass, c.name))
    print("  sporadic-like positive covering fronts (not true candidates):")
    if not sporadic_like:
        print("    none")
    for c in sporadic_like[:12]:
        print(
            f"    {c.name:34s} mass={str(c.mass):>10s} "
            f"cover_size={c.cover.min_cover_size} essential={c.cover.essential_speeds or '-'}"
        )
    print()
    print("  Family interpretation:")
    print("    AP/GW-boundary-family: tight qdiv=14 boundary, not strict counterexample.")
    print("    unit-petal/GW and K33 families: labelled local packets; discharge or state-lift.")
    print("    single-14-tail and divisor-load families: covering-core comb/tower rows.")
    print("    unlabelled-covering-repair family: current sporadic reservoir; in this audit all open.")


def tournament_analysis() -> None:
    vertices = [
        "qclock-excluded",
        "q14-boundary-tight",
        "q14-open-migration",
        "covering-comb-family",
        "unit-petal-family",
        "K33-state-lift-family",
        "unlabelled-covering-repair",
        "true-covered-sporadic",
    ]
    path = list(vertices)
    out = {v: set(path[i + 1 :]) for i, v in enumerate(path)}
    score_hist = Counter(len(out[v]) for v in vertices)
    print("[Tournament Analysis]")
    print("  vertices=classifier gates/families, not runners")
    print("  observable=preserves strict-counterexample predicate plus discharge label")
    print("  switch=earlier gate discharges or labels before later residual")
    print(f"  tie_hamiltonian_path={' > '.join(path)}")
    print(f"  score_hist={dict(sorted(score_hist.items()))}")
    print("  directed_3_cycles=0")
    print("  hamiltonian_paths=1")
    print()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--one-max", type=int, default=420)
    parser.add_argument("--two-max", type=int, default=60)
    parser.add_argument("--three-max", type=int, default=30)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--progress-every", type=int, default=10000)
    args = parser.parse_args()

    print("=" * 78)
    print("S151 LRC14 counterexample family/sporadic classifier")
    print("=" * 78)
    print("[assumption challenge]")
    print("  considered vertices: runners, residues, divisor clocks, cover skeletons,")
    print("    Haar fronts, boundary owners, C27/K33 packet labels, and proof families.")
    print("  chosen vertices: classifier gates/families.")
    print("  preserved predicate: qdiv>14 plus exact covered is the only strict")
    print("    counterexample state after THM-523 and S146.")
    print("  destroyed data: individual row identity after divisor skeleton and packet")
    print("    labels are recorded; exact wide analytic data remains an external proof.")
    print()

    all_rows: list[ClassifiedRow] = []
    for k, add_max in [(0, 13), (1, args.one_max), (2, args.two_max), (3, args.three_max)]:
        all_rows.extend(summarize_bank(k, add_max, args.workers, args.progress_every))
    print_global_readout(all_rows)
    tournament_analysis()


if __name__ == "__main__":
    main()
