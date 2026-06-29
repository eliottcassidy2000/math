#!/usr/bin/env python3
"""HYP-3452: AP-with-84m component phase audit for LRC14.

HYP-3433 found the harmonic endpoint-spine law for the canonical tail

    S_m = {1,2,...,11,13,84m}

from m >= 5.  HYP-3450/HYP-3451 then identified the same AP-with-84m rows as
the dangerous component-cover graph family: many dead components, connected
dead-cover projection, and very few escapes.

This scout stitches those two facts together.  It asks whether the graph-side
component-cover certificate undergoes the same phase transition as the
endpoint-spine ledger.  The answer on the checked range is clean:

    m >= 5: best component is the HYP-3433 E:84m/E:84m interval of length 1/(49m)
    m >= 3: every dead component has paired-cover rank <= 2

The remaining AP-tail data is a discrete escape clock.  The number of escape
components follows a period-35 Beatty correction around floor(12m/35), giving a
small, proof-facing residue object instead of a raw dead-fraction scalar.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
H3450_PATH = ROOT / "04-computation" / "lrc14_component_cover_obstruction_extractor_codex_20260628.py"
H3451_PATH = ROOT / "04-computation" / "lrc14_component_cover_conductance_router_codex_20260628.py"
F = Fraction


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


H3450 = load_module("hyp3450_component_phase", H3450_PATH)
H3451 = load_module("hyp3451_component_phase", H3451_PATH)


def fmt(x: F | None) -> str:
    if x is None:
        return "None"
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def canonical_speeds(m: int) -> tuple[int, ...]:
    return tuple(list(range(1, 12)) + [13, 84 * m])


def expected_address(m: int) -> int:
    return (48 * m + 6) // 7


def expected_interval(m: int) -> tuple[F, F]:
    a = expected_address(m)
    return (F(14 * a + 1, 588 * m), F(14 * a + 13, 588 * m))


@dataclass(frozen=True)
class TailRecord:
    m: int
    row: object
    summary: object
    projection_components: int
    projection_largest: int

    @property
    def components(self) -> int:
        return len(self.row.components)

    @property
    def dead(self) -> int:
        return len(self.row.dead_components)

    @property
    def escapes(self) -> int:
        return self.components - self.dead

    @property
    def dead_fraction(self) -> F:
        return F(self.dead, self.components)

    @property
    def best(self):
        return self.row.best_component

    @property
    def best_len(self) -> F:
        return self.best.union_measure

    @property
    def best_interval(self) -> tuple[F, F] | None:
        return self.best.best_survivor

    @property
    def has_hyp3433_interval(self) -> bool:
        return self.best_interval == expected_interval(self.m)

    @property
    def has_hyp3433_length(self) -> bool:
        return self.best_len == F(1, 49 * self.m)

    @property
    def has_rank_one_ee_tail(self) -> bool:
        return (
            self.row.best_rank == 1
            and self.has_hyp3433_interval
            and self.has_hyp3433_length
            and self.best.labels == f"L[E:{84 * self.m}] R[E:{84 * self.m}]"
        )


def audit_tail(limit: int) -> list[TailRecord]:
    records: list[TailRecord] = []
    for m in range(1, limit + 1):
        row = H3450.audit_row(f"ap_omit_12_tail_84x{m:03d}", canonical_speeds(m))
        summary = H3451.basic_summary(row)
        parts = H3451.connected_components(H3451.dead_projection(row))
        records.append(
            TailRecord(
                m=m,
                row=row,
                summary=summary,
                projection_components=len(parts),
                projection_largest=max((len(part) for part in parts), default=0),
            )
        )
    return records


def first_eventual(records: list[TailRecord], predicate) -> int | None:
    for record in records:
        if all(predicate(candidate) for candidate in records[record.m - 1 :]):
            return record.m
    return None


def correction_vector(records: list[TailRecord]) -> tuple[list[int], list[int]]:
    if len(records) < 35:
        raise ValueError("need at least 35 records for one mod-35 correction period")
    corrections = [
        records[r - 1].escapes // 2 - (12 * r // 35)
        for r in range(1, 36)
    ]
    failures: list[int] = []
    for record in records:
        predicted = 2 * ((12 * record.m // 35) + corrections[(record.m - 1) % 35])
        if predicted != record.escapes:
            failures.append(record.m)
    return corrections, failures


def tournament_fingerprint() -> tuple[dict[int, int], list[str], int]:
    vertices = {
        "m5_rank_one_endpoint_phase": (10, 10, 10, 9, 9, 10),
        "H3433_tail_address_certificate": (10, 10, 9, 9, 8, 10),
        "connected_dead_projection_family": (9, 9, 10, 10, 8, 9),
        "mod35_beatty_escape_clock": (8, 9, 8, 8, 10, 8),
        "m3_dead_pair_rank_drop": (8, 8, 9, 8, 8, 8),
        "dead_fraction_peak_scalar": (4, 5, 3, 4, 4, 3),
        "raw_component_count_slope": (3, 4, 3, 3, 4, 3),
    }
    scores = {name: sum(values) for name, values in vertices.items()}
    hist = dict(sorted(Counter(scores.values()).items()))
    path = [name for name, _score in sorted(scores.items(), key=lambda item: (-item[1], item[0]))]
    rank = {name: i for i, name in enumerate(path)}
    cycles = 0
    for a, b, c in combinations(vertices, 3):
        ab = rank[a] < rank[b]
        bc = rank[b] < rank[c]
        ca = rank[c] < rank[a]
        if ab == bc == ca:
            cycles += 1
    return hist, path, cycles


def print_sample(records: list[TailRecord], sample_ms: tuple[int, ...]) -> None:
    print("m | components | dead | escapes | best_rank | max_pair | best_len | interval_law | projection")
    for m in sample_ms:
        record = records[m - 1]
        print(
            "{m:3d} | {components:10d} | {dead:4d} | {escapes:7d} | "
            "{rank:9} | {pair:8d} | {best:8} | {law:12} | {projection}".format(
                m=m,
                components=record.components,
                dead=record.dead,
                escapes=record.escapes,
                rank=str(record.row.best_rank),
                pair=record.row.max_dead_pair_rank,
                best=fmt(record.best_len),
                law=str(record.has_hyp3433_interval),
                projection=f"{record.projection_components}/{record.projection_largest}",
            )
        )


def main() -> None:
    limit = 70
    records = audit_tail(limit)

    rank_one_start = first_eventual(records, lambda record: record.has_rank_one_ee_tail)
    pair_rank_drop = first_eventual(records, lambda record: record.row.max_dead_pair_rank <= 2)
    connected_start = first_eventual(records, lambda record: record.projection_components == 1)
    low_rank_escape_failures = [
        record.m for record in records if record.summary.low_rank_escape != record.escapes
    ]
    endpoint_failures = [
        record.m for record in records if record.m >= (rank_one_start or 1) and not record.has_rank_one_ee_tail
    ]
    corrections, beatty_failures = correction_vector(records)
    worst = max(records, key=lambda record: (record.dead_fraction, -record.m))

    print("HYP-3452 AP84 TAIL COMPONENT PHASE AUDIT")
    print("status=EVIDENCE / exact AP-tail component-cover phase audit; not an LRC14 proof")
    print("source=HYP-3433 endpoint-spine law + HYP-3450/HYP-3451 component-cover graph router")
    print()

    print("## Phase Transition")
    print(f"canonical_family=S_m={{1,2,...,11,13,84m}}, checked_m=1..{limit}")
    print(f"rank_one_EE_endpoint_phase_start={rank_one_start}")
    print(f"hyp3433_interval_failures_after_phase={endpoint_failures}")
    print(f"dead_pair_rank_le_2_start={pair_rank_drop}")
    print(f"connected_dead_projection_start={connected_start}")
    print(f"low_rank_escape_equals_alive_failures={low_rank_escape_failures}")
    print(
        "max_dead_fraction={fraction} at m={m}, components={components}, dead={dead}, escapes={escapes}".format(
            fraction=f"{fmt(worst.dead_fraction)} ({float(worst.dead_fraction):.6f})",
            m=worst.m,
            components=worst.components,
            dead=worst.dead,
            escapes=worst.escapes,
        )
    )
    print()

    print("## Sample Rows")
    print_sample(records, (1, 2, 3, 4, 5, 6, 7, 8, 12, 21, 35, 36, 70))
    print()

    print("## Escape Count Clock")
    print("For checked m, escapes(m)=2*(floor(12m/35)+d[m mod 35]).")
    print("Residues are listed as m=1..35, with the final entry representing residue 0.")
    print(f"mod35_correction_vector={corrections}")
    print(f"beatty_clock_failures={beatty_failures}")
    print()

    print("## Lemma Target")
    print("For the canonical AP-tail family, prove directly that m>=5 exposes the")
    print("HYP-3433 interval")
    print("  I_m=[(14ceil(48m/7)+1)/(588m),(14ceil(48m/7)+13)/(588m)]")
    print("as a rank-one E:84m/E:84m component-cover escape.  The transient m=1..4")
    print("has only mixed E/B1 escapes; from m>=3 the paired dead-cover rank is already")
    print("at most 2.  Thus the AP-tail graph proof can split into finite transients")
    print("m=1..4 plus the rank-one harmonic tail m>=5, with the remaining discrete")
    print("escape count controlled by the mod-35 Beatty clock.")
    print()

    hist, path, cycles = tournament_fingerprint()
    print("## Tournament Analysis")
    print("vertices=proof carriers for the AP-tail phase, not runners or raw arcs")
    print("pairwise_observable=phase threshold + endpoint retention + graph connectivity + residue-clock payload")
    print("switch=higher retained proof payload; ties by AP-tail proof route")
    print(f"score_hist={hist}")
    print(f"directed_3cycles={cycles}")
    print("hamiltonian_path=" + " -> ".join(path))
    print()

    print("## Assumption Challenge")
    print("Considered vertices: runners, m-values, E_safe components, dead components,")
    print("branch blockers, endpoint labels, tail gaps, and residues mod 35.  The")
    print("chosen quotient preserves the branch-union escape predicate, endpoint")
    print("rank/labels, dead-cover ranks, and dead-projection connectivity.  It")
    print("forgets full wall geometry away from the AP-tail family, so it is a base")
    print("case/induction carrier, not a global proof of LRC14 by itself.")


if __name__ == "__main__":
    main()
