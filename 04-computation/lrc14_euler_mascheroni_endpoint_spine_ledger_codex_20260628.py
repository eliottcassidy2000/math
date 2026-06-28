#!/usr/bin/env python3
"""HYP-3433: Euler-Mascheroni finite parts for the LRC14 endpoint spine.

HYP-3429 found that the two-adic relocation survivor can often be certified by
one or two endpoint walls.  HYP-3431 then closed the canonical 84m tower by a
corridor-fence argument, while HYP-3432 showed that harmonic wall budgets are
sidecars rather than certificates.  The user's Euler-Mascheroni prompt suggests
a complementary test: if an infinite packet family emits endpoint windows whose
lengths are harmonic in the family parameter, then gamma should enter only as
the renormalized finite part after the logarithmic tail is removed.

The canonical covering tower

    S_m = {1,2,...,11,13,84m}

is an ideal stress test.  This script reuses the exact HYP-3429 interval
machinery and audits the endpoint spine for m <= N.  It proves no theorem, but
it does produce an exact certificate law on the checked tail:

    best endpoint-spine length = 1/(49m), labelled only by E:84m.

Consequently the cumulative best-spine capacity has finite part

    (gamma - H_4) / 49

after subtracting log(M)/49.  The proof-facing message is a guardrail:
Euler-Mascheroni is useful as a labelled tail normalizer, not as an unlabeled
scalar replacement for endpoint, branch, owner, or component data.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib import util
from itertools import combinations
from math import log
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
H3429_PATH = ROOT / "04-computation" / "lrc14_component_spine_certificate_codex_20260628.py"
F = Fraction
EULER_GAMMA = 0.577215664901532860606512090082402431


def load_h3429():
    spec = util.spec_from_file_location("h3429_endpoint_spine_tools", H3429_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {H3429_PATH}")
    module = util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


h3429 = load_h3429()


def fmt(x: F | None) -> str:
    if x is None:
        return "n/a"
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def harmonic(n: int) -> F:
    return sum((F(1, m) for m in range(1, n + 1)), F(0))


def canonical_speeds(m: int) -> tuple[int, ...]:
    return tuple(list(range(1, 12)) + [13, 84 * m])


@dataclass(frozen=True)
class CanonicalRecord:
    m: int
    row: object

    @property
    def expected_best_length(self) -> F:
        return F(1, 49 * self.m)

    @property
    def expected_label(self) -> tuple[tuple[str, int], ...]:
        return (("E", 84 * self.m),)

    @property
    def expected_address(self) -> int:
        return (48 * self.m + 6) // 7

    @property
    def expected_interval(self) -> tuple[F, F]:
        address = self.expected_address
        return (F(14 * address + 1, 588 * self.m), F(14 * address + 13, 588 * self.m))

    @property
    def has_tail_law(self) -> bool:
        return (
            self.row.best.length == self.expected_best_length
            and self.row.best.labels == self.expected_label
        )

    @property
    def has_address_law(self) -> bool:
        return self.row.best.interval == self.expected_interval

    @property
    def endpoint_weight(self) -> F:
        return F(1, 84 * self.m)


def audit_canonical(limit: int) -> list[CanonicalRecord]:
    return [
        CanonicalRecord(
            m,
            h3429.audit_row(f"canonical_84m_gamma_{m:03d}", canonical_speeds(m)),
        )
        for m in range(1, limit + 1)
    ]


def first_eventual_law(records: list[CanonicalRecord]) -> int | None:
    for start in range(1, len(records) + 1):
        if all(record.has_tail_law for record in records[start - 1 :]):
            return start
    return None


def branch_blocks(records: list[CanonicalRecord]) -> list[tuple[int, int, tuple[int, ...]]]:
    blocks: list[tuple[int, int, tuple[int, ...]]] = []
    start = records[0].m
    active = records[0].row.best.branches
    for record in records[1:]:
        if record.row.best.branches != active:
            blocks.append((start, record.m - 1, active))
            start = record.m
            active = record.row.best.branches
    blocks.append((start, records[-1].m, active))
    return blocks


def label_text(labels: tuple[tuple[str, int], ...]) -> str:
    return ",".join(f"{kind}:{value}" for kind, value in labels) if labels else "none"


AXES = (
    "predicate_retention",
    "endpoint_exactness",
    "tail_compression",
    "branch_sheet_retention",
    "owner_sidecar_compatibility",
    "finite_lemma_usefulness",
    "scalar_risk_control",
)


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: tuple[int, ...]

    @property
    def total(self) -> int:
        return sum(self.scores)


CARRIERS = (
    Carrier("eventual_E_endpoint_harmonic_law", (10, 10, 10, 8, 7, 9, 10)),
    Carrier("gamma_finite_part_sidecar", (8, 8, 10, 6, 7, 8, 10)),
    Carrier("H3429_endpoint_spine_certificate", (10, 10, 8, 9, 8, 10, 9)),
    Carrier("Mertens_loglog_tail_normalizer", (7, 7, 9, 6, 8, 7, 9)),
    Carrier("endpoint_weight_exchange_rate", (9, 9, 8, 7, 7, 8, 9)),
    Carrier("total_good_measure_density_probe", (8, 6, 5, 7, 6, 5, 6)),
    Carrier("raw_window_count_growth", (5, 5, 4, 4, 4, 4, 4)),
    Carrier("raw_gamma_constant_scalar", (3, 2, 7, 1, 1, 2, 1)),
)


def tournament() -> tuple[dict[int, int], list[str], int]:
    hist = dict(sorted(Counter(carrier.total for carrier in CARRIERS).items()))
    ordered = sorted(CARRIERS, key=lambda carrier: (-carrier.total, CARRIERS.index(carrier)))
    rank = {carrier.name: index for index, carrier in enumerate(ordered)}
    cycles = 0
    for a, b, c in combinations(CARRIERS, 3):
        ab = rank[a.name] < rank[b.name]
        bc = rank[b.name] < rank[c.name]
        ca = rank[c.name] < rank[a.name]
        if ab == bc == ca:
            cycles += 1
    return hist, [carrier.name for carrier in ordered], cycles


def main() -> None:
    limit = 180
    records = audit_canonical(limit)
    law_start = first_eventual_law(records)
    if law_start is None:
        raise AssertionError("no eventual harmonic law found in checked range")

    tail = [record for record in records if record.m >= law_start]
    failures = [record.m for record in tail if not record.has_tail_law]
    address_failures = [record.m for record in tail if not record.has_address_law]
    best_tail_sum = sum((record.row.best.length for record in tail), F(0))
    expected_tail_sum = (harmonic(limit) - harmonic(law_start - 1)) / 49
    finite_part = float(best_tail_sum) - log(limit) / 49
    expected_finite_part = (EULER_GAMMA - float(harmonic(law_start - 1))) / 49
    endpoint_weight_sum = sum((record.endpoint_weight for record in tail), F(0))
    exchange_failures = [
        record.m
        for record in tail
        if record.row.best.length != F(12, 7) * record.endpoint_weight
    ]

    total_measures = [record.row.total_good_measure for record in records]
    tail_total_measures = [record.row.total_good_measure for record in tail]
    sample_ms = (1, 2, 3, 4, 5, 6, 12, 24, 48, 96, 144, 180)

    print("HYP-3433 EULER-MASCHERONI ENDPOINT-SPINE LEDGER")
    print("=" * 72)
    print("Identity:")
    print("  Gamma is tested as a finite-part tail normalizer, not as a raw")
    print("  LRC floor scalar.  The packet remains labelled by endpoint walls.")
    print()

    print("A. Canonical 84m exact audit")
    print(f"  rows audited:                         {len(records)}")
    print(f"  canonical family:                     {{1,...,11,13,84m}}")
    print(f"  eventual law begins at m:             {law_start}")
    print(f"  failures on checked tail:             {failures}")
    print(f"  tail label law:                       best label = E:84m")
    print(f"  tail length law:                      best_len = 1/(49m)")
    print(f"  tail address law:                     a_m=ceil(48m/7)")
    print("                                        best=[(14a_m+1)/(588m),(14a_m+13)/(588m)]")
    print(f"  tail address failures:                {address_failures}")
    print(f"  exact scaled law residual:            {fmt(49 * best_tail_sum - (harmonic(limit) - harmonic(law_start - 1)))}")
    print(f"  endpoint-weight exchange residuals:   {exchange_failures}")
    print(f"  branch blocks for best spine:         {branch_blocks(records)}")
    print()

    print("B. Euler-Mascheroni finite part")
    harmonic_tail_name = f"(H_{limit}-H_{law_start - 1})"
    print(
        f"  sum_{{m={law_start}}}^{limit} best_len:       "
        f"{harmonic_tail_name}/49 = {float(best_tail_sum):.12f}"
    )
    print(f"  expected harmonic sum residual:       {fmt(best_tail_sum - expected_tail_sum)}")
    print(f"  sum minus log(M)/49 at M={limit}:      {finite_part:.12f}")
    print(f"  limiting finite part (gamma-H_4)/49:  {expected_finite_part:.12f}")
    print(f"  finite-part error at M={limit}:        {finite_part - expected_finite_part:.12e}")
    print(
        f"  sum endpoint weights:                 "
        f"{harmonic_tail_name}/84 = {float(endpoint_weight_sum):.12f}"
    )
    print("  exchange rate:                        best_len = (12/7) * 1/(84m)")
    print()

    print("C. Separate density signal")
    print(f"  total_good_measure range, all rows:   {fmt(min(total_measures))} .. {fmt(max(total_measures))}")
    print(f"  total_good_measure range, tail rows:  {fmt(min(tail_total_measures))} .. {fmt(max(tail_total_measures))}")
    print("  samples: m | windows | total_good_measure | best_len | branch | label")
    for m in sample_ms:
        record = records[m - 1]
        print(
            f"  {m:3d} | {record.row.window_count:7d} | "
            f"{fmt(record.row.total_good_measure):>18} | "
            f"{fmt(record.row.best.length):>8} | "
            f"{record.row.best.branches!s:>6} | {label_text(record.row.best.labels)}"
        )
    print()

    print("D. Proof interpretation")
    print("  The canonical tail separates two clocks.  The survivor-density clock stays")
    print("  bounded and mildly oscillatory, while the selected endpoint-spine clock is")
    print("  exactly harmonic after m=4.  Summing selected spine lengths therefore")
    print("  produces a logarithmic divergence whose finite part is Euler's constant.")
    print("  This does not prove LRC14; it says that any infinite family/packet proof")
    print("  should normalize harmonic endpoint tails only after retaining the label")
    print("  E:84m and the branch/component sidecar.")
    print()

    hist, path, cycles = tournament()
    print("E. Tournament Analysis")
    print("  vertices are proof carriers and tail normalizers, not runners.")
    print(f"  axes={','.join(AXES)}")
    print(f"  score_hist={hist}")
    print(f"  directed_3cycles={cycles}")
    print("  hamiltonian_path=" + " -> ".join(path))
    print()

    print("F. Assumption challenge")
    print("  Considered vertices: runners, the m-parameter, survivor components,")
    print("  endpoint labels, branch masks, harmonic tail atoms, total measures,")
    print("  Mertens/gamma constants, and proof obligations.  The chosen quotient")
    print("  keeps endpoint labels plus the actual survivor window.  It destroys")
    print("  raw runner identities and most component geometry, so it is legal only")
    print("  for tail normalization and not for the final relocation predicate.")


if __name__ == "__main__":
    main()
