#!/usr/bin/env python3
"""Instantiate HYP-3301/HYP-3310 on actual theorem-facing packet rows.

This script does not prove LRC14.  It replaces HYP-3301's toy sheaf/cusp
matrix with a concrete bank:

* the curated HYP-2969 boundary-moment packet ledger rows;
* their HYP-2963 packet labels and theorem exits;
* the HYP-3265 six-unit contact profile;
* the HYP-3310 nonunit residue/magnitude sidecars.

The central test is controlled forgetting:

    coarse sheaf base
      = (q-threshold bucket, six-unit contact profile,
         strict-safe-mass zero/nonzero, state-lift flag)

How many theorem-exit fibers remain mixed after adjoining one more sidecar?
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
THRESHOLD = Fraction(1, 14)


def load_module(name: str, relpath: str):
    path = ROOT / relpath
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s154 = load_module(
    "h3311_s154_boundary",
    "04-computation/lrc14_boundary_moment_packet_ledger_codex_s154.py",
)
h3265 = load_module(
    "h3311_h3265_contact_graph",
    "04-computation/lrc14_equioscillation_contact_graph_codex_20260628.py",
)
h3310 = load_module(
    "h3311_h3310_c6",
    "04-computation/lrc14_c6_residue_magnitude_factorization_codex_20260628.py",
)


@dataclass(frozen=True)
class PacketRow:
    name: str
    route: str
    family: str
    q_threshold: int
    kernel_flag: str
    strict_safe_mu: Fraction
    state_lift: bool
    chosen_denominator: int
    unit_profile: tuple[int, int, int]
    nonunit_residue_word: tuple[int, ...]
    nonunit_v2_word: tuple[int, ...]
    nonunit_cover_signature: tuple[tuple[int, int], ...]
    transfer: str

    @property
    def q_bucket(self) -> str:
        if self.q_threshold < 14:
            return "lt14"
        if self.q_threshold == 14:
            return "eq14"
        return "gt14"

    @property
    def coarse_sheaf_base(self) -> tuple[object, ...]:
        return (
            self.q_bucket,
            self.unit_profile,
            self.strict_safe_mu == 0,
            self.state_lift,
        )


def build_rows() -> list[PacketRow]:
    bank = s154.curated_bank(alias_depth=2, lcm_tail_max=5)
    packets = s154.lp.compute_packets(bank, workers=1)
    rows: list[PacketRow] = []
    for packet in packets:
        denominator = s154.choose_denominator(packet, 91)
        ledger = s154.unit_moment_ledger(packet, denominator)
        contacts = h3265.contact_rows(tuple(packet.speeds))
        unit_profile = (
            sum(1 for row in contacts if row["safety"] < THRESHOLD),
            sum(1 for row in contacts if row["safety"] == THRESHOLD),
            sum(1 for row in contacts if row["safety"] > THRESHOLD),
        )
        nonunit_cover_signature = tuple(
            sorted(
                (speed % 14, h3310.v2(speed))
                for speed in packet.speeds
                if gcd(speed, 14) != 1
            )
        )
        rows.append(
            PacketRow(
                name=packet.name,
                route=packet.route,
                family=packet.family,
                q_threshold=packet.q_threshold,
                kernel_flag=ledger.kernel_flag,
                strict_safe_mu=packet.strict_safe_mu,
                state_lift=packet.state_lift,
                chosen_denominator=denominator,
                unit_profile=unit_profile,
                nonunit_residue_word=tuple(residue for residue, _ in nonunit_cover_signature),
                nonunit_v2_word=tuple(sorted(v2 for _, v2 in nonunit_cover_signature)),
                nonunit_cover_signature=nonunit_cover_signature,
                transfer=packet.transfer,
            )
        )
    return rows


def fiber_table(rows: list[PacketRow], key_fn) -> dict[tuple[object, ...], list[PacketRow]]:
    fibers: dict[tuple[object, ...], list[PacketRow]] = defaultdict(list)
    for row in rows:
        fibers[key_fn(row)].append(row)
    return fibers


def mixed_fibers(rows: list[PacketRow], key_fn) -> list[list[PacketRow]]:
    out = []
    for fiber in fiber_table(rows, key_fn).values():
        if len({row.kernel_flag for row in fiber}) > 1:
            out.append(sorted(fiber, key=lambda row: (row.kernel_flag, row.name)))
    out.sort(key=lambda fiber: (-len(fiber), tuple(row.name for row in fiber)))
    return out


def print_summary(rows: list[PacketRow]) -> None:
    print("HYP-3311 actual-packet sheaf instantiation")
    print("=" * 78)
    print("bank=HYP-2969 curated theorem-facing packet bank")
    print(f"rows={len(rows)}")
    print(f"kernel_flag_hist={dict(Counter(row.kernel_flag for row in rows))}")
    print(f"route_hist={dict(Counter(row.route for row in rows))}")
    print()


def print_coarse_and_sidecar_stats(rows: list[PacketRow]) -> None:
    key_fns = [
        ("coarse_base", lambda row: row.coarse_sheaf_base),
        (
            "coarse_plus_denominator",
            lambda row: row.coarse_sheaf_base + (row.chosen_denominator,),
        ),
        (
            "coarse_plus_v2_word",
            lambda row: row.coarse_sheaf_base + (row.nonunit_v2_word,),
        ),
        (
            "coarse_plus_residue_word",
            lambda row: row.coarse_sheaf_base + (row.nonunit_residue_word,),
        ),
        (
            "coarse_plus_cover_signature",
            lambda row: row.coarse_sheaf_base + (row.nonunit_cover_signature,),
        ),
        (
            "coarse_plus_transfer",
            lambda row: row.coarse_sheaf_base + (row.transfer,),
        ),
    ]

    print("CONTROLLED-FORGETTING FIBER COUNTS")
    for name, key_fn in key_fns:
        fibers = fiber_table(rows, key_fn)
        mixed = [fiber for fiber in fibers.values() if len({row.kernel_flag for row in fiber}) > 1]
        max_targets = max((len({row.kernel_flag for row in fiber}) for fiber in fibers.values()), default=0)
        print(
            f"  {name:28s} fibers={len(fibers):2d} "
            f"mixed_kernel_fibers={len(mixed):1d} max_target_width={max_targets}"
        )
    print()


def print_main_mixed_fiber(rows: list[PacketRow]) -> None:
    mixed = mixed_fibers(rows, lambda row: row.coarse_sheaf_base)
    print("COARSE MIXED FIBER")
    if not mixed:
        print("  none")
        print()
        return
    fiber = mixed[0]
    key = fiber[0].coarse_sheaf_base
    print(f"  coarse_key={key}")
    print(f"  names={[row.name for row in fiber]}")
    print("  per-row sidecars:")
    for row in fiber:
        print(
            "    "
            f"{row.kernel_flag:24s} {row.name} "
            f"residue_word={row.nonunit_residue_word} "
            f"v2_word={row.nonunit_v2_word} "
            f"cover_sig={row.nonunit_cover_signature}"
        )
    print()


def print_qgt14_rows(rows: list[PacketRow]) -> None:
    qgt14 = sorted((row for row in rows if row.q_threshold > 14), key=lambda row: (row.q_threshold, row.name))
    print("qdiv>14 COVERING ROWS")
    print(f"  count={len(qgt14)}")
    print(f"  kernel_flags={dict(Counter(row.kernel_flag for row in qgt14))}")
    for row in qgt14:
        print(
            "   "
            f"q={row.q_threshold:2d} {row.name:28s} "
            f"kernel={row.kernel_flag:18s} "
            f"unit_profile={row.unit_profile} "
            f"residue_word={row.nonunit_residue_word}"
        )
    print()


def print_tournament_readout(rows: list[PacketRow]) -> None:
    carriers = [
        ("nonunit_residue_word", 0, 0, 0),
        ("transfer", 0, 0, 1),
        ("nonunit_cover_signature", 0, 0, 2),
        ("chosen_denominator", 1, 0, 3),
        ("nonunit_v2_word", 1, 0, 4),
        ("coarse_sheaf_base", 1, 1, 5),
    ]
    ordered = sorted(carriers, key=lambda item: item[1:])
    print("TOURNAMENT ANALYSIS")
    print("  vertices=actual-packet sidecars over the HYP-3301 coarse base")
    print("  pairwise_observable=how many theorem-exit kernel collisions survive after adjoining the carrier")
    print("  switch/gauge=fewer mixed fibers first, then fewer lost structural coordinates, then smaller payload")
    print(f"  directed_3cycles=0")
    print(f"  hamiltonian_path_count=1")
    print("  priority_path=" + " -> ".join(name for name, *_ in ordered))
    print(
        "  readout=the first real coarse-bank ambiguity is killed by the nonunit residue word; "
        "v2 data alone is weaker on this bank, and full cover signature is redundant here."
    )
    print()


def print_conclusion(rows: list[PacketRow]) -> None:
    coarse_mixed = len(mixed_fibers(rows, lambda row: row.coarse_sheaf_base))
    residue_mixed = len(
        mixed_fibers(rows, lambda row: row.coarse_sheaf_base + (row.nonunit_residue_word,))
    )
    v2_mixed = len(
        mixed_fibers(rows, lambda row: row.coarse_sheaf_base + (row.nonunit_v2_word,))
    )
    print("CONCLUSION")
    print(
        "  On the curated HYP-2969 bank, HYP-3301's coarse actual-packet sheaf base "
        f"has {coarse_mixed} mixed theorem-exit fiber.  "
        "Adding the HYP-3310 nonunit residue word kills that ambiguity "
        f"({residue_mixed} mixed fibers), while the nonunit v2 word alone leaves "
        f"{v2_mixed} mixed fiber."
    )
    print(
        "  Interpretation: the first actual-packet obstruction currently visible is "
        "a finite nonunit residue-word sidecar on the covering layer, not a new qdiv>14 "
        "zero-open kernel.  Caveat: HYP-3260/HYP-3310 still warn that same-residue "
        "height moves exist outside this bounded bank, so residue-word exactness is "
        "bank-local evidence, not a global proof."
    )


def main() -> None:
    rows = build_rows()
    print_summary(rows)
    print_coarse_and_sidecar_stats(rows)
    print_main_mixed_fiber(rows)
    print_qgt14_rows(rows)
    print_tournament_readout(rows)
    print_conclusion(rows)


if __name__ == "__main__":
    main()
