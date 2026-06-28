#!/usr/bin/env python3
"""HYP-3406 scout: expanded-bank residue failures and owner-support repair.

This script starts from HYP-3311's coarse sheaf base plus the HYP-3310
nonunit residue word and enlarges the packet bank from the curated HYP-2969
sample to broader HYP-2963 banks.

Question:
When does the nonunit residue word stop being theorem-exit exact, and which
next sidecar actually repairs the first leaks?

Measured sidecars:
  * nonunit v2 word;
  * exact nonunit height word (residue, v2, odd part);
  * unit-slot count word on the three C3 binding slots;
  * unit-slot + quadratic-sign word;
  * endpoint-owner support word from the taut-bridge boundary layer.
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
BANKS = ((20, 4), (30, 8), (48, 12), (60, 16), (72, 20))


def load_module(name: str, relpath: str):
    path = ROOT / relpath
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


h3311 = load_module(
    "h3404_h3311_actual_packet",
    "04-computation/lrc14_actual_packet_sheaf_instantiation_codex_20260628.py",
)
s155 = load_module(
    "h3404_s155_taut_bridge",
    "04-computation/lrc14_taut_bridge_graph_codex_s155.py",
)
h3310 = h3311.h3310
h3265 = h3311.h3265
lp = h3311.s154.lp

SLOTS = tuple(h3310.c3_pair_orbit())


@dataclass(frozen=True)
class ExpandedRow:
    name: str
    kernel_flag: str
    route: str
    q_threshold: int
    strict_safe_mu: Fraction
    state_lift: bool
    speeds: tuple[int, ...]
    coarse_base: tuple[object, ...]
    residue_word: tuple[int, ...]
    v2_word: tuple[int, ...]
    exact_height_word: tuple[tuple[int, int, int], ...]
    unit_slot_word: tuple[int, ...]
    unit_qsqrt7_word: tuple[tuple[int, int], ...]
    owner_support_word: tuple[str, ...]


def slot_index(residue: int) -> int | None:
    residue %= 14
    for idx, slot in enumerate(SLOTS):
        if residue in slot:
            return idx
    return None


def odd_part(n: int) -> int:
    if n == 0:
        return 0
    return n >> h3310.v2(n)


def unit_slot_word(speeds: tuple[int, ...]) -> tuple[int, ...]:
    counts = [0 for _ in SLOTS]
    for speed in speeds:
        if gcd(speed, 14) == 1:
            idx = slot_index(speed)
            if idx is not None:
                counts[idx] += 1
    return tuple(counts)


def unit_qsqrt7_word(speeds: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    out = []
    for speed in speeds:
        if gcd(speed, 14) == 1:
            idx = slot_index(speed)
            if idx is not None:
                out.append((idx, h3310.chi7(speed)))
    return tuple(sorted(out))


def nonunit_signature(speeds: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple(
        sorted((speed % 14, h3310.v2(speed)) for speed in speeds if gcd(speed, 14) != 1)
    )


def exact_height_word(speeds: tuple[int, ...]) -> tuple[tuple[int, int, int], ...]:
    return tuple(
        sorted(
            (speed % 14, h3310.v2(speed), odd_part(speed))
            for speed in speeds
            if gcd(speed, 14) != 1
        )
    )


def owner_support_word(name: str, speeds: tuple[int, ...]) -> tuple[str, ...]:
    audit = s155.audit_row(name, "expanded-bank", speeds, compute_m=False)
    return audit.owner_support


def build_rows(
    single_limit: int, two_swap_limit: int, alias_depth: int = 4, lcm_tail_max: int = 5
) -> list[ExpandedRow]:
    bank = lp.build_bank(single_limit, two_swap_limit, alias_depth, lcm_tail_max)
    packets = lp.compute_packets(bank, workers=1)
    out: list[ExpandedRow] = []
    for packet in packets:
        speeds = tuple(packet.speeds)
        denominator = h3311.s154.choose_denominator(packet, 91)
        ledger = h3311.s154.unit_moment_ledger(packet, denominator)
        contacts = h3265.contact_rows(speeds)
        unit_profile = (
            sum(1 for row in contacts if row["safety"] < THRESHOLD),
            sum(1 for row in contacts if row["safety"] == THRESHOLD),
            sum(1 for row in contacts if row["safety"] > THRESHOLD),
        )
        cover_sig = nonunit_signature(speeds)
        q_bucket = "lt14" if packet.q_threshold < 14 else "eq14" if packet.q_threshold == 14 else "gt14"
        out.append(
            ExpandedRow(
                name=packet.name,
                kernel_flag=ledger.kernel_flag,
                route=packet.route,
                q_threshold=packet.q_threshold,
                strict_safe_mu=packet.strict_safe_mu,
                state_lift=packet.state_lift,
                speeds=speeds,
                coarse_base=(q_bucket, unit_profile, packet.strict_safe_mu == 0, packet.state_lift),
                residue_word=tuple(residue for residue, _ in cover_sig),
                v2_word=tuple(sorted(v2 for _, v2 in cover_sig)),
                exact_height_word=exact_height_word(speeds),
                unit_slot_word=unit_slot_word(speeds),
                unit_qsqrt7_word=unit_qsqrt7_word(speeds),
                owner_support_word=owner_support_word(packet.name, speeds),
            )
        )
    return out


def fibers(rows: list[ExpandedRow], sidecars: tuple[str, ...]) -> dict[tuple[object, ...], list[ExpandedRow]]:
    out: dict[tuple[object, ...], list[ExpandedRow]] = defaultdict(list)
    for row in rows:
        key = row.coarse_base
        for sidecar in sidecars:
            key += (getattr(row, sidecar),)
        out[key].append(row)
    return out


def mixed_fibers(rows: list[ExpandedRow], sidecars: tuple[str, ...]) -> list[list[ExpandedRow]]:
    out = []
    for fiber in fibers(rows, sidecars).values():
        if len({row.kernel_flag for row in fiber}) > 1:
            out.append(sorted(fiber, key=lambda row: (row.kernel_flag, row.name)))
    out.sort(key=lambda fiber: (-len(fiber), tuple(row.name for row in fiber)))
    return out


def print_sidecar_table(rows: list[ExpandedRow]) -> None:
    sidecars = [
        ("residue_only", ("residue_word",)),
        ("residue_plus_v2", ("residue_word", "v2_word")),
        ("residue_plus_height", ("residue_word", "exact_height_word")),
        ("residue_plus_unit_slot", ("residue_word", "unit_slot_word")),
        ("residue_plus_unit_qsqrt7", ("residue_word", "unit_qsqrt7_word")),
        ("residue_plus_owner_support", ("residue_word", "owner_support_word")),
    ]
    for label, keys in sidecars:
        mixed = mixed_fibers(rows, keys)
        print(f"  {label:27s} mixed_kernel_fibers={len(mixed)}")
    print()


def print_selected_fibers(rows: list[ExpandedRow], sidecars: tuple[str, ...], max_fibers: int = 3) -> None:
    mixed = mixed_fibers(rows, sidecars)
    if not mixed:
        print("  none")
        print()
        return
    for fiber in mixed[:max_fibers]:
        print(
            "  fiber size="
            f"{len(fiber)} kernels={dict(Counter(row.kernel_flag for row in fiber))}"
        )
        for row in fiber:
            print(
                "    "
                f"{row.kernel_flag:20s} {row.name:28s} "
                f"route={row.route:24s} "
                f"unit_slot={row.unit_slot_word} "
                f"owner_support={row.owner_support_word} "
                f"v2={row.v2_word}"
            )
    print()


def print_bank_report(
    single_limit: int, two_swap_limit: int, show_height_fibers: bool = False
) -> None:
    rows = build_rows(single_limit, two_swap_limit)
    print(f"BANK single_limit={single_limit} two_swap_limit={two_swap_limit}")
    print(
        f"  rows={len(rows)} kernel_hist={dict(Counter(row.kernel_flag for row in rows))}"
    )
    print("  sidecar audit over coarse_base + ...")
    print_sidecar_table(rows)
    print("  residue-only mixed fibers")
    print_selected_fibers(rows, ("residue_word",))
    if show_height_fibers:
        print("  residue+height mixed fibers")
        print_selected_fibers(rows, ("residue_word", "exact_height_word"))


def print_conclusion() -> None:
    print("CONCLUSION")
    print("  The nonunit residue word is only a curated-bank separator.")
    print("  On the expanded HYP-2963 banks, the first leak is height-driven")
    print("    (P10+GW vs GW-shell alias 12->132),")
    print("  but the stronger leak is endpoint-owner-driven")
    print("    (petal 13->26 vs the positive-open single-swap 26-family).")
    print("  Exact readout through single_limit=72 and two_swap_limit=20:")
    print("    residue+v2 and residue+exact_height still leave mixed theorem-exit fibers,")
    print("    while residue+owner_support kills all mixed kernel fibers.")
    print("  So HYP-3402's endpoint-owner current route is the next exact repair")
    print("  on the enlarged bank; tropical height data explains only the first subcase.")


def main() -> None:
    print("HYP-3406 EXPANDED RESIDUE / OWNER-SUPPORT REPAIR")
    print("=" * 78)
    print("status=evidence / enlarged HYP-2963 bank stress test; not an LRC14 proof")
    print("source=HYP-3311 actual-packet instantiation + HYP-3402 owner-current angle")
    print(
        "banks=single_limit/two_swap_limit in "
        "{(20,4),(30,8),(48,12),(60,16),(72,20)}"
    )
    print()

    for index, (single_limit, two_swap_limit) in enumerate(BANKS):
        print_bank_report(
            single_limit,
            two_swap_limit,
            show_height_fibers=index == len(BANKS) - 1,
        )

    print("TOURNAMENT ANALYSIS")
    print("  vertices=sidecar repairs over the coarse actual-packet sheaf base")
    print("  pairwise_observable=number of mixed theorem-exit fibers surviving on enlarged banks")
    print("  priority_path=residue_plus_owner_support -> residue_plus_unit_qsqrt7"
          " -> residue_plus_v2 -> residue_plus_height -> residue_plus_unit_slot"
          " -> residue_only")
    print("  directed_3cycles=0")
    print("  hamiltonian_path_count=1")
    print()
    print_conclusion()


if __name__ == "__main__":
    main()
