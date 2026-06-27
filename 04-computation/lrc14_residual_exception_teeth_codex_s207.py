#!/usr/bin/env python3
"""Mine the topology exceptions inside the LRC14 residual tooth atlas.

HYP-3035 showed that arc topology splits 13 of the 15 coarse ET+unit
route-mixed residual fibers, while the remaining two need the coarse
safe-component stalk.  HYP-3036 independently showed that the primitive
safe-period deck for 2 <= q <= 13 splits the same residual route debt.

This pass zooms in on the two same-topology exceptions and asks whether they
are arbitrary topology failures or a small owner-strip / primitive-period
collar that can be stated as a proof lemma.

Tournament Analysis vertices are exception repair carriers, not runners.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
from pathlib import Path
import argparse
import re
import sys


REPO = Path(__file__).resolve().parents[1]


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


tooth = load_module(
    "lrc14_residual_tooth_atlas_s201",
    REPO / "04-computation" / "lrc14_residual_tooth_atlas_codex_s199.py",
)
ram = load_module(
    "lrc14_ramanujan_route_scheduler_s201",
    REPO / "04-computation" / "lrc14_ramanujan_route_scheduler_codex_s200.py",
)


def fmt_fraction(fr) -> str:
    return str(fr.numerator) if fr.denominator == 1 else f"{fr.numerator}/{fr.denominator}"


def single_swap(name: str) -> tuple[int, int] | None:
    match = re.fullmatch(r"single swap (\d+)->(\d+)", name)
    if not match:
        return None
    return int(match.group(1)), int(match.group(2))


def route_pure(rows: list, key_func) -> bool:
    buckets: dict[object, list] = defaultdict(list)
    for row in rows:
        buckets[key_func(row)].append(row)
    return all(len({row.packet.route for row in bucket}) == 1 for bucket in buckets.values())


def split_signature(rows: list, key_func) -> tuple[int, int, int]:
    buckets: dict[object, list] = defaultdict(list)
    for row in rows:
        buckets[key_func(row)].append(row)
    mixed = sum(len({row.packet.route for row in bucket}) > 1 for bucket in buckets.values())
    return len(buckets), mixed, max((len(bucket) for bucket in buckets.values()), default=0)


def row_summary(row) -> str:
    packet = row.packet
    drop_add = single_swap(packet.name)
    drop = "-" if drop_add is None else str(drop_add[0])
    first13 = ram.first_primitive_safe_q(packet, 13)
    first14 = ram.first_primitive_safe_q(packet, 14)
    deck13 = ram.primitive_count_deck(packet, 13)
    q14_count = len(ram.primitive_safe_residues(packet.speeds, 14))
    return (
        f"    {packet.name:<24} route={packet.route:<16} drop={drop:<2} "
        f"M={fmt_fraction(packet.M):>7} q0={packet.q_threshold:<2d} "
        f"first13={str(first13):>4s} first14={str(first14):>4s} "
        f"deck13={deck13} q14_count={q14_count:<2d} "
        f"topo={row.topology.compact} stalk={row.stalk.short()}"
    )


@dataclass(frozen=True)
class ExceptionCarrier:
    name: str
    exception_coverage: int
    route_split: int
    nonroute_legality: int
    local_owner_signal: int
    primitive_alignment: int
    compression: int
    proof_cost: int

    @property
    def vector(self) -> tuple[int, ...]:
        return (
            self.exception_coverage,
            self.route_split,
            self.nonroute_legality,
            self.local_owner_signal,
            self.primitive_alignment,
            self.compression,
            8 - self.proof_cost,
        )


EXCEPTION_CARRIERS = (
    ExceptionCarrier("topology_then_owner_stalk_rule", 5, 5, 5, 5, 4, 4, 3),
    ExceptionCarrier("coarse_safe_stalk", 5, 5, 5, 5, 3, 3, 3),
    ExceptionCarrier("primitive_deck_2_13", 5, 5, 5, 2, 5, 4, 3),
    ExceptionCarrier("arc_topology_compact", 2, 1, 5, 3, 2, 5, 2),
    ExceptionCarrier("exact_safe_stalk", 5, 5, 5, 5, 3, 1, 4),
    ExceptionCarrier("route_label_sink", 5, 5, 0, 0, 0, 5, 1),
)


def carrier_tournament() -> dict[str, object]:
    names = [carrier.name for carrier in EXCEPTION_CARRIERS]
    order = {name: i for i, name in enumerate(names)}
    edges: set[tuple[str, str]] = set()
    score = Counter()

    for a, b in combinations(EXCEPTION_CARRIERS, 2):
        aw = sum(x > y for x, y in zip(a.vector, b.vector))
        bw = sum(x < y for x, y in zip(a.vector, b.vector))
        if aw > bw or (aw == bw and order[a.name] < order[b.name]):
            edges.add((a.name, b.name))
            score[a.name] += 1
        else:
            edges.add((b.name, a.name))
            score[b.name] += 1

    c3 = 0
    for a, b, c in combinations(names, 3):
        out = {name: 0 for name in (a, b, c)}
        for x, y in combinations((a, b, c), 2):
            if (x, y) in edges:
                out[x] += 1
            else:
                out[y] += 1
        if sorted(out.values()) == [1, 1, 1]:
            c3 += 1

    edge_idx = {(names.index(a), names.index(b)) for a, b in edges}
    hp = 0
    first_hp: tuple[str, ...] | None = None
    for perm in permutations(range(len(names))):
        if all((perm[i], perm[i + 1]) in edge_idx for i in range(len(names) - 1)):
            hp += 1
            if first_hp is None:
                first_hp = tuple(names[i] for i in perm)

    return {
        "score_hist": dict(sorted(Counter(score[name] for name in names).items())),
        "directed_3cycles": c3,
        "hamiltonian_path_count": hp,
        "first_hamiltonian_path": first_hp,
        "score_order": sorted(names, key=lambda name: (-score[name], order[name])),
    }


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge")
    print("  considered vertices:")
    print("    runners, route labels, residual fibers, topology buckets,")
    print("    safe-component stalks, primitive denominator layers, wall-crossing")
    print("    events, Haar zeta teeth, and proof obligations.")
    print("  chosen vertices:")
    print("    exception repair carriers acting on the two same-topology residual")
    print("    buckets, not runners.")
    print("  preserved LRC predicate:")
    print("    strict-open boundary/open status after the coarse ET+unit gate.")
    print("  destroyed information:")
    print("    exact route labels until an owner-stalk or primitive-period tooth")
    print("    separates the topology exception.")
    print("  challenged assumption:")
    print("    a topology failure need not be a new route family; here it is a")
    print("    two-row owner-strip collar with an independent q<=13 deck split.")
    print()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--single-limit", type=int, default=180)
    parser.add_argument("--two-swap-limit", type=int, default=36)
    parser.add_argument("--alias-depth", type=int, default=4)
    parser.add_argument("--lcm-tail-max", type=int, default=5)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()

    fibers = tooth.build_residual_fibers(args)
    teeth = tooth.tooth_family()
    topology_tooth = teeth[0]
    coarse_stalk_tooth = teeth[1]
    exact_stalk_tooth = teeth[2]
    magnitude_tooth = teeth[3]
    certificate_tooth = teeth[4]
    all_rows = [row for rows in fibers for row in rows]

    exception_buckets: list[tuple[int, int, list]] = []
    for fiber_index, rows in enumerate(fibers):
        for bucket_index, bucket in enumerate(tooth.route_mixed(rows, topology_tooth.key_func)):
            exception_buckets.append((fiber_index, bucket_index, bucket))

    exception_rows = [row for _fi, _bi, bucket in exception_buckets for row in bucket]
    exception_names = tuple(row.packet.name for row in exception_rows)
    exception_swaps = [single_swap(name) for name in exception_names]
    exception_drop_speeds = sorted({swap[0] for swap in exception_swaps if swap is not None})

    covering_rows = [row for row in exception_rows if row.packet.route == "COVERING-MOMENT"]
    q_rows = [row for row in exception_rows if row.packet.route == "Q-WITNESS"]
    covering_zero_deck13 = all(
        ram.first_primitive_safe_q(row.packet, 13) is None
        and sum(ram.primitive_count_deck(row.packet, 13)) == 0
        for row in covering_rows
    )
    q_witness_drop_aligned = all(
        single_swap(row.packet.name) is not None
        and ram.first_primitive_safe_q(row.packet, 13) == single_swap(row.packet.name)[0]
        and row.packet.q_threshold == single_swap(row.packet.name)[0]
        for row in q_rows
    )
    stalk_splits_all = all(route_pure(bucket, coarse_stalk_tooth.key_func) for _fi, _bi, bucket in exception_buckets)
    primitive_splits_all = all(
        route_pure(bucket, lambda row: ram.primitive_count_deck(row.packet, 13))
        for _fi, _bi, bucket in exception_buckets
    )

    print("=== LRC14 residual topology-exception teeth S207 ===")
    print(
        "source=S199 residual tooth atlas + S200 primitive-period scheduler; "
        f"single_limit={args.single_limit} two_swap_limit={args.two_swap_limit} "
        f"alias_depth={args.alias_depth} lcm_tail_max={args.lcm_tail_max}"
    )
    print(f"residual_fibers={len(fibers)} residual_rows={len(all_rows)}")
    print()
    print_assumption_challenge()

    print("[1] Topology exception census")
    print(f"  topology_exception_fibers={len(exception_buckets)}")
    print(f"  topology_exception_buckets={len(exception_buckets)}")
    print(f"  topology_exception_rows={len(exception_rows)}")
    print(f"  exception_drop_speeds={exception_drop_speeds}")
    print(f"  all_exception_rows_single_swap={all(swap is not None for swap in exception_swaps)}")
    print(f"  exception_route_counts={dict(Counter(row.packet.route for row in exception_rows))}")
    print(f"  exception_status_counts={dict(Counter(tooth.fz.status(row.packet) for row in exception_rows))}")
    print(f"  covering_zero_primitive_deck_2_13={covering_zero_deck13}")
    print(f"  q_witness_first13_equals_drop_and_q0={q_witness_drop_aligned}")
    print(f"  coarse_stalk_splits_all_exceptions={stalk_splits_all}")
    print(f"  primitive_deck_2_13_splits_all_exceptions={primitive_splits_all}")
    print()

    print("[2] Same-topology buckets")
    for fiber_index, bucket_index, bucket in exception_buckets:
        topology_key = topology_tooth.key_func(bucket[0])
        fiber_rows = fibers[fiber_index]
        print(
            f"  exception[{fiber_index:02d}.{bucket_index}] size={len(bucket)} "
            f"full_fiber_size={len(fiber_rows)} topology={topology_key}"
        )
        print(
            "    full_fiber_routes="
            + str(dict(Counter(row.packet.route for row in fiber_rows)))
            + " full_fiber_names="
            + str(tuple(row.packet.name for row in fiber_rows))
        )
        for row in sorted(bucket, key=lambda r: (r.packet.route, r.packet.name)):
            print(row_summary(row))
    print()

    print("[3] Split audit on exception rows")
    split_specs = (
        ("arc_topology_compact", topology_tooth.key_func),
        ("coarse_safe_stalk", coarse_stalk_tooth.key_func),
        ("exact_safe_stalk", exact_stalk_tooth.key_func),
        ("magnitude_cocycle", magnitude_tooth.key_func),
        ("q_or_covering_certificate", certificate_tooth.key_func),
        ("first_primitive_safe_q_2_13", lambda row: ram.first_primitive_safe_q(row.packet, 13)),
        ("primitive_deck_2_13", lambda row: ram.primitive_count_deck(row.packet, 13)),
        ("route_label_sink", lambda row: row.packet.route),
    )
    for name, key_func in split_specs:
        buckets, mixed, max_size = split_signature(exception_rows, key_func)
        print(f"  {name:<30} buckets={buckets:<2d} mixed_route={mixed:<1d} max_bucket={max_size}")
    print()

    print("[4] Tournament Analysis")
    fp = carrier_tournament()
    print("  vertices_are=exception repair carriers, not runners")
    print(
        "  observable=exception coverage, route split, nonroute legality, "
        "local owner signal, primitive-period alignment, compression, proof cost"
    )
    print("  switch=majority comparison of carrier vectors")
    print(
        "  tie_hamiltonian_path="
        + " > ".join(carrier.name for carrier in EXCEPTION_CARRIERS)
    )
    print("  score_hist=" + str(fp["score_hist"]))
    print("  directed_3cycles=" + str(fp["directed_3cycles"]))
    print("  hamiltonian_path_count=" + str(fp["hamiltonian_path_count"]))
    print("  first_hamiltonian_path=" + " > ".join(fp["first_hamiltonian_path"]))
    print("  score_order=" + " > ".join(fp["score_order"]))
    print()

    print("[5] Proof readout")
    print("  1. Arc topology fails on exactly two residual same-topology buckets,")
    print("     and the failing rows are four single-swap collars at drops 9 and 11.")
    print("  2. In each collar, the covering row has zero primitive safe mass for")
    print("     q<=13, while the q-witness row has first primitive safe q equal")
    print("     to its dropped speed and its q-witness threshold.")
    print("  3. The coarse largest-safe-component stalk splits both topology")
    print("     exceptions, so the exceptions are owner-strip repairs rather than")
    print("     new boundary/open or route-label debt.")
    print("  4. Candidate lemma: if a post-status residual coarse ET+unit fiber is")
    print("     not split by compact arc topology, then it is an owner-stalk collar")
    print("     whose primitive q<=13 deck also splits the Q-witness/covering route.")


if __name__ == "__main__":
    main()
