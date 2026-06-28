#!/usr/bin/env python3
"""HYP-3403 scout: shadow-charge packet gluing for LRC14.

This is a proof-route stress test, not a proof.  It takes the actual-packet
bank from HYP-3311 and asks which of the current creative reframes are
observable on theorem-facing rows:

* ambient index / C3 unit-slot data;
* Q(sqrt(-7)) quadratic sign data;
* nonunit covering residue data;
* 2-adic magnitude and exact height data;
* transfer/state-lift exits.

The controlled-forgetting test is deliberately conservative.  A sidecar is
useful only if it separates theorem exits without becoming an exact row hash.
Tournament Analysis uses proof sidecars as vertices, not runners or arcs.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]


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
    "h3401_actual_packet_sheaf",
    "04-computation/lrc14_actual_packet_sheaf_instantiation_codex_20260628.py",
)
h3310 = h3311.h3310


SLOTS = tuple(h3310.c3_pair_orbit())


@dataclass(frozen=True)
class ChargeRow:
    name: str
    route: str
    kernel_flag: str
    family: str
    q_threshold: int
    transfer: str
    state_lift: bool
    strict_safe_positive: bool
    coarse_sheaf_base: tuple[object, ...]
    unit_profile: tuple[int, int, int]
    chosen_denominator: int
    speeds: tuple[int, ...]
    unit_c3_slot_word: tuple[int, ...]
    unit_qsqrt7_word: tuple[tuple[int, int], ...]
    qsqrt7_balance: tuple[int, int, int]
    nonunit_residue_word: tuple[int, ...]
    nonunit_v2_word: tuple[int, ...]
    nonunit_cover_signature: tuple[tuple[int, int], ...]
    nonunit_exact_height_word: tuple[tuple[int, int, int], ...]

    @property
    def q_bucket(self) -> str:
        if self.q_threshold < 14:
            return "lt14"
        if self.q_threshold == 14:
            return "eq14"
        return "gt14"


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


def qsqrt_balance(speeds: tuple[int, ...]) -> tuple[int, int, int]:
    counts = Counter(h3310.chi7(speed) for speed in speeds)
    return (counts[1], counts[-1], counts[0])


def unit_slot_word(speeds: tuple[int, ...]) -> tuple[int, ...]:
    counts = [0 for _ in SLOTS]
    for speed in speeds:
        if gcd(speed, 14) == 1:
            idx = slot_index(speed)
            if idx is not None:
                counts[idx] += 1
    return tuple(counts)


def unit_qsqrt7_word(speeds: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    rows = []
    for speed in speeds:
        if gcd(speed, 14) == 1:
            idx = slot_index(speed)
            if idx is not None:
                rows.append((idx, h3310.chi7(speed)))
    return tuple(sorted(rows))


def nonunit_height_word(speeds: tuple[int, ...]) -> tuple[tuple[int, int, int], ...]:
    return tuple(
        sorted(
            (speed % 14, h3310.v2(speed), odd_part(speed))
            for speed in speeds
            if gcd(speed, 14) != 1
        )
    )


def build_charge_rows() -> list[ChargeRow]:
    actual_rows = {row.name: row for row in h3311.build_rows()}
    bank = h3311.s154.curated_bank(alias_depth=2, lcm_tail_max=5)
    packets = h3311.s154.lp.compute_packets(bank, workers=1)
    out: list[ChargeRow] = []
    for packet in packets:
        row = actual_rows[packet.name]
        speeds = tuple(packet.speeds)
        out.append(
            ChargeRow(
                name=row.name,
                route=row.route,
                kernel_flag=row.kernel_flag,
                family=row.family,
                q_threshold=row.q_threshold,
                transfer=row.transfer,
                state_lift=row.state_lift,
                strict_safe_positive=row.strict_safe_mu > 0,
                coarse_sheaf_base=row.coarse_sheaf_base,
                unit_profile=row.unit_profile,
                chosen_denominator=row.chosen_denominator,
                speeds=speeds,
                unit_c3_slot_word=unit_slot_word(speeds),
                unit_qsqrt7_word=unit_qsqrt7_word(speeds),
                qsqrt7_balance=qsqrt_balance(speeds),
                nonunit_residue_word=row.nonunit_residue_word,
                nonunit_v2_word=row.nonunit_v2_word,
                nonunit_cover_signature=row.nonunit_cover_signature,
                nonunit_exact_height_word=nonunit_height_word(speeds),
            )
        )
    return out


def fiber_table(rows: list[ChargeRow], key_fn) -> dict[tuple[object, ...], list[ChargeRow]]:
    fibers: dict[tuple[object, ...], list[ChargeRow]] = defaultdict(list)
    for row in rows:
        fibers[key_fn(row)].append(row)
    return fibers


@dataclass(frozen=True)
class FiberStats:
    fibers: int
    mixed_kernel: int
    mixed_route: int
    mixed_transfer: int
    max_kernel_width: int
    max_route_width: int
    max_fiber_size: int


def fiber_stats(rows: list[ChargeRow], key_fn) -> FiberStats:
    fibers = fiber_table(rows, key_fn)
    mixed_kernel = 0
    mixed_route = 0
    mixed_transfer = 0
    max_kernel_width = 0
    max_route_width = 0
    max_fiber_size = 0
    for fiber in fibers.values():
        kernels = {row.kernel_flag for row in fiber}
        routes = {row.route for row in fiber}
        transfers = {row.transfer for row in fiber}
        mixed_kernel += int(len(kernels) > 1)
        mixed_route += int(len(routes) > 1)
        mixed_transfer += int(len(transfers) > 1)
        max_kernel_width = max(max_kernel_width, len(kernels))
        max_route_width = max(max_route_width, len(routes))
        max_fiber_size = max(max_fiber_size, len(fiber))
    return FiberStats(
        fibers=len(fibers),
        mixed_kernel=mixed_kernel,
        mixed_route=mixed_route,
        mixed_transfer=mixed_transfer,
        max_kernel_width=max_kernel_width,
        max_route_width=max_route_width,
        max_fiber_size=max_fiber_size,
    )


@dataclass(frozen=True)
class Carrier:
    name: str
    key_label: str
    payload_cost: int
    preserves: frozenset[str]
    destroys: frozenset[str]
    priority: int


PAYLOAD_WEIGHTS = {
    "kernel_exit": 20,
    "route_exit": 12,
    "unit_c3_slot": 8,
    "quadratic_sign": 6,
    "nonunit_residue": 10,
    "nonunit_v2": 7,
    "exact_height": 4,
    "transfer_exit": 9,
    "qdiv_cusp": 5,
    "strict_open": 5,
    "state_lift": 5,
}


CARRIERS = [
    Carrier(
        "raw_coarse_sheaf_base",
        "coarse_base",
        1,
        frozenset({"qdiv_cusp", "strict_open", "state_lift"}),
        frozenset({"unit_c3_slot", "quadratic_sign", "nonunit_residue", "nonunit_v2", "exact_height"}),
        90,
    ),
    Carrier(
        "ambient_index_unit_c3",
        "coarse_plus_unit_c3",
        2,
        frozenset({"unit_c3_slot", "qdiv_cusp", "strict_open", "state_lift"}),
        frozenset({"quadratic_sign", "nonunit_residue", "nonunit_v2", "exact_height"}),
        80,
    ),
    Carrier(
        "imaginary_quadratic_balance",
        "coarse_plus_qsqrt_balance",
        2,
        frozenset({"quadratic_sign", "qdiv_cusp", "strict_open", "state_lift"}),
        frozenset({"unit_c3_slot", "nonunit_residue", "nonunit_v2", "exact_height"}),
        70,
    ),
    Carrier(
        "c3_plus_qsqrt_binding_packet",
        "coarse_plus_c3_qsqrt",
        3,
        frozenset({"unit_c3_slot", "quadratic_sign", "qdiv_cusp", "strict_open", "state_lift"}),
        frozenset({"nonunit_residue", "nonunit_v2", "exact_height"}),
        60,
    ),
    Carrier(
        "covering_residue_sheaf",
        "coarse_plus_nonunit_residue",
        3,
        frozenset({"kernel_exit", "route_exit", "nonunit_residue", "qdiv_cusp", "strict_open", "state_lift"}),
        frozenset({"unit_c3_slot", "quadratic_sign", "nonunit_v2", "exact_height"}),
        10,
    ),
    Carrier(
        "two_adic_magnitude_gate",
        "coarse_plus_nonunit_v2",
        3,
        frozenset({"nonunit_v2", "qdiv_cusp", "strict_open", "state_lift"}),
        frozenset({"unit_c3_slot", "quadratic_sign", "nonunit_residue", "exact_height"}),
        50,
    ),
    Carrier(
        "residue_magnitude_cover_signature",
        "coarse_plus_cover_signature",
        4,
        frozenset(
            {
                "kernel_exit",
                "route_exit",
                "nonunit_residue",
                "nonunit_v2",
                "qdiv_cusp",
                "strict_open",
                "state_lift",
            }
        ),
        frozenset({"unit_c3_slot", "quadratic_sign", "exact_height"}),
        20,
    ),
    Carrier(
        "three_shadow_packet",
        "coarse_plus_c3_qsqrt_cover",
        6,
        frozenset(
            {
                "kernel_exit",
                "route_exit",
                "unit_c3_slot",
                "quadratic_sign",
                "nonunit_residue",
                "nonunit_v2",
                "qdiv_cusp",
                "strict_open",
                "state_lift",
            }
        ),
        frozenset({"exact_height", "transfer_exit"}),
        30,
    ),
    Carrier(
        "exact_height_hinge_oracle",
        "coarse_plus_exact_height",
        9,
        frozenset({"kernel_exit", "route_exit", "exact_height", "qdiv_cusp", "strict_open", "state_lift"}),
        frozenset({"unit_c3_slot", "quadratic_sign", "nonunit_residue", "nonunit_v2"}),
        100,
    ),
    Carrier(
        "analytic_lifting_ledger",
        "coarse_plus_c3_qsqrt_cover_transfer",
        7,
        frozenset(
            {
                "kernel_exit",
                "route_exit",
                "transfer_exit",
                "unit_c3_slot",
                "quadratic_sign",
                "nonunit_residue",
                "nonunit_v2",
                "qdiv_cusp",
                "strict_open",
                "state_lift",
            }
        ),
        frozenset({"exact_height"}),
        40,
    ),
]


def key_fn(label: str):
    if label == "coarse_base":
        return lambda row: row.coarse_sheaf_base
    if label == "coarse_plus_unit_c3":
        return lambda row: row.coarse_sheaf_base + (row.unit_c3_slot_word,)
    if label == "coarse_plus_qsqrt_balance":
        return lambda row: row.coarse_sheaf_base + (row.qsqrt7_balance,)
    if label == "coarse_plus_c3_qsqrt":
        return lambda row: row.coarse_sheaf_base + (row.unit_c3_slot_word, row.qsqrt7_balance)
    if label == "coarse_plus_nonunit_residue":
        return lambda row: row.coarse_sheaf_base + (row.nonunit_residue_word,)
    if label == "coarse_plus_nonunit_v2":
        return lambda row: row.coarse_sheaf_base + (row.nonunit_v2_word,)
    if label == "coarse_plus_cover_signature":
        return lambda row: row.coarse_sheaf_base + (row.nonunit_cover_signature,)
    if label == "coarse_plus_c3_qsqrt_cover":
        return lambda row: row.coarse_sheaf_base + (
            row.unit_c3_slot_word,
            row.qsqrt7_balance,
            row.nonunit_cover_signature,
        )
    if label == "coarse_plus_exact_height":
        return lambda row: row.coarse_sheaf_base + (row.nonunit_exact_height_word,)
    if label == "coarse_plus_c3_qsqrt_cover_transfer":
        return lambda row: row.coarse_sheaf_base + (
            row.unit_c3_slot_word,
            row.qsqrt7_balance,
            row.nonunit_cover_signature,
            row.transfer,
        )
    raise ValueError(label)


def carrier_stats(rows: list[ChargeRow]) -> dict[str, FiberStats]:
    return {carrier.name: fiber_stats(rows, key_fn(carrier.key_label)) for carrier in CARRIERS}


def preserve_score(carrier: Carrier) -> int:
    kept = sum(PAYLOAD_WEIGHTS[item] for item in carrier.preserves)
    lost = sum(max(1, PAYLOAD_WEIGHTS[item] // 3) for item in carrier.destroys)
    return kept - lost


def carrier_rank_key(carrier: Carrier, stats: FiberStats) -> tuple[int, int, int, int, int, int]:
    return (
        stats.mixed_kernel,
        stats.mixed_route,
        carrier.payload_cost,
        stats.mixed_transfer,
        -preserve_score(carrier),
        carrier.priority,
    )


def adjacency(rows: list[ChargeRow]) -> dict[tuple[str, str], bool]:
    stats = carrier_stats(rows)
    out: dict[tuple[str, str], bool] = {}
    for a, b in combinations(CARRIERS, 2):
        ak = carrier_rank_key(a, stats[a.name])
        bk = carrier_rank_key(b, stats[b.name])
        if ak < bk:
            out[(a.name, b.name)] = True
            out[(b.name, a.name)] = False
        else:
            out[(b.name, a.name)] = True
            out[(a.name, b.name)] = False
    return out


def directed_3cycles(adj: dict[tuple[str, str], bool]) -> list[tuple[str, str, str]]:
    names = [carrier.name for carrier in CARRIERS]
    cycles = []
    for a, b, c in combinations(names, 3):
        if adj[(a, b)] and adj[(b, c)] and adj[(c, a)]:
            cycles.append((a, b, c))
        elif adj[(a, c)] and adj[(c, b)] and adj[(b, a)]:
            cycles.append((a, c, b))
    return cycles


def hamiltonian_path_count(adj: dict[tuple[str, str], bool]) -> int:
    names = [carrier.name for carrier in CARRIERS]
    count = 0

    def extend(path: list[str], remaining: set[str]) -> None:
        nonlocal count
        if not remaining:
            count += 1
            return
        for nxt in list(remaining):
            if adj[(path[-1], nxt)]:
                remaining.remove(nxt)
                path.append(nxt)
                extend(path, remaining)
                path.pop()
                remaining.add(nxt)

    for start in names:
        rest = set(names)
        rest.remove(start)
        extend([start], rest)
    return count


def scc_sizes(adj: dict[tuple[str, str], bool]) -> list[int]:
    names = [carrier.name for carrier in CARRIERS]
    unseen = set(names)
    sizes = []
    while unseen:
        start = next(iter(unseen))
        forward = reachable(start, names, adj)
        backward = {name for name in names if start in reachable(name, names, adj)}
        comp = forward & backward
        sizes.append(len(comp))
        unseen -= comp
    return sorted(sizes, reverse=True)


def reachable(start: str, names: list[str], adj: dict[tuple[str, str], bool]) -> set[str]:
    seen = {start}
    stack = [start]
    while stack:
        cur = stack.pop()
        for nxt in names:
            if nxt != cur and nxt not in seen and adj[(cur, nxt)]:
                seen.add(nxt)
                stack.append(nxt)
    return seen


def mixed_fibers(rows: list[ChargeRow], label: str, target: str) -> list[list[ChargeRow]]:
    get_key = key_fn(label)
    out = []
    for fiber in fiber_table(rows, get_key).values():
        values = {getattr(row, target) for row in fiber}
        if len(values) > 1:
            out.append(sorted(fiber, key=lambda row: (getattr(row, target), row.name)))
    out.sort(key=lambda fiber: (-len(fiber), tuple(row.name for row in fiber)))
    return out


def print_summary(rows: list[ChargeRow]) -> None:
    print("HYP-3403 shadow-charge packet gluing")
    print("=" * 78)
    print("bank=HYP-3311 actual-packet rows over the curated HYP-2969 ledger")
    print(f"rows={len(rows)}")
    print(f"kernel_flag_hist={dict(Counter(row.kernel_flag for row in rows))}")
    print(f"route_hist={dict(Counter(row.route for row in rows))}")
    print(f"q_bucket_hist={dict(Counter(row.q_bucket for row in rows))}")
    print()


def print_sidecar_table(rows: list[ChargeRow]) -> None:
    print("CONTROLLED-FORGETTING SIDECARE TABLE")
    print(
        "  sidecar                              fibers kernel route transfer "
        "maxK maxR maxFiber"
    )
    stats = carrier_stats(rows)
    for carrier in CARRIERS:
        s = stats[carrier.name]
        print(
            f"  {carrier.name:36s} "
            f"{s.fibers:6d} {s.mixed_kernel:6d} {s.mixed_route:5d} "
            f"{s.mixed_transfer:8d} {s.max_kernel_width:4d} {s.max_route_width:4d} "
            f"{s.max_fiber_size:8d}"
        )
    print()


def print_failure_modes(rows: list[ChargeRow]) -> None:
    print("REFRAME FAILURE MODES")
    for label, title in (
        ("coarse_plus_c3_qsqrt", "ambient C3 + Qsqrt(-7) binding packet"),
        ("coarse_plus_nonunit_v2", "2-adic v2 word alone"),
        ("coarse_plus_nonunit_residue", "nonunit residue word"),
    ):
        mixed = mixed_fibers(rows, label, "kernel_flag")
        print(f"  {title}: mixed_kernel_fibers={len(mixed)}")
        if mixed:
            fiber = mixed[0]
            print(f"    first_mixed_names={[row.name for row in fiber]}")
            print(f"    first_mixed_kernels={sorted({row.kernel_flag for row in fiber})}")
    print()


def print_qdiv_and_height(rows: list[ChargeRow]) -> None:
    qgt = [row for row in rows if row.q_threshold > 14]
    same_residue = defaultdict(list)
    for row in rows:
        same_residue[(row.coarse_sheaf_base, row.nonunit_residue_word)].append(row)
    height_debts = [
        fiber
        for fiber in same_residue.values()
        if len({row.nonunit_exact_height_word for row in fiber}) > 1
    ]
    print("HEIGHT AND CUSP CHECKS")
    print(f"  qdiv_gt14_rows={len(qgt)}")
    print(f"  qdiv_gt14_kernel_hist={dict(Counter(row.kernel_flag for row in qgt))}")
    print(f"  same_residue_height_debt_fibers={len(height_debts)}")
    if height_debts:
        fiber = sorted(height_debts, key=lambda f: (-len(f), tuple(row.name for row in f)))[0]
        print(f"  largest_same_residue_height_debt={[row.name for row in fiber]}")
        for row in fiber:
            print(
                f"    {row.name}: kernel={row.kernel_flag}, "
                f"v2={row.nonunit_v2_word}, height={row.nonunit_exact_height_word}"
            )
    print()


def print_tournament(rows: list[ChargeRow]) -> None:
    stats = carrier_stats(rows)
    adj = adjacency(rows)
    cycles = directed_3cycles(adj)
    scc = scc_sizes(adj)
    path = sorted(CARRIERS, key=lambda carrier: carrier_rank_key(carrier, stats[carrier.name]))
    score_hist = Counter()
    for carrier in CARRIERS:
        key = carrier_rank_key(carrier, stats[carrier.name])
        score_hist[key[:4]] += 1
    print("TOURNAMENT ANALYSIS")
    print("  vertices=proof sidecars / shadow-charge carriers, not runners or arcs")
    print(
        "  pairwise_observable=(mixed kernel fibers, mixed route fibers, "
        "payload cost, mixed transfer fibers, preserved charge)"
    )
    print("  switch_gauge=A -> B iff A has fewer mixed theorem-exit fibers, then lower payload debt")
    print(f"  vertex_count={len(CARRIERS)}")
    print(f"  score_hist={dict(sorted(score_hist.items()))}")
    print(f"  directed_3cycles={len(cycles)}")
    print(f"  scc_sizes={scc}")
    print(f"  hamiltonian_path_count={hamiltonian_path_count(adj)}")
    print("  priority_path=")
    for carrier in path:
        s = stats[carrier.name]
        print(
            f"    {carrier.name}: key={carrier_rank_key(carrier, s)}, "
            f"fibers={s.fibers}, preserves={sorted(carrier.preserves)}, "
            f"destroys={sorted(carrier.destroys)}"
        )
    print()


def print_conclusion(rows: list[ChargeRow]) -> None:
    stats = carrier_stats(rows)
    c3q = stats["c3_plus_qsqrt_binding_packet"]
    residue = stats["covering_residue_sheaf"]
    v2 = stats["two_adic_magnitude_gate"]
    three = stats["three_shadow_packet"]
    print("CONCLUSION")
    print(
        "  The most creative-looking ambient reframes are descriptive on this "
        f"bank, not terminal: C3+Qsqrt(-7) still leaves {c3q.mixed_kernel} "
        "mixed kernel fiber."
    )
    print(
        "  The first theorem-facing repair is the covering residue sheaf: "
        f"mixed_kernel={residue.mixed_kernel}, mixed_route={residue.mixed_route}. "
        f"The v2 word alone leaves mixed_kernel={v2.mixed_kernel}."
    )
    print(
        "  The three-shadow packet (unit C3 slot, quadratic balance, and "
        "nonunit residue/v2 cover signature) also separates the bank, but it "
        "is intentionally more expensive.  Use it as the gluing interface; "
        "use the exact height word only as a debt detector for same-residue "
        "moves, not as a proof quotient."
    )
    print(
        "  Assumption challenged: the tournament vertices need not be runners, "
        "arcs, or residues.  Here the useful vertices are proof sidecars and "
        "shadow-charge carriers.  Preserved predicate: theorem-exit separability "
        "for LRC14 packet rows.  Destroyed coordinate: exact height/endpoint "
        "owner/off-grid floor unless retained or routed to named debt."
    )


def main() -> None:
    rows = build_charge_rows()
    print_summary(rows)
    print_sidecar_table(rows)
    print_failure_modes(rows)
    print_qdiv_and_height(rows)
    print_tournament(rows)
    print_conclusion(rows)


if __name__ == "__main__":
    main()
