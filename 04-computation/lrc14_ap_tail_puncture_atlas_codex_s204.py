#!/usr/bin/env python3
"""AP-tail puncture atlas for LRC14.

HYP-3029 left two coarse largest-stalk mixed-route teeth in the hard
automatic fiber `MFCMMCCFFFCCC`:

    single swap 13->159 / single swap 13->117
    single swap 13->118 / single swap 13->104

This script asks whether those teeth are a new residual or an already visible
AP-tail clock.  For the one-tail family

    S_m = {1, 2, ..., 12, m},

there is a very small certificate split:

  * if 13 does not divide m, then t = 1/13 is a strict witness;
  * if m = 13s with s >= 2, then t = s/(13s+1) is a strict witness;
  * m = 13 is the AP boundary atom.

The residual mixed pairs are exactly pairs where a coarse mod-14/stalk quotient
kept the same owner strip but forgot the `13 | m` puncture.  Tournament
Analysis is over repair clocks/proof carriers, not runners.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
from pathlib import Path
import argparse
import os
import sys


REPO = Path(__file__).resolve().parents[1]
N = 14
AP_TAIL_BASE = tuple(range(1, 13))
AP = tuple(range(1, 14))
THRESHOLD = Fraction(1, N)
TARGET_WORD = "MFCMMCCFFFCCC"


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s193 = load_module(
    "lrc14_safe_component_stalk_s193",
    REPO / "04-computation" / "lrc14_safe_component_stalk_codex_s193.py",
)


def fmt(fr: Fraction) -> str:
    return str(fr.numerator) if fr.denominator == 1 else f"{fr.numerator}/{fr.denominator}"


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def circle_distance(x: Fraction) -> Fraction:
    y = frac_part(x)
    return min(y, Fraction(1) - y)


def row_height(row: tuple[int, ...], t: Fraction) -> Fraction:
    return min(circle_distance(Fraction(v) * t) for v in row)


def ap_tail_m(speeds: tuple[int, ...]) -> int | None:
    if tuple(speeds[:12]) != AP_TAIL_BASE:
        return None
    if len(speeds) != 13:
        return None
    m = speeds[-1]
    if m in AP_TAIL_BASE:
        return None
    return m


@dataclass(frozen=True)
class TailCertificate:
    m: int
    kind: str
    witness: Fraction
    height: Fraction
    margin: Fraction
    proof_line: str

    @property
    def strict(self) -> bool:
        return self.height > THRESHOLD


def ap_tail_certificate(m: int) -> TailCertificate:
    row = AP_TAIL_BASE + (m,)
    if m == 13:
        witness = THRESHOLD
        height = row_height(row, witness)
        return TailCertificate(
            m=m,
            kind="ap_boundary",
            witness=witness,
            height=height,
            margin=height - THRESHOLD,
            proof_line="m=13 is the AP equality atom at t=1/14.",
        )

    if m % 13 != 0:
        witness = Fraction(1, 13)
        height = row_height(row, witness)
        return TailCertificate(
            m=m,
            kind="q13_puncture",
            witness=witness,
            height=height,
            margin=height - THRESHOLD,
            proof_line=(
                "13 does not divide m, so every residue v/13 is nonzero; "
                "the AP core and the tail all sit at distance at least 1/13."
            ),
        )

    s = m // 13
    witness = Fraction(s, 13 * s + 1)
    height = row_height(row, witness)
    return TailCertificate(
        m=m,
        kind="fixed_point_tail",
        witness=witness,
        height=height,
        margin=height - THRESHOLD,
        proof_line=(
            "m=13s: at t=s/(13s+1), runners 1 and m bind at height t; "
            "the AP core has distances j*t or 1-j*t, all at least t."
        ),
    )


def tail_proof_check(max_m: int) -> dict[str, object]:
    certs = [ap_tail_certificate(m) for m in range(13, max_m + 1)]
    non_strict = [cert for cert in certs if cert.m != 13 and not cert.strict]
    strict = [cert for cert in certs if cert.strict]
    return {
        "certs": certs,
        "counts": Counter(cert.kind for cert in certs),
        "non_strict": non_strict,
        "min_margin": min((cert.margin for cert in strict), default=Fraction(0)),
        "min_margin_rows": [
            cert for cert in strict
            if cert.margin == min((c.margin for c in strict), default=Fraction(0))
        ],
    }


def residue_terminal(row) -> tuple[object, ...]:
    return s193.residue_terminal(row.packet)


def magnitude_key(row) -> tuple[object, ...]:
    return s193.magnitude_key(row.packet)


def target_rows(single_limit: int, two_swap_limit: int, alias_depth: int, lcm_tail_max: int, workers: int):
    bank = s193.lp.build_bank(single_limit, two_swap_limit, alias_depth, lcm_tail_max)
    target = [row for row in bank if s193.lp.automatic_word(row[2]) == TARGET_WORD]
    packets = s193.lp.compute_packets(target, workers)
    rows = [s193.Row(packet=packet, stalk=s193.largest_component_stalk(packet)) for packet in packets]
    return target, packets, rows


@dataclass(frozen=True)
class Carrier:
    name: str
    description: str
    key_func: object
    exact_period: int
    owner_strip: int
    local_geometry: int
    family_proof: int
    proof_cost: int


def tail_clock(row) -> tuple[object, ...]:
    m = ap_tail_m(tuple(row.packet.speeds))
    if m is None:
        return ("not_ap_tail", row.packet.name)
    cert = ap_tail_certificate(m)
    return (cert.kind, m % 13 == 0)


def tail_certificate_key(row) -> tuple[object, ...]:
    m = ap_tail_m(tuple(row.packet.speeds))
    if m is None:
        return ("not_ap_tail", row.packet.name)
    cert = ap_tail_certificate(m)
    return (cert.kind, cert.witness, cert.height)


def peak_height_key(row) -> tuple[object, ...]:
    return (row.stalk.peak_time, row.stalk.peak_height)


def exact_stalk_key(row) -> tuple[object, ...]:
    return row.stalk.exact_key


def coarse_stalk_key(row) -> tuple[object, ...]:
    return next(carrier for carrier in s193.carriers() if carrier.name == "coarse_component_stalk").key_func(row)


def carriers() -> tuple[Carrier, ...]:
    return (
        Carrier(
            "residue_terminal",
            "automatic word plus mod-14 terminal-state fiber",
            residue_terminal,
            0,
            1,
            0,
            0,
            1,
        ),
        Carrier(
            "coarse_stalk",
            "largest-component coarse owner strip and height bucket",
            coarse_stalk_key,
            0,
            4,
            2,
            0,
            2,
        ),
        Carrier(
            "coarse_plus_q13_bit",
            "coarse stalk plus the hidden 13-divisibility puncture",
            lambda row: (coarse_stalk_key(row), tail_clock(row)),
            5,
            4,
            2,
            4,
            2,
        ),
        Carrier(
            "coarse_plus_tail_certificate",
            "coarse stalk plus explicit AP-tail witness type and height",
            lambda row: (coarse_stalk_key(row), tail_certificate_key(row)),
            5,
            5,
            4,
            5,
            3,
        ),
        Carrier(
            "coarse_plus_peak_height",
            "coarse stalk plus exact largest-stalk peak time and height",
            lambda row: (coarse_stalk_key(row), peak_height_key(row)),
            2,
            5,
            5,
            3,
            3,
        ),
        Carrier(
            "exact_stalk",
            "exact largest safe-component stalk",
            lambda row: (residue_terminal(row), exact_stalk_key(row)),
            2,
            5,
            5,
            3,
            4,
        ),
        Carrier(
            "magnitude_cocycle",
            "exact M/q/Farey/tail magnitude cocycle",
            lambda row: (row.packet.automatic_word, magnitude_key(row)),
            4,
            1,
            1,
            2,
            4,
        ),
    )


@dataclass(frozen=True)
class CarrierReport:
    carrier: Carrier
    fibers: int
    mixed_route: int
    max_mixed: int
    mixed_rows: int

    @property
    def route_purity(self) -> int:
        if self.fibers == 0:
            return 0
        return (1000 * (self.fibers - self.mixed_route)) // self.fibers

    @property
    def vector(self) -> tuple[int, ...]:
        return (
            self.route_purity,
            -self.max_mixed,
            self.carrier.exact_period,
            self.carrier.owner_strip,
            self.carrier.local_geometry,
            self.carrier.family_proof,
            6 - self.carrier.proof_cost,
        )


def carrier_report(rows: list, carrier: Carrier) -> CarrierReport:
    groups: dict[object, list] = defaultdict(list)
    for row in rows:
        groups[carrier.key_func(row)].append(row)
    mixed = [group for group in groups.values() if len({row.packet.route for row in group}) > 1]
    return CarrierReport(
        carrier=carrier,
        fibers=len(groups),
        mixed_route=len(mixed),
        max_mixed=max((len(group) for group in mixed), default=0),
        mixed_rows=sum(len(group) for group in mixed),
    )


def orient(a: CarrierReport, b: CarrierReport) -> str:
    av = a.vector
    bv = b.vector
    awins = sum(x > y for x, y in zip(av, bv))
    bwins = sum(x < y for x, y in zip(av, bv))
    if awins != bwins:
        return a.carrier.name if awins > bwins else b.carrier.name
    names = [carrier.name for carrier in carriers()]
    return a.carrier.name if names.index(a.carrier.name) < names.index(b.carrier.name) else b.carrier.name


def tournament_fingerprint(reports: list[CarrierReport]) -> tuple[dict[int, int], int, list[int], int, list[str]]:
    names = [report.carrier.name for report in reports]
    edges: dict[tuple[str, str], bool] = {}
    score = Counter()
    for a, b in combinations(reports, 2):
        winner = orient(a, b)
        loser = b.carrier.name if winner == a.carrier.name else a.carrier.name
        edges[(winner, loser)] = True
        score[winner] += 1
        score.setdefault(loser, 0)

    def beats(x: str, y: str) -> bool:
        return (x, y) in edges

    c3 = 0
    for a, b, c in combinations(names, 3):
        if (
            beats(a, b)
            and beats(b, c)
            and beats(c, a)
            or beats(a, c)
            and beats(c, b)
            and beats(b, a)
        ):
            c3 += 1

    graph: dict[str, list[str]] = defaultdict(list)
    reverse: dict[str, list[str]] = defaultdict(list)
    for winner, loser in edges:
        graph[winner].append(loser)
        reverse[loser].append(winner)

    def reach(seed: str, graph_map: dict[str, list[str]]) -> set[str]:
        seen = {seed}
        stack = [seed]
        while stack:
            v = stack.pop()
            for w in graph_map[v]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        return seen

    remaining = set(names)
    scc_sizes = []
    while remaining:
        seed = min(remaining)
        comp = reach(seed, graph) & reach(seed, reverse)
        scc_sizes.append(len(comp))
        remaining -= comp

    hp = 0
    first_path: tuple[str, ...] | None = None
    for order in permutations(names):
        if all(beats(order[i], order[i + 1]) for i in range(len(order) - 1)):
            hp += 1
            if first_path is None:
                first_path = order

    ranking = sorted(names, key=lambda name: (-score[name], name))
    if first_path is not None:
        ranking = list(first_path)
    return dict(sorted(Counter(score.values()).items())), c3, sorted(scc_sizes, reverse=True), hp, ranking


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge")
    print("  considered vertices:")
    print("    runners, the AP core, tail parameter m, mod-14 owner strips,")
    print("    q=13 punctures, exact-period witnesses, fixed-point equations,")
    print("    safe-component stalks, route labels, and proof obligations.")
    print("  chosen vertices:")
    print("    AP-tail repair clocks and stalk quotients, not runners.")
    print("  preserved LRC predicate:")
    print("    existence of a strict witness for S_m={1,...,12,m}, and route")
    print("    purity in the target AP/GW automatic fiber.")
    print("  destroyed information:")
    print("    raw runner identity beyond the AP-tail parameter, non-largest")
    print("    components, and non-tail packet labels unless reattached.")
    print("  challenged assumption:")
    print("    a mod-14 owner strip is not a complete clock; it must be paired")
    print("    with the q=13 puncture/fixed-point address.")
    print()


def print_tail_theorem(scan: dict[str, object], max_m: int) -> None:
    print("[1] AP-tail two-clock theorem")
    print(f"  scanned_m=13..{max_m}")
    print(f"  counts={dict(scan['counts'])}")
    print(f"  non_strict_except_AP={len(scan['non_strict'])}")
    print(f"  min_strict_margin={fmt(scan['min_margin'])}")
    for cert in scan["min_margin_rows"][:5]:
        print(
            f"    min_margin_row m={cert.m} kind={cert.kind} "
            f"t={fmt(cert.witness)} height={fmt(cert.height)}"
        )
    print("  proof:")
    print("    If 13 does not divide m, t=1/13 gives every speed a nonzero")
    print("    residue mod 13, hence distance at least 1/13 > 1/14.")
    print("    If m=13s and s>=2, t=s/(13s+1) makes runners 1 and m")
    print("    bind at distance t.  For 2<=j<=12 in the AP core, the")
    print("    distance is j*t or 1-j*t, and both are at least t; moreover")
    print("    t>1/14 exactly when s>1.  The single exception s=1 is AP.")
    print()


def print_residual_atlas(rows: list) -> None:
    coarse = next(carrier for carrier in s193.carriers() if carrier.name == "coarse_component_stalk")
    groups: dict[object, list] = defaultdict(list)
    for row in rows:
        groups[coarse.key_func(row)].append(row)
    residuals = [
        group for group in groups.values()
        if len({row.packet.route for row in group}) > 1
    ]
    residuals.sort(key=lambda group: tuple(row.packet.name for row in group))

    print("[2] HYP-3029 residual tooth atlas")
    print(f"  coarse_stalk_mixed_route_fibers={len(residuals)}")
    for idx, group in enumerate(residuals, 1):
        print(f"  residual[{idx}] size={len(group)} routes={dict(Counter(row.packet.route for row in group))}")
        for row in sorted(group, key=lambda item: item.packet.name):
            packet = row.packet
            m = ap_tail_m(tuple(packet.speeds))
            assert m is not None
            cert = ap_tail_certificate(m)
            s = m // 13 if m % 13 == 0 else None
            print(
                "    {name:24s} route={route:16s} m={m:3d} mod14={r14:2d} "
                "mod13={r13:2d} s={s:>4s} cert={kind:16s} "
                "t={t:>7s} h={h:>7s} M={M:>7s} stalk_len={length:>7s} "
                "peak={peak:>7s}".format(
                    name=packet.name.replace("single swap ", ""),
                    route=packet.route,
                    m=m,
                    r14=m % 14,
                    r13=m % 13,
                    s=str(s) if s is not None else "-",
                    kind=cert.kind,
                    t=fmt(cert.witness),
                    h=fmt(cert.height),
                    M=fmt(packet.M),
                    length=fmt(row.stalk.length),
                    peak=fmt(row.stalk.peak_height),
                )
            )
    print("  readout:")
    print("    The coarse stalk kept the same mod-14 owner strip but erased")
    print("    whether the tail punctures q=13.  Reattaching the q=13 bit")
    print("    or the explicit AP-tail certificate makes these fibers route-pure.")
    print()


def print_repair_ladder(rows: list) -> list[CarrierReport]:
    reports = [carrier_report(rows, carrier) for carrier in carriers()]
    print("[3] Repair-carrier ladder on the target automatic fiber")
    print("  carrier                         fibers mixR maxMix mixRows purity")
    for report in reports:
        print(
            "  {name:31s} {fibers:6d} {mix:4d} {maxm:6d} {rows:7d} {purity:5.1f}%".format(
                name=report.carrier.name,
                fibers=report.fibers,
                mix=report.mixed_route,
                maxm=report.max_mixed,
                rows=report.mixed_rows,
                purity=report.route_purity / 10,
            )
        )
    print()
    return reports


def print_tournament(reports: list[CarrierReport]) -> None:
    hist, c3, scc, hp, ranking = tournament_fingerprint(reports)
    print("[4] Tournament Analysis")
    print("  vertices_are=repair clocks and proof carriers, not runners")
    print("  observable=route_purity,max_mixed,exact_period_retention,owner_strip,")
    print("             local_geometry,family_proof_scope,proof_cost")
    print("  switch=majority comparison of observable vectors")
    print("  tie_hamiltonian_path=" + " > ".join(carrier.name for carrier in carriers()))
    print(f"  score_hist={hist}")
    print(f"  directed_3cycles={c3}")
    print(f"  scc_sizes={scc}")
    print(f"  hamiltonian_path_count={hp}")
    print("  first_hamiltonian_path=" + " > ".join(ranking))
    print()


def print_proof_readout() -> None:
    print("[5] Proof readout")
    print("  1. The AP-tail family S_m={1,...,12,m} is completely discharged")
    print("     by the q=13 puncture clock plus the fixed-point tail clock.")
    print("  2. The two S193 coarse-stalk residuals are not new F7 debt.  They")
    print("     are controlled-kernel examples: mod-14 owner geometry collided,")
    print("     and the missing coordinate was exactly m mod 13.")
    print("  3. HYP-3031's zeta repair class for these teeth is nested")
    print("     refinement with an owner-strip witness: the quotient may forget")
    print("     the mixed coordinate only after proving the AP-tail formula.")
    print("  4. New LRC14 target: identify other large automatic-fiber residuals")
    print("     that are just AP-core plus one or two tail punctures, then prove")
    print("     the same reciprocal fixed-point certificate before invoking")
    print("     heavier Fejer/THM-572 machinery.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-tail", type=int, default=400)
    parser.add_argument("--single-limit", type=int, default=180)
    parser.add_argument("--two-swap-limit", type=int, default=36)
    parser.add_argument("--alias-depth", type=int, default=4)
    parser.add_argument("--lcm-tail-max", type=int, default=5)
    parser.add_argument("--workers", type=int, default=max(1, min(os.cpu_count() or 1, 8)))
    args = parser.parse_args()

    scan = tail_proof_check(args.max_tail)
    target, packets, rows = target_rows(
        args.single_limit,
        args.two_swap_limit,
        args.alias_depth,
        args.lcm_tail_max,
        args.workers,
    )

    print("=== LRC14 AP-tail puncture atlas S204 ===")
    print(
        f"tail_family=S_m={{1,...,12,m}} max_tail={args.max_tail} "
        f"target_word={TARGET_WORD}"
    )
    print(
        f"bank=HYP-2963 target_rows={len(target)} packets={len(packets)} "
        f"single_limit={args.single_limit} two_swap_limit={args.two_swap_limit}"
    )
    print()
    print_assumption_challenge()
    print_tail_theorem(scan, args.max_tail)
    print_residual_atlas(rows)
    reports = print_repair_ladder(rows)
    print_tournament(reports)
    print_proof_readout()


if __name__ == "__main__":
    main()
