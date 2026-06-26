#!/usr/bin/env python3
"""Ramanujan / exact-period route scheduler for LRC14 residual fibers.

HYP-3030 leaves 15 route-mixed coarse ET+Henselian-unit fibers on the full
HYP-2963 packet bank.  Those fibers are all strict-open, so they are harmless
for boundary/open status, but they still need a theorem-route scheduler.

This script tests a small primitive-period deck on exactly those S194 residual
rows.  The deck counts primitive phases a/q, 2 <= q <= 13, at which every
runner is at closed distance at least 1/14 from an integer.  It is a
Ramanujan-style exact-period carrier: primitive residues are the period-q
layers, and the optional trace deck records integer Ramanujan sums c_q(v).

Tournament Analysis declaration:
  vertices: proof carriers and scheduler decks, not runners;
  pairwise observable: status purity, route purity, primitive-period
    certifiability, harmonic trace visibility, compression, and proof cost;
  switch/gauge: majority comparison of observable vectors;
  tie Hamiltonian path: primitive_count_deck_2_13 > first_safe_q_2_13 >
    ramanujan_trace_deck_2_14 > coarse_et_unit_status_gate >
    exact_magnitude_cocycle > raw_residue_terminal_word.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
from math import gcd
from pathlib import Path
import argparse
import re
import sys


REPO = Path(__file__).resolve().parents[1]
S194_RESULT = REPO / "05-knowledge" / "results" / "lrc14_status_topology_gate_codex_s194.out"


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


fz = load_module(
    "lrc14_fiber_zipper_convergence_s200",
    REPO / "04-computation" / "lrc14_fiber_zipper_convergence_codex_s188.py",
)


def fmt(fr: Fraction) -> str:
    return str(fr.numerator) if fr.denominator == 1 else f"{fr.numerator}/{fr.denominator}"


def parse_s194_residual_names(path: Path = S194_RESULT) -> tuple[str, ...]:
    text = path.read_text()
    try:
        section = text.split("[2] Coarse ET+unit route-mixed fibers", 1)[1]
        section = section.split("[3] Tournament Analysis", 1)[0]
    except IndexError as exc:
        raise ValueError(f"could not find S194 residual section in {path}") from exc

    names: list[str] = []
    for line in section.splitlines():
        match = re.match(r"\s{4}(.+?)\s+status=", line)
        if match:
            names.append(match.group(1).strip())
    if not names:
        raise ValueError(f"no residual packet names parsed from {path}")
    return tuple(names)


def packets_for_names(
    names: tuple[str, ...],
    single_limit: int,
    two_swap_limit: int,
    alias_depth: int,
    lcm_tail_max: int,
    workers: int,
) -> list:
    bank = fz.lp.build_bank(single_limit, two_swap_limit, alias_depth, lcm_tail_max)
    rows_by_name = {name: (name, family, speeds) for name, family, speeds in bank}
    missing = [name for name in names if name not in rows_by_name]
    if missing:
        raise ValueError("S194 residual names not present in rebuilt bank: " + ", ".join(missing))
    rows = [rows_by_name[name] for name in names]
    return fz.lp.compute_packets(rows, workers)


def coarse_key(packet):
    split = next(s for s in fz.SPLITS if s.name == "coarse_et_unit_gate")
    return split.key_func(packet)


def group_by(packets: list, key_func) -> dict[object, list]:
    groups: dict[object, list] = defaultdict(list)
    for packet in packets:
        groups[key_func(packet)].append(packet)
    return groups


def mixed_route_groups(groups: dict[object, list]) -> list[list]:
    return [rows for rows in groups.values() if len({p.route for p in rows}) > 1]


def mixed_status_groups(groups: dict[object, list]) -> list[list]:
    return [rows for rows in groups.values() if len({fz.status(p) for p in rows}) > 1]


@lru_cache(maxsize=None)
def closed_safe_at_phase(speeds: tuple[int, ...], q: int, a: int) -> bool:
    return all(14 * min((v * a) % q, q - ((v * a) % q)) >= q for v in speeds)


@lru_cache(maxsize=None)
def primitive_safe_residues(speeds: tuple[int, ...], q: int) -> tuple[int, ...]:
    return tuple(
        a for a in range(1, q)
        if gcd(a, q) == 1 and closed_safe_at_phase(speeds, q, a)
    )


def primitive_count_deck(packet, qmax: int) -> tuple[int, ...]:
    return tuple(len(primitive_safe_residues(packet.speeds, q)) for q in range(2, qmax + 1))


def first_primitive_safe_q(packet, qmax: int) -> int | None:
    for q in range(2, qmax + 1):
        if primitive_safe_residues(packet.speeds, q):
            return q
    return None


def mobius(n: int) -> int:
    x = n
    p = 2
    factors = 0
    while p * p <= x:
        if x % p == 0:
            x //= p
            factors += 1
            if x % p == 0:
                return 0
            while x % p == 0:
                x //= p
        p += 1
    if x > 1:
        factors += 1
    return -1 if factors % 2 else 1


def divisors(n: int) -> tuple[int, ...]:
    return tuple(d for d in range(1, n + 1) if n % d == 0)


@lru_cache(maxsize=None)
def ramanujan_sum(q: int, n: int) -> int:
    g = gcd(q, n)
    return sum(d * mobius(q // d) for d in divisors(g))


def ramanujan_trace_deck(packet, qmax: int) -> tuple[tuple[tuple[int, int], ...], ...]:
    deck = []
    for q in range(2, qmax + 1):
        deck.append(tuple(sorted(Counter(ramanujan_sum(q, v) for v in packet.speeds).items())))
    return tuple(deck)


@dataclass(frozen=True)
class SplitAudit:
    name: str
    groups: int
    mixed_route: int
    mixed_status: int
    max_mixed: int
    mixed_packets: int
    route_counts: Counter


def audit_split(name: str, packets: list, key_func) -> tuple[SplitAudit, list[list]]:
    groups = group_by(packets, key_func)
    route_mixed = mixed_route_groups(groups)
    route_counts = Counter(packet.route for rows in route_mixed for packet in rows)
    audit = SplitAudit(
        name=name,
        groups=len(groups),
        mixed_route=len(route_mixed),
        mixed_status=len(mixed_status_groups(groups)),
        max_mixed=max((len(rows) for rows in route_mixed), default=0),
        mixed_packets=sum(len(rows) for rows in route_mixed),
        route_counts=route_counts,
    )
    return audit, route_mixed


@dataclass(frozen=True)
class Carrier:
    name: str
    status_purity: int
    route_purity: int
    primitive_period: int
    harmonic_trace: int
    avoids_exact_m: int
    compression: int
    proof_cost: int

    @property
    def values(self) -> tuple[int, ...]:
        return (
            self.status_purity,
            self.route_purity,
            self.primitive_period,
            self.harmonic_trace,
            self.avoids_exact_m,
            self.compression,
            8 - self.proof_cost,
        )


CARRIERS = (
    Carrier("primitive_count_deck_2_13", 5, 5, 5, 3, 5, 4, 3),
    Carrier("first_safe_q_2_13", 5, 5, 4, 2, 5, 5, 2),
    Carrier("ramanujan_trace_deck_2_14", 5, 5, 3, 5, 5, 3, 4),
    Carrier("coarse_et_unit_status_gate", 5, 1, 2, 2, 5, 5, 2),
    Carrier("exact_magnitude_cocycle", 5, 5, 2, 1, 0, 2, 4),
    Carrier("raw_residue_terminal_word", 1, 0, 1, 0, 5, 5, 1),
)


def carrier_tournament() -> dict[str, object]:
    names = [carrier.name for carrier in CARRIERS]
    order = {name: i for i, name in enumerate(names)}
    edges: set[tuple[str, str]] = set()
    score = Counter()

    for a, b in combinations(CARRIERS, 2):
        aw = sum(x > y for x, y in zip(a.values, b.values))
        bw = sum(x < y for x, y in zip(a.values, b.values))
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

    n = len(names)
    edge_idx = {(names.index(a), names.index(b)) for a, b in edges}
    hp = 0
    first_hp: tuple[str, ...] | None = None
    for perm in permutations(range(n)):
        if all((perm[i], perm[i + 1]) in edge_idx for i in range(n - 1)):
            hp += 1
            if first_hp is None:
                first_hp = tuple(names[i] for i in perm)

    ranking = sorted(names, key=lambda name: (-score[name], order[name]))
    return {
        "score_hist": dict(sorted(Counter(score[name] for name in names).items())),
        "directed_3cycles": c3,
        "hamiltonian_path_count": hp,
        "first_hamiltonian_path": first_hp,
        "score_order": ranking,
    }


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge")
    print("  considered vertices:")
    print("    runners, gaps, residue words, primitive denominator charts,")
    print("    Ramanujan c_q modes, safe phase residues, coarse zipper fibers,")
    print("    exact magnitude atoms, safe components, topology gates, and proof")
    print("    obligations.")
    print("  chosen vertices:")
    print("    proof carriers: coarse ET+unit gate, first primitive safe q,")
    print("    primitive count deck, Ramanujan trace deck, exact magnitude, and")
    print("    raw residue terminal word.")
    print("  preserved LRC predicate:")
    print("    boundary/open status from HYP-3030, plus q<=13 closed-safe witness")
    print("    availability for scheduling Q-WITNESS versus COVERING-MOMENT rows.")
    print("  destroyed information:")
    print("    exact M, exact safe interval lengths, full safe-component barcode,")
    print("    and arc-Cech topology unless those coordinates are reattached.")
    print()


def print_split_audit(audit: SplitAudit) -> None:
    print(
        f"  {audit.name:<32} fibers={audit.groups:<3d} "
        f"mixed_route={audit.mixed_route:<2d} mixed_status={audit.mixed_status:<2d} "
        f"max_mixed={audit.max_mixed:<2d} mixed_packets={audit.mixed_packets:<2d} "
        f"mixed_routes={dict(audit.route_counts)}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--single-limit", type=int, default=180)
    parser.add_argument("--two-swap-limit", type=int, default=36)
    parser.add_argument("--alias-depth", type=int, default=4)
    parser.add_argument("--lcm-tail-max", type=int, default=5)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--s194-result", type=Path, default=S194_RESULT)
    args = parser.parse_args()

    names = parse_s194_residual_names(args.s194_result)
    packets = packets_for_names(
        names,
        args.single_limit,
        args.two_swap_limit,
        args.alias_depth,
        args.lcm_tail_max,
        args.workers,
    )

    audits: list[SplitAudit] = []
    residual_groups: list[list] = []
    specs = (
        ("coarse_et_unit_gate", coarse_key),
        ("coarse_plus_first_q_2_13", lambda p: (coarse_key(p), first_primitive_safe_q(p, 13))),
        ("coarse_plus_first_q_2_14", lambda p: (coarse_key(p), first_primitive_safe_q(p, 14))),
        ("coarse_plus_primitive_deck_2_13", lambda p: (coarse_key(p), primitive_count_deck(p, 13))),
        ("coarse_plus_primitive_deck_2_14", lambda p: (coarse_key(p), primitive_count_deck(p, 14))),
        ("coarse_plus_ramanujan_trace_2_14", lambda p: (coarse_key(p), ramanujan_trace_deck(p, 14))),
    )
    for name, key_func in specs:
        audit, mixed = audit_split(name, packets, key_func)
        audits.append(audit)
        if name == "coarse_et_unit_gate":
            residual_groups = mixed

    print("=== LRC14 Ramanujan route scheduler S200 ===")
    print(
        "source=S194 stored residual rows from "
        f"{args.s194_result.relative_to(REPO)}"
    )
    print(
        "bank=HYP-2963 default row generator "
        f"single_limit={args.single_limit} two_swap_limit={args.two_swap_limit} "
        f"alias_depth={args.alias_depth} lcm_tail_max={args.lcm_tail_max}"
    )
    print(f"residual_packets={len(packets)} parsed_names={len(names)}")
    print()
    print_assumption_challenge()

    print("[1] Scheduler split audit on S194 residual route fibers")
    for audit in audits:
        print_split_audit(audit)
    print()

    print("[2] Coarse residual fibers with primitive-period labels")
    for index, rows in enumerate(sorted(residual_groups, key=lambda rows: (-len(rows), rows[0].name))):
        rows = sorted(rows, key=lambda p: (p.route, p.name))
        print(
            f"  coarse_route_mixed_fiber[{index}] size={len(rows)} "
            f"routes={dict(Counter(p.route for p in rows))}"
        )
        for packet in rows:
            first13 = first_primitive_safe_q(packet, 13)
            first14 = first_primitive_safe_q(packet, 14)
            deck13 = primitive_count_deck(packet, 13)
            q14_count = len(primitive_safe_residues(packet.speeds, 14))
            print(
                f"    {packet.name:<34} route={packet.route:<16} "
                f"M={fmt(packet.M):>7} q0={packet.q_threshold:<2d} "
                f"first13={str(first13):>4s} first14={str(first14):>4s} "
                f"deck13={deck13} q14_count={q14_count}"
            )
    print()

    print("[3] Deck readout")
    q_witness = [p for p in packets if p.route == "Q-WITNESS"]
    covering = [p for p in packets if p.route == "COVERING-MOMENT"]
    q_decks = Counter(primitive_count_deck(p, 13) for p in q_witness)
    covering_decks = Counter(primitive_count_deck(p, 13) for p in covering)
    print(f"  q_witness_packets={len(q_witness)} unique_deck13={len(q_decks)}")
    print(f"  covering_packets={len(covering)} unique_deck13={len(covering_decks)}")
    print("  q_witness_deck13_counts=" + str(dict(q_decks)))
    print("  covering_deck13_counts=" + str(dict(covering_decks)))
    print("  theorem_hint=Q-WITNESS iff deck13 has positive mass on some primitive period <=13")
    print("  guardrail=q=14 mass appears for many covering rows and is not a q<=13 theorem route")
    print()

    print("[4] Tournament Analysis over scheduler carriers")
    fp = carrier_tournament()
    print("  vertices_are=proof carriers and exact-period decks, not runners")
    print("  score_hist=" + str(fp["score_hist"]))
    print("  directed_3cycles=" + str(fp["directed_3cycles"]))
    print("  hamiltonian_path_count=" + str(fp["hamiltonian_path_count"]))
    print("  first_hamiltonian_path=" + " > ".join(fp["first_hamiltonian_path"]))
    print("  score_order=" + " > ".join(fp["score_order"]))
    print()

    print("[5] Proof readout")
    print("  1. S194 reduced the full-bank route debt to 15 strict-open")
    print("     coarse ET+unit fibers containing 38 packets.")
    print("  2. On those same packets, adding first primitive safe denominator")
    print("     for 2<=q<=13 gives 30 fibers, 0 mixed route, 0 mixed status.")
    print("  3. The full primitive count deck for 2<=q<=13 gives the same")
    print("     route-pure scheduler without exact M, safe interval length, or")
    print("     topology.")
    print("  4. Ramanujan trace profiles c_q(v), 2<=q<=14, also separate the")
    print("     residual route fibers, but they are a diagnostic unless paired")
    print("     with the safe-phase inequality.")
    print("  5. Candidate theorem target: after the coarse ET+unit status gate,")
    print("     every residual Q-WITNESS row has positive primitive safe mass")
    print("     on some q<=13, while every residual COVERING-MOMENT row has zero")
    print("     such mass and must be routed by q=14/covering or by a later")
    print("     boundary-moment certificate.")


if __name__ == "__main__":
    main()
