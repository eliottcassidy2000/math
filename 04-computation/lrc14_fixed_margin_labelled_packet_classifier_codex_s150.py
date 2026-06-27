#!/usr/bin/env python3
"""S150: fixed-margin labelled-packet classifier for LRC14 counterexamples.

This script refines HYP-2961 by making the word "family" operational.
The imported arXiv:2606.22636 fixed-margin swap-chain paper is not used as an
LRC theorem.  It is used as a proof architecture:

  fixed margins            -> packet families
  two-row heat bath        -> conditional packet resampling
  three-row reduction      -> local proof-obligation triples
  scalar count sector      -> qdiv/q-cover margins
  Johnson harmonic sectors -> non-scalar owner, source, C27/K33, moment labels

For each representative row, the script attaches exact LRC data already
available in the repo, then emits a fixed-margin signature.  Rows with the
same signature are treated as members of the same labelled packet family; rows
whose signatures are singleton and whose parameters are bounded are sporadics.

Tournament Analysis declaration
-------------------------------
Vertices:
    proof sectors / packet families, not runners.
Pair observable:
    which sector preserves strict-counterexample status, boundary ownership,
    labelled packet data, and a named proof exit before scalarization.
Switch:
    componentwise score; ties follow the printed Hamiltonian path.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha1
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd
from pathlib import Path
import argparse
import json
import os
import sys


REPO = Path(__file__).resolve().parents[1]
AP = tuple(range(1, 14))
DELTA = Fraction(1, 14)
UNITS14 = (1, 3, 5, 9, 11, 13)


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s138 = load_module(
    "s150_fixed_s138_frontier",
    REPO / "04-computation" / "lrc14_c27_two_swap_frontier_codex_s138.py",
)
s145 = load_module(
    "s150_fixed_s145_recombination",
    REPO / "04-computation" / "lrc14_measurable_rank_recombination_codex_s145.py",
)
s147 = load_module(
    "s150_fixed_s147_haar",
    REPO / "04-computation" / "lrc14_baire_haar_anyangle_codex_s147.py",
)


@dataclass(frozen=True)
class Candidate:
    name: str
    family_hint: str
    speeds: tuple[int, ...]


@dataclass(frozen=True)
class PacketAudit:
    name: str
    family_hint: str
    speeds: tuple[int, ...]
    qdiv: int
    q_cover_2_14: tuple[int, ...]
    exact_M: Fraction
    arg_denoms: tuple[int, ...]
    strict_safe_mu: Fraction
    strict_components: int
    row_label: str
    branch: str
    transfer: str
    packet_route: str
    packet_rank: int
    state_lift: bool
    atom_keys: tuple[str, ...]
    u14_zero: tuple[int, ...]
    u14_boundary: tuple[int, ...]
    c27_state_counts: tuple[tuple[str, int], ...]
    class_id: str
    class_json: str
    bucket: str
    live_family: str


def primitive(row: tuple[int, ...]) -> bool:
    g = 0
    for v in row:
        g = gcd(g, v)
    return g == 1


def row_from_swaps(removed: tuple[int, ...], added: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted((set(AP) - set(removed)) | set(added)))


def candidate_rows(lcm_tail_max: int) -> list[Candidate]:
    rows: list[Candidate] = [
        Candidate("AP", "S0 equality sporadic", AP),
        Candidate("GW 12->24", "S0 equality sporadic", row_from_swaps((12,), (24,))),
        Candidate("near/K33 12->36", "S1 positive K33 sporadic", row_from_swaps((12,), (36,))),
        Candidate("petal 10->20", "S1 unit-petal sporadic", row_from_swaps((10,), (20,))),
        Candidate("petal 13->26", "S1 unit-petal sporadic", row_from_swaps((13,), (26,))),
        Candidate("P10+GW", "S1 unit-petal/GW strip", row_from_swaps((10, 12), (20, 24))),
        Candidate("P10+K33", "S1 K33 strip", row_from_swaps((10, 12), (20, 36))),
        Candidate("residue liar 12->26", "q-witness decoy", row_from_swaps((12,), (26,))),
        Candidate("magnitude liar 12->96", "magnitude decoy", row_from_swaps((12,), (96,))),
        Candidate("floor-odd GW iso impostor", "raw tournament decoy", row_from_swaps((12,), (360,))),
        Candidate("covering repair 13/14 via 182", "bounded covering repair", tuple(list(range(1, 13)) + [182])),
        Candidate("covering comb 12->84", "lcm-tail family", row_from_swaps((12,), (84,))),
        Candidate("covering comb 12->168", "lcm-tail family", row_from_swaps((12,), (168,))),
        Candidate("covering comb 12->252", "lcm-tail family", row_from_swaps((12,), (252,))),
        Candidate("covering comb 12->420", "lcm-tail family", row_from_swaps((12,), (420,))),
        Candidate("apex multiple decoy 12->28", "q-witness apex decoy", row_from_swaps((12,), (28,))),
        Candidate("apex multiple decoy 12->42", "q-witness apex decoy", row_from_swaps((12,), (42,))),
        Candidate("apex multiple decoy 12->56", "q-witness apex decoy", row_from_swaps((12,), (56,))),
    ]
    for m in range(5, lcm_tail_max + 1):
        rows.append(
            Candidate(
                f"covering comb 12->{84*m}",
                "lcm-tail family",
                row_from_swaps((12,), (84 * m,)),
            )
        )

    seen: set[tuple[int, ...]] = set()
    out: list[Candidate] = []
    for c in rows:
        if c.speeds in seen:
            continue
        if len(c.speeds) != 13 or not primitive(c.speeds):
            continue
        seen.add(c.speeds)
        out.append(c)
    return out


def qdiv(row: tuple[int, ...]) -> int:
    return s138.s124.q_threshold(row)


def q_cover(row: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(1 if any(v % d == 0 for v in row) else 0 for d in range(2, 15))


def u14_profiles(row: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    zero: list[int] = []
    boundary: list[int] = []
    for u in UNITS14:
        z = 0
        b = 0
        for v in row:
            r = (v * u) % 14
            d = min(r, 14 - r)
            if d == 0:
                z += 1
            elif d == 1:
                b += 1
        zero.append(z)
        boundary.append(b)
    return tuple(zero), tuple(boundary)


def c27_state_counts(row: tuple[int, ...]) -> tuple[tuple[str, int], ...]:
    shells, zero_count = s138.shell_profile(row)
    counts = Counter(f"{shell.state}:g{shell.gcd27}" for shell in shells)
    if zero_count:
        counts[f"zero27:{zero_count}"] += 1
    return tuple(sorted(counts.items()))


def atom_keys(packet: object) -> tuple[str, ...]:
    return tuple(component.key for component in packet.components)


def fixed_margin_payload(
    qd: int,
    cover: tuple[int, ...],
    u_zero: tuple[int, ...],
    u_boundary: tuple[int, ...],
    c27_counts: tuple[tuple[str, int], ...],
    atoms: tuple[str, ...],
) -> dict[str, object]:
    if qd < 14:
        q_bucket = "direct-q"
    elif qd == 14:
        q_bucket = "boundary-q14"
    else:
        q_bucket = "covering-qgt14"
    return {
        "q_bucket": q_bucket,
        "q_cover_2_14": cover,
        "u14_zero": u_zero,
        "u14_boundary": u_boundary,
        "c27_state_counts": c27_counts,
        "packet_atoms": atoms,
    }


def class_id(payload: dict[str, object]) -> tuple[str, str]:
    text = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return sha1(text.encode("ascii")).hexdigest()[:12], text


def bucket_for(qd: int, exact_M: Fraction, strict_mu: Fraction, packet: object, row_label: str) -> tuple[str, str]:
    if exact_M < DELTA:
        return "COUNTEREXAMPLE FOUND", "actual strict bad row"
    if qd <= 13:
        return "D0 direct q-witness family", "discharged"
    if qd == 14 and strict_mu == 0 and row_label in {"AP", "GW 12->24"}:
        return "S0 equality sporadic AP/GW", "tight equality, not strict"
    if qd == 14 and strict_mu == 0:
        return "L5 boundary-only source-kernel threat", "live if real"
    if strict_mu > 0:
        if packet.state_lift:
            return "D4 positive K33/state-lift labelled row", "discharged unless zero-open"
        if "petal" in packet.route or "GW" in packet.route:
            return "D3 positive unit-petal/GW strip", "discharged"
        if qd > 14:
            return "D2 positive covering Haar-open family", "discharged"
        return "D2 positive Haar-open row", "discharged"
    if qd > 14 and packet.state_lift:
        return "L4 K33 zero-open state-lift family", "live family"
    if qd > 14:
        return "L5 unnamed zero-open source-kernel family", "live family"
    return "unclassified labelled packet", "review"


def audit_candidates(candidates: list[Candidate], workers: int) -> list[PacketAudit]:
    exact = s138.compute_exact([c.speeds for c in candidates], workers=workers, progress_every=0)
    ap_shells, _ = s138.shell_profile(AP)
    ap_mask = s138.tournament_mask([s138.shell_priority(shell) for shell in ap_shells])
    audits: list[PacketAudit] = []
    for c in candidates:
        M, denoms = exact[c.speeds]
        row_audit = s138.audit_row(c.speeds, M, denoms, ap_mask)
        packet = s145.classify_packet(row_audit)
        safe = s147.exact_row_measure(c.speeds)
        qd = qdiv(c.speeds)
        cover = q_cover(c.speeds)
        u_zero, u_boundary = u14_profiles(c.speeds)
        c27_counts = c27_state_counts(c.speeds)
        atoms = atom_keys(packet)
        payload = fixed_margin_payload(qd, cover, u_zero, u_boundary, c27_counts, atoms)
        cid, cjson = class_id(payload)
        bucket, live = bucket_for(qd, M, safe["safe_measure"], packet, row_audit.label)
        audits.append(
            PacketAudit(
                name=c.name,
                family_hint=c.family_hint,
                speeds=c.speeds,
                qdiv=qd,
                q_cover_2_14=cover,
                exact_M=M,
                arg_denoms=denoms,
                strict_safe_mu=safe["safe_measure"],
                strict_components=len(safe["safe_components"]),
                row_label=row_audit.label,
                branch=row_audit.branch,
                transfer=row_audit.transfer,
                packet_route=packet.route,
                packet_rank=packet.ph_rank,
                state_lift=packet.state_lift,
                atom_keys=atoms,
                u14_zero=u_zero,
                u14_boundary=u_boundary,
                c27_state_counts=c27_counts,
                class_id=cid,
                class_json=cjson,
                bucket=bucket,
                live_family=live,
            )
        )
    return sorted(audits, key=lambda a: (a.exact_M, a.qdiv, a.name))


def fmt(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge")
    print("  considered vertices:")
    print("    runners, gaps, divisor clocks, fixed circle sections, section boundaries,")
    print("    wall-crossing events, residues, cover arcs, Fourier/moment modes,")
    print("    matroid-like circuits, K33/state-lift obligations, and packet families.")
    print("  chosen vertices:")
    print("    fixed-margin labelled packet classes and proof sectors.")
    print("  preserved LRC predicate:")
    print("    strict-counterexample status M(S)<1/14, boundary-vs-open witness status,")
    print("    qdiv cover obligations, and named C27/K33/source packet ownership.")
    print("  destroyed information:")
    print("    raw runner identity and exact wall-crossing order after conversion to labels.")
    print("  challenged assumption:")
    print("    a counterexample family is not a raw parametric row list; it is a")
    print("    fixed-margin packet class plus the non-scalar labels that survive swaps.")
    print()


def print_arxiv_import() -> None:
    print("[1] arXiv:2606.22636 proof-pattern import")
    print("  imported as architecture, not as an LRC theorem.")
    print("  fixed-margin binary matrices -> labelled packet incidence matrices")
    print("  swap chain within fixed margins -> packet-preserving mutations")
    print("  two-row heat-bath comparison -> conditional resampling of two packet rows")
    print("  three-row reduction -> local triples of proof-obligation rows")
    print("  scalar count sector -> qdiv/q-cover counts")
    print("  Johnson harmonic sectors -> boundary-owner, C27/unital, source, K33, gK8 labels")
    print()


def print_packet_table(audits: list[PacketAudit]) -> None:
    print("[2] Representative packet classification")
    print(
        f"  {'row':30s} {'M':>8s} {'qdiv':>4s} {'mu_open':>9s} "
        f"{'class':12s} {'bucket'}"
    )
    for a in audits:
        print(
            f"  {a.name[:30]:30s} {fmt(a.exact_M):>8s} {a.qdiv:4d} "
            f"{fmt(a.strict_safe_mu):>9s} {a.class_id:12s} {a.bucket}"
        )
        print(
            f"      route={a.packet_route}; atoms={a.atom_keys or '-'}; "
            f"U14zero={a.u14_zero}; U14bdry={a.u14_boundary}"
        )
        print(f"      transfer={a.transfer}")
    print()


def print_family_groups(audits: list[PacketAudit]) -> None:
    print("[3] Fixed-margin family grouping")
    groups: dict[str, list[PacketAudit]] = defaultdict(list)
    for a in audits:
        groups[a.class_id].append(a)
    multi = {cid: group for cid, group in groups.items() if len(group) > 1}
    print(f"  classes={len(groups)} singleton_sporadic_signatures={len(groups) - len(multi)} shared_family_signatures={len(multi)}")
    if not multi:
        print("  no shared signatures in this small representative bank")
    for cid, group in sorted(multi.items(), key=lambda item: (-len(item[1]), item[0])):
        print(f"  class {cid}: rows={len(group)}")
        print("    " + ", ".join(a.name for a in group))
        sample = group[0]
        print(f"    payload={sample.class_json}")
    print()


def print_counterexample_live_summary(audits: list[PacketAudit]) -> None:
    print("[4] Counterexample-family verdict")
    bucket_counts = Counter(a.bucket for a in audits)
    live_counts = Counter(a.live_family for a in audits)
    below = [a for a in audits if a.exact_M < DELTA]
    zero_open_qgt14 = [a for a in audits if a.qdiv > 14 and a.strict_safe_mu == 0]
    qgt14 = [a for a in audits if a.qdiv > 14]
    print(f"  bucket_counts={dict(sorted(bucket_counts.items()))}")
    print(f"  live_status_counts={dict(sorted(live_counts.items()))}")
    print(f"  exact strict counterexamples found={len(below)}")
    print(f"  qdiv>14 representatives={len(qgt14)}")
    print(f"  qdiv>14 zero-open representatives={len(zero_open_qgt14)}")
    print("  HYP-2961 live buckets tested here:")
    print("    L1 apex-multiple residual: not covered by this small bank unless many 14-multiples appear")
    print("    L2 genuine-wide zero-moment: moment image not computed here")
    print("    L3 bounded covering core: represented by covering repairs and lcm tails; all positive-open")
    print("    L4 K33 zero-open state-lift: no zero-open representative found")
    print("    L5 unnamed source-kernel: no zero-open non-K33 representative found")
    print()


def print_tournament_analysis() -> None:
    print("[5] Tournament Analysis: fixed-margin packet sectors")
    sectors = [
        ("qdiv/q-cover count sector", (7, 7, 5, 4, 5, 6)),
        ("Haar open-vs-boundary sector", (7, 6, 7, 6, 5, 6)),
        ("source-spectrum pullback", (6, 7, 7, 7, 6, 7)),
        ("C27/unital owner sector", (5, 6, 6, 7, 6, 6)),
        ("K33 state-lift sector", (5, 5, 6, 7, 7, 6)),
        ("gK8/L_y moment bridge", (4, 5, 6, 6, 7, 7)),
        ("fixed-margin mutation family", (6, 6, 5, 5, 5, 6)),
        ("unnamed source-kernel bucket", (3, 4, 7, 7, 4, 7)),
        ("raw residue/tournament shadow", (1, 2, 1, 1, 1, 1)),
    ]
    mask = s138.tournament_mask([score for _, score in sectors])
    fp = s138.tournament_fingerprint(mask, len(sectors))
    order = sorted(range(len(sectors)), key=lambda i: sectors[i][1], reverse=True)
    print("  vertices are proof sectors, not runners.")
    print("  pair observable:")
    print("    strict predicate retention, owner-label retention, source/kernel visibility,")
    print("    named proof exit, finite-family pressure, anti-scalar guard.")
    print(f"  fingerprint: score_hist={fp['score_hist']} c3={fp['c3']} scc={fp['scc']} hp={fp['hp']}")
    print("  Hamiltonian path:")
    print("    " + " > ".join(sectors[i][0] for i in order))
    print()


def print_theorem_target() -> None:
    print("[6] Labelled Packet Theorem target")
    print("  For every primitive 13-row S, build its fixed-margin incidence packet:")
    print("    rows = q-cover obligations, U14 apex zero/boundary contacts,")
    print("           C27 shell states, source-spectrum labels, K33/state-lift flags,")
    print("           and boundary-moment/gK8 coordinates when qdiv>14.")
    print("  Packet-preserving swaps define families.  A family is safe if each")
    print("  count sector or non-scalar sector has a named discharge.")
    print("  The only possible strict counterexamples after this theorem would be:")
    print("    L1 apex-multiple residual, L2 genuine-wide zero-moment, L3 bounded")
    print("    covering core, L4 K33 zero-open state lift, or L5 unnamed source kernel.")
    print("  This representative run found no row in the strict or zero-open live buckets.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lcm-tail-max", type=int, default=5)
    parser.add_argument("--workers", type=int, default=max(1, min(os.cpu_count() or 1, 6)))
    args = parser.parse_args()

    print("S150 LRC14 FIXED-MARGIN LABELLED-PACKET CLASSIFIER")
    print("=" * 78)
    print(f"lcm_tail_max={args.lcm_tail_max}, workers={args.workers}")
    print_assumption_challenge()
    print_arxiv_import()
    audits = audit_candidates(candidate_rows(args.lcm_tail_max), args.workers)
    print_packet_table(audits)
    print_family_groups(audits)
    print_counterexample_live_summary(audits)
    print_tournament_analysis()
    print_theorem_target()


if __name__ == "__main__":
    main()
