#!/usr/bin/env python3
"""Automatic-fiber zipper splitter for the LRC14 packet bank.

This is a creative proof-interface scout, not a proof of LRC14.  It follows
HYP-3016 and HYP-3017 by asking which coordinates must be zipped onto a
Moser/fibbinary automatic word before its fibers stop mixing theorem routes.

Tournament Analysis declaration:
  vertices: quotient / sidecar bundles, not runners;
  pairwise observable: route purity on the HYP-2963 bank, retained magnitude
    cocycle, safe-topology data, packet-label data, finite-state checkability,
    and proof cost;
  switch/gauge: majority comparison of the observable vector;
  tie Hamiltonian path: the listed zipper-teeth order.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
from pathlib import Path
import argparse
import sys


REPO = Path(__file__).resolve().parents[1]
TARGET_WORD = "MFCMMCCFFFCCC"


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


lp = load_module(
    "lrc14_labelled_packet_classifier_s184",
    REPO / "04-computation" / "lrc14_labelled_packet_counterexample_classifier_codex_20260624.py",
)


def fmt(fr: Fraction) -> str:
    return str(fr.numerator) if fr.denominator == 1 else f"{fr.numerator}/{fr.denominator}"


def residue_mfc_pairs(packet) -> tuple[tuple[int, str], ...]:
    return tuple(sorted((v % 14, lp.automatic_letter(v)) for v in packet.speeds))


def terminal_word(packet) -> tuple[str, ...]:
    return tuple(
        (
            "Fok" if lp.is_fibbinary(v) else "Fdead",
            "Mok" if lp.is_moser(v) else "Mdead",
        )
        for v in packet.speeds
    )


def magnitude_cocycle(packet) -> tuple[Fraction, int, int, Fraction]:
    return packet.M, packet.q_threshold, packet.farey_excess, packet.lacunary_tail_ratio


def barcode_shadow(packet) -> tuple[str, int, Fraction, Fraction, int]:
    status = "boundary" if packet.strict_safe_mu == 0 else "open"
    return (
        status,
        packet.strict_components,
        packet.strict_safe_mu,
        packet.M - lp.THRESHOLD,
        packet.boundary_count,
    )


def local_packet_shadow(packet) -> tuple[str, str, bool, str, str, str]:
    return (
        packet.transfer,
        packet.packet_route,
        packet.state_lift,
        packet.source_family,
        packet.power_lift_guard,
        packet.q_factorization,
    )


@dataclass(frozen=True)
class Split:
    name: str
    description: str
    key_func: object
    finite_state: int
    magnitude: int
    topology: int
    packet_label: int
    proof_cost: int


SPLITS: tuple[Split, ...] = (
    Split(
        "automatic_word",
        "Moser/fibbinary/carry word only",
        lambda p: (p.automatic_word,),
        5,
        0,
        0,
        0,
        1,
    ),
    Split(
        "residue_terminal_fiber",
        "automatic word plus residue and DFA terminal states",
        lambda p: (p.automatic_word, residue_mfc_pairs(p), terminal_word(p)),
        5,
        1,
        0,
        1,
        2,
    ),
    Split(
        "magnitude_cocycle",
        "adds exact M, q-threshold, Farey excess, and tail ratio",
        lambda p: (p.automatic_word, residue_mfc_pairs(p), magnitude_cocycle(p)),
        4,
        5,
        0,
        2,
        3,
    ),
    Split(
        "barcode_shadow",
        "adds open/boundary status, component count, safe mass, margin, boundary units",
        lambda p: (
            p.automatic_word,
            residue_mfc_pairs(p),
            magnitude_cocycle(p),
            barcode_shadow(p),
        ),
        3,
        5,
        5,
        2,
        4,
    ),
    Split(
        "packet_zipper",
        "adds C27/K33/source/power/factorization local packet labels",
        lambda p: (
            p.automatic_word,
            residue_mfc_pairs(p),
            magnitude_cocycle(p),
            barcode_shadow(p),
            local_packet_shadow(p),
        ),
        3,
        5,
        5,
        5,
        5,
    ),
    Split(
        "route_labelled_packet",
        "adds the theorem route label; this is the current classifier endpoint",
        lambda p: (
            p.automatic_word,
            residue_mfc_pairs(p),
            magnitude_cocycle(p),
            barcode_shadow(p),
            local_packet_shadow(p),
            p.route,
        ),
        2,
        5,
        5,
        6,
        6,
    ),
)


@dataclass(frozen=True)
class SplitReport:
    split: Split
    fibers: int
    mixed_route_fibers: int
    mixed_family_fibers: int
    max_fiber_size: int
    max_mixed_size: int
    largest_mixed_key: object
    largest_mixed_routes: tuple[tuple[str, int], ...]
    largest_mixed_families: tuple[tuple[str, int], ...]

    @property
    def route_purity_score(self) -> int:
        if self.fibers == 0:
            return 0
        return (1000 * (self.fibers - self.mixed_route_fibers)) // self.fibers

    @property
    def vector(self) -> tuple[int, ...]:
        return (
            self.route_purity_score,
            -self.max_mixed_size,
            self.split.magnitude,
            self.split.topology,
            self.split.packet_label,
            self.split.finite_state,
            7 - self.split.proof_cost,
        )


def split_report(packets: list, split: Split) -> SplitReport:
    groups: dict[object, list] = defaultdict(list)
    for packet in packets:
        groups[split.key_func(packet)].append(packet)

    mixed_route = []
    mixed_family = []
    for key, rows in groups.items():
        if len({p.route for p in rows}) > 1:
            mixed_route.append((key, rows))
        if len({p.family for p in rows}) > 1:
            mixed_family.append((key, rows))

    if mixed_route:
        largest_key, largest_rows = max(mixed_route, key=lambda item: (len(item[1]), str(item[0])))
    else:
        largest_key, largest_rows = None, []

    return SplitReport(
        split=split,
        fibers=len(groups),
        mixed_route_fibers=len(mixed_route),
        mixed_family_fibers=len(mixed_family),
        max_fiber_size=max((len(rows) for rows in groups.values()), default=0),
        max_mixed_size=max((len(rows) for _, rows in mixed_route), default=0),
        largest_mixed_key=largest_key,
        largest_mixed_routes=tuple(Counter(p.route for p in largest_rows).most_common()),
        largest_mixed_families=tuple(Counter(p.family for p in largest_rows).most_common()),
    )


def orient(a: SplitReport, b: SplitReport) -> str:
    av = a.vector
    bv = b.vector
    awins = sum(x > y for x, y in zip(av, bv))
    bwins = sum(x < y for x, y in zip(av, bv))
    if awins != bwins:
        return a.split.name if awins > bwins else b.split.name
    if av != bv:
        return a.split.name if av > bv else b.split.name
    names = [s.name for s in SPLITS]
    return a.split.name if names.index(a.split.name) < names.index(b.split.name) else b.split.name


def tournament_fingerprint(reports: list[SplitReport]) -> tuple[Counter, int, list[int], int, list[str]]:
    names = [r.split.name for r in reports]
    edges: dict[tuple[int, int], bool] = {}
    score = Counter()
    for i, j in combinations(range(len(reports)), 2):
        winner = orient(reports[i], reports[j])
        i_wins = winner == names[i]
        edges[(i, j)] = i_wins
        score[i if i_wins else j] += 1

    score_hist = Counter(score[i] for i in range(len(reports)))
    c3 = 0
    for a, b, c in combinations(range(len(reports)), 3):
        ab = edges[(a, b)]
        ac = edges[(a, c)]
        bc = edges[(b, c)]
        # a->b, b->c, c->a or reverse.
        if (ab and bc and not ac) or ((not ab) and (not bc) and ac):
            c3 += 1

    mat = [[False] * len(reports) for _ in reports]
    for i, j in combinations(range(len(reports)), 2):
        if edges[(i, j)]:
            mat[i][j] = True
        else:
            mat[j][i] = True

    def reach(start: int, graph: list[list[bool]]) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            u = stack.pop()
            for v, ok in enumerate(graph[u]):
                if ok and v not in seen:
                    seen.add(v)
                    stack.append(v)
        return seen

    reverse = [[mat[j][i] for j in range(len(mat))] for i in range(len(mat))]
    remaining = set(range(len(mat)))
    scc_sizes: list[int] = []
    while remaining:
        seed = min(remaining)
        comp = reach(seed, mat) & reach(seed, reverse)
        scc_sizes.append(len(comp))
        remaining -= comp

    hp = 0
    for order in permutations(range(len(mat))):
        if all(mat[order[i]][order[i + 1]] for i in range(len(order) - 1)):
            hp += 1

    ranking = sorted(range(len(reports)), key=lambda i: (-score[i], reports[i].split.name))
    return score_hist, c3, sorted(scc_sizes, reverse=True), hp, [names[i] for i in ranking]


def sample_rows(rows: list, limit: int = 8) -> list:
    ordered = sorted(rows, key=lambda p: (p.route, p.M, p.strict_safe_mu, p.name))
    if len(ordered) <= limit:
        return ordered
    return ordered[: limit // 2] + ordered[-(limit // 2) :]


def print_target_word_templates(packets: list) -> None:
    rows = [p for p in packets if p.automatic_word == TARGET_WORD]
    print(f"[3] Target automatic word {TARGET_WORD}")
    print(f"  rows={len(rows)}")
    print(f"  routes={dict(Counter(p.route for p in rows))}")
    print(f"  families={dict(Counter(p.family for p in rows))}")
    print(f"  distinct M values={len({p.M for p in rows})}")
    print(f"  distinct residue-MFC fibers={len({residue_mfc_pairs(p) for p in rows})}")
    print()
    for split in SPLITS[1:5]:
        groups: dict[object, list] = defaultdict(list)
        for packet in rows:
            groups[split.key_func(packet)].append(packet)
        mixed = [
            packet_rows
            for packet_rows in groups.values()
            if len({p.route for p in packet_rows}) > 1
        ]
        print(f"  {split.name}: fibers={len(groups)} mixed_route={len(mixed)}")
        if mixed:
            largest = max(mixed, key=len)
            print(f"    largest_mixed_rows={len(largest)} routes={dict(Counter(p.route for p in largest))}")
            for packet in sample_rows(largest, 6):
                print(
                    "      {name:34s} M={M:>7s} mu={mu:>9s} comps={comp:2d} "
                    "q0={q0:2d} route={route}".format(
                        name=packet.name[:34],
                        M=fmt(packet.M),
                        mu=fmt(packet.strict_safe_mu),
                        comp=packet.strict_components,
                        q0=packet.q_threshold,
                        route=packet.route,
                    )
                )
    print()


def print_reports(reports: list[SplitReport]) -> None:
    print("[1] Zipper split reports")
    print(
        "  {name:24s} {fibers:>7s} {mixed:>7s} {mfam:>7s} {maxf:>7s} "
        "{maxm:>7s} {purity:>7s}"
    .format(
            name="split",
            fibers="fibers",
            mixed="mixR",
            mfam="mixF",
            maxf="maxF",
            maxm="maxMix",
            purity="purity",
        )
    )
    for report in reports:
        print(
            "  {name:24s} {fibers:7d} {mixed:7d} {mfam:7d} {maxf:7d} "
            "{maxm:7d} {purity:6.1f}%".format(
                name=report.split.name,
                fibers=report.fibers,
                mixed=report.mixed_route_fibers,
                mfam=report.mixed_family_fibers,
                maxf=report.max_fiber_size,
                maxm=report.max_mixed_size,
                purity=report.route_purity_score / 10,
            )
        )
    print()
    print("  largest mixed fiber by split:")
    for report in reports:
        print(f"    {report.split.name}: max_mixed={report.max_mixed_size} routes={dict(report.largest_mixed_routes)}")
    print()


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge")
    print("  considered vertices:")
    print("    runners, gaps, fixed circle sections, endpoint owners, residues,")
    print("    residue-automaton fibers, persistence bars, magnitude cocycles,")
    print("    C27/K33 packet labels, and proof obligations.")
    print("  chosen vertices:")
    print("    quotient/sidecar bundles over HYP-2963 packets.")
    print("  preserved LRC predicate:")
    print("    route-purity for the labelled theorem buckets, plus strict-open versus")
    print("    boundary equality at threshold 1/14.")
    print("  challenged assumption:")
    print("    an automatic word is not a proof coordinate until its mixed fibers are")
    print("    split by magnitude, safe topology, and packet labels.")
    print()


def print_tournament(reports: list[SplitReport]) -> None:
    score_hist, c3, scc, hp, ranking = tournament_fingerprint(reports)
    print("[4] Tournament Analysis")
    print("  vertices_are=quotient/sidecar bundles, not runners")
    print("  observable=route purity, max mixed-fiber size, magnitude retention,")
    print("             safe-topology retention, packet-label retention, finite-state")
    print("             checkability, and proof cost")
    print("  switch=majority comparison of observable vectors")
    print("  tie_hamiltonian_path=" + " > ".join(s.name for s in SPLITS))
    print(f"  score_hist={dict(sorted(score_hist.items()))}")
    print(f"  directed_3cycles={c3}")
    print(f"  scc_sizes={scc}")
    print(f"  hamiltonian_path_count={hp}")
    print("  score_order=" + " > ".join(ranking))
    print()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--single-limit", type=int, default=180)
    parser.add_argument("--two-swap-limit", type=int, default=36)
    parser.add_argument("--alias-depth", type=int, default=4)
    parser.add_argument("--lcm-tail-max", type=int, default=5)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()

    rows = lp.build_bank(args.single_limit, args.two_swap_limit, args.alias_depth, args.lcm_tail_max)
    packets = lp.compute_packets(rows, args.workers)
    reports = [split_report(packets, split) for split in SPLITS]

    print("=== LRC14 automatic fiber zipper splitter S187 ===")
    print(
        "bank=HYP-2963 default rows "
        f"single_limit={args.single_limit} two_swap_limit={args.two_swap_limit} "
        f"alias_depth={args.alias_depth} lcm_tail_max={args.lcm_tail_max}"
    )
    print(f"packets={len(packets)}")
    print()
    print_assumption_challenge()
    print_reports(reports)
    print("[2] Proof readout")
    print("  automatic_word is the side-channel warning layer; it has many mixed")
    print("  route fibers.")
    print("  residue_terminal_fiber sharply reduces but does not eliminate mixing.")
    print("  magnitude_cocycle is the first tested non-route coordinate with zero")
    print("  mixed theorem-route fibers on this bank.")
    print("  barcode_shadow and packet_zipper then become certificate-anchor")
    print("  refinements rather than the first purity repair.")
    print("  route_labelled_packet is therefore not a shortcut; it is the endpoint")
    print("  record that the proof should reconstruct by family lemmas.")
    print()
    print_target_word_templates(packets)
    print_tournament(reports)


if __name__ == "__main__":
    main()
