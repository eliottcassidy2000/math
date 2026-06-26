#!/usr/bin/env python3
"""LRC14 zipper-fiber convergence scout.

HYP-3023 found that the automatic/residue fiber zipper becomes route-pure only
after the magnitude cocycle is attached.  HYP-3020 suggested a trident:
Erdos-Turan discrepancy, height, and Hensel data.  This script splices those
ideas on the full HYP-2963 packet bank and asks a narrower question:

    which clocks actually shrink mixed zipper fibers, and which clocks are
    only warnings?

The Hensel coordinate is refined to a unit rule: roots in F_p^* are genuine
local unit clocks, while the forced zero root is recorded separately as scale
or nilpotent debt.  Tournament Analysis is over quotient/proof carriers, not
runners.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
from math import cos, pi, sin, sqrt
from pathlib import Path
import argparse
import sys


REPO = Path(__file__).resolve().parents[1]
TARGET_WORD = "MFCMMCCFFFCCC"
ET_QS = (14, 27, 41)
HENSEL_PS = (2, 3, 7)


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


lp = load_module(
    "lrc14_labelled_packet_classifier_s188",
    REPO / "04-computation" / "lrc14_labelled_packet_counterexample_classifier_codex_20260624.py",
)


def fmt(fr: Fraction) -> str:
    return str(fr.numerator) if fr.denominator == 1 else f"{fr.numerator}/{fr.denominator}"


def residue_mfc_pairs(packet) -> tuple[tuple[int, str], ...]:
    return tuple(sorted((v % 14, lp.automatic_letter(v)) for v in packet.speeds))


def terminal_word(packet) -> tuple[tuple[str, str], ...]:
    return tuple(
        (
            "Fok" if lp.is_fibbinary(v) else "Fdead",
            "Mok" if lp.is_moser(v) else "Mdead",
        )
        for v in packet.speeds
    )


def residue_terminal(packet):
    return packet.automatic_word, residue_mfc_pairs(packet), terminal_word(packet)


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


@lru_cache(maxsize=None)
def residue_l1_speeds(speeds: tuple[int, ...], q: int) -> Fraction:
    counts = Counter(v % q for v in speeds)
    n = len(speeds)
    return sum(abs(Fraction(counts[r], 1) - Fraction(n, q)) for r in range(q))


def residue_l1(packet, q: int) -> Fraction:
    return residue_l1_speeds(packet.speeds, q)


@lru_cache(maxsize=None)
def erdos_turan_proxy_speeds(speeds: tuple[int, ...], q: int, hmax: int = 6) -> tuple[int, int]:
    """Return a scaled ET sum and its dominant harmonic."""

    total = 0.0
    best_h = 0
    best_term = -1.0
    for h in range(1, hmax + 1):
        real = 0.0
        imag = 0.0
        for v in speeds:
            angle = 2.0 * pi * h * (v % q) / q
            real += cos(angle)
            imag += sin(angle)
        term = sqrt(real * real + imag * imag) / h
        total += term
        if term > best_term:
            best_term = term
            best_h = h
    return int(round(1000.0 * total)), best_h


def erdos_turan_proxy(packet, q: int, hmax: int = 6) -> tuple[int, int]:
    return erdos_turan_proxy_speeds(packet.speeds, q, hmax)


@lru_cache(maxsize=None)
def et_clock_speeds(speeds: tuple[int, ...]) -> tuple[tuple[Fraction, int, int], ...]:
    return tuple(
        (residue_l1_speeds(speeds, q), *erdos_turan_proxy_speeds(speeds, q))
        for q in ET_QS
    )


def et_clock(packet) -> tuple[tuple[Fraction, int, int], ...]:
    return et_clock_speeds(packet.speeds)


@lru_cache(maxsize=None)
def et_coarse_clock_speeds(speeds: tuple[int, ...]) -> tuple[tuple[int, int, int], ...]:
    """A deliberately coarse ET gate: L1 numerator buckets plus ET 500-bins."""

    out = []
    for q in ET_QS:
        l1 = residue_l1_speeds(speeds, q)
        et, h = erdos_turan_proxy_speeds(speeds, q)
        out.append((int(20 * l1), et // 500, h))
    return tuple(out)


def et_coarse_clock(packet) -> tuple[tuple[int, int, int], ...]:
    return et_coarse_clock_speeds(packet.speeds)


@lru_cache(maxsize=None)
def poly_value_mod_speeds(speeds: tuple[int, ...], x: int, p: int) -> tuple[int, int]:
    val = sum(pow(x, v, p) for v in speeds) % p
    der = sum((v % p) * pow(x, v - 1, p) for v in speeds) % p
    return val, der


def poly_value_mod(packet, x: int, p: int) -> tuple[int, int]:
    return poly_value_mod_speeds(packet.speeds, x, p)


@lru_cache(maxsize=None)
def hensel_unit_signature_speeds(speeds: tuple[int, ...], p: int) -> tuple[tuple[int, ...], tuple[int, ...], str]:
    """Roots in F_p^* plus a separate status for the forced zero root."""

    unit_roots = []
    singular_units = []
    for x in range(1, p):
        val, der = poly_value_mod_speeds(speeds, x, p)
        if val == 0:
            unit_roots.append(x)
            if der == 0:
                singular_units.append(x)

    zero_val, zero_der = poly_value_mod_speeds(speeds, 0, p)
    assert zero_val == 0
    zero_status = "zero-simple" if zero_der != 0 else "zero-singular"
    return tuple(unit_roots), tuple(singular_units), zero_status


def hensel_unit_signature(packet, p: int) -> tuple[tuple[int, ...], tuple[int, ...], str]:
    return hensel_unit_signature_speeds(packet.speeds, p)


@lru_cache(maxsize=None)
def hensel_unit_rule_speeds(speeds: tuple[int, ...]) -> tuple[tuple[int, tuple[int, ...], tuple[int, ...], str], ...]:
    return tuple(
        (p, *hensel_unit_signature_speeds(speeds, p))
        for p in HENSEL_PS
    )


def hensel_unit_rule(packet) -> tuple[tuple[int, tuple[int, ...], tuple[int, ...], str], ...]:
    return hensel_unit_rule_speeds(packet.speeds)


@lru_cache(maxsize=None)
def hensel_unit_counts_speeds(speeds: tuple[int, ...]) -> tuple[tuple[int, int, int, int], ...]:
    out = []
    for p in HENSEL_PS:
        roots, singular, zero_status = hensel_unit_signature_speeds(speeds, p)
        out.append((p, len(roots), len(singular), int(zero_status == "zero-singular")))
    return tuple(out)


def hensel_unit_counts(packet) -> tuple[tuple[int, int, int, int], ...]:
    return hensel_unit_counts_speeds(packet.speeds)


@dataclass(frozen=True)
class Split:
    name: str
    description: str
    key_func: object
    finite_state: int
    discrepancy: int
    unit_lift: int
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
        0,
        0,
        1,
    ),
    Split(
        "residue_terminal_fiber",
        "automatic word plus mod-14 MFC and terminal-state word",
        residue_terminal,
        5,
        1,
        0,
        0,
        0,
        1,
        2,
    ),
    Split(
        "residue_plus_et",
        "residue-terminal fiber plus Erdos-Turan clocks at 14,27,41",
        lambda p: (residue_terminal(p), et_clock(p)),
        4,
        5,
        0,
        0,
        0,
        1,
        3,
    ),
    Split(
        "residue_plus_unit_hensel",
        "residue-terminal fiber plus Henselian unit-root rule at 2,3,7",
        lambda p: (residue_terminal(p), hensel_unit_rule(p)),
        4,
        1,
        5,
        0,
        0,
        1,
        3,
    ),
    Split(
        "et_unit_convergence",
        "residue-terminal plus ET clocks and Henselian unit rule",
        lambda p: (residue_terminal(p), et_clock(p), hensel_unit_rule(p)),
        3,
        5,
        5,
        0,
        0,
        1,
        4,
    ),
    Split(
        "coarse_et_unit_gate",
        "same as ET+unit but with coarser ET and count-only unit data",
        lambda p: (residue_terminal(p), et_coarse_clock(p), hensel_unit_counts(p)),
        4,
        4,
        4,
        0,
        0,
        1,
        3,
    ),
    Split(
        "magnitude_cocycle",
        "adds exact M, q-threshold, Farey excess, and tail ratio",
        lambda p: (p.automatic_word, residue_mfc_pairs(p), magnitude_cocycle(p)),
        4,
        1,
        0,
        5,
        0,
        2,
        3,
    ),
    Split(
        "magnitude_et_unit",
        "magnitude cocycle plus ET clocks and Henselian unit rule",
        lambda p: (
            p.automatic_word,
            residue_mfc_pairs(p),
            magnitude_cocycle(p),
            et_clock(p),
            hensel_unit_rule(p),
        ),
        3,
        5,
        5,
        5,
        0,
        2,
        5,
    ),
    Split(
        "barcode_packet_zipper",
        "magnitude, barcode, and local labelled packet shadow",
        lambda p: (
            p.automatic_word,
            residue_mfc_pairs(p),
            magnitude_cocycle(p),
            barcode_shadow(p),
            local_packet_shadow(p),
        ),
        3,
        2,
        1,
        5,
        5,
        5,
        5,
    ),
)


@dataclass(frozen=True)
class SplitReport:
    split: Split
    fibers: int
    mixed_route_fibers: int
    mixed_family_fibers: int
    mixed_status_fibers: int
    max_fiber_size: int
    max_mixed_size: int
    largest_mixed_key: object
    largest_mixed_routes: tuple[tuple[str, int], ...]
    largest_mixed_families: tuple[tuple[str, int], ...]
    largest_mixed_status: tuple[tuple[str, int], ...]

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
            self.split.discrepancy,
            self.split.unit_lift,
            self.split.magnitude,
            self.split.topology,
            self.split.packet_label,
            self.split.finite_state,
            8 - self.split.proof_cost,
        )


def status(packet) -> str:
    return "boundary" if packet.strict_safe_mu == 0 else "open"


def split_groups(packets: list, split: Split) -> dict[object, list]:
    groups: dict[object, list] = defaultdict(list)
    for packet in packets:
        groups[split.key_func(packet)].append(packet)
    return groups


def split_report(packets: list, split: Split) -> SplitReport:
    groups = split_groups(packets, split)
    mixed_route = []
    mixed_family = []
    mixed_status = []
    for key, rows in groups.items():
        if len({p.route for p in rows}) > 1:
            mixed_route.append((key, rows))
        if len({p.family for p in rows}) > 1:
            mixed_family.append((key, rows))
        if len({status(p) for p in rows}) > 1:
            mixed_status.append((key, rows))

    if mixed_route:
        largest_key, largest_rows = max(mixed_route, key=lambda item: (len(item[1]), str(item[0])))
    else:
        largest_key, largest_rows = None, []

    return SplitReport(
        split=split,
        fibers=len(groups),
        mixed_route_fibers=len(mixed_route),
        mixed_family_fibers=len(mixed_family),
        mixed_status_fibers=len(mixed_status),
        max_fiber_size=max((len(rows) for rows in groups.values()), default=0),
        max_mixed_size=max((len(rows) for _, rows in mixed_route), default=0),
        largest_mixed_key=largest_key,
        largest_mixed_routes=tuple(Counter(p.route for p in largest_rows).most_common()),
        largest_mixed_families=tuple(Counter(p.family for p in largest_rows).most_common()),
        largest_mixed_status=tuple(Counter(status(p) for p in largest_rows).most_common()),
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
    scc_sizes = []
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


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge")
    print("  considered vertices:")
    print("    runners, gaps, fixed circle sections, section boundaries, residues,")
    print("    Fourier clocks, p-adic unit roots, zero-root scale debt, magnitude")
    print("    cocycles, barcode fibers, packet labels, and proof obligations.")
    print("  chosen vertices:")
    print("    quotient/proof-carrier bundles over HYP-2963 packets.")
    print("  preserved LRC predicate:")
    print("    theorem-route purity and boundary/open safety at threshold 1/14.")
    print("  destroyed information:")
    print("    raw runner identity, exact endpoint owners, and full Fejer atom banks")
    print("    until barcode/packet zipper fields are reattached.")
    print("  challenged assumption:")
    print("    Hensel information should be read through p-adic unit roots, not")
    print("    through the forced nonunit zero root alone.")
    print()


def print_reports(reports: list[SplitReport]) -> None:
    print("[1] Zipper convergence table")
    print(
        "  {name:25s} {fibers:>7s} {mixr:>7s} {mixf:>7s} {mixs:>7s} "
        "{maxf:>7s} {maxm:>7s} {purity:>7s}".format(
            name="split",
            fibers="fibers",
            mixr="mixR",
            mixf="mixF",
            mixs="mixS",
            maxf="maxF",
            maxm="maxMix",
            purity="purity",
        )
    )
    previous = None
    for report in reports:
        print(
            "  {name:25s} {fibers:7d} {mixr:7d} {mixf:7d} {mixs:7d} "
            "{maxf:7d} {maxm:7d} {purity:6.1f}%".format(
                name=report.split.name,
                fibers=report.fibers,
                mixr=report.mixed_route_fibers,
                mixf=report.mixed_family_fibers,
                mixs=report.mixed_status_fibers,
                maxf=report.max_fiber_size,
                maxm=report.max_mixed_size,
                purity=report.route_purity_score / 10,
            )
        )
        if previous:
            print(
                f"    delta_vs_previous: fibers={report.fibers - previous.fibers:+d} "
                f"mixR={report.mixed_route_fibers - previous.mixed_route_fibers:+d} "
                f"maxMix={report.max_mixed_size - previous.max_mixed_size:+d}"
            )
        previous = report
    print()
    print("  largest mixed fiber by split:")
    for report in reports:
        print(
            f"    {report.split.name}: max_mixed={report.max_mixed_size} "
            f"routes={dict(report.largest_mixed_routes)} status={dict(report.largest_mixed_status)}"
        )
    print()


def print_unit_rule_readout(packets: list) -> None:
    print("[2] Erdos-Turan and Henselian unit clocks")
    et_support = Counter(et_coarse_clock(p) for p in packets)
    unit_support = Counter(hensel_unit_counts(p) for p in packets)
    exact_unit_support = Counter(hensel_unit_rule(p) for p in packets)
    singular_unit = sum(
        1
        for p in packets
        if any(len(hensel_unit_signature(p, prime)[1]) for prime in HENSEL_PS)
    )
    zero_singular = sum(
        1
        for p in packets
        if any(hensel_unit_signature(p, prime)[2] == "zero-singular" for prime in HENSEL_PS)
    )
    print(f"  coarse_et_clock_fibers={len(et_support)} largest={et_support.most_common(1)[0][1]}")
    print(f"  hensel_unit_count_fibers={len(unit_support)} largest={unit_support.most_common(1)[0][1]}")
    print(f"  hensel_unit_exact_fibers={len(exact_unit_support)} largest={exact_unit_support.most_common(1)[0][1]}")
    print(f"  packets_with_singular_unit_root={singular_unit}")
    print(f"  packets_with_zero_singular_debt={zero_singular}")
    print("  readout:")
    print("    ET clocks are global phase/discrepancy clocks; they split many")
    print("    residue fibers but remain magnitude-blind.")
    print("    Unit-Hensel roots are local lift clocks; zero-singular status is")
    print("    separate scale debt and should not masquerade as a unit witness.")
    print()


def print_target_word(packets: list, reports: list[SplitReport]) -> None:
    rows = [p for p in packets if p.automatic_word == TARGET_WORD]
    print(f"[3] Target automatic word {TARGET_WORD}")
    print(f"  rows={len(rows)}")
    print(f"  routes={dict(Counter(p.route for p in rows))}")
    print(f"  families={dict(Counter(p.family for p in rows))}")
    print(f"  distinct_M={len({p.M for p in rows})}")
    print(f"  distinct_ET_clocks={len({et_clock(p) for p in rows})}")
    print(f"  distinct_unit_rules={len({hensel_unit_rule(p) for p in rows})}")
    print()
    for report in reports[1:8]:
        split = report.split
        groups: dict[object, list] = defaultdict(list)
        for packet in rows:
            groups[split.key_func(packet)].append(packet)
        mixed = [g for g in groups.values() if len({p.route for p in g}) > 1]
        print(f"  {split.name}: fibers={len(groups)} mixed_route={len(mixed)}")
        if mixed:
            largest = max(mixed, key=len)
            print(
                f"    largest_mixed_rows={len(largest)} "
                f"routes={dict(Counter(p.route for p in largest))}"
            )
            for packet in sample_rows(largest, 6):
                print(
                    "      {name:34s} M={M:>7s} mu={mu:>9s} "
                    "q0={q0:2d} unit={unit} route={route}".format(
                        name=packet.name[:34],
                        M=fmt(packet.M),
                        mu=fmt(packet.strict_safe_mu),
                        q0=packet.q_threshold,
                        unit=hensel_unit_counts(packet),
                        route=packet.route,
                    )
                )
    print()


def print_tournament(reports: list[SplitReport]) -> None:
    score_hist, c3, scc, hp, ranking = tournament_fingerprint(reports)
    print("[4] Tournament Analysis")
    print("  vertices_are=quotient/proof-carrier bundles, not runners")
    print("  observable=route purity, max mixed-fiber size, ET discrepancy,")
    print("             Henselian unit stability, magnitude retention,")
    print("             topology retention, packet-label retention, finite-state")
    print("             checkability, and proof cost")
    print("  switch=majority comparison of observable vectors")
    print("  tie_hamiltonian_path=" + " > ".join(s.name for s in SPLITS))
    print(f"  score_hist={dict(sorted(score_hist.items()))}")
    print(f"  directed_3cycles={c3}")
    print(f"  scc_sizes={scc}")
    print(f"  hamiltonian_path_count={hp}")
    print("  score_order=" + " > ".join(ranking))
    print()


def print_proof_readout() -> None:
    print("[5] Proof readout")
    print("  1. Erdos-Turan clocks improve the residue zipper, but they do not")
    print("     replace exact magnitude; they are best read as a discrepancy")
    print("     certificate scheduler.")
    print("  2. The Henselian unit rule separates genuine unit-root lifting from")
    print("     zero-root scale debt.  This makes p-adic data a routing rule,")
    print("     not a scalar discriminator.")
    print("  3. The magnitude cocycle is still the first tested non-route gate")
    print("     with zero mixed theorem-route fibers on the full bank.  ET+unit")
    print("     data explains which fibers should receive analytic or p-adic")
    print("     certificates before that exact cocycle is discharged.")
    print("  4. The next proof target is a convergence theorem: every residue")
    print("     fiber either has a bounded ET discrepancy certificate, a unit")
    print("     Hensel lift/zero-debt exit, a familywise magnitude formula,")
    print("     or named K33/F7/THM-572 residual debt.")
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

    print("=== LRC14 zipper-fiber convergence scout S188 ===")
    print(
        "bank=HYP-2963 default rows "
        f"single_limit={args.single_limit} two_swap_limit={args.two_swap_limit} "
        f"alias_depth={args.alias_depth} lcm_tail_max={args.lcm_tail_max}"
    )
    print(f"packets={len(packets)}")
    print()
    print_assumption_challenge()
    print_reports(reports)
    print_unit_rule_readout(packets)
    print_target_word(packets, reports)
    print_tournament(reports)
    print_proof_readout()


if __name__ == "__main__":
    main()
