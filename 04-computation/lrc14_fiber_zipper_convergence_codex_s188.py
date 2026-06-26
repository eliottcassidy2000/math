#!/usr/bin/env python3
"""Fiber-zipper convergence scout for LRC14.

This extends HYP-3023 by asking whether the mixed automatic/residue fibers
converge before the exact magnitude cocycle is attached.  The tested middle
teeth are:

* Erdos-Turan low-frequency residue discrepancy bins;
* a Henselian unit rule for A_S(x)=sum_{v in S} x^v on p-adic unit roots;
* a coarse q-threshold / unit-excess rule.

Tournament Analysis declaration:
  vertices: zipper teeth / proof carriers, not runners;
  pairwise observable: route purity, max mixed-fiber size, convergence from the
    previous tooth, retained analytic discrepancy, retained Henselian unit
    data, retained magnitude, packet-label retention, and proof cost;
  switch/gauge: majority comparison of the observable vector;
  tie Hamiltonian path: the listed zipper order.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
from math import cos, gcd, pi, sin, sqrt
from pathlib import Path
import argparse
import os
import sys


REPO = Path(__file__).resolve().parents[1]
TARGET_WORD = "MFCMMCCFFFCCC"
LOCAL_PRIMES = (2, 3, 7)
ET_MODULI = (14, 27, 41)
THRESHOLD_DEN = 14


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


def residue_l1_bucket(speeds: tuple[int, ...], modulus: int) -> int:
    counts = Counter(v % modulus for v in speeds)
    n = len(speeds)
    l1 = sum(abs(Fraction(counts[r], 1) - Fraction(n, modulus)) for r in range(modulus))
    return int(l1 * 28)


def erdos_turan_proxy(speeds: tuple[int, ...], modulus: int, hmax: int = 6) -> int:
    total = 0.0
    for h in range(1, hmax + 1):
        real = 0.0
        imag = 0.0
        for v in speeds:
            angle = 2.0 * pi * h * (v % modulus) / modulus
            real += cos(angle)
            imag += sin(angle)
        total += sqrt(real * real + imag * imag) / h
    return int(round(1000.0 * total))


def erdos_turan_signature(packet) -> tuple[tuple[int, int], ...]:
    """Coarse low-frequency discrepancy signature.

    The first coordinate is an exact L1 residue-discrepancy bucket.  The second
    is a low-frequency Erdos-Turan proxy bucket.  These are intentionally
    coarser than exact magnitude: they model convergence pressure, not identity.
    """

    return tuple(
        (
            residue_l1_bucket(packet.speeds, modulus),
            erdos_turan_proxy(packet.speeds, modulus) // 250,
        )
        for modulus in ET_MODULI
    )


def vp(n: int, p: int) -> int:
    if n == 0:
        return 99
    n = abs(n)
    out = 0
    while n % p == 0:
        n //= p
        out += 1
    return out


def unit_root_counts(speeds: tuple[int, ...], p: int) -> tuple[int, int, tuple[int, ...]]:
    """Return simple/singular unit roots plus exponent counts on F_p^*.

    A unit root x with A_S(x)=0 and A'_S(x) nonzero is the Hensel-stable case:
    it lifts uniquely.  A singular unit root is retained as local lift debt.
    """

    simple = 0
    singular = 0
    period = max(1, p - 1)
    exp_counts = [0] * period
    for v in speeds:
        exp_counts[v % period] += 1
    for x in range(1, p):
        val = sum(pow(x, v, p) for v in speeds) % p
        der = sum((v % p) * pow(x, v - 1, p) for v in speeds) % p
        if val == 0:
            if der == 0:
                singular += 1
            else:
                simple += 1
    return simple, singular, tuple(exp_counts)


def local_denominator_unit(packet, p: int) -> tuple[int, int, int, int]:
    den = packet.M.denominator
    num = packet.M.numerator
    pe2 = p * p
    den_v = vp(den, p)
    num_v = vp(num, p)
    den_unit = 0 if den_v else den % pe2
    excess_unit = packet.farey_excess % pe2
    return den_v, num_v, den_unit, excess_unit


def henselian_unit_rule(packet) -> tuple[tuple[int, int, tuple[int, ...], tuple[int, int, int, int]], ...]:
    return tuple(
        (*unit_root_counts(packet.speeds, p), local_denominator_unit(packet, p))
        for p in LOCAL_PRIMES
    )


def q_unit_excess_rule(packet) -> tuple[str, str, int, int]:
    if packet.q_threshold <= 13:
        qclass = "q<=13"
    elif packet.q_threshold == 14:
        qclass = "q=14"
    else:
        qclass = "q>14"
    excess = 14 * packet.M.numerator - packet.M.denominator
    if excess == 0:
        exclass = "boundary"
    elif excess > 0:
        exclass = "open-excess"
    else:
        exclass = "below"
    return qclass, exclass, packet.farey_excess % 14, packet.M.denominator.bit_length()


@dataclass(frozen=True)
class Split:
    name: str
    description: str
    key_func: object
    finite_state: int
    erdos_turan: int
    hensel_unit: int
    q_rule: int
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
        0,
        1,
    ),
    Split(
        "residue_terminal_fiber",
        "automatic word plus residue and DFA terminal states",
        lambda p: (p.automatic_word, residue_mfc_pairs(p), terminal_word(p)),
        5,
        0,
        0,
        0,
        1,
        0,
        1,
        2,
    ),
    Split(
        "erdos_turan_residue_zipper",
        "adds low-frequency residue discrepancy bins",
        lambda p: (p.automatic_word, residue_mfc_pairs(p), terminal_word(p), erdos_turan_signature(p)),
        5,
        5,
        0,
        0,
        1,
        0,
        1,
        3,
    ),
    Split(
        "henselian_unit_zipper",
        "adds p-adic unit-root lift rule for p=2,3,7",
        lambda p: (p.automatic_word, residue_mfc_pairs(p), terminal_word(p), henselian_unit_rule(p)),
        5,
        0,
        5,
        1,
        1,
        0,
        1,
        3,
    ),
    Split(
        "et_hensel_unit_zipper",
        "zips Erdos-Turan residue bins with the Henselian unit rule",
        lambda p: (
            p.automatic_word,
            residue_mfc_pairs(p),
            terminal_word(p),
            erdos_turan_signature(p),
            henselian_unit_rule(p),
        ),
        5,
        5,
        5,
        1,
        1,
        0,
        1,
        4,
    ),
    Split(
        "et_hensel_qrule_zipper",
        "adds q-threshold and coarse unit-excess lane without exact M",
        lambda p: (
            p.automatic_word,
            residue_mfc_pairs(p),
            terminal_word(p),
            erdos_turan_signature(p),
            henselian_unit_rule(p),
            q_unit_excess_rule(p),
        ),
        5,
        5,
        5,
        5,
        2,
        0,
        2,
        5,
    ),
    Split(
        "magnitude_cocycle",
        "adds exact M, q-threshold, Farey excess, and tail ratio",
        lambda p: (p.automatic_word, residue_mfc_pairs(p), magnitude_cocycle(p)),
        4,
        1,
        1,
        5,
        5,
        0,
        2,
        6,
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
        1,
        1,
        5,
        5,
        5,
        2,
        7,
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
        1,
        1,
        5,
        5,
        5,
        5,
        8,
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
    mixed_rows: int
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
            -self.mixed_rows,
            self.split.erdos_turan,
            self.split.hensel_unit,
            self.split.q_rule,
            self.split.magnitude,
            self.split.topology,
            self.split.packet_label,
            self.split.finite_state,
            10 - self.split.proof_cost,
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
        mixed_rows=sum(len(rows) for _, rows in mixed_route),
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


def target_split_report(packets: list, split: Split):
    rows = [p for p in packets if p.automatic_word == TARGET_WORD]
    groups: dict[object, list] = defaultdict(list)
    for packet in rows:
        groups[split.key_func(packet)].append(packet)
    mixed = [packet_rows for packet_rows in groups.values() if len({p.route for p in packet_rows}) > 1]
    largest = max(mixed, key=len) if mixed else []
    return rows, len(groups), mixed, largest


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge")
    print("  considered vertices:")
    print("    runners, gaps, fixed circle sections, residue classes, Erdos-Turan")
    print("    Fourier modes, Henselian unit roots, denominator units, barcode bars,")
    print("    packet labels, and proof obligations.")
    print("  chosen vertices:")
    print("    zipper teeth / proof-carrier bundles over HYP-2963 packets.")
    print("  preserved LRC predicate:")
    print("    theorem-route purity and boundary-versus-open status at threshold 1/14.")
    print("  challenged assumption:")
    print("    exact magnitude is not the only possible convergence carrier; test")
    print("    analytic discrepancy and p-adic unit-lift rules before accepting that")
    print("    route purity needs near-identity packet data.")
    print()


def print_reports(reports: list[SplitReport]) -> None:
    print("[1] Convergence ladder")
    print(
        "  {name:29s} {fibers:>7s} {mixed:>7s} {mfam:>7s} {maxf:>7s} "
        "{maxm:>7s} {mrows:>8s} {purity:>7s}"
        .format(
            name="split",
            fibers="fibers",
            mixed="mixR",
            mfam="mixF",
            maxf="maxF",
            maxm="maxMix",
            mrows="mixRows",
            purity="purity",
        )
    )
    prior: SplitReport | None = None
    for report in reports:
        print(
            "  {name:29s} {fibers:7d} {mixed:7d} {mfam:7d} {maxf:7d} "
            "{maxm:7d} {mrows:8d} {purity:6.1f}%".format(
                name=report.split.name,
                fibers=report.fibers,
                mixed=report.mixed_route_fibers,
                mfam=report.mixed_family_fibers,
                maxf=report.max_fiber_size,
                maxm=report.max_mixed_size,
                mrows=report.mixed_rows,
                purity=report.route_purity_score / 10,
            )
        )
        if prior:
            print(
                "    delta from previous: mixed_fibers {:+d}, max_mixed {:+d}, mixed_rows {:+d}".format(
                    report.mixed_route_fibers - prior.mixed_route_fibers,
                    report.max_mixed_size - prior.max_mixed_size,
                    report.mixed_rows - prior.mixed_rows,
                )
            )
        prior = report
    print()
    print("  largest mixed fiber by split:")
    for report in reports:
        print(
            f"    {report.split.name}: max_mixed={report.max_mixed_size} "
            f"routes={dict(report.largest_mixed_routes)}"
        )
    print()


def print_target_word_templates(packets: list) -> None:
    rows = [p for p in packets if p.automatic_word == TARGET_WORD]
    print(f"[3] Target automatic word {TARGET_WORD}")
    print(f"  rows={len(rows)}")
    print(f"  routes={dict(Counter(p.route for p in rows))}")
    print(f"  families={dict(Counter(p.family for p in rows))}")
    print(f"  distinct M values={len({p.M for p in rows})}")
    print()
    for split in SPLITS[1:]:
        _, fiber_count, mixed, largest = target_split_report(packets, split)
        print(
            f"  {split.name}: fibers={fiber_count} mixed_route={len(mixed)} "
            f"largest_mixed={len(largest)}"
        )
        if largest:
            print(f"    routes={dict(Counter(p.route for p in largest))}")
            for packet in sample_rows(largest, 6):
                print(
                    "      {name:34s} M={M:>7s} mu={mu:>9s} q0={q0:2d} "
                    "et={et} hU={hu} qrule={qr} route={route}".format(
                        name=packet.name[:34],
                        M=fmt(packet.M),
                        mu=fmt(packet.strict_safe_mu),
                        q0=packet.q_threshold,
                        et=erdos_turan_signature(packet),
                        hu=tuple((x[0], x[1], x[3]) for x in henselian_unit_rule(packet)),
                        qr=q_unit_excess_rule(packet),
                        route=packet.route,
                    )
                )
    print()


def print_rule_readout(reports: list[SplitReport]) -> None:
    by_name = {report.split.name: report for report in reports}
    residue = by_name["residue_terminal_fiber"]
    et = by_name["erdos_turan_residue_zipper"]
    hensel = by_name["henselian_unit_zipper"]
    both = by_name["et_hensel_unit_zipper"]
    qrule = by_name["et_hensel_qrule_zipper"]
    mag = by_name["magnitude_cocycle"]
    print("[2] Proof readout")
    print("  Residue-terminal fibers are close but still mixed.")
    print(
        "  Erdos-Turan bins alone change mixed fibers by {:+d} and max mixed by {:+d}.".format(
            et.mixed_route_fibers - residue.mixed_route_fibers,
            et.max_mixed_size - residue.max_mixed_size,
        )
    )
    print(
        "  Henselian unit-rule bins alone change mixed fibers by {:+d} and max mixed by {:+d}.".format(
            hensel.mixed_route_fibers - residue.mixed_route_fibers,
            hensel.max_mixed_size - residue.max_mixed_size,
        )
    )
    print(
        "  Zipping both changes mixed fibers by {:+d} and max mixed by {:+d}.".format(
            both.mixed_route_fibers - residue.mixed_route_fibers,
            both.max_mixed_size - residue.max_mixed_size,
        )
    )
    print(
        "  Adding the coarse q/unit-excess lane leaves {mixed} mixed route fibers; exact magnitude leaves {mag}.".format(
            mixed=qrule.mixed_route_fibers,
            mag=mag.mixed_route_fibers,
        )
    )
    if hensel.mixed_route_fibers == 0:
        print("  Interpretation: on this target fiber, the Henselian unit rule is")
        print("  already a route-pure convergence carrier, while Erdos-Turan bins explain")
        print("  most of the visible contraction before the p-adic unit rule finishes it.")
        print("  The next lemma should prove this unit-lift split familywise, then stress")
        print("  the same rule on the full HYP-2963 bank before replacing exact magnitude.")
    else:
        print("  Interpretation: Erdos-Turan controls convergence pressure and the")
        print("  Henselian unit rule names local lift debt, but the current coarse bins")
        print("  do not yet replace the exact magnitude cocycle.  The next lemma should")
        print("  prove which q/unit-excess sublanes admit compression without route mixing.")
    print()


def print_tournament(reports: list[SplitReport]) -> None:
    score_hist, c3, scc, hp, ranking = tournament_fingerprint(reports)
    print("[4] Tournament Analysis")
    print("  vertices_are=zipper teeth / proof carriers, not runners")
    print("  observable=route purity, max mixed-fiber size, convergence from prior")
    print("             tooth, Erdos-Turan retention, Henselian-unit retention,")
    print("             magnitude retention, packet-label retention, proof cost")
    print("  switch=majority comparison of observable vectors")
    print("  tie_hamiltonian_path=" + " > ".join(s.name for s in SPLITS))
    print(f"  score_hist={dict(sorted(score_hist.items()))}")
    print(f"  directed_3cycles={c3}")
    print(f"  scc_sizes={scc}")
    print(f"  hamiltonian_path_count={hp}")
    print("  score_order=" + " > ".join(ranking))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--single-limit", type=int, default=180)
    parser.add_argument("--two-swap-limit", type=int, default=36)
    parser.add_argument("--alias-depth", type=int, default=4)
    parser.add_argument("--lcm-tail-max", type=int, default=5)
    parser.add_argument("--workers", type=int, default=max(1, min(os.cpu_count() or 1, 8)))
    parser.add_argument("--target-word", default=TARGET_WORD)
    parser.add_argument(
        "--full-bank",
        action="store_true",
        help="Audit the whole HYP-2963 bank. Default filters to --target-word before exact packet computation.",
    )
    args = parser.parse_args()

    rows = lp.build_bank(args.single_limit, args.two_swap_limit, args.alias_depth, args.lcm_tail_max)
    all_rows = len(rows)
    if not args.full_bank:
        rows = [row for row in rows if lp.automatic_word(row[2]) == args.target_word]
    packets = lp.compute_packets(rows, args.workers)
    reports = [split_report(packets, split) for split in SPLITS]

    print("=== LRC14 fiber-zipper convergence S188 ===")
    print(
        "bank=HYP-2963 "
        + ("full default rows " if args.full_bank else f"target automatic word {args.target_word} ")
        + f"single_limit={args.single_limit} two_swap_limit={args.two_swap_limit} "
        f"alias_depth={args.alias_depth} lcm_tail_max={args.lcm_tail_max}"
    )
    print(f"candidate_rows={len(rows)} of {all_rows}")
    print(f"packets={len(packets)}")
    print()
    print_assumption_challenge()
    print_reports(reports)
    print_rule_readout(reports)
    print_target_word_templates(packets)
    print_tournament(reports)


if __name__ == "__main__":
    main()
