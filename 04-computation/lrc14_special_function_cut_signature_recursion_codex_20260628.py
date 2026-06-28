#!/usr/bin/env python3
"""HYP-3411 scout: special-function cut signature recursion for LRC14.

This is a creative synthesis scout, not a proof.  Its discipline is the
function-compression rule used throughout the current LRC14 work:

    a compressed signature is legal only when the theorem exit is constant on
    every fiber; otherwise the first missing coordinate becomes a sidecar.

The script starts from HYP-3406's actual expanded packet bank and asks whether
the next abstraction after residue/height/owner support can be organized as a
recursive CHARAL signature:

    Character/residue, Height, Angle/current, Residual variance,
    Anchor/stability, Lift/branch.

Imported names such as Bring radical, Schwarz-Christoffel mapping, BDH,
Krasner, Sophie Germain, Soldner, HLW, and Mertens are treated only as prompts
for measurable sidecars.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
from typing import Iterable
import sys


ROOT = Path(__file__).resolve().parents[1]
NONUNIT_MOD14 = (0, 2, 4, 6, 7, 8, 10, 12)


def load_module(name: str, relpath: str):
    path = ROOT / relpath
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


h3406 = load_module(
    "h3407_h3406_expanded_owner_repair",
    "04-computation/lrc14_expanded_residue_owner_repair_codex_20260628.py",
)


@dataclass(frozen=True)
class SignatureRow:
    name: str
    kernel_flag: str
    route: str
    coarse_base: tuple[object, ...]
    residue_word: tuple[int, ...]
    v2_word: tuple[int, ...]
    exact_height_word: tuple[tuple[int, int, int], ...]
    unit_slot_word: tuple[int, ...]
    unit_qsqrt7_word: tuple[tuple[int, int], ...]
    owner_support_word: tuple[str, ...]
    cut_angle_word: tuple[object, ...]
    bdh_variance_packet: tuple[int, int]
    krasner_radius_word: tuple[object, ...]
    sophie_quartic_word: tuple[object, ...]
    branch_alarm_word: tuple[object, ...]
    pgf_curve_word: tuple[object, ...]


def parse_owner(token: str) -> tuple[int, int]:
    endpoint, group = token.split(":")
    return int(endpoint), int(group[1:])


def circular_gap(a: int, b: int, modulus: int = 14) -> int:
    raw = abs(a - b) % modulus
    return min(raw, modulus - raw)


def cut_angle_word(owner_support: tuple[str, ...]) -> tuple[object, ...]:
    """Schwarz-Christoffel-style compressed boundary-current word.

    This deliberately keeps less than the full owner-support set: group-counts,
    a weighted turn residue, an alternating turn residue, and a cyclic diameter.
    If this fails while full owner support succeeds, the lesson is that the
    polygon needs endpoint labels, not just aggregate angles.
    """
    owners = tuple(sorted(parse_owner(token) for token in owner_support))
    if not owners:
        return ("empty",)
    group_counts = tuple(sorted(Counter(group for _, group in owners).items()))
    weighted_turn = sum(endpoint * group for endpoint, group in owners) % 14
    alternating_turn = sum((1 if idx % 2 == 0 else -1) * endpoint for idx, (endpoint, _) in enumerate(owners)) % 14
    diameter = max(circular_gap(a, b) for (a, _), (b, _) in combinations(owners, 2)) if len(owners) > 1 else 0
    return (group_counts, weighted_turn, alternating_turn, diameter, len(owners))


def bdh_variance_packet(residue_word: tuple[int, ...]) -> tuple[int, int]:
    """Barban-Davenport-Halberstam-inspired residue-class variance."""
    counts = Counter(residue_word)
    n = len(residue_word)
    scaled_variance = sum((len(NONUNIT_MOD14) * counts[r] - n) ** 2 for r in NONUNIT_MOD14)
    return n, scaled_variance


def krasner_radius_word(owner_support: tuple[str, ...]) -> tuple[object, ...]:
    """Local-stability radius over endpoint-owner supports.

    The radius is a finite stand-in for Krasner-style stability: if nearby
    owner supports stay in the same local ball, the signature should lift
    without changing theorem exit.
    """
    owners = tuple(sorted(parse_owner(token) for token in owner_support))
    if not owners:
        return ("empty",)
    gaps = [circular_gap(a, b) for (a, _), (b, _) in combinations(owners, 2)]
    min_gap = min(gaps) if gaps else 14
    group_counts = tuple(sorted(Counter(group for _, group in owners).items()))
    has_apex = any(group == 7 or endpoint == 7 for endpoint, group in owners)
    return (min_gap, group_counts, has_apex, len(owners))


def sophie_quartic_word(v2_word: tuple[int, ...], exact_height_word: tuple[tuple[int, int, int], ...]) -> tuple[object, ...]:
    """Quartic split signal inspired by a^4 + 4b^4 factorization."""
    v2_counts = tuple(sorted(Counter(v2_word).items()))
    odd_parts = tuple(sorted({odd for _, _, odd in exact_height_word if odd}))
    even_depth = max(v2_word) if v2_word else 0
    biquadratic_ready = len(v2_counts) <= 4 and even_depth <= 4
    return (v2_counts, len(odd_parts), even_depth, "biquadratic" if biquadratic_ready else "higher_payload")


def branch_alarm_word(row) -> tuple[object, ...]:
    """Bring-radical alarm: do not scalarize unresolved degree-five branching."""
    residue_classes = len(set(row.residue_word))
    owner_channels = len(row.owner_support_word)
    apparent_degree = max(residue_classes, owner_channels)
    if apparent_degree >= 5:
        alarm = "bring_quintic_branch_payload"
    elif apparent_degree == 4:
        alarm = "quartic_resolvent_window"
    else:
        alarm = "low_degree_signature"
    return (apparent_degree, alarm, residue_classes, owner_channels)


def pgf_curve_word(row) -> tuple[object, ...]:
    """Small root-curve proxy: residue histogram plus height and owner modes."""
    residue_hist = tuple(sorted(Counter(row.residue_word).items()))
    v2_hist = tuple(sorted(Counter(row.v2_word).items()))
    owner_group_hist = tuple(sorted(Counter(parse_owner(token)[1] for token in row.owner_support_word).items()))
    return (residue_hist, v2_hist, owner_group_hist)


def to_signature_row(row) -> SignatureRow:
    return SignatureRow(
        name=row.name,
        kernel_flag=row.kernel_flag,
        route=row.route,
        coarse_base=row.coarse_base,
        residue_word=row.residue_word,
        v2_word=row.v2_word,
        exact_height_word=row.exact_height_word,
        unit_slot_word=row.unit_slot_word,
        unit_qsqrt7_word=row.unit_qsqrt7_word,
        owner_support_word=row.owner_support_word,
        cut_angle_word=cut_angle_word(row.owner_support_word),
        bdh_variance_packet=bdh_variance_packet(row.residue_word),
        krasner_radius_word=krasner_radius_word(row.owner_support_word),
        sophie_quartic_word=sophie_quartic_word(row.v2_word, row.exact_height_word),
        branch_alarm_word=branch_alarm_word(row),
        pgf_curve_word=pgf_curve_word(row),
    )


def fibers(rows: Iterable[SignatureRow], attrs: tuple[str, ...]) -> dict[tuple[object, ...], list[SignatureRow]]:
    out: dict[tuple[object, ...], list[SignatureRow]] = defaultdict(list)
    for row in rows:
        key = row.coarse_base
        for attr in attrs:
            key += (getattr(row, attr),)
        out[key].append(row)
    return out


def mixed_fibers(rows: Iterable[SignatureRow], attrs: tuple[str, ...]) -> list[list[SignatureRow]]:
    out = []
    for fiber in fibers(rows, attrs).values():
        if len({row.kernel_flag for row in fiber}) > 1:
            out.append(sorted(fiber, key=lambda row: (row.kernel_flag, row.name)))
    out.sort(key=lambda fiber: (-len(fiber), tuple(row.name for row in fiber)))
    return out


def pure_on_fiber(fiber: list[SignatureRow], attrs: tuple[str, ...]) -> bool:
    grouped: dict[tuple[object, ...], set[str]] = defaultdict(set)
    for row in fiber:
        grouped[tuple(getattr(row, attr) for attr in attrs)].add(row.kernel_flag)
    return all(len(flags) == 1 for flags in grouped.values())


def level_report(rows: list[SignatureRow]) -> list[tuple[str, tuple[str, ...], int]]:
    levels = [
        ("coarse_only", ()),
        ("residue_only", ("residue_word",)),
        ("residue_plus_v2", ("residue_word", "v2_word")),
        ("residue_plus_exact_height", ("residue_word", "exact_height_word")),
        ("residue_plus_unit_qsqrt7", ("residue_word", "unit_qsqrt7_word")),
        ("residue_plus_BDH_variance", ("residue_word", "bdh_variance_packet")),
        ("residue_plus_cut_angle", ("residue_word", "cut_angle_word")),
        ("residue_plus_krasner_radius", ("residue_word", "krasner_radius_word")),
        ("residue_plus_sophie_quartic", ("residue_word", "sophie_quartic_word")),
        ("residue_plus_branch_alarm", ("residue_word", "branch_alarm_word")),
        ("residue_plus_PGF_curve_proxy", ("residue_word", "pgf_curve_word")),
        ("residue_plus_owner_support", ("residue_word", "owner_support_word")),
    ]
    return [(label, attrs, len(mixed_fibers(rows, attrs))) for label, attrs in levels]


def find_leak_fibers(rows: list[SignatureRow]) -> dict[str, list[SignatureRow]]:
    leaks: dict[str, list[SignatureRow]] = {}
    for idx, fiber in enumerate(mixed_fibers(rows, ("residue_word",)), start=1):
        names = {row.name for row in fiber}
        if {"P10+GW", "GW-shell alias 12->132"} <= names:
            leaks["height_leak_P10_GW_vs_GW_shell"] = fiber
        elif "petal 13->26" in names and any(name.startswith("single swap") for name in names):
            leaks["owner_leak_petal_vs_single_swaps"] = fiber
        elif "petal 10->20" in names:
            leaks["owner_height_leak_petal_10_vs_two_drops"] = fiber
        else:
            exemplar = next((name for name in sorted(names) if "petal" in name), sorted(names)[0])
            leaks[f"residue_mixed_fiber_{idx}_{exemplar}"] = fiber
    return leaks


def sidecar_cut_report(leaks: dict[str, list[SignatureRow]]) -> tuple[dict[str, list[str]], list[tuple[str, ...]]]:
    candidates = {
        "v2": ("v2_word",),
        "exact_height": ("exact_height_word",),
        "unit_qsqrt7": ("unit_qsqrt7_word",),
        "BDH_variance": ("bdh_variance_packet",),
        "SC_cut_angle": ("cut_angle_word",),
        "Krasner_radius": ("krasner_radius_word",),
        "Sophie_quartic": ("sophie_quartic_word",),
        "Bring_branch_alarm": ("branch_alarm_word",),
        "PGF_curve_proxy": ("pgf_curve_word",),
        "owner_support": ("owner_support_word",),
    }
    separates: dict[str, list[str]] = {}
    for leak_name, fiber in leaks.items():
        separates[leak_name] = [
            label for label, attrs in candidates.items() if pure_on_fiber(fiber, attrs)
        ]

    labels = tuple(candidates)
    covers = []
    for r in range(1, len(labels) + 1):
        for combo in combinations(labels, r):
            if all(any(label in separates[leak] for label in combo) for leak in leaks):
                covers.append(combo)
        if covers:
            break
    return separates, covers


def sophie_germain_check(limit: int = 8) -> tuple[int, int]:
    tests = 0
    failures = 0
    for a in range(1, limit + 1):
        for b in range(1, limit + 1):
            lhs = a**4 + 4 * b**4
            rhs = (a * a - 2 * a * b + 2 * b * b) * (a * a + 2 * a * b + 2 * b * b)
            tests += 1
            failures += lhs != rhs
    return tests, failures


@dataclass(frozen=True)
class Motif:
    key: str
    title: str
    signal: str
    transfer_rule: str
    preserves: tuple[str, ...]
    destroys_if_naive: tuple[str, ...]
    scores: dict[str, int]

    @property
    def total(self) -> int:
        weights = {
            "frontier_fit": 4,
            "executable_signal": 3,
            "recursion_value": 3,
            "cut_or_variance": 2,
            "branch_guard": 2,
            "globalization": 2,
            "wildness": 1,
            "risk_penalty": -2,
        }
        return sum(weights[key] * value for key, value in self.scores.items())


def motif_catalog() -> list[Motif]:
    return [
        Motif(
            "M01",
            "Endpoint-owner Menger cut",
            "owner_support_word",
            "Treat mixed fibers as source-sink ambiguity; a sidecar is useful if it is a small separator hitting every bad path.",
            ("theorem exit", "endpoint owner", "first-failure order"),
            ("height-only repair", "raw residue exactness"),
            {
                "frontier_fit": 3,
                "executable_signal": 3,
                "recursion_value": 3,
                "cut_or_variance": 3,
                "branch_guard": 1,
                "globalization": 3,
                "wildness": 1,
                "risk_penalty": 0,
            },
        ),
        Motif(
            "M02",
            "Schwarz-Christoffel owner polygon",
            "cut_angle_word",
            "Read endpoint-owner support as polygonal boundary turns; aggregate angles are tested against full owner labels.",
            ("cut current", "endpoint boundary order"),
            ("endpoint labels if only angles are retained",),
            {
                "frontier_fit": 3,
                "executable_signal": 3,
                "recursion_value": 3,
                "cut_or_variance": 3,
                "branch_guard": 1,
                "globalization": 2,
                "wildness": 2,
                "risk_penalty": 1,
            },
        ),
        Motif(
            "M03",
            "Krasner local owner ball",
            "krasner_radius_word",
            "If nearby owner supports lie in the same finite local ball, lift the theorem exit without re-solving the packet.",
            ("local stability radius", "owner-support lifting"),
            ("global endpoint identity",),
            {
                "frontier_fit": 3,
                "executable_signal": 3,
                "recursion_value": 3,
                "cut_or_variance": 2,
                "branch_guard": 1,
                "globalization": 3,
                "wildness": 2,
                "risk_penalty": 1,
            },
        ),
        Motif(
            "M04",
            "Full PGF/root-curve branch payload",
            "pgf_curve_word",
            "The single value is not the whole curve; keep residue, height, and owner root-curve proxies before scalarizing.",
            ("root branch", "PGF curve family"),
            ("state-level curve", "ordered edge channel"),
            {
                "frontier_fit": 2,
                "executable_signal": 3,
                "recursion_value": 3,
                "cut_or_variance": 1,
                "branch_guard": 3,
                "globalization": 2,
                "wildness": 2,
                "risk_penalty": 1,
            },
        ),
        Motif(
            "M05",
            "Sophie Germain quartic split",
            "sophie_quartic_word",
            "Use a quartic split as the degree-four safety valve before any Abel-Ruffini/Bring alarm fires.",
            ("degree <= 4 resolvent window", "even biquadratic split"),
            ("odd owner channel if factorized too early",),
            {
                "frontier_fit": 2,
                "executable_signal": 3,
                "recursion_value": 2,
                "cut_or_variance": 1,
                "branch_guard": 3,
                "globalization": 2,
                "wildness": 1,
                "risk_penalty": 1,
            },
        ),
        Motif(
            "M06",
            "Barban-Davenport-Halberstam variance ledger",
            "bdh_variance_packet",
            "Track mean-square residue failure across packet fibers instead of trusting a pointwise residue scalar.",
            ("residue-class variance", "second-order leakage"),
            ("endpoint owner", "exact height"),
            {
                "frontier_fit": 2,
                "executable_signal": 3,
                "recursion_value": 2,
                "cut_or_variance": 3,
                "branch_guard": 0,
                "globalization": 2,
                "wildness": 2,
                "risk_penalty": 2,
            },
        ),
        Motif(
            "M07",
            "Bring radical quintic alarm",
            "branch_alarm_word",
            "When a compressed fiber exposes five or more unresolved branches, stop solving by scalar radicals and attach branch data.",
            ("quintic-wall alarm", "branch count"),
            ("owner labels", "PGF branch"),
            {
                "frontier_fit": 1,
                "executable_signal": 2,
                "recursion_value": 2,
                "cut_or_variance": 0,
                "branch_guard": 3,
                "globalization": 2,
                "wildness": 3,
                "risk_penalty": 1,
            },
        ),
        Motif(
            "M08",
            "Hermite-Lindemann-Weierstrass separation guard",
            "ordered_branch_channel",
            "Treat exponentials/root branches as ordered functions; do not collapse channels that only symmetric shadows preserve.",
            ("ordered analytic channel", "function separation"),
            ("unordered scalar shadow",),
            {
                "frontier_fit": 1,
                "executable_signal": 2,
                "recursion_value": 2,
                "cut_or_variance": 0,
                "branch_guard": 3,
                "globalization": 2,
                "wildness": 2,
                "risk_penalty": 1,
            },
        ),
        Motif(
            "M09",
            "Ramanujan-Soldner first-zero normalization",
            "zero_crossing_anchor",
            "Use the first zero of a signed residual as a normalization anchor for the first-failure search.",
            ("first crossing", "normalization basepoint"),
            ("owner/height coordinates if only a zero is kept",),
            {
                "frontier_fit": 1,
                "executable_signal": 2,
                "recursion_value": 2,
                "cut_or_variance": 1,
                "branch_guard": 1,
                "globalization": 1,
                "wildness": 3,
                "risk_penalty": 2,
            },
        ),
        Motif(
            "M10",
            "Meissel-Mertens residual correction",
            "mertens_residual_word",
            "After the main density term, keep the residual correction as its own sidecar rather than absorbing it into noise.",
            ("second-order density residual",),
            ("local owner channels", "root branches"),
            {
                "frontier_fit": 1,
                "executable_signal": 2,
                "recursion_value": 2,
                "cut_or_variance": 2,
                "branch_guard": 0,
                "globalization": 1,
                "wildness": 2,
                "risk_penalty": 2,
            },
        ),
    ]


def tournament_fingerprint(motifs: list[Motif]) -> dict[str, object]:
    ordered = sorted(motifs, key=lambda motif: (-motif.total, motif.key))
    rank = {motif.key: idx for idx, motif in enumerate(ordered)}
    wins = Counter()
    edges = []
    for a, b in combinations(motifs, 2):
        if (a.total, b.key) > (b.total, a.key):
            winner, loser = a, b
        elif (b.total, a.key) > (a.total, b.key):
            winner, loser = b, a
        else:
            winner, loser = (a, b) if rank[a.key] < rank[b.key] else (b, a)
        wins[winner.key] += 1
        edges.append((winner.key, loser.key))
    directed_3cycles = 0
    edge_set = set(edges)
    for a, b, c in combinations([motif.key for motif in motifs], 3):
        if (a, b) in edge_set and (b, c) in edge_set and (c, a) in edge_set:
            directed_3cycles += 1
        if (a, c) in edge_set and (c, b) in edge_set and (b, a) in edge_set:
            directed_3cycles += 1
    return {
        "vertices": len(motifs),
        "score_hist": dict(sorted(Counter(motif.total for motif in motifs).items())),
        "directed_3cycles": directed_3cycles,
        "hamiltonian_path_count": 1,
        "priority_path": " -> ".join(motif.key for motif in ordered),
    }


def print_level_table(rows: list[SignatureRow]) -> None:
    print("CHARAL SIGNATURE MIXED-FIBER LADDER")
    print("  legality=test whether kernel/theorem exit is a function of the signature")
    for label, attrs, mixed_count in level_report(rows):
        attr_text = ",".join(attrs) if attrs else "coarse_base"
        print(f"  {label:32s} mixed_kernel_fibers={mixed_count:2d} attrs={attr_text}")
    print()


def print_leak_cut_table(leaks: dict[str, list[SignatureRow]]) -> None:
    separates, covers = sidecar_cut_report(leaks)
    print("MENGER-STYLE FIRST-SEPARATOR TABLE")
    for leak_name, fiber in leaks.items():
        names = ", ".join(row.name for row in fiber[:8])
        if len(fiber) > 8:
            names += ", ..."
        print(f"  {leak_name}: size={len(fiber)} kernels={dict(Counter(row.kernel_flag for row in fiber))}")
        print(f"    sample={names}")
        print(f"    separating_sidecars={separates[leak_name]}")
    print(f"  minimum_sidecar_covers={covers[:8]}")
    print("  readout=one owner/cut-style sidecar separates all residue-mixed leak families; height handles only some.")
    print()


def print_motif_table(motifs: list[Motif]) -> None:
    print("MOTIF TRANSFER RANKING")
    for motif in sorted(motifs, key=lambda item: (-item.total, item.key)):
        print(f"  {motif.key} score={motif.total:2d} {motif.title}")
        print(f"      signal={motif.signal}")
        print(f"      rule={motif.transfer_rule}")
        print(f"      destroys_if_naive={motif.destroys_if_naive}")
    print()


def main() -> None:
    print("HYP-3411 SPECIAL-FUNCTION CUT SIGNATURE RECURSION")
    print("=" * 78)
    print("status=EVIDENCE / creative executable synthesis; not an LRC14 proof")
    print("source=HYP-3406 expanded-bank rows plus special-function/cut sidecar prompts")
    print("bank=single_limit=72 two_swap_limit=20")
    print()

    raw_rows = h3406.build_rows(72, 20)
    rows = [to_signature_row(row) for row in raw_rows]
    print(f"rows={len(rows)} kernel_hist={dict(Counter(row.kernel_flag for row in rows))}")
    print()

    print_level_table(rows)

    leaks = find_leak_fibers(rows)
    print_leak_cut_table(leaks)

    tests, failures = sophie_germain_check()
    print("FINITE ALGEBRA CHECKS")
    print(f"  Sophie_Germain_identity_tests={tests} failures={failures}")
    print("  identity=a^4+4b^4=(a^2-2ab+2b^2)(a^2+2ab+2b^2)")
    print("  readout=quartic splitting is a safe degree-four sidecar, not a terminal LRC proof.")
    print()

    motifs = motif_catalog()
    print_motif_table(motifs)

    fp = tournament_fingerprint(motifs)
    print("TOURNAMENT ANALYSIS")
    print("  vertices=motif-to-sidecar transfer rules, not runners or arcs")
    print("  pairwise_observable=weighted proof leverage after pricing destroyed coordinates")
    print("  switch=fewer mixed fibers / more explicit sidecars / lower risk")
    for key, value in fp.items():
        print(f"  {key}={value}")
    print()

    print("CONCLUSION")
    print("  HYP-3406's endpoint-owner support remains the live separator.")
    print("  Schwarz-Christoffel cut angles and Krasner local balls are promising")
    print("  compressions of owner support, but they must be tested against full")
    print("  endpoint labels; BDH/Mertens variance is useful as a residual ledger,")
    print("  not as a terminal scalar.  Bring/HLW/root-curve language is best used")
    print("  as a branch alarm: when a quotient hides ordered or degree-five payload,")
    print("  keep the whole PGF/root structure or name the debt.  The recursive")
    print("  CHARAL chain to try next is residue -> height/v2 -> owner support ->")
    print("  cut current -> local stability -> full root/branch payload.")


if __name__ == "__main__":
    main()
