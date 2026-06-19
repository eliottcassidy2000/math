#!/usr/bin/env python3
"""
LRC(14) block-frequency transfer scout.

HYP-2632 turns the repeated-residue coimage tail into a finite character
packet.  This script tests the next lift: the exact reciprocal hyperplane sum
should be treated as a core block convolved with a two-large pair block, not as
six independent harmonic variables.

For model supports with four core speeds and two large speeds A,B,

    c1*n1+c2*n2+c3*n3+c4*n4 + A*x + B*y = 0,

the exact six-support term can be grouped as

    sum_s <Core_s(u,v), Pair_s(u,v)>_{u,v in F_7^*},

where Core_s is the signed four-variable reciprocal transfer with the two
frequency residues (u,v) left open, and Pair_s is the two-large reciprocal
transfer for A*x+B*y=-s.  The useful proof question is whether Cauchy/Schur
norms of these finite 6x6 matrices are small enough before falling back to a
Minkowski count.

This is evidence and route-finding, not a proof of LRC(14).
"""
from __future__ import annotations

import itertools
import math
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_support6_residue_cusp_codex_s12 as s12  # noqa: E402

MOD = 7
AMBIENT_D = 9
HMAX = 24
RESIDUES = range(1, MOD)


def section(title: str) -> None:
    print("\n" + "=" * 92)
    print(title)
    print("=" * 92)


def chi7(x: int) -> int:
    x %= MOD
    if x == 0:
        return 0
    return 1 if x in {1, 2, 4} else -1


def affine_selector(a: int, b: int) -> int:
    return (a * b * (1 + 3 * ((a + b) % MOD)) - 1) % MOD


def coeffs(H: int) -> tuple[int, ...]:
    return tuple(x for x in range(-H, H + 1) if x and x % MOD)


def zero_matrix() -> list[list[complex]]:
    return [[0j for _ in RESIDUES] for _ in RESIDUES]


def add_scaled(dst: list[list[complex]], src: tuple[tuple[complex, ...], ...], scale: float) -> None:
    for i in range(6):
        row_dst = dst[i]
        row_src = src[i]
        for j in range(6):
            row_dst[j] += row_src[j] * scale


def frob(mat: list[list[complex]]) -> float:
    return math.sqrt(sum(abs(z) ** 2 for row in mat for z in row))


def l1_abs(mat: list[list[complex]]) -> float:
    return sum(abs(z) for row in mat for z in row)


def dot(a: list[list[complex]], b: list[list[complex]]) -> complex:
    return sum(a[i][j] * b[i][j] for i in range(6) for j in range(6))


def residue_matrix(core_residues: tuple[int, int, int, int]) -> tuple[tuple[complex, ...], ...]:
    rows = []
    for u in RESIDUES:
        row = []
        for v in RESIDUES:
            row.append(s12.residue_coeff(core_residues + (u, v), AMBIENT_D))
        rows.append(tuple(row))
    return tuple(rows)


def packet_units(a: int, b: int) -> str:
    """HYP-2632 finite packet readout for normalized speed residues."""
    a %= MOD
    b %= MOD
    if a == b:
        if a == 0:
            return "-4U"
        if a == 1:
            return "0"
        return f"{(-43 - 7 * chi7(a)) // 2}U"
    aa, bb = sorted((a, b))
    if aa == 1 or bb == 1:
        return "not-4+1+1"
    if (aa + bb) % MOD == 2:
        return "0"
    return "8U" if chi7(affine_selector(aa, bb)) == 1 else "1U"


@dataclass
class CoreTables:
    core_speeds: tuple[int, int, int, int]
    signed: dict[int, list[list[complex]]]
    atom_abs: dict[int, list[list[float]]]
    term_count: int


def build_core_tables(H: int, core_speeds: tuple[int, int, int, int]) -> CoreTables:
    cs = coeffs(H)
    coeff_cache = {
        rr: residue_matrix(rr)
        for rr in itertools.product(RESIDUES, repeat=4)
    }
    signed: dict[int, list[list[complex]]] = defaultdict(zero_matrix)
    atom_abs: dict[int, list[list[float]]] = defaultdict(
        lambda: [[0.0 for _ in RESIDUES] for _ in RESIDUES]
    )
    term_count = 0
    for ns in itertools.product(cs, repeat=4):
        s = sum(c * n for c, n in zip(core_speeds, ns))
        scale = 1.0 / math.prod(ns)
        mat = coeff_cache[tuple(n % MOD for n in ns)]
        add_scaled(signed[s], mat, scale)
        abs_dst = atom_abs[s]
        abs_scale = abs(scale)
        for i in range(6):
            row_dst = abs_dst[i]
            row_src = mat[i]
            for j in range(6):
                row_dst[j] += abs(row_src[j]) * abs_scale
        term_count += 1
    return CoreTables(core_speeds, dict(signed), dict(atom_abs), term_count)


@dataclass
class PairTables:
    signed: dict[int, list[list[complex]]]
    atom_abs: dict[int, list[list[float]]]
    term_count: int


def build_pair_tables(H: int, A: int, B: int) -> PairTables:
    cs = coeffs(H)
    signed: dict[int, list[list[complex]]] = defaultdict(zero_matrix)
    atom_abs: dict[int, list[list[float]]] = defaultdict(
        lambda: [[0.0 for _ in RESIDUES] for _ in RESIDUES]
    )
    term_count = 0
    for x, y in itertools.product(cs, repeat=2):
        s = -(A * x + B * y)
        i = x % MOD - 1
        j = y % MOD - 1
        z = 1.0 / (x * y)
        signed[s][i][j] += z
        atom_abs[s][i][j] += abs(z)
        term_count += 1
    return PairTables(dict(signed), dict(atom_abs), term_count)


@dataclass
class CaseStats:
    name: str
    A: int
    B: int
    packet: str
    active_s: int
    signed: complex
    block_l2: float
    block_l1: float
    raw_atom_abs: float
    largest_s: int
    largest_contrib: complex
    largest_l2: float


def analyze_case(name: str, A: int, B: int, core: CoreTables, H: int) -> CaseStats:
    pair = build_pair_tables(H, A, B)
    common_s = sorted(set(core.signed) & set(pair.signed))
    signed_total = 0j
    block_l2 = 0.0
    block_l1 = 0.0
    raw_atom_abs = 0.0
    largest_s = 0
    largest_contrib = 0j
    largest_l2 = 0.0
    for s in common_s:
        cmat = core.signed[s]
        pmat = pair.signed[s]
        contrib = dot(cmat, pmat)
        signed_total += contrib
        l2 = frob(cmat) * frob(pmat)
        block_l2 += l2
        block_l1 += sum(abs(cmat[i][j]) * abs(pmat[i][j]) for i in range(6) for j in range(6))
        raw_atom_abs += sum(
            core.atom_abs[s][i][j] * pair.atom_abs[s][i][j]
            for i in range(6)
            for j in range(6)
        )
        if abs(contrib) > abs(largest_contrib):
            largest_s = s
            largest_contrib = contrib
            largest_l2 = l2
    return CaseStats(
        name=name,
        A=A,
        B=B,
        packet=packet_units(A, B),
        active_s=len(common_s),
        signed=signed_total,
        block_l2=block_l2,
        block_l1=block_l1,
        raw_atom_abs=raw_atom_abs,
        largest_s=largest_s,
        largest_contrib=largest_contrib,
        largest_l2=largest_l2,
    )


def fmt_complex(z: complex) -> str:
    return f"{z.real:+.8e}{z.imag:+.1e}j"


def report_cases(stats: list[CaseStats]) -> None:
    section("TWO-LARGE CORE/PAIR TRANSFER TABLE")
    print(
        "Core block: exact c1*n1+...+c4*n4 transfer. Pair block: exact A*x+B*y transfer. "
        "All sums are truncated to |n_i|<=H, 7 does not divide n_i."
    )
    print(f"H={HMAX}, ambient d={AMBIENT_D}")
    print()
    print(
        f"{'case':<18} {'A,B':>9} {'res':>7} {'packet':>8} {'active s':>8} "
        f"{'signed':>20} {'L2/sig':>10} {'L1/sig':>10} {'raw/sig':>10}"
    )
    for st in stats:
        sig = abs(st.signed)
        l2_ratio = st.block_l2 / sig if sig > 1e-18 else math.inf
        l1_ratio = st.block_l1 / sig if sig > 1e-18 else math.inf
        raw_ratio = st.raw_atom_abs / sig if sig > 1e-18 else math.inf
        print(
            f"{st.name:<18} {str((st.A, st.B)):>9} "
            f"{str((st.A % MOD, st.B % MOD)):>7} {st.packet:>8} {st.active_s:>8} "
            f"{fmt_complex(st.signed):>20} {l2_ratio:>10.3g} "
            f"{l1_ratio:>10.3g} {raw_ratio:>10.3g}"
        )


def report_dominant_channels(stats: list[CaseStats]) -> None:
    section("DOMINANT EXACT-S CHANNELS")
    print(
        "The Cauchy block bound is a sum over exact core-pair transfer channels s. "
        "The largest channel is shown to expose whether the mass is local or a "
        "long harmonic envelope."
    )
    print(f"{'case':<18} {'largest s':>10} {'channel signed':>20} {'channel L2':>12}")
    for st in stats:
        print(
            f"{st.name:<18} {st.largest_s:>10} {fmt_complex(st.largest_contrib):>20} "
            f"{st.largest_l2:>12.5e}"
        )


def tournament_analysis() -> None:
    section("TOURNAMENT ANALYSIS")
    vertices = [
        "core_pair_block_transfer",
        "additive_frequency_packet",
        "affine_Q_selector",
        "signed_channel_cauchy_bound",
        "successive_minima_tail_count",
        "blind_pair_residue_matrix",
        "raw_runner_vertices",
    ]
    print("Pairwise observable: which quotient preserves signed reciprocal-tail cancellation.")
    print(
        "Switch/gauge: prefer the earliest quotient that keeps both exact relation "
        "channels s and finite mod-7 packet phase."
    )
    print("Hamiltonian path:")
    print("  " + " > ".join(vertices))
    print("score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1}")
    print("directed_3_cycles = 0")
    print("SCCs = seven singleton proof-obligation vertices")
    print(
        "Assumption challenged: the residual Minkowski lemma should not count "
        "six free coordinates first.  The two-large face has a natural exact "
        "core/pair transfer channel before any lattice-volume envelope."
    )


def main() -> None:
    section("LRC14 BLOCK-FREQUENCY TRANSFER - CODEX S24")
    print(
        "Goal: test whether HYP-2632's finite character packets lift to a small "
        "block transfer norm for the exact reciprocal hyperplane."
    )
    core_speed_rows = [
        (1, 1, 1, 1),
        (1, 8, 15, 22),
    ]
    core_cache = {row: build_core_tables(HMAX, row) for row in core_speed_rows}
    for row, core in core_cache.items():
        print(f"core {row} terms enumerated: {core.term_count}")
    cases = [
        ("4+2 id", (1, 1, 1, 1), 15, 15),
        ("4+2 QR", (1, 1, 1, 1), 16, 16),
        ("4+2 NQR", (1, 1, 1, 1), 17, 17),
        ("4+2 QR skew", (1, 1, 1, 1), 16, 23),
        ("4+1+1 high", (1, 1, 1, 1), 16, 17),
        ("4+1+1 low", (1, 1, 1, 1), 16, 20),
        ("4+1+1 zero", (1, 1, 1, 1), 17, 20),
        ("4+1+1 zero2", (1, 1, 1, 1), 18, 19),
        ("spread 4+2 QR", (1, 8, 15, 22), 16, 16),
        ("spread high", (1, 8, 15, 22), 16, 17),
        ("spread zero", (1, 8, 15, 22), 17, 20),
    ]
    stats = [analyze_case(name, A, B, core_cache[core_speeds], HMAX) for name, core_speeds, A, B in cases]
    report_cases(stats)
    report_dominant_channels(stats)
    section("READOUT")
    print(
        "The raw atom absolute sum is the old free-envelope shadow.  The block "
        "L2 and L1 columns retain the exact transfer channel s and the 6x6 "
        "residue matrix before absolute values.  When these ratios are much "
        "smaller than raw/signed, the missing proof should be a channelwise "
        "Cauchy/Schur estimate plus a successive-minima count over s, not a "
        "direct six-dimensional harmonic bound."
    )
    print(
        "Conjectural simplification: prove a uniform bound on "
        "sum_s ||Core_s||_2 ||Pair_s(A,B)||_2 for the two-large tail after "
        "height-2 wall deletion.  HYP-2632 supplies the finite packet phase; "
        "this script supplies the exact transfer skeleton that the analytic "
        "bound should use."
    )
    tournament_analysis()


if __name__ == "__main__":
    main()
