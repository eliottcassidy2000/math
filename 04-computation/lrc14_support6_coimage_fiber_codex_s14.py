#!/usr/bin/env python3
"""
LRC(14) support-6 coimage fiber atlas.

HYP-2614/S12 attached the residue-address coordinate

    K(n_1,...,n_6,0,...) = C_d(n mod 7)/(n_1...n_6).

HYP-2615/S13 interpreted the small signed mass as a coimage phenomenon.  This
script computes that coimage directly.  For a six-speed support with residues
a_i = e_i mod 7, reduce the relation hyperplane modulo 7:

    a_1 r_1 + ... + a_6 r_6 = 0,     r_i in F_7^*.

The leading residue coimage coefficient is

    S_d(a) = sum_{r in (F_7^*)^6, a.r=0} C_d(r).

Zeros are allowed in the speed residues a_i.  This matters: a support can touch
a speed divisible by 7 while its Fourier coefficient remains nonzero and
7-coprime.  That is exactly a coimage degeneration, because the mod-7 relation
forgets that speed coordinate.

This is a finite atlas, not a proof of LRC(14).  It identifies the finite
residue classes that any final reciprocal-tail theorem must control.
"""
from __future__ import annotations

import itertools
import math
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_support6_residue_cusp_codex_s12 as s12  # noqa: E402

MOD = 7
RESIDUES = tuple(itertools.product(range(1, MOD), repeat=6))


def section(title: str) -> None:
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


def canon_support(vals: tuple[int, ...]) -> tuple[int, ...]:
    """Canonical support residue class modulo scalar and coordinate permutation."""
    residues = tuple(v % MOD for v in vals)
    if all(v == 0 for v in residues):
        return tuple(sorted(residues))
    return min(tuple(sorted((lam * v) % MOD for v in residues)) for lam in range(1, MOD))


def support_classes() -> list[tuple[int, ...]]:
    return sorted(
        {
            canon_support(vals)
            for vals in itertools.product(range(MOD), repeat=6)
            if any(vals)
        }
    )


@dataclass
class FiberStats:
    cls: tuple[int, ...]
    d: int
    count: int
    signed: complex
    abs_sum: float
    max_boundary_abs: float
    max_boundary_key: tuple[int, int]
    max_boundary_signed: complex
    max_boundary_abs_sum: float
    max_boundary_count: int

    @property
    def z(self) -> int:
        return self.cls.count(0)

    @property
    def signed_abs(self) -> float:
        return abs(self.signed)

    @property
    def ratio(self) -> float:
        return self.abs_sum / self.signed_abs if self.signed_abs > 1e-18 else math.inf

    @property
    def boundary_ratio(self) -> float:
        b = abs(self.max_boundary_signed)
        return self.max_boundary_abs_sum / b if b > 1e-18 else math.inf


def compute_stats_for_d(d: int, classes: list[tuple[int, ...]]) -> list[FiberStats]:
    coeffs = [(r, s12.residue_coeff(r, d)) for r in RESIDUES]
    rows: list[FiberStats] = []
    for cls in classes:
        signed = 0j
        abs_sum = 0.0
        count = 0
        boundary_signed: dict[tuple[int, int], complex] = defaultdict(complex)
        boundary_abs: dict[tuple[int, int], float] = defaultdict(float)
        boundary_count: dict[tuple[int, int], int] = defaultdict(int)
        for r, c in coeffs:
            if sum(a * ri for a, ri in zip(cls, r)) % MOD:
                continue
            signed += c
            abs_sum += abs(c)
            count += 1
            for i, u in enumerate(r):
                key = (i, u)
                boundary_signed[key] += c
                boundary_abs[key] += abs(c)
                boundary_count[key] += 1
        if boundary_signed:
            best_key, best_signed = max(
                boundary_signed.items(), key=lambda item: abs(item[1])
            )
            best_abs = boundary_abs[best_key]
            best_count = boundary_count[best_key]
        else:
            best_key = (5, 6)
            best_signed = 0j
            best_abs = 0.0
            best_count = 0
        rows.append(
            FiberStats(
                cls=cls,
                d=d,
                count=count,
                signed=signed,
                abs_sum=abs_sum,
                max_boundary_abs=abs(best_signed),
                max_boundary_key=best_key,
                max_boundary_signed=best_signed,
                max_boundary_abs_sum=best_abs,
                max_boundary_count=best_count,
            )
        )
    return rows


def fmt_complex(z: complex) -> str:
    return f"{z.real:.8g}{z.imag:+.1g}j"


NAMED_SUPPORTS = {
    "AP core / dissociated 211": ((1, 2, 3, 4, 5, 6), 7),
    "resonant 21 support": ((1, 2, 3, 4, 5, 21), 7),
    "k=9 wide 68 support": ((2, 3, 4, 5, 6, 68), 8),
    "k=10 wall 22 support": ((1, 2, 4, 7, 8, 22), 9),
}


def class_inventory(classes: list[tuple[int, ...]]) -> None:
    section("PROJECTIVE SPEED-RESIDUE COIMAGE CLASSES")
    hist = Counter(cls.count(0) for cls in classes)
    print(f"projective support classes modulo F_7^* and S_6: {len(classes)}")
    print(f"zero-speed-residue histogram: {dict(sorted(hist.items()))}")
    print(
        "The 159 classes are the finite coimage table for exact support 6. "
        "A zero here means a support speed is divisible by 7; that coordinate is "
        "still Fourier-live, but disappears from the mod-7 relation."
    )


def summary_by_d(all_stats: dict[int, list[FiberStats]]) -> None:
    section("COIMAGE FIBER SUMMARY BY AMBIENT DIMENSION")
    print(
        f"{'d':>3} {'max |S_d|':>12} {'class':>22} {'z':>2} "
        f"{'zero-ish classes':>16} {'median |S_d|':>15} {'max boundary':>13}"
    )
    for d, rows in all_stats.items():
        max_row = max(rows, key=lambda r: r.signed_abs)
        vals = sorted(r.signed_abs for r in rows)
        median = vals[len(vals) // 2]
        zeroish = sum(r.signed_abs < 1e-12 for r in rows)
        max_boundary = max(rows, key=lambda r: r.max_boundary_abs)
        print(
            f"{d:>3} {max_row.signed_abs:>12.8g} {str(max_row.cls):>22} "
            f"{max_row.z:>2} {zeroish:>16} {median:>15.8g} "
            f"{max_boundary.max_boundary_abs:>13.8g}"
        )


def zero_count_summary(all_stats: dict[int, list[FiberStats]]) -> None:
    section("MAX COIMAGE MASS BY NUMBER OF ZERO SPEED RESIDUES")
    for d in (7, 8, 9, 10, 13):
        rows = all_stats[d]
        print(f"\nd={d}")
        print(f"{'z':>2} {'classes':>7} {'max |S_d|':>12} {'class':>22} {'count':>7} {'abs/signed':>12}")
        for z in range(6):
            zrows = [r for r in rows if r.z == z]
            if not zrows:
                continue
            best = max(zrows, key=lambda r: r.signed_abs)
            print(
                f"{z:>2} {len(zrows):>7} {best.signed_abs:>12.8g} "
                f"{str(best.cls):>22} {best.count:>7} {best.ratio:>12.6g}"
            )
        print(
            "  read: high zero-speed-residue classes are often coimage-null or "
            "coimage-small; divisibility by 7 is a degeneration to ledger, not "
            "automatically an analytic tail threat."
        )


def named_support_table(all_stats: dict[int, list[FiberStats]]) -> None:
    section("NAMED S12 SUPPORTS IN THE COIMAGE ATLAS")
    lookup = {(r.d, r.cls): r for rows in all_stats.values() for r in rows}
    print(
        f"{'support':<28} {'d':>3} {'class':>22} {'z':>2} {'fiber count':>11} "
        f"{'S_d':>18} {'|S_d|':>12} {'abs/signed':>12} {'max boundary':>13}"
    )
    for name, (support, d) in NAMED_SUPPORTS.items():
        cls = canon_support(support)
        row = lookup[(d, cls)]
        print(
            f"{name:<28} {d:>3} {str(cls):>22} {row.z:>2} {row.count:>11} "
            f"{fmt_complex(row.signed):>18} {row.signed_abs:>12.8g} "
            f"{row.ratio:>12.6g} {row.max_boundary_abs:>13.8g}"
        )
    print(
        "The k=10 height-one wall is especially revealing: in its relevant "
        "ambient dimension d=9, its leading coimage fiber is numerically zero "
        "even though its absolute fiber mass is large.  This is the coimage "
        "version of the S12 cusp cancellation signal."
    )


def top_classes(all_stats: dict[int, list[FiberStats]]) -> None:
    section("TOP COIMAGE FIBERS")
    for d in (7, 9, 13):
        print(f"\nd={d}")
        print(f"{'rank':>4} {'class':>22} {'z':>2} {'|S_d|':>12} {'abs total':>12} {'ratio':>12} {'boundary |sum|':>15}")
        rows = sorted(all_stats[d], key=lambda r: r.signed_abs, reverse=True)
        for i, row in enumerate(rows[:10], 1):
            print(
                f"{i:>4} {str(row.cls):>22} {row.z:>2} {row.signed_abs:>12.8g} "
                f"{row.abs_sum:>12.8g} {row.ratio:>12.6g} "
                f"{row.max_boundary_abs:>15.8g}"
            )


def coimage_null_table(all_stats: dict[int, list[FiberStats]]) -> None:
    section("COIMAGE-NULL AND COIMAGE-SMALL CLASSES")
    thresholds = [1e-12, 1e-3, 1e-2, 1e-1]
    print(f"{'d':>3} " + " ".join(f"<{t:g}".rjust(10) for t in thresholds))
    for d, rows in all_stats.items():
        print(
            f"{d:>3} "
            + " ".join(f"{sum(r.signed_abs < t for r in rows):>10}" for t in thresholds)
        )
    classes = [r.cls for r in all_stats[6]]
    null_all = [
        cls
        for cls in classes
        if all(next(r for r in all_stats[d] if r.cls == cls).signed_abs < 1e-12 for d in all_stats)
    ]
    print(f"classes coimage-null for every d=6..13: {len(null_all)}")
    print("first null classes:", null_all[:12])


def tournament_analysis() -> None:
    section("TOURNAMENT ANALYSIS: COIMAGE PROOF QUOTIENTS")
    vertices = [
        "raw_relation_volume",
        "speed_residue_projectivization",
        "coimage_fiber_sum",
        "fixed_boundary_residue",
        "named_wall_nullity",
        "low_height_wall_ledger",
        "reciprocal_tail_estimate",
    ]
    # Higher score means better proof leverage / less ghost mass.
    score = {
        "raw_relation_volume": 0,
        "speed_residue_projectivization": 2,
        "coimage_fiber_sum": 5,
        "fixed_boundary_residue": 3,
        "named_wall_nullity": 6,
        "low_height_wall_ledger": 4,
        "reciprocal_tail_estimate": 1,
    }
    outscore = {v: 0 for v in vertices}
    cycles = 0
    for i, a in enumerate(vertices):
        for b in vertices[i + 1 :]:
            if score[a] > score[b] or (score[a] == score[b] and i < vertices.index(b)):
                outscore[a] += 1
            else:
                outscore[b] += 1
    for a, b, c in itertools.combinations(range(len(vertices)), 3):
        tri = [vertices[a], vertices[b], vertices[c]]
        wins = 0
        for x, y in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            if score[x] > score[y]:
                wins += 1
        if wins in (0, 3):
            cycles += 1
    hist = Counter(outscore.values())
    path = sorted(vertices, key=lambda v: (-score[v], vertices.index(v)))
    print(f"Hamiltonian proof path: {path}")
    print(f"score histogram: {dict(sorted(hist.items()))}")
    print(f"directed 3-cycles: {cycles}")
    print(
        "Assumption challenged: the vertices are neither runners nor residue tuples. "
        "They are quotient stages in the proof obligation.  The map from raw "
        "relation volume to speed-residue projective class to coimage fiber "
        "preserves the analytic predicate while discarding witness-time data."
    )


def main() -> None:
    section("LRC(14) SUPPORT-6 COIMAGE FIBER ATLAS S14")
    print(
        "Goal: compute the finite mod-7 coimage of the support-6 relation "
        "hyperplane, with speed residues allowed to be 0."
    )
    classes = support_classes()
    class_inventory(classes)
    all_stats = {d: compute_stats_for_d(d, classes) for d in range(6, 14)}
    summary_by_d(all_stats)
    zero_count_summary(all_stats)
    named_support_table(all_stats)
    top_classes(all_stats)
    coimage_null_table(all_stats)
    tournament_analysis()
    section("S14 READING")
    print(
        "1. The support-6 residue coimage has only 159 projective speed-residue "
        "classes.  This is the finite address table sitting between S13's "
        "sequence spine and the infinite reciprocal tail."
    )
    print(
        "2. Allowing speed residues 0 mod 7 is essential.  Those coordinates are "
        "Fourier-live but invisible to the mod-7 relation, which is exactly a "
        "coimage degeneration."
    )
    print(
        "3. The leading coimage fiber often kills or shrinks the apparent cusp "
        "mass.  The k=10 height-one wall support is coimage-null in its relevant "
        "ambient dimension d=9, so its danger belongs in the finite wall ledger, "
        "not in the analytic tail."
    )
    print(
        "4. LRC(14) remains open.  The next theorem target is now more concrete: "
        "control the finite list of non-null coimage fiber classes by a signed "
        "reciprocal-tail estimate after deleting low-height walls."
    )


if __name__ == "__main__":
    main()
