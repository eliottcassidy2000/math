#!/usr/bin/env python3
"""
LRC(14) alternating cusp sequence atlas.

The prompt asks for the integer/fractional sequences behind the key empirical
fact:

    large absolute mass on boundary/cusp faces,
    but tiny signed mass after the seven-sector kernel.

This script makes the "alternating series" mechanism explicit.  It separates
three layers:

  1. residue coefficient signs C_d(r), where all-positive summation would be
     the wrong divergent/bosonic object;
  2. shell-by-shell signed reciprocal sums on relation hyperplanes, where
     boundary/cusp absolute mass collapses to a much smaller signed shell;
  3. projective mod-7 coimage fibers, where many support classes are null or
     small after quotienting by scalar and coordinate permutation.

This is route-finding, not a proof of LRC(14).  The new proof target is a
two-stage conditional convergence theorem: finite low-height wall ledger, then
signed reciprocal-tail summation over non-null coimage fibers.

Tournament Analysis declaration.
  Vertex set: proof-obligation sequence families, not runners:
    coefficient_sign_inventory, shell_alternation_ladder, boundary_cusp_faces,
    projective_coimage_nullity, height_one_wall_ledger, dedekind_tail_bound,
    raw_absolute_volume.
  Pairwise observable: proof leverage, exactness, portability in d/k, and how
    much ghost absolute mass the quotient removes.  The switch is lexicographic
    on that tuple.

Assumption challenge.
  I considered runner vertices, residue tuples, shell heights, boundary faces,
  Fourier modes, coimage classes, and proof obligations.  Runner vertices hide
  the support-floor cancellation.  Residue tuples are too fine.  The chosen
  quotient preserves the analytic LRC predicate "the support-6 correction stays
  below cap margin after wall deletion" while discarding witness-time geometry.
"""
from __future__ import annotations

import itertools
import math
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_support6_coimage_fiber_codex_s14 as s14  # noqa: E402
import lrc14_support6_residue_cusp_codex_s12 as s12  # noqa: E402

MOD = 7
RESIDUES = tuple(itertools.product(range(1, MOD), repeat=6))
TOL = 1e-12


def section(title: str) -> None:
    print("\n" + "=" * 92)
    print(title)
    print("=" * 92)


def real_part(z: complex) -> float:
    if abs(z.imag) > 1e-8:
        raise ValueError(f"unexpected imaginary residue {z}")
    return z.real


def neg_residue_tuple(r: tuple[int, ...]) -> tuple[int, ...]:
    return tuple((-x) % MOD for x in r)


def sign_symbol(x: float) -> str:
    if x > TOL:
        return "+"
    if x < -TOL:
        return "-"
    return "0"


def fmt_float(x: float) -> str:
    if math.isinf(x):
        return "inf"
    return f"{x:.8g}"


def coefficient_sign_inventory(max_d: int = 16) -> None:
    section("RESIDUE COEFFICIENT SIGN INVENTORY")
    print(
        "C_d(r) is the finite mod-7 coefficient in "
        "K(n_1,...,n_6,0,...)=C_d(n mod 7)/(n_1...n_6)."
    )
    print(
        "Individual C_d(r) are complex, so the signed real atoms are conjugate "
        "pairs r and -r.  The table compares the signed paired total with the "
        "all-positive paired total.  Large abs/net is the finite residue version "
        "of an alternating series that would be unusable with all signs made "
        "positive."
    )
    print(
        f"{'d':>3} {'+ pairs':>8} {'- pairs':>8} {'0':>5} {'sum+':>13} "
        f"{'sum- abs':>13} {'net':>13} {'abs total':>13} {'abs/net':>11} "
        f"{'max |pair|':>12} {'max |C|':>10}"
    )
    max_seq: list[float] = []
    ratio_seq: list[float] = []
    balance_seq: list[float] = []
    for d in range(6, max_d + 1):
        vals = []
        max_raw = 0.0
        seen: set[tuple[int, ...]] = set()
        for r in RESIDUES:
            if r in seen:
                continue
            nr = neg_residue_tuple(r)
            c = s12.residue_coeff(r, d)
            cn = s12.residue_coeff(nr, d)
            max_raw = max(max_raw, abs(c), abs(cn))
            vals.append(real_part(c + cn))
            seen.add(r)
            seen.add(nr)
        pos = [x for x in vals if x > TOL]
        neg = [x for x in vals if x < -TOL]
        zero = len(vals) - len(pos) - len(neg)
        sum_pos = sum(pos)
        sum_neg_abs = -sum(neg)
        net = sum(vals)
        abs_total = sum(abs(x) for x in vals)
        ratio = abs_total / abs(net) if abs(net) > TOL else math.inf
        max_abs = max(abs(x) for x in vals)
        max_seq.append(max_abs)
        ratio_seq.append(ratio)
        balance_seq.append(sum_pos / sum_neg_abs if sum_neg_abs else math.inf)
        print(
            f"{d:>3} {len(pos):>8} {len(neg):>8} {zero:>5} "
            f"{sum_pos:>13.7g} {sum_neg_abs:>13.7g} {net:>13.7g} "
            f"{abs_total:>13.7g} {fmt_float(ratio):>11} {max_abs:>12.7g} "
            f"{max_raw:>10.7g}"
        )
    print("max |paired C_d| sequence:", ", ".join(f"{x:.9g}" for x in max_seq))
    print("abs/net sequence:", ", ".join(fmt_float(x) for x in ratio_seq))
    print("positive/negative balance sequence:", ", ".join(fmt_float(x) for x in balance_seq))
    print(
        "Readout: the residue coefficient is already an alternating object.  "
        "The positive and negative piles stay comparable while the max coefficient "
        "decays; this is the finite shadow of the signed theta tail."
    )


@dataclass
class ShellRow:
    h: int
    count: int
    signed: complex
    abs_mass: float

    @property
    def signed_abs(self) -> float:
        return abs(self.signed)

    @property
    def ratio(self) -> float:
        return self.abs_mass / self.signed_abs if self.signed_abs > TOL else math.inf

    @property
    def sign(self) -> str:
        return sign_symbol(real_part(self.signed))


CASES = [
    ("AP core support", (1, 2, 3, 4, 5, 6), 7, 28),
    ("resonant 21 support", (1, 2, 3, 4, 5, 21), 7, 28),
    ("dissociated 211 support", (2, 3, 4, 5, 6, 211), 7, 34),
    ("k=9 wide 68 support", (2, 3, 4, 5, 6, 68), 8, 30),
    ("k=10 wall 22 support", (1, 2, 4, 7, 8, 22), 9, 24),
]


@lru_cache(None)
def cached_shell_data(
    vals: tuple[int, ...], d: int, hmax: int
) -> tuple[dict[int, s12.ShellStats], dict[tuple[int, int], s12.ShellStats]]:
    return s12.six_support_shells(vals, d, hmax)


def shell_rows(vals: tuple[int, ...], d: int, hmax: int) -> list[ShellRow]:
    by_shell, _ = cached_shell_data(vals, d, hmax)
    rows = []
    for h in range(1, hmax + 1):
        stat = by_shell.get(h)
        if stat and stat.count:
            rows.append(ShellRow(h, stat.count, stat.signed, stat.absK))
    return rows


def run_lengths(signs: list[str]) -> list[int]:
    nonzero = [s for s in signs if s != "0"]
    if not nonzero:
        return []
    out = []
    cur = nonzero[0]
    length = 1
    for s in nonzero[1:]:
        if s == cur:
            length += 1
        else:
            out.append(length)
            cur = s
            length = 1
    out.append(length)
    return out


def shell_alternation_ladders() -> None:
    section("RELATION-HYPERPLANE SHELL ALTERNATION LADDERS")
    print(
        "Each case is an exact support-6 relation hyperplane.  raw abs is the "
        "all-positive/bosonic mass; signed variation is sum_h |shell signed_h|; "
        "net is |sum_h shell signed_h|.  Thus raw/net decomposes into "
        "(raw/signed-variation) times (signed-variation/net)."
    )
    print(
        f"{'case':<28} {'shells':>6} {'sign flips':>10} {'run lengths':>18} "
        f"{'raw abs':>11} {'signed var':>11} {'net':>11} "
        f"{'raw/var':>10} {'var/net':>10} {'raw/net':>10}"
    )
    summary: list[tuple[str, float, float, float]] = []
    for name, vals, d, hmax in CASES:
        rows = shell_rows(vals, d, hmax)
        raw_abs = sum(r.abs_mass for r in rows)
        signed_var = sum(r.signed_abs for r in rows)
        net = abs(sum((r.signed for r in rows), 0j))
        signs = [r.sign for r in rows]
        flips = sum(1 for a, b in zip(signs, signs[1:]) if a != "0" and b != "0" and a != b)
        runs = run_lengths(signs)
        raw_var = raw_abs / signed_var if signed_var > TOL else math.inf
        var_net = signed_var / net if net > TOL else math.inf
        raw_net = raw_abs / net if net > TOL else math.inf
        summary.append((name, raw_var, var_net, raw_net))
        print(
            f"{name:<28} {len(rows):>6} {flips:>10} {str(runs[:8]):>18} "
            f"{raw_abs:>11.6g} {signed_var:>11.6g} {net:>11.6g} "
            f"{fmt_float(raw_var):>10} {fmt_float(var_net):>10} {fmt_float(raw_net):>10}"
        )
    print("\nCompact shell sign words, by height:")
    for name, vals, d, hmax in CASES:
        rows = shell_rows(vals, d, hmax)
        word = "".join(r.sign for r in rows)
        heights = [r.h for r in rows]
        print(f"{name:<28} heights={heights}")
        print(f"{'':<28} signs  ={word}")
    print("\nTwo-stage cancellation factor sequences:")
    print("raw/signed-variation:", ", ".join(f"{x[1]:.6g}" for x in summary))
    print("signed-variation/net:", ", ".join(fmt_float(x[2]) for x in summary))
    print("raw/net:", ", ".join(fmt_float(x[3]) for x in summary))
    print(
        "Readout: the proof should not try to make raw abs small.  It should "
        "prove a cusp-to-signed-shell collapse and then an alternating-shell "
        "summation bound."
    )


def boundary_face_sequences() -> None:
    section("BOUNDARY-FACE INTEGER/FRACTIONAL SEQUENCES")
    print(
        "The final nonempty shell is where the max-norm cusp is most visible.  "
        "touch=t means t coordinates sit on the max boundary.  Counts are integer "
        "mass; ratios are the fractional cancellation address."
    )
    final_counts: list[int] = []
    final_ratios: list[float] = []
    one_touch_counts: list[int] = []
    one_touch_ratios: list[float] = []
    for name, vals, d, hmax in CASES:
        by_shell, by_touch = cached_shell_data(vals, d, hmax)
        final_h = max(h for h, stat in by_shell.items() if stat.count)
        total = by_shell[final_h]
        total_ratio = total.absK / abs(total.signed) if abs(total.signed) > TOL else math.inf
        one = by_touch.get((final_h, 1), s12.ShellStats())
        one_ratio = one.absK / abs(one.signed) if abs(one.signed) > TOL else math.inf
        final_counts.append(total.count)
        final_ratios.append(total_ratio)
        one_touch_counts.append(one.count)
        one_touch_ratios.append(one_ratio)
        print(f"\n{name} final shell h={final_h}")
        print(f"{'touch':>5} {'count':>12} {'signed':>14} {'abs':>11} {'abs/signed':>12}")
        for touch in range(1, 7):
            stat = by_touch.get((final_h, touch), s12.ShellStats())
            ratio = stat.absK / abs(stat.signed) if abs(stat.signed) > TOL else math.inf
            print(
                f"{touch:>5} {stat.count:>12} {real_part(stat.signed):>14.6g} "
                f"{stat.absK:>11.6g} {fmt_float(ratio):>12}"
            )
    print("\nfinal-shell relation count sequence:", final_counts)
    print("final-shell abs/signed sequence:", ", ".join(fmt_float(x) for x in final_ratios))
    print("one-face relation count sequence:", one_touch_counts)
    print("one-face abs/signed sequence:", ", ".join(fmt_float(x) for x in one_touch_ratios))
    print(
        "Readout: the large integers are not evidence of analytic danger by "
        "themselves.  They are the pre-quotient cusp count; the signed fraction "
        "is the proof-relevant address."
    )


def precompute_fibers(classes: list[tuple[int, ...]]) -> dict[tuple[int, ...], list[int]]:
    fibers: dict[tuple[int, ...], list[int]] = {}
    for cls in classes:
        idxs = []
        for i, r in enumerate(RESIDUES):
            if sum(a * ri for a, ri in zip(cls, r)) % MOD == 0:
                idxs.append(i)
        fibers[cls] = idxs
    return fibers


def coimage_sequence_extension(max_d: int = 16) -> None:
    section("PROJECTIVE COIMAGE SEQUENCE EXTENSION")
    classes = s14.support_classes()
    fibers = precompute_fibers(classes)
    print(
        f"Projective classes: {len(classes)}.  Zero-speed-residue histogram: "
        f"{dict(sorted(Counter(cls.count(0) for cls in classes).items()))}."
    )
    print(
        f"{'d':>3} {'max |S_d|':>12} {'argmax class':>22} {'null':>6} "
        f"{'<.001':>6} {'<.01':>6} {'<.1':>6} {'median':>11} {'same sign top?':>14}"
    )
    max_seq: list[float] = []
    null_seq: list[int] = []
    small01_seq: list[int] = []
    for d in range(6, max_d + 1):
        coeffs = [s12.residue_coeff(r, d) for r in RESIDUES]
        rows = []
        for cls in classes:
            signed = sum((coeffs[i] for i in fibers[cls]), 0j)
            rows.append((cls, abs(signed), real_part(signed)))
        rows.sort(key=lambda x: x[1], reverse=True)
        vals = sorted(x[1] for x in rows)
        max_seq.append(rows[0][1])
        null_count = sum(v < 1e-12 for _, v, _ in rows)
        null_seq.append(null_count)
        small01 = sum(v < 0.01 for _, v, _ in rows)
        small01_seq.append(small01)
        top_signs = "".join(sign_symbol(x[2]) for x in rows[:5])
        print(
            f"{d:>3} {rows[0][1]:>12.8g} {str(rows[0][0]):>22} "
            f"{null_count:>6} {sum(v < 0.001 for _, v, _ in rows):>6} "
            f"{small01:>6} {sum(v < 0.1 for _, v, _ in rows):>6} "
            f"{vals[len(vals)//2]:>11.6g} {top_signs:>14}"
        )
    print("max coimage fiber sequence:", ", ".join(f"{x:.9g}" for x in max_seq))
    print("coimage-null count sequence:", null_seq)
    print("coimage <0.01 count sequence:", small01_seq)
    print(
        "Readout: after d=11 the exact-null floor stabilizes at the three pure "
        "degenerate classes, while most classes stay <0.1.  The tail theorem can "
        "be class-by-class rather than volume-by-volume."
    )


def tournament_analysis() -> None:
    section("TOURNAMENT ANALYSIS: SEQUENCE QUOTIENTS")
    names = [
        "raw_absolute_volume",
        "coefficient_sign_inventory",
        "shell_alternation_ladder",
        "boundary_cusp_faces",
        "projective_coimage_nullity",
        "height_one_wall_ledger",
        "dedekind_tail_bound",
    ]
    metrics = {
        # proof leverage, exactness, portability, removed ghost mass
        "raw_absolute_volume": (0, 4, 4, 0),
        "coefficient_sign_inventory": (3, 5, 5, 2),
        "shell_alternation_ladder": (5, 4, 4, 5),
        "boundary_cusp_faces": (4, 4, 4, 4),
        "projective_coimage_nullity": (5, 5, 5, 5),
        "height_one_wall_ledger": (4, 5, 3, 5),
        "dedekind_tail_bound": (6, 1, 5, 6),
    }
    adj = [[False] * len(names) for _ in names]
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if i == j:
                continue
            adj[i][j] = metrics[a] > metrics[b] or (metrics[a] == metrics[b] and i < j)
    scores = {names[i]: sum(adj[i]) for i in range(len(names))}
    cycles = 0
    for a, b, c in itertools.combinations(range(len(names)), 3):
        if adj[a][b] and adj[b][c] and adj[c][a]:
            cycles += 1
        if adj[a][c] and adj[c][b] and adj[b][a]:
            cycles += 1
    path = sorted(names, key=lambda n: scores[n], reverse=True)
    print("Hamiltonian proof path:", path)
    print("score histogram:", dict(sorted(Counter(scores.values()).items())))
    print("directed 3-cycles:", cycles)
    print(
        "Assumption challenged: the vertices are proof quotients.  This preserves "
        "the signed-tail predicate and the recursive d/k sequences while throwing "
        "away runner-time geometry."
    )


def final_reading() -> None:
    section("S15 READING")
    print(
        "1. The all-positive residue and shell sums are the wrong object, just as "
        "an alternating conditionally convergent series becomes divergent when "
        "all signs are made positive."
    )
    print(
        "2. The useful bound should factor into two inequalities: cusp absolute "
        "mass -> signed shell variation, then alternating shell variation -> net."
    )
    print(
        "3. The persistent sequences are now visible: support floor 0..0,1; "
        "residue sign-balance by d; shell sign words; boundary count/ratio "
        "pairs; coimage null-count and max-fiber ladders."
    )
    print(
        "4. LRC(14) is still open.  The next proof move is a class-by-class "
        "Dedekind/cotangent summation theorem over the non-null coimage fibers, "
        "after HYP-2616-style wall ledger deletion."
    )


def main() -> None:
    section("LRC(14) ALTERNATING CUSP SEQUENCE ATLAS S15")
    print(
        "Goal: expose the recursive integer/fractional sequences behind large "
        "absolute cusp mass but tiny signed mass."
    )
    coefficient_sign_inventory()
    shell_alternation_ladders()
    boundary_face_sequences()
    coimage_sequence_extension()
    tournament_analysis()
    final_reading()


if __name__ == "__main__":
    main()
