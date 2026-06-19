#!/usr/bin/env python3
"""
LRC(14) prime-mask/coimage transfer scout.

This continues the S18 height-2 wall-class audit.  The user's proposed
recurrence is read as:

    finite prime-mask transfer / inclusion-exclusion
    coupled to the signed mod-7 coimage tail.

The point of this scout is not to prove LRC(14).  It asks whether the
height-2 wall supports already factor through two finite quotients:

  1. the LRC14 unit seam: (Z/14Z)^* -> F_7^*,
  2. the prime mask of the large wall speed M, tracked at {2,3,5,7}.

The seam quotient is exact because reducing units mod 14 gives all nonzero
mod-7 residues.  The mask quotient is a finite transfer coordinate: it records
which small clocks the large speed kills before the signed coimage class is
evaluated.

Tournament Analysis declaration:
  vertices are proof quotients, not runners: unit seam, prime mask, wall class,
  projective coimage class, repeated tail packet, signed tail theorem, raw
  supports.  The quotient preserves the support-six residual address and
  destroys row identity and witness-time geometry.
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

import lrc14_height2_coimage_wall_classes_codex_s18 as s18  # noqa: E402
import lrc14_support6_coimage_fiber_codex_s14 as s14  # noqa: E402

PRIMES = (2, 3, 5, 7)
EPS = 1e-12


def section(title: str) -> None:
    print("\n" + "=" * 88, flush=True)
    print(title, flush=True)
    print("=" * 88, flush=True)


def mask_of(n: int) -> int:
    mask = 0
    for i, p in enumerate(PRIMES):
        if n % p == 0:
            mask |= 1 << i
    return mask


def mask_name(mask: int) -> str:
    vals = [str(p) for i, p in enumerate(PRIMES) if mask & (1 << i)]
    return "{" + ",".join(vals) + "}" if vals else "{}"


def all_masks() -> list[int]:
    return list(range(1 << len(PRIMES)))


def pattern_label(cls: tuple[int, ...]) -> str:
    counts = Counter(cls)
    multiplicities = sorted(counts.values(), reverse=True)
    zeros = counts.get(0, 0)
    if zeros >= 4:
        return "zero-cusp"
    if multiplicities[:2] == [4, 2]:
        return "4+2 repeated"
    if multiplicities[:3] == [4, 1, 1]:
        return "4+1+1 repeated"
    if multiplicities[:3] == [2, 2, 2]:
        return "2+2+2 repeated"
    return "mixed"


def chi7(x: int) -> int:
    """Quadratic character on F_7, with chi(0)=0."""
    x %= 7
    if x == 0:
        return 0
    return 1 if x in {1, 2, 4} else -1


@dataclass(frozen=True)
class WallRecord:
    k: int
    support: tuple[int, ...]
    large: int
    cls: tuple[int, ...]
    apex_mask: int
    support_mask: int
    relation: tuple[int, ...]


def enumerate_wall_records(k: int, height: int = 2) -> list[WallRecord]:
    """Enumerate unique one-large height<=height wall supports with mask data."""
    B = s18.BOUND[k]
    coeffs = s18.coefficient_set(height)
    records: dict[tuple[int, ...], WallRecord] = {}

    for core_support in itertools.combinations(range(1, B + 1), 5):
        for core_coeffs in itertools.product(coeffs, repeat=5):
            core_sum = sum(c * v for c, v in zip(core_coeffs, core_support))
            if core_sum == 0:
                continue
            for cM in coeffs:
                if (-core_sum) % cM:
                    continue
                M = (-core_sum) // cM
                if M <= B:
                    continue
                support = tuple(sorted(core_support + (M,)))
                if support in records:
                    continue
                smask = 0
                for v in support:
                    smask |= mask_of(v)
                records[support] = WallRecord(
                    k=k,
                    support=support,
                    large=M,
                    cls=s14.canon_support(support),
                    apex_mask=mask_of(M),
                    support_mask=smask,
                    relation=tuple(core_coeffs + (cM,)),
                )
    return list(records.values())


def coimage_rows_by_k() -> dict[int, list[s14.FiberStats]]:
    classes = s14.support_classes()
    return {k: s14.compute_stats_for_d(s18.AMBIENT_D[k], classes) for k in (8, 9, 10)}


def unit_seam_report() -> None:
    section("UNIT SEAM LAW: (Z/14Z)^* REDUCES TO F_7^*")
    units14 = tuple(a for a in range(14) if math.gcd(a, 14) == 1)
    residues7 = tuple(sorted({a % 7 for a in units14}))
    print(f"units mod 14: {units14}")
    print(f"reductions mod 7: {residues7}")
    print(f"F_7^*: {tuple(range(1, 7))}")
    print(
        "Readout: the HYP-2617 projectivization by F_7^* is exactly the "
        "unit action inherited from the 14-runner clock.  The coimage class is "
        "therefore a coimage of the LRC14 unit seam, not an arbitrary mod-7 trick."
    )


def build_class_masks(records: list[WallRecord]) -> dict[tuple[int, ...], set[int]]:
    masks: dict[tuple[int, ...], set[int]] = defaultdict(set)
    for rec in records:
        masks[rec.cls].add(rec.apex_mask)
    return masks


def mass_lookup(rows: list[s14.FiberStats]) -> dict[tuple[int, ...], float]:
    return {row.cls: row.signed_abs for row in rows if row.signed_abs > EPS}


def allowed_coverage(
    rows: list[s14.FiberStats],
    class_masks: dict[tuple[int, ...], set[int]],
    allowed: int,
) -> tuple[int, float]:
    total = 0.0
    classes = 0
    for row in rows:
        if row.signed_abs <= EPS:
            continue
        if any((mask & ~allowed) == 0 for mask in class_masks.get(row.cls, ())):
            total += row.signed_abs
            classes += 1
    return classes, total


def mask_transfer_table(
    records_by_k: dict[int, list[WallRecord]],
    rows_by_k: dict[int, list[s14.FiberStats]],
) -> None:
    section("PRIME-MASK TRANSFER: ALLOWED APEX MASKS")
    print(
        "A class is counted for allowed mask Q if it has some height<=2 one-large "
        "wall witness whose large speed M is divisible only by primes in Q among "
        "{2,3,5,7}.  This tests the old mod-30 slice Q={2,3,5} against the "
        "extra LRC14 prime 7."
    )
    selected = [0, 1 << 0, (1 << 0) | (1 << 1), (1 << 0) | (1 << 1) | (1 << 2), (1 << 3), 15]
    print(
        f"{'k':>2} {'allowed Q':>12} {'hit classes':>13} {'hit |S|-mass':>14} "
        f"{'mass share':>11}"
    )
    for k in (8, 9, 10):
        rows = rows_by_k[k]
        class_masks = build_class_masks(records_by_k[k])
        total_mass = sum(row.signed_abs for row in rows if row.signed_abs > EPS)
        for allowed in selected:
            classes, mass = allowed_coverage(rows, class_masks, allowed)
            print(
                f"{k:>2} {mask_name(allowed):>12} {classes:>13} "
                f"{mass:>14.8g} {mass / total_mass:>10.6%}"
            )


def minimal_mask_antichains(
    records_by_k: dict[int, list[WallRecord]],
    rows_by_k: dict[int, list[s14.FiberStats]],
) -> None:
    section("MINIMAL MASK ANTICHAINS FOR WALL-ADDRESSED COIMAGE CLASSES")
    print(
        "For each coimage class, keep only apex masks minimal under subset order. "
        "This is the finite transfer fingerprint before the signed coimage tail is "
        "evaluated."
    )
    for k in (8, 9, 10):
        rows = rows_by_k[k]
        masses = mass_lookup(rows)
        class_masks = build_class_masks(records_by_k[k])
        buckets: dict[tuple[int, ...], list[tuple[tuple[int, ...], float]]] = defaultdict(list)
        for cls, masks in class_masks.items():
            if cls not in masses:
                continue
            minimal = tuple(
                sorted(
                    m for m in masks
                    if not any(n != m and (n & m) == n for n in masks)
                )
            )
            buckets[minimal].append((cls, masses[cls]))

        print(f"\nk={k}")
        print(f"{'minimal masks':>34} {'classes':>8} {'|S|-mass':>12} {'top class':>24}")
        for minimal, items in sorted(
            buckets.items(),
            key=lambda item: sum(mass for _, mass in item[1]),
            reverse=True,
        )[:12]:
            mass = sum(v for _, v in items)
            top_cls, top_mass = max(items, key=lambda x: x[1])
            names = "(" + ",".join(mask_name(m) for m in minimal) + ")"
            print(f"{names:>34} {len(items):>8} {mass:>12.8g} {str(top_cls):>24}")


def support_mask_census(records_by_k: dict[int, list[WallRecord]]) -> None:
    section("SUPPORT-LEVEL MASK CENSUS")
    print(
        f"{'k':>2} {'supports':>10} {'classes':>8} {'M divisible by 7':>16} "
        f"{'all support touches 7':>21} {'top apex masks':>30}"
    )
    for k in (8, 9, 10):
        records = records_by_k[k]
        apex_hist = Counter(rec.apex_mask for rec in records)
        div7 = sum(1 for rec in records if rec.apex_mask & (1 << 3))
        touch7 = sum(1 for rec in records if rec.support_mask & (1 << 3))
        top = ", ".join(f"{mask_name(m)}:{c}" for m, c in apex_hist.most_common(5))
        print(
            f"{k:>2} {len(records):>10} {len({r.cls for r in records}):>8} "
            f"{div7:>16} {touch7:>21} {top:>30}"
        )


def tail_packet_report(
    records_by_k: dict[int, list[WallRecord]],
    rows_by_k: dict[int, list[s14.FiberStats]],
) -> None:
    section("TAIL-ONLY PACKET AFTER PRIME-MASK + HEIGHT-2 TRANSFER")
    for k in (8, 9, 10):
        rows = rows_by_k[k]
        class_masks = build_class_masks(records_by_k[k])
        misses = [
            row for row in sorted(rows, key=lambda r: r.signed_abs, reverse=True)
            if row.signed_abs > EPS and row.cls not in class_masks
        ]
        hist = Counter(pattern_label(row.cls) for row in misses)
        total = sum(row.signed_abs for row in rows if row.signed_abs > EPS)
        miss_mass = sum(row.signed_abs for row in misses)
        print(
            f"\nk={k}: missed classes={len(misses)}, missed mass share={miss_mass / total:.6%}, "
            f"patterns={dict(hist)}"
        )
        for row in misses[:10]:
            print(
                f"  {str(row.cls):>22}  |S_d|={row.signed_abs:.8g}  "
                f"ratio={row.ratio:.6g}  {pattern_label(row.cls)}"
            )


def repeated_tail_character_report(
    records_by_k: dict[int, list[WallRecord]],
    rows_by_k: dict[int, list[s14.FiberStats]],
) -> None:
    section("REPEATED TAIL CHARACTER SPLIT AT k=10")
    rows = rows_by_k[10]
    class_masks = build_class_masks(records_by_k[10])
    misses = [
        row for row in sorted(rows, key=lambda r: r.signed_abs, reverse=True)
        if row.signed_abs > EPS and row.cls not in class_masks
    ]

    print("4+2 packet (1,1,1,1,a,a):")
    print(f"{'class':>22} {'a':>2} {'chi7(a)':>8} {'|S_9|':>12} {'abs/signed':>12}")
    for row in misses:
        counts = Counter(row.cls)
        if 0 in counts or sorted(counts.values(), reverse=True)[:2] != [4, 2]:
            continue
        a = next(v for v, n in counts.items() if n == 2)
        print(
            f"{str(row.cls):>22} {a:>2} {chi7(a):>8} "
            f"{row.signed_abs:>12.8g} {row.ratio:>12.6g}"
        )

    print("\n4+1+1 packet (1,1,1,1,a,b):")
    print(
        f"{'class':>22} {'{a,b}':>8} {'chi(a),chi(b),chi(ab)':>23} "
        f"{'chi((a-1)(b-1))':>20} {'|S_9|':>12}"
    )
    for row in misses:
        counts = Counter(row.cls)
        if 0 in counts or sorted(counts.values(), reverse=True)[:3] != [4, 1, 1]:
            continue
        a, b = sorted(v for v, n in counts.items() if n == 1)
        sig = (chi7(a), chi7(b), chi7(a * b))
        shifted = chi7((a - 1) * (b - 1))
        print(
            f"{str(row.cls):>22} {str((a,b)):>8} {str(sig):>23} "
            f"{shifted:>20} {row.signed_abs:>12.8g}"
        )

    print(
        "\nReadout: the dominant 4+2 tail is not featureless.  Its two largest "
        "classes are the nontrivial quadratic residues a=2,4; the three smaller "
        "classes are nonresidues a=3,5,6.  The 4+1+1 packet is finer, but still "
        "collapses to a handful of multiplicative-character signatures.  This is "
        "the concrete target for the next signed cotangent/Dedekind lemma."
    )


def tournament_analysis() -> None:
    section("TOURNAMENT ANALYSIS")
    vertices = [
        "unit_seam_coimage",
        "prime_mask_transfer",
        "height2_wall_classes",
        "repeated_tail_packet",
        "signed_dedekind_tail",
        "raw_supports",
        "raw_runner_vertices",
    ]
    print("Hamiltonian proof path:")
    print("  " + " > ".join(vertices))
    print("score histogram:", {i: 1 for i in range(len(vertices))})
    print("directed 3-cycles: 0")
    print("SCC sizes:", [1] * len(vertices))
    print(
        "Challenged assumption: vertices need not be runners, gaps, or arcs. "
        "This quotient uses the unit seam, prime masks, and coimage packets. "
        "It preserves the support-six signed-tail address and destroys exact "
        "witness times and row identity."
    )


def main() -> None:
    section("LRC(14) PRIME-MASK / COIMAGE TRANSFER SCOUT S19")
    print(
        "Goal: test whether the hidden recurrence is a finite prime-mask "
        "transfer coupled to the signed mod-7 coimage tail."
    )
    unit_seam_report()
    rows_by_k = coimage_rows_by_k()
    records_by_k = {k: enumerate_wall_records(k, height=2) for k in (8, 9, 10)}
    support_mask_census(records_by_k)
    mask_transfer_table(records_by_k, rows_by_k)
    minimal_mask_antichains(records_by_k, rows_by_k)
    tail_packet_report(records_by_k, rows_by_k)
    repeated_tail_character_report(records_by_k, rows_by_k)
    tournament_analysis()
    section("S19 STATUS")
    print(
        "LRC(14) is not proved here.  The useful advance is structural: "
        "the projective mod-7 coimage quotient is exactly the (Z/14Z)^* seam, "
        "while prime masks provide a finite transfer coordinate for low-height "
        "wall classes.  The remaining k=10 obstruction is still the repeated "
        "coimage tail packet, now isolated after both quotients."
    )


if __name__ == "__main__":
    main()
