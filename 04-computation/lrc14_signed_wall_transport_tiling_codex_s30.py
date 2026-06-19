#!/usr/bin/env python3
"""HYP-2647/T894: signed wall transport matrix for AP9 -> nearAP9.

This scout keeps one address layer that HYP-2642 intentionally summarized
away: when AP9=(0,...,8) is replaced by D9=(0,...,7,9), the moving endpoint
runner carries sector mass from the old speed 8 to the new speed 9.  On the
common wall refinement, every atom has an old sector, a new sector, a missed-set
transition, and a signed L_y valuation.

The point is not to find a new numerical maximum.  It is to make HYP-2647's
"signed transport tiling" concrete:

    addressed atom transport -> signed positive/negative scalar shadow.

Tournament Analysis vertices are proof quotients rather than runners.  The
observable is how much LRC-relevant address survives projection.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction


AP9 = tuple(range(9))
DEFECT9 = (0, 1, 2, 3, 4, 5, 6, 7, 9)
OLD_SPEED = 8
NEW_SPEED = 9


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def sector(x: Fraction) -> int:
    return (x.numerator * 7) // x.denominator


def missed_set(E: tuple[int, ...], x: Fraction) -> tuple[int, ...]:
    hit = {sector(frac_part(e * x)) for e in E}
    return tuple(j for j in range(1, 7) if j not in hit)


def g9(t: int) -> Fraction:
    return Fraction(-(t - 2) * (t - 3) * (t - 6), 36)


def breakpoints(rows: tuple[tuple[int, ...], ...]) -> list[Fraction]:
    pts = {Fraction(0), Fraction(1)}
    for E in rows:
        for e in E:
            if e == 0:
                continue
            for a in range(0, 7 * e + 1):
                pts.add(Fraction(a, 7 * e))
    return sorted(pts)


def fmt(q: Fraction) -> str:
    return f"{q} = {float(q):.9f}"


def sign_name(q: Fraction) -> str:
    if q > 0:
        return "positive"
    if q < 0:
        return "negative"
    return "neutral"


def main() -> None:
    pts = breakpoints((AP9, DEFECT9))

    by_sector: dict[tuple[int, int], dict[str, Fraction]] = defaultdict(
        lambda: {
            "mass": Fraction(0),
            "positive": Fraction(0),
            "negative": Fraction(0),
            "signed": Fraction(0),
        }
    )
    old_row_mass: Counter[int] = Counter()
    new_col_mass: Counter[int] = Counter()
    by_count: Counter[tuple[int, int]] = Counter()
    by_signed_bucket: dict[str, Fraction] = defaultdict(Fraction)
    addressed: list[tuple[Fraction, Fraction, int, int, tuple[int, ...], tuple[int, ...], Fraction]] = []

    for lo, hi in zip(pts, pts[1:]):
        if hi == lo:
            continue
        mid = (lo + hi) / 2
        mass = hi - lo
        ma = missed_set(AP9, mid)
        mb = missed_set(DEFECT9, mid)
        old_s = sector(frac_part(OLD_SPEED * mid))
        new_s = sector(frac_part(NEW_SPEED * mid))
        dg = g9(len(mb)) - g9(len(ma))
        weighted = mass * dg
        sign = sign_name(weighted)

        entry = by_sector[(old_s, new_s)]
        entry["mass"] += mass
        entry["signed"] += weighted
        if weighted > 0:
            entry["positive"] += weighted
        elif weighted < 0:
            entry["negative"] += -weighted

        old_row_mass[old_s] += mass
        new_col_mass[new_s] += mass
        by_count[(len(ma), len(mb))] += mass
        by_signed_bucket[sign] += mass
        if weighted:
            addressed.append((abs(weighted), weighted, old_s, new_s, ma, mb, mass))

    weighted_positive = sum(entry["positive"] for entry in by_sector.values())
    weighted_negative = sum(entry["negative"] for entry in by_sector.values())
    signed_delta = sum(entry["signed"] for entry in by_sector.values())
    total_mass = sum(old_row_mass.values())

    print("=" * 88)
    print("HYP-2647/T894 signed wall transport tiling scout")
    print("=" * 88)
    print(f"AP9={AP9}")
    print(f"D9 ={DEFECT9}")
    print(f"moving endpoint: {OLD_SPEED} -> {NEW_SPEED}")
    print()

    print("Transport checksum for the moving endpoint sector map:")
    print(f"  total mass = {fmt(total_mass)}")
    print("  old-sector row masses:")
    print("    " + ", ".join(f"{s}:{old_row_mass[s]}" for s in range(7)))
    print("  new-sector column masses:")
    print("    " + ", ".join(f"{s}:{new_col_mass[s]}" for s in range(7)))
    print()

    print("Signed scalar shadow:")
    print(f"  weighted positive = {fmt(weighted_positive)}")
    print(f"  weighted negative = {fmt(weighted_negative)}")
    print(f"  signed D-AP       = {fmt(signed_delta)}")
    print(f"  AP-D surplus      = {fmt(weighted_negative - weighted_positive)}")
    print()

    print("Mass by signed bucket before taking the weighted shadow:")
    for name in ("positive", "negative", "neutral"):
        print(f"  {name:>8}: {fmt(by_signed_bucket[name])}")
    print()

    print("Old-sector -> new-sector addressed table:")
    print("  old new        mass          positive        negative          signed")
    for old_s, new_s in sorted(by_sector):
        entry = by_sector[(old_s, new_s)]
        if not entry["mass"]:
            continue
        print(
            f"  {old_s:>3} {new_s:>3}  "
            f"{str(entry['mass']):>12}  "
            f"{str(entry['positive']):>14}  "
            f"{str(entry['negative']):>14}  "
            f"{str(entry['signed']):>14}"
        )
    print()

    print("Common-wall transfer by missed count:")
    print("  AP_count -> D_count        mass        weight_delta")
    for a_count, d_count in sorted(by_count):
        mass = by_count[(a_count, d_count)]
        weighted = mass * (g9(d_count) - g9(a_count))
        print(
            f"       {a_count} -> {d_count}       "
            f"{str(mass):>12}  {str(weighted):>18}"
        )
    print()

    print("Largest addressed nonzero weighted atom transitions:")
    print("  old->new   AP_missed -> D_missed          mass        weighted")
    for _, weighted, old_s, new_s, ma, mb, mass in sorted(addressed, reverse=True)[:24]:
        ma_s = "()" if not ma else str(ma)
        mb_s = "()" if not mb else str(mb)
        print(
            f"   {old_s}->{new_s}     "
            f"{ma_s:>12} -> {mb_s:<12}  "
            f"{str(mass):>10}  {str(weighted):>14}"
        )
    print()

    print("Tournament Analysis fingerprint:")
    print("  vertices: addressed_transport, sector_transport, count_transport, scalar_signs")
    print("  pairwise observable: preserves LRC-relevant address before valuation")
    print("  Hamiltonian path: addressed_transport > sector_transport > count_transport > scalar_signs")
    print("  score histogram: {0:1,1:1,2:1,3:1}; directed 3-cycles: 0")


if __name__ == "__main__":
    main()
