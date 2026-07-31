#!/usr/bin/env python3
"""Exact open-component witness for a 12-chart THM-2640 simplex.

The fixed configuration is the first witness from the exhaustive common-
configuration probe.  All interval intersections are exact on the canonical
T-grid.  The final delayed/carry pullback is kept as rational endpoints with
denominator R=13^6, so the output is an actual positive open component rather
than only a positive integral.
"""

from fractions import Fraction
from math import floor
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))
import lrc14_predecessor_carry_private_root_atlas_thm2640 as m


def shift_union(intervals, shift, modulus):
    out = []
    for left, right in intervals:
        length = right - left
        start = (left - shift) % modulus
        stop = start + length
        if stop <= modulus:
            out.append((start, stop))
        else:
            out.append((start, modulus))
            out.append((0, stop - modulus))
    out.sort()
    merged = []
    for left, right in out:
        if merged and left <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], right))
        else:
            merged.append((left, right))
    return merged


def first_delayed_component(base, prefix, carry):
    """First component of base intersect carry/delayed pullback.

    If y={R*x} lies in [a/T,b/T) and the predecessor carry is c, then
    {S*x} lies in [(c+a/T)/13,(c+b/T)/13).  In T-grid coordinates X,
    the pullback components have R-scaled endpoints

        k*(13*T) + c*T + a,  k*(13*T) + c*T + b.
    """
    starts, lens, _ = prefix
    period = 13 * m.T
    best = None
    for left, right, weight in base:
        if not weight:
            continue
        left_scaled = left * m.R
        right_scaled = right * m.R
        for a, length in zip(starts, lens):
            b = a + length
            A = carry * m.T + a
            B = carry * m.T + b
            k = (left_scaled - B) // period + 1
            while k * period + A < right_scaled:
                lo = max(left_scaled, k * period + A)
                hi = min(right_scaled, k * period + B)
                if lo < hi:
                    candidate = (lo, hi)
                    if best is None or candidate < best:
                        best = candidate
                    break
                k += 1
    return best


def frac_part(value):
    return value - floor(value)


def in_interval_union(value, intervals):
    """Exact strict-interior membership; the witness avoids every endpoint."""
    return any(Fraction(left, 1) < value < Fraction(right, 1)
               for left, right in intervals)


def in_weighted_support(value, pieces):
    return any(weight and Fraction(left, 1) < value < Fraction(right, 1)
               for left, right, weight in pieces)


def main():
    cell = (1, 0)
    j = 0
    sector, edge, kappa, h = 0, 0, 0, 6
    carry0 = 0
    carries = tuple(c for c in range(m.P) if c != 6)
    deltas = tuple((2 * (c - carry0)) % m.P for c in carries)

    # Cheap coefficient check on the single selected rail.  The global
    # primitive content 26 is already canonical in THM-2640.
    result = m.shard((j, j + 2))
    metadata, rows = result[5], result[6]
    m.require(metadata == ((1, 0, 12), (1, 0, 0)),
              "selected cell rail metadata changed")
    flags = []
    for c in range(m.P):
        root = (2 * c + (2 * h + kappa) // m.P + (edge == 0)) % m.P
        flags.append(m.is_unit(rows[0][sector][edge][c][kappa][h], root, 26))
    m.require(tuple(c for c, flag in enumerate(flags) if flag) == carries,
              "selected twelve-carry unit configuration changed")

    module, _, _, _, rails, present, starts = m.core.build_carrier_data()
    prefixes = m.build_pair_prefixes(module)
    pieces = rails[j][3]
    rail_support = [(left, right) for left, right, weight in pieces if weight]
    common_rail = list(pieces)
    for delta in deltas:
        if delta == 0:
            continue
        shift = 7 * delta * m.T // m.R
        common_rail = m.old.intersect_weighted_union(
            common_rail, shift_union(rail_support, shift, m.T)
        )
    m.require(common_rail, "twelve translated rail supports became disjoint")

    root0 = 1
    witnesses = []
    base_by_clock = {}
    for ell5 in range(m.Q7):
        common = m.old.intersect_weighted_union(
            common_rail, present[ell5, (-h) % m.P],
            starts[ell5, (-h) % m.P]
        )
        source_present = present[ell5, (-h) % m.P]
        for delta in deltas:
            if delta == 0 or not common:
                continue
            shift = 7 * delta * m.T // m.R
            common = m.old.intersect_weighted_union(
                common, shift_union(source_present, shift, m.T)
            )
        if not common:
            continue
        half = m.old.intersect_weighted_comb(
            common, module.C3, 182, 14 * root0 - 13, 14 * root0
        )
        base_by_clock[ell5] = half
        component = first_delayed_component(
            half, prefixes[sector][ell5][h][kappa], carry0
        )
        if component is not None:
            witnesses.append((ell5, component, len(half)))

    m.require(witnesses, "no rational open component survived")
    m.require(len(witnesses) == 2, "positive delayed-clock component count changed")
    ell5, (lo, hi), base_count = min(witnesses, key=lambda item: item[1])
    left = Fraction(lo, m.R * m.T)
    right = Fraction(hi, m.R * m.T)
    midpoint = (left + right) / 2
    m.require(left < midpoint < right, "component is not positive")
    m.require((ell5, base_count, left, right, right - left, midpoint) == (
        1,
        545,
        Fraction(4855397, 10396204),
        Fraction(23436073938185, 50180491033036),
        Fraction(3, 12545122758259),
        Fraction(3348010562597, 7168641576148),
    ), "canonical twelve-chart component changed")

    # Independent pointwise replay of all physical factors at the rational
    # midpoint.  This does not reuse the multi-intersection routine above.
    direct = []
    ell_prefix = prefixes[sector][ell5][h][kappa]
    delayed_intervals = tuple(
        (start, start + length)
        for start, length in zip(ell_prefix[0], ell_prefix[1])
    )
    for delta in range(m.P):
        translated = frac_part(midpoint + Fraction(7 * delta, m.R))
        grid_x = translated * m.T
        carry = floor(m.P * frac_part(m.S * translated))
        y = frac_part(m.R * translated)
        digit = floor(2 * m.P * y)
        observed_h, observed_kappa = divmod(digit, 2)
        root = (root0 + delta) % m.P
        rail_ok = in_weighted_support(grid_x, pieces)
        present_ok = in_interval_union(
            grid_x, present[ell5, (-h) % m.P]
        )
        delayed_ok = in_interval_union(y * m.T, delayed_intervals)
        if root:
            lo, hi = 14 * root - 13, 14 * root
            deep = frac_part(module.C3 * translated)
            deep_ok = Fraction(lo, 182) < deep < Fraction(hi, 182)
        else:
            deep_ok = False
        expected_carry = (carry0 + 7 * delta) % m.P
        unit_ok = root != 0 and flags[expected_carry]
        ok = (carry == expected_carry and observed_h == h
              and observed_kappa == kappa and rail_ok and present_ok
              and delayed_ok and deep_ok and unit_ok)
        direct.append((delta, carry, root, ok))
    m.require(tuple(delta for delta, _, _, ok in direct if ok)
              == tuple(sorted(deltas)),
              f"independent midpoint chart-label replay changed: {direct}")

    # Test the omitted thirteenth label against the FULL configuration union
    # on the same base cell.  This is exactly where the fixed-configuration
    # root-zero argument can fail: adjacent deep half-tooths overlap, so the
    # missing label may switch edge/root and/or rail/sector.
    missing_delta = 12
    translated = frac_part(midpoint + Fraction(7 * missing_delta, m.R))
    grid_x = translated * m.T
    carry = floor(m.P * frac_part(m.S * translated))
    y = frac_part(m.R * translated)
    digit = floor(2 * m.P * y)
    observed_h, observed_kappa = divmod(digit, 2)
    m.require((carry, observed_h, observed_kappa) == (6, h, kappa),
              "missing-label physical digits changed")
    union_witnesses = []
    for jj in (0, 1):
        jj_pieces = rails[jj][3]
        rail_ok = in_weighted_support(grid_x, jj_pieces)
        if not rail_ok:
            continue
        for alt_sector in range(2):
            for alt_edge in range(2):
                alt_root = (2 * carry + (2 * h + kappa) // m.P
                            + (alt_edge == 0)) % m.P
                if not alt_root:
                    continue
                unit_ok = m.is_unit(
                    rows[jj][alt_sector][alt_edge][carry][kappa][h],
                    alt_root, 26,
                )
                if not unit_ok:
                    continue
                if alt_edge == 0:
                    lo, hi = 14 * alt_root - 13, 14 * alt_root
                else:
                    lo, hi = 14 * alt_root, 14 * alt_root + 13
                deep = frac_part(module.C3 * translated)
                deep_ok = Fraction(lo, 182) < deep < Fraction(hi, 182)
                if not deep_ok:
                    continue
                for alt_ell5 in range(m.Q7):
                    present_ok = in_interval_union(
                        grid_x, present[alt_ell5, (-h) % m.P]
                    )
                    alt_prefix = prefixes[alt_sector][alt_ell5][h][kappa]
                    alt_delayed = tuple(
                        (start, start + length)
                        for start, length in zip(alt_prefix[0], alt_prefix[1])
                    )
                    delayed_ok = in_interval_union(y * m.T, alt_delayed)
                    if present_ok and delayed_ok:
                        union_witnesses.append(
                            (jj, alt_sector, alt_edge, kappa, h,
                             alt_ell5, carry, alt_root)
                        )

    # Exhaust every open component of the same 12-fold packet against every
    # allowed rail/sector/edge/clock realization of the omitted label.  This
    # is stronger than the midpoint check and is still a finite exact sweep.
    extension_witnesses = []
    missing_shift = 7 * missing_delta * m.T // m.R
    shifted_alt_rails = {
        jj: shift_union(
            [(a, b) for a, b, weight in rails[jj][3] if weight],
            missing_shift, m.T,
        )
        for jj in (0, 1)
    }
    for source_ell5, source_base in sorted(base_by_clock.items()):
        source_prefix = prefixes[sector][source_ell5][h][kappa]
        source_delayed = tuple(
            (start, start + length)
            for start, length in zip(source_prefix[0], source_prefix[1])
        )
        for jj in (0, 1):
            alt_rail = m.old.intersect_weighted_union(
                source_base, shifted_alt_rails[jj]
            )
            if not alt_rail:
                continue
            for alt_sector in range(2):
                for alt_edge in range(2):
                    alt_root = (2 * carry + (2 * h + kappa) // m.P
                                + (alt_edge == 0)) % m.P
                    if not alt_root or not m.is_unit(
                        rows[jj][alt_sector][alt_edge][carry][kappa][h],
                        alt_root, 26,
                    ):
                        continue
                    base_root = (alt_root - missing_delta) % m.P
                    if alt_edge == 0:
                        lo, hi = 14 * base_root - 13, 14 * base_root
                    else:
                        lo, hi = 14 * base_root, 14 * base_root + 13
                    alt_deep = module.make_comb(module.C3, 182, lo, hi)
                    deep_base = m.old.intersect_weighted_union(
                        alt_rail, alt_deep
                    )
                    if not deep_base:
                        continue
                    for alt_ell5 in range(m.Q7):
                        shifted_alt_present = shift_union(
                            present[alt_ell5, (-h) % m.P],
                            missing_shift, m.T,
                        )
                        full_base = m.old.intersect_weighted_union(
                            deep_base, shifted_alt_present
                        )
                        if not full_base:
                            continue
                        alt_prefix = prefixes[alt_sector][alt_ell5][h][kappa]
                        alt_delayed = tuple(
                            (start, start + length)
                            for start, length in zip(alt_prefix[0], alt_prefix[1])
                        )
                        common_delayed = m.old.intersect_sorted(
                            list(source_delayed), list(alt_delayed)
                        )
                        if not common_delayed:
                            continue
                        common_prefix = module.make_prefix(common_delayed)
                        extension = first_delayed_component(
                            full_base, common_prefix, carry0
                        )
                        if extension is not None:
                            extension_witnesses.append(
                                (source_ell5, jj, alt_sector, alt_edge,
                                 alt_ell5, carry, alt_root, base_root,
                                 extension)
                            )

    m.require(not union_witnesses,
              "selected midpoint gained a missing-label configuration")
    m.require(not extension_witnesses,
              "selected open facet gained a missing-label extension")

    print("THM-2640 twelve-chart open-component witness")
    print(f"cell={cell} rail={j} sector_edge_kappa_h="
          f"{(sector, edge, kappa, h)}")
    print(f"carry_labels={carries} missing_carry=6")
    print(f"target_deltas={deltas}")
    print(f"unit_flags={tuple(flags)}")
    print(f"positive_clock_components={len(witnesses)} first_clock={ell5}")
    print(f"first_clock_base_half_components={base_count}")
    print(f"component_left={left}")
    print(f"component_right={right}")
    print(f"component_length={right-left}")
    print(f"component_midpoint={midpoint}")
    print(f"direct_midpoint_delta_carry_root_ok={tuple(direct)}")
    print(f"missing_delta_full_union_witnesses={tuple(union_witnesses)}")
    print(f"missing_delta_full_component_extension_witnesses="
          f"{tuple(extension_witnesses[:20])}; count={len(extension_witnesses)}")
    print("scope_boundary=this proves one explicit nonextendable open "
          "12-fold component inside cell (1,0); it does not exclude a "
          "13-fold intersection elsewhere after configuration switching")


if __name__ == "__main__":
    main()
