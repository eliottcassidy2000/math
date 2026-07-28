#!/usr/bin/env python3
"""Exact scout for one D-edge followed by a slope-seven configuration switch.

The first leg is the physical THM-2680 handoff ``D(x)={13x}``.  At the
following coordinate ``z=D(x)`` we retain one of THM-2672's honest
same-configuration slope-seven facets.  The script asks for a rational point
strictly inside

    E_current(x) AND E_following(Dx)
      AND intersection_delta T_delta^-1 P_(e,c+7 delta)(Dx).

Every factor is evaluated before integration.  A positive answer is only a
support witness for a mixed correspondence: the slope-seven leg is not the
canonical THM-2365 target action and supplies no owner/endpoint transport.
"""

from fractions import Fraction
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_dilation_reversed_clock_fibre_product_probe as d
import lrc14_predecessor_carry_private_root_atlas_thm2640 as m
import lrc14_slope7_rebase_facet_torsor_thm2672 as facet


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def strict_in_union(value, intervals):
    return any(Fraction(left) < value < Fraction(right)
               for left, right in intervals)


def strict_in_weighted(value, pieces):
    return any(weight and Fraction(left) < value < Fraction(right)
               for left, right, weight in pieces)


def frac(value):
    return value - (value.numerator // value.denominator)


def point_in_atom(atom, x, prefixes):
    """Strict-interior membership in one unintegrated THM-2680 atom."""
    base = frac(x) * m.T
    if not strict_in_weighted(base, atom["pieces"]):
        return False
    y = frac(m.R * x) * m.T
    prefix = prefixes[atom["future"]][atom["h"]]
    delayed = tuple((start, start + length)
                    for start, length in zip(prefix[0], prefix[1]))
    return strict_in_union(y, delayed)


def delayed_components(base, prefix, speed, grid):
    """All exact open components of base AND Q({speed*x})."""
    starts, lengths, _ = prefix
    out = []
    for left, right, weight in base:
        if not weight:
            continue
        base_left = left * speed
        base_right = right * speed
        for start, length in zip(starts, lengths):
            stop = start + length
            branch = (base_left - stop) // grid + 1
            while branch * grid + start < base_right:
                lo = max(base_left, branch * grid + start)
                hi = min(base_right, branch * grid + stop)
                if lo < hi:
                    out.append((Fraction(lo, speed * grid),
                                Fraction(hi, speed * grid)))
                branch += 1
    return tuple(sorted(set(out)))


def selected_two_edge_components(module, prefixes, rails, present, starts):
    """Exact components of THM-2680's first displayed positive atom pair."""
    current = d.build_atom(
        module, prefixes, present, starts, rails[2], 3, 2, 0, 0
    )
    following = d.build_atom(
        module, prefixes, present, starts, rails[0], 1, 6, 1, 1
    )
    require((current["j"], current["h"], following["j"], following["h"])
            == (5, 2, 2, 6), "selected D-edge labels changed")
    base = d.intersect_weighted(
        d.scale_weighted(current["pieces"], m.P),
        d.pullback_dilation_weighted(following["pieces"], m.T),
    )
    current_q = d.prefix_intervals(prefixes[3][2])
    following_q = d.preimage_times_13(
        d.prefix_intervals(prefixes[1][6]), m.T
    )
    joint_q = m.old.intersect_sorted(current_q, following_q)
    joint_prefix = module.make_prefix(
        [(m.P * left, m.P * right) for left, right in joint_q]
    )
    components = delayed_components(
        base, joint_prefix, m.R, m.P * m.T
    )
    require(components, "selected positive D edge lost every open component")
    return current, following, components


def point_in_private_configuration(module, pair_prefixes, rails, present,
                                   rail_index, sector, edge, local_kappa, h,
                                   shallow, carry, root, z):
    """Strict membership in the full THM-2640 packet with displayed labels."""
    base = frac(z) * m.T
    if not strict_in_weighted(base, rails[rail_index][3]):
        return False
    if not strict_in_union(base, present[shallow, (-h) % m.P]):
        return False
    deep = frac(module.C3 * z) * 182
    lo, hi = ((14 * root - 13, 14 * root)
              if edge == 0 else (14 * root, 14 * root + 13))
    if not (Fraction(lo) < deep < Fraction(hi)):
        return False
    y = frac(m.R * z)
    carry_actual = ((m.R * z).numerator // (m.R * z).denominator) % m.P
    if carry_actual != carry:
        return False
    prefix = pair_prefixes[sector][shallow][h][local_kappa]
    delayed = tuple((start, start + length)
                    for start, length in zip(prefix[0], prefix[1]))
    return strict_in_union(y * m.T, delayed)


def find_mixed_triangle(module, prefixes, rails, present, starts):
    """Search one proved D-edge component for a two-chart slope switch."""
    current, following, components = selected_two_edge_components(
        module, prefixes, rails, present, starts
    )
    pair_prefixes = m.build_pair_prefixes(module)
    shard = m.shard((0, 1))
    require(shard[1] == 26 and shard[5] == ((1, 0, 12),),
            "selected THM-2640 rail/content changed")
    rows = shard[6][0]
    tested = 0
    for left, right in components:
        x = (left + right) / 2
        z = frac(m.P * x)
        require(point_in_atom(current, x, prefixes)
                and point_in_atom(following, z, prefixes),
                "D-edge component midpoint left one of its endpoint atoms")
        carry = ((m.R * z).numerator // (m.R * z).denominator) % m.P
        y = frac(m.R * z)
        h = min(m.P - 1, (m.P * y.numerator) // y.denominator)
        local_kappa = min(1, (2 * m.P * y.numerator) // y.denominator - 2 * h)
        shallow = following["future"]
        for sector in range(2):
            for edge in range(2):
                root = (2 * carry + (2 * h + local_kappa) // m.P
                        + (edge == 0)) % m.P
                if not root or not m.is_unit(
                        rows[sector][edge][carry][local_kappa][h], root, 26):
                    continue
                if not point_in_private_configuration(
                        module, pair_prefixes, rails, present, 0, sector, edge,
                        local_kappa, h, shallow, carry, root, z):
                    continue
                successful_deltas = []
                for delta in range(1, m.P):
                    tested += 1
                    shifted = frac(z + Fraction(7 * delta, m.R))
                    carry2 = (carry + 7 * delta) % m.P
                    root2 = (root + delta) % m.P
                    if not root2 or not m.is_unit(
                            rows[sector][edge][carry2][local_kappa][h],
                            root2, 26):
                        continue
                    if point_in_private_configuration(
                            module, pair_prefixes, rails, present, 0,
                            sector, edge, local_kappa, h, shallow,
                            carry2, root2, shifted):
                        successful_deltas.append((delta, carry2, root2, shifted))
                if successful_deltas:
                    return {
                            "component": (left, right),
                            "x": x,
                            "z": z,
                            "D_labels": (3, 1, 0),
                            "D_current": (1, 2, 5, 2, 0, 0),
                            "D_following": (1, 0, 2, 6, 1, 1),
                            "private_config": (
                                0, sector, edge, local_kappa, h, shallow,
                            ),
                            "carry_root": (carry, root),
                            "successful_slope_switches": tuple(successful_deltas),
                            "components_scanned": components.index((left, right)) + 1,
                            "slope_candidates_tested": tested,
                        }
    return None


def build_selected_facets():
    """Return one exact first component for every source carry."""
    rail_index = 0
    sector, edge, local_kappa, h = 0, 0, 0, 6
    module, _, _, _, rails, present, _ = m.core.build_carrier_data()
    pair_prefixes = m.build_pair_prefixes(module)
    pieces = rails[rail_index][3]
    rail_support = [(left, right) for left, right, weight in pieces if weight]
    roots = tuple(
        (2 * c + (2 * h + local_kappa) // m.P + (edge == 0)) % m.P
        for c in range(m.P)
    )
    missing_carry = 6
    rows = []
    for source_carry in range(m.P):
        missing_delta = 2 * (missing_carry - source_carry) % m.P
        deltas = tuple(delta for delta in range(m.P)
                       if delta != missing_delta)
        anchor = deltas[0]
        anchor_shift = 7 * anchor * m.T // m.R
        common_rail = facet.shift_weighted(pieces, anchor_shift, m.T)
        for delta in deltas[1:]:
            shift = 7 * delta * m.T // m.R
            common_rail = m.old.intersect_weighted_union(
                common_rail, facet.shift_union(rail_support, shift, m.T)
            )
            if not common_rail:
                break
        require(common_rail, "selected facet lost its common rail")

        for shallow in range(m.Q7):
            source_present = present[shallow, (-h) % m.P]
            common = m.old.intersect_weighted_union(
                common_rail,
                facet.shift_union(source_present, anchor_shift, m.T),
            )
            for delta in deltas[1:]:
                if not common:
                    break
                shift = 7 * delta * m.T // m.R
                common = m.old.intersect_weighted_union(
                    common, facet.shift_union(source_present, shift, m.T)
                )
            if not common:
                continue
            half = facet.intersect_root_half(
                common, module.C3, edge, roots[source_carry]
            )
            component = facet.first_delayed_component(
                half,
                pair_prefixes[sector][shallow][h][local_kappa],
                source_carry,
            )
            if component is None:
                continue
            left_raw, right_raw = component
            rows.append({
                "source_carry": source_carry,
                "missing_delta": missing_delta,
                "deltas": deltas,
                "root": roots[source_carry],
                "shallow": shallow,
                "left": Fraction(left_raw, m.R * m.T),
                "right": Fraction(right_raw, m.R * m.T),
            })
            break
    require(len(rows) == m.P, "did not recover all thirteen selected facets")
    return rows


def main():
    module, prefixes, _, _, rails, present, starts = m.core.build_carrier_data()
    rails_by_cell = {}
    for rail_index, rail in enumerate(rails):
        rails_by_cell.setdefault(rail[:2], []).append((rail_index, rail))

    mixed_triangle = find_mixed_triangle(
        module, prefixes, rails, present, starts
    )
    require(mixed_triangle is not None,
            "no two-chart mixed D/slope-seven physical simplex found")

    facets = build_selected_facets()
    witnesses = []
    for row in facets:
        z = (row["left"] + row["right"]) / 2
        following_atoms = []
        # THM-2672's selected rail metadata is (source,owner,deep)=(1,0,12).
        # Search the sharp THM-2680 refinements of that same rail and shallow
        # clock, then retain only factors containing the exact facet midpoint.
        for global_kappa in (0, 1):
            atom = d.build_atom(
                module, prefixes, present, starts, rails[0],
                row["shallow"], 6, 1, global_kappa,
            )
            if atom["value"] and point_in_atom(atom, z, prefixes):
                atom = dict(atom)
                atom["rail_index"] = 0
                following_atoms.append(atom)

        for following in following_atoms:
            middle_clock = row["shallow"]
            next_clock = 0
            for branch in range(m.P):
                x = (z + branch) / m.P
                require(frac(m.P * x) == z,
                        "chosen inverse branch does not map to the facet")
                for previous_clock in range(m.Q7):
                    if previous_clock == middle_clock:
                        continue
                    for source in range(1, m.P):
                        for rail_index, rail in rails_by_cell.get(
                                (source, middle_clock), ()):
                            for epsilon in (0, 1):
                                for global_kappa in (0, 1):
                                    current = d.build_atom(
                                        module, prefixes, present, starts, rail,
                                        previous_clock, following["j"],
                                        epsilon, global_kappa,
                                    )
                                    if (current["value"]
                                            and point_in_atom(
                                                current, x, prefixes)):
                                        witnesses.append({
                                            "source_carry": row["source_carry"],
                                            "missing_delta": row["missing_delta"],
                                            "slope_label_count": len(row["deltas"]),
                                            "root": row["root"],
                                            "clock_path": (
                                                previous_clock,
                                                middle_clock,
                                                next_clock,
                                            ),
                                            "branch": branch,
                                            "x": x,
                                            "z": z,
                                            "facet": (row["left"], row["right"]),
                                            "current": (
                                                source, rail_index,
                                                current["j"], current["h"],
                                                epsilon, global_kappa,
                                            ),
                                            "following": (
                                                1, 0, following["j"],
                                                following["h"], 1,
                                                global_kappa,
                                            ),
                                        })
                                        break
                                if witnesses:
                                    break
                            if witnesses:
                                break
                        if witnesses:
                            break
                    if witnesses:
                        break
                if witnesses:
                    break
            if witnesses:
                break
        if witnesses:
            break

    # The stronger twelve-label spot-check is allowed to fail.  It tests only
    # the selected first-component midpoint in each carry stratum, not every
    # component of the canonical facet bank.  The scoped question here needs
    # only one nonzero slope-switch edge.
    witness = witnesses[0] if witnesses else None
    print("LRC14 MIXED D/SLOPE-SEVEN CONFIGURATION-SWITCH SCOUT")
    print("object=E_current(x) AND E_following(Dx) AND two THM2640 slope-seven charts at Dx")
    print(f"two_chart_witness={mixed_triangle}")
    print(f"selected_twelve_chart_first_midpoint_spotcheck={witness}")
    print("strict_membership=all factors hold at rational interior points; hence a positive open component exists")
    print("preserved=physical common point, D owner/shallow covariance, D digit covariance, rail/source labels, predecessor carry, future half-digit, and private root across one nonzero slope label")
    print("lost=slope leg is a translated packet correspondence, not THM2365 target action or canonical owner/endpoint transition")
    print("verdict=PASS: mixed support bypass exists; D-chain nilpotence is grammar-specific, not a universal two-simplex obstruction")


if __name__ == "__main__":
    main()
