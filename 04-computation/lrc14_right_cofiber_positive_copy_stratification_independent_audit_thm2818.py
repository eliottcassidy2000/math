#!/usr/bin/env python3
"""Independent hostile audit of the THM-2771 all-cell copy decomposition.

This companion imports only pinned canonical computation modules.  It does
not import the earlier ``.scratch/lrc_cofiber_copy_bridge`` companions and
does not encode their cell list, copy counts, block lengths, weights,
ancestry digests, or Bockstein values.
"""

from __future__ import annotations

from bisect import bisect_right
from collections import Counter
from hashlib import sha256
from math import gcd
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


PINNED = {
    "lrc14_fully_marked_root_zero_target_profile_thm2749.py":
        "d67c852c52f88feaadb2fcaa0a9a07a212f2e47018040b455855df886200595e",
    "lrc14_root_zero_clutch_mayer_vietoris_wing_shear_thm2751.py":
        "25cbed38026d61891173c687006250a69fe38aea56d67439406bd8bb60fa2552",
    "lrc14_root_zero_full_target_semantic_clutch_20260728.py":
        "208f71020efa19fa47f66d2da061ab03fa7bc87beeb077b4008c069f499736d8",
    "lrc14_full_arm_orbit_path_sheet_audit_thm2791.py":
        "1e00b6711db0d878fca70047f5f1532518084dbf6928551cd28fe51283fde543",
}
for name, expected in PINNED.items():
    actual = sha256(
        (COMP / name).read_bytes().replace(b"\r\n", b"\n")
    ).hexdigest()
    require(actual == expected, f"pinned dependency changed: {name}")


import lrc14_fully_marked_root_zero_target_profile_thm2749 as marked
import lrc14_full_arm_orbit_path_sheet_audit_thm2791 as ancestry
import lrc14_root_zero_clutch_mayer_vietoris_wing_shear_thm2751 as wing
import lrc14_root_zero_full_target_semantic_clutch_20260728 as physical


P = physical.P
T = physical.T
SHIFT = physical.SHIFT
DEPTH = P**5
RAIL_INDEX = 8
SOURCE_ROLE = 3
SOURCE_LABEL = 0


def weighted_intersection(weighted, intervals):
    return tuple(
        physical.relative.private.old.intersect_weighted_union(
            weighted, intervals
        )
    )


def physical_setup():
    (
        module, _prefixes, _whole, _masses, rails, present, _starts
    ) = physical.relative.lift.m.core.build_carrier_data()
    pair_prefixes = physical.relative.private.build_pair_prefixes(module)
    _sv, _tv, rail_pairs, details = physical.overlap.overlap_vectors(
        module, pair_prefixes, rails, present, rail_index=RAIL_INDEX
    )
    require(
        all(left == right for _a, _b, left, right in rail_pairs),
        "selected rail stopped having aligned weights",
    )

    full_module = physical.target.load_present_module()
    require(full_module.T == module.T, "physical grids diverged")
    e3 = physical.target.exclusive_source(full_module, SOURCE_ROLE)
    fork = physical.target.deepest_fork(full_module)
    clocks = tuple(
        full_module.make_comb(
            full_module.C1,
            182,
            26 * clock - 13,
            26 * clock + 13,
        )
        for clock in range(7)
    )
    q_pairs = physical.q_restricted_pair_prefixes(
        full_module, pair_prefixes, fork
    )

    marked_pairs = marked.private.build_pair_prefixes(module)
    delayed = marked.marked_prefixes(
        module,
        marked_pairs,
        marked.two.deepest_fork(module),
    )
    _source_weight, target_weight, rail_common = marked.rail_data(
        rails, marked.RAIL
    )
    return (
        module,
        rails,
        present,
        details,
        full_module,
        e3,
        clocks,
        q_pairs,
        delayed,
        target_weight,
        rail_common,
    )


def physical_right(details, full_module, e3, clocks, clock, target):
    section = physical.target.source_present_section(
        full_module,
        e3,
        clock,
        SOURCE_LABEL,
        target,
        clocks,
    )
    source_base, target_base = details[clock]
    source = weighted_intersection(source_base, section)
    target_native = weighted_intersection(target_base, section)
    target_pullback = physical.overlap.shift_weighted(
        target_native, -SHIFT
    )
    aligned = physical.overlap.intersect_weighted_profiles(
        source, target_pullback
    )
    common = tuple(
        (left, right, source_weight)
        for left, right, source_weight, target_weight in aligned
        if source_weight == target_weight
    )
    require(
        all(
            source_weight == target_weight
            for _left, _right, source_weight, target_weight in aligned
        ),
        f"unequal overlap weights at {(clock, target)}",
    )
    return physical.subtract_weighted(target_pullback, common)


def raw_right(module, rail_common, clocks, clock, target):
    source = marked.two.exclusive_source(module, SOURCE_ROLE)
    raw = tuple(marked.two.intersect_sorted(source, clocks[clock]))
    raw = tuple(
        module.subtract_comb(raw, module.W[1], 182, -13, 13)
    )
    raw = tuple(
        module.subtract_comb(raw, module.C2, 182, -13, 13)
    )
    raw = tuple(
        module.subtract_comb(
            raw,
            module.W[2],
            182,
            -14 * target - 13,
            -14 * target + 13,
        )
    )
    raw = tuple(
        module.subtract_comb(
            raw,
            module.C3,
            182,
            14 * target - 13,
            14 * target + 13,
        )
    )
    left = marked.intersect(rail_common, raw)
    right = marked.intersect(
        rail_common, marked.shift_union(raw, SHIFT)
    )
    return wing.difference(right, left)


def semantic_partition(pieces, q_pair):
    source_cache = {}
    target_cache = {}
    live = []
    dead = []
    hostile = []
    for piece in pieces:
        left, right, _weight = piece
        unit = ((left, right, 1),)
        target_unit = physical.overlap.shift_weighted(unit, SHIFT)
        source_value = physical.relative.private.delayed_carry_pair(
            unit, q_pair, source_cache
        )[12][1]
        target_value = physical.relative.private.delayed_carry_pair(
            target_unit, q_pair, target_cache
        )[6][1]
        if source_value == target_value and source_value != 0:
            live.append(piece)
        elif source_value == 0 and target_value == 0:
            dead.append(piece)
        else:
            hostile.append((piece, source_value, target_value))
    require(
        not hostile,
        f"semantic source/target mismatch: {hostile[:2]}",
    )
    return tuple(live), tuple(dead)


def ancestry_setup():
    e_intervals = tuple(
        ancestry.base.build_set(
            ancestry.base.PAT_E, ancestry.base.ZELL
        )
    )
    q_intervals = tuple(
        ancestry.base.build_set(
            ancestry.host.PAT_QB, ancestry.base.ZELL
        )
    )
    q_depth, q_starts = ancestry.scaled_intervals(
        q_intervals, DEPTH
    )
    e_depth_pack, e_depth_pack_starts = ancestry.scaled_intervals(
        e_intervals, DEPTH * P**2
    )
    e_depth, e_depth_starts = ancestry.scaled_intervals(
        e_intervals, DEPTH
    )
    arguments = (
        q_depth,
        q_starts,
        e_depth_pack,
        e_depth_pack_starts,
        e_depth,
        e_depth_starts,
    )
    u_events = tuple(
        sorted(
            set(
                ancestry.mapped_events(q_intervals, DEPTH)
                + ancestry.mapped_events(
                    e_intervals, DEPTH * P**2
                )
            )
        )
    )
    v_events = ancestry.mapped_events(
        e_intervals,
        DEPTH,
        ancestry.RAIL_DISPLACEMENT * T // P,
    )
    return arguments, u_events, v_events


def wall_chamber(events, left, right):
    index = bisect_right(events, left)
    require(
        index > 0 and index < len(events),
        "ancestry chamber wrapped",
    )
    require(
        events[index] > right,
        f"ancestry wall crossed hull [{left},{right})",
    )
    return events[index - 1], events[index]


def ancestry_record(pieces, arguments, u_events, v_events, cache):
    require(
        tuple(sorted(pieces)) == pieces,
        "physical pieces stopped being ordered",
    )
    hull = pieces[0][0], pieces[-1][1]
    chambers = (
        wall_chamber(u_events, *hull),
        wall_chamber(v_events, *hull),
    )
    if chambers not in cache:
        coordinate = (pieces[0][0] + pieces[0][1]) // 2
        u_labels, v_labels = ancestry.contributor_sets(
            coordinate, *arguments
        )
        cache[chambers] = (u_labels, v_labels)
    u_labels, v_labels = cache[chambers]
    path_a, path_b, path_e = ancestry.SUPPLIED_PATH
    path_active = (
        path_a * P**2 + path_b in u_labels
        and path_e in v_labels
    )
    require(path_active, "supplied ancestry path became inactive")
    return (
        len(u_labels),
        len(v_labels),
        ancestry.path_digest(u_labels, v_labels),
        chambers,
        u_labels,
        v_labels,
    )


def exceptional_chain_profile(pieces, live, half_step):
    live_support = {piece[:2] for piece in live}
    blocks = []
    gap_steps = []
    current = []
    for index, piece in enumerate(pieces):
        if index:
            step = piece[0] - pieces[index - 1][0]
            if step != half_step:
                blocks.append(tuple(current))
                current = []
                gap_steps.append(step)
        current.append(piece[:2] in live_support)
    blocks.append(tuple(current))
    require(
        all(
            value == (index % 2 == 0)
            for block in blocks
            for index, value in enumerate(block)
        ),
        "a mixed block stopped alternating 1,0 from live",
    )
    profiles = tuple(
        (len(block), sum(block), len(block) - sum(block))
        for block in blocks
    )
    require(
        all(live_count == (length + 1) // 2
            for length, live_count, _dead_count in profiles),
        "block live count stopped being a ceiling half",
    )
    return profiles, tuple(gap_steps)


def valuation(value, prime):
    answer = 0
    while value and value % prime == 0:
        answer += 1
        value //= prime
    return answer


def main():
    require(
        T % (2 * DEPTH) == 0,
        "half-depth step left the integer grid",
    )
    half_step = T // (2 * DEPTH)
    (
        module,
        _rails,
        present,
        details,
        full_module,
        e3,
        clocks,
        q_pairs,
        delayed,
        target_weight,
        rail_common,
    ) = physical_setup()
    ancestry_arguments, u_events, v_events = ancestry_setup()

    cells = []
    semantic_contents = set()
    ancestry_cache = {}
    exceptional_profiles = {}
    exceptional_gaps = {}
    for clock in range(len(clocks)):
        for target in range(P):
            pieces = physical_right(
                details,
                full_module,
                e3,
                clocks,
                clock,
                target,
            )
            if not pieces:
                continue
            require(
                all(weight > 0 for _left, _right, weight in pieces),
                f"nonpositive physical piece at {(clock, target)}",
            )
            require(
                len({right - left for left, right, _weight in pieces})
                == 1,
                f"unequal piece lengths at {(clock, target)}",
            )
            require(
                len({weight for _left, _right, weight in pieces}) == 1,
                f"unequal piece weights at {(clock, target)}",
            )
            live, dead = semantic_partition(
                pieces, q_pairs[clock]
            )
            require(live, f"raw-only cell at {(clock, target)}")
            for left, right, _weight in live:
                unit = ((left, right, 1),)
                semantic_contents.add(
                    physical.relative.private.delayed_carry_pair(
                        unit, q_pairs[clock], {}
                    )[12][1]
                )

            ancestry_data = ancestry_record(
                pieces,
                ancestry_arguments,
                u_events,
                v_events,
                ancestry_cache,
            )
            if dead:
                profile, gaps = exceptional_chain_profile(
                    pieces, live, half_step
                )
                exceptional_profiles[clock, target] = profile
                exceptional_gaps[clock, target] = gaps

            raw_support = raw_right(
                module, rail_common, clocks, clock, target
            )
            total_vector, _total_masses = wing.coefficient_vector(
                module,
                delayed,
                present,
                target_weight,
                raw_support,
                12,
                [{} for _ in range(7)],
            )
            copy_vector, _copy_masses = wing.coefficient_vector(
                module,
                delayed,
                present,
                target_weight,
                (live[0][:2],),
                12,
                [{} for _ in range(7)],
            )
            require(
                copy_vector[0] == 0
                and len(set(copy_vector[1:])) == 1
                and copy_vector[1] > 0,
                f"live copy vector lost its scalar form at {(clock, target)}",
            )
            require(
                total_vector
                == (0,) + (len(live) * copy_vector[1],) * 6,
                f"raw coefficient is not the live-copy sum at "
                f"{(clock, target)}",
            )
            if dead:
                dead_vector, _dead_masses = wing.coefficient_vector(
                    module,
                    delayed,
                    present,
                    target_weight,
                    (dead[0][:2],),
                    12,
                    [{} for _ in range(7)],
                )
                require(
                    dead_vector == (0,) * 7,
                    f"dead copy contributed at {(clock, target)}",
                )

            weight = pieces[0][2]
            copy_value = copy_vector[1]
            cells.append(
                {
                    "cell": (clock, target),
                    "raw": len(pieces),
                    "live": len(live),
                    "dead": len(dead),
                    "length": pieces[0][1] - pieces[0][0],
                    "weight": weight,
                    "U": ancestry_data[0],
                    "V": ancestry_data[1],
                    "digest": ancestry_data[2],
                    "chambers": ancestry_data[3],
                    "u_labels": ancestry_data[4],
                    "v_labels": ancestry_data[5],
                    "copy_value": copy_value,
                    "total_value": total_vector[1],
                }
            )

    require(len(cells) == 28, "discovered nonzero-cell count is not 28")
    require(
        len(semantic_contents) == 1 and 0 not in semantic_contents,
        "live semantic content stopped being globally uniform",
    )
    semantic_content = next(iter(semantic_contents))
    require(
        all(
            row["copy_value"] == row["weight"] * semantic_content
            for row in cells
        ),
        "copy coefficient stopped factoring as weight times content",
    )
    require(
        len({row["length"] for row in cells}) == 1,
        "copy length stopped being global",
    )

    copy_values = tuple(sorted({row["copy_value"] for row in cells}))
    common_gcd = 0
    for value in copy_values:
        common_gcd = gcd(common_gcd, value)
    require(common_gcd > 0, "copy-value gcd vanished")
    primitive_copy_factors = tuple(
        value // common_gcd for value in copy_values
    )
    require(
        gcd(*primitive_copy_factors) == 1,
        "derived copy factors retain a common divisor",
    )
    require(
        valuation(common_gcd, P) == 1
        and all(valuation(value, P) == 1 for value in copy_values),
        "copy coefficients lost their simple 13-divisibility",
    )

    for row in cells:
        row["primitive_copy_factor"] = (
            row["copy_value"] // common_gcd
        )
        row["copy_beta"] = (row["copy_value"] // P) % P
        row["total_beta"] = (row["total_value"] // P) % P
        require(
            row["total_beta"]
            == row["live"] * row["copy_beta"] % P,
            f"Bockstein failed copywise additivity at {row['cell']}",
        )
        require(
            row["total_value"] // row["copy_value"] == row["live"],
            f"primitive multiplier is not the live count at {row['cell']}",
        )
        require(
            row["weight"] == row["U"] * row["V"],
            f"ancestry cardinality product missed weight at {row['cell']}",
        )

    cells.sort(key=lambda row: row["cell"])
    baseline = cells[0]
    first_global_change = next(
        row for row in cells if row["digest"] != baseline["digest"]
    )
    baseline_u = set(baseline["u_labels"])
    changed_u = set(first_global_change["u_labels"])
    baseline_v = set(baseline["v_labels"])
    changed_v = set(first_global_change["v_labels"])
    prototype_change = {
        "cell": first_global_change["cell"],
        "U_lost": len(baseline_u - changed_u),
        "U_gained": len(changed_u - baseline_u),
        "V_lost": len(baseline_v - changed_v),
        "V_gained": len(changed_v - baseline_v),
        "path_active_before": True,
        "path_active_after": True,
    }
    require(
        baseline_v == changed_v,
        "first prototype change unexpectedly altered V",
    )

    cell_table = tuple(
        (
            row["cell"],
            row["raw"],
            row["live"],
            row["dead"],
            row["length"],
            row["weight"],
            row["U"],
            row["V"],
            row["digest"],
            row["copy_value"],
            row["primitive_copy_factor"],
            row["copy_beta"],
            row["total_beta"],
        )
        for row in cells
    )
    partition_histogram = Counter(
        (row["raw"], row["live"], row["dead"]) for row in cells
    )
    prototype_groups = {}
    for row in cells:
        prototype_groups.setdefault(row["digest"], []).append(row["cell"])
    prototype_groups = tuple(
        (digest, tuple(group))
        for digest, group in prototype_groups.items()
    )
    exceptional_gap_ratios = {
        cell: tuple(step // half_step for step in steps)
        for cell, steps in exceptional_gaps.items()
    }

    print("THM-2771 ALL-CELL HOSTILE RECONSTRUCTION")
    print(
        "status=FINITE-EXACT independent audit; imports=pinned canonical modules only;"
        "earlier cofiber scratch not imported"
    )
    print(f"dependency_hashes={PINNED}")
    print(
        "columns=(cell,raw,live,dead,length,weight,U,V,digest,"
        "copy_value,primitive_copy_factor,copy_beta,total_beta)"
    )
    for row in cell_table:
        print(row)
    print(f"nonzero_cell_count={len(cell_table)}")
    print(f"partition_histogram={dict(partition_histogram)}")
    print(f"global_copy_length={cells[0]['length']}")
    print(f"semantic_content={semantic_content}")
    print(f"half_step=T/(2*13^5)={half_step}")
    print(f"exceptional_chain_profiles={exceptional_profiles}")
    print(f"exceptional_interblock_steps={exceptional_gaps}")
    print(
        "exceptional_interblock_half_step_ratios="
        f"{exceptional_gap_ratios}"
    )
    print(
        "exceptional_formula="
        "within each discovered block live positions are exactly even "
        "indices, hence live=sum((block_length+1)//2)"
    )
    print(f"ancestry_wall_free_chamber_count={len(ancestry_cache)}")
    print(f"prototype_groups={prototype_groups}")
    print(f"first_global_prototype_change={prototype_change}")
    print(
        "ancestry_scope=identity is proved within each raw cell by absence "
        "of U/V walls across its full hull; it is not a global prototype "
        "claim"
    )
    print(f"derived_copy_values={copy_values}")
    print(f"derived_common_gcd={common_gcd}")
    print(f"derived_primitive_copy_factors={primitive_copy_factors}")
    print(
        "derived_primitive_multiplier_values="
        f"{tuple(sorted({row['live'] for row in cells}))}"
    )
    print(
        "bockstein_audit=every copy value and their gcd have v13=1;"
        "total_beta=(live_count*copy_beta) mod 13 in 28/28 cells"
    )
    print(
        "verdict=the 28-cell table, exceptional alternating half-step "
        "chains, local ancestry constancy, global U-prototype change, and "
        "copywise Bockstein multipliers are independently reproduced"
    )


if __name__ == "__main__":
    main()
