#!/usr/bin/env python3
"""Exact endpoint-arrival guard-danger and D-return bypass scout.

This is a standalone, noncanonical scout downstream of THM-2584/2586 and
THM-2623.  It asks why the endpoint arrival digits ``v=0,12`` disappeared
from the guard-safe carrier that eventually became THM-2616, and tests the
lawful complementary guard-danger leg on both the present and delayed words.

The computation first retains the physical THM-2584 rail label
``(s,ell,v,t)``, the absolute target section ``q=0``, the following shallow
clock ``ell5``, and the nonzero deep probe ``r``.  It then audits all 162
positive endpoint rails against all thirteen target sections, divides by one
global content, and computes every primitive seven-clock multiplication
determinant and owner Bockstein.  Positive integer profile weights are never
cancelled or reprimitivized slice by slice.

It also proves the elementary base-thirteen return-clock formula for every
arrival digit and for the four endpoint-to-endpoint digit pairs.  No
semantic old head, adjacent-event gluing, iterated physical D-chain, row
exclusion, or LRC(14) conclusion is tested.
"""

from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from math import gcd
import hashlib

import lrc14_base_only_bridge_opus_20260728 as base
import lrc14_b_r5_owner_clock_host_thm2581 as host
import lrc14_b_r5_theta_target_tensor_thm2584 as tensor
import lrc14_cross_time_target_future_diagonal_thm2616 as cross
import lrc14_fallback_rail_digit_diagonal_pullback_thm2592 as old
import lrc14_guard_cospan_successor_private_clock_collapse as guard


P = 13
Q7 = 7
T = old.T
ENDPOINTS = ((0, 0), (12, 12))
EXPECTED_ENDPOINT_ZEROS = {
    (0, 0): ((7, 4), (7, 5), (7, 6)),
    (12, 12): ((6, 1), (6, 2), (6, 3)),
}
EXPECTED_FIRST_NUMERATORS = {
    0: 14_955_985_505_560_745_933_829_360,
    12: 637_612_774_685_803_357_216_776,
}
PRIMITIVE_SHARDS = ((0, 41), (41, 81), (81, 122), (122, 162))
EXPECTED_PRIMITIVE_DIGEST = (
    "5d90d8c4b1f9a6537ce3383537a698c0383700646aaa9643f4e250f524f94a9f"
)
EXPECTED_BOCKSTEIN_DIGEST = (
    "5945bf59b4d90ca9098b78b2582d74fe281b2a2385bfad54129205308fbed370"
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def merge_adjacent(intervals):
    """Normalize a sorted half-open union, merging overlaps and contacts."""
    out = []
    for left, right in sorted(intervals):
        if out and left <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], right))
        else:
            out.append((left, right))
    return out


def interval_mass(intervals):
    return sum(right - left for left, right in intervals)


def build_present_omit_guard(module, ell, shift):
    """The THM-2592 present packet with only its guard factor omitted."""
    intervals = module.make_comb(
        module.C1, 182, 26 * ell - 13, 26 * ell + 13
    )
    intervals = module.subtract_comb(
        intervals, module.W[1], 182,
        -14 * shift - 13, -14 * shift + 13,
    )
    for index in (2, 3, 4, 5):
        intervals = module.subtract_comb(
            intervals, module.W[index], 182, -13, 13
        )
    intervals = module.subtract_comb(
        intervals, module.C2, 182,
        14 * shift - 13, 14 * shift + 13,
    )
    intervals = module.subtract_comb(
        intervals, module.C3, 182, -13, 13
    )
    return intervals


def build_present_danger(module, ell, shift):
    """Replace the present packet's guard-safe factor by guard danger."""
    intervals = module.make_comb(
        module.C1, 182, 26 * ell - 13, 26 * ell + 13
    )
    intervals = module.intersect_comb(
        intervals, module.W[module.GUARD], 91, -13, 13
    )
    intervals = module.subtract_comb(
        intervals, module.W[1], 182,
        -14 * shift - 13, -14 * shift + 13,
    )
    for index in (2, 3, 4, 5):
        intervals = module.subtract_comb(
            intervals, module.W[index], 182, -13, 13
        )
    intervals = module.subtract_comb(
        intervals, module.C2, 182,
        14 * shift - 13, 14 * shift + 13,
    )
    intervals = module.subtract_comb(
        intervals, module.C3, 182, -13, 13
    )
    return intervals


def build_endpoint_rails():
    """Rebuild all positive THM-2584 endpoint rails with literal pieces."""
    e_set = base.build_set(base.PAT_E, base.ZELL)
    q_b = base.build_set(host.PAT_QB, base.ZELL)
    ust, uv, vst, vv = tensor.packet_profiles(e_set, q_b)
    _, _, k_tensor = tensor.exact_tensor(e_set, q_b)
    owner = base.clock_cells(P, Q7, T, P * P)
    deep = tensor.deep_cells()
    arrival = [
        [(v * (T // P), (v + 1) * (T // P))]
        for v in range(P)
    ]

    rails = {}
    zero_sets = {}
    for v, t in ENDPOINTS:
        bank = []
        zeros = []
        for s in range(1, P):
            rotated_starts, rotated_values = base.rotate_profile(
                vst, vv, s * (T // P), T
            )
            profile_starts, profile_values, _, _ = base.product_cum(
                ust, uv, rotated_starts, rotated_values, T
            )
            for ell in range(Q7):
                mass = k_tensor[s][v][t][ell]
                if not mass:
                    zeros.append((s, ell))
                    continue
                cell = old.intersect_sorted(
                    old.intersect_sorted(owner[ell], deep[t]),
                    arrival[v],
                )
                pieces = old.profile_on_intervals(
                    cell, profile_starts, profile_values
                )
                rebuilt_mass = P * sum(
                    weight * (right - left)
                    for left, right, weight in pieces
                )
                require(rebuilt_mass == mass > 0,
                        "endpoint rail mass reconstruction failed")
                bank.append({
                    "s": s,
                    "ell": ell,
                    "v": v,
                    "t": t,
                    "pieces": pieces,
                    "mass": mass,
                })
        require(tuple(zeros) == EXPECTED_ENDPOINT_ZEROS[v, t],
                "endpoint rail zero set changed")
        require(len(bank) == 81,
                "endpoint rail bank is no longer 81 of 84 cells")
        rails[v] = bank
        zero_sets[v] = tuple(zeros)
    return rails, zero_sets


def clock_label(x):
    """Half-open seven-clock label, centred at ell/7."""
    return int((7 * x + Fraction(1, 2)) // 1) % 7


def return_clock_profile(first_digit, second_digit):
    """Masses of (clock(z),clock(Dz)) on z in I_v, Dz in I_w."""
    left = Fraction(second_digit, P)
    right = Fraction(second_digit + 1, P)
    cuts = {left, right}
    boundaries = tuple(Fraction(2 * ell + 1, 14) for ell in range(Q7))
    for boundary in boundaries:
        if left < boundary < right:
            cuts.add(boundary)
        pulled = P * boundary - first_digit
        if left < pulled < right:
            cuts.add(pulled)
    cuts = sorted(cuts)
    masses = Counter()
    for a, b in zip(cuts, cuts[1:]):
        middle = (a + b) / 2
        z = (first_digit + middle) / P
        masses[clock_label(z), clock_label(middle)] += (b - a) / P
    require(sum(masses.values(), Fraction(0)) == Fraction(1, P * P),
            "return-clock profile has wrong total mass")
    return dict(sorted(masses.items()))


def pure_expected_profile(v):
    if v % 2 == 1:
        label = (v + 1) // 2
        return {(label, label): Fraction(1, 169)}
    if v == 6:
        return {
            (3, 3): Fraction(1, 338),
            (4, 4): Fraction(1, 338),
        }
    if v < 6:
        label = v // 2
        return {
            (label, label): Fraction(13 - v, 2366),
            (label, label + 1): Fraction(v + 1, 2366),
        }
    label = (v // 2 + 1) % Q7
    return {
        (label, v // 2): Fraction(13 - v, 2366),
        (label, label): Fraction(v + 1, 2366),
    }


def format_profile(profile):
    return tuple((pair, str(mass)) for pair, mass in profile.items())


def build_primitive_data():
    """Rebuild the complete endpoint-danger coefficient input."""
    (module, inherited_safe, _, _, _, _, _) = cross.build_carrier_data()
    delayed_prefixes, _, _, _ = guard.build_guard_cospan(
        module, inherited_safe
    )
    rails, _ = build_endpoint_rails()
    rail_bank = tuple(rail for v in (0, 12) for rail in rails[v])
    require(len(rail_bank) == 162,
            "endpoint primitive rail-bank census changed")

    present = {}
    starts = {}
    for q in range(P):
        shift = (-q) % P
        for ell in range(Q7):
            intervals = build_present_danger(module, ell, shift)
            present[q, ell] = intervals
            starts[q, ell] = [left for left, _ in intervals]
    return module, delayed_prefixes["danger"], rail_bank, present, starts


def primitive_shard(bounds):
    """Compute raw numerators on one disjoint rail shard."""
    start, stop = bounds
    module, delayed, rail_bank, present, starts = build_primitive_data()
    content = 0
    positive = 0
    controls = 0
    metadata = []
    rows = []
    for j in range(start, stop):
        rail = rail_bank[j]
        v = rail["v"]
        metadata.append((rail["s"], rail["ell"], v, rail["t"]))
        raw = [[[0] * P for _ in range(Q7)] for _ in range(P)]
        for q in range(P):
            for ell5 in range(Q7):
                overlap = old.intersect_weighted_union(
                    rail["pieces"], present[q, ell5], starts[q, ell5]
                )
                if not overlap:
                    continue
                prefix = delayed[ell5][v]
                for r in range(1, P):
                    probed = old.intersect_weighted_comb(
                        overlap, module.C3, 182,
                        14 * r - 13, 14 * r + 13,
                    )
                    value = old.delayed_weighted_numerator(probed, prefix)
                    if value and controls < 2:
                        require(
                            value == old.grouped_original_numerator(
                                module, probed, prefix
                            ),
                            "weighted and grouped delayed routes disagree",
                        )
                        controls += 1
                    raw[q][ell5][r] = value
                    if value:
                        positive += 1
                        content = gcd(content, value)
        rows.append(tuple(
            tuple(tuple(root_row) for root_row in clock_rows)
            for clock_rows in raw
        ))
    return bounds, content, positive, controls, tuple(metadata), tuple(rows)


def primitive_unit_data(raw, content):
    """Return the seven coefficients, determinant, and owner evaluations."""
    y = tuple(
        sum(
            (raw[ell][r] // content) * pow(r, -1, P)
            for r in range(1, P)
        ) % P
        for ell in range(Q7)
    )
    reduced = tuple((y[ell] - y[-1]) % P for ell in range(Q7 - 1))
    determinant = old.sat.multiplication_determinant_7(reduced)
    owner_nonzero = tuple(
        any(old.sat.owner_factor(y, kappa)) for kappa in range(1, Q7)
    )
    return y, determinant, owner_nonzero


def run_primitive_atlas():
    """Audit all 162 endpoint rails and all thirteen target sections."""
    with ProcessPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(primitive_shard, PRIMITIVE_SHARDS))

    content = 0
    positive = 0
    controls = 0
    metadata = []
    joint = []
    for bounds, local_content, local_positive, local_controls, meta, rows in results:
        require(bounds in PRIMITIVE_SHARDS,
                "worker returned an unknown primitive shard")
        content = gcd(content, local_content)
        positive += local_positive
        controls += local_controls
        metadata.extend(meta)
        joint.extend(rows)
    require(len(joint) == len(metadata) == 162,
            "joined primitive endpoint bank census changed")
    require(content > 0, "primitive endpoint bank is identically zero")

    valuation13 = 0
    temp = content
    while temp % P == 0:
        valuation13 += 1
        temp //= P

    pair_supports = []
    pair_totals = []
    determinant_hist = Counter()
    unit_flags = [[False] * P for _ in joint]
    beta_flags = [[False] * P for _ in joint]
    y_rows = []
    zero_y_slices = []
    zero_beta_slices = []
    nonunit_slices = []
    for j, rows in enumerate(joint):
        for q in range(P):
            support = sum(
                bool(rows[q][ell][r])
                for ell in range(Q7) for r in range(1, P)
            )
            total = sum(
                rows[q][ell][r]
                for ell in range(Q7) for r in range(1, P)
            )
            pair_supports.append(support)
            pair_totals.append(total)
            y, determinant, owner_nonzero = primitive_unit_data(
                rows[q], content
            )
            y_rows.append(y)
            determinant_hist[determinant] += 1
            unit_flags[j][q] = determinant != 0
            beta_flags[j][q] = all(owner_nonzero)
            if not any(y):
                zero_y_slices.append((metadata[j], q))
            if not all(owner_nonzero):
                zero_beta_slices.append((metadata[j], q, y, owner_nonzero))
            if determinant == 0:
                nonunit_slices.append((metadata[j], q, y))

    primitive_serialization = ";".join(
        ",".join(
            str(joint[j][q][ell][r] // content)
            for ell in range(Q7) for r in range(P)
        )
        for j in range(len(joint)) for q in range(P)
    )
    primitive_digest = hashlib.sha256(
        primitive_serialization.encode("ascii")
    ).hexdigest()
    y_digest = hashlib.sha256(str(y_rows).encode("ascii")).hexdigest()

    by_cell = {}
    for j, (s, ell, v, _) in enumerate(metadata):
        by_cell.setdefault((s, ell), []).append((j, v))
    require(len(by_cell) == 84, "primitive endpoint bank lost a base cell")

    fixed_q_positive_cells = []
    fixed_q_beta_cells = []
    fixed_q_unit_cells = []
    for q in range(P):
        fixed_q_positive_cells.append(sum(
            any(pair_totals[j * P + q] > 0 for j, _ in edges)
            for edges in by_cell.values()
        ))
        fixed_q_beta_cells.append(sum(
            any(beta_flags[j][q] for j, _ in edges)
            for edges in by_cell.values()
        ))
        fixed_q_unit_cells.append(sum(
            any(unit_flags[j][q] for j, _ in edges)
            for edges in by_cell.values()
        ))

    full_q = tuple(q for q, count in enumerate(fixed_q_unit_cells)
                   if count == 84)
    selected_contents = {}
    selected_patterns = {}
    for q in full_q:
        chosen = []
        for cell, edges in sorted(by_cell.items()):
            candidates = [(j, v) for j, v in edges if unit_flags[j][q]]
            require(candidates, "declared full target section lost a unit")
            j, v = min(candidates, key=lambda pair: pair[1])
            chosen.append((cell, j, v))
        selected_content = 0
        for _, j, _ in chosen:
            for ell in range(Q7):
                for r in range(1, P):
                    selected_content = gcd(
                        selected_content, joint[j][q][ell][r]
                    )
        selected_contents[q] = selected_content
        selected_patterns[q] = tuple(
            cell for cell, _, v in chosen if v == 12
        )

    fixed_q_unit_missing = tuple(
        (
            q,
            tuple(
                cell for cell, edges in sorted(by_cell.items())
                if not any(unit_flags[j][q] for j, _ in edges)
            ),
        )
        for q in range(P) if fixed_q_unit_cells[q] < 84
    )
    cell_unit_choice_histogram = Counter(
        sum(unit_flags[j][q]
            for j, _ in edges for q in range(P))
        for edges in by_cell.values()
    )
    q0_expected_v12 = tuple(
        (s, ell) for s in (5, 7) for ell in range(Q7)
    )

    require(positive == 68_904,
            "primitive positive fine-entry census changed")
    require(sum(total > 0 for total in pair_totals) == 2_002,
            "primitive positive rail/section census changed")
    require(Counter(pair_supports) == Counter({
        0: 104, 12: 352, 24: 330, 36: 550, 48: 770,
    }), "primitive rail/section support histogram changed")
    require(content == 86_814 and valuation13 == 1,
            "endpoint danger global content changed")
    require(controls == 8,
            "independent grouped-route control census changed")
    require(sum(map(sum, unit_flags)) == 1_770,
            "primitive unit-slice census changed")
    require(sum(map(sum, beta_flags)) == 1_786,
            "primitive owner-Bockstein census changed")
    require(len(zero_y_slices) == len(zero_beta_slices) == 320
            and len(nonunit_slices) == 336,
            "primitive zero/nonunit census changed")
    require(determinant_hist == Counter({
        0: 336, 1: 360, 2: 68, 3: 202, 4: 128, 5: 122, 6: 82,
        7: 54, 8: 54, 9: 76, 10: 150, 11: 40, 12: 434,
    }), "primitive determinant histogram changed")
    require(fixed_q_positive_cells
            == [84, 84, 83, 84, 84, 84, 84, 84, 84, 84, 84, 83, 84],
            "fixed-section positive coverage changed")
    require(fixed_q_beta_cells
            == [84, 84, 82, 82, 84, 82, 84, 84, 82, 84, 82, 82, 84],
            "fixed-section owner-Bockstein coverage changed")
    require(fixed_q_unit_cells
            == [84, 84, 82, 82, 84, 81, 84, 84, 81, 84, 82, 82, 84],
            "fixed-section unit coverage changed")
    require(full_q == (0, 1, 4, 6, 7, 9, 12),
            "uniform unit target-section set changed")
    require(all(selected_contents[q] == content for q in full_q),
            "a uniform selector acquired extra primitive content")
    require(selected_patterns[0] == q0_expected_v12,
            "clean q=0 endpoint selector changed")
    require(primitive_digest == EXPECTED_PRIMITIVE_DIGEST,
            "primitive tensor digest changed")
    require(y_digest == EXPECTED_BOCKSTEIN_DIGEST,
            "primitive Bockstein-row digest changed")

    print("primitive_atlas=all 162 positive endpoint rails x 13 target sections")
    print(
        f"primitive_positive_fine={positive}/176904 "
        f"positive_pairs={sum(total > 0 for total in pair_totals)}/2106"
    )
    print(
        "primitive_pair_support_histogram="
        + str(tuple(sorted(Counter(pair_supports).items())))
    )
    print(f"primitive_global_content={content} v13={valuation13}")
    print(
        f"primitive_units={sum(map(sum, unit_flags))}/2106 "
        f"good_owner_bocksteins={sum(map(sum, beta_flags))}/2106"
    )
    print(f"primitive_tensor_digest={primitive_digest}")
    print(f"primitive_bockstein_row_digest={y_digest}")
    print(f"fixed_q_positive_cells={tuple(fixed_q_positive_cells)}")
    print(f"fixed_q_unit_cells={tuple(fixed_q_unit_cells)}")
    print(f"fixed_q_unit_missing={fixed_q_unit_missing}")
    print(f"uniform_unit_target_sections={full_q}")
    print(
        "uniform_selector_contents="
        + str(tuple((q, selected_contents[q]) for q in full_q))
    )
    print(
        "uniform_selector_v12_counts="
        + str(tuple((q, len(selected_patterns[q])) for q in full_q))
    )
    print(
        "q0_unit_selector=v12 iff source s in {5,7}, otherwise v0; "
        "selected_units=84/84; selected_content=86814"
    )
    print(
        "cell_unit_choice_histogram="
        + str(tuple(sorted(cell_unit_choice_histogram.items())))
    )


def main():
    print("LRC14 endpoint-arrival guard-danger D-return bypass scout")
    print("status=VERIFIED-EXACT SCOUT; not a theorem dependency")
    print(
        "scope=THM-2584 endpoint rails; present and delayed guard-danger "
        "cospan legs; exact q=0 support and globally primitive all-q atlas"
    )
    print(
        "not_tested=semantic old head; adjacent-event gluing; iterated "
        "physical D-chain; row consequence"
    )

    (module, inherited_safe, _, _, _, _, _) = cross.build_carrier_data()
    delayed_prefixes, _, _, _ = guard.build_guard_cospan(
        module, inherited_safe
    )
    rails, zero_sets = build_endpoint_rails()

    guard_danger = module.make_comb(
        module.W[module.GUARD], 91, -13, 13
    )
    require(guard_danger == [
        (0, T // 7),
        (6 * T // 7, T),
    ], "present guard-danger interval changed")
    endpoint_windows = {
        0: (0, T // 26),
        12: (25 * T // 26, T),
    }
    for v in (0, 12):
        lo, hi = endpoint_windows[v]
        require(all(
            lo <= left < right <= hi
            for rail in rails[v]
            for left, right, _ in rail["pieces"]
        ), "endpoint rail escaped its arrival/deep interval")
        require(all(
            any(a <= left and right <= b for a, b in guard_danger)
            for rail in rails[v]
            for left, right, _ in rail["pieces"]
        ), "endpoint rail is not structurally guard-danger")

    present_safe = {}
    present_danger = {}
    shallow = {}
    present_partition_checks = 0
    for ell in range(Q7):
        shallow[ell] = module.make_comb(
            module.C1, 182, 26 * ell - 13, 26 * ell + 13
        )
        for shift in range(P):
            safe = module.build_F(ell, shift)
            danger = build_present_danger(module, ell, shift)
            omitted = build_present_omit_guard(module, ell, shift)
            require(not old.intersect_sorted(safe, danger),
                    "present guard cospan legs overlap")
            require(interval_mass(safe) + interval_mass(danger)
                    == interval_mass(omitted),
                    "present guard cospan does not partition omitted word")
            require(merge_adjacent(safe + danger) == omitted,
                    "present guard cospan union changed")
            present_safe[ell, shift] = safe
            present_danger[ell, shift] = danger
            present_partition_checks += 1
    require(present_partition_checks == 91,
            "present guard partition-check universe changed")

    # Locate the exact first killed factor on the endpoint rails.  The
    # shallow C1 clock still meets them, but guard-safety kills every pair.
    safe_gate_counts = {}
    for v in (0, 12):
        shallow_positive = 0
        safe_positive = 0
        for rail in rails[v]:
            for ell in range(Q7):
                if old.intersect_weighted_union(
                    rail["pieces"], shallow[ell],
                    [left for left, _ in shallow[ell]],
                ):
                    shallow_positive += P
                for shift in range(P):
                    safe = present_safe[ell, shift]
                    if old.intersect_weighted_union(
                        rail["pieces"], safe,
                        [left for left, _ in safe],
                    ):
                        safe_positive += 1
        require(shallow_positive == 4_056,
                "endpoint/shallow positive census changed")
        require(safe_positive == 0,
                "an endpoint rail survived present guard-safety")
        safe_gate_counts[v] = (shallow_positive, safe_positive)

    danger_starts = {
        key: [left for left, _ in intervals]
        for key, intervals in present_danger.items()
    }

    # Fixed q=0.  Retain each endpoint option separately, then select v=12
    # exactly on source s=7 and v=0 otherwise.
    options = {}
    first_numerators = {}
    for v in (0, 12):
        for rail in rails[v]:
            positive_entries = []
            for ell5 in range(Q7):
                present = old.intersect_weighted_union(
                    rail["pieces"], present_danger[ell5, 0],
                    danger_starts[ell5, 0],
                )
                if not present:
                    continue
                for r in range(1, P):
                    probed = old.intersect_weighted_comb(
                        present, module.C3, 182,
                        14 * r - 13, 14 * r + 13,
                    )
                    if not probed:
                        continue
                    value = old.delayed_weighted_numerator(
                        probed, delayed_prefixes["danger"][ell5][v]
                    )
                    if value:
                        positive_entries.append((ell5, r, value))
                        first_numerators.setdefault(
                            v,
                            (rail["s"], rail["ell"], 0, ell5, r, value),
                        )
            if positive_entries:
                options.setdefault((rail["s"], rail["ell"]), {})[v] = (
                    tuple(positive_entries)
                )

    option_histogram = Counter(len(options.get((s, ell), {}))
                               for s in range(1, P)
                               for ell in range(Q7))
    require(option_histogram == Counter({2: 70, 1: 14}),
            "endpoint danger-option histogram changed")
    require(all(options.get((s, ell)) for s in range(1, P)
                for ell in range(Q7)),
            "endpoint danger atlas lost a base cell")
    require(all(set(options[s, ell]) == {0} for s in (6,)
                for ell in range(Q7)),
            "source six endpoint singleton changed")
    require(all(set(options[s, ell]) == {12} for s in (7,)
                for ell in range(Q7)),
            "source seven endpoint singleton changed")

    selected = []
    selected_support_histogram = Counter()
    selected_total_numerators = []
    for s in range(1, P):
        for ell in range(Q7):
            v = 12 if s == 7 else 0
            require(v in options[s, ell],
                    "uniform endpoint source selector is not positive")
            entries = options[s, ell][v]
            selected.append((s, ell, v, len(entries)))
            selected_support_histogram[len(entries)] += 1
            selected_total_numerators.append(sum(value for _, _, value in entries))
    require(Counter(v for _, _, v, _ in selected)
            == Counter({0: 77, 12: 7}),
            "uniform endpoint selector census changed")
    require(first_numerators[0][-1] == EXPECTED_FIRST_NUMERATORS[0]
            and first_numerators[12][-1] == EXPECTED_FIRST_NUMERATORS[12],
            "first endpoint delayed numerator changed")

    print(
        "endpoint_zero_sets="
        + str(tuple((v, zero_sets[v]) for v in (0, 12)))
    )
    print(
        "endpoint_support_intervals="
        "v0:[0,1/26) v12:[25/26,1); both lie in present guard danger"
    )
    print(
        "present_guard_partition_checks=91 "
        f"safe_gate_counts={tuple((v, safe_gate_counts[v]) for v in (0,12))}"
    )
    print(
        "safe_gate_columns=(arrival_v,shallow_positive_pairs,"
        "guard_safe_positive_pairs)"
    )
    print(
        "danger_q0_option_histogram=" + str(tuple(sorted(option_histogram.items())))
    )
    print(
        "danger_q0_support_selector=v12 iff source s=7 else v0; "
        "selected_counts=(v0:77,v12:7)"
    )
    print(
        "danger_q0_selected_fine_support_histogram="
        + str(tuple(sorted(selected_support_histogram.items())))
    )
    print(
        "danger_q0_selected_min_total_raw_numerator="
        + str(min(selected_total_numerators))
    )
    print(
        "first_raw_numerators="
        + str(tuple((v, first_numerators[v]) for v in (0, 12)))
    )

    run_primitive_atlas()

    diagonal_digits = []
    bypass_digits = []
    for v in range(P):
        profile = return_clock_profile(v, v)
        expected = pure_expected_profile(v)
        require(profile == expected,
                f"pure arrival-{v} D-return profile changed")
        cylinder = (Fraction(14 * v, 169),
                    Fraction(14 * v + 1, 169))
        diagonal = all(left == right for left, right in profile)
        (diagonal_digits if diagonal else bypass_digits).append(v)
        print(
            f"pure_return_v={v} cylinder=({cylinder[0]},{cylinder[1]}) "
            f"clock_profile={format_profile(profile)}"
        )
    require(tuple(diagonal_digits) == (1, 3, 5, 6, 7, 9, 11),
            "diagonal-trap arrival digits changed")
    require(tuple(bypass_digits) == (0, 2, 4, 8, 10, 12),
            "off-diagonal bypass arrival digits changed")

    mixed_expected = {
        (0, 0): pure_expected_profile(0),
        (0, 12): {
            (0, 6): Fraction(1, 2366),
            (1, 0): Fraction(1, 182),
        },
        (12, 0): {
            (6, 0): Fraction(1, 182),
            (0, 1): Fraction(1, 2366),
        },
        (12, 12): pure_expected_profile(12),
    }
    for pair in ((0, 0), (0, 12), (12, 0), (12, 12)):
        profile = return_clock_profile(*pair)
        require(profile == mixed_expected[pair],
                "mixed endpoint return profile changed")
        off_diagonal = sum(
            mass for (left, right), mass in profile.items()
            if left != right
        )
        print(
            f"endpoint_return={pair} clock_profile={format_profile(profile)} "
            f"off_diagonal_mass={off_diagonal}"
        )
    print(
        "return_formula=C_(v,w)=[(13v+w)/169,(13v+w+1)/169); "
        "pure diagonal iff v odd or v=6"
    )
    print(
        "verdict=PASS: endpoint arrivals have off-diagonal return-clock mass "
        "and admit a globally primitive q0 all-84 unit carrier"
    )
    print(
        "boundary=coefficient common-x atlas only; no semantic endpoint, "
        "adjacent D-event gluing, iterated chain, row exclusion, or LRC14; "
        "no conflict with full rail-envelope D nilpotence"
    )


if __name__ == "__main__":
    main()
