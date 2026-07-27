#!/usr/bin/env python3
"""Exact companion for THM-2600: attach all cells through arrival digit six.

This reuses THM-2592's common-x arithmetic, but replaces the theta-zero
priority selector by the middle vertex of THM-2584's four-edge path.  The
two edges (v,t)=(6,12),(6,0) have disjoint zero sets, so they cover all 84
nonzero-displacement/owner-clock cells.  The THM-2585 delayed word has
positive digit-six mass in every owner phase.
"""

from collections import Counter
from fractions import Fraction
from math import gcd
import hashlib

import lrc14_fallback_rail_digit_diagonal_pullback_thm2592 as old


P = old.P
Q7 = old.Q7
T = old.T
R = old.R
base = old.base
host = old.host
rail = old.rail
sat = old.sat


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    module = old.load_present_module()
    require(module.W == base.W, "row mismatch")

    # Split the THM-2585 delayed word by its literal old/future digit.
    word = module.build_word_Ta()
    word_prefix = [[None] * P for _ in range(Q7)]
    for ell in range(Q7):
        q_ell = module.subtract_comb(
            word, module.C1, 182, 26 * ell - 13, 26 * ell + 13
        )
        for h in range(P):
            digit = sat.intersect_interval(
                q_ell, h * T // P, (h + 1) * T // P
            )
            word_prefix[ell][h] = module.make_prefix(digit)
    digit_lengths = [[prefix[2][-1] for prefix in row] for row in word_prefix]
    require(all(digit_lengths[ell][0] == 0 for ell in range(Q7)),
            "digit-zero hostile disappeared")
    require(all(digit_lengths[ell][12] == 0 for ell in range(Q7)),
            "digit-twelve hostile disappeared")
    require(all(digit_lengths[ell][6] > 0 for ell in range(Q7)),
            "digit-six support disappeared")
    old.validate_delayed_formula(module, word_prefix)

    # Rebuild the fine THM-2584 tensor and retain every positive middle edge.
    e4 = base.build_set(base.PAT_E, base.ZELL)
    qb = base.build_set(host.PAT_QB, base.ZELL)
    ust, uv, vst, vv = rail.packet_profiles(e4, qb)
    _, _, k_tensor = rail.exact_tensor(e4, qb)
    owner = base.clock_cells(P, Q7, T, P * P)
    deep = rail.deep_cells()
    arrival = [[(v * (T // P), (v + 1) * (T // P))] for v in range(P)]

    zero_sets = {}
    for edge in ((0, 0), (6, 0), (6, 12), (12, 12)):
        v, t = edge
        zero_sets[edge] = tuple(
            (s, ell) for s in range(1, P) for ell in range(Q7)
            if k_tensor[s][v][t][ell] == 0
        )
    require(set(zero_sets[(6, 0)]).isdisjoint(zero_sets[(6, 12)]),
            "middle-edge zero sets overlap")

    rail_bank = []
    selected_mass = []
    theta_hist = Counter()
    for s4 in range(1, P):
        rvst, rvv = base.rotate_profile(vst, vv, s4 * (T // P), T)
        ps, pv, _, _ = base.product_cum(ust, uv, rvst, rvv, T)
        for ell4 in range(Q7):
            choices = [(6, t) for t in (12, 0)
                       if k_tensor[s4][6][t][ell4] > 0]
            require(choices, "constant-six middle selector failed")
            for v, t in choices:
                theta = (t - 2 * v) % P
                require(theta in (0, 1), "selected edge left the two rails")
                theta_hist[theta] += 1
                cell = old.intersect_sorted(
                    old.intersect_sorted(owner[ell4], deep[t]), arrival[v]
                )
                pieces = old.profile_on_intervals(cell, ps, pv)
                numerator = P * sum(w * (b - a) for a, b, w in pieces)
                require(numerator == k_tensor[s4][v][t][ell4],
                        "rail piece mass mismatch")
                rail_bank.append((s4, ell4, v, t, pieces))
                selected_mass.append(numerator)
    require(len(rail_bank) == 162 and theta_hist == Counter({0: 81, 1: 81}),
            "complete middle-edge census changed")
    nrails = len(rail_bank)

    # Present-side target/source packets.
    f_cache = {}
    f_starts = {}
    for s5 in range(P):
        for ell5 in range(Q7):
            f = module.build_F(ell5, s5)
            f_cache[ell5, s5] = f
            f_starts[ell5, s5] = [a for a, _ in f]

    joint = [[[[0] * P for _ in range(Q7)] for _ in range(P)]
             for _ in rail_bank]
    content = 0
    positive_entries = 0
    grouped_zero_controls = 0
    grouped_positive_controls = 0
    comb_controls = 0
    for j, (_, _, v, _, rail_pieces) in enumerate(rail_bank):
        require(v == 6, "arrival digit drifted")
        for q in range(P):
            s5 = (-q) % P
            for ell5 in range(Q7):
                present = old.intersect_weighted_union(
                    rail_pieces, f_cache[ell5, s5], f_starts[ell5, s5]
                )
                if not present:
                    continue
                prefix = word_prefix[ell5][v]
                for r in range(1, P):
                    probed = old.intersect_weighted_comb(
                        present, module.C3, 182, 14 * r - 13, 14 * r + 13
                    )
                    value = old.delayed_weighted_numerator(probed, prefix)
                    if comb_controls < 3:
                        grouped = {}
                        for a, b, w in present:
                            grouped.setdefault(w, []).append((a, b))
                        expected_len = 0
                        for w, intervals in grouped.items():
                            expected = module.intersect_comb(
                                intervals, module.C3, 182,
                                14 * r - 13, 14 * r + 13,
                            )
                            expected_len += w * sat.interval_length(expected)
                        require(
                            expected_len == sum(w * (b - a)
                                                for a, b, w in probed),
                            "direct comb and complement routes disagree",
                        )
                        comb_controls += 1
                    need_grouped = (
                        probed
                        and ((value == 0 and grouped_zero_controls < 8)
                             or (value > 0 and grouped_positive_controls < 8))
                    )
                    if need_grouped:
                        require(
                            value == old.grouped_original_numerator(
                                module, probed, prefix
                            ),
                            "weighted and grouped delayed sweeps disagree",
                        )
                        if value:
                            grouped_positive_controls += 1
                        else:
                            grouped_zero_controls += 1
                    joint[j][q][ell5][r] = value
                    if value:
                        positive_entries += 1
                        content = gcd(content, value)
        if (j + 1) % 14 == 0 or j + 1 == nrails:
            print(f"joint rail progress {j + 1}/{nrails}", flush=True)

    require(content > 0, "joint bank is identically zero")
    require(comb_controls == 3 and grouped_positive_controls == 8,
            "independent route-control census changed")

    pair_totals = [
        sum(joint[j][q][ell][r]
            for ell in range(Q7) for r in range(1, P))
        for j in range(nrails) for q in range(P)
    ]
    zero_pairs = [
        (rail_bank[j][:4], q)
        for j in range(nrails) for q in range(P)
        if pair_totals[j * P + q] == 0
    ]
    pair_supports = [
        sum(joint[j][q][ell][r] != 0
            for ell in range(Q7) for r in range(1, P))
        for j in range(nrails) for q in range(P)
    ]

    # One global primitive reduction, then all first Bocksteins.
    valuation13 = 0
    temp = content
    while temp % P == 0:
        valuation13 += 1
        temp //= P
    beta_nonzero = 0
    unit_slices = 0
    determinant_hist = Counter()
    zero_beta_slices = []
    nonunit_slices = []
    y_rows = []
    good_beta_by_rail = [0] * nrails
    unit_by_rail = [0] * nrails
    good_beta_by_q = [0] * P
    unit_by_q = [0] * P
    is_beta_good = [[False] * P for _ in range(nrails)]
    is_unit = [[False] * P for _ in range(nrails)]
    for j in range(nrails):
        for q in range(P):
            y = []
            for ell in range(Q7):
                y.append(sum(
                    (joint[j][q][ell][r] // content) * pow(r, -1, P)
                    for r in range(1, P)
                ) % P)
            y = tuple(y)
            reduced = tuple((y[i] - y[-1]) % P for i in range(Q7 - 1))
            det = sat.multiplication_determinant_7(reduced)
            determinant_hist[det] += 1
            unit_slices += int(det != 0)
            owner_nonzero = []
            for kappa in range(1, Q7):
                nz = any(sat.owner_factor(y, kappa))
                owner_nonzero.append(nz)
                beta_nonzero += int(nz)
            if not all(owner_nonzero):
                zero_beta_slices.append((rail_bank[j][:4], q, y,
                                         tuple(owner_nonzero)))
            else:
                good_beta_by_rail[j] += 1
                good_beta_by_q[q] += 1
                is_beta_good[j][q] = True
            if det == 0:
                nonunit_slices.append((rail_bank[j][:4], q, y))
            else:
                unit_by_rail[j] += 1
                unit_by_q[q] += 1
                is_unit[j][q] = True
            y_rows.append(y)

    serialization = ";".join(
        ",".join(str(joint[j][q][ell][r] // content)
                 for ell in range(Q7) for r in range(P))
        for j in range(nrails) for q in range(P)
    )
    digest = hashlib.sha256(serialization.encode("ascii")).hexdigest()
    y_digest = hashlib.sha256(str(y_rows).encode("ascii")).hexdigest()

    actual_content = P * content
    actual_denominator = 169 * (P**5) ** 2 * R * T
    primitive_denominator = Fraction(actual_denominator, actual_content)

    by_cell = {}
    for j, (s, ell, _, t, _) in enumerate(rail_bank):
        by_cell.setdefault((s, ell), []).append((j, t))
    require(len(by_cell) == 84, "middle edges lost a base cell")
    cell_positive = {}
    cell_beta = {}
    cell_units = {}
    for cell, edges in by_cell.items():
        positive = []
        beta = []
        units = []
        for j, t in edges:
            for q in range(P):
                if pair_totals[j * P + q] > 0:
                    positive.append((t, q))
                if is_beta_good[j][q]:
                    beta.append((t, q))
                if is_unit[j][q]:
                    units.append((t, q))
        cell_positive[cell] = positive
        cell_beta[cell] = beta
        cell_units[cell] = units
    fixed_q_unit_cells = [
        sum(any(is_unit[j][q] for j, _ in edges) for edges in by_cell.values())
        for q in range(P)
    ]

    # Freeze one uniform target section q=0 and prefer the theta-zero edge
    # t=12 whenever its slice is a unit.  The selected 84-slice bank is then
    # independently reduced by its own single global content.
    q0_selector = {}
    selected_indices = []
    for cell, edges in sorted(by_cell.items()):
        candidates = [(j, t) for j, t in edges if is_unit[j][0]]
        require(candidates, "q=0 lost a unit middle edge")
        j, t = min(candidates, key=lambda pair: (pair[1] != 12, pair[1]))
        require(pair_totals[j * P] > 0,
                "q=0 unit selector lost positive physical mass")
        q0_selector[cell] = t
        selected_indices.append(j)
    selected_content = 0
    for j in selected_indices:
        for ell in range(Q7):
            for r in range(1, P):
                selected_content = gcd(selected_content, joint[j][0][ell][r])
    require(selected_content > 0, "selected q=0 bank is empty")
    selected_unit_count = 0
    selected_y_rows = []
    for j in selected_indices:
        y = tuple(sum(
            (joint[j][0][ell][r] // selected_content) * pow(r, -1, P)
            for r in range(1, P)
        ) % P for ell in range(Q7))
        selected_y_rows.append(y)
        reduced = tuple((y[i] - y[-1]) % P for i in range(Q7 - 1))
        selected_unit_count += int(sat.multiplication_determinant_7(reduced) != 0)
    selected_serialization = ";".join(
        ",".join(str(joint[j][0][ell][r] // selected_content)
                 for ell in range(Q7) for r in range(P))
        for j in selected_indices
    )
    selected_digest = hashlib.sha256(
        selected_serialization.encode("ascii")
    ).hexdigest()

    require(zero_sets[(6, 0)] == ((7, 1), (7, 2), (7, 3)),
            "theta-one middle-edge zero set changed")
    require(zero_sets[(6, 12)] == ((6, 4), (6, 5), (6, 6)),
            "theta-zero middle-edge zero set changed")
    require(len(pair_totals) - len(zero_pairs) == 2024 and len(zero_pairs) == 82,
            "positive pair census changed")
    require(positive_entries == 61248,
            "fine positive-support census changed")
    require(Counter(pair_supports) == Counter({0: 82, 12: 154, 24: 814,
                                               36: 902, 48: 154}),
            "pair-support histogram changed")
    require(content == 4244240 and valuation13 == 1,
            "global primitive content changed")
    require(actual_content == 55175120,
            "physical raw content changed")
    require(primitive_denominator == 607037547933467614874742741,
            "primitive denominator changed")
    require(digest == "058e45bd039e8235add3d6b0eabbad0e57a90bfe9f2e4eebfc07e033c0f3fe8c",
            "primitive tensor digest changed")
    require(y_digest == "5dd45705d4bd5e3654aa8152e6da75753e8b8060ae44ecc668a42a0e7e1fa191",
            "Bockstein-row digest changed")
    require(beta_nonzero == 10452 and len(zero_beta_slices) == 364,
            "Bockstein census changed")
    require(unit_slices == 1740,
            "unit-slice census changed")
    require(fixed_q_unit_cells == [84, 84, 81, 83, 83, 83, 84,
                                   84, 83, 83, 83, 81, 84],
            "fixed-section unit coverage changed")
    require(all(cell_positive.values()) and all(cell_beta.values())
            and all(cell_units.values()),
            "one base cell lost every middle-edge attachment")
    expected_theta_one = (
        {(6, ell) for ell in range(Q7)}
        | {(11, ell) for ell in range(Q7)}
        | {(8, 2)}
    )
    require({cell for cell, t in q0_selector.items() if t == 0}
            == expected_theta_one,
            "uniform-q selector pattern changed")
    require(selected_content == content and selected_unit_count == 84,
            "selected q=0 bank lost global primitivity or a unit")
    require(selected_digest
            == "ab59b257dd3a9788d6108a46f1e01705caa2598fc9c47cee95f26165a0242110",
            "selected q=0 bank digest changed")

    # Coefficient-level projective owner-cycle support test.  In the affine
    # root chart of THM-2603, the owner map has two seven-cycles; the second
    # passes through the projective boundary infinity, for which this bank has
    # no target section.  Rotate each cycle against the seven clock positions
    # and ask only whether every finite vertex admits a positive unit rail.
    # This is deliberately not a composability or physical-intertwiner test.
    projective_cycles = {
        "O0": (0, 4, 6, 10, 7, 5, 3),
        "O1_finite": (1, 8, None, 2, 9, 12, 11),
    }
    cycle_min_theta_one = {}
    cycle_best_phases = {}
    for name, cycle in projective_cycles.items():
        minima = []
        phases = []
        for s in range(1, P):
            phase_candidates = []
            for phase in range(Q7):
                forced_theta_one = 0
                admissible = True
                for ell in range(Q7):
                    q = cycle[(ell + phase) % Q7]
                    if q is None:
                        continue
                    choices = [
                        t for j, t in by_cell[(s, ell)]
                        if pair_totals[j * P + q] > 0 and is_unit[j][q]
                    ]
                    if not choices:
                        admissible = False
                        break
                    forced_theta_one += int(12 not in choices)
                if admissible:
                    phase_candidates.append((forced_theta_one, phase))
            require(phase_candidates,
                    f"projective owner cycle lost finite unit support: {name}, s={s}")
            best = min(phase_candidates)
            minima.append(best[0])
            phases.append(best[1])
        cycle_min_theta_one[name] = tuple(minima)
        cycle_best_phases[name] = tuple(phases)

    require(cycle_min_theta_one["O0"]
            == (0, 0, 0, 0, 0, 5, 0, 1, 0, 0, 7, 0),
            "first projective owner-cycle theta invoice changed")
    require(cycle_min_theta_one["O1_finite"]
            == (0, 0, 0, 0, 0, 4, 0, 1, 0, 0, 6, 0),
            "second projective owner-cycle theta invoice changed")
    require({s for s in range(1, P)
             if cycle_min_theta_one["O0"][s - 1]}
            == {6, 8, 11},
            "projective owner-cycle exceptional displacement set changed")

    print("== constant-v=6 middle-rail common-x pullback ==")
    print(f"common row/grid: {module.W} / {T}")
    print("middle-edge zero sets:", zero_sets[(6, 0)], zero_sets[(6, 12)])
    print("selector theta histogram:", dict(sorted(theta_hist.items())))
    print("digit controls: h=0 and h=12 empty; h=6 positive in all seven phases")
    print(f"selected rails: {len(rail_bank)}; min numerator {min(selected_mass)}")
    print(f"positive rail x target-slice pairs: {len(pair_totals)-len(zero_pairs)}/{len(pair_totals)}")
    zero_pair_groups = {}
    for label, q in zero_pairs:
        zero_pair_groups.setdefault(label, []).append(q)
    print("zero-pair groups:", sorted(zero_pair_groups.items()))
    positive_by_rail = [sum(pair_totals[j * P + q] > 0 for q in range(P))
                        for j in range(nrails)]
    positive_by_q = [sum(pair_totals[j * P + q] > 0 for j in range(nrails))
                     for q in range(P)]
    print("positive target sections per rail histogram:",
          sorted(Counter(positive_by_rail).items()))
    print("positive rails by fixed q:", positive_by_q)
    print(f"positive fine entries: {positive_entries}/{nrails*P*Q7*(P-1)}")
    print("pair-support histogram:", sorted(Counter(pair_supports).items()))
    print(f"pre-route numerator global gcd: {content}; v13={valuation13}")
    print(f"actual raw global content: {actual_content}")
    print(f"primitive common denominator: {primitive_denominator}")
    print(f"primitive tensor digest: {digest}")
    print(f"primitive Bockstein-row digest: {y_digest}")
    print(f"Bockstein owner profiles nonzero: {beta_nonzero}/{nrails*P*(Q7-1)}")
    print("good-Bockstein sections per rail histogram:",
          sorted(Counter(good_beta_by_rail).items()))
    print("good-Bockstein rails by fixed q:", good_beta_by_q)
    print("zero-Bockstein slice count:", len(zero_beta_slices))
    print(f"unit slices: {unit_slices}/{nrails*P}")
    print("unit sections per rail histogram:",
          sorted(Counter(unit_by_rail).items()))
    print("unit rails by fixed q:", unit_by_q)
    nonzero_nonunits = [x for x in nonunit_slices if any(x[2])]
    print("nonzero nonunit slices:", nonzero_nonunits)
    print("cell-level positive-choice histogram:",
          sorted(Counter(map(len, cell_positive.values())).items()))
    print("cell-level good-Bockstein-choice histogram:",
          sorted(Counter(map(len, cell_beta.values())).items()))
    print("cell-level unit-choice histogram:",
          sorted(Counter(map(len, cell_units.values())).items()))
    print("cells with no positive choice:",
          sorted(cell for cell, choices in cell_positive.items() if not choices))
    print("cells with no good-Bockstein choice:",
          sorted(cell for cell, choices in cell_beta.items() if not choices))
    print("cells with no unit choice:",
          sorted(cell for cell, choices in cell_units.items() if not choices))
    print("cells admitting a unit by fixed q:", fixed_q_unit_cells)
    q0_zero = sorted(cell for cell, t in q0_selector.items() if t == 0)
    print("q=0 selector: theta one exactly on s in {6,11} plus (s,ell)=(8,2);")
    print("  theta zero on the other 69 cells")
    print("q=0 theta-one t=0 cells:", q0_zero)
    print(f"q=0 selected-bank gcd: {selected_content}; ratio to full {selected_content // content}")
    print(f"q=0 selected-bank units: {selected_unit_count}/84")
    print(f"q=0 selected-bank digest: {selected_digest}")
    print("projective owner-cycle coefficient support: O0 all 7 finite; "
          "O1 all 6 finite, with infinity absent from the target-section atlas")
    print("projective cycle minimum theta-one switches, s=1..12:",
          cycle_min_theta_one)
    print("projective cycle first minimizing phases, s=1..12:",
          cycle_best_phases)
    print("projective cycle exceptional displacements: {6,8,11}; "
          "coefficient support only, no ordered gluing")
    print("determinant histogram:", sorted(determinant_hist.items()))
    print("independent controls: "
          f"comb={comb_controls}, grouped zero/positive="
          f"{grouped_zero_controls}/{grouped_positive_controls}")


if __name__ == "__main__":
    main()
