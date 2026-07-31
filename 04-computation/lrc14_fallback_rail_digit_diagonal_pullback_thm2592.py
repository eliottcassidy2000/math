#!/usr/bin/env python3
"""Exact companion for THM-2592's THM-2584/2586 x THM-2585 pullback.

It keeps the route-two physical coordinate ``x`` before either theorem
marginalizes.  On the THM-2586-selected theta=0 rail it intersects

* the b-word depth-five Perron packet and its (ell4,v,t) cell; and
* the THM-2585 literal target slice s5=-q, source phase ell5, deep probe r,
  and the old=future digit diagonal at the fixed clock 13^6.

Because v=floor(13x), the THM-2585 digit restriction uses h=v and enforces
floor(13 {13^6 x})=v.  The final hostile pass replaces h=v by every affine
graph h=v+c and, independently, tests each primary rail against the entire
unconditioned delayed word.  All arithmetic is integer/Fraction arithmetic.
One global content is removed only after the complete joint bank is built.
"""

from bisect import bisect_right
from collections import Counter
from fractions import Fraction
from math import gcd
from pathlib import Path
import importlib.util
import hashlib
import sys


HERE = Path(__file__).resolve().parents[1] / "04-computation"
sys.path.insert(0, str(HERE))

import lrc14_base_only_bridge_opus_20260728 as base
import lrc14_b_r5_owner_clock_host_thm2581 as host
import lrc14_b_r5_theta_target_tensor_thm2584 as rail
import lrc14_saturated_target_projector_bockstein_thm2585 as sat


P = 13
Q7 = 7
T = base.T_DEN
R = P**6


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def load_present_module():
    path = HERE / "lrc14_replica_dichotomy_typed_row_opus_20260727.py"
    spec = importlib.util.spec_from_file_location("joint_2584_2585_present", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    require(module.T == T, "base-grid mismatch")
    return module


def intersect_sorted(left, right):
    out = []
    i = j = 0
    while i < len(left) and j < len(right):
        a = max(left[i][0], right[j][0])
        b = min(left[i][1], right[j][1])
        if a < b:
            out.append((a, b))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return out


def profile_on_intervals(intervals, starts, values):
    """Restrict a nonnegative step profile to a sorted interval union."""
    out = []
    for a, b in intervals:
        i = bisect_right(starts, a) - 1
        x = a
        while x < b:
            y = min(b, starts[i + 1] if i + 1 < len(starts) else T)
            w = values[i]
            if w:
                if out and out[-1][1] == x and out[-1][2] == w:
                    out[-1] = (out[-1][0], y, w)
                else:
                    out.append((x, y, w))
            x = y
            if i + 1 < len(starts) and x == starts[i + 1]:
                i += 1
    return out


def intersect_weighted_union(pieces, intervals, starts=None):
    """Intersect sorted weighted pieces with a sorted interval union."""
    if starts is None:
        starts = [a for a, _ in intervals]
    out = []
    for a, b, w in pieces:
        j = max(0, bisect_right(starts, a) - 1)
        while j < len(intervals) and intervals[j][0] < b:
            c = max(a, intervals[j][0])
            d = min(b, intervals[j][1])
            if c < d:
                out.append((c, d, w))
            j += 1
    return out


def intersect_weighted_comb(pieces, speed, pd, lo, hi):
    """Intersect with nonwrapping comb windows, retaining integer weights."""
    require(0 <= lo < hi <= pd, "this probe only uses nonwrapping r!=0 windows")
    require(T % (pd * speed) == 0, "comb is not resolved by the base grid")
    unit = T // (pd * speed)
    step = pd * unit
    width = (hi - lo) * unit
    origin = lo * unit
    out = []
    for a, b, w in pieces:
        # First p=origin+k*step whose right endpoint is strictly after a.
        k = (a - width - origin) // step + 1
        p = origin + k * step
        while p < b:
            c = max(a, p)
            d = min(b, p + width)
            if c < d:
                out.append((c, d, w))
            p += step
    return out


def phi_at(x, starts, lens, pref):
    i = bisect_right(starts, x) - 1
    if i < 0:
        return 0
    d = x - starts[i]
    return pref[i] + min(d, lens[i])


def delayed_weighted_numerator(pieces, prefix):
    """Numerator over R*T of integral weight(x) 1_Q(Rx) dx."""
    starts, lens, pref = prefix
    length_q = pref[-1]
    weighted_len = 0
    acc_r = 0
    acc_phi = 0
    rred = R % T
    for a, b, w in pieces:
        ra = a * rred % T
        rb = b * rred % T
        weighted_len += w * (b - a)
        acc_r += w * (rb - ra)
        acc_phi += w * (phi_at(rb, starts, lens, pref) - phi_at(ra, starts, lens, pref))
    floor_numerator = R * weighted_len - acc_r
    require(floor_numerator % T == 0, "weighted floor count is not integral")
    result = length_q * (floor_numerator // T) + acc_phi
    require(result >= 0, "negative overlap numerator")
    return result


def grouped_original_numerator(module, pieces, prefix):
    """Independent linear route through THM-2585's unweighted sweep."""
    starts, lens, pref = prefix
    by_weight = {}
    for a, b, w in pieces:
        by_weight.setdefault(w, []).append((a, b))
    total = 0
    for w, intervals in by_weight.items():
        ar, ap = module.sweep_acc(intervals, R % T, starts, lens, pref)
        value = module.IR_from_acc(
            R, sat.interval_length(intervals), pref[-1], ar, ap
        )
        scaled = value * R * T
        require(scaled.denominator == 1, "original sweep did not clear")
        total += w * scaled.numerator
    return total


def validate_delayed_formula(module, word_prefix):
    """Recover selected raw THM-2585 entries before adding the rail weight."""
    for ell, s, r in ((0, 0, 1), (3, 7, 8), (6, 12, 12)):
        f = module.build_F(ell, s)
        e = module.intersect_comb(f, module.C3, 182, 14 * r - 13, 14 * r + 13)
        total = 0
        for h in range(P):
            eh = sat.intersect_interval(e, h * T // P, (h + 1) * T // P)
            total += delayed_weighted_numerator([(a, b, 1) for a, b in eh], word_prefix[ell][h])
        starts, lens, pref = word_prefix[ell][0]
        # Directly compare with the theorem's original sweep, digit by digit.
        original = 0
        for h in range(P):
            eh = sat.intersect_interval(e, h * T // P, (h + 1) * T // P)
            starts, lens, pref = word_prefix[ell][h]
            ar, ap = module.sweep_acc(eh, R % T, starts, lens, pref)
            frac = module.IR_from_acc(R, sat.interval_length(eh), pref[-1], ar, ap)
            original += frac * R * T
        require(original.denominator == 1 and total == original.numerator,
                "weighted delayed-overlap formula does not recover THM-2585")


def main():
    module = load_present_module()
    require(module.W == base.W, "row mismatch")

    # THM-2585 delayed word, split by its old/future root digit.
    word = module.build_word_Ta()
    word_prefix = [[None] * P for _ in range(Q7)]
    whole_word_prefix = [None] * Q7
    for ell in range(Q7):
        q_ell = module.subtract_comb(word, module.C1, 182, 26 * ell - 13, 26 * ell + 13)
        whole_word_prefix[ell] = module.make_prefix(q_ell)
        for h in range(P):
            digit = sat.intersect_interval(q_ell, h * T // P, (h + 1) * T // P)
            word_prefix[ell][h] = module.make_prefix(digit)
    digit_lengths = [[prefix[2][-1] for prefix in row] for row in word_prefix]
    require(
        all(sum(digit_lengths[ell]) == whole_word_prefix[ell][2][-1]
            for ell in range(Q7)),
        "future-digit pieces do not partition the delayed word",
    )
    require(all(digit_lengths[ell][0] == 0 for ell in range(Q7)),
            "the exact future-digit-zero hole disappeared")
    require(all(digit_lengths[ell][6] > 0 for ell in range(Q7)),
            "the exact future-digit-six support disappeared")
    validate_delayed_formula(module, word_prefix)

    # THM-2584 route-two profile and exact selected theta=0 rails.
    e4 = base.build_set(base.PAT_E, base.ZELL)
    qb = base.build_set(host.PAT_QB, base.ZELL)
    ust, uv, vst, vv = rail.packet_profiles(e4, qb)
    _, _, k_tensor = rail.exact_tensor(e4, qb)
    owner = base.clock_cells(P, Q7, T, P * P)
    deep = rail.deep_cells()
    arrival = [[(v * (T // P), (v + 1) * (T // P))] for v in range(P)]

    rail_bank = []
    selected_mass = []
    for s4 in range(1, P):
        rvst, rvv = base.rotate_profile(vst, vv, s4 * (T // P), T)
        ps, pv, _, _ = base.product_cum(ust, uv, rvst, rvv, T)
        for ell4 in range(Q7):
            choices = [(v, t) for v, t in ((0, 0), (6, 12)) if k_tensor[s4][v][t][ell4] > 0]
            require(len(choices) >= 1, "selected theta-zero rail disappeared")
            v, t = choices[0]
            cell = intersect_sorted(intersect_sorted(owner[ell4], deep[t]), arrival[v])
            pieces = profile_on_intervals(cell, ps, pv)
            numerator = P * sum(w * (b - a) for a, b, w in pieces)
            require(numerator == k_tensor[s4][v][t][ell4], "rail piece mass mismatch")
            rail_bank.append((s4, ell4, v, t, pieces))
            selected_mass.append(numerator)
    require(len(rail_bank) == 84, "rail census changed")

    # Cache the 91 present-side target/source packets.
    f_cache = {}
    f_starts = {}
    for s5 in range(P):
        for ell5 in range(Q7):
            f = module.build_F(ell5, s5)
            f_cache[ell5, s5] = f
            f_starts[ell5, s5] = [a for a, _ in f]

    # joint[rail_index][q][ell5][r], with r=0 identically zero because the
    # present packet retains c3-safe.  The stored integer is the numerator
    # over R*T before the global route-two scalar P/(169*13^10).
    joint = [[[[0] * P for _ in range(Q7)] for _ in range(P)] for _ in rail_bank]
    positive_pairs = 0
    pair_totals = []
    content = 0
    positive_entries = 0
    grouped_zero_controls = 0
    grouped_positive_controls = 0
    comb_controls = 0
    primary_whole_totals = {}
    fallback_digit_totals = {}
    for j, (s4, ell4, v, t, rail_pieces) in enumerate(rail_bank):
        for q in range(P):
            s5 = (-q) % P
            total = 0
            primary_whole_total = 0
            fallback_by_h = [0] * P
            for ell5 in range(Q7):
                present = intersect_weighted_union(
                    rail_pieces, f_cache[ell5, s5], f_starts[ell5, s5]
                )
                if not present:
                    continue
                prefix = word_prefix[ell5][v]
                for r in range(1, P):
                    probed = intersect_weighted_comb(
                        present, module.C3, 182, 14 * r - 13, 14 * r + 13
                    )
                    if v == 0:
                        primary_whole_total += delayed_weighted_numerator(
                            probed, whole_word_prefix[ell5]
                        )
                        value = delayed_weighted_numerator(probed, prefix)
                    else:
                        values_by_h = [
                            delayed_weighted_numerator(probed, word_prefix[ell5][h])
                            for h in range(P)
                        ]
                        for h, h_value in enumerate(values_by_h):
                            fallback_by_h[h] += h_value
                        value = values_by_h[v]
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
                        require(expected_len == sum(w * (b - a) for a, b, w in probed),
                                "direct comb route disagrees with complement route")
                        comb_controls += 1
                    need_grouped_control = (
                        probed
                        and ((value == 0 and grouped_zero_controls < 8)
                             or (value > 0 and grouped_positive_controls < 8))
                    )
                    if need_grouped_control:
                        require(
                            value == grouped_original_numerator(module, probed, prefix),
                            "weighted sweep disagrees with grouped original route",
                        )
                        if value:
                            grouped_positive_controls += 1
                        else:
                            grouped_zero_controls += 1
                    joint[j][q][ell5][r] = value
                    if value:
                        total += value
                        positive_entries += 1
                        content = gcd(content, value)
            label = (s4, ell4, v, t, q)
            if v == 0:
                primary_whole_totals[label] = primary_whole_total
            else:
                fallback_digit_totals[label] = tuple(fallback_by_h)
            pair_totals.append(total)
            positive_pairs += int(total > 0)
        if (j + 1) % 7 == 0:
            print(f"joint rail progress {j + 1}/84", flush=True)

    require(content > 0, "joint bank is identically zero")
    require(comb_controls == 3 and grouped_positive_controls == 8,
            "independent positive route-control census changed")
    zero_pairs = [(rail_bank[j][:4], q) for j in range(84) for q in range(P)
                  if sum(joint[j][q][ell][r] for ell in range(Q7) for r in range(1, P)) == 0]
    positive_pair_data = [
        (rail_bank[j][:4], q,
         sum(joint[j][q][ell][r] for ell in range(Q7) for r in range(1, P)))
        for j in range(84) for q in range(P)
        if sum(joint[j][q][ell][r] for ell in range(Q7) for r in range(1, P)) > 0
    ]
    pair_supports = [
        sum(joint[j][q][ell][r] != 0 for ell in range(Q7) for r in range(1, P))
        for j in range(84) for q in range(P)
        if sum(joint[j][q][ell][r] for ell in range(Q7) for r in range(1, P)) > 0
    ]

    # Union-first affine-offset hostile.  The 81 primary rails are orthogonal
    # to the whole delayed word, so every one of its thirteen nonnegative
    # digit pieces (and every selector subordinate to them) vanishes there.
    require(len(primary_whole_totals) == 81 * P,
            "primary whole-word pair census changed")
    require(not any(primary_whole_totals.values()),
            "a primary rail meets the unconditioned delayed word")
    require(len(fallback_digit_totals) == 3 * P,
            "fallback digit pair census changed")
    fallback_positive_by_h = Counter()
    for totals in fallback_digit_totals.values():
        for h, h_value in enumerate(totals):
            fallback_positive_by_h[h] += int(h_value > 0)
    require(
        fallback_positive_by_h == Counter({h: 39 for h in range(1, 12)}),
        "fallback future-digit support changed",
    )
    offset_positive_pairs = []
    for c in range(P):
        # Every fallback rail has v=6.  Primary rails vanish before splitting.
        offset_positive_pairs.append(fallback_positive_by_h[(6 + c) % P])
    expected_offset_pairs = [39] * P
    expected_offset_pairs[6] = expected_offset_pairs[7] = 0
    require(offset_positive_pairs == expected_offset_pairs,
            "affine future-digit offset spectrum changed")

    # Remove exactly one global content, then recompute every first Bockstein.
    valuation13 = 0
    temp = content
    while temp % P == 0:
        valuation13 += 1
        temp //= P
    beta_nonzero = 0
    beta_total = 84 * P * (Q7 - 1)
    unit_slices = 0
    determinant_hist = Counter()
    zero_beta_profiles = []
    y_digest_rows = []
    positive_slice_data = []
    fallback_y_rows = {4: {}, 5: {}, 6: {}}
    for j in range(84):
        for q in range(P):
            y = []
            for ell in range(Q7):
                value = sum(
                    (joint[j][q][ell][r] // content) * pow(r, -1, P)
                    for r in range(1, P)
                ) % P
                y.append(value)
            reduced = tuple((y[i] - y[-1]) % P for i in range(Q7 - 1))
            det = sat.multiplication_determinant_7(reduced)
            determinant_hist[det] += 1
            unit_slices += int(det != 0)
            if any(y):
                positive_slice_data.append((rail_bank[j][:4], q, tuple(y), det))
            s4, ell4, v, t = rail_bank[j][:4]
            if (s4, v, t) == (7, 6, 12) and ell4 in fallback_y_rows:
                fallback_y_rows[ell4][q] = tuple(y)
            for kappa in range(1, Q7):
                nz = any(sat.owner_factor(y, kappa))
                beta_nonzero += int(nz)
                if not nz:
                    zero_beta_profiles.append((rail_bank[j][:4], q, kappa, tuple(y)))
            y_digest_rows.append(tuple(y))

    # Boolean q-section aggregates on each actual fallback carrier.  These
    # are sums after the one global primitive-content reduction, never a
    # slice-by-slice reprimitivization.
    boolean_atlas = {}
    for ell4, rows in fallback_y_rows.items():
        require(set(rows) == set(range(P)), "fallback q-row atlas incomplete")
        fibres = Counter()
        zero_subsets = []
        unit_subsets = 0
        owner_nonzero = 0
        for mask in range(1 << P):
            raw = tuple(
                sum(rows[q][ell] for q in range(P) if mask >> q & 1) % P
                for ell in range(Q7)
            )
            reduced = tuple((raw[ell] - raw[-1]) % P for ell in range(Q7 - 1))
            fibres[reduced] += 1
            if reduced == (0,) * (Q7 - 1):
                zero_subsets.append(tuple(q for q in range(P) if mask >> q & 1))
            unit_subsets += int(sat.multiplication_determinant_7(reduced) != 0)
            if reduced != (0,) * (Q7 - 1):
                for kappa in range(1, Q7):
                    require(any(sat.owner_factor(raw, kappa)),
                            "nonzero section aggregate lost an owner colour")
                    owner_nonzero += 1
        boolean_atlas[ell4] = (
            len(fibres), tuple(zero_subsets), max(fibres.values()),
            unit_subsets, (1 << P) - unit_subsets, owner_nonzero,
        )
    require(
        boolean_atlas == {
            4: (864, ((),), 40, 8118, 74, 8191 * 6),
            5: (516, ((),), 88, 8029, 163, 8191 * 6),
            6: (328, ((), (10,)), 110, 8114, 78, 8190 * 6),
        },
        "fallback Boolean q-section atlas changed",
    )

    # Cross-chart transversals: choose exactly one q-section on each selected
    # fallback chart, always using the same global primitive normalization.
    pair_transversals = {}
    for left, right in ((4, 5), (4, 6), (5, 6)):
        fibres = Counter()
        nonzero = 0
        for q_left in range(P):
            for q_right in range(P):
                raw = tuple(
                    (fallback_y_rows[left][q_left][ell]
                     + fallback_y_rows[right][q_right][ell]) % P
                    for ell in range(Q7)
                )
                reduced = tuple((raw[ell] - raw[-1]) % P
                                for ell in range(Q7 - 1))
                fibres[reduced] += 1
                nonzero += int(reduced != (0,) * (Q7 - 1))
        pair_transversals[left, right] = (
            nonzero, len(fibres), max(fibres.values())
        )
    require(all(data[0] == P * P for data in pair_transversals.values()),
            "a two-chart singleton transversal cancelled")

    triple_fibres = Counter()
    triple_units = 0
    triple_owner_nonzero = 0
    for q4 in range(P):
        for q5 in range(P):
            for q6 in range(P):
                raw = tuple(
                    (fallback_y_rows[4][q4][ell]
                     + fallback_y_rows[5][q5][ell]
                     + fallback_y_rows[6][q6][ell]) % P
                    for ell in range(Q7)
                )
                reduced = tuple((raw[ell] - raw[-1]) % P
                                for ell in range(Q7 - 1))
                require(reduced != (0,) * (Q7 - 1),
                        "a three-chart singleton transversal cancelled")
                triple_fibres[reduced] += 1
                triple_units += int(sat.multiplication_determinant_7(reduced) != 0)
                for kappa in range(1, Q7):
                    require(any(sat.owner_factor(raw, kappa)),
                            "triple transversal lost an owner colour")
                    triple_owner_nonzero += 1
    triple_transversal = (
        sum(triple_fibres.values()), len(triple_fibres),
        max(triple_fibres.values()), triple_units,
        sum(triple_fibres.values()) - triple_units,
        triple_owner_nonzero,
    )
    require(triple_transversal == (2197, 295, 103, 2173, 24, 2197 * 6),
            "three-chart singleton transversal atlas changed")

    serialization = ";".join(
        ",".join(str(joint[j][q][ell][r] // content)
                 for ell in range(Q7) for r in range(P))
        for j in range(84) for q in range(P)
    )
    digest = hashlib.sha256(serialization.encode("ascii")).hexdigest()

    # The actual joint raw numerator has the extra route-two factor 13, so its
    # one global content is 13*content; the primitive tensor is unchanged.
    actual_content = P * content
    actual_denominator = 169 * (P**5) ** 2 * R * T
    primitive_denominator = Fraction(actual_denominator, actual_content)

    expected_positive_labels = {
        ((7, ell4, 6, 12), q) for ell4 in (4, 5, 6) for q in range(P)
    }
    require({(lab, q) for lab, q, _ in positive_pair_data} == expected_positive_labels,
            "joint support is not exactly the three fallback rails")
    require(positive_pairs == 39 and len(zero_pairs) == 1053,
            "joint pair census changed")
    require(positive_entries == 1404,
            "joint fine-support census changed")
    require(Counter(pair_supports) == Counter({36: 25, 24: 7, 48: 7}),
            "joint fine-support histogram changed")
    require(content == 25465440 and valuation13 == 1,
            "global pre-route content changed")
    require(actual_content == 331050720,
            "global physical numerator content changed")
    require(primitive_denominator == Fraction(202345849311155871624914247, 2),
            "primitive scalar changed")
    require(digest == "1286b1f2baa1f05299d93a1074db19afc3f99eb3d44500fd1de39c97c62c300c",
            "primitive joint digest changed")
    require(unit_slices == 37 and beta_nonzero == 228,
            "primitive Bockstein census changed")

    print("== exact THM-2584/2586 x THM-2585 joint pullback probe ==")
    print(f"common row/grid: {module.W} / {T}")
    print("common coordinate: route-two x; old v=floor(13x), future floor(13{13^6x})=v")
    print("exact word-digit gate: all seven Q_a owner phases have digit-0 mass 0 and digit-6 mass >0")
    print(f"selected theta-zero rails: {len(rail_bank)}, min rail numerator {min(selected_mass)}")
    print(f"positive rail x target-slice pairs: {positive_pairs}/{84 * P}")
    print(f"zero pair count: {len(zero_pairs)}")
    print("positive pair labels: s4=7, ell4 in {4,5,6}, (v,t)=(6,12), all 13 q")
    print(f"positive pair total numerator range: {min(x[2] for x in positive_pair_data)}..{max(x[2] for x in positive_pair_data)}")
    print(f"positive joint (rail,q,ell5,r) entries: {positive_entries}/{84 * P * Q7 * (P - 1)}")
    print("positive-pair support histogram: " + ",".join(
        f"{k}:{v}" for k, v in sorted(Counter(pair_supports).items())
    ))
    print("primary whole-word overlap: 0/1053; hence all digit selectors vanish")
    print("fallback positive pairs by future digit: " + ",".join(
        f"{h}:{fallback_positive_by_h[h]}" for h in range(P)
    ))
    print("affine h=v+c positive pairs: " + ",".join(
        f"{c}:{offset_positive_pairs[c]}" for c in range(P)
    ))
    print("fallback Boolean q-section atlas: " + "; ".join(
        f"ell4={ell}:images={data[0]},zeros={data[1]},max_fibre={data[2]},"
        f"units={data[3]},nonunits={data[4]},owner_nonzero={data[5]}"
        for ell, data in sorted(boolean_atlas.items())
    ))
    print("two-chart singleton transversals: " + "; ".join(
        f"{left}{right}:nonzero={data[0]},images={data[1]},max_fibre={data[2]}"
        for (left, right), data in sorted(pair_transversals.items())
    ))
    print("three-chart singleton transversals: "
          f"total={triple_transversal[0]},images={triple_transversal[1]},"
          f"max_fibre={triple_transversal[2]},units={triple_transversal[3]},"
          f"nonunits={triple_transversal[4]},owner_nonzero={triple_transversal[5]}")
    print(f"pre-route numerator global gcd: {content}; v13={valuation13}")
    print(f"actual raw global content (including route-two 13): {actual_content}")
    print(f"primitive common denominator: {primitive_denominator}")
    print(f"primitive joint digest: {digest}")
    print(f"unit slice polynomials: {unit_slices}/{84 * P}")
    print(f"Bockstein owner profiles nonzero: {beta_nonzero}/{beta_total}")
    positive_labels = {(lab, q) for lab, q, _ in positive_pair_data}
    positive_zero_beta = [x for x in zero_beta_profiles if (x[0], x[1]) in positive_labels]
    collapsed_zero_slices = sorted({(x[0], x[1], x[3]) for x in positive_zero_beta})
    require(collapsed_zero_slices == [((7, 6, 6, 12), 10, (0,) * 7)],
            "positive zero-Bockstein exception changed")
    nonunit_positive = sorted(
        (lab, q, y) for lab, q, y, det in positive_slice_data if det == 0
    )
    require(nonunit_positive == [((7, 4, 6, 12), 11, (9, 11, 0, 10, 0, 0, 0))],
            "nonunit nonzero exception changed")
    print(f"Bockstein-nonzero positive slices: {len(positive_slice_data)}/39; owner profiles 228/234")
    print(f"positive slices with zero Bockstein: {collapsed_zero_slices}")
    print("unit positive slices: 37/39; nonunit nonzero exception: " + str(nonunit_positive))
    print("independent exact route controls: "
          f"comb {comb_controls}, weighted sweep zero/positive "
          f"{grouped_zero_controls}/{grouped_positive_controls}")


if __name__ == "__main__":
    main()
