#!/usr/bin/env python3
"""Exact carrier and cohomology controls for THM-2607."""

from collections import Counter
from itertools import product

import lrc14_fallback_rail_digit_diagonal_pullback_thm2592 as old


P = 13
N = 7
V = 6
GLOBAL_CONTENT = 4244240  # THM-2600's one global primitive reduction


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def selected_t(s, ell):
    """THM-2600 (16), with t=0 the theta-one middle rail."""
    return 0 if s in (6, 11) or (s, ell) == (8, 2) else 12


def rebuild_q0_unit_choices():
    """Rebuild only the q=0 face of THM-2600's exact common-x bank."""
    module = old.load_present_module()
    base = old.base
    host = old.host
    rail = old.rail
    sat = old.sat
    top = old.T

    require(module.W == base.W, "row mismatch")

    word = module.build_word_Ta()
    word_prefix = []
    for ell in range(N):
        q_ell = module.subtract_comb(
            word, module.C1, 182, 26 * ell - 13, 26 * ell + 13
        )
        digit = sat.intersect_interval(
            q_ell, V * top // P, (V + 1) * top // P
        )
        prefix = module.make_prefix(digit)
        require(prefix[2][-1] > 0, "future digit six disappeared")
        word_prefix.append(prefix)

    e4 = base.build_set(base.PAT_E, base.ZELL)
    qb = base.build_set(host.PAT_QB, base.ZELL)
    ust, uv, vst, vv = rail.packet_profiles(e4, qb)
    _, _, k_tensor = rail.exact_tensor(e4, qb)
    owner = base.clock_cells(P, N, top, P * P)
    deep = rail.deep_cells()
    arrival = [[(V * (top // P), (V + 1) * (top // P))]]

    rail_bank = []
    for s in range(1, P):
        rvst, rvv = base.rotate_profile(vst, vv, s * (top // P), top)
        ps, pv, _, _ = base.product_cum(ust, uv, rvst, rvv, top)
        for ell in range(N):
            for t in (12, 0):
                if k_tensor[s][V][t][ell] == 0:
                    continue
                cell = old.intersect_sorted(
                    old.intersect_sorted(owner[ell], deep[t]), arrival[0]
                )
                pieces = old.profile_on_intervals(cell, ps, pv)
                numerator = P * sum(weight * (right - left)
                                    for left, right, weight in pieces)
                require(numerator == k_tensor[s][V][t][ell] > 0,
                        "middle-rail mass mismatch")
                rail_bank.append((s, ell, t, pieces))
    require(len(rail_bank) == 162, "middle-rail census changed")

    # q=0 means s5=0.  No other target section is built here.
    present = {}
    present_starts = {}
    for ell5 in range(N):
        f = module.build_F(ell5, 0)
        present[ell5] = f
        present_starts[ell5] = [left for left, _ in f]

    choices = {(s, ell): [] for s in range(1, P) for ell in range(N)}
    unit_edges = 0
    divisibility_checks = 0
    for s, ell, t, pieces in rail_bank:
        y = []
        total = 0
        for ell5 in range(N):
            overlap = old.intersect_weighted_union(
                pieces, present[ell5], present_starts[ell5]
            )
            value_by_r = [0] * P
            if overlap:
                for r in range(1, P):
                    probed = old.intersect_weighted_comb(
                        overlap, module.C3, 182, 14 * r - 13, 14 * r + 13
                    )
                    value = old.delayed_weighted_numerator(
                        probed, word_prefix[ell5]
                    )
                    require(value % GLOBAL_CONTENT == 0,
                            "THM-2600 global content no longer divides q=0")
                    divisibility_checks += 1
                    value_by_r[r] = value
                    total += value
            y.append(sum(
                (value_by_r[r] // GLOBAL_CONTENT) * pow(r, -1, P)
                for r in range(1, P)
            ) % P)
        reduced = tuple((y[i] - y[-1]) % P for i in range(N - 1))
        determinant = old.sat.multiplication_determinant_7(reduced)
        if determinant != 0:
            require(total > 0, "unit Bockstein has zero physical mass")
            choices[s, ell].append(t)
            unit_edges += 1

    table = {key: tuple(sorted(value)) for key, value in choices.items()}
    require(unit_edges == 138, "q=0 unit-edge census changed")
    require(divisibility_checks == 4752,
            "q=0 divisibility-check census changed")
    return table, unit_edges, divisibility_checks


def endpoint_data(s):
    ts = tuple(selected_t(s, ell) for ell in range(N))
    ws = tuple((7 * t) % P for t in ts)
    cs = tuple((w - V) % P for w in ws)
    return ts, ws, cs, sum(cs) % P


def gauges(a, cs):
    """All h with a+c_ell+h_ell-h_(ell-1)=0, edge ell-1 -> ell."""
    out = []
    for base in range(P):
        h = [None] * N
        previous = base  # h_6
        for ell in range(N - 1):
            h[ell] = (previous - a - cs[ell]) % P
            previous = h[ell]
        h[N - 1] = base
        for ell in range(N):
            require((a + cs[ell] + h[ell] - h[(ell - 1) % N]) % P == 0,
                    "gauge equation failed")
        out.append(tuple(h))
    require(len(set(out)) == P, "gauge torsor does not have size thirteen")
    return tuple(out)


def shifted_digit_support(m):
    """Digits met by h=6 after x -> x+m/91 at response R=13^6."""
    # On the denominator-91 phase circle, h=6 is [42,49).  Since
    # 13^5 == -1 mod 7, the response phase is shifted by -13m/91.
    start = (42 - 13 * m) % 91
    pieces = ((start, start + 7),) if start + 7 <= 91 else (
        (start, 91), (0, start + 7 - 91)
    )
    hit = []
    for h in range(P):
        left, right = 7 * h, 7 * (h + 1)
        if any(max(left, a) < min(right, b) for a, b in pieces):
            hit.append(h)
    return tuple(hit)


def koopman_translation_image(n, numerator):
    """Image of x -> x+numerator/91 under T^n, in denominator seven."""
    require(n >= 1, "Koopman delay must be positive")
    # 13^n*numerator/91 = numerator*13^(n-1)/7.  Returning the
    # numerator modulo seven makes the C91 -> C7 rebasing explicit.
    return (numerator * pow(13, n - 1, 7)) % 7


def main():
    unit_table, unit_edges, divisibility_checks = rebuild_q0_unit_choices()

    both = (0, 12)
    only_zero = (0,)
    only_twelve = (12,)
    expected_table = {}
    for s in (1, 3, 4, 9, 10, 12):
        for ell in range(N):
            expected_table[s, ell] = both
    for s in (2, 7):
        for ell in range(N):
            expected_table[s, ell] = only_twelve
    for ell in range(N):
        expected_table[5, ell] = only_twelve if ell == 5 else both
        expected_table[6, ell] = only_zero
        expected_table[8, ell] = only_zero if ell == 2 else both
        expected_table[11, ell] = only_zero
    require(unit_table == expected_table, "q=0 unit-choice table changed")

    correction_sums = {}
    theta_one_counts = {}
    matched = []
    reverse_matched = []
    gauge_count = 0

    for s in range(1, P):
        ts, ws, cs, total = endpoint_data(s)
        require(all(t in (0, 12) for t in ts), "selector left middle rails")
        require(all(w in (0, 6) for w in ws), "unexpected rail endpoint")
        require(all(c in (0, 7) for c in cs), "unexpected rail correction")
        theta_one_counts[s] = ts.count(0)
        correction_sums[s] = total
        for a in range(1, P):
            if total == (-N * a) % P:
                matched.append((s, a))
                gs = gauges(a, cs)
                gauge_count += len(gs)
                if s in (6, 11):
                    require(all((a + c) % P == 0 for c in cs),
                            "constant-seven lane should cancel edgewise")
            if (-total) % P == (-N * a) % P:
                reverse_matched.append((s, a))

    require(theta_one_counts == {
        1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 7,
        7: 0, 8: 1, 9: 0, 10: 0, 11: 7, 12: 0,
    }, "theta-one count profile changed")
    require(correction_sums == {
        1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 10,
        7: 0, 8: 7, 9: 0, 10: 0, 11: 10, 12: 0,
    }, "rail correction profile changed")
    require(matched == [(6, 6), (8, 12), (11, 6)],
            "forward marker matching changed")
    require(reverse_matched == [(6, 7), (8, 1), (11, 7)],
            "reverse marker matching changed")
    require(gauge_count == 39, "wrong total gauge count")

    attainable_counts = {}
    selector_count = 0
    nontrivial_selectors = 0
    marker_hist = Counter()
    count_level_lanes = []
    all_selector_gauges = 0
    for s in range(1, P):
        counts = set()
        local_choices = [unit_table[s, ell] for ell in range(N)]
        for ts in product(*local_choices):
            selector_count += 1
            n = ts.count(0)
            counts.add(n)
            if n:
                nontrivial_selectors += 1
                a = (-n) % P
                marker_hist[a] += 1
                cs = tuple(7 if t == 0 else 0 for t in ts)
                require(sum(cs) % P == (-N * a) % P,
                        "full selector failed marker invoice")
                all_selector_gauges += len(gauges(a, cs))
        attainable_counts[s] = tuple(sorted(counts))
        for n in sorted(counts - {0}):
            count_level_lanes.append((s, (-n) % P, n))
    require(selector_count == 900 and nontrivial_selectors == 891,
            "q=0 selector-polytope census changed")
    require(len(count_level_lanes) == 57,
            "q=0 count-level invoice census changed")
    require(marker_hist == Counter({
        6: 9, 7: 49, 8: 147, 9: 245,
        10: 245, 11: 147, 12: 49,
    }), "marker multiplicity histogram changed")
    require(all_selector_gauges == 891 * P,
            "full selector gauge count changed")

    expected_digit_pairs = {
        1: (4, 5), 2: (2, 3), 3: (0, 1),
        4: (11, 12), 5: (9, 10), 6: (7, 8),
    }
    digit_pairs = {m: shifted_digit_support(m) for m in range(1, N)}
    require(digit_pairs == expected_digit_pairs,
            "natural C91 digit-transport table changed")
    require(all(V not in pair for pair in digit_pairs.values()),
            "natural chart rotation unexpectedly preserves digit six")
    arrival_vertices = {0, 6, 12}
    require({m: tuple(sorted(set(pair) & arrival_vertices))
             for m, pair in digit_pairs.items()} == {
                 1: (), 2: (), 3: (0,), 4: (12,), 5: (), 6: ()
             }, "arrival/future toothpick match table changed")

    # Chronological transport does not preserve the root deck.  A root
    # translation q/13 is the C91 numerator 7q and is killed by every
    # positive iterate.  A general C91 shift becomes the alternating C7
    # owner phase (-1)^(n-1)m/7.
    for n in range(1, 15):
        require(all(koopman_translation_image(n, 7 * q) == 0
                    for q in range(P)),
                "positive Koopman time did not erase the C13 deck")
        require(tuple(koopman_translation_image(n, m)
                      for m in range(1, N)) == tuple(
                          (((-1) ** (n - 1)) * m) % N
                          for m in range(1, N)),
                "C91 chart shift did not rebase to alternating C7 phase")

    # At (s,ell)=(1,0), both t=12 (c=0) and t=0 (c=7) are positive
    # open rail cells.  Any nonempty open interval maps onto the circle
    # under all sufficiently large powers of T, so the same future point
    # has antecedents with both rail classes.  The finite assertion needed
    # here is that the two hostile cells really occur in the exact bank.
    require(unit_table[1, 0] == (0, 12),
            "two-class chronological hostile lost a positive rail")

    # The s=8 lane is globally but not edgewise flat.  Its combined
    # transition has one value 6 and six values 12, of total zero.
    _, _, c8, _ = endpoint_data(8)
    combined8 = tuple((12 + c) % P for c in c8)
    require(combined8 == (12, 12, 6, 12, 12, 12, 12),
            "s=8 combined transition profile changed")
    require(sum(combined8) % P == 0 and any(combined8),
            "s=8 should be exact but not pointwise zero")

    print("THM-2607 exact constant-six rail-boundary controls")
    print(f"q0_unit_edges={unit_edges}/162 "
          f"global_content_divisibility_checks={divisibility_checks}")
    print("q0_attainable_theta_one_counts=" + str(tuple(
        attainable_counts[s] for s in range(1, P))))
    print(f"q0_selector_polytope={selector_count} "
          f"nontrivial_invoices={nontrivial_selectors} "
          f"count_level_lanes={len(count_level_lanes)}")
    print("forward_marker_selector_hist=" + str(tuple(
        marker_hist[a] for a in range(1, P))))
    print("theta_one_counts=" + str(tuple(theta_one_counts[s]
                                         for s in range(1, P))))
    print("correction_sums=" + str(tuple(correction_sums[s]
                                       for s in range(1, P))))
    print(f"forward_invoice_matches={matched}")
    print(f"reverse_orientation_matches={reverse_matched}")
    print(f"forward_gauges_checked={gauge_count} gauges_per_lane=13")
    print(f"all_selector_gauges_checked={all_selector_gauges}")
    print(f"s8_combined_transition={combined8} sum={sum(combined8) % P}")
    print(f"natural_rotation_digit_pairs={digit_pairs}")
    print("koopman_root_deck_images=" + str(tuple(
        koopman_translation_image(n, 7) for n in range(1, 8))) +
        " c91_generator_images=" + str(tuple(
            koopman_translation_image(n, 1) for n in range(1, 8))))
    print("verdict=PASS: physical rail cochain pays seven forward marker "
          "classes abstractly; natural C91/chart-time transport fails")


if __name__ == "__main__":
    main()
