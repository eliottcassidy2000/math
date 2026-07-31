#!/usr/bin/env python3
"""Exact scout for the smallest changed delayed-base affine handoff.

Integer odometer lifts leave ``y={13^6 x}`` on the autonomous map y->{13y}.
The half lift

    A_k(x)={13x+(k+1/2)/13^6}

instead gives y->{13y+1/2}.  It is the smallest nonintegral phase shift for
which the high-speed root displacement remains on the 13-grid, because
``2(k+1/2)=2k+1`` is integral.  This scratch scout searches the two exact
fixed delayed phases 11/24 and 13/24, then asks whether two source-one
THM-2640 packet points can be joined into a physical two-cycle.
"""

from bisect import bisect_right
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_cross_time_target_future_diagonal_thm2616 as core
import lrc14_predecessor_carry_private_root_atlas_thm2640 as private
import lrc14_successor_halfcell_carry_no_go_thm2623 as prior


P = 13
R = P**6
T = core.T
ALPHA = Fraction(1, 2)
FIXED_PHASES = (Fraction(11, 24), Fraction(13, 24))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def floor_fraction(value):
    return value.numerator // value.denominator


def frac(value):
    return value - floor_fraction(value)


def clock(value):
    return floor_fraction(7 * frac(value) + Fraction(1, 2)) % 7


def shallow(value):
    return clock(P * value)


def owner(value):
    return clock(P * P * value)


def half_handoff(value, k):
    return frac(P * value + Fraction(2 * k + 1, 2 * R))


def half_translate(value, k):
    return frac(value + Fraction(2 * k + 1, 2 * R))


def doubled_digit(value, edge):
    rx = R * frac(value)
    carry = floor_fraction(rx) % P
    y = frac(rx)
    digit = floor_fraction(2 * P * y)
    upper = digit // P
    edge_term = int(edge == 0)
    absolute = (2 * carry + upper + edge_term) % P
    return carry, y, digit, upper, absolute


def audit_half_digit_law(value, k, edge):
    before = doubled_digit(value, edge)
    after = doubled_digit(half_translate(value, k), edge)
    carry, _, digit, upper, absolute = before
    carry2, _, digit2, upper2, absolute2 = after
    require(carry2 == (carry + k + upper) % P,
            "C_(2R) half translation carry law changed")
    require(digit2 == (digit + P) % (2 * P) and upper2 == 1 - upper,
            "C_(2R) half translation future-half law changed")
    require(absolute2 == (absolute + 2 * k + 1) % P,
            "C_(2R) absolute doubled-digit law changed")
    return before, after


def prefix_intervals(prefix):
    starts, lengths = prefix[:2]
    return tuple((Fraction(start, T), Fraction(start + length, T))
                 for start, length in zip(starts, lengths))


def strict_interval_member(value, intervals, starts=None):
    if starts is None:
        starts = [left for left, *_ in intervals]
    index = bisect_right(starts, value) - 1
    if index < 0:
        return False
    left, right = intervals[index][:2]
    return Fraction(left) < value < Fraction(right)


def strict_weighted_member(value, pieces, starts):
    index = bisect_right(starts, value) - 1
    if index < 0:
        return False
    left, right, weight = pieces[index]
    return bool(weight) and Fraction(left) < value < Fraction(right)


def interval_slack(value, intervals, starts=None):
    if starts is None:
        starts = [left for left, *_ in intervals]
    index = bisect_right(starts, value) - 1
    require(index >= 0, "slack point lies before support")
    left, right = intervals[index][:2]
    require(Fraction(left) < value < Fraction(right),
            "slack point is not strict in support")
    return min(value - Fraction(left), Fraction(right) - value)


def clock_slack(value, speed):
    phase = frac(speed * value)
    coordinate = 7 * phase + Fraction(1, 2)
    local = frac(coordinate)
    return min(local, 1 - local) / (7 * speed)


def packet_slack(module, rails, rail_starts, present, present_starts,
                 pair_prefixes, candidate):
    """Exact symmetric physical-x radius preserving every tested factor."""
    x = candidate["x"]
    ell = candidate["shallow"]
    h = candidate["h"]
    kappa = candidate["kappa"]
    edge = candidate["edge"]
    root = candidate["root"]
    coordinate = x * T
    rail_margin = interval_slack(
        coordinate, rails[candidate["j"]][3],
        rail_starts[candidate["j"]]
    ) / T
    present_margin = interval_slack(
        coordinate, present[ell, (-h) % P],
        present_starts[ell, (-h) % P]
    ) / T
    y = frac(R * x)
    delayed_margin = interval_slack(
        y, prefix_intervals(pair_prefixes[0][ell][h][kappa])
    ) / R
    digit = 2 * h + kappa
    half_digit_margin = min(
        y - Fraction(digit, 2 * P),
        Fraction(digit + 1, 2 * P) - y,
    ) / R
    carry_margin = min(y, 1 - y) / R
    deep = frac(module.C3 * x) * 182
    low = 14 * root - 13 if edge == 0 else 14 * root
    high = low + 13
    half_tooth_margin = min(deep - low, high - deep) / (182 * module.C3)
    margins = {
        "rail": rail_margin,
        "present": present_margin,
        "delayed": delayed_margin,
        "half_digit": half_digit_margin,
        "carry": carry_margin,
        "half_tooth": half_tooth_margin,
        "shallow_clock": clock_slack(x, P),
        "owner_clock": clock_slack(x, P * P),
    }
    require(all(value > 0 for value in margins.values()),
            "witness packet has a non-strict factor")
    return min(margins.values()), tuple(sorted(margins.items()))


def half_membership(module, value, edge, root):
    phase = frac(module.C3 * value) * 182
    low = 14 * root - 13 if edge == 0 else 14 * root
    high = low + 13
    return Fraction(low) < phase < Fraction(high)


def half_roots(module, value, edge):
    return tuple(root for root in range(1, P)
                 if half_membership(module, value, edge, root))


def main():
    module, _, _, _, rails, present, present_starts = core.build_carrier_data()
    pair_prefixes = private.build_pair_prefixes(module)
    source_rows = private.shard((0, 14))
    require(source_rows[1] == 26
            and all(meta[0] == 1 for meta in source_rows[5]),
            "source-one THM-2640 row bank changed")
    rows = source_rows[6]

    sector_words = prior.sector_words(module)
    sector_starts = [[left for left, _ in word] for word in sector_words]
    delayed_fixed = {}
    for y in FIXED_PHASES:
        require(frac(P * y + ALPHA) == y,
                "chosen half-shift delayed phase stopped being fixed")
        sectors = tuple(sector for sector in range(2)
                        if strict_interval_member(
                            y * T, sector_words[sector], sector_starts[sector]
                        ))
        clocks = tuple(ell for ell in range(7)
                       if strict_interval_member(
                           y, prefix_intervals(
                               pair_prefixes[0][ell]
                               [floor_fraction(P * y)]
                               [floor_fraction(2 * P * y)
                                - 2 * floor_fraction(P * y)]
                           )
                       ))
        require(sectors == (0,) and clocks == (1, 2, 3, 4, 5, 6),
                "half-shift fixed delayed phases changed")
        delayed_fixed[y] = (sectors, clocks)

    rail_starts = [[left for left, _, _ in rail[3]] for rail in rails[:14]]
    candidates = []
    candidate_counts = defaultdict(int)
    for y in FIXED_PHASES:
        h = floor_fraction(P * y)
        kappa = floor_fraction(2 * P * y) - 2 * h
        require((y, h, kappa) in (
            (Fraction(11, 24), 5, 1),
            (Fraction(13, 24), 7, 0),
        ), "half-shift future labels changed")
        first_n = floor_fraction(Fraction(6 * R, P) - y) + 1
        last_n = floor_fraction(Fraction(7 * R, P) - y)
        for n in range(first_n, last_n + 1):
            x = Fraction(n, R) + Fraction(y, R)
            ell = shallow(x)
            if ell == 0:
                continue
            coordinate = x * T
            if not strict_interval_member(
                    coordinate, present[ell, (-h) % P],
                    present_starts[ell, (-h) % P]):
                continue
            rail_indices = [
                j for j in range(14)
                if strict_weighted_member(
                    coordinate, rails[j][3], rail_starts[j]
                )
            ]
            if not rail_indices:
                continue
            carry = n % P
            for j in rail_indices:
                require(rails[j][0] == 1 and rails[j][1] == owner(x),
                        "source-one rail owner typing changed")
                for edge in (0, 1):
                    root = (2 * carry + (2 * h + kappa) // P
                            + (edge == 0)) % P
                    if root == 0 or not half_membership(
                            module, x, edge, root):
                        continue
                    values = rows[j][0][edge][carry][kappa][h]
                    if not private.is_unit(values, root, 26):
                        continue
                    candidate = {
                        "x": x, "N": n, "y": y, "j": j,
                        "shallow": ell, "owner": owner(x),
                        "carry": carry, "h": h, "kappa": kappa,
                        "edge": edge, "root": root,
                    }
                    candidates.append(candidate)
                    candidate_counts[
                        (y, ell, owner(x), edge)
                    ] += 1

    require(candidates, "half-shift fixed delayed fibres have no full local packet")

    by_clock = defaultdict(list)
    for candidate in candidates:
        by_clock[candidate["shallow"], candidate["owner"]].append(candidate)

    witness = None
    pair_checks = 0
    for first in candidates:
        if first["shallow"] == first["owner"]:
            continue
        compatible = by_clock[first["owner"], first["shallow"]]
        for second in compatible:
            # Keeping one delayed fixed phase makes both half-lifts act on the
            # same recurrent delayed carrier.  The physical fibre may change.
            if second["y"] != first["y"]:
                continue
            if second["x"] == first["x"]:
                continue
            a = floor_fraction(P * first["y"] + ALPHA)
            require(a in (6, 7), "half-shift delayed carry branch changed")
            k0 = (second["N"] - P * first["N"] - a) % R
            k1 = (first["N"] - P * second["N"] - a) % R
            require(half_handoff(first["x"], k0) == second["x"]
                    and half_handoff(second["x"], k1) == first["x"],
                    "half-odometer fibre-closing equation failed")
            delta0 = (2 * k0 + 1) % P
            delta1 = (2 * k1 + 1) % P
            if not delta0 or not delta1:
                continue
            d0 = frac(P * first["x"])
            d1 = frac(P * second["x"])
            roots_d0 = half_roots(module, d0, second["edge"])
            roots_d1 = half_roots(module, d1, first["edge"])
            pair_checks += 1
            if not any((root + delta0) % P == second["root"]
                       for root in roots_d0):
                continue
            if not any((root + delta1) % P == first["root"]
                       for root in roots_d1):
                continue
            witness = (first, second, k0, k1, delta0, delta1,
                       roots_d0, roots_d1)
            break
        if witness is not None:
            break

    require(witness is not None,
            "no clock-glued half-odometer two-cycle passed the root law")
    first, second, k0, k1, delta0, delta1, roots_d0, roots_d1 = witness
    doubled_laws = (
        audit_half_digit_law(frac(P * first["x"]), k0, second["edge"]),
        audit_half_digit_law(frac(P * second["x"]), k1, first["edge"]),
    )
    direct_target_carries = (
        (k0 + floor_fraction(P * first["y"] + ALPHA)) % P,
        (k1 + floor_fraction(P * second["y"] + ALPHA)) % P,
    )
    require(direct_target_carries
            == (second["carry"], first["carry"]),
            "factorized D-then-half-translation carry law changed")
    first_slack, first_margins = packet_slack(
        module, rails, rail_starts, present, present_starts,
        pair_prefixes, first
    )
    second_slack, second_margins = packet_slack(
        module, rails, rail_starts, present, present_starts,
        pair_prefixes, second
    )
    three_state_radius = min(first_slack / P**2, second_slack / P)
    require(three_state_radius > 0,
            "strict half-odometer witness lost its open cylinder")
    require(len(candidates) == 332668,
            "half-odometer full local candidate census changed")
    require(first == {
        "x": Fraction(55232507, 115843416),
        "N": 2301354, "y": Fraction(11, 24), "j": 8,
        "shallow": 1, "owner": 4, "carry": 3,
        "h": 5, "kappa": 1, "edge": 0, "root": 7,
    }, "first half-odometer witness changed")
    require(second == {
        "x": Fraction(58313459, 115843416),
        "N": 2429727, "y": Fraction(11, 24), "j": 3,
        "shallow": 4, "owner": 1, "carry": 1,
        "h": 5, "kappa": 1, "edge": 0, "root": 3,
    }, "second half-odometer witness changed")
    require((pair_checks, k0, k1, delta0, delta1,
             roots_d0, roots_d1)
            == (1, 1472973, 4502560, 4, 8, (12,), (12,)),
            "half-odometer handoff/root law changed")
    require(tuple((after[-1] - before[-1]) % P
                  for before, after in doubled_laws) == (4, 8),
            "doubled-digit and geometric root steps disagree")
    phase_counts = tuple(
        (y, sum(count for (phase, _, _, _), count
                in candidate_counts.items() if phase == y))
        for y in FIXED_PHASES
    )
    nonconstant_counts = tuple(
        (y, sum(count for (phase, left, right, _), count
                in candidate_counts.items()
                if phase == y and left != right))
        for y in FIXED_PHASES
    )
    clock_edge_support = tuple(sorted(set(
        (phase, left, right)
        for phase, left, right, _ in candidate_counts
        if left != right
    )))
    require(len(clock_edge_support) == 51
            and all((y, 1, 4) in clock_edge_support
                    and (y, 4, 1) in clock_edge_support
                    for y in FIXED_PHASES),
            "half-odometer nonconstant clock-edge support changed")

    print("LRC14 HALF-ODOMETER CHANGED DELAYED CARRIER EXACT SCRATCH SCOUT")
    print(f"scales=(R={R},alpha={ALPHA})")
    print(f"fixed_delayed_phases={tuple((y, delayed_fixed[y]) for y in FIXED_PHASES)}")
    print(f"candidate_count={len(candidates)}")
    print(f"candidate_phase_counts={phase_counts}; "
          f"nonconstant_phase_counts={nonconstant_counts}")
    print(f"nonconstant_clock_edge_support_count={len(clock_edge_support)}; "
          "reciprocal_1_4_present_at_both_phases=True")
    print(f"pair_checks={pair_checks}")
    print(f"first={first}")
    print(f"second={second}")
    print(f"half_lifts=(k0={k0},k1={k1}); root_steps=({delta0},{delta1})")
    print(f"intermediate_D_roots=(first={roots_d0},second={roots_d1})")
    print(f"doubled_digit_laws={doubled_laws}")
    print(f"direct_target_carries={direct_target_carries}; "
          f"factorization=A_half=tau_odd_after_D")
    print(f"packet_slacks=(first={first_slack},second={second_slack}); "
          f"three_state_open_radius={three_state_radius}; "
          f"open_length={2*three_state_radius}")
    print(f"binding_margin_ledgers=(first={first_margins},second={second_margins})")
    print("delayed_base_cycle=y -> {13y+1/2}=y; inherited sector0 and clocks1..6 recur")
    print("typing_boundary=C_(2R)=C_2 x C_R; H(x)=x+1/2 is central and commutes with D, every odd lift is H times an integral odometer lift, and the filtration retains residual C2; existing (carry,h,kappa,edge) data type the pure half-translation root step through a_edge=2carry+floor((2h+kappa)/13)+1_{edge=left}")
    print("verdict=POSITIVE SCOUT: a changed delayed carrier has a full local rail/present/delayed/private-unit two-cycle; semantic cospan, endpoint, and LRC consequence remain unproved")


if __name__ == "__main__":
    main()
