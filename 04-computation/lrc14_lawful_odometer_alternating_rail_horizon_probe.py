#!/usr/bin/env python3
"""Exact finite-horizon THM-2584 rail fibre for lawful odometer lifts.

Set ``P=13``, ``R=P^6``, ``S=P^5``, and ``K=S+1``.  The two lawful affine
handoffs

    T_-(x)={P*x-K/R},       T_+(x)={P*x+K/R}

exchange

    x_+=1/2+K/(14R),        x_-=1/2-K/(14R).

Their intrinsic stored clock edges are respectively ``4->3`` and ``3->4``;
each owner is the following event's shallow clock.  The THM-2584 route-two
rails with keys ``(6,1,3,0)`` and ``(6,1,4,12)`` contain ``x_+`` and ``x_-``
strictly.  Their local support intervals are exchanged by reflection.

If the initial point is ``x_++e``, its nth alternating state is
``x_+ + P^n e`` for even n and ``x_- + P^n e`` for odd n.  The local rail
margins nest strictly, so the exact H-event cylinder is the final pulled-back
rail interval.  It has length

    1/(7*13^(H+1))

for every finite ``H>=1``.  The intersection over all H is the repelling
cycle point only; no uniform positive infinite-horizon interval is claimed.

The companion also constructs the literal three-event weighted fibre product
on the common ``13^2*T`` grid as an independent check.  It transports the
THM-2584 b-word depth-five product profile, source displacement, arrival,
owner, and absolute deep-root labels, plus an explicit shallow-clock cut.
It does not yet transport THM-2640's present packet, delayed sector/clock
word, predecessor carry, private edge/root, or primitive unit.
"""

from fractions import Fraction

import lrc14_alternate_arrival_physical_rail_handoff as rail
import lrc14_dilation_reversed_clock_fibre_product_probe as two


P = 13
Q = 7
R = P**6
S = P**5
K = S + 1
LIFTS = (-K, K)
PLUS_KEY = (6, 1, 3, 0)
MINUS_KEY = (6, 1, 4, 12)
KEYS = (PLUS_KEY, MINUS_KEY)
SHALLOW_CLOCKS = (4, 3)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def circle(value):
    return value % 1


def handoff(value, lift):
    return circle(P * value + Fraction(lift, R))


def clock(value):
    return int((Q * circle(value) + Fraction(1, 2)) // 1) % Q


def shallow(value):
    return clock(P * value)


def owner(value):
    return clock(P * P * value)


def strict_piece_at(pieces, point_on_grid):
    hits = tuple(piece for piece in pieces
                 if piece[0] < point_on_grid < piece[1])
    require(len(hits) == 1, "cycle point is not in one strict rail piece")
    return hits[0]


def strict_interval_at(intervals, point_on_grid):
    hits = tuple(interval for interval in intervals
                 if interval[0] < point_on_grid < interval[1])
    require(len(hits) == 1,
            "cycle point is not in one strict shallow interval")
    return hits[0]


def iterate_translation_coefficients(horizon):
    """c_i for x_i={P^i*x+c_i/R} under alternating lifts."""
    coefficients = [0]
    for step in range(horizon - 1):
        coefficients.append(P * coefficients[-1] + LIFTS[step % 2])
    return tuple(coefficients)


def pullback_iterate_weighted(pieces, step, horizon, coefficient, grid):
    """Pull one iterate to the common ``P^(H-1)*grid`` integer grid."""
    degree = P**step
    common_scale = P**(horizon - 1 - step)
    require(grid % R == 0, "THM-2584 grid lost the odometer lattice")
    shift = (coefficient * (grid // R)) % grid
    out = []
    for left, right, weight in pieces:
        for branch in range(degree + 1):
            lo = max(0, left - shift + branch * grid)
            hi = min(degree * grid, right - shift + branch * grid)
            if lo < hi:
                out.append((common_scale * lo, common_scale * hi, weight))
    return sorted(out)


def weighted_grid_mass(pieces):
    return sum((right - left) * weight
               for left, right, weight in pieces)


def main():
    bank = rail.build_rail_bank()
    shallow_cells = rail.build_shallow_cells()
    grid = rail.old.T
    require(all(key in bank for key in KEYS),
            "one alternating THM-2584 rail disappeared")

    amplitude = Fraction(K, (P + 1) * R)
    plus = Fraction(1, 2) + amplitude
    minus = Fraction(1, 2) - amplitude
    require(handoff(plus, -K) == minus
            and handoff(minus, K) == plus,
            "lawful affine two-cycle changed")
    require(((shallow(plus), owner(plus)),
             (shallow(minus), owner(minus))) == ((4, 3), (3, 4)),
            "intrinsic stored clock edges changed")
    require(owner(plus) == shallow(minus)
            and owner(minus) == shallow(plus),
            "owner-to-following-shallow interface stopped gluing")
    require(K % P == 1 and tuple((2 * lift) % P for lift in LIFTS)
            == (11, 2),
            "THM-2657 unit/quotient labels changed")

    centres = (plus, minus)
    local_rails = tuple(
        strict_piece_at(bank[key], centre * grid)
        for key, centre in zip(KEYS, centres)
    )
    local_shallows = tuple(
        strict_interval_at(shallow_cells[label], centre * grid)
        for label, centre in zip(SHALLOW_CLOCKS, centres)
    )
    local_rail_intervals = tuple(
        (Fraction(left, grid), Fraction(right, grid))
        for left, right, _ in local_rails
    )
    local_shallow_intervals = tuple(
        (Fraction(left, grid), Fraction(right, grid))
        for left, right in local_shallows
    )
    require(local_rail_intervals == (
        (Fraction(1195, 2366), Fraction(171, 338)),
        (Fraction(167, 338), Fraction(1171, 2366)),
    ), "alternating local rail intervals changed")
    require(local_shallow_intervals == (
        (Fraction(1, 2), Fraction(93, 182)),
        (Fraction(89, 182), Fraction(1, 2)),
    ), "alternating local shallow intervals changed")
    require(all(sleft <= rleft < rright <= sright
                for (rleft, rright), (sleft, sright)
                in zip(local_rail_intervals, local_shallow_intervals)),
            "a local rail stopped forcing its shallow clock")
    weights = tuple(piece[2] for piece in local_rails)
    require(weights == (27_580_222_516, 27_580_222_516),
            "alternating local rail weights changed")

    # The two margins differ by the single refined endpoint atom 1/(7R).
    outer_margin = plus - local_rail_intervals[0][0]
    inner_margin = local_rail_intervals[0][1] - plus
    require((outer_margin, inner_margin) == (
        Fraction(14_281, 7 * R), Fraction(2_040, R)
    ), "alternating rail margins changed")
    require((minus - local_rail_intervals[1][0],
             local_rail_intervals[1][1] - minus)
            == (inner_margin, outer_margin),
            "reflected rail margins changed")
    require(outer_margin + inner_margin == Fraction(1, 7 * P**2)
            and P * inner_margin > outer_margin,
            "nested finite-horizon margin law changed")

    horizon_rows = []
    for horizon in range(1, 129):
        constraints = []
        for step in range(horizon):
            centre = centres[step % 2]
            left, right = local_rail_intervals[step % 2]
            constraints.append((
                (left - centre) / P**step,
                (right - centre) / P**step,
            ))
        perturbation = (
            max(interval[0] for interval in constraints),
            min(interval[1] for interval in constraints),
        )
        if (horizon - 1) % 2 == 0:
            expected_perturbation = (
                -outer_margin / P**(horizon - 1),
                inner_margin / P**(horizon - 1),
            )
        else:
            expected_perturbation = (
                -inner_margin / P**(horizon - 1),
                outer_margin / P**(horizon - 1),
            )
        require(perturbation == expected_perturbation,
                "last local rail stopped controlling the horizon cylinder")
        interval = (plus + perturbation[0], plus + perturbation[1])
        length = interval[1] - interval[0]
        require(length == Fraction(1, 7 * P**(horizon + 1)),
                "finite-horizon cylinder length changed")
        horizon_rows.append((horizon, interval, length))

    # Independent direct common-grid fibre product for H=3.
    direct_horizon = 3
    coefficients = iterate_translation_coefficients(direct_horizon)
    direct = None
    common_grid = P**(direct_horizon - 1) * grid
    for step in range(direct_horizon):
        pulled_rail = pullback_iterate_weighted(
            bank[KEYS[step % 2]], step, direct_horizon,
            coefficients[step], grid,
        )
        direct = (pulled_rail if direct is None
                  else two.intersect_weighted(direct, pulled_rail))
        pulled_shallow = pullback_iterate_weighted(
            [(left, right, 1)
             for left, right in shallow_cells[SHALLOW_CLOCKS[step % 2]]],
            step, direct_horizon, coefficients[step], grid,
        )
        direct = two.intersect_weighted(direct, pulled_shallow)
    require(len(direct) == 1,
            "direct three-event alternating fibre acquired extra components")
    expected_h3 = horizon_rows[direct_horizon - 1][1]
    require((Fraction(direct[0][0], common_grid),
             Fraction(direct[0][1], common_grid), direct[0][2]) == (
        expected_h3[0], expected_h3[1], weights[0]**direct_horizon,
    ), "direct common-grid route disagrees with the horizon formula")
    direct_numerator = weighted_grid_mass(direct)

    sampled = tuple(
        (horizon, interval, length)
        for horizon, interval, length in horizon_rows
        if horizon in (1, 2, 3, 7, 12)
    )
    print("LRC14 lawful odometer alternating THM2584 rail-horizon probe")
    print(f"scales=P:{P} R:{R} S:{S} lift_magnitude:{K} quotient_steps={(11,2)}")
    print(f"cycle=(plus:{plus},minus:{minus}) amplitude={amplitude}")
    print("cycle_stored_edges=((4,3),(3,4)) interfaces_glue=True")
    print(f"selected_keys={KEYS} selected_source=1 local_weights={weights}")
    print(f"local_rail_intervals={local_rail_intervals}")
    print(f"local_shallow_intervals={local_shallow_intervals}")
    print(f"one_sided_margins=(outer:{outer_margin},inner:{inner_margin})")
    print("arbitrary_horizon_law=H-event component length 1/(7*13^(H+1)); last pulled-back rail is binding")
    print(f"sample_horizon_components={sampled}")
    print(f"direct_H3_coefficients={coefficients} components={len(direct)} weighted_grid_numerator={direct_numerator}")
    print("independent_routes=128 exact local inverse-cylinder intersections plus literal H3 weighted common-grid fibre agree")
    print("infinite_horizon_boundary=repelling cycle singleton only; no uniform positive open interval")
    print("transported=THM2584 b-word depth-five product profile,source displacement,arrival,owner,deep label; explicit shallow clock")
    print("not_transported=THM2640 present factor,delayed sector/clock word,predecessor carry,private edge/root,primitive unit,semantic transition,row exclusion,LRC14")


if __name__ == "__main__":
    main()
