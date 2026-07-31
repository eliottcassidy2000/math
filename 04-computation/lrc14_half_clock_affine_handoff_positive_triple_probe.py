#!/usr/bin/env python3
"""Exact positive three-rail escape for the half-clock affine handoff.

For ``q>=3`` put

    p=2*q-1,  delta=1/(4*q),  beta=1/(2*q),
    F(x)={p*x+beta},  D(x)={p*x},
    x*=delta-(p-3)*delta/p**3.

Then ``F^2=D^2`` because ``(p+1)*beta=1``.  Direct algebra gives

    F(x*) = 1/2+delta-(p-3)*delta/p**2,
    F^2(x*) = 3*delta/p.

Thus the orbit visits the low endpoint tooth, the central tooth, and the
low endpoint tooth.  With the intrinsic event clock

    shallow(u)=c_q(Du),  owner(u)=c_q(D^2u),
    c_q(v)=floor(q*{v}+1/2) mod q,

its first three stored edges begin

    floor(q/2) -> 0 -> 1.

At q=7, the exact THM-2584 route-two rail keys

    (arrival,source,owner,deep)
      (0,1,0,0), (6,1,1,0), (0,1,3,0)

support the three events at ``x,F(x),F^2(x)``.  This companion forms their
literal weighted fibre product, including shallow clocks 3,0,1, on the
common ``13^2*T`` grid.  It proves a positive open component and independently
recovers that component by pulling back the three local rational intervals.

This is a boundary-scale affine handoff scout.  The translation ``1/14`` is
not a THM-2657 lift ``k/13^6`` and does not preserve its predecessor-carry,
root, or delayed-prefix packet.  No semantic transition, unit, row exclusion,
or LRC(14) conclusion is asserted.
"""

from fractions import Fraction

import lrc14_alternate_arrival_physical_rail_handoff as rail
import lrc14_dilation_reversed_clock_fibre_product_probe as two


P = 13
Q = 7
BETA = Fraction(1, 2 * Q)
KEYS = ((0, 1, 0, 0), (6, 1, 1, 0), (0, 1, 3, 0))
SHALLOW_CLOCKS = (3, 0, 1)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def circle(value):
    return value % 1


def dilation(q, value):
    return circle((2 * q - 1) * value)


def affine_handoff(q, value):
    return circle((2 * q - 1) * value + Fraction(1, 2 * q))


def event_clock(q, value):
    return int((q * circle(value) + Fraction(1, 2)) // 1) % q


def shallow_clock(q, value):
    return event_clock(q, dilation(q, value))


def owner_clock(q, value):
    return event_clock(q, dilation(q, dilation(q, value)))


def general_orbit_audit():
    """Hostile finite replay of identities proved algebraically above."""
    checks = 0
    for q in tuple(range(3, 129)) + (257, 509, 1021):
        p = 2 * q - 1
        delta = Fraction(1, 4 * q)
        beta = Fraction(1, 2 * q)
        x = delta - Fraction(p - 3, p**3) * delta
        y = affine_handoff(q, x)
        z = affine_handoff(q, y)
        require(
            y == Fraction(1, 2) + delta
                 - Fraction(p - 3, p**2) * delta,
            "affine first iterate formula failed",
        )
        require(z == Fraction(3, p) * delta,
                "affine second iterate formula failed")
        require(z == dilation(q, dilation(q, x)),
                "F^2=D^2 failed")
        require(0 < x < delta,
                "source left the low endpoint tooth")
        require(Fraction(1, 2) - delta < y
                < Fraction(1, 2) + delta,
                "first image left the central tooth")
        require(0 < z < delta,
                "second image left the low endpoint tooth")
        require(
            (shallow_clock(q, x), owner_clock(q, x),
             shallow_clock(q, y), owner_clock(q, y),
             shallow_clock(q, z))
            == (q // 2, 0, 0, 1, 1),
            "intrinsic half-clock chain changed",
        )
        checks += 1
    return checks


def pullback_affine_weighted(pieces, grid):
    """Pull ``pieces(Fx)`` to the common ``P^2*grid`` integer grid."""
    require(grid % (2 * Q) == 0,
            "base grid does not resolve the half-clock translation")
    shift = grid // (2 * Q)
    refined = []
    # On X=P*grid*x, F(x)={X/grid+shift/grid}.  The shifted coordinate can
    # cross branches 0,...,P, after which we scale once more to P^2*grid.
    for left, right, weight in pieces:
        for branch in range(P + 1):
            lo = max(0, left - shift + branch * grid)
            hi = min(P * grid, right - shift + branch * grid)
            if lo < hi:
                refined.append((P * lo, P * hi, weight))
    return sorted(refined)


def pullback_affine_intervals(intervals, grid):
    return [
        (left, right)
        for left, right, _ in pullback_affine_weighted(
            [(left, right, 1) for left, right in intervals], grid
        )
    ]


def weighted_grid_mass(pieces):
    return sum((right - left) * weight
               for left, right, weight in pieces)


def strict_piece_at(pieces, point_on_grid):
    hits = tuple(piece for piece in pieces
                 if piece[0] < point_on_grid < piece[1])
    require(len(hits) == 1, "control point is not in one strict rail piece")
    return hits[0]


def strict_interval_at(intervals, point_on_grid):
    hits = tuple(interval for interval in intervals
                 if interval[0] < point_on_grid < interval[1])
    require(len(hits) == 1,
            "control point is not in one strict shallow interval")
    return hits[0]


def main():
    formula_checks = general_orbit_audit()
    bank = rail.build_rail_bank()
    shallow = rail.build_shallow_cells()
    grid = rail.old.T
    grid2 = P * P * grid
    require(P == 2 * Q - 1 and grid % (2 * Q) == 0,
            "canonical 13/7 half-clock scales changed")
    require(all(key in bank for key in KEYS),
            "one selected THM-2584 rail key disappeared")

    delta = Fraction(1, 4 * Q)
    x = delta - Fraction(P - 3, P**3) * delta
    orbit = (x, affine_handoff(Q, x),
             affine_handoff(Q, affine_handoff(Q, x)))
    require(orbit == (
        Fraction(2187, 61516),
        Fraction(2525, 4732),
        Fraction(3, 364),
    ), "q=7 affine orbit changed")
    require(tuple(shallow_clock(Q, value) for value in orbit)
            == SHALLOW_CLOCKS,
            "q=7 shallow clock word changed")
    require((owner_clock(Q, orbit[0]), owner_clock(Q, orbit[1]))
            == SHALLOW_CLOCKS[1:],
            "q=7 owner/shallow covariance failed")

    # Direct common-grid route through the complete weighted profiles.
    current = two.scale_weighted(bank[KEYS[0]], P * P)
    following = pullback_affine_weighted(bank[KEYS[1]], grid)
    third = rail.pullback_d2_weighted(bank[KEYS[2]])
    require(weighted_grid_mass(following)
            == P * P * weighted_grid_mass(bank[KEYS[1]]),
            "affine pullback did not preserve normalized weighted mass")
    joint = two.intersect_weighted(
        two.intersect_weighted(current, following), third
    )

    shallow_current = [
        (P * P * left, P * P * right)
        for left, right in shallow[SHALLOW_CLOCKS[0]]
    ]
    shallow_following = pullback_affine_intervals(
        shallow[SHALLOW_CLOCKS[1]], grid
    )
    shallow_third = rail.pullback_d2_intervals(
        shallow[SHALLOW_CLOCKS[2]]
    )
    for intervals in (shallow_current, shallow_following, shallow_third):
        joint = rail.intersect_weighted_union(joint, intervals)
    require(joint, "the affine three-rail clock fibre product is empty")

    x_grid = x * grid2
    require(x_grid.denominator == 1,
            "canonical point left the common integer grid")
    witness = strict_piece_at(joint, x_grid.numerator)
    expected_left = Fraction(14215, 399854)
    expected_right = Fraction(2031, 57122)
    expected_weight = (
        27_397_688_484 * 27_581_135_604 * 27_396_781_668
    )
    require((Fraction(witness[0], grid2),
             Fraction(witness[1], grid2), witness[2])
            == (expected_left, expected_right, expected_weight),
            "canonical positive component changed")
    require(expected_right - expected_left == Fraction(1, 199927),
            "canonical positive component length changed")

    # Independent local route: find one strict rail and shallow component at
    # each orbit point, then invert the branch used by F and F^2 rationally.
    local_rails = tuple(
        strict_piece_at(bank[key], value * grid)
        for key, value in zip(KEYS, orbit)
    )
    local_shallows = tuple(
        strict_interval_at(shallow[clock], value * grid)
        for clock, value in zip(SHALLOW_CLOCKS, orbit)
    )
    require(tuple((Fraction(a, grid), Fraction(b, grid))
                  for a, b, _ in local_rails) == (
        (Fraction(83, 2366), Fraction(1, 28)),
        (Fraction(97, 182), Fraction(1263, 2366)),
        (Fraction(19, 2366), Fraction(3, 338)),
    ), "local rail intervals changed")
    require(tuple((Fraction(a, grid), Fraction(b, grid))
                  for a, b in local_shallows) == (
        (Fraction(5, 182), Fraction(1, 26)),
        (Fraction(97, 182), Fraction(99, 182)),
        (Fraction(1, 182), Fraction(3, 182)),
    ), "local shallow intervals changed")
    local_constraints = (
        tuple(Fraction(endpoint, grid) for endpoint in local_rails[0][:2]),
        tuple((Fraction(endpoint, grid) - BETA) / P
              for endpoint in local_rails[1][:2]),
        tuple((Fraction(6) + Fraction(endpoint, grid)) / (P * P)
              for endpoint in local_rails[2][:2]),
    )
    local_intersection = (
        max(interval[0] for interval in local_constraints),
        min(interval[1] for interval in local_constraints),
    )
    require(local_intersection == (expected_left, expected_right),
            "independent rational pullback missed the direct component")
    require(all(a < x < b for a, b in local_constraints),
            "canonical point left one local pullback interval")

    # Positive length is the actual support certificate.  The product weight
    # is recorded only as a reproducibility sidecar, not a canonical mass.
    component_grid_length = witness[1] - witness[0]
    component_numerator = component_grid_length * witness[2]
    owner_word = tuple(key[2] for key in KEYS)
    stored_edges = tuple(zip(SHALLOW_CLOCKS, owner_word))
    require(stored_edges == ((3, 0), (0, 1), (1, 3))
            and all(left != right for left, right in stored_edges),
            "selected stored clock word stopped being nonconstant")

    print("LRC14 half-clock affine-handoff positive triple probe")
    print(f"general_q_formula_controls={formula_checks} q_range=3..128_plus_257_509_1021")
    print("identity=F(x)={p*x+1/(2q)} has F^2=D_p^2 because (p+1)/(2q)=1")
    print(f"q=7_orbit=(x,Fx,F2x)={orbit}")
    print("tooth_word=endpoint-low,central,endpoint-low")
    print(f"intrinsic_shallow_word={SHALLOW_CLOCKS} owner_prefix={owner_word[:2]}")
    print(f"thm2584_rail_keys={KEYS} stored_clock_edges={stored_edges}")
    print(f"positive_component=({expected_left},{expected_right}) length={expected_right-expected_left}")
    print(f"component_grid_length={component_grid_length} component_weight={witness[2]}")
    print(f"component_weighted_grid_numerator={component_numerator}")
    print(f"all_selected_fibre_components={len(joint)} total_weighted_grid_mass={weighted_grid_mass(joint)}")
    print(f"local_inverse_constraints={local_constraints}")
    print("independent_routes=full weighted common-grid pullback and local rational branch intersection agree")
    print("scope=positive THM2584 rail-clock fibre product for a boundary-scale affine handoff only")
    print("not_lawful_odometer=1/14 is not k/13^6; carry,root,delayed-prefix,unit,semantic transition,row exclusion,and LRC14 remain untested")


if __name__ == "__main__":
    main()
