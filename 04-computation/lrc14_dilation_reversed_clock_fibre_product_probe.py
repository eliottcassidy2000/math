#!/usr/bin/env python3
"""Exact two-edge LRC(14) fibre products under the dilation handoff.

The common-x matrices downstream of THM-2616/2629/2630 are indexed by an
owner (deep) clock ``ell4`` and a present (shallow) clock ``ell5``.  Literal
same-x composition of cells with different owner clocks is vacuous, since
the owner-clock cells partition the circle.  The naturally typed handoff is
instead

    D(x) = {13*x}.

It has two exact covariance laws:

    shallow_ell(Dx) = owner_ell(x),
    j(Dx) = h(x).

Consequently a D-compatible candidate chronology traverses the printed clock
arrows backwards.  For three clock labels ``(ell0,ell1,ell2)`` this script
forms the genuine support fibre product

    E_(s0; owner=ell1, shallow=ell0; j,h)(x)
      AND
    E_(s1; owner=ell2, shallow=ell1; h,k)(D x).

The two source labels are allowed to drift independently.  Each E is the
exact nonnegative sharp-graph atom before integration: a THM-2616 middle
rail, present packet, sharp deep root, tooth half, predecessor digit, carry
half, and delayed guard-safe word digit.  Rail profile weights are multiplied
only as a positivity-preserving support witness; the printed magnitudes are
not asserted to be canonical transition masses.

The pullback by D is represented exactly on the refined grid 13*T.  Its two
delayed-word factors reduce to one set in y={13^6*x}:

    Q_A(y) AND Q_B({13*y}).

No floating point arithmetic, Bockstein unit test, endpoint quotient, or
Perron-sheet choice occurs here.  The finite atlas is a typed feasibility
probe, not a chronology theorem and not a proof of LRC(14).
"""

from bisect import bisect_right
from collections import Counter

import lrc14_cross_time_target_future_diagonal_thm2616 as cross
import lrc14_fallback_rail_digit_diagonal_pullback_thm2592 as old


P = 13
Q7 = 7
INV2 = 7
PREDECESSOR_SPEED = P**5
SUCCESSOR_SPEED = P**6
FUTURE_DIGITS = tuple(range(1, 12))

EXPECTED_GENERIC_CLOCK_SUPPORT = (
    (3, 0, 3, 6),
    (3, 0, 4, 10),
    (3, 0, 5, 10),
    (3, 0, 6, 12),
    (3, 1, 0, 20),
    (4, 0, 1, 14),
    (4, 0, 2, 14),
    (4, 0, 3, 14),
    (4, 0, 4, 12),
    (4, 6, 0, 18),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def merge_adjacent(intervals):
    """Normalize a half-open interval union without changing its set."""
    out = []
    for left, right in sorted(intervals):
        if out and out[-1][1] == left:
            out[-1] = (out[-1][0], right)
        else:
            out.append((left, right))
    return out


def prefix_intervals(prefix):
    starts, lengths, _ = prefix
    return [(left, left + length) for left, length in zip(starts, lengths)]


def preimage_times_13(intervals, grid):
    """Pull a set back by y -> {13y}, retaining the original grid."""
    out = []
    for left, right in intervals:
        require(left % P == right % P == 0,
                "delayed-word endpoint left the 13-divisible grid")
        for branch in range(P):
            out.append(((left + branch * grid) // P,
                        (right + branch * grid) // P))
    return sorted(out)


def scale_weighted(pieces, factor):
    return [(factor * left, factor * right, weight)
            for left, right, weight in pieces]


def pullback_dilation_weighted(pieces, grid):
    """Represent x -> {13x} pullback on the refined grid 13*grid."""
    return sorted(
        (left + branch * grid, right + branch * grid, weight)
        for left, right, weight in pieces
        for branch in range(P)
    )


def intersect_weighted(left, right):
    """Intersect two sorted weighted step supports, multiplying weights."""
    out = []
    i = j = 0
    while i < len(left) and j < len(right):
        a = max(left[i][0], right[j][0])
        b = min(left[i][1], right[j][1])
        if a < b:
            out.append((a, b, left[i][2] * right[j][2]))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return out


def phi_at(x, starts, lengths, prefix_sums):
    index = bisect_right(starts, x) - 1
    if index < 0:
        return 0
    return prefix_sums[index] + min(x - starts[index], lengths[index])


def delayed_weighted_numerator(pieces, prefix, speed, grid):
    """Exact numerator for a weighted E(x) AND Q({speed*x}) overlap."""
    starts, lengths, prefix_sums = prefix
    q_length = prefix_sums[-1]
    weighted_length = 0
    residue_accumulator = 0
    prefix_accumulator = 0
    reduced_speed = speed % grid
    for left, right, weight in pieces:
        rleft = left * reduced_speed % grid
        rright = right * reduced_speed % grid
        weighted_length += weight * (right - left)
        residue_accumulator += weight * (rright - rleft)
        prefix_accumulator += weight * (
            phi_at(rright, starts, lengths, prefix_sums)
            - phi_at(rleft, starts, lengths, prefix_sums)
        )
    floor_numerator = speed * weighted_length - residue_accumulator
    require(floor_numerator % grid == 0,
            "generalized delayed floor count is not integral")
    value = q_length * (floor_numerator // grid) + prefix_accumulator
    require(value >= 0, "negative generalized delayed overlap")
    return value


def build_atom(module, prefixes, present, present_starts, rail,
               future_clock, h, epsilon, kappa):
    """Build one unintegrated sharp-graph atom and its one-edge mass."""
    _, _, _, rail_pieces = rail
    r = (-h - 1) % P
    require(r != 0, "sharp graph entered the absent deep root")
    pieces = old.intersect_weighted_union(
        rail_pieces,
        present[future_clock, (-h) % P],
        present_starts[future_clock, (-h) % P],
    )
    left, right = (
        (14 * r, 14 * r + 13)
        if epsilon == 0
        else (14 * r - 13, 14 * r)
    )
    pieces = old.intersect_weighted_comb(
        pieces, module.C3, 182, left, right
    )
    j = INV2 * (r - epsilon - kappa) % P
    pieces = old.intersect_weighted_comb(
        pieces, PREDECESSOR_SPEED, P, j, j + 1
    )
    pieces = old.intersect_weighted_comb(
        pieces, SUCCESSOR_SPEED, 2, kappa, kappa + 1
    )
    value = old.delayed_weighted_numerator(
        pieces, prefixes[future_clock][h]
    )
    return {
        "future": future_clock,
        "j": j,
        "h": h,
        "epsilon": epsilon,
        "kappa": kappa,
        "pieces": pieces,
        "value": value,
    }


def main():
    (module, prefixes, _, _, rails,
     present, present_starts) = cross.build_carrier_data()
    grid = old.T
    refined_grid = P * grid
    require(module.C1 == P and module.C3 == 2 * PREDECESSOR_SPEED,
            "inherited LRC scale constants changed")

    # The literal same-x product is the wrong type: distinct owner cells are
    # disjoint.  D instead identifies the current owner with the next shallow
    # clock.  Check that identity as exact interval unions, including ell=0's
    # wraparound presentation.
    owner_cells = old.base.clock_cells(P, Q7, grid, P * P)
    direct_owner_disjoint_checks = 0
    dilation_clock_checks = 0
    for ell in range(Q7):
        for other in range(Q7):
            if ell == other:
                continue
            require(not old.intersect_sorted(owner_cells[ell],
                                             owner_cells[other]),
                    "distinct owner-clock cells overlap")
            direct_owner_disjoint_checks += 1
        shallow = module.make_comb(
            module.C1, 182, 26 * ell - 13, 26 * ell + 13
        )
        scaled_owner = merge_adjacent(
            (P * left, P * right) for left, right in owner_cells[ell]
        )
        pulled_shallow = merge_adjacent(
            (left + branch * grid, right + branch * grid)
            for left, right in shallow
            for branch in range(P)
        )
        require(scaled_owner == pulled_shallow,
                "D did not transport an owner clock to the shallow clock")
        dilation_clock_checks += 1

    rails_by_cell = {}
    for index, rail in enumerate(rails):
        rails_by_cell.setdefault(rail[:2], []).append((index, rail))

    atom_bank = {}

    def positive_atoms(s, owner_clock, future_clock):
        key = (s, owner_clock, future_clock)
        if key in atom_bank:
            return atom_bank[key]
        atoms = []
        for rail_index, rail in rails_by_cell[s, owner_clock]:
            for h in FUTURE_DIGITS:
                for epsilon in (0, 1):
                    for kappa in (0, 1):
                        atom = build_atom(
                            module, prefixes, present, present_starts,
                            rail, future_clock, h, epsilon, kappa,
                        )
                        if atom["value"]:
                            atom["rail_index"] = rail_index
                            atoms.append(atom)
        atom_bank[key] = atoms
        return atoms

    def atom_relation(s, owner_clock, future_clock):
        """Boolean ``j -> h`` relation after forgetting rail/binary labels."""
        rows = [0] * P
        for atom in positive_atoms(s, owner_clock, future_clock):
            rows[atom["j"]] |= 1 << atom["h"]
        return tuple(rows)

    def boolean_product(left, right):
        """Formal label product; this deliberately performs no D-pullback."""
        result = []
        for middle_mask in left:
            row = 0
            while middle_mask:
                bit = middle_mask & -middle_mask
                row |= right[bit.bit_length() - 1]
                middle_mask -= bit
            result.append(row)
        return tuple(result)

    def relation_support(relation):
        return sum(row.bit_count() for row in relation)

    joint_prefix_cache = {}

    def fibre_product_value(current, following):
        """Product-weight support numerator for E(x) AND E'(D x)."""
        base = intersect_weighted(
            scale_weighted(current["pieces"], P),
            pullback_dilation_weighted(following["pieces"], grid),
        )
        key = (
            current["future"], current["h"],
            following["future"], following["h"],
        )
        if key not in joint_prefix_cache:
            current_q = prefix_intervals(
                prefixes[current["future"]][current["h"]]
            )
            following_q_pullback = preimage_times_13(
                prefix_intervals(
                    prefixes[following["future"]][following["h"]]
                ),
                grid,
            )
            joint_q = old.intersect_sorted(
                current_q, following_q_pullback
            )
            joint_prefix_cache[key] = module.make_prefix([
                (P * left, P * right) for left, right in joint_q
            ])
        return delayed_weighted_numerator(
            base,
            joint_prefix_cache[key],
            SUCCESSOR_SPEED,
            refined_grid,
        )

    # Independent formula controls: scaling both the base and Q grids by 13
    # multiplies the old numerator by 13 and changes no positivity decision.
    controls = 0
    for atom in positive_atoms(1, 1, 3):
        scaled_prefix = module.make_prefix([
            (P * left, P * right)
            for left, right in prefix_intervals(
                prefixes[atom["future"]][atom["h"]]
            )
        ])
        scaled_value = delayed_weighted_numerator(
            scale_weighted(atom["pieces"], P),
            scaled_prefix,
            SUCCESSOR_SPEED,
            refined_grid,
        )
        require(scaled_value == P * atom["value"],
                "refined-grid delayed formula disagrees with inherited route")
        controls += 1
        if controls == 8:
            break
    require(controls == 8, "not enough positive independent route controls")

    def classify_triple(clock_triple):
        """Exhaust all 12^2 source-label pairs for one reversed clock path."""
        ell0, ell1, ell2 = clock_triple
        require(ell0 != ell1 and ell1 != ell2,
                "a tested clock edge repeats a label")
        histogram = Counter()
        zero_sources = []
        example = None
        for s0 in range(1, P):
            for s1 in range(1, P):
                current_bank = positive_atoms(s0, ell1, ell0)
                following_bank = positive_atoms(s1, ell2, ell1)
                positives = 0
                for current in current_bank:
                    for following in following_bank:
                        # This equality is forced by j(Dx)=h(x).  Imposing it
                        # before the exact interval sweep only removes pairs
                        # lying in disjoint digit cells of y={13^6 x}.
                        if current["h"] != following["j"]:
                            continue
                        value = fibre_product_value(current, following)
                        if value:
                            positives += 1
                            if example is None:
                                example = (
                                    s0, s1, value,
                                    tuple(current[name] for name in (
                                        "rail_index", "j", "h",
                                        "epsilon", "kappa",
                                    )),
                                    tuple(following[name] for name in (
                                        "rail_index", "j", "h",
                                        "epsilon", "kappa",
                                    )),
                                )
                histogram[positives] += 1
                if positives == 0:
                    zero_sources.append((s0, s1))
        return histogram, tuple(zero_sources), example

    first_hist, first_zeros, first_example = classify_triple((3, 1, 0))
    require(first_hist == Counter({20: 132, 0: 12}),
            "(3,1,0) fibre-product histogram changed")
    require(first_zeros == tuple((s, 6) for s in range(1, P)),
            "(3,1,0) inherited following-source exception changed")
    require(first_example == (
        1, 1, 126816337986097204341478787325120,
        (2, 5, 2, 0, 0),
        (0, 2, 6, 1, 1),
    ), "(3,1,0) positive numerator certificate changed")
    require(
        tuple(s for s in range(1, P) if not positive_atoms(s, 0, 1)) == (6,),
        "(3,1,0) inherited following-edge zero source changed",
    )

    second_hist, second_zeros, second_example = classify_triple((4, 0, 1))
    require(second_hist == Counter({14: 132, 0: 12}),
            "(4,0,1) fibre-product histogram changed")
    require(second_zeros == tuple((s, 7) for s in range(1, P)),
            "(4,0,1) inherited following-source exception changed")
    require(second_example == (
        1, 1, 90407221954729505049108909024000,
        (1, 11, 3, 0, 0),
        (3, 3, 5, 1, 0),
    ), "(4,0,1) positive numerator certificate changed")
    require(
        tuple(s for s in range(1, P) if not positive_atoms(s, 1, 0)) == (7,),
        "(4,0,1) inherited following-edge zero source changed",
    )

    # Preserve the old zero control but type its mechanism honestly: the
    # following cell (owner=2, shallow=1) is already a universally empty
    # one-edge incidence cell.  It is not a hostile to D-gluing itself.
    inherited_zero_hist, inherited_zero_pairs, inherited_zero_example = (
        classify_triple((3, 1, 2))
    )
    require(
        inherited_zero_hist == Counter({0: 144})
        and len(inherited_zero_pairs) == 144
        and inherited_zero_example is None
        and all(not positive_atoms(s, 2, 1) for s in range(1, P)),
        "inherited universal-zero clock control changed",
    )

    # Sharp D-hostile: at source pair (1,1), both one-edge relations and
    # their formal Boolean label product are positive, yet the actual
    # D-pullback fibre product is empty.  The physical zero persists after
    # exhausting all 144 independently retained source pairs.
    sharp_hostile_hist, sharp_hostile_pairs, sharp_hostile_example = (
        classify_triple((0, 1, 0))
    )
    sharp_current_relation = atom_relation(1, 1, 0)
    sharp_following_relation = atom_relation(1, 0, 1)
    sharp_formal_support = relation_support(boolean_product(
        sharp_current_relation, sharp_following_relation
    ))
    require(
        sharp_hostile_hist == Counter({0: 144})
        and len(sharp_hostile_pairs) == 144
        and sharp_hostile_example is None
        and relation_support(sharp_current_relation) > 0
        and relation_support(sharp_following_relation) > 0
        and sharp_formal_support > 0,
        "sharp formal-positive/D-empty hostile changed",
    )

    # Cheap broader hostile/control: at the one generic source pair (1,1),
    # exhaust all 7*6*6 reversed two-edge clock triples.  This is not asserted
    # to be the union over source drift; the full 12^2 exhaustion above is the
    # scoped result for the four displayed triples.
    generic_clock_support = []
    generic_edge_positive_triples = 0
    generic_formal_positive_triples = 0
    for ell0 in range(Q7):
        for ell1 in range(Q7):
            if ell0 == ell1:
                continue
            for ell2 in range(Q7):
                if ell1 == ell2:
                    continue
                current_relation = atom_relation(1, ell1, ell0)
                following_relation = atom_relation(1, ell2, ell1)
                if (relation_support(current_relation) > 0
                        and relation_support(following_relation) > 0):
                    generic_edge_positive_triples += 1
                if relation_support(boolean_product(
                        current_relation, following_relation)) > 0:
                    generic_formal_positive_triples += 1
                positives = 0
                for current in positive_atoms(1, ell1, ell0):
                    for following in positive_atoms(1, ell2, ell1):
                        if (current["h"] == following["j"]
                                and fibre_product_value(current, following)):
                            positives += 1
                if positives:
                    generic_clock_support.append(
                        (ell0, ell1, ell2, positives)
                    )
    require(tuple(generic_clock_support) == EXPECTED_GENERIC_CLOCK_SUPPORT,
            "generic reversed-clock support atlas changed")
    require(
        generic_edge_positive_triples == 146
        and generic_formal_positive_triples == 146
        and len(generic_clock_support) == 10,
        "formal-versus-physical two-edge clock census changed",
    )

    print("LRC14 exact dilation-handoff reversed-clock fibre-product probe")
    print(
        "scope=THM-2616 raw nonnegative guard-safe carrier; sharp graph "
        "and predecessor/carry labels"
    )
    print("handoff=D(x)={13x}; owner_current=shallow_following and h_current=j_following")
    print(
        "clock_orientation=D-compatible candidate chronology traverses printed "
        "(owner,shallow) edges backwards"
    )
    print(
        f"direct_owner_disjoint_checks={direct_owner_disjoint_checks} "
        f"dilation_clock_covariance_checks={dilation_clock_checks}"
    )
    print(f"refined_grid_formula_controls={controls}")
    print(
        "triple_(shallow,current_owner,next_owner)=(3,1,0) "
        "source_pair_positive_rail_labelled_atom_pair_count_hist="
        f"{tuple(sorted(first_hist.items()))} "
        "zero_pairs=((s,6):s=1..12)"
    )
    print(f"triple_310_first_positive_numerator_certificate={first_example}")
    print(
        "triple_(shallow,current_owner,next_owner)=(4,0,1) "
        "source_pair_positive_rail_labelled_atom_pair_count_hist="
        f"{tuple(sorted(second_hist.items()))} "
        "zero_pairs=((s,7):s=1..12)"
    )
    print(f"triple_401_first_positive_numerator_certificate={second_example}")
    print(
        "inherited_zero_edge_control_"
        "(shallow,current_owner,next_owner)=(3,1,2) "
        "source_pair_positive_rail_labelled_atom_pair_count_hist=((0,144),)"
    )
    print(
        "sharp_formal_positive_D_empty_hostile_"
        "(shallow,current_owner,next_owner)=(0,1,0) "
        f"formal_support_at_source_pair_11={sharp_formal_support} "
        "physical_source_pair_hist=((0,144),)"
    )
    print(
        "generic_source_pair_11_two_edge_clock_census="
        f"(edge_positive={generic_edge_positive_triples},"
        f"formal_positive={generic_formal_positive_triples},"
        f"physical_positive={len(generic_clock_support)},total=252)"
    )
    print(
        "generic_source_pair_11_physical_clock_support="
        f"{tuple(generic_clock_support)}"
    )
    print(
        "verdict=PASS: a correctly typed two-edge fibre product is nonempty, "
        "but sparse and source-sensitive"
    )
    print(
        "semantics=the former fixed-offset Boolean zero is not a physical "
        "chronology no-go; the next gate is iterated D-coherence with source drift"
    )


if __name__ == "__main__":
    main()
