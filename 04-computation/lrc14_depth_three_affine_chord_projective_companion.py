#!/usr/bin/env python3
"""Exact affine-chord projective companion to THM-2250.

The scalar valuation profile is (3,4,5).  Put

    U_k = union_j D_(13^(lambda_j-k) u_j),  0 <= k <= 3.

Then U_1 is the second forward image of U_3, so

    Y = (169 1_(U_1) - 1_(U_3)) / 168

has zero signed response against the residual R.  It is one on the negative
carrier U_1 intersect U_3.  Thus the residual mass is at most

    E[1_(C_H intersect U_0 intersect U_2) max(Y,0)].

The usual arbitrary-coupling Bellman forgets the affine root chord.  After
normalizing the thirteen roots by H, a danger process with unit core u has
the fixed nonzero chord direction

    a = -H u^(-1) mod 13.

If its future bit is zero, the two danger roots are
``{beta,beta+a}``; if its future bit is one, the unique danger root is one
of those endpoints.  The affine anchor ``beta`` depends on the parent and
is therefore allowed to range independently over all thirteen roots at
every Bellman history.  The guard roots are {2,...,10}, with optional root
1 or 11 when the future guard bit is zero.

This companion fixes (a_1,a_2,a_3) throughout the Bellman recursion, allows
every affine anchor and phase orientation independently at every history (a
safe enlargement), and enumerates all 12^3 oriented direction triples.  Since
an affine two-point chord does not distinguish ``a`` from ``-a``, the actual
state space is ``F_13^*/{+-1}``.  Every transition is a literal thirteen-root
count, so the dynamic program uses integers.  Only the terminal arbitrary
coupling requires an LP; all exact vertices are enumerated and every distinct
terminal objective receives an exact primal/dual check.

An independent direct Fraction root enumeration checks the advertised guard
and danger root sets for every label and phase class.  Ordinary and ``-O``
executions are intended to be byte-identical.

Canonical status: this is a finite-exact explanatory companion, not a
current proof-graph theorem.  THM-2250 subsequently proved the stronger
literal equality of all three normalized cores, and THM-2257/THM-2258
independently closed the whole profile.
"""

from argparse import ArgumentParser
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
from multiprocessing import get_context

import numpy as np


P = 13
PROFILE = (3, 4, 5)
TARGET = Fraction(961, 6930)
SCORE_DENOMINATOR = 168
TERMINAL_VERTEX_DENOMINATOR = 42

BIT_VECTORS = tuple(product((0, 1), repeat=4))
BIT_INDEX = {bits: index for index, bits in enumerate(BIT_VECTORS)}
STATES = tuple(range(len(BIT_VECTORS)))
COUPLING_COLUMNS = tuple(
    (Fraction(1),) + tuple(Fraction(bit) for bit in bits)
    for bits in BIT_VECTORS
)
RANK = 5
TERMINAL_RHS = (
    Fraction(1),
    Fraction(5, 7),
    Fraction(1, 7),
    Fraction(1, 7),
    Fraction(1, 7),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def solve_square(matrix, rhs):
    """Solve a square rational system, or return None when singular."""
    size = len(rhs)
    augmented = [list(matrix[row]) + [rhs[row]] for row in range(size)]
    for column in range(size):
        pivot = next(
            (
                row
                for row in range(column, size)
                if augmented[row][column]
            ),
            None,
        )
        if pivot is None:
            return None
        augmented[column], augmented[pivot] = (
            augmented[pivot],
            augmented[column],
        )
        scale = augmented[column][column]
        augmented[column] = [
            entry / scale for entry in augmented[column]
        ]
        for row in range(size):
            if row == column:
                continue
            scale = augmented[row][column]
            if scale:
                augmented[row] = [
                    augmented[row][entry]
                    - scale * augmented[column][entry]
                    for entry in range(size + 1)
                ]
    return tuple(augmented[row][-1] for row in range(size))


def basis_matrix(basis):
    return tuple(
        tuple(
            COUPLING_COLUMNS[basis[column]][row]
            for column in range(RANK)
        )
        for row in range(RANK)
    )


INVERTIBLE_BASES = tuple(
    basis
    for basis in combinations(STATES, RANK)
    if solve_square(
        basis_matrix(basis),
        (Fraction(1),) + (Fraction(0),) * 4,
    )
    is not None
)


def terminal_vertices():
    """Return exact terminal vertices together with all witnessing bases."""
    by_distribution = {}
    for basis in INVERTIBLE_BASES:
        weights = solve_square(basis_matrix(basis), TERMINAL_RHS)
        require(weights is not None, "terminal basis became singular")
        if any(weight < 0 for weight in weights):
            continue
        distribution = [Fraction() for _ in STATES]
        for state, weight in zip(basis, weights):
            distribution[state] = weight
        distribution = tuple(distribution)
        by_distribution.setdefault(distribution, []).append(basis)
    require(by_distribution, "terminal coupling polytope has no vertex")
    return tuple(
        (distribution, tuple(bases))
        for distribution, bases in sorted(by_distribution.items())
    )


TERMINAL_VERTICES = terminal_vertices()
TERMINAL_INTEGER_WEIGHTS = tuple(
    tuple(
        int(weight * TERMINAL_VERTEX_DENOMINATOR)
        for weight in distribution
    )
    for distribution, _ in TERMINAL_VERTICES
)
require(
    all(
        Fraction(integer, TERMINAL_VERTEX_DENOMINATOR) == weight
        for (distribution, _), integers in zip(
            TERMINAL_VERTICES,
            TERMINAL_INTEGER_WEIGHTS,
        )
        for weight, integer in zip(distribution, integers)
    ),
    "terminal common denominator drift",
)


TERMINAL_LP_CHECKS = 0


@lru_cache(maxsize=None)
def maximize_terminal(values):
    """Exact terminal primal maximum plus a matching exact basic dual."""
    global TERMINAL_LP_CHECKS
    TERMINAL_LP_CHECKS += 1
    scored = tuple(
        (
            sum(
                weight * value
                for weight, value in zip(distribution, values)
            ),
            distribution,
            bases,
        )
        for distribution, bases in TERMINAL_VERTICES
    )
    primal = max(row[0] for row in scored)
    candidate_bases = tuple(
        basis
        for value, _, bases in scored
        if value == primal
        for basis in bases
    )
    for basis in candidate_bases:
        dual = solve_square(
            tuple(COUPLING_COLUMNS[state] for state in basis),
            tuple(values[state] for state in basis),
        )
        require(dual is not None, "terminal dual basis became singular")
        if any(
            sum(
                feature * multiplier
                for feature, multiplier in zip(column, dual)
            )
            < value
            for column, value in zip(COUPLING_COLUMNS, values)
        ):
            continue
        dual_value = sum(
            rhs * multiplier
            for rhs, multiplier in zip(TERMINAL_RHS, dual)
        )
        require(dual_value == primal, "terminal primal/dual value drift")
        return primal
    raise RuntimeError("terminal optimum has no exact basic dual")


def fractional_part(value):
    return value - value.numerator // value.denominator


def circle_distance(value):
    value = fractional_part(value)
    return min(value, 1 - value)


def inverse_mod(value):
    return pow(value, -1, P)


def danger_direct_roots(label, anchor, phase):
    """Directly enumerate danger roots in guard-normalized root labels."""
    unit = (-inverse_mod(label)) % P
    require(unit, "danger unit residue vanished")
    integer_part = (anchor * inverse_mod(label)) % P
    roots = set()
    for root in range(P):
        current_phase = (
            Fraction(integer_part + unit * root, P) + phase / P
        )
        if circle_distance(current_phase) < Fraction(1, 14):
            roots.add(root)
    return frozenset(roots)


def guard_direct_roots(phase):
    """Directly enumerate guard roots after setting H=1 mod 13."""
    roots = set()
    for root in range(P):
        current_phase = Fraction(root, P) + phase / P
        if circle_distance(current_phase) > Fraction(1, 7):
            roots.add(root)
    return frozenset(roots)


DANGER_PHASES = {
    "future_zero": Fraction(1, 2),
    "future_one_anchor_root": Fraction(1, 28),
    "future_one_shifted_root": Fraction(27, 28),
}
GUARD_PHASES = {
    "future_one": Fraction(1, 2),
    "future_zero_root_11": Fraction(1, 28),
    "future_zero_root_1": Fraction(27, 28),
}


def direct_root_audit():
    """Independent exact check of every advertised one-coordinate root set."""
    checks = 0
    for label in range(1, P):
        for anchor in range(P):
            expected = {
                "future_zero": frozenset(
                    (anchor, (anchor + label) % P)
                ),
                "future_one_anchor_root": frozenset((anchor,)),
                "future_one_shifted_root": frozenset(
                    ((anchor + label) % P,)
                ),
            }
            for name, phase in DANGER_PHASES.items():
                require(
                    danger_direct_roots(label, anchor, phase)
                    == expected[name],
                    (
                        "direct danger root formula failed at "
                        f"label={label} anchor={anchor} {name}"
                    ),
                )
                checks += 1
    guard_expected = {
        "future_one": frozenset(range(2, 11)),
        "future_zero_root_11": frozenset(range(2, 12)),
        "future_zero_root_1": frozenset(range(1, 11)),
    }
    for name, phase in GUARD_PHASES.items():
        require(
            guard_direct_roots(phase) == guard_expected[name],
            f"direct guard root formula failed at {name}",
        )
        checks += 1
    return checks


DIRECT_ROOT_CHECKS = direct_root_audit()


def guard_root_options(future_bit):
    if future_bit:
        return (frozenset(range(2, 11)),)
    return (
        frozenset(range(2, 12)),
        frozenset(range(1, 11)),
    )


def danger_root_options(label, future_bit):
    if future_bit:
        return tuple(frozenset((root,)) for root in range(P))
    return tuple(
        frozenset((anchor, (anchor + label) % P))
        for anchor in range(P)
    )


def sign_quotient_audit():
    """Check directly that every transition option forgets chord orientation."""
    checks = 0
    for label in range(1, (P + 1) // 2):
        for future_bit in (0, 1):
            positive = frozenset(danger_root_options(label, future_bit))
            negative = frozenset(
                danger_root_options(P - label, future_bit)
            )
            require(
                positive == negative,
                (
                    "orientation quotient failed at "
                    f"label={label} future_bit={future_bit}"
                ),
            )
            checks += 1
    return checks


SIGN_QUOTIENT_CHECKS = sign_quotient_audit()


def root_policies(labels, future_state):
    """All safe enlarged thirteen-root incidence policies for one state."""
    future_bits = BIT_VECTORS[future_state]
    options = (
        guard_root_options(future_bits[0]),
        *(
            danger_root_options(label, future_bits[index + 1])
            for index, label in enumerate(labels)
        ),
    )
    rows = []
    raw_count = 1
    for choices in options:
        raw_count *= len(choices)
    for guard_roots, danger_1, danger_2, danger_3 in product(*options):
        counts = np.zeros(len(STATES), dtype=np.int64)
        for root in range(P):
            bits = (
                int(root in guard_roots),
                int(root in danger_1),
                int(root in danger_2),
                int(root in danger_3),
            )
            counts[BIT_INDEX[bits]] += 1
        require(sum(counts) == P, "root policy lost a root")
        require(
            sum(
                counts[state] * BIT_VECTORS[state][0]
                for state in STATES
            )
            == 10 - future_bits[0],
            "guard root marginal drift",
        )
        for core in range(3):
            require(
                sum(
                    counts[state] * BIT_VECTORS[state][core + 1]
                    for state in STATES
                )
                == 2 - future_bits[core + 1],
                "danger root marginal drift",
            )
        rows.append(counts)
    unique = np.unique(np.asarray(rows, dtype=np.int64), axis=0)
    return unique, raw_count


HIT_MASKS = []
for time in range(PROFILE[-1] + 1):
    row = []
    for bits in BIT_VECTORS:
        hit = int(time == 0 and bits[0])
        for checkpoint in range(PROFILE[0] + 1):
            if any(
                bits[core + 1]
                and time == PROFILE[core] - checkpoint
                for core in range(3)
            ):
                hit |= 1 << (checkpoint + 1)
        row.append(hit)
    HIT_MASKS.append(tuple(row))
HIT_MASKS = tuple(HIT_MASKS)


FULL_DENOMINATOR = (
    SCORE_DENOMINATOR
    * (P ** (PROFILE[-1] + 1))
    * TERMINAL_VERTEX_DENOMINATOR
)


def profile_row(labels):
    """Return one exact integer Bellman row and its terminal objective."""
    policy_data = tuple(
        root_policies(labels, future_state)
        for future_state in STATES
    )
    policies = tuple(row[0] for row in policy_data)
    raw_policy_count = sum(row[1] for row in policy_data)
    unique_policy_count = sum(len(row) for row in policies)

    # Summary bits are guard,U_0,U_1,U_2,U_3.
    value = np.zeros((32, len(STATES)), dtype=np.int64)
    for summary in range(32):
        guard = summary & 1
        u_0 = (summary >> 1) & 1
        u_1 = (summary >> 2) & 1
        u_2 = (summary >> 3) & 1
        u_3 = (summary >> 4) & 1
        payoff_numerator = 0
        if guard and u_0 and u_2:
            payoff_numerator = max(169 * u_1 - u_3, 0)
        value[summary, :] = payoff_numerator

    for time, hit_masks in enumerate(HIT_MASKS):
        hit_masks = np.asarray(hit_masks, dtype=np.int64)
        updated = np.zeros_like(value)
        for summary in range(32):
            costs = value[
                np.bitwise_or(summary, hit_masks),
                np.arange(len(STATES)),
            ]
            for future_state in STATES:
                updated[summary, future_state] = int(
                    np.max(policies[future_state] @ costs)
                )
        value = updated

    terminal_values = tuple(int(value) for value in value[0])
    integer_terminal = max(
        sum(weight * cost for weight, cost in zip(vertex, terminal_values))
        for vertex in TERMINAL_INTEGER_WEIGHTS
    )
    return (
        tuple(labels),
        integer_terminal,
        terminal_values,
        raw_policy_count,
        unique_policy_count,
    )


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "--workers",
        type=int,
        default=8,
        help="parallel exact label rows (default: 8; use 1 for serial)",
    )
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    require(P == 13 and PROFILE == (3, 4, 5), "frozen universe drift")
    require(len(INVERTIBLE_BASES) == 3008, "terminal basis count drift")
    require(len(TERMINAL_VERTICES) == 34, "terminal vertex count drift")
    require(DIRECT_ROOT_CHECKS == 471, "direct root audit count drift")
    require(SIGN_QUOTIENT_CHECKS == 12, "sign quotient audit drift")

    # Each affine chord family is unchanged by reversing its direction.  Run
    # only the 6^3 projective triples, then expand them back to the complete
    # 12^3 oriented census so the scope count remains explicit.
    label_universe = tuple(
        product(range(1, (P + 1) // 2), repeat=3)
    )
    if args.workers == 1:
        projective_rows = tuple(
            profile_row(labels) for labels in label_universe
        )
    else:
        with get_context("fork").Pool(processes=args.workers) as pool:
            projective_rows = tuple(
                pool.imap(profile_row, label_universe, chunksize=2)
            )
    projective_rows = tuple(sorted(projective_rows))
    require(len(projective_rows) == 6**3, "projective census drift")
    projective_lookup = {row[0]: row for row in projective_rows}
    detailed_rows = tuple(
        (
            labels,
            *projective_lookup[
                tuple(min(label, P - label) for label in labels)
            ][1:],
        )
        for labels in product(range(1, P), repeat=3)
    )
    detailed_rows = tuple(sorted(detailed_rows))
    rows = tuple((row[0], row[1]) for row in detailed_rows)
    raw_policy_count = sum(row[3] for row in detailed_rows)
    unique_policy_count = sum(row[4] for row in detailed_rows)
    projective_raw_policy_count = sum(row[3] for row in projective_rows)
    projective_unique_policy_count = sum(
        row[4] for row in projective_rows
    )

    terminal_objectives = {
        row[2]: row[1] for row in detailed_rows
    }
    for values, integer_terminal in terminal_objectives.items():
        exact_terminal = maximize_terminal(
            tuple(Fraction(value) for value in values)
        )
        require(
            exact_terminal
            == Fraction(
                integer_terminal,
                TERMINAL_VERTEX_DENOMINATOR,
            ),
            "integer/exact terminal maximization drift",
        )

    require(len(rows) == 12**3, "projective label census drift")
    same_direction_rows = tuple(
        row
        for row in rows
        if len({min(label, P - label) for label in row[0]}) == 1
    )
    mixed_direction_rows = tuple(
        row
        for row in rows
        if len({min(label, P - label) for label in row[0]}) > 1
    )
    require(len(same_direction_rows) == 48, "same-direction census drift")
    require(len(mixed_direction_rows) == 1680, "mixed-direction census drift")

    def exact_value(row):
        return Fraction(row[1], FULL_DENOMINATOR)

    worst_all = max(rows, key=lambda row: row[1])
    worst_same_direction = max(
        same_direction_rows, key=lambda row: row[1]
    )
    best_same_direction = min(
        same_direction_rows, key=lambda row: row[1]
    )
    worst_mixed_direction = max(
        mixed_direction_rows, key=lambda row: row[1]
    )
    all_maximizers = tuple(
        labels
        for labels, numerator in rows
        if numerator == worst_all[1]
    )
    mixed_direction_maximizers = tuple(
        labels
        for labels, numerator in mixed_direction_rows
        if numerator == worst_mixed_direction[1]
    )
    same_direction_minimizers = tuple(
        labels
        for labels, numerator in same_direction_rows
        if numerator == best_same_direction[1]
    )
    passing_mixed_direction = tuple(
        row for row in mixed_direction_rows if exact_value(row) < TARGET
    )
    failing_same_direction = tuple(
        row for row in same_direction_rows if exact_value(row) >= TARGET
    )

    require(
        len(passing_mixed_direction) == 1680,
        "mixed-direction closure drift",
    )
    require(
        len(failing_same_direction) == 48,
        "same-direction failure count drift",
    )
    require(
        exact_value(worst_mixed_direction)
        == Fraction(567556891, 5676327384),
        "mixed-direction worst bound drift",
    )
    require(
        TARGET - exact_value(worst_mixed_direction)
        == Fraction(36232889557, 936594018360),
        "mixed-direction target margin drift",
    )
    require(
        exact_value(worst_same_direction)
        == Fraction(8907541, 62377224),
        "same-direction hostile maximum drift",
    )
    require(
        exact_value(worst_same_direction) - TARGET
        == Fraction(42493973, 10292241960),
        "same-direction hostile excess drift",
    )
    require(
        exact_value(best_same_direction)
        == Fraction(5168763, 36386714),
        "same-direction closest hostile bound drift",
    )
    require(
        exact_value(best_same_direction) - TARGET
        == Fraction(30424837, 9005711715),
        "same-direction closest hostile excess drift",
    )
    require(len(mixed_direction_maximizers) == 80, "mixed maximizer drift")
    require(len(all_maximizers) == 16, "all maximizer drift")
    require(
        len(same_direction_minimizers) == 24,
        "same-direction minimizer drift",
    )
    require(
        raw_policy_count == 91113984
        and unique_policy_count == 1338240,
        "oriented transition policy census drift",
    )
    require(
        projective_raw_policy_count == 11389248
        and projective_unique_policy_count == 167280,
        "projective transition policy census drift",
    )

    transcript = "\n".join(
        f"{labels[0]},{labels[1]},{labels[2]}|{numerator}"
        for labels, numerator in rows
    )
    row_digest = sha256(transcript.encode("ascii")).hexdigest()
    require(
        row_digest
        == "39d7332e4d8bbeb565dad8d9ebcf2429dd8da8fb44662b4c97f8fd7c38121960",
        "label-row digest drift",
    )

    print("DEPTH-THREE AFFINE-CHORD PROJECTIVE COMPANION")
    print(f"profile={PROFILE} target={TARGET}")
    print("score=(169*1_U1-1_U3)/168")
    print("positive_carrier=C_H intersect U0 intersect U2")
    print("negative_carrier=U1 intersect U3 score_on_negative=1")
    print("chord_direction=a_j=-H*u_j^(-1) mod 13, modulo sign")
    print("affine_anchor=beta_j arbitrary per parent/history")
    print(f"direct_fraction_root_checks={DIRECT_ROOT_CHECKS}")
    print(f"direct_sign_quotient_checks={SIGN_QUOTIENT_CHECKS}")
    print(
        f"terminal_invertible_bases={len(INVERTIBLE_BASES)} "
        f"terminal_vertices={len(TERMINAL_VERTICES)}"
    )
    print(
        f"raw_transition_policies={raw_policy_count} "
        f"unique_transition_policies={unique_policy_count}"
    )
    print(
        f"projective_raw_transition_policies={projective_raw_policy_count} "
        "projective_unique_transition_policies="
        f"{projective_unique_policy_count}"
    )
    print("label_universe=12^3=1728")
    print("mixed_direction_labels=1680 same_direction_labels=48")
    print(
        "mixed_direction_closed=1680/1680 "
        f"worst_labels={worst_mixed_direction[0]} "
        f"worst_bound={exact_value(worst_mixed_direction)} "
        f"target_margin={TARGET-exact_value(worst_mixed_direction)}"
    )
    print(
        f"mixed_direction_worst_maximizers={mixed_direction_maximizers}"
    )
    print(
        "same_direction_failing=48/48 "
        f"worst_labels={worst_same_direction[0]} "
        f"worst_bound={exact_value(worst_same_direction)} "
        f"bound_minus_target={exact_value(worst_same_direction)-TARGET}"
    )
    print(
        "same_direction_closest_to_target "
        f"labels={best_same_direction[0]} "
        f"bound={exact_value(best_same_direction)} "
        f"bound_minus_target={exact_value(best_same_direction)-TARGET}"
    )
    print(f"same_direction_minimizers={same_direction_minimizers}")
    print(f"all_label_maximizers={all_maximizers}")
    print(f"row_digest={row_digest}")
    require(TERMINAL_LP_CHECKS == 23, "terminal LP check census drift")
    print(f"distinct_terminal_primal_dual_checks={TERMINAL_LP_CHECKS}")
    print("integer_terminal_vertex_parity=PASS")
    print("independent_direct_root_incidence=PASS")
    print(
        "finite_exact_companion_consequence=every_hypothetical_"
        "(3,4,5)_survivor_has_u1=+-u2=+-u3_mod_13"
    )
    print(
        "canonical_status=SUPERSEDED_IN_CONSEQUENCE_BY_"
        "THM-2250_exact_core_equality"
    )
    print(
        "current_scope=profile_EMPTY_independently_by_"
        "THM-2257_and_THM-2258_LRC14_OPEN"
    )


if __name__ == "__main__":
    main()
