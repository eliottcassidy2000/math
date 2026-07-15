#!/usr/bin/env python3
"""Exact verifier for the shell-five universal endpoint-grid obstruction.

This continuation concerns only a *single numerator chosen before the seven
free signed lifts are known*.  Put

    P_d = {u <= B : u mod 13 in {1,2,4,5,8,11,12}}
          union {B-3,B-2,B},       B=(13d+5)/2.

A point which is deep for every shell-five lift must be deep for every member
of P_d.  For d mod 52 in {15,37,41}, the proof in THM-836 localizes that
universal deep set to four explicit intervals, and none contains a unit point
on either endpoint grid q=5d or q=13d.

This does NOT say that a numerator selected after a particular ten-speed set
U is known cannot work.  It is an obstruction to the uniform one-column
template that succeeded for d=11 mod 52.

Every verdict below uses integers or fractions.Fraction.  The finite audits
are replays of all-size inequalities stated in THM-836, not extrapolations.
"""

from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
from math import gcd


FREE_RESIDUES = (1, 2, 4, 5, 8, 11, 12)
OPEN_CLASSES_MOD_52 = (15, 37, 41)
DEEP = F(1, 11)


def balanced(value: int, modulus: int) -> int:
    residue = value % modulus
    if 2 * residue > modulus:
        residue -= modulus
    return residue


def norm(value: F) -> F:
    residue = value % 1
    return min(residue, 1 - residue)


def shell_B(d: int) -> int:
    assert d > 0 and d % 2 == 1
    return (13 * d + 5) // 2


def candidate_pool(d: int) -> tuple[int, ...]:
    B = shell_B(d)
    free = [u for u in range(1, B + 1) if u % 13 in FREE_RESIDUES]
    return tuple(sorted(set(free + [B - 3, B - 2, B])))


def interval_intersection(
    left: tuple[tuple[F, F], ...], right: tuple[tuple[F, F], ...]
) -> tuple[tuple[F, F], ...]:
    answer: list[tuple[F, F]] = []
    i = j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo <= hi:
            answer.append((lo, hi))
        if left[i][1] < right[j][1]:
            i += 1
        else:
            j += 1
    return tuple(answer)


def deep_intervals(speeds: tuple[int, ...]) -> tuple[tuple[F, F], ...]:
    """Closed components of {t in [0,1]: ||ut||>=1/11 for every u}."""

    answer = ((F(0), F(1)),)
    for speed in speeds:
        safe = tuple(
            (F(11 * k + 1, 11 * speed), F(11 * k + 10, 11 * speed))
            for k in range(speed)
        )
        answer = interval_intersection(answer, safe)
    return answer


SKELETON = tuple(u for u in range(1, 65) if u % 13 in FREE_RESIDUES)
C2 = (F(89, 583), F(109, 704))
C4 = (F(47, 154), F(219, 704))
C9 = (F(243, 704), F(164, 473))


def reflect(interval: tuple[F, F]) -> tuple[F, F]:
    return (1 - interval[1], 1 - interval[0])


SKELETON_COMPONENTS = (C2, C4, C9, reflect(C9), reflect(C4), reflect(C2))


def skeleton_audit() -> tuple[object, ...]:
    components = deep_intervals(SKELETON)
    assert components == SKELETON_COMPONENTS

    # Within each skeleton arc, successive phases in the progressions used
    # by the no-jump lemma move by strictly less than the forbidden width
    # 2/11.  The initial phase remains on the indicated signed lift.
    centers = (F(2, 13), F(4, 13), F(9, 26))
    steps = (13, 13, 26)
    for center, (lo, hi), step in zip(centers, (C2, C4, C9), steps):
        assert step * max(center - lo, hi - center) < F(2, 11)

    # Positive carriers for left motion and negative carriers for right
    # motion: raw 1/12 at 2/13, raw 4/12 at 4/13, and mod-26 states 12/17
    # at 9/26.
    signed_initial_phases = (
        (F(2, 13) + C2[0], -F(2, 13) + C2[1]),
        (F(3, 13) + 4 * (C4[0] - F(4, 13)), -F(4, 13) + 12 * (C4[1] - F(4, 13))),
        (F(2, 13) + 12 * (C9[0] - F(9, 26)), -F(3, 26) + 17 * (C9[1] - F(9, 26))),
    )
    # The first entry in the first row is only a coarse positivity check;
    # the exact progression starts at phase 2/13+x for speed 1.
    assert F(2, 13) + (C2[0] - F(2, 13)) > 0
    assert -F(2, 13) + 12 * (C2[1] - F(2, 13)) < 0
    for positive, negative in signed_initial_phases[1:]:
        assert positive > 0
        assert negative < 0

    return (
        len(SKELETON),
        SKELETON,
        components,
        tuple(
            (center, step * max(center - lo, hi - center))
            for center, (lo, hi), step in zip(centers, (C2, C4, C9), steps)
        ),
    )


def residue_parity_center_audit() -> tuple[object, ...]:
    states = tuple(
        next(u for u in range(1, 27) if u % 13 == residue and u % 2 == parity)
        for residue in FREE_RESIDUES
        for parity in (0, 1)
    )
    free_safe = tuple(
        a
        for a in range(1, 26)
        if all(abs(balanced(a * u, 26)) >= 3 for u in states)
    )
    assert free_safe == (4, 8, 9, 17, 18, 22)

    rows = []
    for residue_class in (11, 15, 37, 41):
        B = shell_B(residue_class)
        forced_states = (B - 3, B - 2, B)
        safe = tuple(
            a
            for a in free_safe
            if all(abs(balanced(a * u, 26)) >= 3 for u in forced_states)
        )
        assert safe == (8, 18)
        rows.append(
            (
                residue_class,
                B % 26,
                tuple(u % 26 for u in forced_states),
                safe,
            )
        )
    return (states, free_safe, tuple(rows))


def localization_intervals(d: int) -> tuple[tuple[F, F], tuple[F, F]]:
    """The two positive-half arcs allowed by the all-size localization proof."""

    B = shell_B(d)
    around_4_over_13 = (
        F(4, 13) - F(9, 143 * (B - 2)),
        F(4, 13) + F(9, 143 * (B - 3)),
    )
    if d % 52 == 15:
        around_9_over_26 = (
            F(9, 26) - F(9, 143 * (B - 10)),
            F(9, 26) - F(2, 143 * (B - 2)),
        )
    else:
        assert d % 52 in (37, 41)
        around_9_over_26 = (
            F(9, 26) + F(2, 143 * (B - 3)),
            F(9, 26) + F(7, 286 * (B - 18)),
        )
    return around_4_over_13, around_9_over_26


def all_size_inequality_audit() -> tuple[object, ...]:
    """Integer reductions used to keep forced phases on one signed lift."""

    rows = []
    for residue_class in OPEN_CLASSES_MOD_52:
        B = shell_B(residue_class)
        assert B >= 100
        assert B % 13 == 9

        # 2/13 branch: free raw 1 and 12 give
        # -9/[143(B-8)] <= x <= 9/[143(B-10)].  The inequalities below
        # exclude the far-side safe arcs of forced B-3 and B-2.
        assert 9 * (B - 3) < 24 * (B - 10)
        assert 9 * (B - 2) < 24 * (B - 8)

        # 4/13 branch: raw 4 and 12 first localize x; the forced pair then
        # gives the displayed I_4 interval.
        assert 20 * (B - 2) < 35 * (B - 5)
        assert 31 * (B - 3) < 35 * (B - 10)
        assert F(2, 13) + F(31 * (B - 2), 143 * (B - 10)) < F(1, 2)
        assert F(2, 13) + F(20 * (B - 3), 143 * (B - 5)) < F(1, 2)

        if residue_class == 15:
            assert B % 26 == 22
            # At 9/26, the top +4/26 free state is B-10 and forced B-2
            # has phase -2/26.  The free upper bound also prevents a jump
            # to the positive far-side safe arc.
            assert 7 * (B - 2) < 48 * (B - 5)
        else:
            assert B % 26 == 9
            # The top -3/26 free state is B-18 and forced B-3 has +2/26.
            assert 9 * (B - 3) < 24 * (B - 23)

        intervals = localization_intervals(residue_class)
        rows.append((residue_class, B, intervals))

    # Every displayed comparison becomes only stronger as B increases by
    # 338 when d increases by 52.  Record their positive linear margins at
    # the least B in each class.
    margins = {
        "c2_left": 24 * (100 - 10) - 9 * (100 - 3),
        "c2_right": 24 * (100 - 8) - 9 * (100 - 2),
        "c4_left": 35 * (100 - 5) - 20 * (100 - 2),
        "c4_right": 35 * (100 - 10) - 31 * (100 - 3),
        "c9_left_class15": 48 * (100 - 5) - 7 * (100 - 2),
        "c9_right_classes37_41": 24 * (243 - 23) - 9 * (243 - 3),
    }
    assert all(value > 0 for value in margins.values())
    return (tuple(rows), tuple(sorted(margins.items())))


def interval_contains(interval: tuple[F, F], value: F) -> bool:
    return interval[0] <= value <= interval[1]


def unit_grid_points_in_intervals(
    d: int, denominator: int, intervals: tuple[tuple[F, F], ...]
) -> tuple[int, ...]:
    return tuple(
        p
        for p in range(1, denominator)
        if gcd(p, denominator) == 1
        and any(interval_contains(interval, F(p, denominator)) for interval in intervals)
    )


def q13_direct_inverse_audit(d: int, p: int) -> tuple[int, int]:
    """Return the raw +/-1 danger speed from the direct q=13d proof."""

    q = 13 * d
    assert gcd(p, q) == 1
    j = balanced(p, 13)
    assert 1 <= abs(j) <= 6
    inverse = pow(p, -1, q)
    v = balanced(inverse * j, q)
    assert v != 0
    u = abs(v)
    assert u <= q // 2 < shell_B(d)
    assert u % 13 in (1, 12)
    assert norm(F(p * u, q)) == F(abs(j), q) < DEEP
    return u, j


def endpoint_grid_audit(limit_k: int = 24) -> tuple[object, ...]:
    rows = []
    inverse_hasher = sha256()
    for residue_class in OPEN_CLASSES_MOD_52:
        for k in range(limit_k + 1):
            d = residue_class + 52 * k
            B = shell_B(d)
            I4, I9 = localization_intervals(d)
            intervals = (I4, I9, reflect(I9), reflect(I4))

            q5 = 5 * d
            q13 = 13 * d
            assert not unit_grid_points_in_intervals(d, q5, intervals)
            assert not unit_grid_points_in_intervals(d, q13, intervals)

            # Exact gap reductions used in the all-size proof.
            assert F(45 * d, 143 * (B - 3)) < F(1, 13)
            assert F(117 * d, 143 * (B - 3)) < 1
            if residue_class == 15:
                assert (45 * d) % 26 == 25
                assert F(45 * d, 143 * (B - 10)) < F(25, 26)
                assert F(117 * d, 143 * (B - 10)) < F(1, 2)
            else:
                assert (45 * d) % 26 in (1, 25)
                assert F(35 * d, 286 * (B - 18)) < F(1, 26)
                assert F(91 * d, 286 * (B - 18)) < F(1, 2)

            # Directly replay every q=13d unit numerator at the first row of
            # each class; the displayed inverse argument is all-size.
            if k == 0:
                inverse_rows = tuple(
                    (p, *q13_direct_inverse_audit(d, p))
                    for p in range(1, q13)
                    if gcd(p, q13) == 1
                )
                inverse_hasher.update(repr(inverse_rows).encode())

            rows.append((d, B, I4, I9))

    return (
        limit_k,
        len(rows),
        rows[0],
        rows[-1],
        sha256(repr(rows).encode()).hexdigest(),
        inverse_hasher.hexdigest(),
    )


def finite_full_pool_replay() -> tuple[object, ...]:
    """Reconnaissance only: scan several representatives without extrapolation."""

    rows = []
    for residue_class in OPEN_CLASSES_MOD_52:
        for k in range(3):
            d = residue_class + 52 * k
            pool = candidate_pool(d)
            counts = []
            for denominator in (5 * d, 13 * d):
                universal = tuple(
                    p
                    for p in range(1, denominator)
                    if gcd(p, denominator) == 1
                    and all(norm(F(p * u, denominator)) >= DEEP for u in pool)
                )
                assert not universal
                counts.append((denominator, len(universal)))
            rows.append((d, shell_B(d), len(pool), tuple(counts)))
    return (tuple(rows), sha256(repr(rows).encode()).hexdigest())


def tournament_telemetry() -> tuple[object, ...]:
    """Telemetry on proof obligations, not a theorem-bearing runner graph."""

    vertices = ("skeleton", "free-progression", "forced-owner", "grid-gap")
    raw_order = vertices
    # Switch from logical dependency to the size of the exact certificate.
    switched_order = ("grid-gap", "forced-owner", "free-progression", "skeleton")
    flips = 6
    return (
        vertices,
        raw_order,
        switched_order,
        flips,
        (0, 1, 2, 3),
        0,
        4,
        1,
    )


def main() -> None:
    skeleton = skeleton_audit()
    centers = residue_parity_center_audit()
    inequalities = all_size_inequality_audit()
    grids = endpoint_grid_audit()
    finite = finite_full_pool_replay()
    tournament = tournament_telemetry()
    certificate = sha256(
        repr((skeleton, centers, inequalities, grids, finite, tournament)).encode()
    ).hexdigest()

    print("LRC14 SHELL-FIVE UNIVERSAL ENDPOINT-GRID TEMPLATE OBSTRUCTION")
    print("arithmetic=integer+fractions.Fraction optimizer=none SAT=none floating_point=none sampled_circle=none")
    print("scope=one_unit_numerator_chosen_uniformly_before_the_seven_free_lifts")
    print("not_claimed=absence_of_U_dependent_endpoint_grid_witnesses")
    print()
    print("LOW_SPEED_SKELETON")
    print(f"speed_count={skeleton[0]} speeds={skeleton[1]}")
    print(f"exact_components={skeleton[2]}")
    print(f"no_jump_step_bounds={skeleton[3]} forbidden_width=2/11")
    print("centers_positive_half=(2/13,4/13,9/26)")
    print()
    print("RESIDUE_PARITY_CENTER_VIEW")
    print(f"free_mod26_states={centers[0]}")
    print(f"free_safe_center_numerators_mod26={centers[1]}")
    print(f"after_forced_owner_states={centers[2]}")
    print("interpretation=exact_center_robust_slopes_reduce_to_plus_or_minus_4/13")
    print()
    print("ALL_SIZE_LOCALIZATION")
    print("candidate_pool=all_free_raw_lifts_union_{B-3,B-2,B}")
    print("2_over_13_branch=empty_by_oppositely_signed_forced_B-3_and_B-2")
    print("I4=[4/13-9/(143(B-2)),4/13+9/(143(B-3))]")
    print("I9_class15=[9/26-9/(143(B-10)),9/26-2/(143(B-2))]")
    print("I9_classes37_41=[9/26+2/(143(B-3)),9/26+7/(286(B-18))]")
    print("plus_reflections_under_t_to_1-t")
    print(f"least_class_rows={inequalities[0]}")
    print(f"positive_linear_margins={inequalities[1]}")
    print()
    print("UNIT_ENDPOINT_GRIDS")
    print("q=5d:no_unit_grid_point_in_localized_intervals")
    print("q=13d:no_unit_grid_point_in_localized_intervals")
    print("q13_direct_obstruction=balanced_j_[p]_13_and_u=abs([p_inverse*j]_(13d))")
    print("q13_direct_properties=raw_u_in_{1,12},u<=13d/2<B,phase<=6/(13d)<1/11")
    print(f"grid_audit_k_limit={grids[0]} rows={grids[1]} first={grids[2]} last={grids[3]}")
    print(f"grid_rows_sha256={grids[4]} inverse_rows_sha256={grids[5]}")
    print()
    print("FINITE_FULL_POOL_REPLAY")
    print(f"rows={finite[0]}")
    print(f"rows_sha256={finite[1]}")
    print("scope_guard=finite_rows_are_replay_only_not_the_source_of_all_size_quantifiers")
    print()
    print("TOURNAMENT_ANALYSIS")
    print(f"vertices={tournament[0]}")
    print(f"observable_order={tournament[1]} switch_order={tournament[2]} edge_flips={tournament[3]}")
    print(f"score_histogram={tournament[4]} directed_cycles={tournament[5]} SCCs={tournament[6]} Hamiltonian_paths={tournament[7]}")
    print("tie_Hamiltonian_path=logical_dependency_order")
    print("telemetry_only=binary_order_forgets_signed_phase_margins_and_grid_arithmetic")
    print("challenged_assumption=vertices_are_proof_obligations_not_runners")
    print()
    print("VERDICT")
    print("classes_mod52=(15,37,41):uniform_single_unit_numerator_template_fails_for_q=5d_and_q=13d")
    print("classes_not_closed=U_dependent_or_multi_numerator_divisor_arguments_remain_possible")
    print(f"certificate_sha256={certificate}")
    print("PASS")


if __name__ == "__main__":
    main()
