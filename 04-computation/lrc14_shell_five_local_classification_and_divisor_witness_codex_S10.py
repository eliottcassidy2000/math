#!/usr/bin/env python3
"""Exact shell-five continuation probe beyond THM-836.

THM-836 leaves the local owner shell

    s = 2B - 13d = 5

because its two coefficient-nine owner inequalities are feasible at
``(d,B)=(11,74)``.  This script separates three logically different facts.

1. An elementary all-size argument classifies *every* locally feasible
   shell-five owner row.  Necessarily

       d mod 52 in {11,15,37,41},
       missing signed class = {3,10},
       owner speeds = {B-3,B-2}.

   Conversely those owner speeds satisfy the two local THM-836 inequalities
   in each of the four congruence classes.  This is only a classification of
   the local owner shell, not of ten-speed cores.

2. The class ``d=11 mod 52`` is uniformly impossible after adjoining the
   cheapest exception-divisor obligation, the raw grid of denominator
   ``q=5d``.  The explicit unit numerator

       p=(45d-1)/26

   lies in ``E_U`` for every signed-complement lift allowed by the shell-five
   classification, while ``Q_d(p/(5d))=2/5<11/13``.  This contradicts even
   ``E_U subset H_d``, hence also THM-836's guarded containment.

3. A finite exact census at the next class ``d=15`` joins the two endpoint
   grids ``q=5d,13d`` to THM-772 divisor completeness and THM-803 parity
   support.  The conjunction is empty at this one height.  This finite row is
   reconnaissance, not an extrapolated uniform theorem.

All verdicts use integers or fractions.Fraction.  No optimizer, SAT solver,
floating point, or sampled circle is used.  Tournament Analysis is telemetry
only: the proof carrier is the incidence among unit-grid numerators, signed
residue/lift choices, and divisor/parity obligations.
"""

from __future__ import annotations

from fractions import Fraction as F
from functools import reduce
from hashlib import sha256
from itertools import combinations, product
from math import gcd


THRESHOLD = F(11, 13)
LOCAL_CLASSES_MOD_52 = (11, 15, 37, 41)
SIGNED_COMPLEMENT = (1, 2, 4, 5, 6, 7, 8, 9, 11, 12)
FREE_RESIDUES = (1, 2, 4, 5, 8, 11, 12)


def balanced(value: int, modulus: int) -> int:
    residue = value % modulus
    if 2 * residue > modulus:
        residue -= modulus
    return residue


def least_absolute(value: int, modulus: int) -> int:
    return abs(balanced(value, modulus))


def norm(value: F) -> F:
    residue = value % 1
    return min(residue, 1 - residue)


def shell_B(d: int) -> int:
    assert d > 0 and d % 2 == 1
    return (13 * d + 5) // 2


def packing_width(d: int) -> F:
    """G=B-L in THM-836 for s=5."""

    return F(55 * (13 * d + 5), 2 * (117 * d + 55))


def largest_with_residue(B: int, residue: int) -> int:
    residue %= 13
    assert residue != 0
    return B - ((B - residue) % 13)


def locally_feasible(d: int) -> tuple[bool, int, int]:
    """Whether the two coefficient-nine owner inequalities can both hold."""

    B = shell_B(d)
    denominator = 117 * d + 55
    threshold_numerator = 117 * d * B
    left = largest_with_residue(B, 3 * d)
    right = largest_with_residue(B, -3 * d)
    feasible = (
        left * denominator >= threshold_numerator
        and right * denominator >= threshold_numerator
    )
    return feasible, min(left, right), max(left, right)


def local_classification_audit(limit: int = 9_999) -> tuple[object, ...]:
    # Symbolic comparisons used in the proof:
    # G<4 because 4(234d+110)-(715d+275)=221d+165>0.
    # G>3 exactly when 13d>55; every feasible congruence has d>=11.
    for d in (1, 3, 5, 11, 15, 37, 41):
        G = packing_width(d)
        assert G < 4
        assert (G > 3) == (13 * d > 55)

    # B=9 mod 13.  The only signed pair wholly contained in the top four
    # classes {9,8,7,6} is {6,7}; solving +/-3d={6,7} gives d=+/-2 mod 13.
    top_classes = {9, 8, 7, 6}
    signed_pairs = {
        frozenset((r, (-r) % 13))
        for r in range(1, 13)
    }
    pairs_inside = tuple(sorted(tuple(sorted(pair)) for pair in signed_pairs if pair <= top_classes))
    assert pairs_inside == ((6, 7),)
    assert tuple(r for r in range(52) if r % 2 and r % 13 in (2, 11)) == LOCAL_CLASSES_MOD_52

    rows = []
    for d in range(1, limit + 1, 2):
        if d % 13 == 0:
            continue
        feasible, owner_lo, owner_hi = locally_feasible(d)
        predicted = d % 52 in LOCAL_CLASSES_MOD_52
        assert feasible == predicted
        if feasible:
            B = shell_B(d)
            assert (owner_lo, owner_hi) == (B - 3, B - 2)
            assert {owner_lo % 13, owner_hi % 13} == {6, 7}
            assert {(5 * d) % 13, (-5 * d) % 13} == {3, 10}
            rows.append((d, B, owner_lo, owner_hi, packing_width(d)))

    return (
        limit,
        len(rows),
        tuple(sorted({d % 52 for d, *_rest in rows})),
        rows[0],
        rows[-1],
        sha256(repr(rows).encode()).hexdigest(),
    )


def residue_parity_phase_table() -> tuple[tuple[int, str, int], ...]:
    rows = []
    for residue in range(1, 13):
        for parity in (0, 1):
            lift = next(u for u in range(1, 27) if u % 13 == residue and u % 2 == parity)
            rows.append((residue, "odd" if parity else "even", balanced(9 * lift, 26)))
    return tuple(rows)


def uniform_11_mod_52_witness_audit(limit_k: int = 99) -> tuple[object, ...]:
    table = residue_parity_phase_table()
    phase = {(residue, parity): value for residue, parity, value in table}

    # In d=11 mod 52, B is even.  The seven free signed classes have base
    # phase a/26 with either a>=4 or a<=-3.  The three forced classes are
    # B-3 (residue 6, odd, a=-11), B-2 (residue 7, even, a=-2), and
    # B (residue 9, even, a=-10).
    free_values = tuple(
        (residue, phase[(residue, parity)])
        for residue in FREE_RESIDUES
        for parity in ("even", "odd")
    )
    assert all(value >= 4 or value <= -3 for _residue, value in free_values)
    forced_values = (
        phase[(6, "odd")],
        phase[(7, "even")],
        phase[(9, "even")],
    )
    assert forced_values == (-11, -2, -10)

    # The following integer differences prove the three phase lower bounds:
    #   a>=4:  (27d-5)/(260d) - 1/11 = (37d-55)/(2860d)>0;
    #   a<=-3: 3/26 > 1/11, with no half-circle wrap because B<10d;
    #   a=-2:  (33d+1)/(260d) - 1/11 = (103d+11)/(2860d)>0.
    rows = []
    for k in range(limit_k + 1):
        d = 52 * k + 11
        B = shell_B(d)
        p = (45 * d - 1) // 26
        q = 5 * d
        assert 26 * p == 45 * d - 1
        assert gcd(p, q) == 1
        assert p % 5 == 4
        assert B % 2 == 0
        assert B < 10 * d
        assert 37 * d - 55 > 0
        assert 103 * d + 11 > 0

        t = F(p, q)
        assert norm(9 * d * t) + norm(4 * d * t) == F(2, 5) < THRESHOLD

        # A finite direct replay of every allowed lift at these audit rows.
        forced = {6: B - 3, 7: B - 2, 9: B}
        for residue in SIGNED_COMPLEMENT:
            candidates = (
                (forced[residue],)
                if residue in forced
                else tuple(range(residue, B + 1, 13))
            )
            assert all(norm(u * t) >= F(1, 11) for u in candidates)
        rows.append((d, B, p, q))

    return (
        limit_k,
        len(rows),
        table,
        rows[0],
        rows[-1],
        sha256(repr(rows).encode()).hexdigest(),
    )


def parity_twisted_class(speed: int) -> int:
    residue = (speed if speed % 2 else 7 * speed) % 13
    if residue == 0:
        return 0
    return min(residue, 13 - residue)


def folded_good(d: int, numerator: int, denominator: int) -> bool:
    total = least_absolute(9 * d * numerator, denominator)
    total += least_absolute(4 * d * numerator, denominator)
    return 13 * total >= 11 * denominator


def grid_data(d: int, denominator: int, B: int) -> tuple[tuple[int, ...], int, dict[int, int]]:
    forbidden = tuple(
        p
        for p in range(1, denominator)
        if gcd(p, denominator) == 1 and not folded_good(d, p, denominator)
    )
    full = (1 << len(forbidden)) - 1
    masks = {
        speed: sum(
            1 << index
            for index, p in enumerate(forbidden)
            if 11 * least_absolute(speed * p, denominator) < denominator
        )
        for speed in range(1, B + 1)
    }
    return forbidden, full, masks


def d15_endpoint_grid_census() -> tuple[object, ...]:
    """Pure exhaustive finite census of the first remaining congruence class."""

    d = 15
    B = shell_B(d)
    options = {residue: list(range(residue, B + 1, 13)) for residue in SIGNED_COMPLEMENT}
    options[6] = [B - 3]
    options[7] = [B - 2]
    options[9] = [B]
    fixed = (B - 3, B - 2, B)
    assert fixed == (97, 98, 100)

    divisor_full = (1 << 11) - 1
    divisor_mask = {
        speed: sum(1 << (modulus - 2) for modulus in range(2, 13) if speed % modulus == 0)
        for speed in range(1, B + 1)
    }
    parity_full = (1 << 6) - 1
    parity_mask = {
        speed: (
            0
            if parity_twisted_class(speed) == 0
            else 1 << (parity_twisted_class(speed) - 1)
        )
        for speed in range(1, B + 1)
    }
    grids = {
        denominator: grid_data(d, denominator, B)
        for denominator in (5 * d, 13 * d)
    }

    total = divisor_complete = structural = 0
    pass_5d = pass_13d = pass_both = 0
    pass_13d_rows = []
    hasher = sha256()
    for values in product(*(options[residue] for residue in FREE_RESIDUES)):
        total += 1
        speeds = values + fixed
        divisors = reduce(int.__or__, (divisor_mask[u] for u in speeds), 0)
        if divisors != divisor_full:
            continue
        divisor_complete += 1
        parities = reduce(int.__or__, (parity_mask[u] for u in speeds), 0)
        if parities != parity_full:
            continue
        structural += 1

        verdicts = {}
        for denominator, (_forbidden, full, masks) in grids.items():
            covered = reduce(int.__or__, (masks[u] for u in speeds), 0)
            verdicts[denominator] = covered == full
        pass_5d += int(verdicts[5 * d])
        pass_13d += int(verdicts[13 * d])
        pass_both += int(all(verdicts.values()))
        if verdicts[13 * d]:
            row = tuple(sorted(speeds))
            pass_13d_rows.append(row)
            hasher.update(repr(row).encode())

    assert total == 1_605_632
    assert divisor_complete == 121_352
    assert structural == 71_644
    assert pass_5d == 3_004
    assert pass_13d == 4
    assert pass_both == 0
    assert len(pass_13d_rows) == 4
    forbidden_5d, _full_5d, masks_5d = grids[5 * d]
    q5d_uncovered = []
    for row in pass_13d_rows:
        covered = reduce(int.__or__, (masks_5d[u] for u in row), 0)
        missed = tuple(
            p
            for index, p in enumerate(forbidden_5d)
            if not (covered >> index) & 1
        )
        assert missed
        q5d_uncovered.append((row, missed))
    return (
        d,
        B,
        total,
        divisor_complete,
        structural,
        pass_5d,
        pass_13d,
        pass_both,
        tuple(pass_13d_rows),
        tuple(q5d_uncovered),
        tuple((denominator, len(data[0])) for denominator, data in sorted(grids.items())),
        hasher.hexdigest(),
    )


def total_order_edges(order: tuple[int, ...]) -> set[tuple[int, int]]:
    rank = {vertex: index for index, vertex in enumerate(order)}
    return {
        (left, right) if rank[left] < rank[right] else (right, left)
        for left, right in combinations(order, 2)
    }


def tournament_telemetry() -> tuple[object, ...]:
    vertices = LOCAL_CLASSES_MOD_52
    local_slack = {d: packing_width(d) - 3 for d in vertices}
    raw_order = tuple(sorted(vertices, key=lambda d: (local_slack[d], d)))
    # Switch to proof-search priority: the uniformly eliminated class first,
    # followed by decreasing local slack and the numerical tie path.
    switched_order = tuple(sorted(vertices, key=lambda d: (d != 11, -local_slack[d], d)))
    raw_edges = total_order_edges(raw_order)
    switched_edges = total_order_edges(switched_order)
    flips = len(raw_edges.symmetric_difference(switched_edges)) // 2
    assert raw_order == (11, 15, 37, 41)
    assert switched_order == (11, 41, 37, 15)
    assert flips == 3
    return (
        vertices,
        tuple((d, local_slack[d]) for d in vertices),
        raw_order,
        switched_order,
        flips,
        (0, 1, 2, 3),
        0,
        4,
        1,
    )


def main() -> None:
    local = local_classification_audit()
    witness = uniform_11_mod_52_witness_audit()
    d15 = d15_endpoint_grid_census()
    tournament = tournament_telemetry()
    certificate = sha256(repr((local, witness, d15, tournament)).encode()).hexdigest()

    print("LRC14 SHELL-FIVE LOCAL CLASSIFICATION AND DIVISOR WITNESS")
    print("arithmetic=integer+fractions.Fraction optimizer=none SAT=none floating_point=none sampled_circle=none")
    print("base_hypotheses=THM-836 exact ten-speed signed complement plus guarded containment")
    print()
    print("ALL_SIZE_LOCAL_CLASSIFICATION")
    print("s=5 B=(13d+5)/2 B_mod13=9 packing_width_G=55(13d+5)/(2(117d+55))<4")
    print("top_four_classes={9,8,7,6} unique_signed_pair={6,7}")
    print("local_feasible_iff=d_mod13_in_{2,11}_and_d_odd")
    print(f"equivalently_d_mod52={LOCAL_CLASSES_MOD_52}")
    print("missing_signed_class={3,10} owner_speeds={B-3,B-2}")
    print(f"finite_audit_odd_d_limit={local[0]} feasible_rows={local[1]} residues={local[2]}")
    print(f"first_row={local[3]} last_row={local[4]} audit_sha256={local[5]}")
    print()
    print("UNIFORM_DIVISOR_EXCLUSION")
    print("class=d_mod52_11 denominator=q=5d numerator=p=(45d-1)/26")
    print("p_is_unit=26p=45d-1_and_p_mod5=4")
    print("all_allowed_core_phases_at_least=1/11 folded_target_value=2/5<11/13")
    print("consequence=E_U_not_subset_H_d_hence_guarded_containment_fails")
    print("remaining_local_classes_mod52=(15,37,41)")
    print(f"witness_audit_k_limit={witness[0]} rows={witness[1]} first={witness[3]} last={witness[4]}")
    print(f"witness_audit_sha256={witness[5]}")
    print()
    print("D15_FINITE_ENDPOINT_GRID_CENSUS")
    print(f"d={d15[0]} B={d15[1]} signed_lift_rows={d15[2]}")
    print(f"divisor_complete={d15[3]} plus_parity_support={d15[4]}")
    print(f"pass_q_5d={d15[5]} pass_q_13d={d15[6]} pass_both={d15[7]}")
    print(f"q13d_survivors={d15[8]}")
    print(f"q5d_uncovered_on_q13d_survivors={d15[9]}")
    print(f"forbidden_unit_counts={d15[10]} survivor_sha256={d15[11]}")
    print("scope=d15_only_no_uniform_extrapolation")
    print()
    print("TOURNAMENT_ANALYSIS")
    print(f"vertices_mod52={tournament[0]} local_slacks={tournament[1]}")
    print(f"observable_order={tournament[2]} switch_order={tournament[3]} edge_flips={tournament[4]}")
    print(f"score_histogram={tournament[5]} directed_cycles={tournament[6]} SCCs={tournament[7]} Hamiltonian_paths={tournament[8]}")
    print("tie_Hamiltonian_path=numerical_mod52_order")
    print("telemetry_only=both_orders_destroy_the_unit_numerator_by_signed_lift_incidence")
    print("proof_carrier=(divisor_grid,numerator)_to_(signed_class,lift)_to_(divisor,parity_obligation)")
    print("challenged_assumption=vertices_need_not_be_runners_or_target_components")
    print()
    print("SCOPE_GUARDRAIL")
    print("proved=exact_local_s5_classification_plus_uniform_elimination_of_d_mod52_11")
    print("not_proved=existence_or_nonexistence_of_full_cores_in_remaining_classes_15_37_41")
    print("local_survivor_means_owner_inequalities_only_not_guarded_containment")
    print(f"certificate_sha256={certificate}")
    print("PASS")


if __name__ == "__main__":
    main()
