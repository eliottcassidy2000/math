#!/usr/bin/env python3
"""Exact carry/selector product-extension and recombination audit.

This scratch companion starts after THM-2882.  It classifies the smallest
joint carrier for the independent order-thirteen carry and Boolean
origin-selector parity, then tests every natural quotient against the
recombined two-Prony-channel coefficient.

It proves no physical selector construction, row exclusion, or LRC(14).
"""

from __future__ import annotations

from hashlib import sha256
from itertools import product
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
COMP = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge/results"
THEOREMS = ROOT / "01-canon/theorems"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_hash(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


PINNED = {
    THEOREMS / "THM-2882-event-twisted-all-q-coefficient-carry-lift.md":
        "a0f722ff9982b3011e0db23895b0949c9d6ae61ac837cb096a5aecf5d6882486",
    COMP / "lrc14_event_twisted_all_q_coefficient_carry_lift_thm2882.py":
        "3ed346e0c631b34bd61f0c4d27d7f161e8d35b70decfb95f5207c5f57893d005",
    RESULTS / "lrc14_event_twisted_all_q_coefficient_carry_lift_thm2882.out":
        "0faa0a24f6ba8b6c88b6bbfc4f225e38667097b1a937d977741453499884901c",
    COMP / "lrc14_q0_q3_one_fibre_selector_provenance_obstruction_thm2880.py":
        "7d379279f08f4df8f16d1fc699f4c6f7a9b657fa151e54d2794ca883b2ceee24",
    RESULTS / "lrc14_q0_q3_one_fibre_selector_provenance_obstruction_thm2880.out":
        "7c7a74d38eac97155b55178d08b317b4658f52cff9651e93113f6b556d19f2e1",
    COMP / "lrc14_endpoint_factor_exit_carry_transducer_thm2878.py":
        "b379b9278f6c0d0864908bbc2da2123f4d208eb83c35738d12f651119e7a3366",
    RESULTS / "lrc14_endpoint_factor_exit_carry_transducer_thm2878.out":
        "35bdec6bc5b2bb3c0287bd5aee26c66e8485876e066bf423e2fadb3a94727224",
    COMP / "lrc14_q3_q11_q7_bockstein_holonomy_thm2851.py":
        "2227f59c717095da0f2042096ada145de4e3661530c9aa2cc9020f42c8237a8b",
    RESULTS / "lrc14_q3_q11_q7_bockstein_holonomy_thm2851.out":
        "424fd2e83a618f862a5ee1b5f073a93fe236d10cdc5412eab1b54dec5e537eac",
}
for path, expected in PINNED.items():
    require(lf_hash(path) == expected, f"pinned dependency changed: {path.name}")


P = 13
F2 = (0, 1)
V4 = tuple(product(F2, repeat=2))
FIELD = 53
OMEGA = pow(2, 4, FIELD)
Z = pow(OMEGA, 3, FIELD)


def q_next(q, h):
    return (q + h) % P


def kappa(q, h):
    return (q + h) // P


def truth(q):
    """(zero-origin E3 truth, stepped-origin E3 truth)."""
    return (
        int(q in (0, 3, 11)),
        int(q in (0, 11)),
    )


def truth_parity(q):
    left, right = truth(q)
    return left ^ right


def positive_selector(q):
    """Zero-origin-full, stepped-origin-empty Boolean selector."""
    left, right = truth(q)
    return left, right ^ 1


def selector_parity(q):
    left, right = positive_selector(q)
    return left ^ right


def selector_edge(q, h):
    return selector_parity(q) ^ selector_parity(q_next(q, h))


def truth_edge(q, h):
    source = truth(q)
    target = truth(q_next(q, h))
    return source[0] ^ target[0], source[1] ^ target[1]


def add_v4(left, right):
    return left[0] ^ right[0], left[1] ^ right[1]


def lift169(state, h):
    a, q = state
    return (a + kappa(q, h)) % P, q_next(q, h)


def lift338(state, h):
    a, parity, q = state
    return (
        (a + kappa(q, h)) % P,
        parity ^ selector_edge(q, h),
        q_next(q, h),
    )


def lift676(state, h):
    a, selector0, selector1, q = state
    selector = selector0, selector1
    moved_selector = add_v4(selector, truth_edge(q, h))
    return (
        (a + kappa(q, h)) % P,
        moved_selector[0],
        moved_selector[1],
        q_next(q, h),
    )


def repeat_unit(state, steps, lift):
    result = state
    for _index in range(steps):
        result = lift(result, 1)
    return result


def orbit(state, lift):
    rows = []
    current = state
    while current not in rows:
        rows.append(current)
        current = lift(current, 1)
    require(current == state, "unit orbit returned to the wrong state")
    return tuple(rows)


def rank_mod(matrix, prime):
    work = [[entry % prime for entry in row] for row in matrix]
    rows = len(work)
    columns = len(work[0]) if rows else 0
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (
                row for row in range(pivot_row, rows)
                if work[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        inverse = pow(work[pivot_row][column], -1, prime)
        work[pivot_row] = [
            value * inverse % prime for value in work[pivot_row]
        ]
        for row in range(rows):
            if row == pivot_row or work[row][column] == 0:
                continue
            scalar = work[row][column]
            work[row] = [
                (left - scalar * right) % prime
                for left, right in zip(work[row], work[pivot_row])
            ]
        pivot_row += 1
    return pivot_row


def v4_h2_dimensions():
    nonzero = (1, 2, 3)
    one_index = {value: index for index, value in enumerate(nonzero)}
    two_keys = tuple((left, right) for left in nonzero for right in nonzero)
    two_index = {value: index for index, value in enumerate(two_keys)}

    delta_one = []
    for left, right in two_keys:
        row = [0] * len(nonzero)
        row[one_index[left]] += 1
        row[one_index[right]] += 1
        if left ^ right:
            row[one_index[left ^ right]] -= 1
        delta_one.append(row)

    delta_two = []
    for left in range(4):
        for middle in range(4):
            for right in range(4):
                row = [0] * len(two_keys)
                terms = (
                    (+1, (middle, right)),
                    (-1, (left ^ middle, right)),
                    (+1, (left, middle ^ right)),
                    (-1, (left, middle)),
                )
                for sign, key in terms:
                    if 0 not in key:
                        row[two_index[key]] += sign
                delta_two.append(row)

    b2_rank = rank_mod(delta_one, P)
    z2_dimension = len(two_keys) - rank_mod(delta_two, P)
    return len(nonzero), len(two_keys), b2_rank, z2_dimension


def carry_coboundary_ranks():
    nonzero = tuple(range(1, P))
    index = {value: position for position, value in enumerate(nonzero)}
    delta_one = []
    carry = []
    for h in range(P):
        for k in range(P):
            row = [0] * len(nonzero)
            if h:
                row[index[h]] += 1
            if k:
                row[index[k]] += 1
            total = (h + k) % P
            if total:
                row[index[total]] -= 1
            delta_one.append(row)
            carry.append((h + k) // P)
    rank = rank_mod(delta_one, P)
    augmented = rank_mod(
        [row + [carry[position]] for position, row in enumerate(delta_one)],
        P,
    )
    return rank, augmented


def matrix_apply(matrix, vector):
    return (
        (matrix[0] * vector[0] + matrix[1] * vector[1]) % 2,
        (matrix[2] * vector[0] + matrix[3] * vector[1]) % 2,
    )


def compose_matrix(left, right):
    # left after right
    return (
        (left[0] * right[0] + left[1] * right[2]) % 2,
        (left[0] * right[1] + left[1] * right[3]) % 2,
        (left[2] * right[0] + left[3] * right[2]) % 2,
        (left[2] * right[1] + left[3] * right[3]) % 2,
    )


def matrix_order(matrix):
    identity = (1, 0, 0, 1)
    current = identity
    for order in range(1, 7):
        current = compose_matrix(matrix, current)
        if current == identity:
            return order
    raise RuntimeError("GL2(F2) element order exceeded six")


def main() -> None:
    require(
        pow(OMEGA, P, FIELD) == 1
        and OMEGA != 1
        and pow(Z, P, FIELD) == 1
        and Z != 1,
        "order-thirteen field control changed",
    )

    truth_rows = tuple(truth(q) for q in range(P))
    selector_rows = tuple(positive_selector(q) for q in range(P))
    parity_rows = tuple(selector_parity(q) for q in range(P))
    require(
        truth_rows
        == (
            (1, 1), (0, 0), (0, 0), (1, 0), (0, 0), (0, 0),
            (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (1, 1),
            (0, 0),
        )
        and selector_rows
        == tuple((left, right ^ 1) for left, right in truth_rows)
        and parity_rows == tuple(1 ^ int(q == 3) for q in range(P))
        and tuple(q for q in range(P) if parity_rows[q] != parity_rows[(q + 1) % P])
        == (2, 3),
        "two-origin truth or positive-selector section changed",
    )

    # The selector edge is the exact coboundary delta p.  Its two unit
    # edges are disjoint from the unique carry unit edge q12->q0.
    selector_cocycle_checks = 0
    carry_composition_checks = 0
    joint_composition_checks = 0
    for q in range(P):
        for h in range(P):
            for k in range(P):
                qh = q_next(q, h)
                reduced = (h + k) % P
                central = (h + k) // P
                require(
                    selector_edge(q, h) ^ selector_edge(qh, k)
                    == selector_edge(q, reduced),
                    "selector edge stopped being an exact F2 coboundary",
                )
                selector_cocycle_checks += 1
                require(
                    kappa(q, h) + kappa(qh, k)
                    == kappa(q, reduced) + central,
                    "carry composition changed",
                )
                carry_composition_checks += 1
                for a in range(P):
                    for parity in F2:
                        left = lift338(lift338((a, parity, q), h), k)
                        direct = lift338((a, parity, q), reduced)
                        direct = (
                            (direct[0] + central) % P,
                            direct[1],
                            direct[2],
                        )
                        require(left == direct, "joint lift composition changed")
                        joint_composition_checks += 1
    require(
        selector_cocycle_checks == P**3
        and carry_composition_checks == P**3
        and joint_composition_checks == 2 * P**4,
        "composition census changed",
    )

    unit_selector_edges = tuple(q for q in range(P) if selector_edge(q, 1))
    unit_carry_edges = tuple(q for q in range(P) if kappa(q, 1))
    require(
        unit_selector_edges == (2, 3)
        and unit_carry_edges == (12,)
        and selector_edge(0, 3) == 1
        and kappa(0, 3) == 0
        and selector_edge(12, 1) == 0
        and kappa(12, 1) == 1,
        "selector/carry independence witnesses changed",
    )

    # Gauge c=b XOR p(q) trivializes the F2 connection.
    gauge_checks = 0
    for a in range(P):
        for parity in F2:
            for q in range(P):
                for h in range(P):
                    moved = lift338((a, parity, q), h)
                    require(
                        parity ^ selector_parity(q)
                        == moved[1] ^ selector_parity(moved[2]),
                        "selector parity gauge ceased to be constant",
                    )
                    gauge_checks += 1
    require(gauge_checks == 2 * P**3, "selector gauge census changed")

    parity_orbits = tuple(
        orbit((0, parity, 0), lift338) for parity in F2
    )
    require(
        all(len(rows) == 169 for rows in parity_orbits)
        and not set(parity_orbits[0]) & set(parity_orbits[1])
        and set(parity_orbits[0]) | set(parity_orbits[1])
        == {
            (a, parity, q)
            for a in range(P) for parity in F2 for q in range(P)
        },
        "338-state unit action stopped being two C169 cycles",
    )

    def parity_flip(state):
        a, parity, q = state
        return a, parity ^ 1, q

    commute_checks = 0
    for state in set(parity_orbits[0]) | set(parity_orbits[1]):
        require(
            parity_flip(lift338(state, 1))
            == lift338(parity_flip(state), 1),
            "selector flip stopped commuting with positive carry",
        )
        commute_checks += 1
    generated338 = {
        parity_flip(repeat_unit((0, 0, 0), n, lift338)) if bit
        else repeat_unit((0, 0, 0), n, lift338)
        for n in range(169) for bit in F2
    }
    require(
        commute_checks == 338 and len(generated338) == 338,
        "C169 x C2 regular action changed",
    )

    # Full two-origin selector transport is V4, not only its parity
    # quotient.  Regrading by truth makes mismatch x=s XOR t constant.
    selector_transport_checks = 0
    for a in range(P):
        for selector in V4:
            for q in range(P):
                for h in range(P):
                    moved = lift676((a, selector[0], selector[1], q), h)
                    source_mismatch = add_v4(selector, truth(q))
                    target_selector = moved[1], moved[2]
                    target_mismatch = add_v4(target_selector, truth(moved[3]))
                    require(
                        source_mismatch == target_mismatch,
                        "full selector mismatch gauge changed",
                    )
                    selector_transport_checks += 1
    require(
        selector_transport_checks == 4 * P**3,
        "full-selector transport census changed",
    )
    positive_mismatch = tuple(
        add_v4(positive_selector(q), truth(q)) for q in range(P)
    )
    require(
        positive_mismatch == ((0, 1),) * P,
        "positive selector stopped being a constant mismatch section",
    )

    full_orbits = tuple(
        orbit((0, selector[0], selector[1], 0), lift676)
        for selector in V4
    )
    require(
        all(len(rows) == 169 for rows in full_orbits)
        and len(set().union(*(set(rows) for rows in full_orbits))) == 676,
        "full-selector carrier stopped being four C169 cycles",
    )

    # Selector quotient anatomy.  Parity preserves support but identifies
    # the + and - signed currents.  The three-value quotient preserves the
    # scalar value but is not a congruence for the V4 translation action.
    mismatch_amplitude = {
        mismatch: mismatch[1] - mismatch[0] for mismatch in V4
    }
    mismatch_parity = {
        mismatch: mismatch[0] ^ mismatch[1] for mismatch in V4
    }
    parity_classes = tuple(
        tuple(mismatch for mismatch in V4 if mismatch_parity[mismatch] == bit)
        for bit in F2
    )
    require(
        parity_classes == (
            ((0, 0), (1, 1)),
            ((0, 1), (1, 0)),
        )
        and all(
            bool(mismatch_amplitude[mismatch])
            == bool(mismatch_parity[mismatch])
            for mismatch in V4
        )
        and {
            mismatch_amplitude[mismatch] for mismatch in parity_classes[1]
        } == {-1, 1},
        "V4 parity quotient anatomy changed",
    )
    three_value_failure_witness = None
    for left in V4:
        for right in V4:
            if mismatch_amplitude[left] != mismatch_amplitude[right]:
                continue
            for shift in V4:
                moved_left = add_v4(left, shift)
                moved_right = add_v4(right, shift)
                if (
                    mismatch_amplitude[moved_left]
                    != mismatch_amplitude[moved_right]
                ):
                    three_value_failure_witness = (
                        left, right, shift, moved_left, moved_right
                    )
                    break
            if three_value_failure_witness is not None:
                break
        if three_value_failure_witness is not None:
            break
    require(
        three_value_failure_witness
        == ((0, 0), (1, 1), (0, 1), (0, 1), (1, 0)),
        "three-value selector quotient acquired a V4 action",
    )

    # Group-action classification.  C169 cannot act nontrivially on V4
    # because GL2(F2)=S3 has no element of order 13.  A C2 action on C169
    # is identity or inversion; inversion reverses the positive generator
    # and sends exceptional truth q3 to ordinary q10.
    gl2 = tuple(
        matrix
        for matrix in product(F2, repeat=4)
        if (matrix[0] * matrix[3] - matrix[1] * matrix[2]) % 2 == 1
    )
    gl2_orders = tuple(sorted(matrix_order(matrix) for matrix in gl2))
    c169_eligible_actions = tuple(
        matrix for matrix in gl2 if 169 % matrix_order(matrix) == 0
    )
    c169_involutive_units = tuple(
        unit for unit in range(169)
        if gcd(unit, 169) == 1 and unit * unit % 169 == 1
    )
    v4_to_inversion_characters = tuple(
        (
            alpha,
            beta,
            tuple(
                168 if (alpha * left + beta * right) % 2 else 1
                for left, right in V4
            ),
        )
        for alpha, beta in V4
    )
    orientation_preserving_v4_actions = tuple(
        character
        for character in v4_to_inversion_characters
        if set(character[2]) == {1}
    )
    inversion_truth_witness = (
        3,
        (-3) % P,
        truth(3),
        truth((-3) % P),
    )
    require(
        len(gl2) == 6
        and gl2_orders == (1, 2, 2, 2, 3, 3)
        and c169_eligible_actions == ((1, 0, 0, 1),)
        and c169_involutive_units == (1, 168)
        and len(v4_to_inversion_characters) == 4
        and len(orientation_preserving_v4_actions) == 1
        and inversion_truth_witness == (3, 10, (1, 0), (0, 0)),
        "semidirect-action classification changed",
    )

    # V4 carries no intrinsic F13-valued H2 class, whereas the C13 carry
    # cocycle is not a coboundary.
    v4_dimensions = v4_h2_dimensions()
    carry_ranks = carry_coboundary_ranks()
    require(
        v4_dimensions == (3, 9, 3, 3)
        and carry_ranks == (11, 12),
        "V4/carry cohomology separation changed",
    )

    # Recombination no-go.  On the vertical carry T, selector and truth are
    # fixed while (U,V) transforms by diag(Z,1).  The sum map R=(1,1)
    # cannot intertwine this with any scalar.  The selector C2 cannot
    # cancel an order-thirteen character.
    scalar_intertwiners = tuple(
        scalar for scalar in range(FIELD)
        if (Z, 1) == (scalar, scalar)
    )
    mu2 = tuple(value for value in range(1, FIELD) if value * value % FIELD == 1)
    mu13 = tuple(value for value in range(1, FIELD) if pow(value, P, FIELD) == 1)
    require(
        scalar_intertwiners == ()
        and mu2 == (1, FIELD - 1)
        and set(mu2) & set(mu13) == {1},
        "selector character unexpectedly cancelled the carry character",
    )
    kernel_vector = (1, FIELD - 1)
    moved_kernel = (Z, FIELD - 1)
    kernel_invariant = (
        kernel_vector[0] * moved_kernel[1]
        - kernel_vector[1] * moved_kernel[0]
    ) % FIELD == 0
    require(not kernel_invariant, "recombination kernel became T-invariant")

    recombined_failure_count = 0
    recombined_survivor_count = 0
    for q in range(P):
        for h in range(P):
            carry = kappa(q, h)
            edge = pow(Z, h + carry, FIELD)
            source_r = (q - 3) % P
            target_r = (carry + q_next(q, h) - 3) % P
            source_u = pow(Z, source_r, FIELD)
            target_u = pow(Z, target_r, FIELD)
            require(
                target_u == edge * source_u % FIELD,
                "charged edge model changed",
            )
            source_s = (source_u + 1) % FIELD
            target_s = (target_u + 1) % FIELD
            if target_s == edge * source_s % FIELD:
                recombined_survivor_count += 1
            else:
                recombined_failure_count += 1
    require(
        (recombined_failure_count, recombined_survivor_count) == (144, 25),
        "recombined edge hostile changed",
    )

    smallest_fixed_section = P**2
    smallest_independent_parity = P**2 * 2
    smallest_full_selector_torsor = P**2 * 4
    require(
        (
            smallest_fixed_section,
            smallest_independent_parity,
            smallest_full_selector_torsor,
        ) == (169, 338, 676),
        "minimal joint state invoices changed",
    )

    thm2880_output = (
        RESULTS / "lrc14_q0_q3_one_fibre_selector_provenance_obstruction_thm2880.out"
    ).read_text(encoding="utf-8").replace("\r\n", "\n")
    thm2882_output = (
        RESULTS / "lrc14_event_twisted_all_q_coefficient_carry_lift_thm2882.out"
    ).read_text(encoding="utf-8").replace("\r\n", "\n")
    require(
        "unique_D_to_S=q12_to_q0;q0_to_q3_events=0;"
        "uniform_over_169_addresses=1" in thm2880_output
        and "q0_to_q3_changes_selector_parity_but_kappa(0,3)=0"
        in thm2882_output
        and "two-node_pair_transport=diag(E,1)" in thm2882_output,
        "pinned selector/carry/current boundary changed",
    )

    print("CARRY / SELECTOR PRODUCT-EXTENSION AND QUOTIENT AUDIT")
    print(
        f"truth_rows={truth_rows}; positive_selector={selector_rows}; "
        f"selector_parity={parity_rows}; selector_unit_edges={unit_selector_edges}; "
        f"carry_unit_edges={unit_carry_edges}"
    )
    print(
        "independence_witnesses=q0--3->q3:(selector_delta=1,kappa=0),"
        "q12--1->q0:(selector_delta=0,kappa=1); "
        "selector_edge=delta(p)_is_exact; carry_is_nonsplit"
    )
    print(
        f"composition_checks=(selector={selector_cocycle_checks},"
        f"carry={carry_composition_checks},joint={joint_composition_checks}); "
        f"selector_gauge_checks={gauge_checks}"
    )
    print(
        "fixed_positive_section=truth_mismatch_(0,1)_constant; "
        "no_extra_state_if_this_external_section_is_frozen; "
        f"carrier_size={smallest_fixed_section}"
    )
    print(
        "independent_selector_parity_lift=C169xC2_direct_product_torsor; "
        f"carrier_size={smallest_independent_parity}; "
        "unit_action=two_disjoint_C169_cycles; adjoining_free_flip_is_regular"
    )
    print(
        "full_two_origin_selector_lift=C169xV4_direct_product_torsor; "
        f"carrier_size={smallest_full_selector_torsor}; "
        "unit_action=four_disjoint_C169_cycles; mismatch=s_XOR_truth_is_constant"
    )
    print(
        f"semidirect_classification=GL2F2_orders={gl2_orders}; "
        f"C169_to_AutV4_eligible={len(c169_eligible_actions)}_identity_only; "
        f"C2_to_AutC169_involutions={c169_involutive_units}; "
        "V4_to_AutC169_actions=4_with_only_1_orientation-preserving; "
        "inversion_reverses_positive_orientation_and_moves_truth_q3_to_q10; "
        "therefore_joint_action_is_direct_not_dihedral"
    )
    print(
        f"cohomology=V4_(C1,C2,rankB2,dimZ2)={v4_dimensions},H2=0; "
        f"C13_carry_(coboundary_rank,augmented_rank)={carry_ranks}; "
        "selector_F2_coboundary_cannot_host_Bockstein"
    )
    print(
        f"parity_quotient_classes={parity_classes}; "
        "preserves_nonzero_support_but_identifies_amplitudes_+1_and_-1; "
        f"three_value_quotient_not_V4_equivariant_witness="
        f"{three_value_failure_witness}"
    )
    print(
        f"vertical_T=diag(omega^3,1); scalar_sum_intertwiners="
        f"{scalar_intertwiners}; mu2_intersection_mu13={set(mu2) & set(mu13)}; "
        "any_quotient_retaining_T_has_no_scalar_action_on_U+V"
    )
    print(
        f"recombined_edges=(fail={recombined_failure_count},"
        f"survive={recombined_survivor_count}); "
        "a_quotient_killing_T_allows_recombination_only_by_losing_the_"
        "character-three_carry; split_pair_or_charged_U_remains_lawful"
    )
    print(
        "THM2880_connection=q0_to_q3_selector_provenance_changes_with_"
        "zero_carry; independent-origin_F2_is_exactly_its_unclosed_"
        "broader_signed_coupling_sidecar"
    )
    print(
        "VERDICT=smallest_faithful_joint_parity_extension_is_338-state_"
        "C169xC2_direct_product_groupoid/torsor; full_signed_selector_"
        "transport_requires_676-state_C169xV4; no selector quotient_"
        "both_preserves_sign_and_descends_the_recombined_current_while_"
        "retaining_Bockstein_T"
    )
    print(
        "scope=exact abstract residue/selector/coefficient representation; "
        "no constructed physical origin coupling, QA/QAB attachment, row "
        "exclusion, or LRC14 conclusion"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
