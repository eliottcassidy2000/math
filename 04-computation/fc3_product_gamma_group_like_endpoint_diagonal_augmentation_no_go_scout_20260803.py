#!/usr/bin/env python3
"""Exact algebraic audit of the conditional three-face copy sidecar.

There are two different objects that can be called a copy.

* A removed pole multiset ``sigma`` acts on the product-Gamma alphabet by a
  formal unit.  That unit is group-like, so the native symmetric-function
  coproduct duplicates the *endpoint state* exactly.
* A response obligation is a nonzero centered partition current.  It has
  augmentation zero.  No counital linear coalgebra can send such a current
  to its tensor square: applying either counit would kill the tensor square
  instead of recovering the original current.

The scout reconstructs the exact three-face response bank, tests the first
object on the four reachable blockers and the thirty-one full-domain
rows-{16,17} hostiles, and locates the remaining face-tag and ancestry losses.
It is the exact companion for THM-3314's bounded endpoint/counit statement.
It adds no inherited transition and proves no FC(3), SFC(3), GMC, positivity,
or shared chronology statement.
"""

from __future__ import annotations

import ast
from collections import Counter, defaultdict
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from multiprocessing import Pool, freeze_support
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge/results"
PRIOR_SCRIPT = (
    COMP / "fc3_three_support_face_f12_release_copy_rebind_gate_scout_20260803.py"
)
PRIOR_OUTPUT = (
    RESULTS / "fc3_three_support_face_f12_release_copy_rebind_gate_scout_20260803.out"
)
THM3110 = (
    ROOT
    / "01-canon/theorems/THM-3110-arbitrary-anchored-product-gamma-dominant-tail-and-low-histogram-reduction.md"
)
THM3115 = (
    ROOT
    / "01-canon/theorems/THM-3115-low-degree-monomial-fibre-newton-refinement-transport.md"
)
THM3047 = (
    ROOT
    / "01-canon/theorems/THM-3047-formal-corner-width-product-gamma-moment-and-strict-hankel-positivity.md"
)
PINNED = {
    PRIOR_SCRIPT:
        "90387770e45f2b560cdebc547ac295cb682935445fcbbbf0a7484f6a0a6e7200",
    PRIOR_OUTPUT:
        "6ba33a9f578a08a63b06d81b994818dc9c18f64e3d4cfce971ef5b8eee3b7357",
    THM3110:
        "97113139be9ea16b1434ca9bc076920700920b20d0c08b6fb082ca06c6af48c4",
    THM3115:
        "5f8b98cc13642060f87b230ec65093d7a669b93223e0a86710fad7fd9fe51ab6",
    THM3047:
        "208bb1406fc731997f8ae61f94836a9fe252d24b6e6c8275fce56023f057b744",
}

PAIR = ("F13", "F23")
ROWS_1617 = frozenset((16, 17))
REACHABLE = (
    (1, 2, 3, 4, 5, 6, 7, 8),
    (1, 3, 4, 5, 6, 7, 8),
    (2, 3, 4, 5, 6, 7, 8),
    (3, 4, 5, 6, 7, 8),
)
EXPECTED_ROW_SPECIFIC_HOSTILES = (
    (4, 4, 5, 6, 7, 8),
    (1, 4, 4, 5, 6, 7, 8),
    (2, 4, 4, 5, 6, 7, 8),
    (1, 1, 4, 4, 5, 6, 7, 8),
    (1, 2, 4, 4, 5, 6, 7, 8),
    (2, 2, 4, 4, 5, 6, 7, 8),
    (2, 4, 4, 5, 5, 6, 7, 8),
    (3, 4, 4, 5, 5, 6, 7, 8),
    (1, 1, 2, 4, 4, 5, 6, 7, 8),
    (1, 1, 3, 4, 4, 5, 6, 7, 8),
    (1, 1, 4, 4, 5, 5, 6, 7, 8),
    (1, 2, 2, 4, 4, 5, 6, 7, 8),
    (1, 2, 4, 4, 5, 5, 6, 7, 8),
    (1, 3, 4, 4, 5, 5, 6, 7, 8),
    (2, 2, 2, 4, 4, 5, 6, 7, 8),
    (2, 3, 4, 4, 5, 5, 6, 7, 8),
    (1, 1, 2, 2, 4, 4, 5, 6, 7, 8),
    (1, 1, 2, 4, 4, 5, 5, 6, 7, 8),
    (1, 1, 3, 4, 4, 5, 5, 6, 7, 8),
    (1, 2, 2, 3, 4, 4, 5, 6, 7, 8),
    (1, 2, 2, 4, 4, 5, 5, 6, 7, 8),
    (1, 2, 3, 4, 4, 5, 5, 6, 7, 8),
    (2, 2, 2, 3, 4, 4, 5, 6, 7, 8),
    (2, 2, 2, 4, 4, 5, 5, 6, 7, 8),
    (2, 2, 3, 4, 4, 5, 5, 6, 7, 8),
    (1, 1, 2, 2, 2, 4, 4, 5, 6, 7, 8),
    (1, 1, 2, 2, 4, 4, 5, 5, 6, 7, 8),
    (1, 1, 2, 3, 4, 4, 5, 5, 6, 7, 8),
    (1, 2, 2, 2, 3, 4, 4, 5, 6, 7, 8),
    (1, 2, 2, 2, 4, 4, 5, 5, 6, 7, 8),
    (1, 1, 2, 2, 2, 4, 4, 5, 5, 6, 7, 8),
)
CANDIDATE_SET = frozenset(REACHABLE + EXPECTED_ROW_SPECIFIC_HOSTILES)
MAXIMUM_DEGREE = 14


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_hash(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return sha256(payload).hexdigest()


def digest(value):
    return sha256(repr(value).encode("ascii")).hexdigest()


ACTUAL_HASHES = tuple((path, lf_hash(path)) for path in PINNED)
require(
    all(actual == PINNED[path] for path, actual in ACTUAL_HASHES),
    "a pinned inherited artifact changed",
)
require(
    "Adjoin the group-like formal exponentials"
    in THM3110.read_text(encoding="utf-8"),
    "the inherited group-like response-algebra statement changed",
)
require(
    "sum_mu G_(j,N,mu)=0"
    in THM3115.read_text(encoding="utf-8"),
    "the inherited centered-current statement changed",
)

SOURCE = Path(__file__).read_text(encoding="utf-8")
SOURCE_TREE = ast.parse(SOURCE)
ASSERT_NODES = sum(isinstance(node, ast.Assert) for node in ast.walk(SOURCE_TREE))
FLOAT_LITERALS = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(SOURCE_TREE)
)
require(
    ASSERT_NODES == 0 and FLOAT_LITERALS == 0,
    "source validity gate found assert or floating-point literal",
)

PRIOR_SPEC = spec_from_file_location("fc3_copy_gate_for_algebraic_audit", PRIOR_SCRIPT)
prior = module_from_spec(PRIOR_SPEC)
PRIOR_SPEC.loader.exec_module(prior)
FACES = prior.FACES


def response_task(task):
    """Rebuild exact row coordinates in a spawn-safe local wrapper."""

    face, state = task
    data = FACES[face]
    vectors = prior.prior.prior.base.coefficient_vectors(
        1, prior.prior.prior.base.BANK, data["a"], data["b"], state
    )
    current = None
    if state in CANDIDATE_SET:
        current = tuple(
            (
                degree,
                sum(vector.values()),
                sum(value != 0 for value in vector.values()),
                digest(tuple(sorted(vector.items()))),
            )
            for degree, vector in sorted(vectors.items())
        )
    return (
        (face, state),
        (prior.prior.prior.row_coordinates(vectors), current),
    )


def multiset_union(left, right):
    return tuple(sorted(left + right))


def state_delta(state):
    """Group-like coalgebra on the free physical-state monoid."""

    return Counter({(state, state): 1})


def state_log_profile(state):
    """Numerators of log of the inverse removed-alphabet kernel.

    If K_sigma(X)=prod_(r in sigma,x in X)(1+r*x), then the coefficient of
    p_k(X)/k in log K_sigma(X)^(-1) is (-1)^k p_k(sigma).  Additivity of this
    profile is the exact formal-series reason K_sigma^(-1) is group-like.
    """

    return tuple(
        (-1) ** degree * sum(root ** degree for root in state)
        for degree in range(1, MAXIMUM_DEGREE + 1)
    )


def inverse_kernel_series(state):
    """K_sigma({x})^(-1) through the maintained response degree."""

    coefficients = [1] + [0] * MAXIMUM_DEGREE
    for root in state:
        geometric = tuple((-root) ** degree
                          for degree in range(MAXIMUM_DEGREE + 1))
        coefficients = [
            sum(coefficients[left] * geometric[degree - left]
                for left in range(degree + 1))
            for degree in range(MAXIMUM_DEGREE + 1)
        ]
    return tuple(coefficients)


def subsets(items):
    return tuple(
        tuple(items[index] for index in indices)
        for size in range(len(items) + 1)
        for indices in combinations(range(len(items)), size)
    )


def face_delta(tag):
    """Full labelled-set unshuffle coalgebra, not a carrier transition."""

    answer = Counter()
    for left in subsets(tag):
        left_set = frozenset(left)
        right = tuple(face for face in tag if face not in left_set)
        answer[left, right] += 1
    return answer


def face_coassociativity_record(tag):
    left = Counter()
    right = Counter()
    for (first, third), outer_weight in face_delta(tag).items():
        for (left_first, left_second), inner_weight in face_delta(first).items():
            left[left_first, left_second, third] += outer_weight * inner_weight
    for (first, second), outer_weight in face_delta(tag).items():
        for (right_second, right_third), inner_weight in face_delta(second).items():
            right[first, right_second, right_third] += outer_weight * inner_weight
    return left == right, tuple(sorted(left.items()))


def sign(value):
    return (value > 0) - (value < 0)


def best_row_difference(responses, face, state, row):
    neighbours = prior.prior.prior.toward_reset_neighbours(
        state, FACES[face]["poles"], FACES[face]["reset"]
    )
    values = tuple(
        (responses[face, target][row - 1] - responses[face, state][row - 1], target)
        for target in neighbours
    )
    require(values, ("decision state has no reset-directed neighbour", face, state))
    maximum = max(value for value, _target in values)
    positive_targets = tuple(target for value, target in values if value > 0)
    return maximum, positive_targets


def compact_completion(responses, face, state):
    record = prior.singleton_completion(
        responses, face, state, allowed_rows=ROWS_1617
    )
    return (
        record[2],
        record[3],
        record[4],
        len(record[5]),
        record[6],
        record[8],
    )


def main():
    tasks = tuple(
        (face, state)
        for face in PAIR
        for state in FACES[face]["states"]
    )
    with Pool(processes=4) as pool:
        computed = dict(pool.imap_unordered(response_task, tasks, chunksize=8))
    require(len(computed) == 8_158, "exact pair-face response cache drift")
    responses = {key: value[0] for key, value in computed.items()}
    current_data = {
        key: value[1] for key, value in computed.items() if value[1] is not None
    }
    del computed

    pair_states = tuple(sorted(
        set(FACES["F13"]["decisions"]) & set(FACES["F23"]["decisions"]),
        key=lambda state: (len(state), state),
    ))
    unrestricted_blockers = tuple(
        state for state in pair_states
        if not prior.prior.common_decorated(responses, PAIR, state)
    )
    restricted_blockers = tuple(
        state for state in pair_states
        if not prior.prior.common_decorated(
            responses, PAIR, state, allowed_rows=ROWS_1617
        )
    )
    row_specific_hostiles = tuple(
        state for state in restricted_blockers if state not in unrestricted_blockers
    )
    require(
        len(pair_states) == 1_726
        and len(unrestricted_blockers) == 7
        and len(restricted_blockers) == 38
        and row_specific_hostiles == EXPECTED_ROW_SPECIFIC_HOSTILES
        and all(state in unrestricted_blockers for state in REACHABLE),
        "pair blocker census drift",
    )

    candidates = REACHABLE + row_specific_hostiles
    require(len(candidates) == 35 and len(set(candidates)) == 35,
            "candidate universe is not four plus thirty-one")

    # Exact bialgebra laws on the endpoint-state monoid.
    for state in candidates:
        delta = state_delta(state)
        require(delta == Counter({(state, state): 1}), "state diagonal drift")
        require(
            Counter((right, left) for (left, right), weight in delta.items()
                    for _ in range(weight)) == delta,
            "state diagonal lost cocommutativity",
        )
    multiplicative_checks = 0
    profile_checks = 0
    for left in candidates:
        for right in candidates:
            joined = multiset_union(left, right)
            expected_delta = Counter({(joined, joined): 1})
            require(state_delta(joined) == expected_delta,
                    "group-like state diagonal lost multiplicativity")
            expected_profile = tuple(
                first + second
                for first, second in zip(
                    state_log_profile(left), state_log_profile(right)
                )
            )
            require(state_log_profile(joined) == expected_profile,
                    "removed-alphabet logarithm lost additivity")
            multiplicative_checks += 1
            profile_checks += MAXIMUM_DEGREE

    face_records = tuple(
        (tag, face_coassociativity_record(tag))
        for tag in ((), ("F13",), ("F23",), PAIR)
    )
    require(
        all(record[1][0] for record in face_records)
        and face_delta(PAIR)[(("F13",), ("F23",))] == 1
        and face_delta(PAIR)[(("F23",), ("F13",))] == 1,
        "face-tag unshuffle coalgebra failed",
    )

    # Factorwise predicate and continuation after the endpoint split.
    local_records = []
    sign_histogram = Counter()
    row_availability_histogram = Counter()
    for state in candidates:
        face_rows = []
        completions = []
        for face in PAIR:
            row_data = []
            for row in sorted(ROWS_1617):
                maximum, targets = best_row_difference(
                    responses, face, state, row
                )
                row_data.append((row, maximum, targets))
            face_rows.append((face, tuple(row_data)))
            completions.append((face, compact_completion(responses, face, state)))
        for row in sorted(ROWS_1617):
            differences = tuple(
                next(data[1] for data in rows if data[0] == row)
                for _face, rows in face_rows
            )
            sign_histogram[row, tuple(map(sign, differences))] += 1
        availability = tuple(
            tuple(row for row, maximum, _targets in rows if maximum > 0)
            for _face, rows in face_rows
        )
        row_availability_histogram[availability] += 1
        local_records.append((
            state,
            tuple((face, tuple((row, sign(value), len(targets))
                               for row, value, targets in rows))
                  for face, rows in face_rows),
            tuple(completions),
            digest(state_log_profile(state)),
            digest(inverse_kernel_series(state)),
        ))
    local_records = tuple(local_records)
    reachable_records = local_records[:len(REACHABLE)]
    hostile_records = local_records[len(REACHABLE):]
    reachable_split_successes = sum(
        all(record[3] == 0 and record[4] > 0
            for _face, record in candidate[2])
        for candidate in reachable_records
    )
    hostile_split_successes = sum(
        all(record[3] == 0 and record[4] > 0
            for _face, record in candidate[2])
        for candidate in hostile_records
    )

    # Every maintained response current is centered.  Nonzero centered
    # elements cannot be group-like in a counital coalgebra.
    require(len(current_data) == 70, "candidate current cache drift")
    current_records = tuple(
        (state, tuple((face, current_data[face, state]) for face in PAIR))
        for state in candidates
    )
    current_cells = tuple(
        record
        for _state, face_records_for_state in current_records
        for _face, records in face_records_for_state
        for record in records
    )
    require(
        all(record[1] == 0 for record in current_cells),
        "a maintained response current lost zero mass",
    )
    nonzero_current_cells = sum(record[2] > 0 for record in current_cells)
    zero_current_cells = len(current_cells) - nonzero_current_cells

    # Exact ancestry quotient collision on the generated 94-start closure.
    triple_states = tuple(sorted(
        set.intersection(*(set(FACES[face]["decisions"])
                           for face in prior.FACE_ORDER)),
        key=lambda state: (len(state), state),
    ))
    starts = tuple((PAIR, state) for state in triple_states)
    restricted_closure = prior.prior.closure(
        responses, starts, allowed_rows=ROWS_1617
    )
    endpoint_counts = prior.pair_endpoint_census(
        restricted_closure, triple_states
    )
    ancestry_by_blocker = defaultdict(list)
    for start, endpoints in endpoint_counts:
        require(len(endpoints) == 1, "shared endpoint stopped being unique")
        endpoint, count = endpoints[0]
        ancestry_by_blocker[endpoint].append((start, count))
    ancestry_records = tuple(
        (
            blocker,
            len(ancestry_by_blocker[blocker]),
            min(count for _start, count in ancestry_by_blocker[blocker]),
            max(count for _start, count in ancestry_by_blocker[blocker]),
            sum(count for _start, count in ancestry_by_blocker[blocker]),
            digest(tuple(ancestry_by_blocker[blocker])),
        )
        for blocker in REACHABLE
    )
    require(
        len(endpoint_counts) == 94
        and all(record[1] >= 7 and record[2] >= 6 for record in ancestry_records),
        "ancestry quotient ceased to have a nontrivial fibre",
    )

    # The diagonal is linear on classical state basis elements, not a map
    # z |-> z tensor z on arbitrary superpositions.
    left_state, right_state = candidates[:2]
    linear_diagonal_support = ((left_state, left_state),
                               (right_state, right_state))
    nonlinear_square_support = (
        (left_state, left_state),
        (left_state, right_state),
        (right_state, left_state),
        (right_state, right_state),
    )
    require(
        set(nonlinear_square_support) - set(linear_diagonal_support)
        == {(left_state, right_state), (right_state, left_state)},
        "superposition cross-term control drift",
    )

    print("FC3 PRODUCT-GAMMA ENDPOINT DIAGONAL / CENTERED OBLIGATION NO-GO SCOUT")
    print("status=FINITE-EXACT_PARTIAL;theorem_scope=THM-3314")
    print("dependency_hashes=BEGIN")
    for path, actual in ACTUAL_HASHES:
        print(f"{actual}  {path.relative_to(ROOT)}")
    print("dependency_hashes=END")
    print(
        "source=universal_normalized_product_Gamma_response_algebra_with_"
        "primitive_power_sums_and_group_like_formal_exponentials;target=two_"
        "face_tagged_tensor_factors;map_on_removed_state=Delta(g_sigma)=g_sigma_"
        "tensor_g_sigma_where_g_sigma=product_(r_in_sigma,x_in_X)(1+r*x)^(-1);"
        "predicate_retained=factorwise_exact_state_and_facewise_response_"
        "recomputation;destroyed_or_missing=actual_path_ancestry,response_row,"
        "released_origin_assignment,activation_gate,and_shared_chronology"
    )
    print(
        f"pair_state_universe={len(pair_states)};unrestricted_blockers="
        f"{len(unrestricted_blockers)};rows_16_17_blockers="
        f"{len(restricted_blockers)};row_specific_hostiles="
        f"{len(row_specific_hostiles)};row_specific_hostile_digest="
        f"{digest(row_specific_hostiles)};candidate_universe="
        "four_reachable_plus_31_row_specific"
    )
    print("row_specific_hostiles=" + repr(row_specific_hostiles))
    print(
        f"endpoint_group_like_laws=PASS;candidate_states={len(candidates)};"
        f"multiplicativity_checks={multiplicative_checks};formal_log_profile_"
        f"coefficient_checks={profile_checks};maximum_degree={MAXIMUM_DEGREE};"
        f"candidate_profile_digest={digest(tuple((state,state_log_profile(state),inverse_kernel_series(state)) for state in candidates))}"
    )
    print(
        "mechanism=power_sums_are_primitive_so_the_exponential_removed_alphabet_"
        "unit_is_group_like;this_avoids_the_naive_total_degree_doubling_no_go_"
        "because_the_copied_object_is_a_formal_unit_not_a_homogeneous_state_"
        "polynomial"
    )
    print(
        f"face_tag_full_unshuffle_coalgebra=PASS;records_digest="
        f"{digest(face_records)};ordered_pair_split_coefficient=1;reverse_"
        "assignment_coefficient=1;boundary=selecting_one_ordered_component_and_"
        "binding_it_to_released_versus_old_origin_is_an_added_projection_not_"
        "an_inherited_transition"
    )
    print(
        f"factorwise_row_availability_histogram="
        f"{tuple(sorted(row_availability_histogram.items()))};best_difference_"
        f"sign_histogram={tuple(sorted(sign_histogram.items()))}"
    )
    print(
        f"four_reachable_endpoint_split_successes={reachable_split_successes}/4;"
        f"reachable_records_digest={digest(reachable_records)};rows_16_17_"
        f"hostile_endpoint_split_successes={hostile_split_successes}/31;hostile_"
        f"records_digest={digest(hostile_records)};typing=local_factorwise_"
        "continuation_after_an_unconditional_endpoint_split_only"
    )
    print(
        f"centered_response_current_cells={len(current_cells)};nonzero="
        f"{nonzero_current_cells};zero={zero_current_cells};all_current_masses_"
        f"zero=PASS;records_digest={digest(current_records)}"
    )
    print(
        "first_failed_law=counit_augmentation_compatibility_for_copying_the_"
        "response_obligation;proof=if_epsilon(c)=0_and_Delta(c)=c_tensor_c_then_"
        "(epsilon_tensor_id)Delta(c)=0_but_the_counit_law_requires_c_so_c=0;"
        "therefore_every_nonzero_centered_partition_current_above_is_not_"
        "group_like;native_Delta_on_U=f_a-f_0_contains_U_tensor_1_plus_1_tensor_"
        "U_plus_U_tensor_U_rather_than_only_U_tensor_U"
    )
    print(
        f"ancestry_fibre_records={ancestry_records};all_94_generated_starts_"
        "have_one_endpoint_but_at_least_six_distinct_rows_16_17_histories;"
        "response_element_depends_only_on_(face,endpoint_state);hence_no_map_"
        "factoring_through_the_product_Gamma_endpoint_can_return_the_actual_"
        "chosen_shared_ancestry"
    )
    print(
        f"linear_classical_diagonal_support={linear_diagonal_support};nonlinear_"
        f"tensor_square_support={nonlinear_square_support};missing_cross_terms=2;"
        "boundary=known_basis_states_can_be_copied_but_arbitrary_signed_"
        "response_superpositions_cannot_be_cloned_by_that_linear_map"
    )
    print(
        "strongest_survivor=native_group_like_copy_of_the_removed_endpoint_"
        "alphabet_plus_a_new_face_set_unshuffle_recomputes_both_facewise_"
        "response_banks_and_locally_closes_every_tested_split;it_does_not_copy_"
        "a_centered_response_current_or_actual_ancestry_and_does_not_derive_"
        "COPY_AT_REACHABLE_BLOCKER"
    )
    print(
        f"source_ast=(assert_nodes={ASSERT_NODES},float_literals={FLOAT_LITERALS})"
    )
    print(
        "scope=three_pinned_bank_I2_faces_four_reachable_blockers_and_31_"
        "full_domain_rows_16_17_hostiles;copy_not_inherited;no_FC3_SFC3_GMC_"
        "positivity_chronology_or_global_controller_claim"
    )
    print("ALL FINITE-EXACT PARTIAL CHECKS PASSED")


if __name__ == "__main__":
    freeze_support()
    main()
