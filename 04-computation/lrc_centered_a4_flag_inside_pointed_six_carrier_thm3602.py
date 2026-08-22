#!/usr/bin/env python3
"""Exact THM-3602 comparison of the THM-3593 flag with pointed P6.

Reconstruct the THM-3593 relation flag from its pinned source tensor and
independently reconstruct the six pointed relation rows from the clean-room
four-way audit.  The only comparison map is the typed isomorphism

    im(P_r) -> Fun_0(F_13,F_p),
    h(state,relation) |-> (1/4) sum_state h(state,relation),

followed by relation centering C_13 on the raw pointed carrier P6.  The script
certifies the centered flag, every raw/centered intersection used by THM-3602,
the rank-one augmentation graph recovering P6 from Q6=C_13(P6), and all 156
affine target reindexings t |-> a*t+b.

This is a static exact computation over one pinned finite field.  It supplies
no pointed-difference-to-state lift, chronology, physical current, row
exclusion, characteristic-zero statement, or LRC(14) conclusion.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
import importlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMPUTATION = ROOT / "04-computation"
THM3593_PATH = COMPUTATION / "lrc_common_a4_anova_graph_flag_thm3593.py"
P6_CLEAN_PATH = COMPUTATION / (
    "lrc_r5_owner_node_inverse_branch_root_difference_four_way_"
    "independent_audit_20260816.py"
)

EXPECTED_PARENT_SHA256 = {
    "thm3593": "3ee8486833d539599a4c8add304172b72929c69f3859e9cc1ce22e3018199516",
    "p6_clean": "cf2efd6317dd6ba2a05865ad11d115f1ab239e8ebafb2b8b3764447d1e0e67a3",
}
EXPECTED_PRIME = 755373809845391722745761
EXPECTED_P = 13
EXPECTED_V = 4
EXPECTED_ZETA13 = 266737884585332483769981
EXPECTED_POINTED_STATES = ((0, 0), (1, 0), (1, 6), (3, 6), (3, 12), (2, 12))

EXPECTED_SPACE_DIGESTS = {
    "K2": "193c82c6337f88c9b6c1bb2464198808336f42cc7036f5c3ae30af449e147357",
    "R4": "28bca0c67b0d94431fe3c43181a46d7b55e6cfce98fdc232f837acfde0ccec8c",
    "P6": "6e9083f15408f6d2d85fb3f2747ba0bd1f987e83ce4b836cb7298aaccc84e0c4",
    "Q6": "ae86ad2a3fd03bea95c823b2454b78f244581aa048cb2da63a03a75f484cc596",
    "R4_inter_P6": "ad92a4ba1888bb52d2e9cbcb19fd5967fe10f048c16705817e79201b7b2359bf",
    "R4_inter_Q6": "28bca0c67b0d94431fe3c43181a46d7b55e6cfce98fdc232f837acfde0ccec8c",
    "K2_inter_P6": "193c82c6337f88c9b6c1bb2464198808336f42cc7036f5c3ae30af449e147357",
    "K2_inter_Q6": "193c82c6337f88c9b6c1bb2464198808336f42cc7036f5c3ae30af449e147357",
    "P6_inter_Q6": "739439d079f8639399af6f3737345c6c175af24b51c5c159b46045e52c3420e1",
}
EXPECTED_AFFINE_SHA256 = (
    "f074cff98c370dbebfefc998f96ef761641ba28657a21c1b5c20f95b790678a6"
)
EXPECTED_SEMANTIC_SHA256 = (
    "fdce64910b2eb0c200ef9761cd6518f65c6b3295083c161d00c76a5c8c83d095"
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def import_pinned(path: Path, expected_sha256: str):
    observed = lf_sha256(path)
    require(observed == expected_sha256, (path.name, observed, expected_sha256))
    if str(path.parent) not in sys.path:
        sys.path.insert(0, str(path.parent))
    module = importlib.import_module(path.stem)
    require(Path(module.__file__).resolve() == path.resolve(), (path, module.__file__))
    return module


A = import_pinned(THM3593_PATH, EXPECTED_PARENT_SHA256["thm3593"])
F = import_pinned(P6_CLEAN_PATH, EXPECTED_PARENT_SHA256["p6_clean"])
P, V, PRIME = A.P, A.V, A.PRIME
require((P, V, PRIME) == (EXPECTED_P, EXPECTED_V, EXPECTED_PRIME),
        ("THM-3593 ambient", P, V, PRIME))
require((F.P, F.V, F.PRIME) == (P, V, PRIME),
        ("P6 ambient", F.P, F.V, F.PRIME))


def freeze_rows(rows, width: int | None = None):
    frozen = tuple(tuple(value % PRIME for value in row) for row in rows)
    if width is None and frozen:
        width = len(frozen[0])
    if width is not None:
        require(all(len(row) == width for row in frozen), ("row width", width))
    return frozen


def canonical_row_basis(rows, width: int | None = None):
    matrix = [list(row) for row in freeze_rows(rows, width)]
    if not matrix:
        return ()
    columns = len(matrix[0])
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(pivot_row, len(matrix)) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        inverse = pow(matrix[pivot_row][column], -1, PRIME)
        matrix[pivot_row] = [value * inverse % PRIME for value in matrix[pivot_row]]
        for row in range(len(matrix)):
            if row == pivot_row or matrix[row][column] == 0:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                (left - factor * right) % PRIME
                for left, right in zip(matrix[row], matrix[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return tuple(tuple(row) for row in matrix[:pivot_row])


def rank(rows, width: int | None = None) -> int:
    return len(canonical_row_basis(rows, width))


def rowspace_digest(rows, width: int | None = None) -> str:
    return digest_json(canonical_row_basis(rows, width))


def nullspace(matrix, columns: int):
    rows = canonical_row_basis(matrix, columns)
    pivots = []
    for row in rows:
        pivot = next((column for column, value in enumerate(row) if value), None)
        require(pivot is not None, "zero row in RREF")
        pivots.append(pivot)
    free = tuple(column for column in range(columns) if column not in pivots)
    answer = []
    for free_column in free:
        vector = [0] * columns
        vector[free_column] = 1
        for row, pivot in zip(rows, pivots):
            vector[pivot] = -row[free_column] % PRIME
        answer.append(tuple(vector))
    return tuple(answer)


def rowspace_intersection(left, right, width: int):
    left_basis = canonical_row_basis(left, width)
    right_basis = canonical_row_basis(right, width)
    if not left_basis or not right_basis:
        return ()
    equations = tuple(
        tuple(row[column] for row in left_basis)
        + tuple(-row[column] % PRIME for row in right_basis)
        for column in range(width)
    )
    kernel = nullspace(equations, len(left_basis) + len(right_basis))
    candidates = tuple(
        tuple(
            sum(coefficients[index] * left_basis[index][column]
                for index in range(len(left_basis))) % PRIME
            for column in range(width)
        )
        for coefficients in kernel
    )
    answer = canonical_row_basis(candidates, width)
    require(rank(answer + left_basis, width) == len(left_basis), "intersection-left")
    require(rank(answer + right_basis, width) == len(right_basis), "intersection-right")
    return answer


def coordinates_in_row_basis(basis, target):
    basis = freeze_rows(basis)
    target = tuple(value % PRIME for value in target)
    require(basis and all(len(row) == len(target) for row in basis), "coordinate width")
    variables = len(basis)
    augmented = [
        [basis[index][column] for index in range(variables)] + [target[column]]
        for column in range(len(target))
    ]
    pivot_row = 0
    pivots = {}
    for variable in range(variables):
        pivot = next(
            (row for row in range(pivot_row, len(augmented))
             if augmented[row][variable]),
            None,
        )
        require(pivot is not None, ("dependent basis", variable))
        augmented[pivot_row], augmented[pivot] = augmented[pivot], augmented[pivot_row]
        inverse = pow(augmented[pivot_row][variable], -1, PRIME)
        augmented[pivot_row] = [value * inverse % PRIME
                                for value in augmented[pivot_row]]
        for row in range(len(augmented)):
            if row == pivot_row or augmented[row][variable] == 0:
                continue
            factor = augmented[row][variable]
            augmented[row] = [
                (left - factor * right) % PRIME
                for left, right in zip(augmented[row], augmented[pivot_row])
            ]
        pivots[variable] = pivot_row
        pivot_row += 1
    require(all(any(row[:variables]) or row[-1] == 0 for row in augmented),
            "target outside row span")
    answer = tuple(augmented[pivots[variable]][-1] for variable in range(variables))
    rebuilt = tuple(
        sum(answer[index] * basis[index][column] for index in range(variables)) % PRIME
        for column in range(len(target))
    )
    require(rebuilt == target, "coordinate reconstruction")
    return answer


def relation_profile(row):
    row = tuple(value % PRIME for value in row)
    require(len(row) == V * P, ("relation-profile width", len(row)))
    blocks = tuple(row[state * P:(state + 1) * P] for state in range(V))
    require(all(block == blocks[0] for block in blocks), "not relation-only")
    require(sum(blocks[0]) % PRIME == 0, "not relation-zero-sum")
    # This equals state averaging on im(P_r), because all four slices agree.
    averaged = tuple(
        sum(blocks[state][relation] for state in range(V)) * pow(V, -1, PRIME)
        % PRIME
        for relation in range(P)
    )
    require(averaged == blocks[0], "state-average transport")
    return averaged


def centre_relation(rows):
    inverse_p = pow(P, -1, PRIME)
    return tuple(
        tuple((value - sum(row) * inverse_p) % PRIME for value in row)
        for row in freeze_rows(rows, P)
    )


def affine_reindex(rows, multiplier: int, shift: int):
    require(1 <= multiplier < P and 0 <= shift < P, "affine parameters")
    return tuple(
        tuple(row[(multiplier * relation + shift) % P] for relation in range(P))
        for row in freeze_rows(rows, P)
    )


def combine_rows(coefficients, rows):
    rows = freeze_rows(rows)
    return tuple(
        sum(coefficients[index] * rows[index][column]
            for index in range(len(rows))) % PRIME
        for column in range(len(rows[0]))
    )


def reconstruct_relation_flag():
    tensor, reconstruction = A.M.reconstruct_source_current()
    require(reconstruction == A.M.EXPECTED_PARENT_RECONSTRUCTION_SHA256[2:],
            ("THM-3593 source reconstruction", reconstruction))
    raw = A.M.flatten_source_current(tensor)
    constant, state, relation, _interaction = A.component_families(raw)
    additive = A.add_rows(constant, state, relation)
    _pure_state, _state_only, pure_relation, _relation_only = A.standard_spaces()

    r4 = canonical_row_basis(tuple(relation_profile(row) for row in relation), P)
    k2_ambient = rowspace_intersection(additive, pure_relation, V * P)
    k2 = canonical_row_basis(tuple(relation_profile(row) for row in k2_ambient), P)
    require((len(k2), len(r4)) == (2, 4), ("relation flag ranks", len(k2), len(r4)))
    require(rank(k2 + r4, P) == 4, "K2 not contained in R4")
    return k2, r4, reconstruction


def reconstruct_pointed_carrier():
    (
        _gamma_actual, _gamma_support, gamma_pointed,
        point_cores, _support_cores, work_counts, state_counts, branch_counts,
    ) = F.build_banks()
    require(work_counts == F.EXPECTED_WORK_COUNTS, ("P6 work counts", work_counts))
    require(state_counts == F.EXPECTED_STATE_COUNTS, ("P6 state counts", state_counts))
    require(branch_counts == F.EXPECTED_BRANCH_COUNTS,
            ("P6 branch counts", branch_counts))
    require(F.digest_json(gamma_pointed) == F.PS.EXPECTED_DIGESTS[0][0],
            "P6 pointed gamma drift")

    zeta = pow(F.R.JOINT_ROOT, F.R.JOINT_ORDER // P, PRIME)
    point_branch = F.invert_point_cores(point_cores, zeta)
    pointed = F.pointed_marginal(point_branch)
    require(F.digest_json(pointed) == F.PS.EXPECTED_DIGESTS[1][0],
            "P6 pointed tensor drift")
    p6_rows = tuple(
        tuple(
            sum(pointed[point][difference][relation] for difference in range(P))
            % PRIME
            for relation in range(P)
        )
        for point in range(len(F.POINTS))
    )
    p6 = canonical_row_basis(p6_rows, P)
    q6 = canonical_row_basis(centre_relation(p6_rows), P)
    return p6, q6, zeta, tuple(F.POINTS), F.digest_json(pointed)


def main() -> None:
    k2, r4, reconstruction = reconstruct_relation_flag()
    p6, q6, zeta_p6, pointed_states, pointed_tensor_digest = (
        reconstruct_pointed_carrier()
    )
    zeta_a4 = A.M.SOURCE.B.context()["zeta"]
    require(zeta_a4 == zeta_p6 == EXPECTED_ZETA13,
            ("zeta alignment", zeta_a4, zeta_p6))
    require(pointed_states == EXPECTED_POINTED_STATES,
            ("pointed state order", pointed_states))

    zero_sum = canonical_row_basis(
        tuple(
            tuple((int(relation == target) - int(relation == P - 1)) % PRIME
                  for relation in range(P))
            for target in range(P - 1)
        ),
        P,
    )
    ones = ((1,) * P,)
    require(len(zero_sum) == 12, "zero-sum ambient rank")
    require(all(sum(row) % PRIME == 0 for row in r4 + k2 + q6),
            "centered space misses Fun_0")

    spaces = {"K2": k2, "R4": r4, "P6": p6, "Q6": q6}
    space_ranks = {name: len(rows) for name, rows in spaces.items()}
    space_digests = {name: rowspace_digest(rows, P) for name, rows in spaces.items()}
    require(space_ranks == {"K2": 2, "R4": 4, "P6": 6, "Q6": 6},
            ("space ranks", space_ranks))
    require(space_digests == {name: EXPECTED_SPACE_DIGESTS[name] for name in spaces},
            ("space digests", space_digests))

    stack_ranks = {
        "K2_R4": rank(k2 + r4, P),
        "R4_P6": rank(r4 + p6, P),
        "R4_Q6": rank(r4 + q6, P),
        "K2_P6": rank(k2 + p6, P),
        "K2_Q6": rank(k2 + q6, P),
        "P6_Q6": rank(p6 + q6, P),
        "P6_constants": rank(p6 + ones, P),
        "Q6_constants": rank(q6 + ones, P),
        "Q6_zero_sum": rank(q6 + zero_sum, P),
    }
    require(stack_ranks == {
        "K2_R4": 4, "R4_P6": 7, "R4_Q6": 6,
        "K2_P6": 6, "K2_Q6": 6, "P6_Q6": 7,
        "P6_constants": 7, "Q6_constants": 7, "Q6_zero_sum": 12,
    }, ("stack ranks", stack_ranks))

    intersection_pairs = {
        "R4_inter_P6": (r4, p6),
        "R4_inter_Q6": (r4, q6),
        "K2_inter_P6": (k2, p6),
        "K2_inter_Q6": (k2, q6),
        "P6_inter_Q6": (p6, q6),
    }
    intersections = {
        name: rowspace_intersection(left, right, P)
        for name, (left, right) in intersection_pairs.items()
    }
    intersection_dimensions = {name: len(rows) for name, rows in intersections.items()}
    intersection_digests = {
        name: rowspace_digest(rows, P) for name, rows in intersections.items()
    }
    require(intersection_dimensions == {
        "R4_inter_P6": 3, "R4_inter_Q6": 4,
        "K2_inter_P6": 2, "K2_inter_Q6": 2, "P6_inter_Q6": 5,
    }, ("intersection dimensions", intersection_dimensions))
    require(intersection_digests == {
        name: EXPECTED_SPACE_DIGESTS[name] for name in intersection_pairs
    }, ("intersection digests", intersection_digests))

    # Since centering is injective on P6, use a raw P6 basis and its centered
    # image to exhibit P6 as the graph q |-> q + lambda(q)*1 over Q6.
    q_images = centre_relation(p6)
    require(rank(q_images, P) == 6, "centering is not injective on P6")
    lambda_values = []
    for raw_row, centered_row in zip(p6, q_images):
        difference = tuple((left - right) % PRIME
                           for left, right in zip(raw_row, centered_row))
        require(all(value == difference[0] for value in difference),
                "augmentation difference is not constant")
        lambda_values.append(difference[0])
    lambda_values = tuple(lambda_values)
    require(any(lambda_values), "augmentation functional has rank zero")
    require(all(
        tuple((q_images[index][column] + lambda_values[index]) % PRIME
              for column in range(P)) == p6[index]
        for index in range(len(p6))
    ), "augmentation graph reconstruction")

    lambda_on_r4 = tuple(
        sum(coefficient * value for coefficient, value in
            zip(coordinates_in_row_basis(q_images, row), lambda_values)) % PRIME
        for row in r4
    )
    lambda_on_k2 = tuple(
        sum(coefficient * value for coefficient, value in
            zip(coordinates_in_row_basis(q_images, row), lambda_values)) % PRIME
        for row in k2
    )
    require(any(lambda_on_r4), "lambda restricted to R4 has rank zero")
    require(not any(lambda_on_k2), "lambda does not vanish on K2")

    global_kernel_coefficients = nullspace((lambda_values,), len(p6))
    global_kernel = canonical_row_basis(
        tuple(combine_rows(coefficients, q_images)
              for coefficients in global_kernel_coefficients),
        P,
    )
    r4_kernel_coefficients = nullspace((lambda_on_r4,), len(r4))
    r4_kernel = canonical_row_basis(
        tuple(combine_rows(coefficients, r4)
              for coefficients in r4_kernel_coefficients),
        P,
    )
    require(global_kernel == intersections["P6_inter_Q6"],
            "global lambda kernel mismatch")
    require(r4_kernel == intersections["R4_inter_P6"],
            "R4 lambda kernel mismatch")

    augmentation = {
        "centering_domain_rank": 6,
        "centering_image_rank": rank(q_images, P),
        "centering_kernel_dimension": 6 - rank(q_images, P),
        "lambda_rank": int(any(lambda_values)),
        "lambda_kernel_dimension": len(global_kernel),
        "lambda_on_R4_rank": int(any(lambda_on_r4)),
        "lambda_on_R4_kernel_dimension": len(r4_kernel),
        "lambda_on_K2_rank": int(any(lambda_on_k2)),
        "lambda_basis_values_sha256": digest_json(lambda_values),
        "lambda_R4_values_sha256": digest_json(lambda_on_r4),
    }
    require(tuple(augmentation[key] for key in (
        "centering_domain_rank", "centering_image_rank",
        "centering_kernel_dimension", "lambda_rank", "lambda_kernel_dimension",
        "lambda_on_R4_rank", "lambda_on_R4_kernel_dimension",
        "lambda_on_K2_rank",
    )) == (6, 6, 0, 1, 5, 1, 3, 0), ("augmentation ranks", augmentation))

    projection_ranks = {
        "R4_to_F13_mod_P6": rank(r4 + p6, P) - len(p6),
        "R4_to_Fun0_mod_Q6": rank(r4 + q6, P) - len(q6),
        "K2_to_Fun0_mod_Q6": rank(k2 + q6, P) - len(q6),
    }
    require(projection_ranks == {
        "R4_to_F13_mod_P6": 1,
        "R4_to_Fun0_mod_Q6": 0,
        "K2_to_Fun0_mod_Q6": 0,
    }, ("projection ranks", projection_ranks))

    affine_atlas = []
    for multiplier in range(1, P):
        for shift in range(P):
            image = affine_reindex(q6, multiplier, shift)
            r_intersection = len(rowspace_intersection(r4, image, P))
            k_intersection = len(rowspace_intersection(k2, image, P))
            affine_atlas.append((multiplier, shift, r_intersection, k_intersection))
    affine_atlas = tuple(affine_atlas)
    affine_histogram = tuple(sorted(Counter(
        (r_dimension, k_dimension)
        for _multiplier, _shift, r_dimension, k_dimension in affine_atlas
    ).items()))
    affine_exceptions = tuple(
        record for record in affine_atlas if record[2:] != (0, 0)
    )
    affine_digest = digest_json(affine_atlas)
    require(len(affine_atlas) == 156, ("affine atlas size", len(affine_atlas)))
    require(affine_histogram == (((0, 0), 155), ((4, 2), 1)),
            ("affine histogram", affine_histogram))
    require(affine_exceptions == ((1, 0, 4, 2),),
            ("affine exceptions", affine_exceptions))
    require(affine_digest == EXPECTED_AFFINE_SHA256,
            ("affine digest", affine_digest, EXPECTED_AFFINE_SHA256))

    semantic_record = (
        PRIME, V, P, zeta_p6, pointed_states,
        tuple(sorted(EXPECTED_PARENT_SHA256.items())), reconstruction,
        pointed_tensor_digest, tuple(sorted(space_ranks.items())),
        tuple(sorted(space_digests.items())), tuple(sorted(stack_ranks.items())),
        tuple(sorted(intersection_dimensions.items())),
        tuple(sorted(intersection_digests.items())),
        tuple(sorted(augmentation.items())), tuple(sorted(projection_ranks.items())),
        affine_histogram, affine_exceptions, affine_digest,
        "static finite-field relation target; difference coordinate marginalized",
    )
    semantic = digest_json(semantic_record)
    require(semantic == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== THM-3602 centered A4 flag inside pointed six carrier ==")
    print(f"field=(prime={PRIME},state={V},relation={P},H_dimension={V*P},Fun0_dimension=12)")
    print(f"parent_sha256_lf={EXPECTED_PARENT_SHA256}")
    print(f"shared_relation_order={tuple(range(P))};shared_zeta13={zeta_p6}")
    print(f"pointed_states={pointed_states}")
    print(f"space_ranks={space_ranks}")
    print(f"space_digests={space_digests}")
    print(f"stack_ranks={stack_ranks}")
    print(f"intersection_dimensions={intersection_dimensions}")
    print(f"intersection_digests={intersection_digests}")
    print(f"augmentation_graph={augmentation}")
    print(f"projection_ranks={projection_ranks}")
    print(f"flag=K2(2)<R4(4)<Q6(6)<Fun0(F13)(12)")
    print(f"refined_flag=K2(2)<(R4_inter_P6)(3)<R4(4)<Q6(6)")
    print(f"affine_reindexing_count={len(affine_atlas)}")
    print(f"affine_intersection_histogram={affine_histogram}")
    print(f"affine_nonzero_exceptions={affine_exceptions}")
    print(f"affine_atlas_sha256={affine_digest}")
    print(f"semantic_sha256={semantic}")
    print(f"script_sha256_lf={lf_sha256(Path(__file__).resolve())}")
    print("status=FINITE-EXACT companion candidate;independent audit still required")
    print("scope=static finite-field relation target;not difference-fibre lift;not chronology/current/entry/characteristic-zero/LRC14")
    print("PASS")


if __name__ == "__main__":
    main()
