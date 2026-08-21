#!/usr/bin/env python3
"""Exact pointed root-difference lift of the THM-3602 relation flag.

The clean-room four-way parent has coordinates

    point x branch_digit x root_difference x relation.

Sum the branch digit but retain root difference.  The resulting six rows span
P6sharp in Fun(F13_difference x F13_relation,Fp).  This companion proves that
difference marginalization identifies P6sharp with THM-3602's P6 and, after
relation centering, identifies Q6sharp with Q6.  It then lifts K2<R4<Q6,
certifies the rank-two sharp augmentation and its one-dimensional lost
difference mode, and runs exact slice, Fourier, reversal, affine, and
coordinate-swap hostiles.

The branch-resolved rows and their branch sums have the same target rowspace,
but no coefficient/address intertwiner is constructed.  The theorem is static
and finite-field only: no chronology, current, row exclusion, characteristic-
zero transfer, or LRC(14) conclusion follows.

The rank-two augmentation target A2sharp is abstract.  In particular, the
tempting physical map Delta(c)=(c1-c2,c3-c4) on the two duplicate-state point
pairs neither has kernel K2sharp nor factors LambdaSharp; an exact hostile
record below prevents that false identification from being reintroduced.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
import importlib
import json
from itertools import combinations
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
CANON = ROOT / "01-canon/theorems"
COMPUTATION = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge/results"

THM3602_THEOREM_PATH = (
    CANON / "THM-3602-lrc-centered-a4-flag-inside-pointed-six-carrier.md"
)
THM3602_SCRIPT_PATH = (
    COMPUTATION / "lrc_centered_a4_flag_inside_pointed_six_carrier_thm3602.py"
)
THM3602_OUTPUT_PATH = (
    RESULTS / "lrc_centered_a4_flag_inside_pointed_six_carrier_thm3602.out"
)
CLEAN_PARENT_PATH = COMPUTATION / (
    "lrc_r5_owner_node_inverse_branch_root_difference_four_way_"
    "independent_audit_20260816.py"
)

EXPECTED_PARENT_SHA256 = {
    "thm3602_theorem": "457d6d422a9cf5f7ff4068264b9a23862fa5aea18e558792754c16dc3c5931a6",
    "thm3602_script": "c6954206f6a3187a98772fdf2313f0fd91d7df12f38ec7438d39be7bad7235e8",
    "thm3602_output": "bdb1a7cb2b85c1dce1ad1dfce54c9e02b4c935ff810101768fd08030d6531d56",
    "clean_four_way": "cf2efd6317dd6ba2a05865ad11d115f1ab239e8ebafb2b8b3764447d1e0e67a3",
}

EXPECTED_PRIME = 755373809845391722745761
EXPECTED_P = 13
EXPECTED_V = 4
EXPECTED_ZETA13 = 266737884585332483769981
EXPECTED_POINTS = ((0, 0), (1, 0), (1, 6), (3, 6), (3, 12), (2, 12))
EXPECTED_POINTED_TENSOR_SHA256 = (
    "9c5e227d9c142373973a562a54c6a67cac60a82da1121a028b9658920d155a19"
)

EXPECTED_MARGINAL_RECORDS = {
    "P6sharp_to_P6": (
        6, 6, 0,
        "09c232a61bd3bd308cc4b66e5c22b8ff4c177cade263b6fea12b519dd355b2c2",
    ),
    "Q6sharp_to_Q6": (
        6, 6, 0,
        "0291a20c81074a28474f986e526522397d0e252c2bdb0ad8600462a4aeb59d26",
    ),
}

EXPECTED_SPACE_DIGESTS = {
    "K2sharp": "7b73626f894ea6730fde769fbe57abbc3dca13819c3ab38e3c62aa6fbe2a95af",
    "R4sharp": "f0ca0f513a106e883962ea97f1df6ed56ddf584472a5186d02275c92f3570936",
    "P6sharp": "6014f99bee469da11163238341e815f6901d1460ec2f3d1ac62515b65470fda5",
    "Q6sharp": "167fb13d3fea42a269881e4be9bd7d2e7405bc5a19c50a7d01c0f9c5ec0c0240",
}

EXPECTED_INTERSECTION_DIMENSIONS = {
    "K2sharp_inter_R4sharp": 2,
    "K2sharp_inter_P6sharp": 2,
    "K2sharp_inter_Q6sharp": 2,
    "R4sharp_inter_P6sharp": 2,
    "R4sharp_inter_Q6sharp": 4,
    "P6sharp_inter_Q6sharp": 4,
}

EXPECTED_INTERSECTION_DIGESTS = {
    "K2sharp_inter_R4sharp": EXPECTED_SPACE_DIGESTS["K2sharp"],
    "K2sharp_inter_P6sharp": EXPECTED_SPACE_DIGESTS["K2sharp"],
    "K2sharp_inter_Q6sharp": EXPECTED_SPACE_DIGESTS["K2sharp"],
    "R4sharp_inter_P6sharp": EXPECTED_SPACE_DIGESTS["K2sharp"],
    "R4sharp_inter_Q6sharp": EXPECTED_SPACE_DIGESTS["R4sharp"],
    "P6sharp_inter_Q6sharp": (
        "bd56007794e27feb2800d4a45ecde11b1b0621dc9629ef8f466ffd949bccc6cb"
    ),
}

EXPECTED_AUGMENTATION = {
    "sharp_rank": 2,
    "sharp_sha256": "33697ba84b3d4a3d6f64ea4577bcef00f5a3be511a3de0cbe9ff3b1bdf6914d3",
    "difference_profile_rank": 2,
    "difference_profile_sha256": (
        "b61977c3665a790c561151aacfb17abde5efd013429651b32c80631120985d70"
    ),
    "marginal_rank": 1,
    "marginal_sha256": "fa0409594f7a0aa164fbb64068ab145feede6c5769baaf14964c07ffb8785cb1",
    "marginal_kernel_dimension": 1,
    "marginal_kernel_sha256": (
        "bbf00a18f8ec98a9a07661f5df8c9413d3d493faac29ed5d459710332a9c994e"
    ),
    "lambda_on_R4sharp_rank": 2,
    "lambda_on_R4sharp_image_sha256": (
        "33697ba84b3d4a3d6f64ea4577bcef00f5a3be511a3de0cbe9ff3b1bdf6914d3"
    ),
    "lambda_on_R4sharp_kernel_dimension": 2,
    "lambda_on_R4sharp_kernel_sha256": EXPECTED_SPACE_DIGESTS["K2sharp"],
    "lambda_on_K2sharp_rank": 0,
    "marginal_lambda_on_R4sharp_rank": 1,
    "marginal_lambda_on_R4sharp_kernel_dimension": 3,
    "marginal_lambda_on_R4sharp_kernel_sha256": (
        "e04ab1f05892cdf142fd8ff3d16801e5ae5181246521f1852a1351b0c87699d2"
    ),
    "reversal": {
        "sharp_plus_minus_ranks": (1, 1),
        "marginal_plus_minus_ranks": (1, 1),
        "plus_sha256": "979d31a886c99c7bac225a87f0a53135e23f07a4ae494253442b57504cec92af",
        "minus_sha256": "bfd775153cae2f9129491ec006765576dba939855fef06a3e2ae6b230181d24a",
    },
}

EXPECTED_DUPLICATE_STATE_HOSTILE = {
    "point_order": EXPECTED_POINTS,
    "delta_definition": "(c1-c2,c3-c4)",
    "point_correction_pattern": "(0,+H1,-H1,+H2,-H2,0)",
    "point_correction_pattern_flags": (True, False, False, True),
    "point_correction_pattern_holds": False,
    "point_corrections_ordered_sha256": (
        "ca0c79202edd5034d971a1b8f73bb53f0998eac66bbc476f7d09a4226c94b09c"
    ),
    "pattern_residual_rank": 2,
    "pattern_residual_ordered_sha256": (
        "11f6e5f424811fc4d8dc97bc6bccc0b1ee8d3dbcb37e73b9dd392dea8e35b0d2"
    ),
    "pattern_residual_rowspace_sha256": (
        "33697ba84b3d4a3d6f64ea4577bcef00f5a3be511a3de0cbe9ff3b1bdf6914d3"
    ),
    "delta_on_R4sharp_rank": 2,
    "delta_on_R4sharp_sha256": (
        "5c57ba69b9d73197b874254e5b67eb9a987366ced21f1f2d40e3f5a15b864dfb"
    ),
    "delta_graph_sha256": (
        "ac4be69cd44d49955d24cd8a197a7b35f847e97c89a08c163af1ef6007b34b29"
    ),
    "delta_kernel_dimension": 2,
    "delta_kernel_sha256": (
        "52ef24befd72009f398e4ccf182f4ebf66cfe369bd22163d06c591b213f01604"
    ),
    "delta_on_K2sharp_rank": 2,
    "delta_kernel_inter_K2sharp_dimension": 0,
    "delta_kernel_inter_K2sharp_sha256": (
        "4f53cda18c2baa0c0354bb5f9a3ecbe5ed12ab4d8e11ba873c2f11161202b945"
    ),
    "delta_kernel_equals_K2sharp": False,
    "lambda_coordinate_rank": 2,
    "delta_lambda_joint_rank": 4,
    "lambda_factors_through_delta": False,
    "identification_scope": "A2sharp is abstract, not duplicate-state physical",
}

EXPECTED_DIFFERENCE_SLICE_RANKS = {
    "K2sharp": (0,) + (1,) * 12,
    "R4sharp": (0, 3, 3, 3, 3, 3, 2, 2, 3, 3, 3, 3, 3),
    "P6sharp": (0, 3, 3, 3, 3, 3, 2, 2, 3, 3, 3, 3, 3),
    "Q6sharp": (0, 3, 3, 3, 3, 3, 2, 2, 3, 3, 3, 3, 3),
}

EXPECTED_RELATION_SLICE_RANKS = {
    "K2sharp": (2,) * 13,
    "R4sharp": (4,) * 13,
    "P6sharp": (6,) * 13,
    "Q6sharp": (6,) * 13,
}

EXPECTED_FOURIER = {
    "K2sharp": {
        "nonzero_count": 130,
        "difference_slice_ranks": (2,) * 13,
        "difference_frequency_support": tuple(range(13)),
        "support_rows": (10,) * 13,
        "support_columns": (0, 0) + (13,) * 10 + (0,),
        "support_sha256": "a95110de44c6daff17595d42892b6dc380a3ba658cac0ab3bf585b7102fee460",
    },
    "R4sharp": {
        "nonzero_count": 156,
        "difference_slice_ranks": (4,) * 13,
        "difference_frequency_support": tuple(range(13)),
        "support_rows": (12,) * 13,
        "support_columns": (0,) + (13,) * 12,
        "support_sha256": "71faefe862d1fd6473469c7af281e679ab9f7d0b976d8e1a8c4cbebc6688626d",
    },
    "P6sharp": {
        "nonzero_count": 169,
        "difference_slice_ranks": (6,) * 13,
        "difference_frequency_support": tuple(range(13)),
        "support_rows": (13,) * 13,
        "support_columns": (13,) * 13,
        "support_sha256": "236e8613c6ad9fb99b020a4dacd73ffb04b4f34bceb6a1ca9741485ee18d7188",
    },
    "Q6sharp": {
        "nonzero_count": 156,
        "difference_slice_ranks": (6,) * 13,
        "difference_frequency_support": tuple(range(13)),
        "support_rows": (12,) * 13,
        "support_columns": (0,) + (13,) * 12,
        "support_sha256": "71faefe862d1fd6473469c7af281e679ab9f7d0b976d8e1a8c4cbebc6688626d",
    },
}

EXPECTED_REVERSAL = {
    "point_reversal": (5, 4, 3, 2, 1, 0),
    "sharp_graph_rank": 6,
    "sharp_graph_sha256": "039753708e18bdc32f27485bceae10a87ce7ce7667da9e082b8ee274f45f30c9",
    "marginal_graph_rank": 6,
    "marginal_graph_sha256": (
        "67850c2586ac9a4f242ada7173eba2fe1b3c5e6d1ae20a5e2b6d7ad0dc7d3f88"
    ),
    "sharp_plus_minus_ranks": (3, 3),
    "target_difference_negation_relation_affine_stabilizers": (),
    "pointwise_target_covariance_solutions": (),
    "mu_graph_equivariant": True,
}

EXPECTED_AFFINE = {
    "relation_flag_histogram": (((0, 0, 0), 155), ((4, 2, 6), 1)),
    "relation_flag_nonzero": ((1, 0, 4, 2, 6),),
    "relation_flag_sha256": "b90a1366ecfb55134816025974c931f2555a2b6cf538e57c674e6b46bdb593c2",
    "difference_flag_histogram": (((0, 0, 0), 155), ((4, 2, 6), 1)),
    "difference_flag_sha256": "b90a1366ecfb55134816025974c931f2555a2b6cf538e57c674e6b46bdb593c2",
    "relation_self_histogram": (((0, 0), 155), ((6, 6), 1)),
    "relation_self_stabilizers": ((1, 0),),
    "relation_self_sha256": "52e94fe22f0b18a4e3a25ace96c10ee08aa453e00d0442402445e64c4288d9e9",
    "difference_self_histogram": (((0, 0), 155), ((6, 6), 1)),
    "difference_self_stabilizers": ((1, 0),),
    "difference_self_sha256": "52e94fe22f0b18a4e3a25ace96c10ee08aa453e00d0442402445e64c4288d9e9",
    "relation_and_difference_flag_atlases_equal": True,
    "relation_and_difference_self_atlases_equal": True,
    "coordinate_swap": {
        "K2sharp": (
            0, 4, "84931f695994c9c3b92763dba6433e52f0d4f7f5c506a271642dcb6272e8b5a4"
        ),
        "R4sharp": (
            0, 8, "ff6260dae4b6924e39ab20f87016facb9a2dd65d8d3c59ebe521ebc0d5bbb506"
        ),
        "P6sharp": (
            0, 12, "35a55b4461462cd756c22f08a504ed8876b92a59f3a096b3a53716f9fb07a1b5"
        ),
        "Q6sharp": (
            0, 12, "d06b861b472e67fd3445de7173d642de428ff0c2a2f13397607e9b16ca9a09a2"
        ),
    },
}

EXPECTED_SEMANTIC_SHA256 = (
    "7994b4ff0428249ac299e6a90517de7e36d0329603f7473d648e199677b1e9c9"
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


for parent_name, parent_path in (
    ("thm3602_theorem", THM3602_THEOREM_PATH),
    ("thm3602_script", THM3602_SCRIPT_PATH),
    ("thm3602_output", THM3602_OUTPUT_PATH),
    ("clean_four_way", CLEAN_PARENT_PATH),
):
    observed = lf_sha256(parent_path)
    require(observed == EXPECTED_PARENT_SHA256[parent_name],
            (parent_name, observed, EXPECTED_PARENT_SHA256[parent_name]))

require("id: THM-3602" in THM3602_THEOREM_PATH.read_text(encoding="utf-8"),
        "THM-3602 theorem identity")
require(THM3602_OUTPUT_PATH.read_text(encoding="utf-8").rstrip().endswith("PASS"),
        "THM-3602 stored output boundary")

if str(COMPUTATION) not in sys.path:
    sys.path.insert(0, str(COMPUTATION))
C = importlib.import_module(THM3602_SCRIPT_PATH.stem)
F = C.F
require(Path(C.__file__).resolve() == THM3602_SCRIPT_PATH.resolve(), "THM-3602 import")
require(Path(F.__file__).resolve() == CLEAN_PARENT_PATH.resolve(), "clean parent import")

P, V, PRIME = C.P, C.V, C.PRIME
WIDTH = P * P
require((P, V, PRIME) == (EXPECTED_P, EXPECTED_V, EXPECTED_PRIME),
        ("ambient", P, V, PRIME))


def mu_row(row):
    require(len(row) == WIDTH, ("mu width", len(row)))
    return tuple(
        sum(row[difference * P + relation] for difference in range(P)) % PRIME
        for relation in range(P)
    )


def mu_rows(rows):
    return tuple(mu_row(row) for row in rows)


def flatten_pointed(pointed):
    return tuple(
        tuple(
            pointed[point][difference][relation]
            for difference in range(P)
            for relation in range(P)
        )
        for point in range(len(F.POINTS))
    )


def flatten_point_branch(point_branch):
    return tuple(
        tuple(
            point_branch[point][branch][difference][relation]
            for difference in range(P)
            for relation in range(P)
        )
        for point in range(len(F.POINTS))
        for branch in range(P)
    )


def centre_relation_sharp(rows):
    inverse_p = pow(P, -1, PRIME)
    answer = []
    for row in rows:
        require(len(row) == WIDTH, ("sharp width", len(row)))
        answer.append(tuple(
            (
                row[difference * P + relation]
                - sum(row[difference * P:(difference + 1) * P]) * inverse_p
            ) % PRIME
            for difference in range(P)
            for relation in range(P)
        ))
    return tuple(answer)


def lift_basis(target_basis, sharp_basis, marginal_basis):
    require(len(sharp_basis) == len(marginal_basis), "lift basis arity")
    answer = tuple(
        C.combine_rows(C.coordinates_in_row_basis(marginal_basis, target), sharp_basis)
        for target in target_basis
    )
    require(mu_rows(answer) == tuple(target_basis), "lift marginal")
    return C.canonical_row_basis(answer, WIDTH)


def physical_slice_ranks(rows):
    return tuple(
        C.rank(
            tuple(row[difference * P:(difference + 1) * P] for row in rows), P
        )
        for difference in range(P)
    )


def relation_slice_ranks(rows):
    return tuple(
        C.rank(
            tuple(
                tuple(row[difference * P + relation] for difference in range(P))
                for row in rows
            ),
            P,
        )
        for relation in range(P)
    )


def physical_support(rows):
    return tuple(int(any(row[index] for row in rows)) for index in range(WIDTH))


def relation_transform(rows, multiplier: int, shift: int):
    return tuple(
        tuple(
            row[difference * P + (multiplier * relation + shift) % P]
            for difference in range(P)
            for relation in range(P)
        )
        for row in rows
    )


def difference_transform(rows, multiplier: int, shift: int):
    return tuple(
        tuple(
            row[((multiplier * difference + shift) % P) * P + relation]
            for difference in range(P)
            for relation in range(P)
        )
        for row in rows
    )


def transpose_coordinates(rows):
    return tuple(
        tuple(
            row[relation * P + difference]
            for difference in range(P)
            for relation in range(P)
        )
        for row in rows
    )


def scalar_rows(rows, scalar: int):
    return tuple(
        tuple(scalar * value % PRIME for value in row) for row in rows
    )


def intersection_dimension(left, right, width=WIDTH):
    return C.rank(left, width) + C.rank(right, width) - C.rank(left + right, width)


def linear_map_kernel(domain_basis, image_rows, image_width: int):
    require(len(domain_basis) == len(image_rows), "linear map arity")
    coefficient_kernel = C.nullspace(
        tuple(
            tuple(image_rows[index][column] for index in range(len(image_rows)))
            for column in range(image_width)
        ),
        len(image_rows),
    )
    return C.canonical_row_basis(
        tuple(C.combine_rows(coefficients, domain_basis)
              for coefficients in coefficient_kernel),
        len(domain_basis[0]),
    )


def fourier_record(rows, zeta: int):
    roots = tuple(
        tuple(pow(zeta, frequency * coordinate % P, PRIME)
              for coordinate in range(P))
        for frequency in range(P)
    )
    support = []
    transformed_by_difference = []
    for difference_frequency in range(P):
        transformed_rows = tuple(
            tuple(
                sum(
                    row[difference * P + relation]
                    * roots[difference_frequency][difference]
                    for difference in range(P)
                ) % PRIME
                for relation in range(P)
            )
            for row in rows
        )
        transformed_by_difference.append(transformed_rows)
        for relation_frequency in range(P):
            support.append(int(any(
                sum(
                    transformed[relation] * roots[relation_frequency][relation]
                    for relation in range(P)
                ) % PRIME
                for transformed in transformed_rows
            )))
    support = tuple(support)
    difference_ranks = tuple(C.rank(rows_at_frequency, P)
                             for rows_at_frequency in transformed_by_difference)
    return {
        "nonzero_count": sum(support),
        "difference_slice_ranks": difference_ranks,
        "difference_frequency_support": tuple(
            frequency for frequency, value in enumerate(difference_ranks) if value
        ),
        "support_rows": tuple(
            sum(support[frequency * P:(frequency + 1) * P])
            for frequency in range(P)
        ),
        "support_columns": tuple(
            sum(support[difference_frequency * P + relation_frequency]
                for difference_frequency in range(P))
            for relation_frequency in range(P)
        ),
        "support_sha256": digest_json(support),
    }


def affine_atlas(flag_r, flag_k, carrier, coordinate: str):
    answer = []
    for multiplier in range(1, P):
        for shift in range(P):
            image = (
                relation_transform(carrier, multiplier, shift)
                if coordinate == "relation"
                else difference_transform(carrier, multiplier, shift)
            )
            answer.append((
                multiplier,
                shift,
                intersection_dimension(flag_r, image),
                intersection_dimension(flag_k, image),
                intersection_dimension(carrier, image),
            ))
    return tuple(answer)


def self_affine_atlas(raw, centered, coordinate: str):
    answer = []
    for multiplier in range(1, P):
        for shift in range(P):
            transform = relation_transform if coordinate == "relation" else difference_transform
            raw_image = transform(raw, multiplier, shift)
            centered_image = transform(centered, multiplier, shift)
            answer.append((
                multiplier,
                shift,
                intersection_dimension(raw, raw_image),
                intersection_dimension(centered, centered_image),
            ))
    return tuple(answer)


def histogram(atlas, start: int):
    return tuple(sorted(Counter(record[start:] for record in atlas).items()))


def reversal_record(generators, marginal_generators):
    lookup = {point: index for index, point in enumerate(F.POINTS)}
    reversal = tuple(
        lookup[(state ^ 2, P - 1 - root)] for state, root in F.POINTS
    )
    reflected = tuple(generators[index] for index in reversal)
    marginal_reflected = tuple(marginal_generators[index] for index in reversal)
    graph = tuple(tuple(left) + tuple(right)
                  for left, right in zip(generators, reflected))
    marginal_graph = tuple(tuple(left) + tuple(right)
                           for left, right in zip(marginal_generators,
                                                 marginal_reflected))
    plus = tuple(
        tuple((left + right) % PRIME for left, right in zip(row, image))
        for row, image in zip(generators, reflected)
    )
    minus = tuple(
        tuple((left - right) % PRIME for left, right in zip(row, image))
        for row, image in zip(generators, reflected)
    )
    target_stabilizers = []
    pointwise_solutions = []
    for multiplier in range(1, P):
        for shift in range(P):
            target = relation_transform(
                difference_transform(generators, P - 1, 0), multiplier, shift
            )
            if C.rank(generators + target, WIDTH) == C.rank(generators, WIDTH):
                target_stabilizers.append((multiplier, shift))
            for sign in (1, PRIME - 1):
                if reflected == scalar_rows(target, sign):
                    pointwise_solutions.append((
                        multiplier, shift, 1 if sign == 1 else -1
                    ))
    return {
        "point_reversal": reversal,
        "sharp_graph_rank": C.rank(graph, 2 * WIDTH),
        "sharp_graph_sha256": C.rowspace_digest(graph, 2 * WIDTH),
        "marginal_graph_rank": C.rank(marginal_graph, 2 * P),
        "marginal_graph_sha256": C.rowspace_digest(marginal_graph, 2 * P),
        "sharp_plus_minus_ranks": (C.rank(plus, WIDTH), C.rank(minus, WIDTH)),
        "target_difference_negation_relation_affine_stabilizers": (
            tuple(target_stabilizers)
        ),
        "pointwise_target_covariance_solutions": tuple(pointwise_solutions),
        "mu_graph_equivariant": mu_rows(reflected) == marginal_reflected,
    }


def main() -> None:
    k2, r4, reconstruction = C.reconstruct_relation_flag()
    (
        _gamma_actual, _gamma_support, gamma_pointed,
        point_cores, _support_cores, work_counts, state_counts, branch_counts,
    ) = F.build_banks()
    require(work_counts == F.EXPECTED_WORK_COUNTS, ("work counts", work_counts))
    require(state_counts == F.EXPECTED_STATE_COUNTS, ("state counts", state_counts))
    require(branch_counts == F.EXPECTED_BRANCH_COUNTS,
            ("branch counts", branch_counts))
    require(F.digest_json(gamma_pointed) == F.PS.EXPECTED_DIGESTS[0][0],
            "pointed gamma drift")

    zeta = pow(F.R.JOINT_ROOT, F.R.JOINT_ORDER // P, PRIME)
    require(zeta == EXPECTED_ZETA13, ("zeta", zeta))
    require(tuple(F.POINTS) == EXPECTED_POINTS, ("point order", F.POINTS))
    point_branch = F.invert_point_cores(point_cores, zeta)
    pointed = F.pointed_marginal(point_branch)
    pointed_digest = F.digest_json(pointed)
    require(pointed_digest == EXPECTED_POINTED_TENSOR_SHA256,
            ("pointed tensor", pointed_digest))

    generators = flatten_pointed(pointed)
    branch_resolved_rows = flatten_point_branch(point_branch)
    psharp = C.canonical_row_basis(generators, WIDTH)
    branch_basis = C.canonical_row_basis(branch_resolved_rows, WIDTH)
    require(branch_basis == psharp, "branch-resolved target image differs from branch sum")
    branch_boundary = {
        "branch_labeled_rows": len(branch_resolved_rows),
        "branch_resolved_rank": len(branch_basis),
        "branch_summed_rows": len(generators),
        "branch_summed_rank": len(psharp),
        "union_rank": C.rank(branch_resolved_rows + generators, WIDTH),
        "common_rowspace_sha256": C.rowspace_digest(branch_basis, WIDTH),
        "operation": "sum_branch_digit_before_P6sharp",
        "destroyed_sidecar": "point_by_branch coefficient/address labels",
    }

    marginal_generators = mu_rows(generators)
    p6 = C.canonical_row_basis(marginal_generators, P)
    qsharp_generators = centre_relation_sharp(generators)
    qsharp = C.canonical_row_basis(qsharp_generators, WIDTH)
    q6 = C.canonical_row_basis(mu_rows(qsharp_generators), P)
    require(C.canonical_row_basis(C.centre_relation(marginal_generators), P) == q6,
            "centering/marginal commutator")
    marginal_records = {
        "P6sharp_to_P6": (
            len(psharp), len(p6), len(psharp) - len(p6),
            C.rowspace_digest(tuple(row + mu_row(row) for row in psharp),
                              WIDTH + P),
        ),
        "Q6sharp_to_Q6": (
            len(qsharp), len(q6), len(qsharp) - len(q6),
            C.rowspace_digest(tuple(row + mu_row(row) for row in qsharp),
                              WIDTH + P),
        ),
    }
    require(marginal_records == EXPECTED_MARGINAL_RECORDS,
            ("marginal isomorphisms", marginal_records))

    rsharp = lift_basis(r4, qsharp, mu_rows(qsharp))
    ksharp = lift_basis(k2, qsharp, mu_rows(qsharp))
    spaces = {
        "K2sharp": ksharp,
        "R4sharp": rsharp,
        "P6sharp": psharp,
        "Q6sharp": qsharp,
    }
    space_ranks = {name: len(rows) for name, rows in spaces.items()}
    space_digests = {
        name: C.rowspace_digest(rows, WIDTH) for name, rows in spaces.items()
    }
    require(space_ranks == {"K2sharp": 2, "R4sharp": 4,
                            "P6sharp": 6, "Q6sharp": 6},
            ("sharp ranks", space_ranks))
    require(space_digests == EXPECTED_SPACE_DIGESTS,
            ("sharp digests", space_digests))
    require(C.rank(ksharp + rsharp, WIDTH) == 4, "K2sharp not in R4sharp")
    require(C.rank(rsharp + qsharp, WIDTH) == 6, "R4sharp not in Q6sharp")

    intersection_dimensions = {}
    intersection_digests = {}
    intersection_bases = {}
    for left_name, right_name in combinations(spaces, 2):
        key = left_name + "_inter_" + right_name
        basis = C.rowspace_intersection(spaces[left_name], spaces[right_name], WIDTH)
        intersection_bases[key] = basis
        intersection_dimensions[key] = len(basis)
        intersection_digests[key] = C.rowspace_digest(basis, WIDTH)
    require(intersection_dimensions == EXPECTED_INTERSECTION_DIMENSIONS,
            ("intersection dimensions", intersection_dimensions))
    require(intersection_digests == EXPECTED_INTERSECTION_DIGESTS,
            ("intersection digests", intersection_digests))
    require(intersection_bases["R4sharp_inter_P6sharp"] == ksharp,
            "R4sharp intersect P6sharp is not K2sharp")

    centered_p_basis = centre_relation_sharp(psharp)
    augmentation_images = tuple(
        tuple((raw - centered) % PRIME
              for raw, centered in zip(raw_row, centered_row))
        for raw_row, centered_row in zip(psharp, centered_p_basis)
    )
    require(all(
        all(len(set(row[difference * P:(difference + 1) * P])) == 1
            for difference in range(P))
        for row in augmentation_images
    ), "augmentation is not relation-constant")
    augmentation_basis = C.canonical_row_basis(augmentation_images, WIDTH)
    augmentation_difference_profiles = C.canonical_row_basis(
        tuple(tuple(row[difference * P] for difference in range(P))
              for row in augmentation_basis), P
    )
    marginal_augmentation = mu_rows(augmentation_basis)
    augmentation_marginal_kernel = linear_map_kernel(
        augmentation_basis, marginal_augmentation, P
    )
    lambda_on_rsharp = tuple(
        C.combine_rows(C.coordinates_in_row_basis(centered_p_basis, row),
                       augmentation_images)
        for row in rsharp
    )
    lambda_on_ksharp = tuple(
        C.combine_rows(C.coordinates_in_row_basis(centered_p_basis, row),
                       augmentation_images)
        for row in ksharp
    )
    lambda_r_kernel = linear_map_kernel(rsharp, lambda_on_rsharp, WIDTH)
    marginal_lambda_on_rsharp = mu_rows(lambda_on_rsharp)
    marginal_lambda_r_kernel = linear_map_kernel(
        rsharp, marginal_lambda_on_rsharp, P
    )

    point_lookup = {point: index for index, point in enumerate(F.POINTS)}
    point_reversal = tuple(
        point_lookup[(state ^ 2, P - 1 - root)] for state, root in F.POINTS
    )
    centered_generators = centre_relation_sharp(generators)
    augmentation_generators = tuple(
        tuple((raw - centered) % PRIME
              for raw, centered in zip(raw_row, centered_row))
        for raw_row, centered_row in zip(generators, centered_generators)
    )
    augmentation_reflected = tuple(
        augmentation_generators[index] for index in point_reversal
    )
    augmentation_plus = tuple(
        tuple((left + right) % PRIME for left, right in zip(row, reflected))
        for row, reflected in zip(augmentation_generators,
                                  augmentation_reflected)
    )
    augmentation_minus = tuple(
        tuple((left - right) % PRIME for left, right in zip(row, reflected))
        for row, reflected in zip(augmentation_generators,
                                  augmentation_reflected)
    )
    augmentation = {
        "sharp_rank": len(augmentation_basis),
        "sharp_sha256": C.rowspace_digest(augmentation_basis, WIDTH),
        "difference_profile_rank": len(augmentation_difference_profiles),
        "difference_profile_sha256": C.rowspace_digest(
            augmentation_difference_profiles, P
        ),
        "marginal_rank": C.rank(marginal_augmentation, P),
        "marginal_sha256": C.rowspace_digest(marginal_augmentation, P),
        "marginal_kernel_dimension": len(augmentation_marginal_kernel),
        "marginal_kernel_sha256": C.rowspace_digest(
            augmentation_marginal_kernel, WIDTH
        ),
        "lambda_on_R4sharp_rank": C.rank(lambda_on_rsharp, WIDTH),
        "lambda_on_R4sharp_image_sha256": C.rowspace_digest(
            lambda_on_rsharp, WIDTH
        ),
        "lambda_on_R4sharp_kernel_dimension": len(lambda_r_kernel),
        "lambda_on_R4sharp_kernel_sha256": C.rowspace_digest(
            lambda_r_kernel, WIDTH
        ),
        "lambda_on_K2sharp_rank": C.rank(lambda_on_ksharp, WIDTH),
        "marginal_lambda_on_R4sharp_rank": C.rank(
            marginal_lambda_on_rsharp, P
        ),
        "marginal_lambda_on_R4sharp_kernel_dimension": len(
            marginal_lambda_r_kernel
        ),
        "marginal_lambda_on_R4sharp_kernel_sha256": C.rowspace_digest(
            marginal_lambda_r_kernel, WIDTH
        ),
        "reversal": {
            "sharp_plus_minus_ranks": (
                C.rank(augmentation_plus, WIDTH),
                C.rank(augmentation_minus, WIDTH),
            ),
            "marginal_plus_minus_ranks": (
                C.rank(mu_rows(augmentation_plus), P),
                C.rank(mu_rows(augmentation_minus), P),
            ),
            "plus_sha256": C.rowspace_digest(augmentation_plus, WIDTH),
            "minus_sha256": C.rowspace_digest(augmentation_minus, WIDTH),
        },
    }
    require(augmentation == EXPECTED_AUGMENTATION,
            ("sharp augmentation", augmentation))
    require(lambda_r_kernel == ksharp, "sharp augmentation kernel is not K2sharp")

    # Hostile control: the duplicate-state coefficient differences are a
    # concrete rank-two map, but they are not the abstract augmentation map.
    # Coefficients use the pinned six-point order and the centered point rows.
    require(C.rank(centered_generators, WIDTH) == 6,
            "centered point generators are not a coefficient frame")
    rsharp_point_coefficients = tuple(
        C.coordinates_in_row_basis(centered_generators, row) for row in rsharp
    )
    ksharp_point_coefficients = tuple(
        C.coordinates_in_row_basis(centered_generators, row) for row in ksharp
    )

    def duplicate_state_delta(coefficients):
        require(len(coefficients) == 6, "duplicate-state coefficient arity")
        return (
            (coefficients[1] - coefficients[2]) % PRIME,
            (coefficients[3] - coefficients[4]) % PRIME,
        )

    delta_on_rsharp = tuple(
        duplicate_state_delta(coefficients)
        for coefficients in rsharp_point_coefficients
    )
    delta_on_ksharp = tuple(
        duplicate_state_delta(coefficients)
        for coefficients in ksharp_point_coefficients
    )
    delta_kernel = linear_map_kernel(rsharp, delta_on_rsharp, 2)
    delta_kernel_inter_ksharp = C.rowspace_intersection(
        delta_kernel, ksharp, WIDTH
    )
    lambda_coordinates = tuple(
        C.coordinates_in_row_basis(augmentation_basis, row)
        for row in lambda_on_rsharp
    )
    delta_lambda_joint = tuple(
        tuple(delta) + tuple(lambda_coordinate)
        for delta, lambda_coordinate in zip(delta_on_rsharp, lambda_coordinates)
    )
    pattern_residuals = (
        augmentation_generators[0],
        augmentation_generators[5],
        tuple((left + right) % PRIME for left, right in zip(
            augmentation_generators[1], augmentation_generators[2]
        )),
        tuple((left + right) % PRIME for left, right in zip(
            augmentation_generators[3], augmentation_generators[4]
        )),
    )
    pattern_flags = (
        not any(augmentation_generators[0]),
        all((left + right) % PRIME == 0 for left, right in zip(
            augmentation_generators[1], augmentation_generators[2]
        )),
        all((left + right) % PRIME == 0 for left, right in zip(
            augmentation_generators[3], augmentation_generators[4]
        )),
        not any(augmentation_generators[5]),
    )
    duplicate_state_hostile = {
        "point_order": tuple(F.POINTS),
        "delta_definition": "(c1-c2,c3-c4)",
        "point_correction_pattern": "(0,+H1,-H1,+H2,-H2,0)",
        "point_correction_pattern_flags": pattern_flags,
        "point_correction_pattern_holds": all(pattern_flags),
        "point_corrections_ordered_sha256": digest_json(augmentation_generators),
        "pattern_residual_rank": C.rank(pattern_residuals, WIDTH),
        "pattern_residual_ordered_sha256": digest_json(pattern_residuals),
        "pattern_residual_rowspace_sha256": C.rowspace_digest(
            pattern_residuals, WIDTH
        ),
        "delta_on_R4sharp_rank": C.rank(delta_on_rsharp, 2),
        "delta_on_R4sharp_sha256": digest_json(delta_on_rsharp),
        "delta_graph_sha256": C.rowspace_digest(
            tuple(tuple(row) + tuple(image)
                  for row, image in zip(rsharp, delta_on_rsharp)),
            WIDTH + 2,
        ),
        "delta_kernel_dimension": len(delta_kernel),
        "delta_kernel_sha256": C.rowspace_digest(delta_kernel, WIDTH),
        "delta_on_K2sharp_rank": C.rank(delta_on_ksharp, 2),
        "delta_kernel_inter_K2sharp_dimension": len(delta_kernel_inter_ksharp),
        "delta_kernel_inter_K2sharp_sha256": C.rowspace_digest(
            delta_kernel_inter_ksharp, WIDTH
        ),
        "delta_kernel_equals_K2sharp": delta_kernel == ksharp,
        "lambda_coordinate_rank": C.rank(lambda_coordinates, 2),
        "delta_lambda_joint_rank": C.rank(delta_lambda_joint, 4),
        "lambda_factors_through_delta": (
            C.rank(delta_lambda_joint, 4) == C.rank(delta_on_rsharp, 2)
        ),
        "identification_scope": "A2sharp is abstract, not duplicate-state physical",
    }
    require(duplicate_state_hostile["delta_on_R4sharp_rank"] == 2,
            ("Delta rank", duplicate_state_hostile))
    require(duplicate_state_hostile["delta_kernel_dimension"] == 2,
            ("Delta kernel dimension", duplicate_state_hostile))
    require(not duplicate_state_hostile["delta_kernel_equals_K2sharp"],
            ("false Delta kernel identification", duplicate_state_hostile))
    require(not duplicate_state_hostile["lambda_factors_through_delta"],
            ("false LambdaSharp factorization", duplicate_state_hostile))
    require(not duplicate_state_hostile["point_correction_pattern_holds"],
            ("false point-correction pattern", duplicate_state_hostile))
    require(duplicate_state_hostile == EXPECTED_DUPLICATE_STATE_HOSTILE,
            ("duplicate-state hostile", duplicate_state_hostile))

    difference_slice_ranks = {
        name: physical_slice_ranks(rows) for name, rows in spaces.items()
    }
    relation_slice_rank_record = {
        name: relation_slice_ranks(rows) for name, rows in spaces.items()
    }
    require(difference_slice_ranks == EXPECTED_DIFFERENCE_SLICE_RANKS,
            ("difference slices", difference_slice_ranks))
    require(relation_slice_rank_record == EXPECTED_RELATION_SLICE_RANKS,
            ("relation slices", relation_slice_rank_record))
    same_root = {
        name: (
            difference_slice_ranks[name][0],
            int(any(row[relation] for row in rows for relation in range(P))),
        )
        for name, rows in spaces.items()
    }
    require(set(same_root.values()) == {(0, 0)}, ("same root", same_root))
    physical_support_records = {
        name: (sum(physical_support(rows)), digest_json(physical_support(rows)))
        for name, rows in spaces.items()
    }
    require(set(physical_support_records.values()) == {
        (156, "9aac9a5917414fdaf88ee5fd53c564109ab357fc6b41dd21404feec1cceb9b31")
    }, ("physical support", physical_support_records))

    fourier = {name: fourier_record(rows, zeta) for name, rows in spaces.items()}
    require(fourier == EXPECTED_FOURIER, ("Fourier ledger", fourier))
    reversal = reversal_record(generators, marginal_generators)
    require(reversal == EXPECTED_REVERSAL, ("reversal ledger", reversal))

    relation_flag_atlas = affine_atlas(rsharp, ksharp, qsharp, "relation")
    difference_flag_atlas = affine_atlas(rsharp, ksharp, qsharp, "difference")
    relation_self_atlas = self_affine_atlas(psharp, qsharp, "relation")
    difference_self_atlas = self_affine_atlas(psharp, qsharp, "difference")
    affine = {
        "relation_flag_histogram": histogram(relation_flag_atlas, 2),
        "relation_flag_nonzero": tuple(
            record for record in relation_flag_atlas if record[2:4] != (0, 0)
        ),
        "relation_flag_sha256": digest_json(relation_flag_atlas),
        "difference_flag_histogram": histogram(difference_flag_atlas, 2),
        "difference_flag_sha256": digest_json(difference_flag_atlas),
        "relation_self_histogram": histogram(relation_self_atlas, 2),
        "relation_self_stabilizers": tuple(
            record[:2] for record in relation_self_atlas if record[2:] == (6, 6)
        ),
        "relation_self_sha256": digest_json(relation_self_atlas),
        "difference_self_histogram": histogram(difference_self_atlas, 2),
        "difference_self_stabilizers": tuple(
            record[:2] for record in difference_self_atlas if record[2:] == (6, 6)
        ),
        "difference_self_sha256": digest_json(difference_self_atlas),
        "relation_and_difference_flag_atlases_equal": (
            relation_flag_atlas == difference_flag_atlas
        ),
        "relation_and_difference_self_atlases_equal": (
            relation_self_atlas == difference_self_atlas
        ),
        "coordinate_swap": {
            name: (
                intersection_dimension(rows, transpose_coordinates(rows)),
                C.rank(rows + transpose_coordinates(rows), WIDTH),
                C.rowspace_digest(transpose_coordinates(rows), WIDTH),
            )
            for name, rows in spaces.items()
        },
    }
    require(affine == EXPECTED_AFFINE, ("affine ledger", affine))

    semantic_record = (
        PRIME, P, V, zeta, tuple(F.POINTS),
        tuple(sorted(EXPECTED_PARENT_SHA256.items())), reconstruction,
        pointed_digest, tuple(sorted(branch_boundary.items())),
        tuple(sorted(marginal_records.items())),
        tuple(sorted(space_ranks.items())), tuple(sorted(space_digests.items())),
        tuple(sorted(intersection_dimensions.items())),
        tuple(sorted(intersection_digests.items())),
        tuple(sorted(augmentation.items())),
        tuple(sorted(duplicate_state_hostile.items())),
        tuple(sorted(difference_slice_ranks.items())),
        tuple(sorted(relation_slice_rank_record.items())),
        tuple(sorted(same_root.items())),
        tuple(sorted(physical_support_records.items())),
        tuple((name, tuple(sorted(record.items())))
              for name, record in sorted(fourier.items())),
        tuple(sorted(reversal.items())), tuple(sorted(affine.items())),
        "branch digit summed; target image retained; coefficient/address labels lost",
        "static finite-field sharp carrier; no chronology/current/row/LRC14",
    )
    semantic = digest_json(semantic_record)
    require(semantic == EXPECTED_SEMANTIC_SHA256,
            ("semantic", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== THM-3615 pointed root-difference lift flag ==")
    print(f"field=(prime={PRIME},difference={P},relation={P},sharp_ambient={WIDTH})")
    print(f"parent_sha256_lf={EXPECTED_PARENT_SHA256}")
    print(f"points={tuple(F.POINTS)};zeta13={zeta}")
    print(f"pointed_tensor_sha256={pointed_digest}")
    print(f"branch_summed_boundary={branch_boundary}")
    print(f"marginal_isomorphisms={marginal_records}")
    print(f"space_ranks={space_ranks}")
    print(f"space_digests={space_digests}")
    print(f"intersection_dimensions={intersection_dimensions}")
    print(f"intersection_digests={intersection_digests}")
    print(f"sharp_augmentation={augmentation}")
    print(f"duplicate_state_physical_hostile={duplicate_state_hostile}")
    print(f"difference_slice_ranks={difference_slice_ranks}")
    print(f"relation_slice_ranks={relation_slice_rank_record}")
    print(f"same_root_s0_(slice_rank,nonzero)={same_root}")
    print(f"physical_support_(count,digest)={physical_support_records}")
    print("fourier_convention=sum_(s,t) f(s,t) zeta^(a*s+b*t)")
    print(f"fourier_records={fourier}")
    print(f"reversal_records={reversal}")
    print(f"affine_records={affine}")
    print("flag=K2sharp(2)<R4sharp(4)<Q6sharp(6)<Fun(F13xF13)(169)")
    print("raw_intersection=R4sharp_inter_P6sharp=K2sharp;dimension=2")
    print(f"semantic_sha256={semantic}")
    print(f"script_sha256_lf={lf_sha256(Path(__file__).resolve())}")
    print("status=FINITE-EXACT+VERIFIED-EXACT companion;independent hostile audit pending")
    print("boundary=branch digit summed before P6sharp;branch target image survives but point-by-branch coefficient/address intertwiner is absent")
    print("augmentation_scope=A2sharp is abstract;Delta(c)=(c1-c2,c3-c4) is not LambdaSharp")
    print("scope=static pinned finite field;not chronology/current/entry/row/characteristic-zero/LRC14")
    print("PASS")


if __name__ == "__main__":
    main()
