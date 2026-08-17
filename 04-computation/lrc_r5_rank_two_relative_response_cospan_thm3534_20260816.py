#!/usr/bin/env python3
"""Audit the rank-two r5 middle-response cospan and its H1 boundary.

Frozen inputs:

* the exact r5 tent-location C4 hostile; and
* the compact exact K2/K3 bank used by the common ten-dimensional quotient.

The expensive endpoint integration is not repeated.  The script reconstructs
the six 13 by 13 parent blocks from the hash-pinned bank, isolates the two
middle oriented blocks, and performs all subsequent linear algebra exactly in
the split prime field.

The main object is deliberately relative.  The two exceptional source rows
and the dual middle-Q10 defect contract to different two-planes.  They become
the same representation only after quotienting the unique endpoint-supported
line at r0=6.  The chamber involution then gives a 1+1 representation, whose
natural twisted C4 H1 has dimension one rather than two.
"""

from __future__ import annotations

from hashlib import sha256
import json
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE_PATH = (
    ROOT / "05-knowledge/results/"
    "lrc_r5_third_digit_78_state_exact_bank_20260816.json"
)

EXPECTED_CERTIFICATE_SHA256 = (
    "472925c638de7ac90b1a7880184766f2e8acec0f8291e2c46cf747b37cd46712"
)
EXPECTED_PARENT_RAW_SHA256 = (
    "a227bc2f385d8a2eaecb27f317fa5ed66623c70938d8a97aba620298a8a7b61b"
)
EXPECTED_PARENT_SEMANTIC_SHA256 = (
    "3d1527fb4ce4931680e50d7135b9d1129c1816e3a9158645523e2728ddc71ec2"
)
EXPECTED_BANK_DIGESTS = (
    "185b2fb843d37a6e6f73be48375e40e93029ff507392df16f0058892184a1db2",
    "b5b0b3bceb926ab8176c0fb5c9fd57c4ad5fa87379b8e29731d7e7094f552263",
    "1a3a8c73c62a9a7293ed9b80337df7211b45db31362b8bf1a5d223a6f584bec6",
)
EXPECTED_SEMANTIC_SHA256 = (
    "57a44964888fdb0a9ca1c890abbe4950c6fa7130f7ab2800c0faa8a5d6a0212d"
)

P = 13
ARCS = 6
ARC_REVERSAL = (1, 0, 3, 2, 5, 4)
CHAMBER_REVERSAL = (5, 4, 3, 2, 1, 0)
POINTS = ((0, 0), (1, 0), (1, 6), (3, 6), (3, 12), (2, 12))
TENT = (12, 12, 9, 3, 0, 0)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":")
    ).encode("ascii")
    return sha256(encoded).hexdigest()


def rref(rows, prime: int, columns: int | None = None):
    matrix = [list(row) for row in rows]
    if columns is None:
        columns = len(matrix[0]) if matrix else 0
    pivot_row = 0
    pivots = []
    for column in range(columns):
        pivot = next(
            (row for row in range(pivot_row, len(matrix))
             if matrix[row][column] % prime),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        inverse = pow(matrix[pivot_row][column] % prime, -1, prime)
        matrix[pivot_row] = [value * inverse % prime
                             for value in matrix[pivot_row]]
        for row in range(len(matrix)):
            if row == pivot_row:
                continue
            factor = matrix[row][column] % prime
            if factor:
                matrix[row] = [
                    (left - factor * right) % prime
                    for left, right in zip(matrix[row], matrix[pivot_row])
                ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    basis = tuple(tuple(value % prime for value in row)
                  for row in matrix[:pivot_row])
    return basis, tuple(pivots)


def rank(rows, prime: int, columns: int | None = None) -> int:
    return len(rref(rows, prime, columns)[0])


def nullspace(rows, prime: int, columns: int | None = None):
    if columns is None:
        columns = len(rows[0]) if rows else 0
    basis, pivots = rref(rows, prime, columns)
    free_columns = [column for column in range(columns)
                    if column not in pivots]
    result = []
    for free in free_columns:
        vector = [0] * columns
        vector[free] = 1
        for row, pivot in enumerate(pivots):
            vector[pivot] = -basis[row][free] % prime
        result.append(tuple(vector))
    return tuple(result)


def transpose(matrix):
    if not matrix:
        return ()
    return tuple(tuple(matrix[row][column] for row in range(len(matrix)))
                 for column in range(len(matrix[0])))


def inverse(matrix, prime: int):
    size = len(matrix)
    require(all(len(row) == size for row in matrix), "inverse square")
    augmented = [
        list(row) + [1 if i == j else 0 for j in range(size)]
        for i, row in enumerate(matrix)
    ]
    reduced, pivots = rref(augmented, prime, size)
    require(pivots == tuple(range(size)), ("inverse singular", pivots))
    return tuple(tuple(row[size:]) for row in reduced)


def row_times_matrix(row, matrix, prime: int):
    return tuple(
        sum(row[index] * matrix[index][column]
            for index in range(len(row))) % prime
        for column in range(len(matrix[0]))
    )


def add(left, right, prime: int):
    return tuple((x + y) % prime for x, y in zip(left, right))


def subtract(left, right, prime: int):
    return tuple((x - y) % prime for x, y in zip(left, right))


def scale(scalar: int, vector, prime: int):
    return tuple(scalar * value % prime for value in vector)


def reverse(vector):
    return tuple(reversed(vector))


def chamber_response(vector, prime: int):
    """Arc-odd chamber action: reverse the digit and reverse orientation."""

    return tuple(-value % prime for value in reversed(vector))


def coordinates(target, basis_rows, prime: int):
    """Coordinates of target in an independent row basis."""

    reduced, pivots = rref(basis_rows, prime, len(target))
    require(len(reduced) == len(basis_rows), "coordinate basis dependent")
    square = tuple(tuple(row[column] for column in pivots)
                   for row in basis_rows)
    square_inverse = inverse(square, prime)
    target_pivots = tuple(target[column] for column in pivots)
    answer = row_times_matrix(target_pivots, square_inverse, prime)
    rebuilt = tuple(
        sum(answer[row] * basis_rows[row][column]
            for row in range(len(basis_rows))) % prime
        for column in range(len(target))
    )
    require(rebuilt == target, ("coordinate reconstruction", target, rebuilt))
    return answer


def load_certificate():
    require(CERTIFICATE_PATH.is_file(), ("missing certificate", CERTIFICATE_PATH))
    certificate_sha = lf_sha256(CERTIFICATE_PATH)
    require(certificate_sha == EXPECTED_CERTIFICATE_SHA256,
            ("certificate drift", certificate_sha))
    certificate = json.loads(CERTIFICATE_PATH.read_text(encoding="ascii"))
    require(certificate["schema"] == "r5-third-digit-78-diagonal-bank-v1",
            "certificate schema")
    require(certificate["parent_raw_sha256"] == EXPECTED_PARENT_RAW_SHA256,
            "parent raw hash")
    require(certificate["parent_semantic_sha256"]
            == EXPECTED_PARENT_SEMANTIC_SHA256, "parent semantic hash")
    require(tuple(certificate["bank_digests"]) == EXPECTED_BANK_DIGESTS,
            "bank digests")
    require(tuple(tuple(point) for point in certificate["points"]) == POINTS,
            "point order")
    return certificate, certificate_sha


def parent_blocks(k2, prime: int):
    """Return blocks with rows r0 and columns r1."""

    return tuple(
        tuple(tuple(k2[arc][r0][r1] % prime for r1 in range(P))
              for r0 in range(P))
        for arc in range(ARCS)
    )


def cycle_primitive(cochain, prime: int):
    """Prefix primitive for a seam-zero edge cochain on oriented C_13."""

    require(sum(cochain) % prime == 0, "cycle seam")
    potential = [0]
    for edge in range(P - 1):
        potential.append((potential[-1] + cochain[edge]) % prime)
    rebuilt = tuple(
        (potential[(edge + 1) % P] - potential[edge]) % prime
        for edge in range(P)
    )
    require(rebuilt == tuple(cochain), ("cycle primitive", rebuilt, cochain))
    return tuple(potential)


def analyse(certificate, certificate_sha: str):
    prime = certificate["prime"]
    require(prime > 13 and gcd(prime, 26) == 1, "coefficient characteristic")
    k2 = certificate["k2"]
    k3 = certificate["k3"]
    blocks = parent_blocks(k2, prime)

    # Frozen tent-location input and its exact arc-odd middle support.
    tent_reversal = tuple(TENT[index] for index in ARC_REVERSAL)
    tent_even = tuple((left + right) // 2
                      for left, right in zip(TENT, tent_reversal))
    tent_odd = tuple((left - right) // 2
                     for left, right in zip(TENT, tent_reversal))
    require(tent_even == (12, 12, 6, 6, 0, 0), tent_even)
    require(tent_odd == (0, 0, 3, -3, 0, 0), tent_odd)

    # Reconstruct the Q10 parent grading and every fixed-r2 surjectivity gate.
    block_ranks = tuple(rank(transpose(block), prime, P) for block in blocks)
    block_codimensions = tuple(P - value for value in block_ranks)
    require(block_ranks == (11, 11, 12, 12, 11, 11), block_ranks)
    require(block_codimensions == (2, 2, 1, 1, 2, 2), block_codimensions)
    require(sum(block_ranks) == 68 and sum(block_codimensions) == 10,
            "Q10 dimensions")
    row_sum_values = tuple(
        tuple(sorted({sum(block[r0]) % prime for r0 in range(P)}))
        for block in blocks
    )
    require(row_sum_values == ((1,),) * ARCS, row_sum_values)
    constant_mode_quotient_ranks = tuple(
        rank(transpose(block) + ((1,) * P,), prime, P)
        - block_ranks[arc]
        for arc, block in enumerate(blocks)
    )
    require(constant_mode_quotient_ranks == (0,) * ARCS,
            constant_mode_quotient_ranks)
    child_quotient_profiles = []
    for r2 in range(P):
        profile = []
        for arc in range(ARCS):
            parent = transpose(blocks[arc])
            child = tuple(
                tuple(
                    k2[arc][r0][r1] * k3[arc][r0][r1][r2] % prime
                    for r0 in range(P)
                )
                for r1 in range(P)
            )
            profile.append(rank(parent + child, prime, P) - block_ranks[arc])
        child_quotient_profiles.append(tuple(profile))
    require(tuple(child_quotient_profiles)
            == (block_codimensions,) * P, "fixed-r2 Q10 surjectivity")

    # The only descended reflection is coupled chamber reflection.  Pure arc
    # reversal adds exactly two dimensions, both on the middle pair.
    pair_defects = []
    for left, right in ((0, 1), (2, 3), (4, 5)):
        left_rows = transpose(blocks[left])
        right_rows = transpose(blocks[right])
        union = rank(left_rows + right_rows, prime, P)
        pair_defects.append((left, right, union,
                             2 * union - block_ranks[left] - block_ranks[right]))
    require(tuple(pair_defects)
            == ((0, 1, 11, 0), (2, 3, 13, 2), (4, 5, 11, 0)),
            pair_defects)
    chamber_exact = all(
        blocks[CHAMBER_REVERSAL[arc]][P - 1 - r0][P - 1 - r1]
        == blocks[arc][r0][r1]
        for arc in range(ARCS) for r0 in range(P) for r1 in range(P)
    )
    require(chamber_exact, "chamber reflection")

    plus = blocks[2]   # B -> C
    minus = blocks[3]  # C -> B
    difference = tuple(
        tuple((plus[r0][r1] - minus[r0][r1]) % prime for r1 in range(P))
        for r0 in range(P)
    )
    expected_support = (
        (0, 11), (0, 12),
        (3, 1), (3, 2), (3, 3), (3, 4), (3, 5),
        (3, 7), (3, 8), (3, 9), (3, 10), (3, 11),
        (6, 5), (6, 7),
        (9, 1), (9, 2), (9, 3), (9, 4), (9, 5),
        (9, 7), (9, 8), (9, 9), (9, 10), (9, 11),
        (12, 0), (12, 1),
    )
    support = tuple((r0, r1) for r0 in range(P) for r1 in range(P)
                    if difference[r0][r1])
    require(support == expected_support, support)
    require(all(
        difference[P - 1 - r0][P - 1 - r1]
        == -difference[r0][r1] % prime
        for r0 in range(P) for r1 in range(P)
    ), "difference chamber anti-symmetry")
    difference_seams = tuple(sum(row) % prime for row in difference)
    require(difference_seams == (0,) * P, difference_seams)
    cycle_primitives = tuple(cycle_primitive(row, prime) for row in difference)
    require(len(cycle_primitives) == P, "digit-cycle primitives")

    response_rank = rank(difference, prime, P)
    source = (difference[3], difference[9])
    endpoint = (difference[0], difference[6], difference[12])
    endpoint_mid = (difference[6],)
    source_rank = rank(source, prime, P)
    endpoint_rank = rank(endpoint, prime, P)
    require((response_rank, source_rank, endpoint_rank,
             rank(source + endpoint, prime, P)) == (5, 2, 3, 5),
            "source/endpoint direct sum")
    nonroot = tuple(digit for digit in range(P) if digit not in (0, 6, 12))
    source_restricted = tuple(tuple(row[digit] for digit in nonroot)
                              for row in source)
    require(rank(source_restricted, prime, len(nonroot)) == 2,
            "r1-blind scalar hostile")
    require(tuple(len(set(row)) for row in source_restricted) == (8, 8),
            "source distinct amplitudes")
    require(all(value for row in source_restricted for value in row),
            "source support")

    # Q_mid^vee consists of the two left-kernel lines of the middle parent
    # blocks.  Link their scalings by the descended chamber reflection.
    left_plus_space = nullspace(transpose(plus), prime, P)
    left_minus_space = nullspace(transpose(minus), prime, P)
    require(len(left_plus_space) == len(left_minus_space) == 1,
            "middle dual quotient lines")
    left_plus = left_plus_space[0]
    left_minus = reverse(left_plus)
    require(rank((left_minus, left_minus_space[0]), prime, P) == 1,
            "chamber-linked minus line")
    require(rank((left_plus, left_minus), prime, P) == 2,
            "distinct Q middle dual lines")

    # Contraction by D gives a second, genuinely two-dimensional response
    # representation.  It is not the source-row plane.
    image_plus = row_times_matrix(left_plus, difference, prime)
    image_minus = row_times_matrix(left_minus, difference, prime)
    image = (image_plus, image_minus)
    require(rank(image, prime, P) == 2, "Q-defect contraction rank")
    require(chamber_response(image_plus, prime) == image_minus
            and chamber_response(image_minus, prime) == image_plus,
            "contracted chamber representation")

    source_image_union_rank = rank(source + image, prime, P)
    source_image_intersection = source_rank + 2 - source_image_union_rank
    require((source_image_union_rank, source_image_intersection) == (3, 1),
            "two distinct response planes")
    require(rank(image + endpoint, prime, P) == 5,
            "image/endpoint direct sum")
    si = source + image
    si_endpoint_intersection = (
        rank(si, prime, P) + endpoint_rank
        - rank(si + endpoint, prime, P)
    )
    require(si_endpoint_intersection == 1, "unique endpoint intersection")
    require(rank(si + endpoint_mid, prime, P) == rank(si, prime, P),
            "middle endpoint line lies in source+image")
    require(rank(source + endpoint_mid, prime, P) == 3
            and rank(image + endpoint_mid, prime, P) == 3,
            "relative quotient injectivity")

    # In the chamber-linked basis, the entire mismatch is one r0=6 endpoint
    # gauge.  The two scalar channels are both nonzero.
    marked_basis = (image_plus, image_minus, difference[6])
    coordinates_3 = coordinates(difference[3], marked_basis, prime)
    coordinates_9 = coordinates(difference[9], marked_basis, prime)
    a, b, c = coordinates_3
    require(coordinates_9 == (b, a, c),
            ("chamber-linked coordinates", coordinates_3, coordinates_9))
    determinant = (a * a - b * b) % prime
    require(a and b and c and determinant and (a + b) % prime
            and (a - b) % prime, "relative transition units")

    source_plus = add(difference[3], difference[9], prime)
    source_minus = subtract(difference[3], difference[9], prime)
    image_even = add(image_plus, image_minus, prime)
    image_odd = subtract(image_plus, image_minus, prime)
    require(source_minus == scale((a - b) % prime, image_odd, prime),
            "literal common anti-line")
    require(source_plus == add(
        scale((a + b) % prime, image_even, prime),
        scale(2 * c % prime, difference[6], prime), prime),
        "invariant endpoint gauge")
    require(rank((source_minus,), prime, P) == 1
            and rank(source + image + (source_minus,), prime, P) == 3,
            "common line nonzero")

    # The natural local-system monodromy is the chamber swap.  In odd
    # characteristic it has one invariant and one anti-invariant line.
    swap = ((0, 1), (1, 0))
    identity = ((1, 0), (0, 1))
    twisted_boundary = tuple(
        tuple((swap[row][column] - identity[row][column]) % prime
              for column in range(2))
        for row in range(2)
    )
    twisted_boundary_rank = rank(twisted_boundary, prime, 2)
    twisted_h0_dimension = 2 - twisted_boundary_rank
    twisted_h1_dimension = 2 - twisted_boundary_rank
    trivial_h1_dimension = 2
    require((twisted_boundary_rank, twisted_h0_dimension,
             twisted_h1_dimension, trivial_h1_dimension) == (1, 1, 1, 2),
            "twisted C4 cohomology")
    anti_vector = (1, -1 % prime)
    plus_vector = (1, 1)
    require(rank(twisted_boundary + (anti_vector,), prime, 2) == 1,
            "anti-line is exact")
    require(rank(twisted_boundary + (plus_vector,), prime, 2) == 2,
            "invariant line survives")

    # Stabilizing Q10 under pure arc reversal fills both middle hyperplanes,
    # so the repaired eight-dimensional quotient kills this dual two-space.
    middle_stabilized_rank = rank(
        transpose(plus) + transpose(minus), prime, P
    )
    require(middle_stabilized_rank == 13, "middle stabilization")
    arc_stable_middle_dual_dimension = 2 * (P - middle_stabilized_rank)
    require(arc_stable_middle_dual_dimension == 0,
            "Q_A8 erases middle dual")

    coefficient_no_go = (
        gcd(prime, 13),
        "Hom_Add(F_p^2,F_13)=Hom_Add(F_13,F_p^2)=0",
        "Hom_Add(F_p^2,char0_flux)=0",
        "no_unital_ring_map_F_p_to_F_13_or_F_2",
    )
    require(coefficient_no_go[0] == 1, coefficient_no_go)

    space_digests = (
        digest_json(rref(source, prime, P)[0]),
        digest_json(rref(endpoint, prime, P)[0]),
        digest_json(rref(image, prime, P)[0]),
        digest_json(rref(difference, prime, P)[0]),
        digest_json(rref((source_minus,), prime, P)[0]),
    )
    record = (
        certificate_sha,
        prime,
        TENT,
        tent_even,
        tent_odd,
        block_ranks,
        block_codimensions,
        row_sum_values,
        constant_mode_quotient_ranks,
        tuple(child_quotient_profiles),
        tuple(pair_defects),
        chamber_exact,
        support,
        difference_seams,
        digest_json(cycle_primitives),
        response_rank,
        source_rank,
        endpoint_rank,
        source_image_union_rank,
        source_image_intersection,
        si_endpoint_intersection,
        coordinates_3,
        coordinates_9,
        determinant,
        (twisted_boundary_rank, twisted_h0_dimension,
         twisted_h1_dimension, trivial_h1_dimension),
        arc_stable_middle_dual_dimension,
        coefficient_no_go,
        space_digests,
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))
    return record, semantic


def main() -> None:
    certificate, certificate_sha = load_certificate()
    record, semantic = analyse(certificate, certificate_sha)
    (
        cert_sha, prime, tent, tent_even, tent_odd,
        block_ranks, block_codimensions, row_sum_values,
        constant_mode_quotient_ranks, child_profiles, pair_defects,
        chamber_exact, support, difference_seams, primitive_digest,
        response_rank, source_rank, endpoint_rank,
        source_image_union_rank, source_image_intersection,
        si_endpoint_intersection, coordinates_3, coordinates_9, determinant,
        cohomology, arc_stable_middle_dual_dimension, coefficient_no_go,
        space_digests,
    ) = record
    print("== r=5 rank-two relative response cospan / twisted-H1 audit ==")
    print(f"frozen_bank=(certificate_sha256={cert_sha},prime={prime},parent_raw={EXPECTED_PARENT_RAW_SHA256},parent_semantic={EXPECTED_PARENT_SEMANTIC_SHA256})")
    print(f"tent=(h={tent},even={tent_even},odd={tent_odd},meaning=exception_location_not_amplitude)")
    print(f"Q10=(block_ranks={block_ranks},block_codimensions={block_codimensions},all_13_child_profiles_equal={len(set(child_profiles)) == 1},profile={child_profiles[0]})")
    print(f"digit_cycle=(parent_row_sums={row_sum_values},Q10_constant_mode_quotient_ranks={constant_mode_quotient_ranks},D_seams={difference_seams},primitive_digest={primitive_digest})")
    print(f"symmetry=(chamber_exact={chamber_exact},arc_reversal_pair_defects={pair_defects})")
    print(f"middle_difference=(support_size={len(support)},support={support},rank={response_rank})")
    print(f"source_endpoint=(source_rows=(3,9),source_rank={source_rank},endpoint_rows=(0,6,12),endpoint_rank={endpoint_rank},direct_sum_rank={response_rank})")
    print(f"two_planes=(source_rank=2,Qmid_dual_contraction_rank=2,union_rank={source_image_union_rank},intersection_rank={source_image_intersection})")
    print(f"minimal_relative_sidecar=(unique_endpoint_intersection_rank={si_endpoint_intersection},line=r0_6,relative_dimension=2)")
    print(f"marked_transition=(row3={coordinates_3},row9={coordinates_9},det={determinant},form=((a,b),(b,a)),endpoint_gauge=c*D6)")
    print(f"twisted_C4=(boundary_rank,H0,H1,trivial_H1)={cohomology};literal_common_line=chamber_anti_and_exact;survivor=chamber_invariant_requires_D6_quotient")
    print(f"arc_stable_QA8=(middle_dual_dimension={arc_stable_middle_dual_dimension},meaning=repair_kills_rank_two_interface)")
    print(f"coefficient_no_go={coefficient_no_go}")
    print(f"space_digests=(source,endpoint,Qcontraction,full_difference,literal_intersection)={space_digests}")
    print(f"semantic_sha256={semantic}")
    print("connection=Loc_{3,9}->source_response->V_rel<-Q10_middle_dual; both arrows are chamber-equivariant isomorphisms only modulo the marked endpoint line k*D_6")
    print("preserved=two reflected source response functions,chamber involution,middle orientation defect,relative zero/nonzero")
    print("lost=r0_6 endpoint amplitude,outer endpoint rows,absolute lift,Q10 child section,r1 chronology,closure edge,word-current semantics")
    print("sidecar=the actual r0_6 endpoint row plus a lawful same-copy D-to-A edge and coefficient/flux realization")
    print("tournament_xor=the intrinsic relation is a symmetric chamber involution (both-way),not a tournament;Boolean parity would erase the nonzero 2*c endpoint gauge and has no coefficient map")
    print("verdict=the scalar hostile has a minimal two-dimensional relative repair,but every response direction is digit-cycle exact,the natural chamber-twisted C4 H1 is only one-dimensional,and every additive route to F13 Kummer or characteristic-zero JC flux is zero")
    print("scope=FINITE-EXACT static good-reduction representation plus elementary local-system lemma;no physical current,no D5 flux map,no row exclusion,no LRC14")
    print("PASS")


if __name__ == "__main__":
    main()
