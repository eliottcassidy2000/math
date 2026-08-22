#!/usr/bin/env python3
"""Exact THM-3636 companion: arc reversal and middle-address quotient.

Starting from the promoted THM-3625 address algebra, this companion restores
the native bidirected-arc coefficient involution S=(01)(23)(45).  It verifies
that S splits the two doubled address characters, while an address-algebra
projector itself identifies R4sharp/K2sharp with the middle-point plane.  It
also retains the refuted endpoint-imbalance factorization as a hostile gate.

All claims are static over one pinned finite field.  There are deliberately no
Python ``assert`` statements; every gate remains active under ``python -O``.
"""

from __future__ import annotations

import ast
from hashlib import sha256
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
CANON = ROOT / "01-canon/theorems"
COMPUTATION = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge/results"
sys.path.insert(0, str(COMPUTATION))

import lrc_point_by_branch_four_character_address_algebra_thm3625 as A
import lrc_pointed_difference_lift_flag_thm3615 as M


MOD = A.MOD
POINTS = A.POINTS
WIDTH = M.WIDTH
INV2 = pow(2, -1, MOD)

PARENT_FILES = (
    (
        "THM3625_theorem",
        CANON / "THM-3625-lrc-point-by-branch-split-four-character-address-algebra.md",
        "7da43a28045b53d0243b0abfeb9451866d2a04af48f71e4758ee92a5e30add7d",
    ),
    (
        "THM3625_script",
        COMPUTATION / "lrc_point_by_branch_four_character_address_algebra_thm3625.py",
        "b32b985e2817f135dcc6d19f8f675d3d174e7f1b2a81bd14848d5424fabaf6ad",
    ),
    (
        "THM3625_output",
        RESULTS / "lrc_point_by_branch_four_character_address_algebra_thm3625.out",
        "8ef35eb1d7dd2f00c67db0171483ac1a9da8851eb21afb88bbdef5976a49b325",
    ),
    (
        "THM3615_theorem",
        CANON / "THM-3615-lrc-pointed-root-difference-lift-flag.md",
        "652b9628a256de9d4970294bc6185e8503464d94e34666ca0b6a971d7b367ec5",
    ),
    (
        "THM3615_script",
        COMPUTATION / "lrc_pointed_difference_lift_flag_thm3615.py",
        "25a37a623f1d955f2e0b6723e0ab103bac32c9ccfabc1d0d9146123be0f6b7c7",
    ),
    (
        "THM3615_output",
        RESULTS / "lrc_pointed_difference_lift_flag_thm3615.out",
        "493cf285d16881fbb7018724b47f44dbcbc183b984aa7602c30043628d0a5f40",
    ),
)

EXPECTED_DIGESTS = {
    "K": "065ca6f0d61df266a5e55aab3d50e03f35ddc37a156373eabf6ec6c04cc8015f",
    "R": "a3ca24b41fb507941de6e274058bf4658918a13d9ca9902e6483fb52be444446",
    "W": "8d2da40cf11f96dd7edfb042ac9d7c87b20eed8596cc3b3de673979e43020936",
    "delta_arc": "beed38852f054455adc87056afa51a4f131584c61dfe144dec3cdc84adc08946",
    "S": "02426ab4b994ba3cb3e8c6ba75fa47c66ba557f67f0be3bddaa3a94fd4ff337d",
    "Pi_W": "a012489eeab9b598118e47cfa78074926fd6c389601d724f718fb1923d872bf6",
    "Pi_K": "584c5dfef50dc58d6661ab1526d9f2bba0c8e1f26f029f57a36a600160bda5eb",
    "Pi_M": "6ed996b3aa6a81a9f774d833948add65d5164224b776344d1a3a128f298535e7",
    "T_profiles": "4af0bf56a78f56ef0af7d3cd6ff6d11b046d43421dc2e6006158b70616e3ede6",
    "T_expanded": "54e969e70b3b42d3ad5ecd4601ee63451bdfbcec18f368ef7f975b744d223ea3",
}

EXPECTED_PI_M_POWER_COORDINATES = (
    231956577553530055552921,
    258341843836556363678220,
    526705381605750340113858,
    221092595255266945805223,
)

EXPECTED_T_MID = (
    (126498113787680818370196, 646561219255993169342961),
    (384727693242765231857657, 452865069702498262026030),
)
EXPECTED_T_DETERMINANT = 275730381649850587765623

EXPECTED_SEMANTIC_SHA256 = "6d230ba65767bdb7080f86d391213c9ef0e6d0fcb7edf3b7e2a7ea6629bd5c73"

CHECKS = 0


def require(condition: bool, payload: object) -> None:
    global CHECKS
    if condition is not True:
        raise RuntimeError(payload)
    CHECKS += 1


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def standard(indices: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(int(index in indices) for index in range(POINTS))


def diagonal(values: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(values[row] if row == column else 0 for column in range(POINTS))
        for row in range(POINTS)
    )


def permutation_matrix(permutation: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(int(column == permutation[row]) for column in range(POINTS))
        for row in range(POINTS)
    )


def canonical(rows: tuple[tuple[int, ...], ...], width: int = POINTS):
    return M.C.canonical_row_basis(rows, width)


def matrix_coordinates(
    basis: tuple[tuple[tuple[int, ...], ...], ...],
    matrix: tuple[tuple[int, ...], ...],
) -> tuple[int, ...]:
    return tuple(
        value % MOD
        for value in M.C.coordinates_in_row_basis(
            tuple(A.flat(item) for item in basis), A.flat(matrix)
        )
    )


def restriction_matrix(
    space: tuple[tuple[int, ...], ...],
    matrix: tuple[tuple[int, ...], ...],
) -> tuple[tuple[int, ...], ...]:
    images = tuple(A.row_times_matrix(row, matrix) for row in space)
    require(
        A.rank(space + images, POINTS) == len(space),
        ("noninvariant restriction", len(space)),
    )
    return tuple(
        tuple(M.C.coordinates_in_row_basis(space, image)) for image in images
    )


def inverse_two(matrix: tuple[tuple[int, int], tuple[int, int]]):
    determinant = (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) % MOD
    require(determinant != 0, ("singular two by two", matrix))
    inverse = pow(determinant, -1, MOD)
    return (
        (matrix[1][1] * inverse % MOD, -matrix[0][1] * inverse % MOD),
        (-matrix[1][0] * inverse % MOD, matrix[0][0] * inverse % MOD),
    )


def multiply_two(left, right):
    return tuple(
        tuple(
            sum(left[row][inner] * right[inner][column] for inner in range(2)) % MOD
            for column in range(2)
        )
        for row in range(2)
    )


def solve_two_map(domain_rows, codomain_rows):
    pair = next(
        (i, j)
        for i in range(len(domain_rows))
        for j in range(i + 1, len(domain_rows))
        if (domain_rows[i][0] * domain_rows[j][1]
            - domain_rows[i][1] * domain_rows[j][0]) % MOD
    )
    domain_minor = (domain_rows[pair[0]], domain_rows[pair[1]])
    codomain_minor = (codomain_rows[pair[0]], codomain_rows[pair[1]])
    answer = multiply_two(inverse_two(domain_minor), codomain_minor)
    require(
        all(
            tuple(
                sum(domain[row] * answer[row][column] for row in range(2)) % MOD
                for column in range(2)
            ) == tuple(codomain)
            for domain, codomain in zip(domain_rows, codomain_rows, strict=True)
        ),
        "two-map factorization",
    )
    return answer


def main() -> None:
    parent_hashes = tuple((label, lf_sha256(path)) for label, path, _ in PARENT_FILES)
    require(
        parent_hashes == tuple((label, expected) for label, _, expected in PARENT_FILES),
        ("parent drift", parent_hashes),
    )

    _zeta, generators, addresses = A.build_addresses()
    identity = A.identity_matrix()
    zero = A.zero_matrix()
    generator = addresses[3]
    powers = tuple(A.matpow(generator, exponent) for exponent in range(4))
    require(A.rank(tuple(A.flat(item) for item in powers), 36) == 4, "power basis")

    # Reconstruct the raw and centered inherited flags in their two typed
    # six-point coefficient frames.
    k2, r4, _reconstruction = M.C.reconstruct_relation_flag()
    psharp = M.C.canonical_row_basis(generators, WIDTH)
    centered_generators = M.centre_relation_sharp(generators)
    qsharp = M.C.canonical_row_basis(centered_generators, WIDTH)
    ksharp = M.lift_basis(k2, qsharp, M.mu_rows(qsharp))
    rsharp = M.lift_basis(r4, qsharp, M.mu_rows(qsharp))
    wsharp = M.C.rowspace_intersection(psharp, qsharp, WIDTH)

    raw_coordinate_rows = lambda rows: tuple(
        tuple(M.C.coordinates_in_row_basis(generators, row)) for row in rows
    )
    centered_coordinate_rows = lambda rows: tuple(
        tuple(M.C.coordinates_in_row_basis(centered_generators, row)) for row in rows
    )

    # Keep the ordered coefficient rows for map computations: they correspond
    # row-for-row to ksharp/rsharp.  Canonicalization is lawful for comparing
    # subspaces, but would silently rebase the domain without applying the same
    # row operations to LambdaSharp's images.
    k_raw_rows = raw_coordinate_rows(ksharp)
    w_raw_rows = raw_coordinate_rows(wsharp)
    k_centered_rows = centered_coordinate_rows(ksharp)
    r_centered_rows = centered_coordinate_rows(rsharp)
    k_raw = canonical(k_raw_rows)
    w_raw = canonical(w_raw_rows)
    k_centered = canonical(k_centered_rows)
    r_centered = canonical(r_centered_rows)

    expected_k = canonical(((1, 1, 0, 0, 0, 0), (0, 0, 0, 0, 1, 1)))
    expected_r = canonical(expected_k + (standard((2,)), standard((3,))))
    expected_w = canonical(tuple(standard((index,)) for index in (0, 1, 4, 5)))
    expected_middle = canonical((standard((2,)), standard((3,))))
    require(k_raw == expected_k == k_centered, ("K coordinate form", k_raw, k_centered))
    require(r_centered == expected_r, ("R coordinate form", r_centered))
    require(w_raw == expected_w, ("W coordinate form", w_raw))
    require(digest(k_raw) == EXPECTED_DIGESTS["K"], "K digest")
    require(digest(r_centered) == EXPECTED_DIGESTS["R"], "R digest")
    require(digest(w_raw) == EXPECTED_DIGESTS["W"], "W digest")

    # The native bidirected-arc coefficient involution and three projectors.
    arc_permutation = (1, 0, 3, 2, 5, 4)
    arc = permutation_matrix(arc_permutation)
    pi_w = diagonal((1, 1, 0, 0, 1, 1))
    pi_m = A.matsub(identity, pi_w)
    pi_k = A.matmul(pi_w, A.scalar(A.matadd(identity, arc), INV2))
    require(A.matmul(arc, arc) == identity, "arc involution")
    require(A.matmul(pi_k, pi_k) == pi_k and A.rank(pi_k, POINTS) == 2, "Pi_K")
    require(canonical(tuple(A.row_times_matrix(row, pi_k) for row in identity)) == k_raw,
            "Pi_K image")
    require(canonical(tuple(A.row_times_matrix(row, pi_w) for row in identity)) == w_raw,
            "Pi_W image")
    require(digest(arc) == EXPECTED_DIGESTS["S"], "S digest")
    require(digest(pi_w) == EXPECTED_DIGESTS["Pi_W"], "Pi_W digest")
    require(digest(pi_k) == EXPECTED_DIGESTS["Pi_K"], "Pi_K digest")
    require(digest(pi_m) == EXPECTED_DIGESTS["Pi_M"], "Pi_M digest")

    pi_w_coordinates = matrix_coordinates(powers, pi_w)
    pi_m_coordinates = matrix_coordinates(powers, pi_m)
    require(pi_m_coordinates == EXPECTED_PI_M_POWER_COORDINATES,
            ("Pi_M power coordinates", pi_m_coordinates))
    require(
        tuple((int(index == 0) - pi_m_coordinates[index]) % MOD for index in range(4))
        == pi_w_coordinates,
        "Pi_W plus Pi_M",
    )

    # Spectral action and crossed-product algebra.
    eigenspaces = tuple(
        canonical(M.C.nullspace(
            A.matsub(generator, A.scalar(identity, eigenvalue)), POINTS
        ))
        for eigenvalue in A.EXPECTED_EIGENVALUES
    )
    require(tuple(len(space) for space in eigenspaces) == (2, 1, 1, 2), "multiplicities")
    arc_block_unions = tuple(
        A.rank(space + tuple(A.row_times_matrix(row, arc) for row in space), POINTS)
        for space in eigenspaces
    )
    require(arc_block_unions == (2, 2, 2, 2), ("arc character action", arc_block_unions))
    require(
        canonical(tuple(A.row_times_matrix(row, arc) for row in eigenspaces[1]))
        == eigenspaces[2]
        and canonical(tuple(A.row_times_matrix(row, arc) for row in eigenspaces[2]))
        == eigenspaces[1],
        "simple-character exchange",
    )
    doubled_sign_ranks = tuple(
        (
            A.rank(tuple(A.row_times_matrix(row, A.matadd(identity, arc)) for row in space), POINTS),
            A.rank(tuple(A.row_times_matrix(row, A.matsub(identity, arc)) for row in space), POINTS),
        )
        for space in (eigenspaces[0], eigenspaces[3])
    )
    require(doubled_sign_ranks == ((1, 1), (1, 1)), "doubled-character split")

    crossed_basis = powers + tuple(A.matmul(power, arc) for power in powers)
    crossed_rank = A.rank(tuple(A.flat(item) for item in crossed_basis), 36)
    require(crossed_rank == 8, "crossed algebra dimension")
    crossed_products = tuple(
        A.flat(A.matmul(left, right)) for left in crossed_basis for right in crossed_basis
    )
    require(
        A.rank(tuple(A.flat(item) for item in crossed_basis) + crossed_products, 36) == 8,
        "crossed product closure",
    )
    commutator_ranks = tuple(
        A.rank(A.matsub(A.matmul(left, right), A.matmul(right, left)), POINTS)
        for left in crossed_basis for right in crossed_basis
    )
    require(any(value for value in commutator_ranks), "crossed algebra noncommutative")
    require(all(A.matmul(arc, A.matmul(address, arc)) in tuple(
        A.matadd(zero, candidate) for candidate in addresses
    ) or A.rank(tuple(A.flat(item) for item in powers) +
                (A.flat(A.matmul(arc, A.matmul(address, arc))),), 36) == 4
                for address in addresses), "arc normalizes address algebra")

    middle_space = canonical(eigenspaces[1] + eigenspaces[2])
    block_span_ranks = (
        A.rank(tuple(A.flat(restriction_matrix(eigenspaces[0], item)) for item in crossed_basis), 4),
        A.rank(tuple(A.flat(restriction_matrix(middle_space, item)) for item in crossed_basis), 4),
        A.rank(tuple(A.flat(restriction_matrix(eigenspaces[3], item)) for item in crossed_basis), 4),
    )
    require(block_span_ranks == (2, 4, 2), ("crossed block algebras", block_span_ranks))
    require(all(A.matmul(pi_k, item) == A.matmul(item, pi_k) for item in crossed_basis),
            "Pi_K central")

    # Corrected quotient: endpoint imbalances vanish on all of R, whereas the
    # middle coordinates and Pi_M have kernel exactly K.
    delta_arc = tuple(
        ((row[0] - row[1]) % MOD, (row[4] - row[5]) % MOD)
        for row in r_centered_rows
    )
    delta_mid = tuple((row[2], row[3]) for row in r_centered_rows)
    require(A.rank(delta_arc, 2) == 0, ("refuted endpoint map", delta_arc))
    require(digest(delta_arc) == EXPECTED_DIGESTS["delta_arc"], "delta_arc digest")
    require(A.rank(delta_mid, 2) == 2, "middle map rank")
    delta_mid_kernel = M.linear_map_kernel(r_centered_rows, delta_mid, 2)
    require(delta_mid_kernel == k_centered, "middle map kernel")
    middle_images = canonical(tuple(A.row_times_matrix(row, pi_m) for row in r_centered_rows))
    require(middle_images == expected_middle, ("Pi_M middle image", middle_images))
    pi_m_images = tuple(A.row_times_matrix(row, pi_m) for row in r_centered_rows)
    pi_m_kernel = M.linear_map_kernel(r_centered_rows, pi_m_images, POINTS)
    require(pi_m_kernel == k_centered, "Pi_M kernel on R")
    require(
        canonical(r_centered + tuple(A.row_times_matrix(row, arc) for row in r_centered))
        == r_centered,
        "arc preserves R",
    )

    # Rebuild LambdaSharp and solve its unique factor through delta_mid.
    centered_p_basis = M.centre_relation_sharp(psharp)
    augmentation_images = tuple(
        tuple((raw - centered) % MOD for raw, centered in zip(raw_row, centered_row))
        for raw_row, centered_row in zip(psharp, centered_p_basis, strict=True)
    )
    augmentation_basis = M.C.canonical_row_basis(augmentation_images, WIDTH)
    lambda_on_r = tuple(
        M.C.combine_rows(
            M.C.coordinates_in_row_basis(centered_p_basis, row), augmentation_images
        )
        for row in rsharp
    )
    lambda_coordinates = tuple(
        tuple(M.C.coordinates_in_row_basis(augmentation_basis, row))
        for row in lambda_on_r
    )
    require(A.rank(lambda_coordinates, 2) == 2, "Lambda rank")
    t_mid = solve_two_map(delta_mid, lambda_coordinates)
    t_determinant = (
        t_mid[0][0] * t_mid[1][1] - t_mid[0][1] * t_mid[1][0]
    ) % MOD
    require(t_mid == EXPECTED_T_MID, ("T_mid", t_mid))
    require(t_determinant == EXPECTED_T_DETERMINANT, ("T determinant", t_determinant))

    expanded_t = tuple(M.C.combine_rows(row, augmentation_basis) for row in t_mid)
    profile_t = tuple(
        tuple(row[difference * M.P] for difference in range(M.P))
        for row in expanded_t
    )
    require(digest(profile_t) == EXPECTED_DIGESTS["T_profiles"],
            ("T profile digest", digest(profile_t)))
    require(digest(expanded_t) == EXPECTED_DIGESTS["T_expanded"],
            ("T expanded digest", digest(expanded_t)))

    semantic_record = {
        "field": (MOD, POINTS, tuple(M.F.POINTS)),
        "parents": parent_hashes,
        "coordinate_digests": {
            "K": digest(k_raw), "R": digest(r_centered), "W": digest(w_raw)
        },
        "operators": {
            "S": digest(arc), "Pi_W": digest(pi_w), "Pi_K": digest(pi_k),
            "Pi_M": digest(pi_m), "Pi_M_power_coordinates": pi_m_coordinates,
        },
        "crossed_algebra": {
            "dimension": crossed_rank,
            "commutator_histogram": A.histogram(commutator_ranks),
            "block_span_ranks": block_span_ranks,
            "arc_block_unions": arc_block_unions,
            "doubled_sign_ranks": doubled_sign_ranks,
        },
        "quotient": {
            "delta_arc": delta_arc,
            "delta_mid": delta_mid,
            "T_mid": t_mid,
            "T_determinant": t_determinant,
            "T_profile_digest": digest(profile_t),
            "T_expanded_digest": digest(expanded_t),
        },
        "scope": "static pinned Fp coefficient/address quotient; no chronology/current/char0/LRC14",
    }
    semantic = digest(semantic_record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic digest", semantic, EXPECTED_SEMANTIC_SHA256))

    source = Path(__file__).resolve()
    source_bytes = source.read_bytes()
    require(b"\r\n" not in source_bytes, "source raw LF")
    require(
        not any(isinstance(node, ast.Assert)
                for node in ast.walk(ast.parse(source_bytes.decode("utf-8")))),
        "Python assert node present",
    )

    print("== THM-3636 LRC arc reversal and middle-address quotient ==")
    print(f"field=(p={MOD},points={tuple(M.F.POINTS)})")
    print(f"parent_sha256_lf={parent_hashes}")
    print(f"coordinate_forms=(K={k_raw},R={r_centered},W={w_raw})")
    print(f"operators=(Pi_M_power={pi_m_coordinates},digests_S_W_K_M={(digest(arc),digest(pi_w),digest(pi_k),digest(pi_m))})")
    print(f"crossed_algebra=(dimension={crossed_rank},commutators={A.histogram(commutator_ranks)},blocks={block_span_ranks},arc_unions={arc_block_unions},sign_ranks={doubled_sign_ranks})")
    print(f"endpoint_hostile=(delta={delta_arc},rank={A.rank(delta_arc,2)},digest={digest(delta_arc)})")
    print(f"middle_quotient=(delta={delta_mid},rank={A.rank(delta_mid,2)},kernel=K,image={middle_images})")
    print(f"augmentation_factor=(T={t_mid},det={t_determinant},profile_digest={digest(profile_t)},expanded_digest={digest(expanded_t)})")
    print(f"semantic_sha256={semantic}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print(f"CHECKS={CHECKS};imported_address_checks={A.CHECKS}")
    print("scope=static pinned finite-field coefficient/address quotient;no chronology/current/characteristic-zero/LRC14")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
