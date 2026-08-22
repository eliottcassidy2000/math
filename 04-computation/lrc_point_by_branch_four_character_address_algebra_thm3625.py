#!/usr/bin/env python3
"""Exact point-by-branch address algebra above THM-3615.

Retain the branch digit in the clean-room four-way tensor used by THM-3615.
For each branch r, express its six point rows in the ordered basis of the six
branch-summed rows.  The resulting 6 by 6 matrices A_r form a commutative
four-dimensional split semisimple algebra over the pinned finite field.

The six-dimensional coefficient module has four characters of multiplicities
2,1,1,2.  The inherited K2sharp plane is one line in each doubled character,
whereas P6sharp intersect Q6sharp is the sum of the two doubled character
spaces.  Consequently branch addresses alone cannot recover K2sharp as the
kernel or image of an algebra element.  This is a finite-field, static
obstruction: it says nothing about chronology, current, characteristic-zero
transfer, row exclusion, or LRC(14).

There are deliberately no Python ``assert`` statements.  Every mathematical
gate remains active under ``python -O``.
"""

from __future__ import annotations

from hashlib import sha256
import json
import os
from pathlib import Path
import sys

# SymPy's python-flint factor sorter rejects this 80-bit modulus on some
# installations; use the exact pure-Python ground domain for degree six.
os.environ["SYMPY_GROUND_TYPES"] = "python"
import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
CANON = ROOT / "01-canon/theorems"
COMPUTATION = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge/results"
sys.path.insert(0, str(COMPUTATION))

import lrc_pointed_difference_lift_flag_thm3615 as M


PARENT_THEOREM = CANON / "THM-3615-lrc-pointed-root-difference-lift-flag.md"
PARENT_SCRIPT = COMPUTATION / "lrc_pointed_difference_lift_flag_thm3615.py"
PARENT_OUTPUT = RESULTS / "lrc_pointed_difference_lift_flag_thm3615.out"

EXPECTED_PARENT_LF_SHA256 = {
    "theorem": "652b9628a256de9d4970294bc6185e8503464d94e34666ca0b6a971d7b367ec5",
    "script": "25a37a623f1d955f2e0b6723e0ab103bac32c9ccfabc1d0d9146123be0f6b7c7",
    "output": "493cf285d16881fbb7018724b47f44dbcbc183b984aa7602c30043628d0a5f40",
}

MOD = 755373809845391722745761
BRANCHES = 13
POINTS = 6
EXPECTED_ZETA13 = 266737884585332483769981
EXPECTED_ADDRESS_RANKS = (6,) * 13
EXPECTED_ADDRESS_SPAN_RANK = 4
EXPECTED_COMMUTATOR_HISTOGRAM = ((0, 169),)
EXPECTED_PRODUCT_EXTENSION_RANK = 4
EXPECTED_REVERSAL_COVARIANCE = ((12, 12, False), (12, 12, True))
EXPECTED_CYCLIC_GENERATOR = ("address", 3)
EXPECTED_MINIMAL_COEFFICIENTS = (
    1,
    533256819274739474231000,
    157450021680728426647686,
    613479327941808322893771,
    64677205047035279561501,
)
EXPECTED_EIGENVALUES = (
    189500145028362740071369,
    380076640875614845827929,
    565142856996729319560976,
    598144967360728788546009,
)
EXPECTED_EIGENSPACE_DIMENSIONS = (2, 1, 1, 2)
EXPECTED_K_EIGEN_INTERSECTIONS = (1, 0, 0, 1)
EXPECTED_PQ_EIGEN_INTERSECTIONS = (2, 0, 0, 2)
EXPECTED_K_COORDINATE_SHA256 = (
    "065ca6f0d61df266a5e55aab3d50e03f35ddc37a156373eabf6ec6c04cc8015f"
)
EXPECTED_PQ_COORDINATE_SHA256 = (
    "8d2da40cf11f96dd7edfb042ac9d7c87b20eed8596cc3b3de673979e43020936"
)
EXPECTED_ADDRESS_SHA256 = (
    "20afa4a3d65e92f7c343710f7e1b9ac163c55fa93ee3d4aa54fde76b347645e0"
)
EXPECTED_FOURIER_SHA256 = (
    "d66be14669f799d2ff227f32ad804c62e002e9faee7e74cb23d53811e8c507e6"
)
EXPECTED_FOURIER_POWER_COORDINATE_SHA256 = (
    "567b9d289ac73e4da5b9cba6cfbd7001ecfd89d744726c238a00dcc0461162af"
)

# Filled after the first deterministic run, then enforced in normal and -O.
EXPECTED_SEMANTIC_SHA256 = (
    "89d69b3eb2dff06f7ebc236d90be73555c397b39542a155e79786370884c06df"
)


CHECKS = 0


def require(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"gate failed: {label}")


def lf_bytes(path):
    return path.read_bytes().replace(b"\r\n", b"\n")


def lf_sha256(path):
    return sha256(lf_bytes(path)).hexdigest()


def digest(obj):
    return sha256(
        json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def rank(rows, width=None):
    if not rows:
        return 0
    if width is None:
        width = len(rows[0])
    return M.C.rank(
        tuple(tuple(value % MOD for value in row) for row in rows), width
    )


def zero_matrix():
    return tuple(tuple(0 for _ in range(POINTS)) for _ in range(POINTS))


def identity_matrix():
    return tuple(
        tuple(int(row == column) for column in range(POINTS))
        for row in range(POINTS)
    )


def matmul(left, right):
    return tuple(
        tuple(
            sum(left[i][k] * right[k][j] for k in range(POINTS)) % MOD
            for j in range(POINTS)
        )
        for i in range(POINTS)
    )


def matadd(left, right):
    return tuple(
        tuple((left[i][j] + right[i][j]) % MOD for j in range(POINTS))
        for i in range(POINTS)
    )


def matsub(left, right):
    return tuple(
        tuple((left[i][j] - right[i][j]) % MOD for j in range(POINTS))
        for i in range(POINTS)
    )


def scalar(matrix, value):
    return tuple(
        tuple(value * matrix[i][j] % MOD for j in range(POINTS))
        for i in range(POINTS)
    )


def transpose(matrix):
    return tuple(
        tuple(matrix[column][row] for column in range(POINTS))
        for row in range(POINTS)
    )


def flat(matrix):
    return tuple(value for row in matrix for value in row)


def matpow(matrix, exponent):
    answer = identity_matrix()
    base = matrix
    while exponent:
        if exponent & 1:
            answer = matmul(answer, base)
        base = matmul(base, base)
        exponent >>= 1
    return answer


def row_times_matrix(row, matrix):
    return tuple(
        sum(row[index] * matrix[index][column] for index in range(POINTS))
        % MOD
        for column in range(POINTS)
    )


def histogram(values):
    return tuple((value, values.count(value)) for value in sorted(set(values)))


def build_addresses():
    (_ga, _gs, gamma_pointed, point_cores, _sc, work, states, branches) = (
        M.F.build_banks()
    )
    require(work == M.F.EXPECTED_WORK_COUNTS, "parent work counts")
    require(states == M.F.EXPECTED_STATE_COUNTS, "parent state counts")
    require(branches == M.F.EXPECTED_BRANCH_COUNTS, "parent branch counts")
    require(
        M.F.digest_json(gamma_pointed) == M.F.PS.EXPECTED_DIGESTS[0][0],
        "parent pointed tensor",
    )
    zeta = pow(M.F.R.JOINT_ROOT, M.F.R.JOINT_ORDER // BRANCHES, MOD)
    require(zeta == EXPECTED_ZETA13, "zeta13")
    point_branch = M.F.invert_point_cores(point_cores, zeta)
    pointed = M.F.pointed_marginal(point_branch)
    generators = M.flatten_pointed(pointed)
    require(rank(generators) == POINTS, "branch-summed basis rank")

    addresses = []
    for branch in range(BRANCHES):
        rows = []
        for point in range(POINTS):
            row = tuple(
                point_branch[point][branch][difference][relation]
                for difference in range(BRANCHES)
                for relation in range(BRANCHES)
            )
            coordinates = M.C.coordinates_in_row_basis(generators, row)
            rows.append(tuple(value % MOD for value in coordinates))
        addresses.append(tuple(rows))
    return zeta, tuple(generators), tuple(addresses)


def main():
    parent_hashes = {
        "theorem": lf_sha256(PARENT_THEOREM),
        "script": lf_sha256(PARENT_SCRIPT),
        "output": lf_sha256(PARENT_OUTPUT),
    }
    require(parent_hashes == EXPECTED_PARENT_LF_SHA256, "parent LF hashes")

    zeta, generators, addresses = build_addresses()
    zero = zero_matrix()
    identity = identity_matrix()

    total = zero
    for address in addresses:
        total = matadd(total, address)
    require(total == identity, "sum_r A_r = identity")

    address_ranks = tuple(rank(address) for address in addresses)
    address_span_rank = rank(tuple(flat(address) for address in addresses), 36)
    commutator_ranks = tuple(
        rank(matsub(matmul(addresses[r], addresses[s]),
                    matmul(addresses[s], addresses[r])))
        for r in range(BRANCHES)
        for s in range(BRANCHES)
    )
    products = tuple(
        flat(matmul(left, right)) for left in addresses for right in addresses
    )
    product_extension_rank = rank(
        tuple(flat(address) for address in addresses) + products, 36
    )
    require(address_ranks == EXPECTED_ADDRESS_RANKS, "all addresses invertible")
    require(address_span_rank == EXPECTED_ADDRESS_SPAN_RANK, "address span")
    require(
        histogram(commutator_ranks) == EXPECTED_COMMUTATOR_HISTOGRAM,
        "address commutators",
    )
    require(
        product_extension_rank == EXPECTED_PRODUCT_EXTENSION_RANK,
        "product closure",
    )
    require(digest(addresses) == EXPECTED_ADDRESS_SHA256, "address digest")

    fourier = []
    for frequency in range(BRANCHES):
        transform = zero
        for branch, address in enumerate(addresses):
            transform = matadd(
                transform,
                scalar(address, pow(zeta, frequency * branch, MOD)),
            )
        fourier.append(transform)
    fourier = tuple(fourier)
    fourier_ranks = tuple(rank(matrix) for matrix in fourier)
    fourier_span_rank = rank(tuple(flat(matrix) for matrix in fourier), 36)
    require(fourier_span_rank == 4, "Fourier span")
    require(digest(fourier) == EXPECTED_FOURIER_SHA256, "Fourier digest")

    reversal = tuple(reversed(range(POINTS)))
    def reversal_conjugate(matrix):
        return tuple(
            tuple(matrix[reversal[i]][reversal[j]] for j in range(POINTS))
            for i in range(POINTS)
        )

    covariance = []
    for multiplier in range(1, BRANCHES):
        for shift in range(BRANCHES):
            for use_transpose in (False, True):
                if all(
                    reversal_conjugate(addresses[r])
                    == (
                        transpose(addresses[(multiplier * r + shift) % BRANCHES])
                        if use_transpose
                        else addresses[(multiplier * r + shift) % BRANCHES]
                    )
                    for r in range(BRANCHES)
                ):
                    covariance.append((multiplier, shift, use_transpose))
    covariance = tuple(covariance)
    all_symmetric = all(matrix == transpose(matrix) for matrix in addresses)
    require(covariance == EXPECTED_REVERSAL_COVARIANCE, "reversal covariance")
    require(all_symmetric, "pinned-gauge symmetry")

    candidates = tuple(
        ("address", index, matrix) for index, matrix in enumerate(addresses)
    ) + tuple(("fourier", index, matrix) for index, matrix in enumerate(fourier))
    generator_record = None
    for family, index, candidate in candidates:
        powers = tuple(matpow(candidate, exponent) for exponent in range(5))
        if rank(tuple(flat(value) for value in powers[:4]), 36) != 4:
            continue
        if rank(tuple(flat(value) for value in powers), 36) != 4:
            continue
        coordinates = tuple(
            M.C.coordinates_in_row_basis(
                tuple(flat(value) for value in powers[:4]), flat(powers[4])
            )
        )
        generator_record = (family, index, candidate, powers[:4], coordinates)
        break
    require(generator_record is not None, "cyclic algebra generator")
    family, generator_index, generator, power_basis, fourth_coordinates = (
        generator_record
    )
    require(
        (family, generator_index) == EXPECTED_CYCLIC_GENERATOR,
        "first cyclic generator",
    )

    variable = sp.symbols("T")
    c0, c1, c2, c3 = fourth_coordinates
    minimal = sp.Poly(
        variable**4 - c3 * variable**3 - c2 * variable**2
        - c1 * variable - c0,
        variable,
        modulus=MOD,
    )
    minimal_coefficients = tuple(int(value) % MOD for value in minimal.all_coeffs())
    require(
        minimal_coefficients == EXPECTED_MINIMAL_COEFFICIENTS,
        "minimal polynomial",
    )
    require(sp.gcd(minimal, minimal.diff()).degree() == 0, "minimal squarefree")
    factorization = sp.factor_list(minimal, modulus=MOD)
    roots = []
    for factor, multiplicity in factorization[1]:
        coefficients = tuple(int(value) % MOD for value in factor.all_coeffs())
        require(
            factor.degree() == 1 and coefficients[0] == 1 and multiplicity == 1,
            "linear squarefree factor",
        )
        roots.append((-coefficients[1]) % MOD)
    eigenvalues = tuple(sorted(roots))
    require(eigenvalues == EXPECTED_EIGENVALUES, "four split characters")

    eigenspaces = tuple(
        M.C.canonical_row_basis(
            M.C.nullspace(matsub(generator, scalar(identity, eigenvalue)), POINTS),
            POINTS,
        )
        for eigenvalue in eigenvalues
    )
    eigenspace_dimensions = tuple(len(space) for space in eigenspaces)
    require(
        eigenspace_dimensions == EXPECTED_EIGENSPACE_DIMENSIONS,
        "character multiplicities",
    )
    require(sum(eigenspace_dimensions) == POINTS, "spectral direct sum")

    address_power_coordinates = tuple(
        tuple(
            M.C.coordinates_in_row_basis(
                tuple(flat(value) for value in power_basis), flat(address)
            )
        )
        for address in addresses
    )
    fourier_power_coordinates = tuple(
        tuple(
            M.C.coordinates_in_row_basis(
                tuple(flat(value) for value in power_basis), flat(matrix)
            )
        )
        for matrix in fourier
    )
    require(
        digest(fourier_power_coordinates)
        == EXPECTED_FOURIER_POWER_COORDINATE_SHA256,
        "Fourier power coordinates",
    )

    k2, _r4, _reconstruction = M.C.reconstruct_relation_flag()
    psharp = M.C.canonical_row_basis(generators, M.WIDTH)
    qsharp = M.C.canonical_row_basis(
        M.centre_relation_sharp(generators), M.WIDTH
    )
    ksharp = M.lift_basis(k2, qsharp, M.mu_rows(qsharp))
    pqsharp = M.C.rowspace_intersection(psharp, qsharp, M.WIDTH)
    k_coordinates = M.C.canonical_row_basis(
        tuple(
            tuple(M.C.coordinates_in_row_basis(generators, row))
            for row in ksharp
        ),
        POINTS,
    )
    pq_coordinates = M.C.canonical_row_basis(
        tuple(
            tuple(M.C.coordinates_in_row_basis(generators, row))
            for row in pqsharp
        ),
        POINTS,
    )
    require(digest(k_coordinates) == EXPECTED_K_COORDINATE_SHA256, "K2sharp")
    require(digest(pq_coordinates) == EXPECTED_PQ_COORDINATE_SHA256, "P cap Q")

    k_invariant = all(
        rank(k_coordinates + tuple(row_times_matrix(row, address)
                                   for row in k_coordinates), POINTS)
        == len(k_coordinates)
        for address in addresses
    )
    pq_invariant = all(
        rank(pq_coordinates + tuple(row_times_matrix(row, address)
                                    for row in pq_coordinates), POINTS)
        == len(pq_coordinates)
        for address in addresses
    )
    require(k_invariant, "K2sharp address invariance")
    require(pq_invariant, "P cap Q address invariance")

    k_intersections = tuple(
        len(M.C.rowspace_intersection(k_coordinates, space, POINTS))
        for space in eigenspaces
    )
    pq_intersections = tuple(
        len(M.C.rowspace_intersection(pq_coordinates, space, POINTS))
        for space in eigenspaces
    )
    require(
        k_intersections == EXPECTED_K_EIGEN_INTERSECTIONS,
        "K2sharp character intersections",
    )
    require(
        pq_intersections == EXPECTED_PQ_EIGEN_INTERSECTIONS,
        "P cap Q character intersections",
    )
    require(
        M.C.canonical_row_basis(eigenspaces[0] + eigenspaces[3], POINTS)
        == pq_coordinates,
        "P cap Q equals doubled-character sum",
    )

    # Every algebra element is scalar on each eigenspace.  Its kernel and image
    # therefore contain either all or none of every eigenspace.  Since K meets
    # each doubled character in one line, it is neither such a kernel nor image.
    kernel_image_dimension_profiles = tuple(
        sorted(
            {
                sum(
                    eigenspace_dimensions[index]
                    for index in range(4)
                    if mask & (1 << index)
                )
                for mask in range(16)
            }
        )
    )
    require(kernel_image_dimension_profiles == (0, 1, 2, 3, 4, 5, 6),
            "possible spectral dimensions")
    multiplicity_debt = (
        k_intersections[0] not in (0, eigenspace_dimensions[0])
        and k_intersections[3] not in (0, eigenspace_dimensions[3])
    )
    require(multiplicity_debt, "K2sharp multiplicity debt")

    semantic_record = {
        "field": (MOD, BRANCHES, POINTS, zeta),
        "parent_hashes": parent_hashes,
        "address_ranks": address_ranks,
        "address_span_rank": address_span_rank,
        "commutator_histogram": histogram(commutator_ranks),
        "product_extension_rank": product_extension_rank,
        "fourier_ranks": fourier_ranks,
        "fourier_span_rank": fourier_span_rank,
        "reversal_covariance": covariance,
        "all_symmetric": all_symmetric,
        "cyclic_generator": (family, generator_index),
        "minimal_coefficients_mod_p": minimal_coefficients,
        "eigenvalues": eigenvalues,
        "eigenspace_dimensions": eigenspace_dimensions,
        "k_invariant": k_invariant,
        "pq_invariant": pq_invariant,
        "k_intersections": k_intersections,
        "pq_intersections": pq_intersections,
        "k_digest": digest(k_coordinates),
        "pq_digest": digest(pq_coordinates),
        "address_power_coordinates_digest": digest(address_power_coordinates),
        "fourier_power_coordinates_digest": digest(fourier_power_coordinates),
        "address_digest": digest(addresses),
        "fourier_digest": digest(fourier),
        "scope": "static pinned Fp; no chronology/current/row/char0/LRC14",
    }
    semantic = digest(semantic_record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256, "semantic digest")

    print("== THM-3625 point-by-branch four-character address algebra ==")
    print(f"field=(p={MOD},branches={BRANCHES},points={POINTS},zeta13={zeta})")
    print(f"parent_sha256_lf={parent_hashes}")
    print(f"address=(ranks={address_ranks},span={address_span_rank},commutators={histogram(commutator_ranks)},product_span={product_extension_rank})")
    print(f"fourier=(ranks={fourier_ranks},span={fourier_span_rank})")
    print(f"symmetry=(all_symmetric={all_symmetric},reversal_covariance={covariance})")
    print(f"algebra=(generator={(family, generator_index)},minimal_coefficients_mod_p={minimal_coefficients},eigenvalues={eigenvalues})")
    print(f"module=(multiplicities={eigenspace_dimensions},K_intersections={k_intersections},PQ_intersections={pq_intersections})")
    print(f"flags=(K_invariant={k_invariant},PQ_invariant={pq_invariant},PQ=E0+E3,K_not_kernel_or_image={multiplicity_debt})")
    print(f"digests=(address={digest(addresses)},fourier={digest(fourier)},K={digest(k_coordinates)},PQ={digest(pq_coordinates)},fourier_power={digest(fourier_power_coordinates)})")
    print(f"semantic_sha256={semantic}")
    print(f"script_sha256_lf={lf_sha256(Path(__file__).resolve())}")
    print(f"CHECKS={CHECKS}")
    print("scope=static pinned finite field; no chronology/current/row/characteristic-zero/LRC14")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
