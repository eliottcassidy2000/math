#!/usr/bin/env python3
"""Exact common-A4 comparison for the two pinned r=5 LRC tensors.

Reconstruct, from their hash-pinned parent programs,

* the two-current-digit tensor T_2(state,r0,r1,relation), and
* the nested source/current tensor T_s(channel,state,relation).

Flatten both into the common target

    H = Fun(V_4(state) x F_13(relation), F_p).

The rows of T_2 are indexed by (r0,r1), while the rows of T_s are indexed
by the seventeen retained source/current character pairs.  We compare their
canonical reduced-row-echelon bases over the pinned prime.  We then apply
double centering on the state and relation coordinates and repeat the test.

This is a static finite-field row-space identity.  Labels such as "source"
and "current" name the parent tensors only: no chronology, physical current,
entry theorem, row exclusion, characteristic-zero lift, or LRC(14) result is
constructed here.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
TWO_CURRENT_PATH = ROOT / (
    "04-computation/"
    "lrc_r5_ufull_owner_node_boolean_square_two_digit_current_ancestry_"
    "probe_20260816.py"
)
SOURCE_CURRENT_PATH = ROOT / (
    "04-computation/"
    "lrc_r5_ufull_owner_node_boolean_square_nested_ancestry_digits_"
    "probe_20260816.py"
)

P = 13
V = 4
EXPECTED_PRIME = 755373809845391722745761

EXPECTED_PARENT_SCRIPT_SHA256 = (
    "3dab580e479e4ba7ac8801c1e5d8523018e0b3dc1c2176c072e7c609033eb6c8",
    "1188df8aa2a7a84c1e8ada5fc3cc8d3b839ece70298b94f1d94c9d440caa88f3",
)
EXPECTED_DEPENDENCY_SHA256 = (
    # two-current one-digit owner parent
    "ae1cf021ea23f325eeded42ff1dea8df903d837b1d7ab551289d62f7ab7a0348",
    # nested source-sheet parent
    "592aa0bce31f2da5d5e2ddff7f3ffe6f1398f3a07b5ce927e0d97c9fe309ae3b",
    # nested inverse-owner parent
    "ae1cf021ea23f325eeded42ff1dea8df903d837b1d7ab551289d62f7ab7a0348",
)
EXPECTED_PARENT_RECONSTRUCTION_SHA256 = (
    # two-current tensor, two-current tau core
    "53eb6e618d0669bdb27841a1800c46e32456bb6c6d3698c590ae0d5e68822033",
    "fd1b837d9de3e4f9e586d29b69ed6726364ce97535d2e48f441c6ccd694250de",
    # source-current gamma bank, source-current tensor
    "9d0aa6823ed9bb83b338350d9da3d86c52e552a844b4334e3571bf2393f04cd1",
    "53e1d8e6de69dcd1cf5740acce7b245f7d9a388cfa2b27da4aefb00f3075df14",
)
EXPECTED_RAW_RREF_SHA256 = (
    "1d9293d05fa3551b785a1537e78bc8be585fcc43dbb5172036c9b32546ca8560"
)
EXPECTED_CENTERED_RREF_SHA256 = (
    "0cfa2e3330f92ab59fd183e5664715c490a702df2ad74491a8180793cae4a21e"
)
EXPECTED_STACK_RREF_SHA256 = (
    "9f4ec33d31337b0100a55871f6284443c6e5cfc0f5133a1493f807d563670821"
)
EXPECTED_SEMANTIC_SHA256 = (
    "861d4d95fd834e62cec842a5fc548e779554b17c79fabcf79cf7338e9db29848"
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


def load_module(path: Path, name: str, expected_sha256: str):
    observed = lf_sha256(path)
    require(observed == expected_sha256, (name, "script hash", observed, expected_sha256))
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, (name, "loader"))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


TWO = load_module(
    TWO_CURRENT_PATH, "thm3585_two_current_parent", EXPECTED_PARENT_SCRIPT_SHA256[0]
)
SOURCE = load_module(
    SOURCE_CURRENT_PATH,
    "thm3585_source_current_parent",
    EXPECTED_PARENT_SCRIPT_SHA256[1],
)

PRIME = TWO.PRIME
require(PRIME == SOURCE.PRIME == EXPECTED_PRIME, "parent prime mismatch")
require(TWO.P == SOURCE.P == P and TWO.V == SOURCE.V == V, "parent shape mismatch")


def two_worker(alpha: int):
    """Return only the fields needed by TWO.inverse_compressed."""
    return TWO.worker(alpha)[:2]


def source_worker(alpha: int):
    """Return only the fields needed by SOURCE.inverse_tensor."""
    return SOURCE.worker(alpha)[:2]


def rank_mod(rows) -> int:
    frozen = tuple(tuple(value % PRIME for value in row) for row in rows)
    return TWO.C.rank_mod(frozen)


def canonical_row_basis(rows):
    """Return the unique reduced-row-echelon row basis over F_p."""
    matrix = [list(value % PRIME for value in row) for row in rows]
    if not matrix:
        return ()
    columns = len(matrix[0])
    require(all(len(row) == columns for row in matrix), "ragged row matrix")
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
            if row == pivot_row or not matrix[row][column]:
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


def rowspace_digest(rows) -> str:
    return digest_json(canonical_row_basis(rows))


def flatten_two_current(tensor):
    return tuple(
        tuple(
            tensor[state][r0][r1][relation]
            for state in range(V)
            for relation in range(P)
        )
        for r0 in range(P)
        for r1 in range(P)
    )


def flatten_source_current(tensor):
    return tuple(
        tuple(
            tensor[channel][state][relation]
            for state in range(V)
            for relation in range(P)
        )
        for channel in range(len(tensor))
    )


def state_relation_center(rows):
    """Apply (I-J_4/4) tensor (I-J_13/13) to each row."""
    inverse_v = pow(V, -1, PRIME)
    inverse_p = pow(P, -1, PRIME)
    inverse_vp = pow(V * P, -1, PRIME)
    answer = []
    for row in rows:
        require(len(row) == V * P, ("row length", len(row), V * P))
        matrix = [list(row[state * P:(state + 1) * P]) for state in range(V)]
        relation_means = tuple(
            sum(matrix[state][relation] for state in range(V)) * inverse_v % PRIME
            for relation in range(P)
        )
        state_means = tuple(
            sum(matrix[state]) * inverse_p % PRIME for state in range(V)
        )
        grand_mean = sum(map(sum, matrix)) * inverse_vp % PRIME
        centered = tuple(
            (
                matrix[state][relation]
                - relation_means[relation]
                - state_means[state]
                + grand_mean
            )
            % PRIME
            for state in range(V)
            for relation in range(P)
        )
        require(
            all(
                sum(centered[state * P + relation] for state in range(V)) % PRIME == 0
                for relation in range(P)
            ),
            "state centering failed",
        )
        require(
            all(
                sum(centered[state * P + relation] for relation in range(P)) % PRIME == 0
                for state in range(V)
            ),
            "relation centering failed",
        )
        answer.append(centered)
    return tuple(answer)


def flatten_two_four_way(tensor):
    interaction = TWO.four_way_interaction(tensor)
    return flatten_two_current(interaction)


def reconstruct_two_current():
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(two_worker, range(P)))
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)), "two worker order")
    tensor, tau_core = TWO.inverse_compressed(chunks, TWO.C.context()["zeta"])
    tensor_digest = TWO.digest_json(tensor)
    tau_core_digest = TWO.digest_json(tau_core)
    require(
        tensor_digest == EXPECTED_PARENT_RECONSTRUCTION_SHA256[0]
        == TWO.EXPECTED_DIGESTS[0],
        ("two-current tensor drift", tensor_digest),
    )
    require(
        tau_core_digest == EXPECTED_PARENT_RECONSTRUCTION_SHA256[1]
        == TWO.EXPECTED_DIGESTS[1],
        ("two-current tau-core drift", tau_core_digest),
    )
    return tensor, (tensor_digest, tau_core_digest)


def reconstruct_source_current():
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(source_worker, range(P)))
    require(
        tuple(chunk[0] for chunk in chunks) == tuple(range(P)),
        "source worker order",
    )
    gamma = tuple(
        row for chunk in chunks for beta_rows in chunk[1] for row in beta_rows
    )
    gamma_digest = SOURCE.digest(gamma)
    tensor = SOURCE.inverse_tensor(gamma, SOURCE.B.context()["zeta"])
    tensor_digest = SOURCE.digest(tensor)
    require(
        gamma_digest == EXPECTED_PARENT_RECONSTRUCTION_SHA256[2],
        ("source-current gamma drift", gamma_digest),
    )
    require(
        tensor_digest == EXPECTED_PARENT_RECONSTRUCTION_SHA256[3],
        ("source-current tensor drift", tensor_digest),
    )
    return tensor, (gamma_digest, tensor_digest)


def main() -> None:
    dependency_hashes = (
        TWO.CURRENT_SHA256,
        SOURCE.SOURCE_SHA256,
        SOURCE.OWNER_SHA256,
    )
    require(
        dependency_hashes == EXPECTED_DEPENDENCY_SHA256,
        ("dependency parent hash drift", dependency_hashes),
    )

    two_tensor, two_reconstruction = reconstruct_two_current()
    source_tensor, source_reconstruction = reconstruct_source_current()
    reconstruction_hashes = two_reconstruction + source_reconstruction
    require(
        reconstruction_hashes == EXPECTED_PARENT_RECONSTRUCTION_SHA256,
        ("parent reconstruction drift", reconstruction_hashes),
    )

    two_raw = flatten_two_current(two_tensor)
    source_raw = flatten_source_current(source_tensor)
    two_centered = state_relation_center(two_raw)
    source_centered = state_relation_center(source_raw)
    two_four_way = flatten_two_four_way(two_tensor)

    raw_bases = (canonical_row_basis(source_raw), canonical_row_basis(two_raw))
    centered_bases = (
        canonical_row_basis(source_centered),
        canonical_row_basis(two_centered),
    )
    four_way_basis = canonical_row_basis(two_four_way)
    stack_basis = canonical_row_basis(source_raw + source_centered)

    raw_ranks = (
        len(raw_bases[0]),
        len(raw_bases[1]),
        rank_mod(source_raw + two_raw),
    )
    centered_ranks = (
        len(centered_bases[0]),
        len(centered_bases[1]),
        rank_mod(source_centered + two_centered),
    )
    complement_ranks = (
        rank_mod(source_raw + source_centered),
        rank_mod(two_raw + two_centered),
        rank_mod(source_raw + two_raw + source_centered + two_centered),
    )
    four_way_ranks = (
        len(four_way_basis),
        rank_mod(two_four_way + source_centered),
    )

    require(raw_bases[0] == raw_bases[1], "raw A4 row spaces differ")
    require(centered_bases[0] == centered_bases[1], "centered A4 row spaces differ")
    require(four_way_basis == centered_bases[0], "four-way image misses centered A4")
    require(raw_ranks == (4, 4, 4), ("raw ranks", raw_ranks))
    require(centered_ranks == (4, 4, 4), ("centered ranks", centered_ranks))
    require(complement_ranks == (8, 8, 8), ("complement ranks", complement_ranks))
    require(four_way_ranks == (4, 4), ("four-way ranks", four_way_ranks))

    raw_digest = digest_json(raw_bases[0])
    centered_digest = digest_json(centered_bases[0])
    stack_digest = digest_json(stack_basis)
    require(raw_digest == EXPECTED_RAW_RREF_SHA256, ("raw RREF drift", raw_digest))
    require(
        centered_digest == EXPECTED_CENTERED_RREF_SHA256,
        ("centered RREF drift", centered_digest),
    )
    require(
        stack_digest == EXPECTED_STACK_RREF_SHA256,
        ("stack RREF drift", stack_digest),
    )

    rank_record = (raw_ranks, centered_ranks, complement_ranks, four_way_ranks)
    rref_record = (raw_digest, centered_digest, stack_digest)
    semantic_surface = (
        EXPECTED_PRIME,
        EXPECTED_PARENT_SCRIPT_SHA256,
        dependency_hashes,
        reconstruction_hashes,
        rank_record,
        rref_record,
    )
    semantic = digest_json(semantic_surface)
    require(
        semantic == EXPECTED_SEMANTIC_SHA256,
        ("semantic drift", semantic, EXPECTED_SEMANTIC_SHA256),
    )

    print("== THM-3585 common A4 channel plane and centering complement ==")
    print(f"field=(prime={PRIME},V={V},relation_modulus={P},target_dimension={V*P})")
    print(f"parent_script_sha256_lf={EXPECTED_PARENT_SCRIPT_SHA256}")
    print(f"dependency_parent_sha256={dependency_hashes}")
    print(f"parent_reconstruction_sha256={reconstruction_hashes}")
    print(
        "raw=(source_rows=17,two_current_rows=169,"
        f"ranks_source_two_stack={raw_ranks},intersection=4,equal=True,"
        f"canonical_rref_sha256={raw_digest})"
    )
    print(
        "centered=(source_rows=17,two_current_rows=169,"
        f"ranks_source_two_stack={centered_ranks},intersection=4,equal=True,"
        f"canonical_rref_sha256={centered_digest})"
    )
    print(
        "raw_plus_centered=(source_stack,two_current_stack,all_stack)="
        f"{complement_ranks};intersection=0;canonical_rref_sha256={stack_digest}"
    )
    print(
        "two_current_four_way_matches_centered="
        f"(four_way_rank,stack_rank)={four_way_ranks}: PASS"
    )
    print(f"semantic_sha256={semantic}")
    print(f"script_sha256_lf={lf_sha256(Path(__file__).resolve())}")
    print("status=PROVED + FINITE-EXACT + VERIFIED-EXACT over the pinned prime")
    print(
        "scope=static finite-field identity only;not chronology;not current;"
        "not physical entry;not row exclusion;not characteristic-zero lift;"
        "not LRC(14)"
    )
    print(
        "commands=python -B 04-computation/"
        "lrc_common_a4_channel_plane_centering_complement_thm3585.py;"
        "python -B -O 04-computation/"
        "lrc_common_a4_channel_plane_centering_complement_thm3585.py"
    )
    print("PASS")


if __name__ == "__main__":
    main()
