#!/usr/bin/env python3
"""Pull actual U_full endpoint pair weights through the r=5 common gauge.

The source dependency constructs a lawful THM-2471 common-ancestry measure
M_(omega,nu,c) on the actual 39 U_full guard atoms.  The endpoint dependency
independently reconstructs the actual AX/BY atom tables and the q_H-q_q5
pair contribution E_(omega,nu) in the certified split field.

This probe performs two distinct tests.

1. Pull the scalar endpoint pair function E back along the atom-address map
   and integrate it against M, retaining every common-offset character.
2. Test whether any fixed 4x4 channel map and one common drift convolution
   maps the common-gauge K4 spectral curve projectively to the full or the
   d!=0-restricted endpoint curve.

A nonzero pullback is an address-weighted auxiliary contraction, not a
physical current.  A failed projective fit is a connection obstruction, not
a failure of the Boolean support realization.  No row exclusion or LRC(14)
claim is made.
"""

from __future__ import annotations

from contextlib import redirect_stdout
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import importlib.util
import io
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
P = 13
SOURCE_PATH = ROOT / "04-computation/lrc_r5_common_ancestry_guard_atom_root_drift_probe_20260816.py"
SOURCE_SHA256 = "83f1fa49ac4d02e21a1d76fed169d101715a6620342714ed05b9172ae967a730"
SOURCE_SEMANTIC = "3d8c88fb7b9762f41ef35c00d980b99fc435c8352baf5dddb9fe412d1baeace0"
TARGET_PATH = ROOT / "04-computation/lrc_ufull_owner_boundary_k4xf13_endpoint_factorization_independent_audit_20260816.py"
TARGET_SHA256 = "f89be10c65bb77270199f9399b155d5a2c82c0da121b3e8589fe3c1f7e9824fc"
TARGET_BUCKET_SHA256 = "553cfe7289b0556a19a8bcd1a0382dc1372545358feac62e5229adca315f8a26"
EXPECTED_SEMANTIC_SHA256 = "6d4f3389f50c98a103c5c6cb2daedea0ed2dd86e267aaa5ab0f9a4585433693d"

WALSH_SIGNS = (
    (1, 1, 1, 1),
    (1, 1, -1, -1),
    (1, -1, 1, -1),
    (1, -1, -1, 1),
)
ACTIVE_PAIRS = ((0, 0), (0, 2), (2, 0), (2, 2))
_TARGET_WORKER_MODULE = None


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    body = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return sha256(body).hexdigest()


def load_module(path: Path, name: str, expected_hash: str):
    require(lf_sha256(path) == expected_hash, (name, "source hash drift"))
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, (name, "module loader"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def target_worker(alpha: int):
    """Process-safe clean target worker with one module cache per child."""
    global _TARGET_WORKER_MODULE
    if _TARGET_WORKER_MODULE is None:
        _TARGET_WORKER_MODULE = load_module(
            TARGET_PATH, "ufull_endpoint_target_worker", TARGET_SHA256
        )
    return _TARGET_WORKER_MODULE.worker(alpha)


def rref_nullspace(matrix: list[list[int]], prime: int) -> tuple[int, tuple[tuple[int, ...], ...]]:
    rows = [[value % prime for value in row] for row in matrix]
    if not rows:
        return 0, ()
    columns = len(rows[0])
    pivot_columns = []
    rank = 0
    for column in range(columns):
        pivot = next((row for row in range(rank, len(rows)) if rows[row][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, prime)
        rows[rank] = [value * inverse % prime for value in rows[rank]]
        for row in range(len(rows)):
            if row == rank or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [
                (value - factor * pivot_value) % prime
                for value, pivot_value in zip(rows[row], rows[rank])
            ]
        pivot_columns.append(column)
        rank += 1
        if rank == len(rows):
            break
    free_columns = [column for column in range(columns) if column not in pivot_columns]
    basis = []
    for free in free_columns:
        vector = [0] * columns
        vector[free] = 1
        for row, pivot in enumerate(pivot_columns):
            vector[pivot] = -rows[row][free] % prime
        basis.append(tuple(vector))
    return rank, tuple(basis)


def rank_mod(rows: list[list[int]], prime: int) -> int:
    return rref_nullspace(rows, prime)[0]


def dft_rows(rows: list[list[int]], zeta: int, prime: int) -> list[list[int]]:
    return [
        [
            sum(value * pow(zeta, -frequency * drift % P, prime)
                for drift, value in enumerate(row)) % prime
            for frequency in range(P)
        ]
        for row in rows
    ]


def walsh(rows: list[list[int]], prime: int) -> list[list[int]]:
    return [
        [
            sum(sign * rows[index][drift] for index, sign in enumerate(signs)) % prime
            for drift in range(P)
        ]
        for signs in WALSH_SIGNS
    ]


def projective_system(
    source_spectra: list[list[int]],
    target_spectra: list[list[int]],
    prime: int,
) -> tuple[int, int, int]:
    """Return augmented rank, nullity, and multiplier-projection rank."""
    equations = []
    for frequency in range(P):
        source = [source_spectra[row][frequency] for row in range(4)]
        target = [target_spectra[row][frequency] for row in range(4)]
        for output in range(4):
            equation = [0] * (16 + P)
            for input_index in range(4):
                equation[4 * output + input_index] = source[input_index]
            equation[16 + frequency] = -target[output] % prime
            equations.append(equation)
    rank, basis = rref_nullspace(equations, prime)
    multiplier_projection = [list(vector[16:]) for vector in basis]
    projection_rank = rank_mod(multiplier_projection, prime) if multiplier_projection else 0
    return rank, 16 + P - rank, projection_rank


def main() -> None:
    source_module = load_module(SOURCE_PATH, "common_gauge_source", SOURCE_SHA256)
    with redirect_stdout(io.StringIO()):
        source = source_module.main()
    require(source["semantic_sha256"] == SOURCE_SEMANTIC, "source semantic drift")

    target = load_module(TARGET_PATH, "ufull_endpoint_target", TARGET_SHA256)
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(target_worker, range(P)))
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)), "target worker order")
    tables = tuple(row for chunk in chunks for row in chunk[1])
    require(len(tables) == P**2, "target table count")

    (
        _word,
        _t_den,
        _nn,
        prime,
        _root,
        zeta,
        _q_intervals,
        _q_starts,
        _embeddings,
        _tabs,
        _atom_intervals,
    ) = target.context()
    require(pow(zeta, P, prime) == 1 and zeta != 1, "target zeta order")

    # Reconstruct the two target guard kernels independently at pair level.
    h_kernel = [[0] * len(target.ATOMS) for _left in target.ATOMS]
    q5_kernel = [[0] * len(target.ATOMS) for _left in target.ATOMS]
    for left_index, (left_sheet, left_chamber) in enumerate(target.ATOMS):
        for right_index, (right_sheet, right_chamber) in enumerate(target.ATOMS):
            products = tuple(
                target.safe(left_chamber, left_sheet + tau)
                * target.safe(right_chamber, right_sheet + tau)
                for tau in range(P)
            )
            q5_kernel[left_index][right_index] = sum(products) % prime
            h_kernel[left_index][right_index] = sum(
                value * pow(zeta, -tau % P, prime)
                for tau, value in enumerate(products)
            ) % prime

    pair_weight = [[0] * len(target.ATOMS) for _left in target.ATOMS]
    table_index = 0
    for alpha in range(P):
        alpha_weight = pow(zeta, -alpha % P, prime)
        for beta in range(P):
            stored_beta, ax_values, by_values = tables[table_index]
            table_index += 1
            require(stored_beta == beta, (alpha, beta, stored_beta))
            phase = pow(zeta, beta, prime) * alpha_weight % prime
            for left_index, left_value in enumerate(ax_values):
                if not left_value:
                    continue
                for right_index, right_value in enumerate(by_values):
                    if not right_value:
                        continue
                    kernel = (h_kernel[left_index][right_index]
                              - q5_kernel[left_index][right_index]) % prime
                    pair_weight[left_index][right_index] = (
                        pair_weight[left_index][right_index]
                        + phase * left_value % prime * right_value % prime * kernel
                    ) % prime
    normalizer = pow(P**3, -1, prime)
    pair_weight = [[value * normalizer % prime for value in row] for row in pair_weight]
    bridge = sum(sum(row) for row in pair_weight) % prime
    require(bridge == target.EXPECTED_BRIDGE, ("target bridge", bridge))

    target_bucket = [[[[0] for _drift in range(P)] for _right in range(3)] for _left in range(3)]
    chamber_index = {name: index for index, name in enumerate(target.CHAMBER_NAMES)}
    flat_buckets = [0] * len(target.BUCKETS)
    for left_index, (left_sheet, left_chamber) in enumerate(target.ATOMS):
        for right_index, (right_sheet, right_chamber) in enumerate(target.ATOMS):
            drift = (right_sheet - left_sheet) % P
            value = pair_weight[left_index][right_index]
            target_bucket[chamber_index[left_chamber]][chamber_index[right_chamber]][drift][0] = (
                target_bucket[chamber_index[left_chamber]][chamber_index[right_chamber]][drift][0]
                + value
            ) % prime
            flat_buckets[target.BUCKET_INDEX[(left_chamber, right_chamber, drift)]] = (
                flat_buckets[target.BUCKET_INDEX[(left_chamber, right_chamber, drift)]] + value
            ) % prime
    require(
        target.digest_json(tuple(zip(target.BUCKETS, tuple(flat_buckets)))) == TARGET_BUCKET_SHA256,
        "target bridge-bucket digest",
    )

    target_full_rows = [
        [target_bucket[left][right][drift][0] for drift in range(P)]
        for left, right in ACTIVE_PAIRS
    ]
    target_restricted_rows = [
        [0 if drift == 0 else row[drift] for drift in range(P)]
        for row in target_full_rows
    ]
    same_sheet = sum(row[0] for row in target_full_rows) % prime
    restricted_bridge = sum(sum(row) for row in target_restricted_rows) % prime
    require(same_sheet == 324498447313453607031, ("same-sheet target", same_sheet))
    require(restricted_bridge == (bridge - same_sheet) % prime, "restricted bridge")
    require(restricted_bridge != 0, "48-bucket endpoint bridge vanished")

    gauge_offset = source["gauge_offset"]
    pair_gauge_offset = source["pair_gauge_offset"]
    source_denominator = source["denominator"]
    require(source_denominator % prime != 0, "source denominator bad at target prime")
    inverse_denominator = pow(source_denominator, -1, prime)

    source_rows_by_owner = []
    projective_full = []
    projective_restricted = []
    pullbacks = []
    target_full_spectra = dft_rows(target_full_rows, zeta, prime)
    target_restricted_spectra = dft_rows(target_restricted_rows, zeta, prime)
    for owner_frequency in range(P):
        source_rows = []
        for left, right in ACTIVE_PAIRS:
            source_rows.append([
                sum(
                    gauge_offset[left][right][drift][offset]
                    * pow(zeta, -owner_frequency * offset % P, prime)
                    for offset in range(P)
                ) % prime * inverse_denominator % prime
                for drift in range(P)
            ])
        source_rows_by_owner.append(source_rows)
        source_spectra = dft_rows(source_rows, zeta, prime)
        projective_full.append(projective_system(source_spectra, target_full_spectra, prime))
        projective_restricted.append(
            projective_system(source_spectra, target_restricted_spectra, prime)
        )

        pullback = 0
        for left_index in range(len(target.ATOMS)):
            for right_index in range(len(target.ATOMS)):
                source_weight = sum(
                    pair_gauge_offset[left_index][right_index][offset]
                    * pow(zeta, -owner_frequency * offset % P, prime)
                    for offset in range(P)
                ) % prime * inverse_denominator % prime
                pullback = (
                    pullback + source_weight * pair_weight[left_index][right_index]
                ) % prime
        pullbacks.append(pullback)

    source_ranks = tuple(rank_mod(rows, prime) for rows in source_rows_by_owner)
    require(source_ranks == (4,) * P, ("source owner ranks", source_ranks))
    require(all(row == (29, 0, 0) for row in projective_full),
            ("full-target projective systems", projective_full))
    require(all(row == (29, 0, 0) for row in projective_restricted),
            ("restricted-target projective systems", projective_restricted))
    pullback_support = tuple(int(value != 0) for value in pullbacks)

    target_full_walsh = walsh(target_full_rows, prime)
    target_restricted_walsh = walsh(target_restricted_rows, prime)
    target_full_fourier_support = tuple(
        sum(value != 0 for value in row)
        for row in dft_rows(target_full_walsh, zeta, prime)
    )
    target_restricted_fourier_support = tuple(
        sum(value != 0 for value in row)
        for row in dft_rows(target_restricted_walsh, zeta, prime)
    )
    require(target_full_fourier_support == (13, 13, 13, 13),
            ("full target spectrum", target_full_fourier_support))
    require(target_restricted_fourier_support == (13, 13, 13, 13),
            ("restricted target spectrum", target_restricted_fourier_support))

    flat_pullback = sum(
        sum(sum(offsets) for offsets in row)
        for row in pair_gauge_offset
    ) % prime * inverse_denominator % prime
    require(flat_pullback != 0, "flat address function lost common-gauge mass")

    record = (
        SOURCE_SHA256,
        SOURCE_SEMANTIC,
        TARGET_SHA256,
        TARGET_BUCKET_SHA256,
        prime,
        zeta,
        bridge,
        same_sheet,
        restricted_bridge,
        digest_json(pair_weight),
        digest_json(target_full_rows),
        digest_json(target_restricted_rows),
        target_full_fourier_support,
        target_restricted_fourier_support,
        source_ranks,
        tuple(projective_full),
        tuple(projective_restricted),
        tuple(pullbacks),
        pullback_support,
        flat_pullback,
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256, (semantic, EXPECTED_SEMANTIC_SHA256))

    print("== r=5 common-gauge / actual U_full endpoint-weight connection probe ==")
    print(f"dependencies=((source,{SOURCE_SHA256},{SOURCE_SEMANTIC}),(target,{TARGET_SHA256},{TARGET_BUCKET_SHA256}))")
    print(f"target_split_field=(prime={prime},zeta13={zeta})")
    print(f"target_pair_weight_sha256={digest_json(pair_weight)}")
    print(f"target_bridge=(full={bridge},same_sheet_d0={same_sheet},common_support_d_nonzero={restricted_bridge})")
    print("target_common_support=48/52; deleting d=0 retains nonzero bridge and full (13,13,13,13) Walsh drift spectrum")
    print(f"target_Walsh_F13_support=(full={target_full_fourier_support},restricted={target_restricted_fourier_support})")
    print(f"source_common_gauge_owner_ranks={source_ranks}")
    print(f"projective_connection_to_full_target={tuple(projective_full)}; each=(rank29,nullity0,lambda_projection0)")
    print(f"projective_connection_to_restricted_target={tuple(projective_restricted)}; each=(rank29,nullity0,lambda_projection0)")
    print(f"endpoint_pair_function_pullback_by_owner_frequency={tuple(pullbacks)}")
    print(f"endpoint_pair_function_pullback_support={sum(pullback_support)}/13 profile={pullback_support}")
    print(f"flat_pair_function_positive_control={flat_pullback}")
    print(f"semantic_sha256={semantic}")
    print("scope=actual endpoint weights pulled through a lawful address measure; auxiliary contraction and sharp connection no-go only; no physical current,row exclusion,LRC(14)")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
