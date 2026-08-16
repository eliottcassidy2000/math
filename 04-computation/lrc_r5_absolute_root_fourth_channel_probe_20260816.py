#!/usr/bin/env python3
"""Exact fourth-channel probe from THM-2594's fixed-root hostile.

The theta-slaved r=5 response has only three nonempty deep windows.  This
probe asks whether the already-lawful fixed-absolute-root reindexing of the
same common-base Boolean table supplies the missing fourth source direction
needed by the rank-four U_full Walsh bank.  The comparison remains a marked
split-field, endpoint-level linear test; it constructs no physical current
or ancestry relation for U_full.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stdout
import hashlib
import importlib.util
import io
import itertools
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PRIMARY_PATH = ROOT / "04-computation/lrc14_stage2_theta_contraction_opus_20260728.py"
TRANSPORTER_PATH = ROOT / "04-computation/lrc_r5_folded_c7_to_ufull_k4_drift_transporter_probe_20260816.py"
PRIMARY_SHA256 = "09c43af0a0a56c7a0833bbfd13ed6a96bc5a7a3718aa1bc6b77a144bde101a06"
TRANSPORTER_SHA256 = "5fe3f696f122869462ff73b1b4ebdb957fa8ca7ee3692c25ef94f0f7efae81cf"
EXPECTED_GAMMA_SHA256 = "1fabc5cfdbaa1455e10cd6bf9264488133616a7b0ff381623d729b4b4bfa9682"
EXPECTED_SEMANTIC_SHA256 = "6ad55eeab7f8e5479194965e819b030f177f4500be37cb189688262dc8af35da"
P = 13


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest(value: object) -> str:
    body = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return hashlib.sha256(body).hexdigest()


def load_module(path: Path, name: str, expected_hash: str):
    actual = lf_sha256(path)
    require(actual == expected_hash, (name, actual, expected_hash))
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, (name, "loader"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


T = load_module(TRANSPORTER_PATH, "lrc_r5_folded_transporter_parent", TRANSPORTER_SHA256)


def ufull_worker(alpha: int) -> tuple[object, ...]:
    return T.worker(alpha)


def primary_joint_worker() -> tuple[tuple[tuple[tuple[tuple[int, ...], ...], ...], ...], int]:
    primary = load_module(PRIMARY_PATH, "lrc_r5_primary_parent", PRIMARY_SHA256)
    buffer = io.StringIO()
    with redirect_stdout(buffer):
        state = primary.main()
    # The primary return packet places the four-dimensional N(u,q,ell,theta)
    # table at slot 12; slots 13 and 14 are its word-weighted companions.
    joint = state[12]
    denominator = state[17]
    return tuple(
        tuple(
            tuple(tuple(int(value) for value in row) for row in ell_rows)
            for ell_rows in q_rows
        )
        for q_rows in joint
    ), int(denominator)


def rank_mod(matrix: tuple[tuple[int, ...], ...], prime: int) -> int:
    rows = [list(row) for row in matrix]
    rank = 0
    for column in range(len(rows[0]) if rows else 0):
        pivot = next(
            (row for row in range(rank, len(rows)) if rows[row][column] % prime),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column] % prime, -1, prime)
        rows[rank] = [entry * inverse % prime for entry in rows[rank]]
        for row in range(len(rows)):
            if row == rank:
                continue
            factor = rows[row][column] % prime
            if factor:
                rows[row] = [
                    (entry - factor * pivot_entry) % prime
                    for entry, pivot_entry in zip(rows[row], rows[rank])
                ]
        rank += 1
    return rank


def projective_system_rank(
    source_rows: tuple[tuple[int, ...], ...],
    target_rows: tuple[tuple[int, ...], ...],
    prime: int,
) -> tuple[int, int, int, int]:
    """Return source rank, equation rank, nullity, and excess nullity.

    For one frequency-independent channel map M and one common circulant C,
    the Fourier vectors must obey Y_a proportional to M X_a.  Eliminating the
    thirteen scalar multipliers gives a homogeneous linear system in the
    4*m entries of M.  Maps annihilating the full source span are automatic
    solutions; excess nullity counts candidate transports beyond them.
    """
    channels = len(source_rows)
    require(channels and all(len(row) == P for row in source_rows), channels)
    require(len(target_rows) == 4 and all(len(row) == P for row in target_rows), "target")
    equations = []
    for frequency in range(P):
        source = tuple(row[frequency] for row in source_rows)
        target = tuple(row[frequency] for row in target_rows)
        require(any(target), ("zero target frequency", frequency))
        for left in range(4):
            for right in range(left + 1, 4):
                equation = [0] * (4 * channels)
                for column in range(channels):
                    equation[right * channels + column] = (
                        equation[right * channels + column]
                        + target[left] * source[column]
                    ) % prime
                    equation[left * channels + column] = (
                        equation[left * channels + column]
                        - target[right] * source[column]
                    ) % prime
                equations.append(tuple(equation))
    source_rank = rank_mod(source_rows, prime)
    equation_rank = rank_mod(tuple(equations), prime)
    nullity = 4 * channels - equation_rank
    annihilator_dimension = 4 * (channels - source_rank)
    require(nullity >= annihilator_dimension, (nullity, annihilator_dimension))
    return source_rank, equation_rank, nullity, nullity - annihilator_dimension


def septimal_profiles(
    table: tuple[tuple[int, ...], ...], zeta7: int, prime: int
) -> tuple[tuple[int, ...], ...]:
    require(len(table) == 7 and all(len(row) == P for row in table), "7x13 table")
    return tuple(
        tuple(
            sum(
                table[ell][root] * pow(zeta7, -septimal * ell % 7, prime)
                for ell in range(7)
            )
            % prime
            for root in range(P)
        )
        for septimal in range(7)
    )


def dilate_spectra(
    rows: tuple[tuple[int, ...], ...], multiplier: int
) -> tuple[tuple[int, ...], ...]:
    """Apply a marked F_13-torsor generator change in Fourier coordinates.

    Translation contributes one common phase and is already absorbed by the
    common circulant.  The remaining affine gauge is multiplication by any
    nonzero field element, which permutes the twelve nontrivial frequencies.
    """
    require(1 <= multiplier < P, multiplier)
    return tuple(
        tuple(row[multiplier * frequency % P] for frequency in range(P))
        for row in rows
    )


def main() -> None:
    with ProcessPoolExecutor(max_workers=5) as pool:
        joint_future = pool.submit(primary_joint_worker)
        chunk_futures = tuple(pool.submit(ufull_worker, alpha) for alpha in range(P))
        joint, denominator = joint_future.result()
        chunks = tuple(future.result() for future in chunk_futures)

    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)), "worker order")
    gamma = tuple(value for chunk in chunks for value in chunk[5])
    require(T.U.digest_integers(gamma) == EXPECTED_GAMMA_SHA256, "gamma digest")

    (_word, _t_den, nn, prime, root, zeta13, *_rest) = T.U.B.context()
    require(nn % 7 == 0, nn)
    zeta7 = pow(root, nn // 7, prime)
    require(pow(zeta7, 7, prime) == 1 and zeta7 != 1, "zeta7 order")
    require(pow(zeta13, P, prime) == 1 and zeta13 != 1, "zeta13 order")
    require(denominator % prime, "bad denominator")

    normalizer = pow(P**3, -1, prime)
    q_buckets = {
        q: tuple(
            sum(chunk[6][q_index][bucket] for chunk in chunks)
            % prime
            * normalizer
            % prime
            for bucket in range(len(T.U.B.BUCKETS))
        )
        for q_index, q in enumerate(T.U.Q_CLASSES)
    }
    q_h = (1, 0, 1)
    q_q5 = (1, 0, 0)
    corners = tuple(
        tuple(
            (
                q_buckets[q_h][T.U.B.BUCKET_INDEX[(left, right, drift)]]
                - q_buckets[q_q5][T.U.B.BUCKET_INDEX[(left, right, drift)]]
            )
            % prime
            for drift in range(P)
        )
        for left, right in T.U.CORNER_PAIRS
    )
    walsh = (
        tuple(sum(corners[row][d] for row in range(4)) % prime for d in range(P)),
        tuple((corners[0][d] + corners[1][d] - corners[2][d] - corners[3][d]) % prime for d in range(P)),
        tuple((corners[0][d] - corners[1][d] + corners[2][d] - corners[3][d]) % prime for d in range(P)),
        tuple((corners[0][d] - corners[1][d] - corners[2][d] + corners[3][d]) % prime for d in range(P)),
    )
    target_spectra = tuple(T.spectrum(row, zeta13, prime) for row in walsh)
    require(rank_mod(walsh, prime) == 4, "target rank")
    require(all(all(row) for row in target_spectra), "target spectral support")

    slaved = [[0] * P for _ in range(7)]
    absolute = [[0] * P for _ in range(7)]
    for u in range(P):
        for q in range(P):
            for ell in range(7):
                for theta in range(3):
                    value = joint[u][q][ell][theta] % prime
                    slaved[ell][theta] = (slaved[ell][theta] + value) % prime
                    absolute_root = (theta + 2 * u) % P
                    absolute[ell][absolute_root] = (
                        absolute[ell][absolute_root] + value
                    ) % prime
    slaved_table = tuple(tuple(row) for row in slaved)
    absolute_table = tuple(tuple(row) for row in absolute)
    require(sum(value != 0 for row in absolute_table for value in row) == 18, "absolute raw support")

    slaved_profiles = septimal_profiles(slaved_table, zeta7, prime)
    absolute_profiles = septimal_profiles(absolute_table, zeta7, prime)
    slaved_spectra = tuple(T.spectrum(row, zeta13, prime) for row in slaved_profiles)
    absolute_spectra = tuple(T.spectrum(row, zeta13, prime) for row in absolute_profiles)
    slaved_rank = rank_mod(slaved_profiles, prime)
    absolute_rank = rank_mod(absolute_profiles, prime)
    union_rank = rank_mod(slaved_profiles + absolute_profiles, prime)
    slaved_support = sum(bool(value) for row in slaved_spectra for value in row)
    absolute_support = sum(bool(value) for row in absolute_spectra for value in row)
    require(slaved_rank == 3 and slaved_support == 91, (slaved_rank, slaved_support))

    whole_systems = {
        "slaved7": projective_system_rank(slaved_spectra, target_spectra, prime),
        "absolute7": projective_system_rank(absolute_spectra, target_spectra, prime),
        "union14": projective_system_rank(slaved_spectra + absolute_spectra, target_spectra, prime),
    }
    affine_torsor_systems = tuple(
        (
            multiplier,
            projective_system_rank(
                dilate_spectra(slaved_spectra + absolute_spectra, multiplier),
                target_spectra,
                prime,
            ),
        )
        for multiplier in range(1, P)
    )

    folded_pairs = ((1, 6), (2, 5), (3, 4))
    one_named_sidecar = []
    for absolute_mode in range(7):
        for signs in itertools.product((0, 1), repeat=3):
            selected = tuple(folded_pairs[index][signs[index]] for index in range(3))
            rows = tuple(slaved_spectra[index] for index in selected) + (
                absolute_spectra[absolute_mode],
            )
            one_named_sidecar.append(
                (absolute_mode, signs, projective_system_rank(rows, target_spectra, prime))
            )
    require(len(one_named_sidecar) == 56, len(one_named_sidecar))

    binary_sidecars = []
    for mask in range(1, 1 << 7):
        sidecar = tuple(
            sum(
                absolute_spectra[index][frequency]
                for index in range(7)
                if mask & (1 << index)
            )
            % prime
            for frequency in range(P)
        )
        rows = tuple(slaved_spectra[index] for index in (1, 2, 3)) + (sidecar,)
        binary_sidecars.append((mask, projective_system_rank(rows, target_spectra, prime)))
    require(len(binary_sidecars) == 127, len(binary_sidecars))

    named_rank4 = sum(record[2][0] == 4 for record in one_named_sidecar)
    named_positive = tuple(record for record in one_named_sidecar if record[2][3] > 0)
    binary_rank4 = sum(record[1][0] == 4 for record in binary_sidecars)
    binary_positive = tuple(record for record in binary_sidecars if record[1][3] > 0)
    binary_rank3_masks = tuple(record[0] for record in binary_sidecars if record[1][0] == 3)
    affine_torsor_positive = tuple(
        record for record in affine_torsor_systems if record[1][3] > 0
    )

    semantic = (
        PRIMARY_SHA256,
        TRANSPORTER_SHA256,
        (prime, zeta7, zeta13),
        denominator,
        slaved_table,
        absolute_table,
        slaved_profiles,
        absolute_profiles,
        target_spectra,
        (slaved_rank, absolute_rank, union_rank, slaved_support, absolute_support),
        whole_systems,
        affine_torsor_systems,
        tuple(one_named_sidecar),
        tuple(binary_sidecars),
        "fixed-root reindex is common-base but endpoint transplant remains formal",
    )
    semantic_hash = hashlib.sha256(repr(semantic).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, (semantic_hash, EXPECTED_SEMANTIC_SHA256))

    print("R5 FIXED-ABSOLUTE-ROOT FOURTH-CHANNEL PROBE")
    print("status=FINITE-EXACT common-base source reindex versus endpoint Walsh target; LRC(14) OPEN")
    print(f"dependencies=(primary={PRIMARY_SHA256},transporter={TRANSPORTER_SHA256})")
    print(f"field=(prime={prime},zeta7={zeta7},zeta13={zeta13})")
    print(f"source_cell_support=(slaved={sum(value != 0 for row in slaved_table for value in row)}/91,absolute={sum(value != 0 for row in absolute_table for value in row)}/91)")
    print(f"source_spectral_support=(slaved={slaved_support}/91,absolute={absolute_support}/91)")
    print(f"source_row_ranks=(slaved={slaved_rank},absolute={absolute_rank},union={union_rank})")
    print(f"source_plane_intersection_dimension={slaved_rank + absolute_rank - union_rank}")
    print(f"whole_projective_systems={whole_systems}")
    print(f"affine_torsor_dilation_systems={affine_torsor_systems}")
    print(f"affine_torsor_positive_cases={len(affine_torsor_positive)}")
    print("projective_tuple_fields=(source_rank,equation_rank,nullity,excess_beyond_source_annihilators)")
    print(f"one_named_absolute_sidecar=(cases=56,rank4={named_rank4},positive_transport_cases={len(named_positive)})")
    print(f"one_named_positive_records={named_positive}")
    print(f"binary_absolute_sidecar=(cases=127,rank4={binary_rank4},positive_transport_cases={len(binary_positive)})")
    print(f"binary_rank3_masks={binary_rank3_masks}")
    print(f"binary_positive_records={binary_positive}")
    print(f"slaved_table_sha256={digest(slaved_table)}")
    print(f"absolute_table_sha256={digest(absolute_table)}")
    print(f"semantic_sha256={semantic_hash}")
    print("nonconsequence=no U_full ancestry relation, physical current, H1 flux, scalar-row exclusion, or LRC(14)")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
