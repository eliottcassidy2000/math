#!/usr/bin/env python3
"""Exact root-difference versus U_full drift alignment probe.

THM-3514 identifies endpoint drift b-a with the difference of two collision
root labels once one common F_13-torsor gauge is chosen.  The earlier folded
transporter instead used the deep coordinate theta as drift.  This probe
returns to THM-2594's joint common-base table and retains the typed root
difference s=u-q before marginalization.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENT_PATH = ROOT / "04-computation/lrc_r5_absolute_root_fourth_channel_probe_20260816.py"
PARENT_SHA256 = "113793a41d39e8ac1ad3d745b9c9238cba92721805a379ee4378a000bacca33c"
EXPECTED_GAMMA_SHA256 = "1fabc5cfdbaa1455e10cd6bf9264488133616a7b0ff381623d729b4b4bfa9682"
EXPECTED_SEMANTIC_SHA256 = "a954aa6ac96fa0a4e7b77e84473971c8444393409d61ca251b2c31de1b27779f"
P = 13


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load_module(path: Path, name: str, expected_hash: str):
    actual = lf_sha256(path)
    require(actual == expected_hash, (name, actual, expected_hash))
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, (name, "loader"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


A = load_module(PARENT_PATH, "lrc_r5_absolute_root_parent", PARENT_SHA256)
T = A.T


def ufull_worker(alpha: int) -> tuple[object, ...]:
    return A.ufull_worker(alpha)


def primary_joint_worker():
    return A.primary_joint_worker()


def digest(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def reconstruct_target(chunks: tuple[tuple[object, ...], ...]):
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)), "worker order")
    gamma = tuple(value for chunk in chunks for value in chunk[5])
    require(T.U.digest_integers(gamma) == EXPECTED_GAMMA_SHA256, "gamma digest")
    (_word, _t_den, nn, prime, root, zeta13, *_rest) = T.U.B.context()
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
    require(A.rank_mod(walsh, prime) == 4, "target rank")
    require(all(all(row) for row in target_spectra), "target spectral support")
    require(nn % 7 == 0, nn)
    zeta7 = pow(root, nn // 7, prime)
    return prime, zeta7, zeta13, walsh, target_spectra


def ratio_equivalent(
    left: tuple[int, ...], right: tuple[int, ...], zeta13: int, prime: int
) -> bool:
    require(all(left) and all(right), "unit ratios required")
    scalar = left[0] * pow(right[0], -1, prime) % prime
    for shift in range(P):
        if all(
            left[frequency]
            == scalar
            * pow(zeta13, frequency * shift % P, prime)
            % prime
            * right[frequency]
            % prime
            for frequency in range(P)
        ):
            return True
    return False


def class_count(
    rows: tuple[tuple[int, ...], ...], zeta13: int, prime: int
) -> int:
    representatives = []
    for row in rows:
        if not any(ratio_equivalent(row, representative, zeta13, prime) for representative in representatives):
            representatives.append(row)
    return len(representatives)


def main() -> None:
    with ProcessPoolExecutor(max_workers=5) as pool:
        joint_future = pool.submit(primary_joint_worker)
        chunk_futures = tuple(pool.submit(ufull_worker, alpha) for alpha in range(P))
        joint, denominator = joint_future.result()
        chunks = tuple(future.result() for future in chunk_futures)

    prime, zeta7, zeta13, walsh, target_spectra = reconstruct_target(chunks)
    require(denominator % prime, "bad denominator")

    # q=u-s is the source root paired with current root u.  Summation over
    # theta retains every deep window while preserving the root difference.
    shift = [[0] * P for _ in range(7)]
    shift_theta = [[[0] * P for _ in range(3)] for _ in range(7)]
    owner_slices = [[[0] * P for _ in range(7)] for _ in range(2)]
    for u in range(P):
        for root_shift in range(P):
            q = (u - root_shift) % P
            for ell in range(7):
                for theta in range(3):
                    value = joint[u][q][ell][theta] % prime
                    shift[ell][root_shift] = (shift[ell][root_shift] + value) % prime
                    shift_theta[ell][theta][root_shift] = (
                        shift_theta[ell][theta][root_shift] + value
                    ) % prime
                    for owner_frequency in range(2):
                        owner_slices[owner_frequency][ell][root_shift] = (
                            owner_slices[owner_frequency][ell][root_shift]
                            + value * pow(zeta13, -owner_frequency * u % P, prime)
                        ) % prime
    shift_table = tuple(tuple(row) for row in shift)
    require(all(shift_table[ell][0] == 0 for ell in range(7)), "same-root slice survives")

    source_profiles = A.septimal_profiles(shift_table, zeta7, prime)
    source_spectra = tuple(T.spectrum(row, zeta13, prime) for row in source_profiles)
    source_rank = A.rank_mod(source_profiles, prime)
    source_spectral_support = sum(bool(value) for row in source_spectra for value in row)
    unit_modes = tuple(index for index, row in enumerate(source_spectra) if all(row))

    # The theta-resolved bank tests whether the three deep windows add
    # channel information once the lawful drift coordinate is s.
    resolved_rows = tuple(
        tuple(
            sum(
                shift_theta[ell][theta][root_shift]
                * pow(zeta7, -septimal * ell % 7, prime)
                for ell in range(7)
            )
            % prime
            for root_shift in range(P)
        )
        for theta in range(3)
        for septimal in range(7)
    )
    resolved_spectra = tuple(T.spectrum(row, zeta13, prime) for row in resolved_rows)
    resolved_rank = A.rank_mod(resolved_rows, prime)

    owner_profiles = tuple(
        A.septimal_profiles(
            tuple(tuple(row) for row in owner_slices[owner_frequency]),
            zeta7,
            prime,
        )
        for owner_frequency in range(2)
    )
    owner_spectra = tuple(
        tuple(T.spectrum(row, zeta13, prime) for row in profile_bank)
        for profile_bank in owner_profiles
    )
    require(owner_profiles[0] == source_profiles, "owner-zero slice mismatch")
    owner_difference_profiles = tuple(
        tuple(
            (owner_profiles[1][septimal][root_shift] - owner_profiles[0][septimal][root_shift])
            % prime
            for root_shift in range(P)
        )
        for septimal in range(7)
    )
    owner_difference_spectra = tuple(
        T.spectrum(row, zeta13, prime) for row in owner_difference_profiles
    )
    owner_union_spectra = owner_spectra[0] + owner_spectra[1]
    owner_ranks = (
        A.rank_mod(owner_profiles[0], prime),
        A.rank_mod(owner_profiles[1], prime),
        A.rank_mod(owner_difference_profiles, prime),
        A.rank_mod(owner_profiles[0] + owner_profiles[1], prime),
    )

    whole_system = A.projective_system_rank(source_spectra, target_spectra, prime)
    resolved_system = A.projective_system_rank(resolved_spectra, target_spectra, prime)
    owner_systems = {
        "k0": A.projective_system_rank(owner_spectra[0], target_spectra, prime),
        "k1": A.projective_system_rank(owner_spectra[1], target_spectra, prime),
        "k1_minus_k0": A.projective_system_rank(
            owner_difference_spectra, target_spectra, prime
        ),
        "k0_union_k1": A.projective_system_rank(
            owner_union_spectra, target_spectra, prime
        ),
    }
    affine_torsor_systems = tuple(
        (
            multiplier,
            A.projective_system_rank(
                A.dilate_spectra(source_spectra, multiplier),
                target_spectra,
                prime,
            ),
        )
        for multiplier in range(1, P)
    )
    resolved_affine_systems = tuple(
        (
            multiplier,
            A.projective_system_rank(
                A.dilate_spectra(resolved_spectra, multiplier),
                target_spectra,
                prime,
            ),
        )
        for multiplier in range(1, P)
    )
    owner_union_affine_systems = tuple(
        (
            multiplier,
            A.projective_system_rank(
                A.dilate_spectra(owner_union_spectra, multiplier),
                target_spectra,
                prime,
            ),
        )
        for multiplier in range(1, P)
    )

    folded_pairs = ((1, 6), (2, 5), (3, 4))
    selected_systems = []
    for signs in itertools.product((0, 1), repeat=3):
        selected = (0,) + tuple(folded_pairs[index][signs[index]] for index in range(3))
        rows = tuple(source_spectra[index] for index in selected)
        selected_systems.append((signs, selected, A.projective_system_rank(rows, target_spectra, prime)))

    names = range(4)
    allocation_census = None
    if len(unit_modes) == 7:
        common_exact = 0
        common_gauge = 0
        observed_classes = []
        for selected in itertools.product(range(7), repeat=4):
            ratios = tuple(
                tuple(
                    target_spectra[name][frequency]
                    * pow(source_spectra[selected[name]][frequency], -1, prime)
                    % prime
                    for frequency in range(P)
                )
                for name in names
            )
            if len(set(ratios)) == 1:
                common_exact += 1
            classes = class_count(ratios, zeta13, prime)
            observed_classes.append(classes)
            if classes == 1:
                common_gauge += 1
        allocation_census = (
            7**4,
            common_exact,
            common_gauge,
            tuple(sorted(set(observed_classes))),
        )

    semantic = (
        PARENT_SHA256,
        (prime, zeta7, zeta13),
        denominator,
        shift_table,
        tuple(tuple(tuple(row) for row in ell_rows) for ell_rows in shift_theta),
        source_profiles,
        source_spectra,
        target_spectra,
        (source_rank, source_spectral_support, unit_modes),
        (resolved_rank, resolved_system),
        owner_profiles,
        owner_spectra,
        owner_difference_profiles,
        owner_difference_spectra,
        owner_ranks,
        owner_systems,
        whole_system,
        affine_torsor_systems,
        resolved_affine_systems,
        owner_union_affine_systems,
        tuple(selected_systems),
        allocation_census,
        "root difference is the equivariant drift candidate; no ancestry conclusion",
    )
    semantic_hash = hashlib.sha256(repr(semantic).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, (semantic_hash, EXPECTED_SEMANTIC_SHA256))

    print("R5 ROOT-DIFFERENCE TO U_FULL DRIFT ALIGNMENT PROBE")
    print("status=FINITE-EXACT typed-coordinate comparison; LRC(14) OPEN")
    print(f"dependency_parent={PARENT_SHA256}")
    print(f"field=(prime={prime},zeta7={zeta7},zeta13={zeta13})")
    print("typed_alignment=one common regular-torsor gauge gives endpoint d=b-a equal to a collision-root difference; source convention s=u-q differs only by an allowed sign")
    print(f"shift_table_support={sum(value != 0 for row in shift_table for value in row)}/91; same_root_column_zero=PASS")
    print(f"source_spectral_support={source_spectral_support}/91; unit_septimal_modes={unit_modes}")
    print(f"source_row_rank={source_rank}; theta_resolved_row_rank={resolved_rank}")
    print(f"owner_phase_row_ranks=(k0={owner_ranks[0]},k1={owner_ranks[1]},k1_minus_k0={owner_ranks[2]},union={owner_ranks[3]})")
    print(f"whole_projective_system={whole_system}")
    print(f"affine_torsor_systems={affine_torsor_systems}")
    print(f"theta_resolved_projective_system={resolved_system}")
    print(f"theta_resolved_affine_systems={resolved_affine_systems}")
    print(f"owner_phase_projective_systems={owner_systems}")
    print(f"owner_phase_union_affine_systems={owner_union_affine_systems}")
    print(f"folded_selected_systems={tuple(selected_systems)}")
    print(f"all_mode_allocation_census={allocation_census}")
    print("projective_tuple_fields=(source_rank,equation_rank,nullity,excess_beyond_source_annihilators)")
    print(f"shift_table_sha256={digest(shift_table)}")
    print(f"semantic_sha256={semantic_hash}")
    print("nonconsequence=abstract root/sheet torsor alignment only; no chamber support map, common U_full ancestry, physical current, H1 flux, row exclusion, or LRC(14)")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
