#!/usr/bin/env python3
"""Exact off-diagonal drift probe for the THM-2594 -> U_full transplant.

The source and target are reconstructed from the independent common-base
connection audit.  The probe asks one narrow question before any ancestry
claim: what happens if the Cartesian U_full endpoint tensor is forced to
obey THM-2471's same-root disjointness by deleting its drift-zero fibre?

This is a representation and support-typing experiment.  The deletion is a
necessary condition for a drift-preserving THM-2471 transplant, not a proved
U_full support relation.  No ancestry, current, H1, bispectrum, scalar-row
exclusion, or LRC(14) consequence is asserted.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
import hashlib
import importlib.util
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "04-computation/lrc_r5_offdiagonal_drift_branch_transplant_probe_20260816.py"
OUTPUT = "05-knowledge/results/lrc_r5_offdiagonal_drift_branch_transplant_probe_20260816.out"
COMMON_PATH = ROOT / (
    "04-computation/"
    "lrc_r5_common_base_connection_obstruction_independent_audit_20260816.py"
)
COMMON_SHA256 = "efbfb738f1901946210e1f438212d24e7cc34ee4ca3f02d436923d61a852bf43"
EXPECTED_SEMANTIC_SHA256: str | None = (
    "84a21447326ac6e86ff7743704653d1c280eee265cfdbe357735f5ad6d80915a"
)

N = 13


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest(value: object) -> str:
    return hashlib.sha256(repr(value).encode("utf-8")).hexdigest()


def load_common():
    require(lf_sha256(COMMON_PATH) == COMMON_SHA256, "common audit source drift")
    name = "common_base_obstruction_for_offdiagonal_probe"
    spec = importlib.util.spec_from_file_location(name, COMMON_PATH)
    require(spec is not None and spec.loader is not None, "common audit loader")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


C = load_common()


def joint_worker() -> tuple[object, ...]:
    return C.joint_worker()


def atom_worker(alpha: int) -> tuple[object, ...]:
    return C.atom_worker(alpha)


def invert_walsh(
    walsh: tuple[tuple[int, ...], ...], prime: int
) -> tuple[tuple[int, ...], ...]:
    require(len(walsh) == 4 and all(len(row) == N for row in walsh), "Walsh shape")
    inv4 = pow(4, -1, prime)
    signs = (
        (1, 1, 1, 1),
        (1, 1, -1, -1),
        (1, -1, 1, -1),
        (1, -1, -1, 1),
    )
    return tuple(
        tuple(
            sum(signs[corner][channel] * walsh[channel][drift] for channel in range(4))
            * inv4
            % prime
            for drift in range(N)
        )
        for corner in range(4)
    )


def delete_drift(
    rows: tuple[tuple[int, ...], ...], drift: int
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(0 if index == drift else value for index, value in enumerate(row))
        for row in rows
    )


def generalized_connection_certificate(
    source_rows: tuple[tuple[int, ...], ...],
    target_rows: tuple[tuple[int, ...], ...],
    prime: int,
) -> tuple[int, ...]:
    """Augmented/wedge certificate for an arbitrary target dimension."""
    channels = len(source_rows)
    outputs = len(target_rows)
    require(channels and outputs, "empty generalized connection")
    require(
        all(len(row) == N for row in source_rows + target_rows),
        "generalized connection row shape",
    )
    equations = []
    wedges = []
    for frequency in range(N):
        source = tuple(row[frequency] for row in source_rows)
        target = tuple(row[frequency] for row in target_rows)
        for output in range(outputs):
            equation = [0] * (outputs * channels + N)
            for column in range(channels):
                equation[output * channels + column] = source[column]
            equation[outputs * channels + frequency] = -target[output] % prime
            equations.append(tuple(equation))
        for left in range(outputs):
            for right in range(left + 1, outputs):
                equation = [0] * (outputs * channels)
                for column in range(channels):
                    equation[right * channels + column] += target[left] * source[column]
                    equation[left * channels + column] -= target[right] * source[column]
                wedges.append(tuple(value % prime for value in equation))

    augmented = tuple(equations)
    basis = C.nullspace(augmented, prime)
    augmented_rank = C.rank_mod(augmented, prime)
    augmented_nullity = outputs * channels + N - augmented_rank
    require(len(basis) == augmented_nullity, "generalized nullity")
    multiplier_rank = C.rank_mod(
        tuple(vector[outputs * channels :] for vector in basis), prime
    )
    source_rank = C.rank_mod(source_rows, prime)
    target_rank = C.rank_mod(target_rows, prime)
    wedge_rank = C.rank_mod(tuple(wedges), prime)
    wedge_nullity = outputs * channels - wedge_rank
    annihilator_dimension = outputs * (channels - source_rank)
    wedge_excess = wedge_nullity - annihilator_dimension
    nonannihilator = 0
    for vector in basis:
        for output in range(outputs):
            for frequency in range(N):
                value = sum(
                    vector[output * channels + column] * source_rows[column][frequency]
                    for column in range(channels)
                ) % prime
                nonannihilator += value != 0
    require((wedge_excess > 0) == (multiplier_rank > 0), (
        wedge_excess,
        multiplier_rank,
    ))
    return (
        source_rank,
        target_rank,
        wedge_rank,
        wedge_nullity,
        wedge_excess,
        augmented_rank,
        augmented_nullity,
        multiplier_rank,
        nonannihilator,
    )


def main() -> None:
    with ProcessPoolExecutor(max_workers=5) as pool:
        joint_future = pool.submit(joint_worker)
        atom_futures = tuple(pool.submit(atom_worker, alpha) for alpha in range(N))
        joint, denominator = joint_future.result()
        chunks = tuple(future.result() for future in atom_futures)

    nn, prime, root, zeta13, walsh, full_target_spectra = C.materialize_target(chunks)
    require(denominator % prime != 0 and nn % 7 == 0, "field compatibility")
    zeta7 = pow(root, nn // 7, prime)
    require(pow(zeta7, 7, prime) == 1 and zeta7 != 1, "zeta7")

    pointwise_same_root = tuple(
        joint[u][u][ell][theta]
        for u in range(N)
        for ell in range(7)
        for theta in range(3)
    )
    require(len(pointwise_same_root) == 273, len(pointwise_same_root))
    require(all(value == 0 for value in pointwise_same_root), "same-root source leak")

    difference = [[0] * N for _ in range(7)]
    for u in range(N):
        for source_root in range(N):
            drift = (u - source_root) % N
            for ell in range(7):
                for theta in range(3):
                    difference[ell][drift] = (
                        difference[ell][drift] + joint[u][source_root][ell][theta]
                    ) % prime
    difference_table = tuple(tuple(row) for row in difference)
    source_row_support = tuple(sum(value != 0 for value in row) for row in difference_table)
    require(source_row_support == (0, 12, 12, 12, 12, 12, 12), source_row_support)
    require(all(row[0] == 0 for row in difference_table), "source diagonal marginal")
    source_reflection = all(
        difference_table[ell][drift]
        == difference_table[(-ell) % 7][(-drift) % N]
        for ell in range(7)
        for drift in range(N)
    )
    require(source_reflection, "source diagonal reflection")

    source_profiles = C.septimal_profiles(difference_table, zeta7, prime)
    source_spectra = tuple(C.spectrum(row, zeta13, prime) for row in source_profiles)
    source_rank = C.rank_mod(source_profiles, prime)
    require(source_rank == 6, source_rank)
    require(all(all(value != 0 for value in row) for row in source_spectra), "source spectrum")

    corners = invert_walsh(walsh, prime)
    require(all(all(value != 0 for value in row) for row in corners), "full corner support")
    diagonal_corners = tuple(row[0] for row in corners)
    same_sheet_bridge = sum(diagonal_corners) % prime
    require(all(diagonal_corners), diagonal_corners)
    require(same_sheet_bridge == 324498447313453607031, same_sheet_bridge)

    off_walsh = delete_drift(walsh, 0)
    off_corners = invert_walsh(off_walsh, prime)
    off_target_spectra = tuple(C.spectrum(row, zeta13, prime) for row in off_walsh)
    off_corner_support = tuple(sum(value != 0 for value in row) for row in off_corners)
    off_walsh_support = tuple(sum(value != 0 for value in row) for row in off_walsh)
    off_fourier_support = tuple(sum(value != 0 for value in row) for row in off_target_spectra)
    off_target_rank = C.rank_mod(off_walsh, prime)
    require(off_corner_support == off_walsh_support == (12, 12, 12, 12), (
        off_corner_support,
        off_walsh_support,
    ))
    require(off_target_rank == 4, off_target_rank)

    edge_source_rows = tuple(difference_table[ell] for ell in range(1, 7))
    edge_source_spectra = tuple(C.spectrum(row, zeta13, prime) for row in edge_source_rows)
    reduced_target_spectra = off_target_spectra[1:]
    edge_to_reduced_certificate = generalized_connection_certificate(
        edge_source_spectra,
        reduced_target_spectra,
        prime,
    )

    full_certificate = C.connection_certificate(source_spectra, full_target_spectra, prime)
    off_certificate = C.connection_certificate(source_spectra, off_target_spectra, prime)
    require(full_certificate == (6, 24, 4, 0, 37, 4, 0, 0), full_certificate)

    dilation_certificates = tuple(
        C.connection_certificate(C.dilate(source_spectra, multiplier), off_target_spectra, prime)
        for multiplier in range(1, N)
    )
    require(all(row[3] == row[6] == row[7] == 0 for row in dilation_certificates), dilation_certificates)

    deletion_census = []
    for deleted in range(N):
        deleted_rows = delete_drift(walsh, deleted)
        deleted_spectra = tuple(C.spectrum(row, zeta13, prime) for row in deleted_rows)
        deletion_census.append(
            (
                deleted,
                C.rank_mod(deleted_rows, prime),
                tuple(sum(value != 0 for value in row) for row in deleted_spectra),
                C.connection_certificate(source_spectra, deleted_spectra, prime),
            )
        )
    deletion_census_t = tuple(deletion_census)

    proof = (
        "THM2471_same_root_disjointness_rebuilt_pointwise_on_273_entries",
        "typed_common_torsor_gauge_preserves_s_zero_iff_d_zero",
        "Cartesian_Ufull_has_four_nonzero_drift_zero_corner_buckets",
        "offdiagonal_projection_matches_punctured_F13_support_but_not_connection",
    )
    boundary = (
        "diagonal_deletion_is_necessary_not_a_proved_Ufull_ancestry_support_map",
        "punctured_source_6x12_and_target_4x12_are_representation_carriers_only",
        "no_ancestry_current_H1_bispectrum_scalar_exclusion_or_LRC14_claim",
    )
    semantic_surface = (
        COMMON_SHA256,
        (prime, zeta7, zeta13),
        digest(difference_table),
        source_row_support,
        source_rank,
        diagonal_corners,
        same_sheet_bridge,
        off_corner_support,
        off_walsh_support,
        off_fourier_support,
        off_target_rank,
        source_reflection,
        edge_to_reduced_certificate,
        full_certificate,
        off_certificate,
        dilation_certificates,
        deletion_census_t,
        proof,
        boundary,
    )
    semantic_sha256 = hashlib.sha256(repr(semantic_surface).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(
            semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
            (semantic_sha256, EXPECTED_SEMANTIC_SHA256),
        )

    source = Path(__file__).resolve()
    print("R5 OFF-DIAGONAL DRIFT BRANCH-TRANSPLANT PROBE")
    print("status=FINITE-EXACT support/representation sidecar;Ufull_ancestry_not_constructed;LRC14_OPEN")
    print(f"dependency={COMMON_PATH.name}:{COMMON_SHA256}")
    print(f"field=(prime={prime},zeta7={zeta7},zeta13={zeta13})")
    print("typed_skeleton=(source_s=u-q,target_d=b-a,common_gauge_u=a+c_and_q=b+c_implies_s=-d;orientation_reversal_is_lawful;zero_fibre_fixed)")
    print(f"source_same_root=(pointwise_zero={sum(value == 0 for value in pointwise_same_root)}/273,row_support={source_row_support},punctured_support={sum(source_row_support)})")
    print(f"source_difference=(sha256={digest(difference_table)},rank={source_rank},Fourier_support={tuple(sum(value != 0 for value in row) for row in source_spectra)})")
    print(f"source_edge_involution=(B_ell_s_equals_B_minusell_minuss={source_reflection},nonzero_ell_channels=6)")
    print(f"target_drift_zero=(corner_values={diagonal_corners},all_nonzero=True,same_sheet_bridge={same_sheet_bridge})")
    print(f"offdiagonal_target=(corner_support={off_corner_support},Walsh_support={off_walsh_support},Fourier_support={off_fourier_support},rank={off_target_rank})")
    print(f"connection_certificates=(full={full_certificate},offdiagonal={off_certificate})")
    print(f"edge_to_reduced_K4_certificate=(fields=source_rank,target_rank,wedge_rank,wedge_nullity,wedge_excess,augmented_rank,augmented_nullity,lambda_projection_rank,nonannihilator_entries;value={edge_to_reduced_certificate})")
    print(f"offdiagonal_dilations=(count={len(dilation_certificates)},unique={tuple(sorted(set(dilation_certificates)))})")
    print(f"single_drift_deletion_census={deletion_census_t}")
    print(f"proof={proof}")
    print(f"boundary={boundary}")
    print(f"semantic_sha256={semantic_sha256}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("reproducibility=no_assert_truth_gates;no_randomness;no_elapsed_fields")
    print(f"commands=python -B {SCRIPT};python -B -O {SCRIPT}")
    print("PASS")


if __name__ == "__main__":
    main()
