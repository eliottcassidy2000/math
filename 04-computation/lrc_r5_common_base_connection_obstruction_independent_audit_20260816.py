#!/usr/bin/env python3
"""Clean-room audit of the R5 common-base connection obstructions.

This script imports neither the fixed-root candidate nor the root-difference
candidate.  It reconstructs their source marginals from the hash-pinned
THM-2594 joint Boolean table and rebuilds the U_full target from the
independently audited THM-3514 endpoint atoms.  The projective connection is
checked by an augmented linear system retaining all thirteen frequency
multipliers; a separate wedge system is kept only as a cross-check.

Scope: finite exact split-field representation statements.  No common U_full
ancestry map, Boolean support transplant, physical current, absolute H1,
bispectrum, scalar exclusion, or LRC(14) theorem is constructed.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stdout
from itertools import combinations, product
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
PRIMARY_PATH = ROOT / "04-computation/lrc14_stage2_theta_contraction_opus_20260728.py"
ATOM_PATH = ROOT / (
    "04-computation/"
    "lrc_ufull_owner_boundary_k4xf13_endpoint_factorization_"
    "independent_audit_20260816.py"
)
PRIMARY_SHA256 = "09c43af0a0a56c7a0833bbfd13ed6a96bc5a7a3718aa1bc6b77a144bde101a06"
ATOM_SHA256 = "f89be10c65bb77270199f9399b155d5a2c82c0da121b3e8589fe3c1f7e9824fc"
EXPECTED_SLAVED_SHA256 = "666200f5e79e48c02e7beaeec90ebfedccb708bd1e383d1a01997d0abf5073f7"
EXPECTED_ABSOLUTE_SHA256 = "9353c6e645b570c659a5d42491a1dd7c79825d4bd338a476aabbd0afe10e69f2"
EXPECTED_DIFFERENCE_SHA256 = "ad63d008da3dfa8525db7ef4e724cf6e3a971b3f6efa3d6da9791fbd8003ed91"
EXPECTED_SEMANTIC_SHA256: str | None = "544b07a84c6806ea63c48f5227b78d74844f466dd0b446b6a320ea8560238895"

N = 13
CHAMBERS = ("left", "middle", "right")
ATOMS = tuple((sheet, chamber) for sheet in range(N) for chamber in CHAMBERS)
BUCKETS = tuple(
    (left, right, drift)
    for left in CHAMBERS
    for right in CHAMBERS
    for drift in range(N)
)
BUCKET_INDEX = {bucket: index for index, bucket in enumerate(BUCKETS)}
CORNERS = (("left", "left"), ("left", "right"), ("right", "left"), ("right", "right"))
Q_H = (1, 0, 1)
Q_Q5 = (1, 0, 0)


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def lf_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def digest(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def load_module(path: Path, name: str, expected_hash: str):
    actual = lf_sha256(path)
    require(actual == expected_hash, (name, actual, expected_hash))
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, (name, "loader"))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


ATOM = load_module(ATOM_PATH, "thm3514_atoms_for_common_base_audit", ATOM_SHA256)


def atom_worker(alpha: int) -> tuple[object, ...]:
    return ATOM.worker(alpha)


def joint_worker() -> tuple[object, ...]:
    primary = load_module(PRIMARY_PATH, "thm2594_joint_for_common_base_audit", PRIMARY_SHA256)
    buffer = io.StringIO()
    with redirect_stdout(buffer):
        state = primary.main()
    joint = state[12]
    denominator = int(state[17])
    frozen = tuple(
        tuple(
            tuple(tuple(int(value) for value in theta_row) for theta_row in ell_rows)
            for ell_rows in q_rows
        )
        for q_rows in joint
    )
    return frozen, denominator


def rref(
    matrix: tuple[tuple[int, ...], ...], prime: int
) -> tuple[tuple[tuple[int, ...], ...], tuple[int, ...]]:
    rows = [list(value % prime for value in row) for row in matrix]
    if not rows:
        return (), ()
    columns = len(rows[0])
    require(all(len(row) == columns for row in rows), "ragged matrix")
    pivots: list[int] = []
    pivot_row = 0
    for column in range(columns):
        pivot = next((row for row in range(pivot_row, len(rows)) if rows[row][column]), None)
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        inverse = pow(rows[pivot_row][column], -1, prime)
        rows[pivot_row] = [value * inverse % prime for value in rows[pivot_row]]
        for row in range(len(rows)):
            if row == pivot_row:
                continue
            factor = rows[row][column]
            if factor:
                rows[row] = [
                    (value - factor * pivot_value) % prime
                    for value, pivot_value in zip(rows[row], rows[pivot_row])
                ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return tuple(tuple(row) for row in rows), tuple(pivots)


def rank_mod(matrix: tuple[tuple[int, ...], ...], prime: int) -> int:
    return len(rref(matrix, prime)[1])


def nullspace(matrix: tuple[tuple[int, ...], ...], prime: int) -> tuple[tuple[int, ...], ...]:
    reduced, pivots = rref(matrix, prime)
    columns = len(reduced[0]) if reduced else 0
    free = tuple(column for column in range(columns) if column not in pivots)
    basis = []
    for free_column in free:
        vector = [0] * columns
        vector[free_column] = 1
        for row, pivot in enumerate(pivots):
            vector[pivot] = -reduced[row][free_column] % prime
        basis.append(tuple(vector))
    return tuple(basis)


def spectrum(values: tuple[int, ...], root: int, prime: int) -> tuple[int, ...]:
    return tuple(
        sum(value * pow(root, (-frequency * index) % len(values), prime) for index, value in enumerate(values))
        % prime
        for frequency in range(len(values))
    )


def septimal_profiles(
    table: tuple[tuple[int, ...], ...], zeta7: int, prime: int
) -> tuple[tuple[int, ...], ...]:
    require(len(table) == 7 and all(len(row) == N for row in table), "7x13 table")
    return tuple(
        tuple(
            sum(table[ell][column] * pow(zeta7, (-mode * ell) % 7, prime) for ell in range(7))
            % prime
            for column in range(N)
        )
        for mode in range(7)
    )


def materialize_target(chunks: tuple[tuple[object, ...], ...]) -> tuple[object, ...]:
    require(tuple(ATOM.ATOMS) == ATOMS and tuple(ATOM.BUCKETS) == BUCKETS, "atom order")
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(N)), "worker order")
    tables = tuple(row for chunk in chunks for row in chunk[1])
    require(len(tables) == N * N, len(tables))
    (_word, _t_den, nn, prime, root, zeta13, *_rest) = ATOM.context()
    raw = {Q_H: [0] * len(BUCKETS), Q_Q5: [0] * len(BUCKETS)}
    table_index = 0
    for alpha in range(N):
        for beta in range(N):
            stored_beta, ax_values, by_values = tables[table_index]
            table_index += 1
            require(stored_beta == beta, (alpha, beta, stored_beta))
            for tau in range(N):
                left_atoms = tuple(
                    (index, value)
                    for index, value in enumerate(ax_values)
                    if value and ATOM.safe(ATOMS[index][1], ATOMS[index][0] + tau)
                )
                right_atoms = tuple(
                    (index, value)
                    for index, value in enumerate(by_values)
                    if value and ATOM.safe(ATOMS[index][1], ATOMS[index][0] + tau)
                )
                slices = [0] * len(BUCKETS)
                for left_index, left_value in left_atoms:
                    left_sheet, left_chamber = ATOMS[left_index]
                    for right_index, right_value in right_atoms:
                        right_sheet, right_chamber = ATOMS[right_index]
                        bucket = BUCKET_INDEX[(left_chamber, right_chamber, (right_sheet - left_sheet) % N)]
                        slices[bucket] = (slices[bucket] + left_value * right_value) % prime
                for q_class in (Q_H, Q_Q5):
                    qa, qb, qt = q_class
                    phase = pow(zeta13, (beta - alpha * qa - beta * qb - tau * qt) % N, prime)
                    for bucket, value in enumerate(slices):
                        raw[q_class][bucket] = (raw[q_class][bucket] + phase * value) % prime
    normalizer = pow(N**3, -1, prime)
    buckets = {
        q_class: tuple(value * normalizer % prime for value in raw[q_class])
        for q_class in (Q_H, Q_Q5)
    }
    bridge = tuple((left - right) % prime for left, right in zip(buckets[Q_H], buckets[Q_Q5]))
    corners = tuple(
        tuple(bridge[BUCKET_INDEX[(left, right, drift)]] for drift in range(N))
        for left, right in CORNERS
    )
    walsh = (
        tuple(sum(corners[row][drift] for row in range(4)) % prime for drift in range(N)),
        tuple((corners[0][d] + corners[1][d] - corners[2][d] - corners[3][d]) % prime for d in range(N)),
        tuple((corners[0][d] - corners[1][d] + corners[2][d] - corners[3][d]) % prime for d in range(N)),
        tuple((corners[0][d] - corners[1][d] - corners[2][d] + corners[3][d]) % prime for d in range(N)),
    )
    target_spectra = tuple(spectrum(row, zeta13, prime) for row in walsh)
    require(rank_mod(walsh, prime) == 4, "target rank")
    require(all(all(row) for row in target_spectra), "target spectral zero")
    require(sum(walsh[0]) % prime == 389266878372286537904, "bridge regression")
    return nn, prime, root, zeta13, walsh, target_spectra


def wedge_system(
    source_rows: tuple[tuple[int, ...], ...],
    target_rows: tuple[tuple[int, ...], ...],
    prime: int,
) -> tuple[tuple[int, ...], ...]:
    channels = len(source_rows)
    equations = []
    for frequency in range(N):
        source = tuple(row[frequency] for row in source_rows)
        target = tuple(row[frequency] for row in target_rows)
        for left, right in combinations(range(4), 2):
            equation = [0] * (4 * channels)
            for column in range(channels):
                equation[right * channels + column] += target[left] * source[column]
                equation[left * channels + column] -= target[right] * source[column]
            equations.append(tuple(value % prime for value in equation))
    return tuple(equations)


def augmented_system(
    source_rows: tuple[tuple[int, ...], ...],
    target_rows: tuple[tuple[int, ...], ...],
    prime: int,
) -> tuple[tuple[int, ...], ...]:
    channels = len(source_rows)
    equations = []
    for frequency in range(N):
        source = tuple(row[frequency] for row in source_rows)
        target = tuple(row[frequency] for row in target_rows)
        for output in range(4):
            equation = [0] * (4 * channels + N)
            for column in range(channels):
                equation[output * channels + column] = source[column]
            equation[4 * channels + frequency] = -target[output] % prime
            equations.append(tuple(equation))
    return tuple(equations)


def connection_certificate(
    source_rows: tuple[tuple[int, ...], ...],
    target_rows: tuple[tuple[int, ...], ...],
    prime: int,
) -> tuple[int, ...]:
    channels = len(source_rows)
    source_rank = rank_mod(source_rows, prime)
    wedge = wedge_system(source_rows, target_rows, prime)
    wedge_rank = rank_mod(wedge, prime)
    wedge_nullity = 4 * channels - wedge_rank
    annihilator_dimension = 4 * (channels - source_rank)
    wedge_excess = wedge_nullity - annihilator_dimension
    require(wedge_excess >= 0, (channels, wedge_nullity, annihilator_dimension))

    augmented = augmented_system(source_rows, target_rows, prime)
    augmented_rank = rank_mod(augmented, prime)
    basis = nullspace(augmented, prime)
    augmented_nullity = 4 * channels + N - augmented_rank
    require(len(basis) == augmented_nullity, (len(basis), augmented_nullity))
    multiplier_rank = rank_mod(tuple(vector[4 * channels :] for vector in basis), prime)
    nonannihilator = 0
    for vector in basis:
        for output in range(4):
            for frequency in range(N):
                value = sum(
                    vector[output * channels + column] * source_rows[column][frequency]
                    for column in range(channels)
                ) % prime
                nonannihilator += value != 0
    require((wedge_excess > 0) == (multiplier_rank > 0), (wedge_excess, multiplier_rank))
    return (
        source_rank,
        wedge_rank,
        wedge_nullity,
        wedge_excess,
        augmented_rank,
        augmented_nullity,
        multiplier_rank,
        nonannihilator,
    )


def dilate(rows: tuple[tuple[int, ...], ...], multiplier: int) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(row[multiplier * frequency % N] for frequency in range(N)) for row in rows)


def canonical_ratio(row: tuple[int, ...], zeta13: int, prime: int) -> tuple[int, ...]:
    require(all(row), "ratio must be a unit vector")
    amplitude = pow(row[0], -1, prime)
    candidates = tuple(
        tuple(
            row[frequency]
            * amplitude
            % prime
            * pow(zeta13, shift * frequency % N, prime)
            % prime
            for frequency in range(N)
        )
        for shift in range(N)
    )
    return min(candidates)


def main() -> None:
    with ProcessPoolExecutor(max_workers=5) as pool:
        joint_future = pool.submit(joint_worker)
        atom_futures = tuple(pool.submit(atom_worker, alpha) for alpha in range(N))
        joint, denominator = joint_future.result()
        chunks = tuple(future.result() for future in atom_futures)

    nn, prime, root, zeta13, walsh, target_spectra = materialize_target(chunks)
    require(nn % 7 == 0 and denominator % prime, (nn, denominator % prime))
    zeta7 = pow(root, nn // 7, prime)
    require(pow(zeta7, 7, prime) == 1 and zeta7 != 1, "zeta7")

    slaved = [[0] * N for _ in range(7)]
    absolute = [[0] * N for _ in range(7)]
    difference = [[0] * N for _ in range(7)]
    resolved = [[[0] * N for _ in range(7)] for _ in range(3)]
    owner = [[[0] * N for _ in range(7)] for _ in range(2)]
    for u in range(N):
        for q in range(N):
            root_difference = (u - q) % N
            for ell in range(7):
                for theta in range(3):
                    value = joint[u][q][ell][theta] % prime
                    slaved[ell][theta] = (slaved[ell][theta] + value) % prime
                    absolute_root = (theta + 2 * u) % N
                    absolute[ell][absolute_root] = (absolute[ell][absolute_root] + value) % prime
                    difference[ell][root_difference] = (difference[ell][root_difference] + value) % prime
                    resolved[theta][ell][root_difference] = (
                        resolved[theta][ell][root_difference] + value
                    ) % prime
                    owner[0][ell][root_difference] = (owner[0][ell][root_difference] + value) % prime
                    owner[1][ell][root_difference] = (
                        owner[1][ell][root_difference]
                        + value * pow(zeta13, -u % N, prime)
                    ) % prime

    slaved_table = tuple(tuple(row) for row in slaved)
    absolute_table = tuple(tuple(row) for row in absolute)
    difference_table = tuple(tuple(row) for row in difference)
    require(digest(slaved_table) == EXPECTED_SLAVED_SHA256, "slaved digest")
    require(digest(absolute_table) == EXPECTED_ABSOLUTE_SHA256, "absolute digest")
    require(digest(difference_table) == EXPECTED_DIFFERENCE_SHA256, "difference digest")
    require(all(difference_table[ell][0] == 0 for ell in range(7)), "same-root column")

    slaved_profiles = septimal_profiles(slaved_table, zeta7, prime)
    absolute_profiles = septimal_profiles(absolute_table, zeta7, prime)
    difference_profiles = septimal_profiles(difference_table, zeta7, prime)
    slaved_spectra = tuple(spectrum(row, zeta13, prime) for row in slaved_profiles)
    absolute_spectra = tuple(spectrum(row, zeta13, prime) for row in absolute_profiles)
    difference_spectra = tuple(spectrum(row, zeta13, prime) for row in difference_profiles)

    resolved_profiles = tuple(
        row
        for theta in range(3)
        for row in septimal_profiles(tuple(tuple(values) for values in resolved[theta]), zeta7, prime)
    )
    resolved_spectra = tuple(spectrum(row, zeta13, prime) for row in resolved_profiles)
    owner_profiles = tuple(
        septimal_profiles(tuple(tuple(values) for values in owner[k]), zeta7, prime)
        for k in range(2)
    )
    owner_spectra = tuple(
        tuple(spectrum(row, zeta13, prime) for row in owner_profiles[k])
        for k in range(2)
    )
    require(owner_profiles[0] == difference_profiles, "owner k0")
    owner_difference_profiles = tuple(
        tuple((owner_profiles[1][mode][column] - owner_profiles[0][mode][column]) % prime for column in range(N))
        for mode in range(7)
    )
    owner_difference_spectra = tuple(spectrum(row, zeta13, prime) for row in owner_difference_profiles)
    owner_union_spectra = owner_spectra[0] + owner_spectra[1]

    raw_supports = tuple(
        sum(value != 0 for row in table for value in row)
        for table in (slaved_table, absolute_table, difference_table)
    )
    spectral_supports = tuple(
        sum(value != 0 for row in bank for value in row)
        for bank in (slaved_spectra, absolute_spectra, difference_spectra)
    )
    ranks = (
        rank_mod(slaved_profiles, prime),
        rank_mod(absolute_profiles, prime),
        rank_mod(slaved_profiles + absolute_profiles, prime),
        rank_mod(difference_profiles, prime),
        rank_mod(resolved_profiles, prime),
        rank_mod(owner_profiles[0], prime),
        rank_mod(owner_profiles[1], prime),
        rank_mod(owner_difference_profiles, prime),
        rank_mod(owner_profiles[0] + owner_profiles[1], prime),
    )
    require(raw_supports == (12, 18, 72), raw_supports)
    require(spectral_supports == (91, 91, 91), spectral_supports)
    require(ranks == (3, 3, 4, 6, 8, 6, 6, 6, 9), ranks)

    banks = {
        "slaved7": slaved_spectra,
        "absolute7": absolute_spectra,
        "union14": slaved_spectra + absolute_spectra,
        "difference7": difference_spectra,
        "resolved21": resolved_spectra,
        "owner_k0": owner_spectra[0],
        "owner_k1": owner_spectra[1],
        "owner_delta": owner_difference_spectra,
        "owner_union14": owner_union_spectra,
    }
    certificates = {name: connection_certificate(bank, target_spectra, prime) for name, bank in banks.items()}
    require(all(row[3] == row[6] == row[7] == 0 for row in certificates.values()), certificates)

    dilation_certificates = {
        name: tuple(connection_certificate(dilate(banks[name], multiplier), target_spectra, prime) for multiplier in range(1, N))
        for name in ("union14", "difference7", "resolved21", "owner_union14")
    }
    require(
        all(all(row[3] == row[6] == row[7] == 0 for row in rows) for rows in dilation_certificates.values()),
        dilation_certificates,
    )

    folded_pairs = ((1, 6), (2, 5), (3, 4))
    difference_selected = tuple(
        (
            bits,
            (0,) + tuple(folded_pairs[index][bits[index]] for index in range(3)),
        )
        for bits in product((0, 1), repeat=3)
    )
    difference_selected_certificates = tuple(
        (bits, modes, connection_certificate(tuple(difference_spectra[mode] for mode in modes), target_spectra, prime))
        for bits, modes in difference_selected
    )
    require(all(row[2] == (4, 16, 0, 0, 29, 0, 0, 0) for row in difference_selected_certificates), difference_selected_certificates)

    named_cases = []
    for absolute_mode in range(7):
        for bits in product((0, 1), repeat=3):
            modes = tuple(folded_pairs[index][bits[index]] for index in range(3))
            bank = tuple(slaved_spectra[mode] for mode in modes) + (absolute_spectra[absolute_mode],)
            named_cases.append(connection_certificate(bank, target_spectra, prime))
    binary_cases = []
    for mask in range(1, 1 << 7):
        sidecar = tuple(
            sum(absolute_spectra[mode][frequency] for mode in range(7) if mask & (1 << mode)) % prime
            for frequency in range(N)
        )
        bank = tuple(slaved_spectra[mode] for mode in (1, 2, 3)) + (sidecar,)
        binary_cases.append((mask, connection_certificate(bank, target_spectra, prime)))
    require(len(named_cases) == 56 and all(row[0] == 4 and row[3] == row[6] == 0 for row in named_cases), "named sidecars")
    require(len(binary_cases) == 127, len(binary_cases))
    require(tuple(mask for mask, row in binary_cases if row[0] == 3) == (127,), "binary rank-three mask")
    require(all(row[3] == row[6] == 0 for _mask, row in binary_cases), "binary projective hostile")

    canonical_classes = []
    exact_common = 0
    gauge_common = 0
    for modes in product(range(7), repeat=4):
        ratios = tuple(
            tuple(
                target_spectra[channel][frequency]
                * pow(difference_spectra[modes[channel]][frequency], -1, prime)
                % prime
                for frequency in range(N)
            )
            for channel in range(4)
        )
        exact_common += len(set(ratios)) == 1
        classes = len({canonical_ratio(row, zeta13, prime) for row in ratios})
        canonical_classes.append(classes)
        gauge_common += classes == 1
    allocation_census = (7**4, exact_common, gauge_common, tuple(sorted(set(canonical_classes))))
    require(allocation_census == (2401, 0, 0, (4,)), allocation_census)

    proof = (
        "source_marginals_rebuilt_from_THM2594_joint_table",
        "target_rebuilt_from_THM3514_independent_endpoint_atoms",
        "augmented_M_lambda_system_retains_all_frequency_multipliers",
        "every_nullspace_multiplier_projection_is_zero",
        "every_augmented_nullvector_annihilates_the_source_curve",
        "rank4_6_8_9_is_necessary_but_not_sufficient_for_common_connection",
    )
    boundary = (
        "failed_class=fixed_frequency_independent_channel_map_plus_one_common_C13_operator",
        "survivor=joint_u_s_ell_theta_to_a_d_C_D_support_relation_before_marginalization",
        "no_ancestry_current_H1_bispectrum_scalar_exclusion_or_LRC14_claim",
    )
    semantic_surface = (
        PRIMARY_SHA256,
        ATOM_SHA256,
        (prime, zeta7, zeta13),
        denominator,
        (digest(slaved_table), digest(absolute_table), digest(difference_table)),
        raw_supports,
        spectral_supports,
        ranks,
        certificates,
        dilation_certificates,
        difference_selected_certificates,
        tuple(named_cases),
        tuple(binary_cases),
        allocation_census,
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
    print("R5 COMMON-BASE CONNECTION OBSTRUCTIONS -- INDEPENDENT AUDIT")
    print("status=ACCEPT_SCOPED_FINITE_EXACT;split_field_representation_only;LRC14_OPEN")
    print(f"dependencies=(thm2594_primary={PRIMARY_SHA256},thm3514_atom_audit={ATOM_SHA256})")
    print(f"field=(prime={prime},zeta7={zeta7},zeta13={zeta13})")
    print(f"source_tables=(digests={(digest(slaved_table),digest(absolute_table),digest(difference_table))},raw_supports={raw_supports},spectral_supports={spectral_supports})")
    print(f"source_ranks=(slaved,absolute,union,difference,resolved,owner0,owner1,owner_delta,owner_union)={ranks}")
    print("certificate_fields=(source_rank,wedge_rank,wedge_nullity,wedge_excess,augmented_rank,augmented_nullity,lambda_projection_rank,nonannihilator_entries)")
    print(f"whole_bank_certificates={certificates}")
    print(f"dilation_certificates={dict((name,(len(rows),tuple(sorted(set(rows))))) for name,rows in dilation_certificates.items())}")
    print(f"difference_folded_sections={difference_selected_certificates}")
    print(f"absolute_sidecar_census=(named_cases={len(named_cases)},named_rank4={sum(row[0]==4 for row in named_cases)},binary_cases={len(binary_cases)},binary_rank4={sum(row[0]==4 for _,row in binary_cases)},binary_rank3_masks={tuple(mask for mask,row in binary_cases if row[0]==3)})")
    print(f"difference_allocation_census={allocation_census}")
    print("augmented_verdict=all_lambda_projections_zero_and_all_nullvectors_source_annihilators")
    print(f"proof={proof}")
    print(f"boundary={boundary}")
    print(f"semantic_sha256={semantic_sha256}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("reproducibility=no_assert_truth_gates;candidate_scripts_not_imported;no_randomness;no_elapsed_fields")
    print("commands=python -B 04-computation/lrc_r5_common_base_connection_obstruction_independent_audit_20260816.py;python -B -O 04-computation/lrc_r5_common_base_connection_obstruction_independent_audit_20260816.py")
    print("PASS")


if __name__ == "__main__":
    main()
