#!/usr/bin/env python3
"""Retain the seven base cells in the source-aligned endpoint pullback.

The source dependency realizes the THM-2471 source-time refinement

    f_omega^src = 1_Q P^K(e P_omega),    e_nu = e P_nu,

on one first-collision base and retains the common F_13 torsor offset c.
This companion inserts THM-2594's seven-cell partition before the final
base integral.  It therefore constructs the finite exact tensor

    M(omega,nu,ell,c),    ell in F_7, c in F_13,

and contracts it against every matrix in THM-3479's refined endpoint
residue bank A_w^-(1,0,t), t in F_13.  The fixed all-unit relation is t=6.

This is a source-time finite simple-kernel 7 x 13 spectrum.  The endpoint
AX/BY values are still preintegrated scalars: no exact-address coefficient,
physical relation current, THM-2512 bridge theorem, row exclusion, or
LRC(14) conclusion is asserted.
"""

from __future__ import annotations

from contextlib import redirect_stdout
from hashlib import sha256
import importlib
import io
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
P = 13
Q = 7
SOURCE_PATH = (
    ROOT
    / "04-computation/lrc_r5_source_aligned_guard_atom_branch_sidecar_probe_20260816.py"
)
SOURCE_SHA256 = "22c5c748392817ccc36889a007c65bd5f44b26c10638df6f6aac48e917547f41"
SOURCE_SEMANTIC = "31e9e90c63053944b590195555be07ccbf84fd4c7abc2101de6a2a3562202de6"
ENDPOINT_PATH = (
    ROOT
    / "04-computation/lrc_r5_source_aligned_actual_endpoint_simple_kernel_probe_20260816.py"
)
ENDPOINT_SHA256 = "fe8fba15e1389711ac569e118a31efbca550c90d69bc95c4d5804f4f2fc73f11"
ENDPOINT_SEMANTIC = "f496b9da968f091b11d27137acdd21c4351785943c45b0e7783fce95b5915df0"
ENDPOINT_PAIR_BANK_SHA256 = "c28119c8b54f47e5b7a46f1508fbba604b0e3997eaadb05b03ad28edd9aed468"
EXPECTED_TENSOR_SHA256 = "39d7a0b4e5b2d8b85631d682ed1967091e44dc41e17b33a77e7184d3dc93e0cf"
EXPECTED_COORDINATE_BANK_SHA256 = "989dafc220a6d09aeacfce4af0e9a4fe13eedacc79fa66032ea39bc107fd8efb"
EXPECTED_SPECTRUM_BANK_SHA256 = "5f173227c5e203309f61bdfd9d47cc64a3b49ae8f14abd0f7bfc469eda278533"
EXPECTED_SEMANTIC_SHA256 = "cd55336bb1dfe5f37f020c242c4bca5b7c6be339ec57e95d69e10bbe68d9dbaa"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    body = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return sha256(body).hexdigest()


def load_endpoint_module():
    require(lf_sha256(ENDPOINT_PATH) == ENDPOINT_SHA256, "endpoint source drift")
    computation = str(ENDPOINT_PATH.parent)
    if computation not in sys.path:
        sys.path.insert(0, computation)
    module = importlib.import_module(ENDPOINT_PATH.stem)
    require(module.EXPECTED_SEMANTIC_SHA256 == ENDPOINT_SEMANTIC,
            "endpoint semantic pin drift")
    return module


def cell_partition_record(cells, grid: int) -> tuple[object, ...]:
    flattened = sorted(interval for cell in cells for interval in cell)
    require(len(cells) == Q, ("cell count", len(cells)))
    require(all(sum(right - left for left, right in cell) == grid // Q
                for cell in cells), "cell masses")
    require(flattened[0][0] == 0 and flattened[-1][1] == grid,
            "cell partition endpoints")
    require(all(left[1] == right[0]
                for left, right in zip(flattened, flattened[1:])),
            "cell partition gap/overlap")
    return (
        tuple(tuple(tuple(interval) for interval in cell) for cell in cells),
        tuple(tuple(interval) for interval in flattened),
        tuple(sum(right - left for left, right in cell) for cell in cells),
    )


def build_cell_tensor(source, source_data):
    """Build M[omega][nu][ell][c] before the final base marginal."""
    base = source.B
    stage = source.M
    grid = source.GRID
    atoms = source.ATOMS

    e_intervals = stage.build_set(stage.PAT_E, stage.ZELL)
    q_intervals = stage.build_set(stage.PAT_QB, stage.ZELL)
    endpoint_groups = base.partition_weighted_pieces(
        [(left, right, 1) for left, right in e_intervals]
    )
    source_groups = tuple(
        source.source_current_group(group, q_intervals) for group in endpoint_groups
    )
    u_profiles = tuple(
        stage.weighted_fold(list(group), stage.DCOLL, grid)
        for group in source_groups
    )
    v_profiles = tuple(
        stage.weighted_fold(list(group), stage.DCOLL, grid)
        for group in endpoint_groups
    )
    u_windows = tuple(
        tuple(stage.extract_window(starts, values, root, P, grid)
              for root in range(P))
        for starts, values in u_profiles
    )
    v_windows = tuple(
        tuple(stage.extract_window(starts, values, root, P, grid)
              for root in range(P))
        for starts, values in v_profiles
    )

    cells = stage.build_cells()
    cell_record = cell_partition_record(cells, grid)
    tensor = [
        [[[0 for _offset in range(P)] for _ell in range(Q)] for _right in atoms]
        for _left in atoms
    ]

    for left_index, (left_sheet, _left_chamber) in enumerate(atoms):
        for right_index, (right_sheet, _right_chamber) in enumerate(atoms):
            drift = (right_sheet - left_sheet) % P
            root_drift = (-drift) % P
            for current_root in range(P):
                source_root = (current_root - root_drift) % P
                profile = stage.product_cum(
                    u_windows[left_index][current_root][0],
                    u_windows[left_index][current_root][1],
                    v_windows[right_index][source_root][0],
                    v_windows[right_index][source_root][1],
                    grid,
                )
                masses = tuple(stage.set_integral(profile, cell) for cell in cells)
                require(sum(masses) == profile[3],
                        ("cell reconstruction", left_index, right_index, current_root))
                offset = (left_sheet - current_root) % P
                require(offset == (right_sheet - source_root) % P,
                        "common torsor offset mismatch")
                for ell, mass in enumerate(masses):
                    tensor[left_index][right_index][ell][offset] += mass

    marginal = [
        [
            [sum(tensor[left][right][ell][offset] for ell in range(Q))
             for offset in range(P)]
            for right in range(len(atoms))
        ]
        for left in range(len(atoms))
    ]
    require(marginal == source_data["pair_gauge_offset"],
            "seven-cell marginal did not recover source tensor")
    pair_support = sum(
        any(tensor[left][right][ell][offset]
            for ell in range(Q) for offset in range(P))
        for left in range(len(atoms)) for right in range(len(atoms))
    )
    entry_support = sum(
        tensor[left][right][ell][offset] != 0
        for left in range(len(atoms)) for right in range(len(atoms))
        for ell in range(Q) for offset in range(P)
    )
    require(pair_support == 362, ("cell pair support", pair_support))
    return tensor, marginal, cell_record, pair_support, entry_support


def contract(
    tensor: list[list[list[list[int]]]],
    pair_function: list[list[int]],
    denominator: int,
    prime: int,
) -> tuple[tuple[int, ...], ...]:
    inverse_denominator = pow(denominator, -1, prime)
    answer = [[0 for _offset in range(P)] for _ell in range(Q)]
    for left in range(len(pair_function)):
        for right in range(len(pair_function)):
            weight = pair_function[left][right]
            if not weight:
                continue
            for ell in range(Q):
                for offset in range(P):
                    answer[ell][offset] = (
                        answer[ell][offset]
                        + tensor[left][right][ell][offset] * weight
                    ) % prime
    return tuple(
        tuple(value * inverse_denominator % prime for value in row)
        for row in answer
    )


def fourier_2d(matrix, eta: int, zeta: int, prime: int):
    eta_phase = tuple(
        tuple(pow(eta, -h * ell % Q, prime) for ell in range(Q))
        for h in range(Q)
    )
    zeta_phase = tuple(
        tuple(pow(zeta, -k * offset % P, prime) for offset in range(P))
        for k in range(P)
    )
    return tuple(
        tuple(
            sum(
                matrix[ell][offset]
                * eta_phase[h][ell]
                * zeta_phase[k][offset]
                for ell in range(Q) for offset in range(P)
            ) % prime
            for k in range(P)
        )
        for h in range(Q)
    )


def support_shape(spectrum) -> tuple[int, int, int, int, int]:
    dc = int(spectrum[0][0] != 0)
    cell_axis = sum(spectrum[h][0] != 0 for h in range(1, Q))
    owner_axis = sum(spectrum[0][k] != 0 for k in range(1, P))
    mixed = sum(
        spectrum[h][k] != 0 for h in range(1, Q) for k in range(1, P)
    )
    return dc + cell_axis + owner_axis + mixed, dc, cell_axis, owner_axis, mixed


def matrix_interaction(matrix, prime: int):
    inv_q = pow(Q, -1, prime)
    inv_p = pow(P, -1, prime)
    inv_qp = pow(Q * P, -1, prime)
    row_sums = tuple(sum(row) % prime for row in matrix)
    col_sums = tuple(
        sum(matrix[ell][offset] for ell in range(Q)) % prime
        for offset in range(P)
    )
    grand = sum(row_sums) % prime
    interaction = tuple(
        tuple(
            (
                matrix[ell][offset]
                - row_sums[ell] * inv_p
                - col_sums[offset] * inv_q
                + grand * inv_qp
            ) % prime
            for offset in range(P)
        )
        for ell in range(Q)
    )
    require(all(sum(row) % prime == 0 for row in interaction),
            "cell-owner interaction row centering")
    require(all(sum(interaction[ell][offset] for ell in range(Q)) % prime == 0
                for offset in range(P)), "cell-owner interaction column centering")
    return interaction


def main() -> None:
    endpoint = load_endpoint_module()
    source = endpoint.load_module(
        SOURCE_PATH, "source_aligned_cell_spectrum_source", SOURCE_SHA256
    )
    with redirect_stdout(io.StringIO()):
        source_data = source.main()
    require(source_data["semantic_sha256"] == SOURCE_SEMANTIC,
            "source semantic drift")

    tensor, marginal, cell_record, pair_support, entry_support = (
        build_cell_tensor(source, source_data)
    )
    tensor_digest = digest_json(tensor)
    require(tensor_digest == EXPECTED_TENSOR_SHA256,
            ("cell tensor", tensor_digest))
    require(digest_json(marginal) == digest_json(source_data["pair_gauge_offset"]),
            "marginal digest mismatch")

    target = endpoint.load_module(
        endpoint.TARGET_PATH,
        "source_aligned_cell_spectrum_target",
        endpoint.TARGET_SHA256,
    )
    _bridge_pair, pair_bank, target_record = endpoint.build_pair_weight(target)
    prime, zeta, _bridge, _support, _pair_digest, endpoint_totals, bank_digest = (
        target_record
    )
    require(bank_digest == ENDPOINT_PAIR_BANK_SHA256, ("pair bank", bank_digest))
    require(source_data["denominator"] % prime != 0, "source denominator mod p")

    (_word, _t_den, nn, context_prime, root, context_zeta, *_rest) = target.context()
    require((prime, zeta) == (context_prime, context_zeta), "target context drift")
    require(nn % (Q * P) == 0, ("91 does not divide embedding order", nn))
    eta = pow(root, nn // Q, prime)
    require(pow(eta, Q, prime) == 1 and eta != 1, "order-seven root")
    require(pow(zeta, P, prime) == 1 and zeta != 1, "order-thirteen root")

    coordinate_bank = tuple(
        contract(tensor, matrix, source_data["denominator"], prime)
        for matrix in pair_bank
    )
    spectrum_bank = tuple(
        fourier_2d(matrix, eta, zeta, prime) for matrix in coordinate_bank
    )
    support_shapes = tuple(support_shape(spectrum) for spectrum in spectrum_bank)
    require(support_shapes == ((91, 1, 6, 12, 72),) * P,
            ("residue 7x13 spectral closure", support_shapes))

    inherited_owner_bank = tuple(
        endpoint.pull_profile(
            source_data["pair_gauge_offset"],
            matrix,
            source_data["denominator"],
            zeta,
            prime,
        )
        for matrix in pair_bank
    )
    require(tuple(spectrum[0] for spectrum in spectrum_bank) == inherited_owner_bank,
            "cell marginal did not recover all thirteen owner profiles")
    require(all(all(value != 0 for value in profile)
                for profile in inherited_owner_bank),
            "inherited residue-owner nonvanishing")

    relation_t = 6
    relation_matrix = coordinate_bank[relation_t]
    relation_spectrum = spectrum_bank[relation_t]
    relation_interaction = matrix_interaction(relation_matrix, prime)
    relation_interaction_spectrum = fourier_2d(
        relation_interaction, eta, zeta, prime
    )
    require(all(relation_interaction_spectrum[0][k] == 0 for k in range(P)),
            "interaction owner axis")
    require(all(relation_interaction_spectrum[h][0] == 0 for h in range(Q)),
            "interaction cell axis")
    require(all(
        relation_interaction_spectrum[h][k] == relation_spectrum[h][k]
        for h in range(1, Q) for k in range(1, P)
    ), "mixed Fourier modes changed under ANOVA centering")

    _left, _right, endpoint_pair_interaction = endpoint.decompose_pair_function(
        pair_bank[relation_t], prime
    )
    endpoint_interaction_matrix = contract(
        tensor, endpoint_pair_interaction, source_data["denominator"], prime
    )
    endpoint_interaction_spectrum = fourier_2d(
        endpoint_interaction_matrix, eta, zeta, prime
    )

    # Hostile controls: erasing either retained coordinate kills every mixed mode.
    inv_q = pow(Q, -1, prime)
    cell_erased = tuple(
        tuple(
            sum(relation_matrix[ell][offset] for ell in range(Q)) * inv_q % prime
            for offset in range(P)
        )
        for _ell in range(Q)
    )
    inv_p = pow(P, -1, prime)
    owner_erased = tuple(
        tuple(sum(relation_matrix[ell]) * inv_p % prime for _offset in range(P))
        for ell in range(Q)
    )
    cell_erased_spectrum = fourier_2d(cell_erased, eta, zeta, prime)
    owner_erased_spectrum = fourier_2d(owner_erased, eta, zeta, prime)
    require(support_shape(cell_erased_spectrum)[4] == 0, "cell-erasure hostile")
    require(support_shape(owner_erased_spectrum)[4] == 0, "owner-erasure hostile")

    relation_shape = support_shape(relation_spectrum)
    relation_interaction_shape = support_shape(relation_interaction_spectrum)
    endpoint_interaction_shape = support_shape(endpoint_interaction_spectrum)
    require(relation_shape == (91, 1, 6, 12, 72),
            ("fixed relation raw spectrum", relation_shape))
    require(relation_interaction_shape == (72, 0, 0, 0, 72),
            ("fixed relation mixed spectrum", relation_interaction_shape))
    require(endpoint_interaction_shape == (91, 1, 6, 12, 72),
            ("endpoint interaction spectrum", endpoint_interaction_shape))

    coordinate_bank_digest = digest_json(coordinate_bank)
    spectrum_bank_digest = digest_json(spectrum_bank)
    require(coordinate_bank_digest == EXPECTED_COORDINATE_BANK_SHA256,
            ("coordinate bank", coordinate_bank_digest))
    require(spectrum_bank_digest == EXPECTED_SPECTRUM_BANK_SHA256,
            ("spectrum bank", spectrum_bank_digest))

    record = (
        SOURCE_SHA256,
        SOURCE_SEMANTIC,
        ENDPOINT_SHA256,
        ENDPOINT_SEMANTIC,
        endpoint.TARGET_SHA256,
        endpoint.TARGET_BUCKET_SHA256,
        tensor_digest,
        pair_support,
        entry_support,
        cell_record,
        source_data["denominator"],
        target_record,
        prime,
        eta,
        zeta,
        endpoint_totals,
        bank_digest,
        coordinate_bank_digest,
        spectrum_bank_digest,
        support_shapes,
        digest_json(inherited_owner_bank),
        relation_t,
        relation_matrix,
        relation_spectrum,
        relation_shape,
        relation_interaction,
        relation_interaction_spectrum,
        relation_interaction_shape,
        digest_json(endpoint_pair_interaction),
        endpoint_interaction_matrix,
        endpoint_interaction_spectrum,
        endpoint_interaction_shape,
        support_shape(cell_erased_spectrum),
        support_shape(owner_erased_spectrum),
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                (semantic, EXPECTED_SEMANTIC_SHA256))

    print("== r=5 source-aligned refined-residue 7 x 13 spectrum probe ==")
    print(f"dependencies=((source,{SOURCE_SHA256},{SOURCE_SEMANTIC}),(endpoint,{ENDPOINT_SHA256},{ENDPOINT_SEMANTIC}))")
    print(f"split_embedding=(prime={prime},zeta7={eta},zeta13={zeta})")
    print(f"seven_cells=exact_half_open_partition; masses={cell_record[2]}")
    print(f"source_cell_tensor=(pair_support={pair_support}/1521,entry_support={entry_support}/138411,sha256={tensor_digest})")
    print("cell_marginal_recovery=all 39x39x13 source-aligned common-offset entries: PASS")
    print(f"target_refined_(1,0,t)_pair_bank_sha256={bank_digest}")
    print("residue_owner_marginal_recovery=all 13x13 prior profiles: PASS")
    print(f"refined_residue_raw_spectrum_shapes_(total,dc,F7axis,F13axis,mixed)={support_shapes}")
    print(f"fixed_relation_refined_class=(1,0,{relation_t}); endpoint_total={endpoint_totals[relation_t]}")
    print(f"fixed_relation_owner_marginal={inherited_owner_bank[relation_t]}")
    print(f"fixed_relation_raw_7x13_spectrum_shape={relation_shape}")
    print(f"fixed_relation_cell_owner_ANOVA_spectrum_shape={relation_interaction_shape}")
    print(f"fixed_relation_endpoint_pair_ANOVA_spectrum_shape={endpoint_interaction_shape}")
    print(f"fixed_relation_mixed_witness_(h=1,k=1)={relation_spectrum[1][1]}")
    print(f"hostile_erasure_shapes=(cell_erased={support_shape(cell_erased_spectrum)},owner_erased={support_shape(owner_erased_spectrum)})")
    print(f"coordinate_bank_sha256={coordinate_bank_digest}")
    print(f"spectrum_bank_sha256={spectrum_bank_digest}")
    print(f"semantic_sha256={semantic}")
    print("scope=finite-exact source-time 7x13 simple-kernel spectrum; preintegrated AX/BY residue pushforward, not exact C(a;X,m),physical current,THM-2512 bridge,row exclusion,LRC(14)")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
