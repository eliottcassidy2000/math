#!/usr/bin/env python3
"""Expand the BY endpoint atom table before the r=5 source contraction.

The source-aligned 7 x 13 probe contracts a common-ancestry tensor

    M(omega, nu, ell, c)

against the U_full endpoint pair bank.  In that bank both endpoint factors
AX and BY were previously represented by preintegrated atom scalars.  This
companion keeps AX preintegrated but rebuilds every BY atom from the literal
half-open intervals used by the endpoint engine.  For an interval [L,R) the
endpoint-numerator summand is

    root^(w_B * 169 * L) - root^(w_B * 169 * R),

where w_B=13+U_full[TARGET_B].  The script checks this direct interval route
against fast_endpoint_sum, then carries the rebuilt values through every
alpha, beta, tau guard, refined residue t, and retained (ell,c) coordinate.

This is an exact finite Fubini certificate on an external product: the
source ancestry variable and the BY endpoint variable remain independent.
It is not a same-stalk temporal transport, a THM-2449 response table, a
THM-2512 bridge, a physical current, or an LRC(14) conclusion.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stdout
from hashlib import sha256
import importlib.util
import io
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
P = 13
Q = 7

SPECTRUM_PATH = (
    ROOT
    / "04-computation/lrc_r5_source_aligned_relation_residue_7x13_spectrum_probe_20260816.py"
)
SPECTRUM_SHA256 = "5f3fbf08bef6f9a61e684f0f7616e80e1dbbda4f6bb2ed4ca3788d3b8b53d65a"
SPECTRUM_SEMANTIC = "cd55336bb1dfe5f37f020c242c4bca5b7c6be339ec57e95d69e10bbe68d9dbaa"
TARGET_PATH = (
    ROOT
    / "04-computation/lrc_ufull_owner_boundary_k4xf13_endpoint_factorization_independent_audit_20260816.py"
)
TARGET_SHA256 = "f89be10c65bb77270199f9399b155d5a2c82c0da121b3e8589fe3c1f7e9824fc"
TARGET_SEMANTIC = "d52c9f0a56c14a83e1e6b175c7b725314c99f09d44509bc8582847a5857f7da6"

EXPECTED_RAW_COUNTS = (
    531398, 541349, 556004, 554331, 546963, 550023, 539135,
    539135, 550023, 546963, 554331, 556004, 541349,
)
EXPECTED_SPLIT_COUNTS = (
    531398, 541470, 556125, 554452, 547084, 550144, 539256,
    539256, 550144, 547084, 554452, 556125, 541470,
)
EXPECTED_PAIR_BANK_SHA256 = "c28119c8b54f47e5b7a46f1508fbba604b0e3997eaadb05b03ad28edd9aed468"
EXPECTED_TENSOR_SHA256 = "39d7a0b4e5b2d8b85631d682ed1967091e44dc41e17b33a77e7184d3dc93e0cf"
EXPECTED_COORDINATE_BANK_SHA256 = "989dafc220a6d09aeacfce4af0e9a4fe13eedacc79fa66032ea39bc107fd8efb"
EXPECTED_SPECTRUM_BANK_SHA256 = "5f173227c5e203309f61bdfd9d47cc64a3b49ae8f14abd0f7bfc469eda278533"
EXPECTED_SEMANTIC_SHA256 = "6b87dc762517773254e3388dffd0bed9a5101685730dd3b93ebf4c639f224e49"

_WORKER_TARGET = None


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


def literal_worker(alpha: int) -> tuple[object, ...]:
    """Return AX, preintegrated BY, and independently summed BY intervals."""
    global _WORKER_TARGET
    if _WORKER_TARGET is None:
        _WORKER_TARGET = load_module(
            TARGET_PATH, "one_leg_fubini_target_worker", TARGET_SHA256
        )
    target = _WORKER_TARGET
    (
        word,
        t_den,
        nn,
        prime,
        root,
        _zeta,
        _q_intervals,
        _q_starts,
        _embeddings,
        _tabs,
        atom_intervals,
    ) = target.context()
    no_guard = dict(target.M.PATTERN_E)
    require(no_guard.pop(target.M.GUARD) == "guard_safe", "removed non-guard")
    y_frequency = target.M.X_FREQUENCY + word[target.M.TARGET_B]

    rows = []
    raw_count = 0
    split_count = 0
    nonempty_count = 0
    total_measure = 0
    for beta in range(P):
        ell = (0, -alpha % P, -beta % P, 0, 0, 0, 0, alpha, beta)
        unguarded = target.M.fast_build_set(word, t_den, no_guard, ell)
        groups = target.partition_two_pointer(unguarded, atom_intervals)
        ax_values, by_preintegrated, _overlap = target.atom_endpoint_tables(groups)

        by_literal = []
        atom_measures = []
        atom_piece_counts = []
        for group in groups:
            literal_total = 0
            measure = 0
            for left, right in group:
                exponent_left = y_frequency * target.M.R_DILATION * left % nn
                exponent_right = y_frequency * target.M.R_DILATION * right % nn
                literal_total = (
                    literal_total
                    + pow(root, exponent_left, prime)
                    - pow(root, exponent_right, prime)
                ) % prime
                measure += right - left
            by_literal.append(literal_total)
            atom_measures.append(measure % prime)
            atom_piece_counts.append(len(group))
            total_measure += measure

        by_literal_tuple = tuple(by_literal)
        require(
            by_literal_tuple == by_preintegrated,
            ("literal BY versus fast endpoint sum", alpha, beta),
        )
        require(len(ax_values) == len(by_literal_tuple) == len(target.ATOMS),
                ("atom widths", alpha, beta))
        rows.append(
            (
                beta,
                ax_values,
                by_preintegrated,
                by_literal_tuple,
                tuple(atom_measures),
                tuple(atom_piece_counts),
            )
        )
        raw_count += len(unguarded)
        split_count += sum(len(group) for group in groups)
        nonempty_count += sum(bool(group) for group in groups)

    return alpha, tuple(rows), raw_count, split_count, nonempty_count, total_measure


def build_guard_kernels(target, zeta: int, prime: int, drop_right: bool = False):
    kernels = [
        [[0 for _right in target.ATOMS] for _left in target.ATOMS]
        for _t in range(P)
    ]
    for left_index, (left_sheet, left_chamber) in enumerate(target.ATOMS):
        for right_index, (right_sheet, right_chamber) in enumerate(target.ATOMS):
            products = tuple(
                target.safe(left_chamber, left_sheet + tau)
                * (
                    1
                    if drop_right
                    else target.safe(right_chamber, right_sheet + tau)
                )
                for tau in range(P)
            )
            for t in range(P):
                kernels[t][left_index][right_index] = sum(
                    value * pow(zeta, -t * tau % P, prime)
                    for tau, value in enumerate(products)
                ) % prime
    return kernels


def build_pair_bank(
    tables: tuple[tuple[object, ...], ...],
    kernels,
    zeta: int,
    prime: int,
    value_slot: int,
) -> tuple[tuple[tuple[int, ...], ...], ...]:
    size = len(tables[0][2])
    bank = [[[0 for _right in range(size)] for _left in range(size)] for _t in range(P)]
    for alpha, beta, ax_values, by_values, _literal, _measures, _counts in tables:
        phase = pow(zeta, (beta - alpha) % P, prime)
        selected = (
            by_values
            if value_slot == 2
            else _literal
            if value_slot == 3
            else _measures
        )
        for left_index, left_value in enumerate(ax_values):
            if not left_value:
                continue
            left_phase = phase * left_value % prime
            for right_index, right_value in enumerate(selected):
                if not right_value:
                    continue
                base = left_phase * right_value % prime
                for t in range(P):
                    bank[t][left_index][right_index] = (
                        bank[t][left_index][right_index]
                        + base * kernels[t][left_index][right_index]
                    ) % prime
    normalizer = pow(P**3, -1, prime)
    return tuple(
        tuple(
            tuple(value * normalizer % prime for value in row)
            for row in matrix
        )
        for matrix in bank
    )


def mismatch_count(left, right) -> int:
    require(len(left) == len(right), "mismatch outer shape")
    return sum(
        x != y
        for left_matrix, right_matrix in zip(left, right)
        for left_row, right_row in zip(left_matrix, right_matrix)
        for x, y in zip(left_row, right_row)
    )


def main() -> None:
    spectrum = load_module(
        SPECTRUM_PATH, "one_leg_fubini_spectrum", SPECTRUM_SHA256
    )
    require(
        spectrum.EXPECTED_SEMANTIC_SHA256 == SPECTRUM_SEMANTIC,
        "spectrum semantic drift",
    )
    endpoint = spectrum.load_endpoint_module()
    target = endpoint.load_module(
        TARGET_PATH, "one_leg_fubini_target", TARGET_SHA256
    )
    require(target.EXPECTED_SEMANTIC_SHA256 == TARGET_SEMANTIC,
            "target semantic drift")

    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(literal_worker, range(P)))
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)), "worker order")
    require(tuple(chunk[2] for chunk in chunks) == EXPECTED_RAW_COUNTS,
            "raw interval census")
    require(tuple(chunk[3] for chunk in chunks) == EXPECTED_SPLIT_COUNTS,
            "split interval census")
    tables = tuple(
        (alpha, *row)
        for alpha, rows, _raw, _split, _nonempty, _measure in chunks
        for row in rows
    )
    require(len(tables) == P**2, "endpoint table count")
    atom_entry_count = len(tables) * len(target.ATOMS)
    literal_piece_count = sum(
        sum(counts) for _a, _b, _ax, _pre, _literal, _measure, counts in tables
    )
    require(literal_piece_count == sum(EXPECTED_SPLIT_COUNTS),
            "literal piece census")
    by_equality_count = sum(
        pre == literal
        for _a, _b, _ax, pre, literal, _measure, _counts in tables
    ) * len(target.ATOMS)
    require(by_equality_count == atom_entry_count, "BY Fubini entry census")

    (
        _word,
        _t_den,
        nn,
        prime,
        root,
        zeta,
        *_rest,
    ) = target.context()
    require(pow(zeta, P, prime) == 1 and zeta != 1, "order-thirteen root")
    eta = pow(root, nn // Q, prime)
    require(pow(eta, Q, prime) == 1 and eta != 1, "order-seven root")

    guarded_kernels = build_guard_kernels(target, zeta, prime)
    right_guard_dropped_kernels = build_guard_kernels(
        target, zeta, prime, drop_right=True
    )
    preintegrated_bank = build_pair_bank(
        tables, guarded_kernels, zeta, prime, value_slot=2
    )
    literal_bank = build_pair_bank(
        tables, guarded_kernels, zeta, prime, value_slot=3
    )
    require(literal_bank == preintegrated_bank, "pair-bank Fubini equality")
    pair_bank_digest = digest_json(literal_bank)
    require(pair_bank_digest == EXPECTED_PAIR_BANK_SHA256,
            ("pinned endpoint bank", pair_bank_digest))

    source = endpoint.load_module(
        spectrum.SOURCE_PATH,
        "one_leg_fubini_source",
        spectrum.SOURCE_SHA256,
    )
    with redirect_stdout(io.StringIO()):
        source_data = source.main()
    require(source_data["semantic_sha256"] == spectrum.SOURCE_SEMANTIC,
            "source semantic drift")
    tensor, _marginal, cell_record, pair_support, entry_support = (
        spectrum.build_cell_tensor(source, source_data)
    )
    tensor_digest = digest_json(tensor)
    require(tensor_digest == EXPECTED_TENSOR_SHA256, ("source tensor", tensor_digest))
    denominator = source_data["denominator"]
    require(denominator % prime != 0, "source denominator at split prime")

    coordinate_bank = tuple(
        spectrum.contract(tensor, matrix, denominator, prime)
        for matrix in literal_bank
    )
    spectrum_bank = tuple(
        spectrum.fourier_2d(matrix, eta, zeta, prime)
        for matrix in coordinate_bank
    )
    coordinate_digest = digest_json(coordinate_bank)
    spectrum_digest = digest_json(spectrum_bank)
    require(coordinate_digest == EXPECTED_COORDINATE_BANK_SHA256,
            ("coordinate Fubini equality", coordinate_digest))
    require(spectrum_digest == EXPECTED_SPECTRUM_BANK_SHA256,
            ("spectrum Fubini equality", spectrum_digest))
    shapes = tuple(spectrum.support_shape(bank) for bank in spectrum_bank)
    require(shapes == ((91, 1, 6, 12, 72),) * P, "full mixed support")

    relation_t = 6
    relation_spectrum = spectrum_bank[relation_t]
    require(relation_spectrum[1][1] == 218019411785559321795,
            "fixed-relation mixed witness")

    # Hostile 1: delete the right-leg guard while retaining all other sums.
    right_guard_dropped_bank = build_pair_bank(
        tables, right_guard_dropped_kernels, zeta, prime, value_slot=3
    )
    right_guard_mismatches = mismatch_count(literal_bank, right_guard_dropped_bank)
    require(right_guard_mismatches > 0, "right-guard deletion was invisible")
    right_guard_relation = spectrum.contract(
        tensor, right_guard_dropped_bank[relation_t], denominator, prime
    )
    require(right_guard_relation != coordinate_bank[relation_t],
            "right-guard hostile reached the lawful relation table")
    right_guard_shape = spectrum.support_shape(
        spectrum.fourier_2d(right_guard_relation, eta, zeta, prime)
    )

    # Hostile 2: retain interval measures but erase the BY character phase.
    measure_only_bank = build_pair_bank(
        tables, guarded_kernels, zeta, prime, value_slot=4
    )
    measure_mismatches = mismatch_count(literal_bank, measure_only_bank)
    require(measure_mismatches > 0, "measure-only BY hostile was invisible")
    measure_relation = spectrum.contract(
        tensor, measure_only_bank[relation_t], denominator, prime
    )
    require(measure_relation != coordinate_bank[relation_t],
            "measure-only hostile reached the lawful relation table")
    measure_shape = spectrum.support_shape(
        spectrum.fourier_2d(measure_relation, eta, zeta, prime)
    )

    # Hostiles 3/4: marginalize either retained common-base coordinate.
    inv_q = pow(Q, -1, prime)
    cell_erased = tuple(
        tuple(
            sum(coordinate_bank[relation_t][ell][offset] for ell in range(Q))
            * inv_q
            % prime
            for offset in range(P)
        )
        for _ell in range(Q)
    )
    inv_p = pow(P, -1, prime)
    owner_erased = tuple(
        tuple(
            sum(coordinate_bank[relation_t][ell]) * inv_p % prime
            for _offset in range(P)
        )
        for ell in range(Q)
    )
    cell_erased_shape = spectrum.support_shape(
        spectrum.fourier_2d(cell_erased, eta, zeta, prime)
    )
    owner_erased_shape = spectrum.support_shape(
        spectrum.fourier_2d(owner_erased, eta, zeta, prime)
    )
    require(cell_erased_shape == (13, 1, 0, 12, 0), "cell erasure hostile")
    require(owner_erased_shape == (7, 1, 6, 0, 0), "owner erasure hostile")

    table_record = tuple(
        (
            alpha,
            beta,
            ax,
            pre,
            literal,
            measures,
            counts,
        )
        for alpha, beta, ax, pre, literal, measures, counts in tables
    )
    record = (
        SPECTRUM_SHA256,
        SPECTRUM_SEMANTIC,
        TARGET_SHA256,
        TARGET_SEMANTIC,
        tuple(chunk[2:] for chunk in chunks),
        digest_json(table_record),
        atom_entry_count,
        literal_piece_count,
        by_equality_count,
        digest_json(guarded_kernels),
        pair_bank_digest,
        tensor_digest,
        pair_support,
        entry_support,
        cell_record,
        coordinate_digest,
        spectrum_digest,
        shapes,
        relation_t,
        relation_spectrum[1][1],
        right_guard_mismatches,
        digest_json(right_guard_dropped_bank),
        right_guard_shape,
        measure_mismatches,
        digest_json(measure_only_bank),
        measure_shape,
        cell_erased_shape,
        owner_erased_shape,
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== r=5 one-leg literal-BY integrand Fubini probe ==")
    print(
        "status=FINITE-EXACT ONE-LEG EXTERNAL-PRODUCT FUBINI; "
        "not same-stalk temporal transport; LRC(14) OPEN"
    )
    print(
        f"dependencies=((spectrum,{SPECTRUM_SHA256},{SPECTRUM_SEMANTIC}),"
        f"(target,{TARGET_SHA256},{TARGET_SEMANTIC}))"
    )
    print(f"split_embedding=(prime={prime},zeta7={eta},zeta13={zeta})")
    print(
        f"literal_BY_interval_census=(pieces={literal_piece_count},"
        f"atom_entries={atom_entry_count},exact_equalities={by_equality_count})"
    )
    print(
        "literal_BY_formula=root^(w_B*169*L)-root^(w_B*169*R) "
        "for each retained half-open [L,R)"
    )
    print("independent_routes=(direct endpoint primitive,fast_endpoint_sum): PASS")
    print(
        f"guarded_pair_bank_Fubini=(entries={P*len(target.ATOMS)**2},"
        f"sha256={pair_bank_digest}): PASS"
    )
    print(
        f"source_common_ancestry=(pair_support={pair_support}/1521,"
        f"entry_support={entry_support}/138411,tensor_sha256={tensor_digest})"
    )
    print(
        f"retained_output_Fubini=(coordinates={P*Q*P},"
        f"coordinate_sha256={coordinate_digest},spectrum_sha256={spectrum_digest}): PASS"
    )
    print("all_13_residue_spectrum_shapes=(91,1,6,12,72)")
    print(
        f"fixed_relation=(1,0,{relation_t});mixed_witness_(1,1)="
        f"{relation_spectrum[1][1]}"
    )
    print(
        f"hostile_right_guard_deleted=(pair_mismatches={right_guard_mismatches},"
        f"relation_shape={right_guard_shape},sha256={digest_json(right_guard_dropped_bank)})"
    )
    print(
        f"hostile_BY_phase_erased=(pair_mismatches={measure_mismatches},"
        f"relation_shape={measure_shape},sha256={digest_json(measure_only_bank)})"
    )
    print(
        f"hostile_common_base_marginals=(cell_erased={cell_erased_shape},"
        f"owner_erased={owner_erased_shape})"
    )
    print(
        "variables=source ancestry y is integrated inside M at fixed (ell,c); "
        "BY endpoint v is separately interval-integrated; AX endpoint x/Q is still "
        "preintegrated; alpha,beta,tau,t are finite sums"
    )
    print(
        "typing=lawful exact iterated/product-space Fubini with every guard and "
        "source coordinate retained; NOT a diagonal/common-stalk temporal map"
    )
    print(
        "nonconsequence=no THM-2449/2512 row,grouped exact address,U_full/U_clock "
        "same-tuple identification,current,row exclusion,or LRC(14) closure"
    )
    print(f"semantic_sha256={semantic}")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
