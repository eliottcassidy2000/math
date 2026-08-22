#!/usr/bin/env python3
"""Exact all-role spectral atlas on the U_full guard-bucket tensor.

The parent guard-bucket theorem decomposes two distinguished endpoint
responses into four active boundary-chamber pairs and thirteen sheet drifts.
This companion retains the same 39-atom construction but evaluates all five
distinct refined role classes.  It then tests the 72 lawful two-K4-and-bridge
role charts in four separately typed representations:

* each of the 52 active chamber-pair/drift buckets;
* the 13 drift marginals after summing the four chamber pairs;
* the 13 Fourier modes of those drift marginals; and
* all 4 x 13 chamber-pair Fourier fibres.

Every edge response in these banks is still obtained from vertex potentials,
so it is a graph coboundary.  The calculation is an endpoint-response and
weighted-tree-spectrum theorem component, not a physical current, ancestry
coupling, absolute H^1 flux, LRC bispectrum, scalar exclusion, or LRC(14)
conclusion.
"""

from __future__ import annotations

import ast
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import hashlib
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "04-computation/lrc_ufull_guard_bucket_all_role_spectral_probe_20260816.py"
BUCKET_PATH = ROOT / "04-computation/lrc_ufull_guard_sheet_drift_bucket_bridge_probe_20260816.py"
ROLE_PATH = ROOT / "04-computation/lrc14_guard_deleted_refined_endpoint_role_probe_20260816.py"
BUCKET_SHA256 = "a1d4b667812949001fc863ba881ff7409bbae3c568a6bf7bc24c9dc88b2766b1"
ROLE_SHA256 = "ee2105742abee578a9c41ff7ec954a07ada324fccc2c643429e7ac6e6e6f8fc2"
EXPECTED_GAMMA_SHA256 = "1fabc5cfdbaa1455e10cd6bf9264488133616a7b0ff381623d729b4b4bfa9682"
EXPECTED_VALUES = {
    (0, 0, 0): 405336876493642499425,
    (0, 1, 0): 518539850465495448196,
    (1, 0, 0): 503604956476841920373,
    (1, 0, 1): 320618948602619577408,
    (1, 12, 0): 15703541686881447885,
}
EXPECTED_SEMANTIC_SHA256 = "2c9495fb8bcb731361ba331d9ca4b84a60f21551dc49b16e6519c1fc4f2e9f97"
P = 13
Q_CLASSES = tuple(sorted(EXPECTED_VALUES))
CORNER_PAIRS = (
    ("left", "left"),
    ("left", "right"),
    ("right", "left"),
    ("right", "right"),
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_integers(values: tuple[int, ...]) -> str:
    return hashlib.sha256(
        ",".join(str(value) for value in values).encode("ascii")
    ).hexdigest()


def digest_json(value: object) -> str:
    body = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return hashlib.sha256(body).hexdigest()


def load_module(path: Path, name: str, expected_hash: str):
    actual_hash = lf_sha256(path)
    if expected_hash != "TO_BE_PINNED":
        require(actual_hash == expected_hash, (name, actual_hash, expected_hash))
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, (name, "loader"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


B = load_module(BUCKET_PATH, "lrc_guard_bucket_parent", BUCKET_SHA256)
R = load_module(ROLE_PATH, "lrc_refined_role_parent", ROLE_SHA256)


def worker(alpha: int) -> tuple[object, ...]:
    (
        word,
        t_den,
        _nn,
        prime,
        _root,
        zeta,
        _q_intervals,
        _q_starts,
        _embeddings,
        _tabs,
        atom_intervals,
        atom_starts,
        tau_h,
        tau_zero,
    ) = B.context()
    no_guard_pattern = dict(B.M.PATTERN_E)
    require(no_guard_pattern.pop(B.M.GUARD) == "guard_safe", "guard factor absent")
    bucket_rows = [[0] * len(B.BUCKETS) for _q in Q_CLASSES]
    gamma_rows: list[int] = []
    raw_count = 0
    split_count = 0
    nonempty_count = 0
    overlap_length = 0
    for beta in range(P):
        ell = (0, -alpha % P, -beta % P, 0, 0, 0, 0, alpha, beta)
        e_unguarded = B.M.fast_build_set(word, t_den, no_guard_pattern, ell)
        groups = B.split_by_guard_atom(e_unguarded, atom_intervals, atom_starts)
        require(
            all(
                not groups[index]
                for index, (_sheet, chamber) in enumerate(B.ATOMS)
                if chamber == "middle"
            ),
            ("middle support", alpha, beta),
        )
        ax_values, by_values, overlap = B.atom_endpoint_values(groups)
        raw_count += len(e_unguarded)
        split_count += sum(len(group) for group in groups)
        nonempty_count += sum(bool(group) for group in groups)
        overlap_length += overlap
        phase = pow(zeta, beta, prime)

        for tau in range(P):
            ax = sum(
                value * B.safe(chamber, sheet + tau)
                for value, (sheet, chamber) in zip(ax_values, B.ATOMS)
            ) % prime
            by = sum(
                value * B.safe(chamber, sheet + tau)
                for value, (sheet, chamber) in zip(by_values, B.ATOMS)
            ) % prime
            gamma_rows.append(phase * ax % prime * by % prime)

        characters = tuple(
            phase
            * pow(zeta, -(alpha * q[0] + beta * q[1]) % P, prime)
            % prime
            for q in Q_CLASSES
        )
        for left_index, (left_sheet, left_chamber) in enumerate(B.ATOMS):
            left_value = ax_values[left_index]
            if left_value == 0:
                continue
            for right_index, (right_sheet, right_chamber) in enumerate(B.ATOMS):
                right_value = by_values[right_index]
                if right_value == 0:
                    continue
                drift = (right_sheet - left_sheet) % P
                bucket_index = B.BUCKET_INDEX[(left_chamber, right_chamber, drift)]
                product = left_value * right_value % prime
                for q_index, q in enumerate(Q_CLASSES):
                    require(q[2] in (0, 1), ("tau character", q))
                    kernel = (
                        tau_h[left_index][right_index]
                        if q[2] == 1
                        else tau_zero[left_index][right_index]
                    )
                    bucket_rows[q_index][bucket_index] = (
                        bucket_rows[q_index][bucket_index]
                        + characters[q_index] * product % prime * kernel
                    ) % prime
    return (
        alpha,
        raw_count,
        split_count,
        nonempty_count,
        overlap_length,
        tuple(gamma_rows),
        tuple(tuple(row) for row in bucket_rows),
    )


def role_values(q_values: dict[tuple[int, int, int], int]) -> dict[str, int]:
    return {label: q_values[q] for label, q in R.ROLE_CLASSES.items()}


def chart_bank(
    entries: tuple[tuple[object, dict[tuple[int, int, int], int]], ...],
    prime: int,
) -> tuple[object, ...]:
    chart_order = R.role_charts()
    factor_zero_totals = [0, 0, 0, 0]
    product_zeros_by_chart = [0] * len(chart_order)
    profiles: Counter[tuple[int, int, int, int]] = Counter()
    payload = []
    all_nonzero_entries = 0
    for key, q_values in entries:
        rows = R.chart_census(role_values(q_values), prime)
        require(tuple(row[0] for row in rows) == chart_order, ("chart order", key))
        factors = tuple(tuple(row[1:]) for row in rows)
        zero_profile = tuple(
            sum(factors[chart][factor] == 0 for chart in range(len(factors)))
            for factor in range(4)
        )
        for factor, count in enumerate(zero_profile):
            factor_zero_totals[factor] += count
        for chart, row in enumerate(factors):
            product_zeros_by_chart[chart] += row[3] == 0
        profiles[zero_profile] += 1
        all_nonzero_entries += zero_profile == (0, 0, 0, 0)
        payload.append((key, factors))
    return (
        len(entries),
        tuple(factor_zero_totals),
        all_nonzero_entries,
        tuple(sorted(profiles.items())),
        tuple(product_zeros_by_chart),
        digest_json(tuple(payload)),
    )


def dft(values: tuple[int, ...], frequency: int, zeta: int, prime: int) -> int:
    return sum(
        value * pow(zeta, -frequency * drift % P, prime)
        for drift, value in enumerate(values)
    ) % prime


def security_certificate(path: Path) -> tuple[int, tuple[str, ...]]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    bad = tuple(type(node).__name__ for node in ast.walk(tree) if isinstance(node, ast.Assert))
    require(not bad, ("security", bad))
    return len(tuple(ast.walk(tree))), bad


def main() -> None:
    require(lf_sha256(ROLE_PATH) == ROLE_SHA256, "role source hash drift")
    if BUCKET_SHA256 != "TO_BE_PINNED":
        require(lf_sha256(BUCKET_PATH) == BUCKET_SHA256, "bucket source hash drift")
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(worker, range(P)))
    require(tuple(row[0] for row in chunks) == tuple(range(P)), "worker order")
    gamma = tuple(value for row in chunks for value in row[5])
    require(len(gamma) == P**3, ("gamma length", len(gamma)))
    gamma_hash = digest_integers(gamma)
    require(gamma_hash == EXPECTED_GAMMA_SHA256, ("gamma hash", gamma_hash))

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
        _atom_starts,
        _tau_h,
        _tau_zero,
    ) = B.context()
    direct_values = {
        q: B.inverse_value(gamma, q, prime, zeta)
        for q in Q_CLASSES
    }
    require(direct_values == EXPECTED_VALUES, ("role values", direct_values))
    normalizer = pow(P**3, -1, prime)
    q_buckets = {
        q: tuple(
            sum(chunk[6][q_index][bucket] for chunk in chunks)
            % prime
            * normalizer
            % prime
            for bucket in range(len(B.BUCKETS))
        )
        for q_index, q in enumerate(Q_CLASSES)
    }
    for q in Q_CLASSES:
        require(sum(q_buckets[q]) % prime == direct_values[q], ("bucket sum", q))

    active_indices = tuple(
        index for index, bucket in enumerate(B.BUCKETS)
        if "middle" not in bucket[:2]
    )
    inactive_indices = tuple(
        index for index, bucket in enumerate(B.BUCKETS)
        if "middle" in bucket[:2]
    )
    require(len(active_indices) == 4 * P and len(inactive_indices) == 5 * P,
            (len(active_indices), len(inactive_indices)))
    require(
        all(q_buckets[q][index] == 0 for q in Q_CLASSES for index in inactive_indices),
        "middle bucket acquired support",
    )

    support_counts = tuple(
        (q, sum(q_buckets[q][index] != 0 for index in active_indices))
        for q in Q_CLASSES
    )
    tensor_rank = B.rank_mod(
        tuple(tuple(q_buckets[q][index] for index in active_indices) for q in Q_CLASSES),
        prime,
    )
    point_entries = tuple(
        (
            B.BUCKETS[index],
            {q: q_buckets[q][index] for q in Q_CLASSES},
        )
        for index in active_indices
    )
    point_bank = chart_bank(point_entries, prime)

    drift_values = {
        q: tuple(
            sum(
                q_buckets[q][B.BUCKET_INDEX[(left, right, drift)]]
                for left, right in CORNER_PAIRS
            ) % prime
            for drift in range(P)
        )
        for q in Q_CLASSES
    }
    for q in Q_CLASSES:
        require(sum(drift_values[q]) % prime == direct_values[q], ("drift sum", q))
    drift_entries = tuple(
        (drift, {q: drift_values[q][drift] for q in Q_CLASSES})
        for drift in range(P)
    )
    drift_bank = chart_bank(drift_entries, prime)

    spectral_values = {
        q: tuple(dft(drift_values[q], frequency, zeta, prime) for frequency in range(P))
        for q in Q_CLASSES
    }
    require(
        all(spectral_values[q][0] == direct_values[q] for q in Q_CLASSES),
        "zero mode does not recover global role value",
    )
    spectral_entries = tuple(
        (frequency, {q: spectral_values[q][frequency] for q in Q_CLASSES})
        for frequency in range(P)
    )
    spectral_bank = chart_bank(spectral_entries, prime)

    chamber_spectral_entries = []
    drift_slice_ranks = []
    spectral_slice_ranks = []
    for drift in range(P):
        drift_slice_ranks.append(
            B.rank_mod(
                tuple(
                    tuple(
                        q_buckets[q][B.BUCKET_INDEX[(left, right, drift)]]
                        for left, right in CORNER_PAIRS
                    )
                    for q in Q_CLASSES
                ),
                prime,
            )
        )
    for frequency in range(P):
        matrix = []
        for q in Q_CLASSES:
            row = []
            for left, right in CORNER_PAIRS:
                values = tuple(
                    q_buckets[q][B.BUCKET_INDEX[(left, right, drift)]]
                    for drift in range(P)
                )
                row.append(dft(values, frequency, zeta, prime))
            matrix.append(tuple(row))
        spectral_slice_ranks.append(B.rank_mod(tuple(matrix), prime))
        for corner_index, (left, right) in enumerate(CORNER_PAIRS):
            chamber_spectral_entries.append(
                (
                    (left, right, frequency),
                    {q: matrix[q_index][corner_index] for q_index, q in enumerate(Q_CLASSES)},
                )
            )
    chamber_spectral_bank = chart_bank(tuple(chamber_spectral_entries), prime)

    global_bank = chart_bank((("global", direct_values),), prime)
    require(global_bank[1] == (0, 0, 0, 0), ("global chart gate", global_bank))
    flat = R.chart_census({label: 1 for label in R.ROLE_CLASSES}, prime)
    require(all(row[1:] == (0, 0, 0, 0) for row in flat), "flat hostile")
    killed_values = role_values(direct_values)
    killed_values["q5"] = killed_values["H"]
    killed = R.chart_census(killed_values, prime)
    require(all(row[1] == 0 and row[4] == 0 for row in killed), "bridge hostile")

    security = security_certificate(ROOT / SCRIPT)
    semantic = (
        lf_sha256(BUCKET_PATH),
        ROLE_SHA256,
        (prime, zeta),
        tuple((q, direct_values[q]) for q in Q_CLASSES),
        tuple((q, q_buckets[q]) for q in Q_CLASSES),
        support_counts,
        tensor_rank,
        point_bank,
        tuple((q, drift_values[q]) for q in Q_CLASSES),
        drift_bank,
        tuple((q, spectral_values[q]) for q in Q_CLASSES),
        spectral_bank,
        tuple(drift_slice_ranks),
        tuple(spectral_slice_ranks),
        chamber_spectral_bank,
        global_bank,
        security,
        "endpoint role potentials only; graph coboundaries; no ancestry/current/H1/LRC14",
    )
    semantic_hash = hashlib.sha256(repr(semantic).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256,
                ("semantic hash", semantic_hash))

    print("LRC U_FULL ALL-ROLE GUARD-BUCKET SPECTRAL ATLAS")
    print("status=FINITE-EXACT endpoint-response theorem component; weighted graph factors; LRC(14) OPEN")
    print(f"dependency_hashes=(bucket={lf_sha256(BUCKET_PATH)},role={ROLE_SHA256})")
    print(f"field=(prime={prime},zeta13={zeta}); gamma_sha256={gamma_hash}")
    print(f"role_classes={tuple((label,R.ROLE_CLASSES[label]) for label in R.ROLE_CLASSES)}")
    print(f"distinct_q_values={tuple((q,direct_values[q]) for q in Q_CLASSES)}")
    print(f"tensor_shape=(q_classes={len(Q_CLASSES)},chamber_pairs={len(CORNER_PAIRS)},drifts={P}); active_entries_per_q={len(active_indices)}")
    print(f"active_bucket_support_counts={support_counts}; five_by_52_rank={tensor_rank}")
    print(f"pointwise_52x72_chart_bank=(entries,factor_zero_totals,all_nonzero_entries,profiles,product_zeros_by_chart,digest)={point_bank}")
    print(f"drift_marginal_13x72_chart_bank={drift_bank}")
    print(f"drift_fourier_13x72_chart_bank={spectral_bank}")
    print(f"chamber_fourier_52x72_chart_bank={chamber_spectral_bank}")
    print(f"slice_ranks=(drift={tuple(drift_slice_ranks)},fourier={tuple(spectral_slice_ranks)})")
    print(f"global_72_chart_bank={global_bank}")
    print("hostiles=(flat role response kills all factors; setting q5=H kills the forced bridge and every product)")
    print("typing=every bank is formed from eight role potentials and hence lands in graph B^1; no absolute H1 flux is inferred")
    print("next_test=insert one lawful character-independent common-stalk ancestry relation before endpoint marginalization, then repeat the same tensor and chart gates")
    print(f"security_ast={security}")
    print(f"semantic_sha256={semantic_hash}")
    print("nonconsequence=no grouped current, all-unit projector, physical chronology, bispectrum theorem, scalar-row exclusion, or LRC(14)")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
