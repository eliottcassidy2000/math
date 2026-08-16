#!/usr/bin/env python3
"""Independent audit of all-role U_full endpoint weighted-tree closure.

This companion imports the independently audited THM-3514 atom-table engine,
not the all-role candidate.  It reconstructs the five refined role-class
bucket rows, writes its own 72-chart generator and matrix-tree evaluator, and
tests pointwise, drift-marginal, drift-Fourier, and chamber-Fourier banks.

Every graph edge remains a difference of vertex potentials, hence a
coboundary.  No ancestry, physical current, absolute H^1 class, grouped
coefficient, scalar-row exclusion, or LRC(14) conclusion is asserted.
"""

from __future__ import annotations

import ast
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import hashlib
import importlib.util
from itertools import permutations
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    "04-computation/"
    "lrc_ufull_all_role_weighted_tree_spectral_closure_"
    "independent_audit_20260816.py"
)
OUTPUT = (
    "05-knowledge/results/"
    "lrc_ufull_all_role_weighted_tree_spectral_closure_"
    "independent_audit_20260816.out"
)
ATOM_AUDIT_PATH = ROOT / (
    "04-computation/"
    "lrc_ufull_owner_boundary_k4xf13_endpoint_factorization_"
    "independent_audit_20260816.py"
)
ATOM_AUDIT_SHA256 = "f89be10c65bb77270199f9399b155d5a2c82c0da121b3e8589fe3c1f7e9824fc"
CANDIDATE_SEMANTIC_SHA256 = "2c9495fb8bcb731361ba331d9ca4b84a60f21551dc49b16e6519c1fc4f2e9f97"
EXPECTED_SEMANTIC_SHA256 = "61db05458acb26cd8d6a9a89991d3bbeed91f95f6bdf035720a3a93a5b4f886b"

P = 13
Q_VALUES = {
    (0, 0, 0): 405336876493642499425,
    (0, 1, 0): 518539850465495448196,
    (1, 0, 0): 503604956476841920373,
    (1, 0, 1): 320618948602619577408,
    (1, 12, 0): 15703541686881447885,
}
Q_CLASSES = tuple(sorted(Q_VALUES))
ROLE_CLASSES = {
    "c1": (0, 0, 0),
    "c2": (1, 0, 0),
    "c3": (0, 1, 0),
    "H": (1, 0, 1),
    "q2": (1, 12, 0),
    "q3": (1, 0, 0),
    "q4": (1, 0, 0),
    "q5": (1, 0, 0),
}
EDGES = (
    (0, 3), (0, 4), (0, 5), (1, 2), (1, 4), (1, 7), (2, 4),
    (2, 7), (3, 4), (3, 5), (4, 5), (4, 6), (4, 7),
)
HUB = 4
LEAF = 6
WINGS = ((0, 3, 5), (1, 2, 7))
BLOCKERS = ("c1", "c2", "c3")
UNITS = ("q2", "q3", "q4")
CORNERS = (
    ("left", "left"),
    ("left", "right"),
    ("right", "left"),
    ("right", "right"),
)
EXPECTED_BANK_DIGESTS = {
    "point": "a1d5183f37cc39deed976876ed91132f662f91137f41a10a79fd3e974ced2dfc",
    "drift": "6bcdd9fa616ba2dc4883d8334dca7358d06375d7342f43c4aac8dea45deb5027",
    "fourier": "e96a3859e6838fa91d04a55ca55fa4e44cafa45bd86b93c7f6b612b1e65f4dc1",
    "chamber_fourier": "c9ca8c9b7ed7000b00e39b496523ad75a44379d84ea35e49163d6b64e4196e73",
}


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
    return hashlib.sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def load_atom_audit():
    require(lf_sha256(ATOM_AUDIT_PATH) == ATOM_AUDIT_SHA256, "THM-3514 audit hash drift")
    spec = importlib.util.spec_from_file_location("thm3514_atom_tables", ATOM_AUDIT_PATH)
    require(spec is not None and spec.loader is not None, "atom audit loader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


A = load_atom_audit()


def worker(alpha: int) -> tuple[object, ...]:
    return A.worker(alpha)


def determinant_mod(matrix: tuple[tuple[int, ...], ...], prime: int) -> int:
    rows = [list(value % prime for value in row) for row in matrix]
    answer = 1
    for column in range(len(rows)):
        pivot = next(
            (row for row in range(column, len(rows)) if rows[row][column]),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            answer = -answer
        pivot_value = rows[column][column] % prime
        answer = answer * pivot_value % prime
        inverse = pow(pivot_value, -1, prime)
        for row in range(column + 1, len(rows)):
            factor = rows[row][column] * inverse % prime
            if factor:
                rows[row] = [
                    (value - factor * pivot_entry) % prime
                    for value, pivot_entry in zip(rows[row], rows[column])
                ]
    return answer % prime


def role_charts() -> tuple[tuple[str, ...], ...]:
    charts = []
    for swap in (0, 1):
        for blockers in permutations(BLOCKERS):
            for units in permutations(UNITS):
                chart = {HUB: "H", LEAF: "q5"}
                chart.update(zip(WINGS[swap], blockers))
                chart.update(zip(WINGS[1 - swap], units))
                charts.append(tuple(chart[vertex] for vertex in range(8)))
    answer = tuple(sorted(set(charts)))
    require(len(answer) == 72, ("chart count", len(answer)))
    return answer


CHARTS = role_charts()


def tree_factor(
    potentials: dict[str, int],
    chart: tuple[str, ...],
    vertices: tuple[int, int, int, int],
    prime: int,
) -> int:
    positions = {vertex: index for index, vertex in enumerate(vertices)}
    laplacian = [[0] * 4 for _ in range(4)]
    for left, right in EDGES:
        if left not in positions or right not in positions:
            continue
        weight = (potentials[chart[left]] - potentials[chart[right]]) % prime
        i = positions[left]
        j = positions[right]
        laplacian[i][i] = (laplacian[i][i] + weight) % prime
        laplacian[j][j] = (laplacian[j][j] + weight) % prime
        laplacian[i][j] = (laplacian[i][j] - weight) % prime
        laplacian[j][i] = (laplacian[j][i] - weight) % prime
    minor = tuple(tuple(row[:3]) for row in laplacian[:3])
    return determinant_mod(minor, prime)


def chart_factors(q_values: dict[tuple[int, int, int], int], prime: int):
    potentials = {label: q_values[q] for label, q in ROLE_CLASSES.items()}
    rows = []
    for chart in CHARTS:
        bridge = (potentials["H"] - potentials["q5"]) % prime
        left = tree_factor(potentials, chart, WINGS[0] + (HUB,), prime)
        right = tree_factor(potentials, chart, WINGS[1] + (HUB,), prime)
        rows.append((bridge, left, right, bridge * left % prime * right % prime))
    return tuple(rows)


def chart_bank(
    entries: tuple[tuple[object, dict[tuple[int, int, int], int]], ...],
    prime: int,
) -> tuple[object, ...]:
    factor_zeros = [0, 0, 0, 0]
    product_zeros_by_chart = [0] * len(CHARTS)
    profiles: Counter[tuple[int, int, int, int]] = Counter()
    payload = []
    all_nonzero_entries = 0
    for key, values in entries:
        factors = chart_factors(values, prime)
        profile = tuple(
            sum(row[factor] == 0 for row in factors)
            for factor in range(4)
        )
        for factor, count in enumerate(profile):
            factor_zeros[factor] += count
        for chart, row in enumerate(factors):
            product_zeros_by_chart[chart] += row[3] == 0
        profiles[profile] += 1
        all_nonzero_entries += profile == (0, 0, 0, 0)
        payload.append((key, factors))
    return (
        len(entries),
        tuple(factor_zeros),
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
    bad = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Assert):
            bad.append("Assert")
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id in {"eval", "exec", "compile", "__import__"}
        ):
            bad.append(node.func.id)
    require(not bad, ("security", bad))
    return len(tuple(ast.walk(tree))), tuple(bad)


def main() -> None:
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(worker, range(P)))
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)), "worker order")
    tables = tuple(row for chunk in chunks for row in chunk[1])
    require(len(tables) == P**2, ("table count", len(tables)))

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
    ) = A.context()

    # Reconstruct the frozen bank independently of the all-role candidate.
    gamma = []
    table_index = 0
    for alpha in range(P):
        for beta in range(P):
            stored_beta, ax_values, by_values = tables[table_index]
            table_index += 1
            require(stored_beta == beta, ("table order", alpha, beta, stored_beta))
            phase = pow(zeta, beta, prime)
            for tau in range(P):
                ax = sum(
                    value * A.safe(chamber, sheet + tau)
                    for value, (sheet, chamber) in zip(ax_values, A.ATOMS)
                ) % prime
                by = sum(
                    value * A.safe(chamber, sheet + tau)
                    for value, (sheet, chamber) in zip(by_values, A.ATOMS)
                ) % prime
                gamma.append(phase * ax % prime * by % prime)
    gamma_tuple = tuple(gamma)
    gamma_hash = digest_integers(gamma_tuple)
    require(gamma_hash == A.EXPECTED_GAMMA_SHA256, ("gamma hash", gamma_hash))

    direct_values = {
        q: A.inverse_value(gamma_tuple, q, prime, zeta)
        for q in Q_CLASSES
    }
    require(direct_values == Q_VALUES, ("role values", direct_values))

    kernels = {}
    for tau_frequency in (0, 1):
        matrix = []
        for left_sheet, left_chamber in A.ATOMS:
            row = []
            for right_sheet, right_chamber in A.ATOMS:
                value = sum(
                    A.safe(left_chamber, left_sheet + tau)
                    * A.safe(right_chamber, right_sheet + tau)
                    * pow(zeta, -tau_frequency * tau % P, prime)
                    for tau in range(P)
                ) % prime
                row.append(value)
            matrix.append(tuple(row))
        kernels[tau_frequency] = tuple(matrix)
        require(
            all(value != 0 for row in kernels[tau_frequency] for value in row),
            ("guard kernel zero", tau_frequency),
        )

    accumulators = {q: [0] * len(A.BUCKETS) for q in Q_CLASSES}
    table_index = 0
    for alpha in range(P):
        for beta in range(P):
            _stored_beta, ax_values, by_values = tables[table_index]
            table_index += 1
            characters = {
                q: pow(zeta, beta - alpha * q[0] - beta * q[1], prime)
                for q in Q_CLASSES
            }
            for left_index, (left_sheet, left_chamber) in enumerate(A.ATOMS):
                left_value = ax_values[left_index]
                if left_value == 0:
                    continue
                for right_index, (right_sheet, right_chamber) in enumerate(A.ATOMS):
                    right_value = by_values[right_index]
                    if right_value == 0:
                        continue
                    drift = (right_sheet - left_sheet) % P
                    bucket = A.BUCKET_INDEX[(left_chamber, right_chamber, drift)]
                    product = left_value * right_value % prime
                    for q in Q_CLASSES:
                        contribution = (
                            characters[q]
                            * product
                            % prime
                            * kernels[q[2]][left_index][right_index]
                        ) % prime
                        accumulators[q][bucket] = (
                            accumulators[q][bucket] + contribution
                        ) % prime
    normalizer = pow(P**3, -1, prime)
    q_buckets = {
        q: tuple(value * normalizer % prime for value in accumulators[q])
        for q in Q_CLASSES
    }
    for q in Q_CLASSES:
        require(sum(q_buckets[q]) % prime == direct_values[q], ("bucket sum", q))

    active = tuple(
        index for index, bucket in enumerate(A.BUCKETS)
        if "middle" not in bucket[:2]
    )
    inactive = tuple(
        index for index, bucket in enumerate(A.BUCKETS)
        if "middle" in bucket[:2]
    )
    require((len(active), len(inactive)) == (52, 65), (len(active), len(inactive)))
    require(
        all(q_buckets[q][index] == 0 for q in Q_CLASSES for index in inactive),
        "middle support",
    )
    support_counts = tuple(
        (q, sum(q_buckets[q][index] != 0 for index in active))
        for q in Q_CLASSES
    )
    require(all(count == 52 for _q, count in support_counts), support_counts)
    tensor_rows = tuple(
        tuple(q_buckets[q][index] for index in active) for q in Q_CLASSES
    )
    tensor_rank = A.rank_mod(tensor_rows, prime)
    require(tensor_rank == 5, ("five-class rank", tensor_rank))

    point_entries = tuple(
        (A.BUCKETS[index], {q: q_buckets[q][index] for q in Q_CLASSES})
        for index in active
    )
    point_bank = chart_bank(point_entries, prime)

    drift_values = {
        q: tuple(
            sum(
                q_buckets[q][A.BUCKET_INDEX[(left, right, drift)]]
                for left, right in CORNERS
            ) % prime
            for drift in range(P)
        )
        for q in Q_CLASSES
    }
    drift_entries = tuple(
        (drift, {q: drift_values[q][drift] for q in Q_CLASSES})
        for drift in range(P)
    )
    drift_bank = chart_bank(drift_entries, prime)

    spectral_values = {
        q: tuple(
            dft(drift_values[q], frequency, zeta, prime)
            for frequency in range(P)
        )
        for q in Q_CLASSES
    }
    require(
        all(spectral_values[q][0] == direct_values[q] for q in Q_CLASSES),
        "zero Fourier mode",
    )
    spectral_entries = tuple(
        (frequency, {q: spectral_values[q][frequency] for q in Q_CLASSES})
        for frequency in range(P)
    )
    spectral_bank = chart_bank(spectral_entries, prime)

    drift_slice_ranks = []
    for drift in range(P):
        matrix = tuple(
            tuple(
                q_buckets[q][A.BUCKET_INDEX[(left, right, drift)]]
                for left, right in CORNERS
            )
            for q in Q_CLASSES
        )
        drift_slice_ranks.append(A.rank_mod(matrix, prime))

    chamber_spectral_entries = []
    spectral_slice_ranks = []
    for frequency in range(P):
        matrix = []
        for q in Q_CLASSES:
            matrix.append(
                tuple(
                    dft(
                        tuple(
                            q_buckets[q][A.BUCKET_INDEX[(left, right, drift)]]
                            for drift in range(P)
                        ),
                        frequency,
                        zeta,
                        prime,
                    )
                    for left, right in CORNERS
                )
            )
        matrix_tuple = tuple(matrix)
        spectral_slice_ranks.append(A.rank_mod(matrix_tuple, prime))
        for corner, (left, right) in enumerate(CORNERS):
            chamber_spectral_entries.append(
                (
                    (left, right, frequency),
                    {q: matrix_tuple[index][corner] for index, q in enumerate(Q_CLASSES)},
                )
            )
    chamber_spectral_bank = chart_bank(tuple(chamber_spectral_entries), prime)

    global_bank = chart_bank((("global", direct_values),), prime)
    banks = {
        "point": point_bank,
        "drift": drift_bank,
        "fourier": spectral_bank,
        "chamber_fourier": chamber_spectral_bank,
    }
    for name, bank in banks.items():
        require(bank[1] == (0, 0, 0, 0), ("factor zero", name, bank[1]))
        require(bank[2] == bank[0], ("not all entries close", name, bank[:3]))
        require(bank[5] == EXPECTED_BANK_DIGESTS[name], ("bank digest", name, bank[5]))
    require(tuple(drift_slice_ranks) == (4,) * P, drift_slice_ranks)
    require(tuple(spectral_slice_ranks) == (4,) * P, spectral_slice_ranks)
    require(global_bank[1] == (0, 0, 0, 0), ("global factors", global_bank))

    flat = {q: 1 for q in Q_CLASSES}
    flat_factors = chart_factors(flat, prime)
    require(all(row == (0, 0, 0, 0) for row in flat_factors), "flat hostile")
    killed = dict(direct_values)
    killed[ROLE_CLASSES["q5"]] = direct_values[ROLE_CLASSES["H"]]
    killed_factors = chart_factors(killed, prime)
    require(all(row[0] == 0 and row[3] == 0 for row in killed_factors), "bridge hostile")

    security = security_certificate(ROOT / SCRIPT)
    semantic = (
        ATOM_AUDIT_SHA256,
        CANDIDATE_SEMANTIC_SHA256,
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
        "all endpoint edge banks are graph coboundaries B^1",
        "no ancestry, current, absolute H1, or LRC14",
    )
    semantic_hash = hashlib.sha256(repr(semantic).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, ("semantic hash", semantic_hash))

    print("LRC U_FULL ALL-ROLE WEIGHTED-TREE SPECTRAL CLOSURE INDEPENDENT AUDIT")
    print("status=FINITE-EXACT INDEPENDENT AUDIT; endpoint B^1 only; LRC(14) OPEN")
    print(f"atom_dependency={ATOM_AUDIT_PATH.name}:{ATOM_AUDIT_SHA256}; imported_all_role_candidate=False")
    print(f"field=(prime={prime},zeta13={zeta}); gamma_sha256={gamma_hash}")
    print(f"role_classes={tuple(ROLE_CLASSES.items())}")
    print(f"distinct_q_values={tuple((q,direct_values[q]) for q in Q_CLASSES)}")
    print(f"five_class_tensor=(shape=(5,4,13),support_counts={support_counts},rank={tensor_rank})")
    print(f"pointwise_52x72=(factor_zero_totals={point_bank[1]},all_nonzero_entries={point_bank[2]}/{point_bank[0]},digest={point_bank[5]})")
    print(f"drift_marginal_13x72=(factor_zero_totals={drift_bank[1]},all_nonzero_entries={drift_bank[2]}/{drift_bank[0]},digest={drift_bank[5]})")
    print(f"drift_fourier_13x72=(factor_zero_totals={spectral_bank[1]},all_nonzero_entries={spectral_bank[2]}/{spectral_bank[0]},digest={spectral_bank[5]})")
    print(f"chamber_fourier_52x72=(factor_zero_totals={chamber_spectral_bank[1]},all_nonzero_entries={chamber_spectral_bank[2]}/{chamber_spectral_bank[0]},digest={chamber_spectral_bank[5]})")
    print(f"slice_ranks=(drift={tuple(drift_slice_ranks)},fourier={tuple(spectral_slice_ranks)})")
    print(f"global_chart_gate=(factor_zero_totals={global_bank[1]},digest={global_bank[5]})")
    print("hostiles=(flat potentials kill every factor; forcing q5=H kills every bridge and product)")
    print("typing=all thirteen edge responses are differences of eight endpoint role potentials, hence lie in graph B^1; no nonzero absolute H1 class")
    print("nonconsequence=no common ancestry, physical current, grouped coefficient, all-unit projector, bispectrum theorem, scalar-row exclusion, or LRC(14)")
    print(f"security_ast={security}")
    print(f"semantic_sha256={semantic_hash}")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
