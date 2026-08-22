#!/usr/bin/env python3
"""Independent audit of the U_full all-role modewise endpoint spectrum.

This audit imports the already promoted THM-3514 atom-table implementation,
not the candidate all-role script or either of the candidate's two helper
modules.  Starting from the materialized unguarded endpoint atom tables, it
reconstructs the five refined role-class bucket tensors in two ways:

1. a direct three-character route that first forms every tau/chamber/drift
   slice and then applies the alpha, beta, and tau inverse characters; and
2. a canonical guard-kernel route that derives and checks the left-sheet
   translation phase before contraction.

The graph route is also independent: all K4 factors are evaluated both by
explicit enumeration of the sixteen spanning trees and by reduced Laplacian
determinants.  The full eight-vertex determinant is checked against the two
K4 factors and forced bridge.  Six fundamental cycle pairings certify that
all resulting edge responses are coboundaries.  Thus this is an endpoint
B^1 audit only, not ancestry, physical current, absolute H^1 flux, a
bispectrum theorem, a row exclusion, or LRC(14).
"""

from __future__ import annotations

import ast
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations, permutations
import hashlib
import importlib.util
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "04-computation/lrc_ufull_all_role_modewise_spectral_independent_audit_20260816.py"
OUTPUT = "05-knowledge/results/lrc_ufull_all_role_modewise_spectral_independent_audit_20260816.out"
ATOM_AUDIT_PATH = ROOT / (
    "04-computation/"
    "lrc_ufull_owner_boundary_k4xf13_endpoint_factorization_"
    "independent_audit_20260816.py"
)
ATOM_AUDIT_SHA256 = "f89be10c65bb77270199f9399b155d5a2c82c0da121b3e8589fe3c1f7e9824fc"
CANDIDATE_PATH = ROOT / "04-computation/lrc_ufull_guard_bucket_all_role_spectral_probe_20260816.py"
CANDIDATE_SHA256 = "98a4cf5c82ca10027302baf2c7fb59acb0f305143e22453d0a8660fef8d90cf0"
CANDIDATE_COMMIT = "b5e18d2b1c84"
CANDIDATE_BUCKET_DEP_SHA256 = "a1d4b667812949001fc863ba881ff7409bbae3c568a6bf7bc24c9dc88b2766b1"
CANDIDATE_ROLE_DEP_SHA256 = "ee2105742abee578a9c41ff7ec954a07ada324fccc2c643429e7ac6e6e6f8fc2"
CANDIDATE_SECURITY = (3458, ())
CANDIDATE_SEMANTIC_SHA256 = "2c9495fb8bcb731361ba331d9ca4b84a60f21551dc49b16e6519c1fc4f2e9f97"
EXPECTED_GAMMA_SHA256 = "1fabc5cfdbaa1455e10cd6bf9264488133616a7b0ff381623d729b4b4bfa9682"
EXPECTED_VALUES = {
    (0, 0, 0): 405336876493642499425,
    (0, 1, 0): 518539850465495448196,
    (1, 0, 0): 503604956476841920373,
    (1, 0, 1): 320618948602619577408,
    (1, 12, 0): 15703541686881447885,
}
EXPECTED_BANK_DIGESTS = (
    "a1d5183f37cc39deed976876ed91132f662f91137f41a10a79fd3e974ced2dfc",
    "6bcdd9fa616ba2dc4883d8334dca7358d06375d7342f43c4aac8dea45deb5027",
    "e96a3859e6838fa91d04a55ca55fa4e44cafa45bd86b93c7f6b612b1e65f4dc1",
    "c9ca8c9b7ed7000b00e39b496523ad75a44379d84ea35e49163d6b64e4196e73",
)
EXPECTED_SEMANTIC_SHA256 = "b745ab5d95936b9a134cfe8b6ffde032e72020d950670ab8b1a4643dd58c8be6"

P = 13
Q_CLASSES = tuple(sorted(EXPECTED_VALUES))
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
CHAMBER_NAMES = ("left", "middle", "right")
CORNER_PAIRS = (
    ("left", "left"),
    ("left", "right"),
    ("right", "left"),
    ("right", "right"),
)
ATOMS = tuple((sheet, chamber) for sheet in range(P) for chamber in CHAMBER_NAMES)
BUCKETS = tuple(
    (left, right, drift)
    for left in CHAMBER_NAMES
    for right in CHAMBER_NAMES
    for drift in range(P)
)
BUCKET_INDEX = {bucket: index for index, bucket in enumerate(BUCKETS)}

# Private-support graph: two K4 blocks sharing the hub, plus one leaf bridge.
EDGES = (
    (0, 3), (0, 4), (0, 5),
    (1, 2), (1, 4), (1, 7),
    (2, 4), (2, 7),
    (3, 4), (3, 5),
    (4, 5), (4, 6), (4, 7),
)
EDGE_INDEX = {edge: index for index, edge in enumerate(EDGES)}
HUB, LEAF = 4, 6
WINGS = ((0, 3, 5), (1, 2, 7))
BLOCKERS = ("c1", "c2", "c3")
UNITS = ("q2", "q3", "q4")


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


def load_atom_audit():
    require(lf_sha256(ATOM_AUDIT_PATH) == ATOM_AUDIT_SHA256, "THM-3514 audit source drift")
    name = "thm3514_atom_table_parent"
    spec = importlib.util.spec_from_file_location(name, ATOM_AUDIT_PATH)
    require(spec is not None and spec.loader is not None, "THM-3514 audit loader")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


A = load_atom_audit()


def atom_worker(alpha: int) -> tuple[object, ...]:
    return A.worker(alpha)


def zpower(zeta: int, exponent: int, prime: int) -> int:
    return pow(zeta, exponent % P, prime)


def inverse_value(
    gamma: tuple[int, ...], q: tuple[int, int, int], prime: int, zeta: int
) -> int:
    total = 0
    index = 0
    for alpha in range(P):
        for beta in range(P):
            for tau in range(P):
                exponent = -(alpha * q[0] + beta * q[1] + tau * q[2])
                total = (total + gamma[index] * zpower(zeta, exponent, prime)) % prime
                index += 1
    return total * pow(P**3, -1, prime) % prime


def det_mod(matrix: tuple[tuple[int, ...], ...], prime: int) -> int:
    rows = [list(value % prime for value in row) for row in matrix]
    require(rows and all(len(row) == len(rows) for row in rows), "square determinant")
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
        pivot_value = rows[column][column]
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


def rank_mod(matrix: tuple[tuple[int, ...], ...], prime: int) -> int:
    require(matrix and matrix[0], "nonempty rank matrix")
    rows = [list(value % prime for value in row) for row in matrix]
    rank = 0
    for column in range(len(rows[0])):
        pivot = next(
            (row for row in range(rank, len(rows)) if rows[row][column]),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, prime)
        rows[rank] = [value * inverse % prime for value in rows[rank]]
        for row in range(len(rows)):
            if row == rank:
                continue
            factor = rows[row][column]
            if factor:
                rows[row] = [
                    (value - factor * pivot_value) % prime
                    for value, pivot_value in zip(rows[row], rows[rank])
                ]
        rank += 1
        if rank == len(rows):
            break
    return rank


def first_full_row_minor(
    matrix: tuple[tuple[int, ...], ...], prime: int
) -> tuple[tuple[int, ...], int]:
    rows = len(matrix)
    columns = len(matrix[0])
    require(rows <= columns, (rows, columns))
    for choice in combinations(range(columns), rows):
        minor = tuple(tuple(row[column] for column in choice) for row in matrix)
        determinant = det_mod(minor, prime)
        if determinant:
            return choice, determinant
    raise RuntimeError("no full-row-rank minor")


def first_full_column_minor(
    matrix: tuple[tuple[int, ...], ...], prime: int
) -> tuple[tuple[int, ...], int]:
    rows = len(matrix)
    columns = len(matrix[0])
    require(columns <= rows, (rows, columns))
    for choice in combinations(range(rows), columns):
        minor = tuple(matrix[row] for row in choice)
        determinant = det_mod(minor, prime)
        if determinant:
            return choice, determinant
    raise RuntimeError("no full-column-rank minor")


def role_charts() -> tuple[tuple[str, ...], ...]:
    charts = []
    for swap in (0, 1):
        for blockers in permutations(BLOCKERS):
            for units in permutations(UNITS):
                chart = {HUB: "H", LEAF: "q5"}
                chart.update(zip(WINGS[swap], blockers))
                chart.update(zip(WINGS[1 - swap], units))
                charts.append(tuple(chart[index] for index in range(8)))
    answer = tuple(sorted(set(charts)))
    require(len(answer) == 72, ("chart count", len(answer)))
    return answer


CHARTS = role_charts()


def edge_weight(potentials: tuple[int, ...], edge: tuple[int, int], prime: int) -> int:
    left, right = edge
    require(left < right, ("orientation", edge))
    return (potentials[left] - potentials[right]) % prime


def k4_spanning_trees(vertices: tuple[int, ...]) -> tuple[tuple[tuple[int, int], ...], ...]:
    vertex_set = frozenset(vertices)
    edges = tuple(edge for edge in EDGES if set(edge) <= vertex_set)
    require(len(edges) == 6, (vertices, edges))
    trees = []
    for chosen in combinations(edges, 3):
        parent = {vertex: vertex for vertex in vertices}

        def find(vertex: int) -> int:
            while parent[vertex] != vertex:
                parent[vertex] = parent[parent[vertex]]
                vertex = parent[vertex]
            return vertex

        acyclic = True
        for left, right in chosen:
            root_left = find(left)
            root_right = find(right)
            if root_left == root_right:
                acyclic = False
                break
            parent[root_left] = root_right
        if acyclic and len({find(vertex) for vertex in vertices}) == 1:
            trees.append(chosen)
    answer = tuple(trees)
    require(len(answer) == 16, ("K4 tree count", vertices, len(answer)))
    return answer


LEFT_VERTICES = WINGS[0] + (HUB,)
RIGHT_VERTICES = WINGS[1] + (HUB,)
LEFT_TREES = k4_spanning_trees(LEFT_VERTICES)
RIGHT_TREES = k4_spanning_trees(RIGHT_VERTICES)


def k4_tree_sum(
    potentials: tuple[int, ...], trees: tuple[tuple[tuple[int, int], ...], ...], prime: int
) -> int:
    total = 0
    for tree in trees:
        product = 1
        for edge in tree:
            product = product * edge_weight(potentials, edge, prime) % prime
        total = (total + product) % prime
    return total


def reduced_laplacian_det(
    potentials: tuple[int, ...], vertices: tuple[int, ...], prime: int
) -> int:
    positions = {vertex: index for index, vertex in enumerate(vertices)}
    laplacian = [[0] * len(vertices) for _vertex in vertices]
    for edge in EDGES:
        left, right = edge
        if left not in positions or right not in positions:
            continue
        weight = edge_weight(potentials, edge, prime)
        i, j = positions[left], positions[right]
        laplacian[i][i] = (laplacian[i][i] + weight) % prime
        laplacian[j][j] = (laplacian[j][j] + weight) % prime
        laplacian[i][j] = (laplacian[i][j] - weight) % prime
        laplacian[j][i] = (laplacian[j][i] - weight) % prime
    reduced = tuple(tuple(row[:-1]) for row in laplacian[:-1])
    return det_mod(reduced, prime)


def full_graph_det(potentials: tuple[int, ...], prime: int) -> int:
    return reduced_laplacian_det(potentials, tuple(range(8)), prime)


def factor_row(
    q_values: dict[tuple[int, int, int], int], chart: tuple[str, ...], prime: int
) -> tuple[int, int, int, int]:
    values = {label: q_values[q] for label, q in ROLE_CLASSES.items()}
    potentials = tuple(values[label] for label in chart)
    bridge = edge_weight(potentials, (HUB, LEAF), prime)
    left_tree = k4_tree_sum(potentials, LEFT_TREES, prime)
    right_tree = k4_tree_sum(potentials, RIGHT_TREES, prime)
    left_det = reduced_laplacian_det(potentials, LEFT_VERTICES, prime)
    right_det = reduced_laplacian_det(potentials, RIGHT_VERTICES, prime)
    require((left_tree, right_tree) == (left_det, right_det), "tree/Laplacian mismatch")
    product = bridge * left_tree % prime * right_tree % prime
    require(full_graph_det(potentials, prime) == product, "cut-vertex determinant mismatch")
    return bridge, left_tree, right_tree, product


def cycle_basis() -> tuple[tuple[int, ...], ...]:
    rows = []
    for wing in WINGS:
        for left, right in combinations(wing, 2):
            path = (HUB, left, right, HUB)
            row = [0] * len(EDGES)
            for start, stop in zip(path, path[1:]):
                edge = (min(start, stop), max(start, stop))
                row[EDGE_INDEX[edge]] += 1 if start < stop else -1
            rows.append(tuple(row))
    answer = tuple(rows)
    require(len(answer) == 6, len(answer))
    return answer


CYCLES = cycle_basis()


def graph_typing_certificate(prime: int) -> tuple[object, ...]:
    incidence = [[0] * len(EDGES) for _vertex in range(8)]
    for column, (left, right) in enumerate(EDGES):
        incidence[left][column] = 1
        incidence[right][column] = -1
    incidence_t = tuple(tuple(row) for row in incidence)
    incidence_rank = rank_mod(incidence_t, prime)
    cycle_rank = rank_mod(CYCLES, prime)
    for cycle in CYCLES:
        require(
            all(sum(incidence[row][edge] * cycle[edge] for edge in range(len(EDGES))) == 0
                for row in range(8)),
            ("cycle boundary", cycle),
        )
    bridge_index = EDGE_INDEX[(HUB, LEAF)]
    require(all(cycle[bridge_index] == 0 for cycle in CYCLES), "bridge entered cycle")
    require((incidence_rank, cycle_rank) == (7, 6), (incidence_rank, cycle_rank))
    return (8, len(EDGES), incidence_rank, cycle_rank, bridge_index)


def check_cycle_pairings(
    entries: tuple[tuple[object, dict[tuple[int, int, int], int]], ...], prime: int
) -> int:
    checked = 0
    for _key, q_values in entries:
        values = {label: q_values[q] for label, q in ROLE_CLASSES.items()}
        for chart in CHARTS:
            potentials = tuple(values[label] for label in chart)
            weights = tuple(edge_weight(potentials, edge, prime) for edge in EDGES)
            for cycle in CYCLES:
                pairing = sum(value * coefficient for value, coefficient in zip(weights, cycle)) % prime
                require(pairing == 0, ("non-coboundary pairing", chart, pairing))
                checked += 1
    return checked


def chart_bank(
    entries: tuple[tuple[object, dict[tuple[int, int, int], int]], ...], prime: int
) -> tuple[object, ...]:
    factor_zero_totals = [0, 0, 0, 0]
    product_zeros_by_chart = [0] * len(CHARTS)
    profiles: Counter[tuple[int, int, int, int]] = Counter()
    payload = []
    all_nonzero_entries = 0
    for key, q_values in entries:
        factors = tuple(factor_row(q_values, chart, prime) for chart in CHARTS)
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


def canonical_guard_kernel(
    left: str, right: str, drift: int, frequency: int, zeta: int, prime: int
) -> int:
    return sum(
        A.safe(left, sheet)
        * A.safe(right, sheet + drift)
        * zpower(zeta, -frequency * sheet, prime)
        for sheet in range(P)
    ) % prime


def guard_phase_certificate(zeta: int, prime: int) -> tuple[object, ...]:
    canonical = {
        (left, right, drift, frequency): canonical_guard_kernel(
            left, right, drift, frequency, zeta, prime
        )
        for left in CHAMBER_NAMES
        for right in CHAMBER_NAMES
        for drift in range(P)
        for frequency in (0, 1)
    }
    checks = 0
    for left_sheet, left in ATOMS:
        for right_sheet, right in ATOMS:
            drift = (right_sheet - left_sheet) % P
            for frequency in (0, 1):
                direct = sum(
                    A.safe(left, left_sheet + tau)
                    * A.safe(right, right_sheet + tau)
                    * zpower(zeta, -frequency * tau, prime)
                    for tau in range(P)
                ) % prime
                translated = (
                    zpower(zeta, frequency * left_sheet, prime)
                    * canonical[(left, right, drift, frequency)]
                ) % prime
                require(direct == translated, ("guard translation phase", left_sheet, right_sheet))
                checks += 1
    return canonical, (checks, digest_json(tuple(sorted(canonical.items()))))


def dft(values: tuple[int, ...], frequency: int, zeta: int, prime: int) -> int:
    return sum(
        value * zpower(zeta, -frequency * drift, prime)
        for drift, value in enumerate(values)
    ) % prime


def inverse_dft(spectrum: tuple[int, ...], drift: int, zeta: int, prime: int) -> int:
    return (
        sum(
            value * zpower(zeta, frequency * drift, prime)
            for frequency, value in enumerate(spectrum)
        )
        * pow(P, -1, prime)
    ) % prime


def security_certificate(path: Path) -> tuple[int, tuple[str, ...]]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    forbidden = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Assert):
            forbidden.append("Assert")
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id in {"eval", "exec", "compile", "__import__"}
        ):
            forbidden.append(node.func.id)
    require(not forbidden, ("security", forbidden))
    return len(tuple(ast.walk(tree))), tuple(forbidden)


def main() -> None:
    require(lf_sha256(CANDIDATE_PATH) == CANDIDATE_SHA256, "candidate source drift")
    require(tuple(A.ATOMS) == ATOMS, "atom order mismatch")
    require(tuple(A.BUCKETS) == BUCKETS, "bucket order mismatch")
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(atom_worker, range(P)))
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)), "worker order")

    tables = tuple(row for chunk in chunks for row in chunk[1])
    require(len(tables) == P**2, ("atom table count", len(tables)))
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

    canonical_kernels, phase_certificate = guard_phase_certificate(zeta, prime)
    direct_raw = {q: [0] * len(BUCKETS) for q in Q_CLASSES}
    wrong_target_raw = {q: [0] * len(BUCKETS) for q in Q_CLASSES}
    wrong_tau_raw = {q: [0] * len(BUCKETS) for q in Q_CLASSES}
    gamma_rows = []

    # Route 1: materialize each tau slice first.  This does not use pair
    # kernels and therefore tests all three inverse-character signs directly.
    table_index = 0
    for alpha in range(P):
        for beta in range(P):
            stored_beta, ax_values, by_values = tables[table_index]
            table_index += 1
            require(stored_beta == beta, ("beta order", alpha, beta, stored_beta))
            active_left = tuple(index for index, value in enumerate(ax_values) if value)
            active_right = tuple(index for index, value in enumerate(by_values) if value)
            for tau in range(P):
                slice_values = [0] * len(BUCKETS)
                left_guarded = tuple(
                    (index, ax_values[index])
                    for index in active_left
                    if A.safe(ATOMS[index][1], ATOMS[index][0] + tau)
                )
                right_guarded = tuple(
                    (index, by_values[index])
                    for index in active_right
                    if A.safe(ATOMS[index][1], ATOMS[index][0] + tau)
                )
                ax_guarded = sum(value for _index, value in left_guarded) % prime
                by_guarded = sum(value for _index, value in right_guarded) % prime
                for left_index, left_value in left_guarded:
                    left_sheet, left_chamber = ATOMS[left_index]
                    for right_index, right_value in right_guarded:
                        right_sheet, right_chamber = ATOMS[right_index]
                        drift = (right_sheet - left_sheet) % P
                        bucket = BUCKET_INDEX[(left_chamber, right_chamber, drift)]
                        slice_values[bucket] = (
                            slice_values[bucket] + left_value * right_value
                        ) % prime
                slice_total = sum(slice_values) % prime
                require(slice_total == ax_guarded * by_guarded % prime, "tau slice product")
                gamma_rows.append(zpower(zeta, beta, prime) * slice_total % prime)
                for q in Q_CLASSES:
                    qa, qb, qt = q
                    correct = beta - alpha * qa - beta * qb - tau * qt
                    wrong_target = beta + alpha * qa + beta * qb - tau * qt
                    wrong_tau = beta - alpha * qa - beta * qb + tau * qt
                    correct_weight = zpower(zeta, correct, prime)
                    wrong_target_weight = zpower(zeta, wrong_target, prime)
                    wrong_tau_weight = zpower(zeta, wrong_tau, prime)
                    for bucket, value in enumerate(slice_values):
                        if value:
                            direct_raw[q][bucket] = (
                                direct_raw[q][bucket] + correct_weight * value
                            ) % prime
                            wrong_target_raw[q][bucket] = (
                                wrong_target_raw[q][bucket] + wrong_target_weight * value
                            ) % prime
                            wrong_tau_raw[q][bucket] = (
                                wrong_tau_raw[q][bucket] + wrong_tau_weight * value
                            ) % prime

    gamma = tuple(gamma_rows)
    require(len(gamma) == P**3, ("gamma size", len(gamma)))
    gamma_hash = digest_integers(gamma)
    require(gamma_hash == EXPECTED_GAMMA_SHA256, ("gamma digest", gamma_hash))
    direct_values = {q: inverse_value(gamma, q, prime, zeta) for q in Q_CLASSES}
    require(direct_values == EXPECTED_VALUES, ("five role values", direct_values))
    normalizer = pow(P**3, -1, prime)
    direct_buckets = {
        q: tuple(value * normalizer % prime for value in direct_raw[q])
        for q in Q_CLASSES
    }

    # Route 2: contract the canonical guard kernels.  The translation identity
    # supplies zeta^(q_t*a_L), whose sign was checked above on all 3042
    # atom-pair/frequency instances.
    kernel_raw = {q: [0] * len(BUCKETS) for q in Q_CLASSES}
    table_index = 0
    for alpha in range(P):
        for beta in range(P):
            stored_beta, ax_values, by_values = tables[table_index]
            table_index += 1
            require(stored_beta == beta, "kernel beta order")
            for left_index, left_value in enumerate(ax_values):
                if left_value == 0:
                    continue
                left_sheet, left_chamber = ATOMS[left_index]
                for right_index, right_value in enumerate(by_values):
                    if right_value == 0:
                        continue
                    right_sheet, right_chamber = ATOMS[right_index]
                    drift = (right_sheet - left_sheet) % P
                    bucket = BUCKET_INDEX[(left_chamber, right_chamber, drift)]
                    product = left_value * right_value % prime
                    for q in Q_CLASSES:
                        qa, qb, qt = q
                        target_weight = zpower(
                            zeta, beta - alpha * qa - beta * qb, prime
                        )
                        translated_kernel = (
                            zpower(zeta, qt * left_sheet, prime)
                            * canonical_kernels[(left_chamber, right_chamber, drift, qt)]
                        ) % prime
                        kernel_raw[q][bucket] = (
                            kernel_raw[q][bucket]
                            + target_weight * product % prime * translated_kernel
                        ) % prime
    kernel_buckets = {
        q: tuple(value * normalizer % prime for value in kernel_raw[q])
        for q in Q_CLASSES
    }
    require(kernel_buckets == direct_buckets, "direct/kernel bucket disagreement")
    q_buckets = direct_buckets
    for q in Q_CLASSES:
        require(sum(q_buckets[q]) % prime == direct_values[q], ("bucket sum", q))

    # Hostiles demonstrate that the signs and P^-3 normalizer are not vacuous.
    wrong_target_buckets = {
        q: tuple(value * normalizer % prime for value in wrong_target_raw[q])
        for q in Q_CLASSES
    }
    wrong_tau_buckets = {
        q: tuple(value * normalizer % prime for value in wrong_tau_raw[q])
        for q in Q_CLASSES
    }
    target_sign_witness = next(
        (q, BUCKETS[index], q_buckets[q][index], wrong_target_buckets[q][index])
        for q in Q_CLASSES
        for index in range(len(BUCKETS))
        if q_buckets[q][index] != wrong_target_buckets[q][index]
    )
    tau_sign_witness = next(
        (q, BUCKETS[index], q_buckets[q][index], wrong_tau_buckets[q][index])
        for q in Q_CLASSES
        for index in range(len(BUCKETS))
        if q_buckets[q][index] != wrong_tau_buckets[q][index]
    )
    first_q = Q_CLASSES[0]
    first_nonzero = next(index for index, value in enumerate(direct_raw[first_q]) if value)
    wrong_normalizer_value = direct_raw[first_q][first_nonzero] * pow(P**2, -1, prime) % prime
    normalizer_witness = (
        first_q,
        BUCKETS[first_nonzero],
        q_buckets[first_q][first_nonzero],
        wrong_normalizer_value,
    )
    require(normalizer_witness[2] != normalizer_witness[3], "normalizer hostile collapsed")

    active_indices = tuple(
        index for index, bucket in enumerate(BUCKETS) if "middle" not in bucket[:2]
    )
    inactive_indices = tuple(
        index for index, bucket in enumerate(BUCKETS) if "middle" in bucket[:2]
    )
    require((len(active_indices), len(inactive_indices)) == (52, 65), "bucket census")
    require(
        all(q_buckets[q][index] == 0 for q in Q_CLASSES for index in inactive_indices),
        "middle support",
    )
    support_counts = tuple(
        (q, sum(q_buckets[q][index] != 0 for index in active_indices))
        for q in Q_CLASSES
    )
    require(all(count == 52 for _q, count in support_counts), support_counts)
    tensor = tuple(
        tuple(q_buckets[q][index] for index in active_indices) for q in Q_CLASSES
    )
    tensor_rank = rank_mod(tensor, prime)
    require(tensor_rank == 5, ("five-by-52 rank", tensor_rank))
    tensor_minor = first_full_row_minor(tensor, prime)

    point_entries = tuple(
        (BUCKETS[index], {q: q_buckets[q][index] for q in Q_CLASSES})
        for index in active_indices
    )
    point_bank = chart_bank(point_entries, prime)

    drift_values = {
        q: tuple(
            sum(
                q_buckets[q][BUCKET_INDEX[(left, right, drift)]]
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
        "zero Fourier mode",
    )
    for q in Q_CLASSES:
        reconstructed = tuple(
            inverse_dft(spectral_values[q], drift, zeta, prime) for drift in range(P)
        )
        require(reconstructed == drift_values[q], ("DFT inverse sign", q))
    spectral_entries = tuple(
        (frequency, {q: spectral_values[q][frequency] for q in Q_CLASSES})
        for frequency in range(P)
    )
    spectral_bank = chart_bank(spectral_entries, prime)

    plus_sign_witness = next(
        (
            q,
            frequency,
            spectral_values[q][frequency],
            sum(
                value * zpower(zeta, frequency * drift, prime)
                for drift, value in enumerate(drift_values[q])
            ) % prime,
        )
        for q in Q_CLASSES
        for frequency in range(1, P)
        if spectral_values[q][frequency]
        != sum(
            value * zpower(zeta, frequency * drift, prime)
            for drift, value in enumerate(drift_values[q])
        ) % prime
    )

    chamber_spectral_entries = []
    drift_slice_ranks = []
    spectral_slice_ranks = []
    drift_minor_certificates = []
    spectral_minor_certificates = []
    for drift in range(P):
        matrix = tuple(
            tuple(
                q_buckets[q][BUCKET_INDEX[(left, right, drift)]]
                for left, right in CORNER_PAIRS
            )
            for q in Q_CLASSES
        )
        rank = rank_mod(matrix, prime)
        drift_slice_ranks.append(rank)
        require(rank == 4, ("drift slice rank", drift, rank))
        drift_minor_certificates.append(first_full_column_minor(matrix, prime))
    for frequency in range(P):
        matrix_rows = []
        for q in Q_CLASSES:
            row = []
            for left, right in CORNER_PAIRS:
                values = tuple(
                    q_buckets[q][BUCKET_INDEX[(left, right, drift)]]
                    for drift in range(P)
                )
                spectrum = tuple(dft(values, mode, zeta, prime) for mode in range(P))
                reconstructed = tuple(
                    inverse_dft(spectrum, drift, zeta, prime) for drift in range(P)
                )
                require(reconstructed == values, ("chamber DFT inverse", q, left, right))
                row.append(spectrum[frequency])
            matrix_rows.append(tuple(row))
        matrix = tuple(matrix_rows)
        rank = rank_mod(matrix, prime)
        spectral_slice_ranks.append(rank)
        require(rank == 4, ("spectral slice rank", frequency, rank))
        spectral_minor_certificates.append(first_full_column_minor(matrix, prime))
        for corner_index, (left, right) in enumerate(CORNER_PAIRS):
            chamber_spectral_entries.append(
                (
                    (left, right, frequency),
                    {q: matrix[q_index][corner_index] for q_index, q in enumerate(Q_CLASSES)},
                )
            )
    chamber_spectral_entries_t = tuple(chamber_spectral_entries)
    chamber_spectral_bank = chart_bank(chamber_spectral_entries_t, prime)

    global_entries = (("global", direct_values),)
    global_bank = chart_bank(global_entries, prime)
    require(global_bank[1] == (0, 0, 0, 0), ("global bank", global_bank))
    flat = {q: 1 for q in Q_CLASSES}
    require(
        all(factor_row(flat, chart, prime) == (0, 0, 0, 0) for chart in CHARTS),
        "flat hostile",
    )
    killed = dict(direct_values)
    killed[(1, 0, 0)] = killed[(1, 0, 1)]
    killed_rows = tuple(factor_row(killed, chart, prime) for chart in CHARTS)
    require(all(row[0] == 0 and row[3] == 0 for row in killed_rows), "bridge hostile")

    bank_digests = (
        point_bank[-1], drift_bank[-1], spectral_bank[-1], chamber_spectral_bank[-1]
    )
    require(bank_digests == EXPECTED_BANK_DIGESTS, ("candidate bank digests", bank_digests))
    require(
        point_bank[1] == drift_bank[1] == spectral_bank[1]
        == chamber_spectral_bank[1] == (0, 0, 0, 0),
        "chart factor zero",
    )

    graph_typing = graph_typing_certificate(prime)
    cycle_checks = (
        check_cycle_pairings(point_entries, prime),
        check_cycle_pairings(drift_entries, prime),
        check_cycle_pairings(spectral_entries, prime),
        check_cycle_pairings(chamber_spectral_entries_t, prime),
        check_cycle_pairings(global_entries, prime),
    )

    # Exact reconstruction of the candidate's semantic schema is a terminal
    # comparison only.  No candidate code, bucket helper, or role helper was
    # imported in deriving any coefficient above.
    candidate_semantic = (
        CANDIDATE_BUCKET_DEP_SHA256,
        CANDIDATE_ROLE_DEP_SHA256,
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
        CANDIDATE_SECURITY,
        "endpoint role potentials only; graph coboundaries; no ancestry/current/H1/LRC14",
    )
    reconstructed_candidate_semantic = hashlib.sha256(
        repr(candidate_semantic).encode("utf-8")
    ).hexdigest()
    require(
        reconstructed_candidate_semantic == CANDIDATE_SEMANTIC_SHA256,
        ("candidate semantic mismatch", reconstructed_candidate_semantic),
    )

    security = security_certificate(ROOT / SCRIPT)
    audit_semantic = (
        CANDIDATE_COMMIT,
        CANDIDATE_SHA256,
        ATOM_AUDIT_SHA256,
        (prime, zeta),
        gamma_hash,
        phase_certificate,
        tuple((q, q_buckets[q]) for q in Q_CLASSES),
        tensor_minor,
        point_bank,
        drift_bank,
        spectral_bank,
        chamber_spectral_bank,
        tuple(drift_minor_certificates),
        tuple(spectral_minor_certificates),
        target_sign_witness,
        tau_sign_witness,
        normalizer_witness,
        plus_sign_witness,
        graph_typing,
        cycle_checks,
        reconstructed_candidate_semantic,
        security,
        "endpoint B1 only; no ancestry/current/absolute H1/bispectrum/row exclusion/LRC14",
    )
    semantic_hash = hashlib.sha256(repr(audit_semantic).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, ("audit semantic", semantic_hash))

    print("LRC U_FULL ALL-ROLE MODEWISE SPECTRAL INDEPENDENT AUDIT")
    print("status=FINITE-EXACT INDEPENDENT AUDIT; endpoint B^1 only; LRC(14) OPEN")
    print(f"candidate=(commit={CANDIDATE_COMMIT},source_sha256={CANDIDATE_SHA256},imported=False)")
    print(f"trusted_parent=(THM-3514 atom audit,source_sha256={ATOM_AUDIT_SHA256})")
    print(f"field=(prime={prime},zeta13={zeta}); gamma_sha256={gamma_hash}")
    print(f"role_classes={tuple((label, ROLE_CLASSES[label]) for label in ROLE_CLASSES)}")
    print(f"five_q_values={tuple((q, direct_values[q]) for q in Q_CLASSES)}")
    print(f"bucket_routes=(direct_tau_slices=PASS,canonical_guard_kernels=PASS,equal=True,phase_checks={phase_certificate[0]})")
    print(f"phase_hostiles=(target_sign={target_sign_witness},tau_sign={tau_sign_witness},drift_dft_sign={plus_sign_witness})")
    print(f"normalizer=(correct=13^-3,witness_wrong_13^-2={normalizer_witness})")
    print(f"active_support={support_counts}; five_by_52_rank={tensor_rank}; rank5_minor={tensor_minor}")
    print(f"pointwise_52x72_bank={point_bank}")
    print(f"drift_13x72_bank={drift_bank}")
    print(f"fourier_13x72_bank={spectral_bank}")
    print(f"chamber_fourier_52x72_bank={chamber_spectral_bank}")
    print(f"slice_ranks=(drift={tuple(drift_slice_ranks)},fourier={tuple(spectral_slice_ranks)})")
    print(f"rank_minor_digests=(drift={digest_json(tuple(drift_minor_certificates))},fourier={digest_json(tuple(spectral_minor_certificates))})")
    print(f"graph_typing=(vertices,edges,incidence_rank,cycle_rank,bridge_edge_index)={graph_typing}")
    print(f"cycle_pairing_checks=(point,drift,fourier,chamber_fourier,global)={cycle_checks}; all_zero=True")
    print("tree_checks=each K4 uses 16 explicit spanning trees; both reduced Laplacians and full 7x7 determinant agree at every bank entry/chart")
    print(f"candidate_semantic_reconstructed={reconstructed_candidate_semantic}")
    print("hostiles=(flat potentials kill all factors; H=q5 kills bridge/product; wrong phase signs and 13^-2 normalization have explicit witnesses)")
    print("typing=edge responses are gradients of eight role potentials, hence lie in B^1; the forced bridge occurs in no cycle and all six cycle pairings vanish")
    print("nonconsequence=no common ancestry, physical or grouped current, absolute H1 flux, all-unit projector, bispectrum theorem, scalar-row exclusion, or LRC(14)")
    print(f"security_ast={security}")
    print(f"semantic_sha256={semantic_hash}")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
