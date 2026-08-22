#!/usr/bin/env python3
"""Independent audit of the folded-C7/K4 drift-transporter obstruction.

This verifier does not import the submitted transporter probe or its U_full
parent.  It takes the proved THM-2594 aggregate response table as source and
reconstructs the THM-3514 H-q5 target directly from the independently audited
unguarded endpoint atoms.  The 48/192 allocation families are obtained by
filtering all 7^4 labelled choices, rather than by the submitted permutation
loops.  Fixed channel mixing plus a common circulant is audited through an
augmented linear system with explicit frequency multipliers.

The result is representation-theoretic and endpoint-only.  It constructs no
common ancestry relation, physical current, absolute H1 flux, bispectrum,
row exclusion, or LRC(14) theorem.
"""

from __future__ import annotations

import ast
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations, product
import hashlib
import importlib.util
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
ATOM_PATH = ROOT / (
    "04-computation/"
    "lrc_ufull_owner_boundary_k4xf13_endpoint_factorization_"
    "independent_audit_20260816.py"
)
R5_OUTPUT = ROOT / "05-knowledge/results/lrc14_stage2_theta_contraction_opus_20260728.out"
ATOM_SHA256 = "f89be10c65bb77270199f9399b155d5a2c82c0da121b3e8589fe3c1f7e9824fc"
R5_OUTPUT_SHA256 = "bef4ee9a18ff3e2f455bad66a95252dd9989b2f60953e26e8ea0c2dc6ae7f5df"
EXPECTED_FULL_KERNEL_SHA256 = "9f03166d1b8d730673aeccf7bff9196d2de42e0ff9fa7c08680d38a772dc041f"
EXPECTED_AUGMENTATION_KERNEL_SHA256 = "351976b6d1ac51915d686b32e1cc0f3f3dfa53d1f33f7c8831f7a74a495e3025"
EXPECTED_SEMANTIC_SHA256 = "a6d6ad7a9891af8bddcba6fe7cf9f2554214731043844098f8556d1003885222"

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
CORNERS = (
    ("left", "left"),
    ("left", "right"),
    ("right", "left"),
    ("right", "right"),
)
CHANNELS = ("constant", "left", "right", "mixed")
Q_H = (1, 0, 1)
Q_Q5 = (1, 0, 0)
EDGES = (
    (0, 3), (0, 4), (0, 5),
    (1, 2), (1, 4), (1, 7),
    (2, 4), (2, 7),
    (3, 4), (3, 5),
    (4, 5), (4, 6), (4, 7),
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def json_digest(value: object) -> str:
    body = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return hashlib.sha256(body).hexdigest()


def load_atom_audit():
    require(lf_sha256(ATOM_PATH) == ATOM_SHA256, "THM-3514 atom audit drift")
    name = "thm3514_atoms_for_folded_transporter_audit"
    spec = importlib.util.spec_from_file_location(name, ATOM_PATH)
    require(spec is not None and spec.loader is not None, "atom loader")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


A = load_atom_audit()


def atom_worker(alpha: int) -> tuple[object, ...]:
    return A.worker(alpha)


def parse_source_table() -> tuple[tuple[tuple[int, ...], ...], int]:
    """Read the proved THM-2594 raw aggregate, without its Fourier helper."""
    require(lf_sha256(R5_OUTPUT) == R5_OUTPUT_SHA256, "THM-2594 output drift")
    lines = R5_OUTPUT.read_text(encoding="utf-8").splitlines()
    denominator_line = next(line for line in lines if "entry denominators: DENC" in line)
    denominator = int(denominator_line.rsplit("=", 1)[1].strip())
    start = next(index for index, line in enumerate(lines) if line.startswith("[8] the response array"))
    stop = next(index for index in range(start + 1, len(lines)) if "fibrewise ANOVA" in lines[index])
    rows: dict[int, tuple[int, ...]] = {}
    for line in lines[start + 1 : stop]:
        stripped = line.strip()
        if not stripped.startswith("ell="):
            continue
        label, payload = stripped.split(": theta=0..2: ", 1)
        values_text = payload.split("  |", 1)[0]
        ell = int(label.split("=", 1)[1])
        values = ast.literal_eval(values_text)
        require(isinstance(values, list) and len(values) == 3, (ell, values))
        rows[ell] = tuple(int(value) for value in values) + (0,) * 10
    require(tuple(sorted(rows)) == tuple(range(7)), sorted(rows))
    return tuple(rows[index] for index in range(7)), denominator


def dft(values: tuple[int, ...], frequency: int, root: int, prime: int) -> int:
    return sum(
        value * pow(root, (-frequency * index) % len(values), prime)
        for index, value in enumerate(values)
    ) % prime


def spectrum(values: tuple[int, ...], root: int, prime: int) -> tuple[int, ...]:
    return tuple(dft(values, frequency, root, prime) for frequency in range(len(values)))


def inverse_spectrum(values: tuple[int, ...], root: int, prime: int) -> tuple[int, ...]:
    inverse = pow(len(values), -1, prime)
    return tuple(
        inverse
        * sum(
            values[frequency] * pow(root, (frequency * index) % len(values), prime)
            for frequency in range(len(values))
        )
        % prime
        for index in range(len(values))
    )


def convolution(left: tuple[int, ...], right: tuple[int, ...], prime: int) -> tuple[int, ...]:
    return tuple(
        sum(left[shift] * right[(index - shift) % len(left)] for shift in range(len(left)))
        % prime
        for index in range(len(left))
    )


def centered(values: tuple[int, ...], prime: int) -> tuple[int, ...]:
    mean = sum(values) * pow(len(values), -1, prime) % prime
    return tuple((value - mean) % prime for value in values)


def rref(matrix: tuple[tuple[int, ...], ...], prime: int) -> tuple[tuple[tuple[int, ...], ...], tuple[int, ...]]:
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
    columns = len(matrix[0])
    free = tuple(column for column in range(columns) if column not in pivots)
    basis = []
    for free_column in free:
        vector = [0] * columns
        vector[free_column] = 1
        for row, pivot in enumerate(pivots):
            vector[pivot] = (-reduced[row][free_column]) % prime
        basis.append(tuple(vector))
    return tuple(basis)


def det_mod(matrix: tuple[tuple[int, ...], ...], prime: int) -> int:
    rows = [list(value % prime for value in row) for row in matrix]
    size = len(rows)
    require(size and all(len(row) == size for row in rows), "square determinant")
    answer = 1
    for column in range(size):
        pivot = next((row for row in range(column, size) if rows[row][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            answer = -answer
        value = rows[column][column]
        answer = answer * value % prime
        inverse = pow(value, -1, prime)
        for row in range(column + 1, size):
            factor = rows[row][column] * inverse % prime
            if factor:
                rows[row] = [
                    (entry - factor * pivot_entry) % prime
                    for entry, pivot_entry in zip(rows[row], rows[column])
                ]
    return answer % prime


def det3_integer(matrix: tuple[tuple[int, int, int], ...]) -> int:
    (a, b, c), (d, e, f), (g, h, i) = matrix
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def first_modular_minor(
    matrix: tuple[tuple[int, ...], ...], size: int, prime: int
) -> tuple[tuple[int, ...], tuple[int, ...], int]:
    for row_choice in combinations(range(len(matrix)), size):
        for column_choice in combinations(range(len(matrix[0])), size):
            minor = tuple(
                tuple(matrix[row][column] for column in column_choice)
                for row in row_choice
            )
            determinant = det_mod(minor, prime)
            if determinant:
                return row_choice, column_choice, determinant
    raise RuntimeError("missing nonzero minor")


def materialize_target() -> tuple[object, ...]:
    """Direct tau-slice reconstruction from the independent THM-3514 atoms."""
    require(tuple(A.ATOMS) == ATOMS, "atom order")
    require(tuple(A.BUCKETS) == BUCKETS, "bucket order")
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(atom_worker, range(N)))
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(N)), "worker order")
    tables = tuple(row for chunk in chunks for row in chunk[1])
    require(len(tables) == N * N, len(tables))
    (_word, _t_den, nn, prime, root, zeta13, *_rest) = A.context()
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
                    if value and A.safe(ATOMS[index][1], ATOMS[index][0] + tau)
                )
                right_atoms = tuple(
                    (index, value)
                    for index, value in enumerate(by_values)
                    if value and A.safe(ATOMS[index][1], ATOMS[index][0] + tau)
                )
                slices = [0] * len(BUCKETS)
                for left_index, left_value in left_atoms:
                    left_sheet, left_chamber = ATOMS[left_index]
                    for right_index, right_value in right_atoms:
                        right_sheet, right_chamber = ATOMS[right_index]
                        address = (
                            left_chamber,
                            right_chamber,
                            (right_sheet - left_sheet) % N,
                        )
                        bucket = BUCKET_INDEX[address]
                        slices[bucket] = (slices[bucket] + left_value * right_value) % prime
                for q in (Q_H, Q_Q5):
                    qa, qb, qt = q
                    phase = pow(zeta13, (beta - alpha * qa - beta * qb - tau * qt) % N, prime)
                    for bucket, value in enumerate(slices):
                        raw[q][bucket] = (raw[q][bucket] + phase * value) % prime
    normalizer = pow(N**3, -1, prime)
    buckets = {
        q: tuple(value * normalizer % prime for value in raw[q])
        for q in (Q_H, Q_Q5)
    }
    bridge = tuple((left - right) % prime for left, right in zip(buckets[Q_H], buckets[Q_Q5]))
    corner_rows = tuple(
        tuple(bridge[BUCKET_INDEX[(left, right, drift)]] for drift in range(N))
        for left, right in CORNERS
    )
    walsh = (
        tuple(sum(corner_rows[row][drift] for row in range(4)) % prime for drift in range(N)),
        tuple((corner_rows[0][d] + corner_rows[1][d] - corner_rows[2][d] - corner_rows[3][d]) % prime for d in range(N)),
        tuple((corner_rows[0][d] - corner_rows[1][d] + corner_rows[2][d] - corner_rows[3][d]) % prime for d in range(N)),
        tuple((corner_rows[0][d] - corner_rows[1][d] - corner_rows[2][d] + corner_rows[3][d]) % prime for d in range(N)),
    )
    target_spectra = tuple(spectrum(row, zeta13, prime) for row in walsh)
    require(all(all(row) for row in walsh), "target physical zero")
    require(all(all(row) for row in target_spectra), "target spectral zero")
    require(rank_mod(walsh, prime) == 4, "target rank")
    require(sum(walsh[0]) % prime == 389266878372286537904, "bridge regression")
    return nn, prime, root, zeta13, buckets, walsh, target_spectra


def primitive_root(prime: int) -> int:
    order = prime - 1
    factors = []
    value = order
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            factors.append(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 1
    if value > 1:
        factors.append(value)
    for candidate in range(2, prime):
        if all(pow(candidate, order // factor, prime) != 1 for factor in factors):
            return candidate
    raise RuntimeError(("primitive root", prime))


def source_at_prime(
    raw: tuple[tuple[int, ...], ...], denominator: int, prime: int
) -> tuple[int, int]:
    require((prime - 1) % 91 == 0 and denominator % prime, (prime, denominator % prime))
    generator = primitive_root(prime)
    zeta91 = pow(generator, (prime - 1) // 91, prime)
    zeta7 = pow(zeta91, 13, prime)
    zeta13 = pow(zeta91, 7, prime)
    require(pow(zeta7, 7, prime) == 1 and zeta7 != 1, (prime, "zeta7"))
    require(pow(zeta13, 13, prime) == 1 and zeta13 != 1, (prime, "zeta13"))
    inverse_denominator = pow(denominator, -1, prime)
    table = tuple(
        tuple(value % prime * inverse_denominator % prime for value in row)
        for row in raw
    )
    transformed = tuple(
        tuple(
            sum(
                table[ell][theta]
                * pow(zeta7, (-septimal * ell) % 7, prime)
                * pow(zeta13, (-frequency * theta) % 13, prime)
                for ell in range(7)
                for theta in range(13)
            )
            % prime
            for frequency in range(13)
        )
        for septimal in range(7)
    )
    return sum(value != 0 for row in transformed for value in row), rank_mod(table, prime)


def folded_algebra(zeta7: int, prime: int) -> tuple[object, ...]:
    residues = {1, 2, 4}
    arcs = frozenset(
        (left, right)
        for left in range(7)
        for right in range(7)
        if left != right and (right - left) % 7 in residues
    )
    orbits = ((0,), (1, 6), (2, 5), (3, 4))
    counts = tuple(
        tuple(sum((left, right) in arcs for left in source for right in target) for target in orbits)
        for source in orbits
    )
    require(counts == ((0, 1, 1, 1), (1, 1, 2, 2), (1, 2, 1, 2), (1, 2, 2, 1)), counts)
    require(all(counts[i][j] == counts[j][i] for i in range(4) for j in range(4)), counts)
    times_two = tuple(
        next(index for index, orbit in enumerate(orbits) if (2 * source[0]) % 7 in orbit)
        for source in orbits[1:]
    )
    require(times_two == (2, 3, 1), times_two)
    vertices = ((0, 0), (0, 1), (1, 0), (1, 1))
    matchings = []
    for character in ((1, 0), (0, 1), (1, 1)):
        matching = tuple(
            (left, right)
            for left in range(4)
            for right in range(left + 1, 4)
            if sum(character[i] * vertices[left][i] for i in range(2)) % 2
            == sum(character[i] * vertices[right][i] for i in range(2)) % 2
        )
        matchings.append(matching)
    require(tuple(matchings) == (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))), matchings)
    discriminant = 1 * 1 * (-2) ** 2 - 4 * (-2) ** 3 - 4 * (-1) - 27 + 18 * 1 * 1 * (-2) * (-1)
    require(discriminant == 49, discriminant)
    require(all(value**3 + value**2 - 2 * value - 1 for value in (-1, 1)), "cubic reducible")
    roots = tuple((pow(zeta7, exponent, prime) + pow(zeta7, -exponent, prime)) % prime for exponent in (1, 2, 3))
    require(len(set(roots)) == 3, roots)
    require(all((value**3 + value**2 - 2 * value - 1) % prime == 0 for value in roots), roots)
    return orbits, counts, times_two, tuple(matchings), discriminant, roots


def fold_class(mode: int) -> int:
    return min(mode, (-mode) % 7)


def ratio_row(
    target: tuple[int, ...], source: tuple[int, ...], prime: int
) -> tuple[int, ...]:
    require(all(target) and all(source), "ratio requires cyclic rows")
    return tuple(left * pow(right, -1, prime) % prime for left, right in zip(target, source))


def gauge_signature(row: tuple[int, ...], zeta13: int, prime: int) -> tuple[int, ...]:
    normalized = tuple(value * pow(row[0], -1, prime) % prime for value in row)
    return min(
        tuple(
            normalized[frequency] * pow(zeta13, (frequency * shift) % N, prime) % prime
            for frequency in range(N)
        )
        for shift in range(N)
    )


def allocation_census(
    allocations: tuple[tuple[int, ...], ...],
    source_spectra: tuple[tuple[int, ...], ...],
    target_spectra: tuple[tuple[int, ...], ...],
    zeta13: int,
    prime: int,
) -> tuple[int, int, tuple[int, ...]]:
    exact = 0
    gauged = 0
    class_counts = set()
    for allocation in allocations:
        ratios = tuple(
            ratio_row(target_spectra[channel], source_spectra[allocation[channel]], prime)
            for channel in range(4)
        )
        exact += len(set(ratios)) == 1
        signatures = tuple(gauge_signature(row, zeta13, prime) for row in ratios)
        classes = len(set(signatures))
        class_counts.add(classes)
        gauged += classes == 1
    return exact, gauged, tuple(sorted(class_counts))


def projective_anchor_system(
    source: tuple[tuple[int, ...], ...],
    target: tuple[tuple[int, ...], ...],
    prime: int,
) -> tuple[tuple[int, ...], ...]:
    equations = []
    for frequency in range(N):
        require(target[0][frequency], ("target anchor", frequency))
        for output in range(1, 4):
            row = [0] * 16
            for channel in range(4):
                row[output * 4 + channel] = target[0][frequency] * source[channel][frequency] % prime
                row[channel] = (row[channel] - target[output][frequency] * source[channel][frequency]) % prime
            equations.append(tuple(row))
    return tuple(equations)


def augmented_transport_system(
    source: tuple[tuple[int, ...], ...],
    target: tuple[tuple[int, ...], ...],
    prime: int,
) -> tuple[tuple[int, ...], ...]:
    equations = []
    for frequency in range(N):
        for output in range(4):
            row = [0] * (16 + N)
            for channel in range(4):
                row[output * 4 + channel] = source[channel][frequency]
            row[16 + frequency] = (-target[output][frequency]) % prime
            equations.append(tuple(row))
    return tuple(equations)


def incidence_cycle_certificate(prime: int) -> tuple[int, int]:
    incidence = [[0] * len(EDGES) for _ in range(8)]
    for column, (left, right) in enumerate(EDGES):
        incidence[left][column] = 1
        incidence[right][column] = -1 % prime
    matrix = tuple(tuple(row) for row in incidence)
    rank = rank_mod(matrix, prime)
    cycles = nullspace(matrix, prime)
    require(rank == 7 and len(cycles) == 6, (rank, len(cycles)))
    require(
        all(
            sum(matrix[vertex][edge] * cycle[edge] for edge in range(len(EDGES))) % prime == 0
            for cycle in cycles
            for vertex in range(8)
        ),
        "cycle basis",
    )
    return rank, len(cycles)


def main() -> None:
    raw_integer, denominator = parse_source_table()
    source_integer_minor = next(
        (rows, det3_integer(tuple(tuple(raw_integer[row][column] for column in range(3)) for row in rows)))
        for rows in combinations(range(7), 3)
        if det3_integer(tuple(tuple(raw_integer[row][column] for column in range(3)) for row in rows))
    )
    require(all(value == 0 for row in raw_integer for value in row[3:]), "source support")

    nn, prime, root, zeta13, target_buckets, walsh, target_spectra = materialize_target()
    require(nn % 7 == 0, nn)
    zeta7 = pow(root, nn // 7, prime)
    require(pow(zeta7, 7, prime) == 1 and zeta7 != 1, "zeta7")
    inverse_denominator = pow(denominator, -1, prime)
    raw = tuple(
        tuple(value % prime * inverse_denominator % prime for value in row)
        for row in raw_integer
    )
    source_profiles = tuple(
        tuple(
            sum(raw[ell][theta] * pow(zeta7, (-septimal * ell) % 7, prime) for ell in range(7)) % prime
            for theta in range(13)
        )
        for septimal in range(7)
    )
    source_spectra = tuple(spectrum(row, zeta13, prime) for row in source_profiles)
    direct_two_dimensional = tuple(
        tuple(
            sum(
                raw[ell][theta]
                * pow(zeta7, (-septimal * ell) % 7, prime)
                * pow(zeta13, (-frequency * theta) % 13, prime)
                for ell in range(7)
                for theta in range(13)
            )
            % prime
            for frequency in range(13)
        )
        for septimal in range(7)
    )
    require(source_spectra == direct_two_dimensional, "separable/direct DFT")
    require(all(all(row) for row in source_spectra), "source spectral zero")
    source_rank = rank_mod(source_profiles, prime)
    target_rank = rank_mod(walsh, prime)
    require((source_rank, target_rank) == (3, 4), (source_rank, target_rank))
    source_minor = first_modular_minor(source_profiles, 3, prime)
    target_minor = first_modular_minor(walsh, 4, prime)
    small_prime_checks = tuple(
        (small_prime,) + source_at_prime(raw_integer, denominator, small_prime)
        for small_prime in (547, 1093, 2003)
    )
    require(all(item[1:] == (91, 3) for item in small_prime_checks), small_prime_checks)

    paley = folded_algebra(zeta7, prime)
    all_choices = tuple(product(range(7), repeat=4))
    preserving = tuple(
        choice
        for choice in all_choices
        if choice[0] == 0 and tuple(sorted(fold_class(mode) for mode in choice)) == (0, 1, 2, 3)
    )
    folded = tuple(
        choice
        for choice in all_choices
        if tuple(sorted(fold_class(mode) for mode in choice)) == (0, 1, 2, 3)
    )
    require((len(preserving), len(folded), len(all_choices)) == (48, 192, 2401), "allocation universes")
    allocation_censuses = (
        allocation_census(preserving, source_spectra, target_spectra, zeta13, prime),
        allocation_census(folded, source_spectra, target_spectra, zeta13, prime),
        allocation_census(all_choices, source_spectra, target_spectra, zeta13, prime),
    )
    require(all(census == (0, 0, (4,)) for census in allocation_censuses), allocation_censuses)

    selected_ranks = []
    projective_ranks = []
    augmented_ranks = []
    augmented_nullities = []
    multiplier_projection_ranks = []
    for signs in product((0, 1), repeat=3):
        selected_modes = (0,) + tuple(((1, 6), (2, 5), (3, 4))[index][signs[index]] for index in range(3))
        selected = tuple(source_spectra[mode] for mode in selected_modes)
        selected_ranks.append(rank_mod(selected, prime))
        projective = projective_anchor_system(selected, target_spectra, prime)
        projective_ranks.append(rank_mod(projective, prime))
        augmented = augmented_transport_system(selected, target_spectra, prime)
        augmented_ranks.append(rank_mod(augmented, prime))
        basis = nullspace(augmented, prime)
        augmented_nullities.append(len(basis))
        multiplier_projection_ranks.append(rank_mod(tuple(vector[16:] for vector in basis), prime))
        require(
            all(
                sum(vector[output * 4 + channel] * selected[channel][frequency] for channel in range(4)) % prime == 0
                for vector in basis
                for output in range(4)
                for frequency in range(N)
            ),
            ("non-annihilator basis", signs),
        )
    require(tuple(selected_ranks) == (3,) * 8, selected_ranks)
    require(tuple(projective_ranks) == (12,) * 8, projective_ranks)
    require(tuple(augmented_ranks) == (25,) * 8, augmented_ranks)
    require(tuple(augmented_nullities) == (4,) * 8, augmented_nullities)
    require(tuple(multiplier_projection_ranks) == (0,) * 8, multiplier_projection_ranks)

    calibration = (0, 1, 2, 3)
    ratios = tuple(
        ratio_row(target_spectra[channel], source_spectra[mode], prime)
        for channel, mode in enumerate(calibration)
    )
    full_kernels = tuple(inverse_spectrum(row, zeta13, prime) for row in ratios)
    require(
        all(convolution(full_kernels[channel], source_profiles[mode], prime) == walsh[channel] for channel, mode in enumerate(calibration)),
        "full circulants",
    )
    full_digest = json_digest(full_kernels)
    require(full_digest == EXPECTED_FULL_KERNEL_SHA256, full_digest)
    augmentation_kernels = []
    margin_multipliers = []
    for channel, mode in enumerate(calibration):
        source_centered = centered(source_profiles[mode], prime)
        target_centered = centered(walsh[channel], prime)
        source_hat = spectrum(source_centered, zeta13, prime)
        target_hat = spectrum(target_centered, zeta13, prime)
        require(source_hat[0] == target_hat[0] == 0, "centred zero")
        multiplier = (0,) + tuple(
            target_hat[frequency] * pow(source_hat[frequency], -1, prime) % prime
            for frequency in range(1, N)
        )
        require(multiplier[1:] == ratios[channel][1:], "centred ratio")
        kernel = inverse_spectrum(multiplier, zeta13, prime)
        require(convolution(kernel, source_centered, prime) == target_centered, "augmentation convolution")
        difference = tuple((full_kernels[channel][index] - kernel[index]) % prime for index in range(N))
        require(len(set(difference)) == 1, "margin is not constant")
        augmentation_kernels.append(kernel)
        margin_multipliers.append(ratios[channel][0])
    augmentation_kernels_tuple = tuple(augmentation_kernels)
    augmentation_digest = json_digest(augmentation_kernels_tuple)
    require(augmentation_digest == EXPECTED_AUGMENTATION_KERNEL_SHA256, augmentation_digest)

    incidence_rank, cycle_rank = incidence_cycle_certificate(prime)
    semantic = (
        ATOM_SHA256,
        R5_OUTPUT_SHA256,
        (prime, root, nn, zeta7, zeta13),
        raw_integer,
        denominator,
        source_integer_minor,
        source_profiles,
        source_spectra,
        target_buckets,
        walsh,
        target_spectra,
        small_prime_checks,
        paley,
        (len(preserving), len(folded), len(all_choices)),
        allocation_censuses,
        tuple(selected_ranks),
        tuple(projective_ranks),
        tuple(augmented_ranks),
        tuple(augmented_nullities),
        tuple(multiplier_projection_ranks),
        source_minor,
        target_minor,
        full_kernels,
        augmentation_kernels_tuple,
        tuple(margin_multipliers),
        (incidence_rank, cycle_rank),
        "THM-3518 phase is target-side; common right actions preserve rank three; cut cycles vanish",
    )
    semantic_hash = hashlib.sha256(repr(semantic).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, (semantic_hash, EXPECTED_SEMANTIC_SHA256))

    print("FOLDED-C7/K4 DRIFT TRANSPORTER -- INDEPENDENT AUDIT")
    print("status=ACCEPT_SCOPED_FINITE_EXACT; endpoint/split-field only; LRC(14) OPEN")
    print(f"dependencies=(thm3514_atom_audit={ATOM_SHA256},thm2594_output={R5_OUTPUT_SHA256})")
    print(f"field=(prime={prime},zeta7={zeta7},zeta13={zeta13})")
    print(f"source_exact_rank=(integer_minor_rows={source_integer_minor[0]},integer_det={source_integer_minor[1]},split_rank={source_rank},split_minor={source_minor})")
    print(f"source_independent_split_primes={small_prime_checks}")
    print(f"source_spectrum=(nonzero=91/91,separable_vs_direct_2d=PASS); target_spectrum=(nonzero=52/52,rank={target_rank},minor={target_minor})")
    print(f"paley_fold=(orbits={paley[0]},arc_counts={paley[1]},times2={paley[2]},walsh_matchings={paley[3]})")
    print(f"rational_descent=(real_cubic_disc={paley[4]},split_roots={paley[5]},QxC3field_not_Q4=True)")
    print(f"allocation_universes=(predicate_filtered_48={len(preserving)},predicate_filtered_192={len(folded)},all_7pow4={len(all_choices)})")
    print(f"allocation_censuses=(preserving={allocation_censuses[0]},folded={allocation_censuses[1]},unrestricted={allocation_censuses[2]})")
    print(f"fixed_mixing=(selected_source_ranks={tuple(selected_ranks)},projective_anchor_ranks={tuple(projective_ranks)},augmented_M_lambda_ranks={tuple(augmented_ranks)},nullities={tuple(augmented_nullities)},lambda_projection_ranks={tuple(multiplier_projection_ranks)})")
    print("rank_obstruction=all augmented solutions have lambda_k=0 and annihilate the rank-three source; arbitrary 4x7 fixed mixing plus any common 13x13 right operator also has rank at most three")
    print(f"formal_survivor=(full_kernel_supports={tuple(sum(value != 0 for value in row) for row in full_kernels)},rank={rank_mod(full_kernels, prime)},sha256={full_digest})")
    print(f"augmentation_survivor=(supports={tuple(sum(value != 0 for value in row) for row in augmentation_kernels_tuple)},rank={rank_mod(augmentation_kernels_tuple, prime)},sha256={augmentation_digest},margin_multipliers={tuple(margin_multipliers)})")
    print(f"thm3518_sidecar=(incidence_rank={incidence_rank},cycle_rank={cycle_rank},cycle_coordinates=all_zero,common_phase_rank_bound={source_rank},lawful_source_channels_added=0)")
    print("typing=the q_t phase is defined on U_full endpoint guard sheets, not the THM-2594 ancestry table; a common phase is rank-preserving, while channel-dependent phases leave the audited transporter class and still supply no common-stalk map")
    print(f"semantic_sha256={semantic_hash}")
    print("nonconsequence=no ancestry transplant, Boolean coupling, physical current, absolute H1, bispectrum, scalar exclusion, or LRC(14)")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
