#!/usr/bin/env python3
"""Primary exact referee for the THM-4099 OCF-curvature candidate.

This program is deliberately self-contained and theorem-ID neutral.  It checks
only finite exact identities and explicit algebraic criteria:

* OCF union-support atoms and Boolean Möbius nonnegativity/parity;
* tensor powers of the directed-triangle face table;
* a strict defect-degree delay;
* explicit Rayleigh and quadratic-Hessian hostiles;
* the directed-cut interpretation of one-ear response and its opposite
  Boolean curvature.

Tournament codes use lexicographic unordered pairs: bit one on ``(i,j)``,
``i<j``, means ``i -> j``.  All executable gates use ``require`` rather than
``assert`` so normal and ``python -O`` runs exercise identical checks.  There
is no floating-point arithmetic and no repository import.
"""

from __future__ import annotations

from hashlib import sha256
from itertools import combinations, permutations
import json


EXPECTED_SEMANTIC_SHA256 = "0eeda799bfc4764678b3aa8c5067b04b7a74944e5b0b6dbf519f427b09dd56b5"


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError("FAILED: " + label)


def decode_raw(code: int, n: int) -> tuple[int, ...]:
    require(code >= 0, "nonnegative tournament code")
    adjacency = [0] * n
    bit_index = 0
    for left in range(n):
        for right in range(left + 1, n):
            if (code >> bit_index) & 1:
                adjacency[left] |= 1 << right
            else:
                adjacency[right] |= 1 << left
            bit_index += 1
    require(code >> bit_index == 0, "tournament code fits declared order")
    require(tournament_well_formed(tuple(adjacency)), "decoded tournament")
    return tuple(adjacency)


def tournament_well_formed(adjacency: tuple[int, ...]) -> bool:
    n = len(adjacency)
    full = (1 << n) - 1
    for vertex in range(n):
        if adjacency[vertex] & ~full:
            return False
        if (adjacency[vertex] >> vertex) & 1:
            return False
    for left in range(n):
        for right in range(left + 1, n):
            if bool((adjacency[left] >> right) & 1) == bool(
                (adjacency[right] >> left) & 1
            ):
                return False
    return True


def arc(adjacency: tuple[int, ...], source: int, target: int) -> bool:
    return bool((adjacency[source] >> target) & 1)


def all_induced_h(adjacency: tuple[int, ...]) -> tuple[int, ...]:
    """Hamilton-path count on every induced mask, with H(empty)=1."""
    require(tournament_well_formed(adjacency), "induced-H tournament input")
    n = len(adjacency)
    size = 1 << n
    dp = [[0] * n for _ in range(size)]
    h = [0] * size
    h[0] = 1
    for vertex in range(n):
        dp[1 << vertex][vertex] = 1
        h[1 << vertex] = 1
    for subset in range(1, size):
        if subset & (subset - 1) == 0:
            continue
        total = 0
        vertices = subset
        while vertices:
            bit = vertices & -vertices
            last = bit.bit_length() - 1
            vertices -= bit
            rest = subset ^ bit
            predecessors = rest
            count = 0
            while predecessors:
                pred_bit = predecessors & -predecessors
                predecessor = pred_bit.bit_length() - 1
                predecessors -= pred_bit
                if arc(adjacency, predecessor, last):
                    count += dp[rest][predecessor]
            dp[subset][last] = count
            total += count
        h[subset] = total
    return tuple(h)


def hamilton_literal(adjacency: tuple[int, ...]) -> int:
    require(tournament_well_formed(adjacency), "literal-H tournament input")
    n = len(adjacency)
    return sum(
        all(arc(adjacency, word[j], word[j + 1]) for j in range(n - 1))
        for word in permutations(range(n))
    )


def directed_odd_cycles(adjacency: tuple[int, ...]) -> tuple[int, ...]:
    """Cycle masks, once per cyclic ordering by rooting at its least vertex."""
    n = len(adjacency)
    cycles: list[int] = []
    for length in range(3, n + 1, 2):
        for vertices in combinations(range(n), length):
            root = vertices[0]
            rest = vertices[1:]
            mask = sum(1 << vertex for vertex in vertices)
            for ordering in permutations(rest):
                cycle = (root,) + ordering
                if all(
                    arc(adjacency, cycle[j], cycle[(j + 1) % length])
                    for j in range(length)
                ):
                    cycles.append(mask)
    return tuple(cycles)


def packing_union_atoms_through_six(
    cycles: tuple[int, ...], n: int
) -> tuple[int, ...]:
    """OCF weight by union mask; at n<=6 at most two cycles are disjoint."""
    require(n <= 6, "two-cycle packing cutoff")
    atoms = [0] * (1 << n)
    atoms[0] = 1
    for cycle in cycles:
        atoms[cycle] += 2
    for first_index, first in enumerate(cycles):
        for second in cycles[first_index + 1 :]:
            if not (first & second):
                atoms[first | second] += 4
    return tuple(atoms)


def boolean_zeta(values: tuple[int, ...], rank: int) -> tuple[int, ...]:
    out = list(values)
    for bit_index in range(rank):
        bit = 1 << bit_index
        for mask in range(1 << rank):
            if mask & bit:
                out[mask] += out[mask ^ bit]
    return tuple(out)


def boolean_mobius(values: tuple[int, ...], rank: int) -> tuple[int, ...]:
    out = list(values)
    for bit_index in range(rank):
        bit = 1 << bit_index
        for mask in range(1 << rank):
            if mask & bit:
                out[mask] -= out[mask ^ bit]
    return tuple(out)


def local_mask_to_global(mask: int, vertices: tuple[int, ...]) -> int:
    return sum(
        1 << vertices[index]
        for index in range(len(vertices))
        if (mask >> index) & 1
    )


def global_mask_to_local(mask: int, vertices: tuple[int, ...]) -> int:
    return sum(
        1 << index
        for index, vertex in enumerate(vertices)
        if (mask >> vertex) & 1
    )


def submasks(mask: int):
    subset = mask
    while True:
        yield subset
        if subset == 0:
            break
        subset = (subset - 1) & mask


def mixed_difference(face: tuple[int, ...], active: int, context: int) -> int:
    require(not (active & context), "disjoint mixed-difference masks")
    total = 0
    active_size = active.bit_count()
    for subset in submasks(active):
        sign = -1 if (active_size - subset.bit_count()) & 1 else 1
        total += sign * face[context | subset]
    return total


def probe_ocf_and_mobius() -> dict[str, object]:
    tournament_count = 0
    induced_mask_count = 0
    cycle_instance_count = 0
    parity_checks = 0
    base_splits = 0
    atom_checks = 0
    mixed_checks = 0
    min_mixed: int | None = None
    max_atom = 0

    for n in range(1, 7):
        code_count = 1 << (n * (n - 1) // 2)
        for code in range(code_count):
            adjacency = decode_raw(code, n)
            h = all_induced_h(adjacency)
            cycles = directed_odd_cycles(adjacency)
            union_atoms = packing_union_atoms_through_six(cycles, n)
            ocf = boolean_zeta(union_atoms, n)
            require(h == ocf, f"OCF induced table n={n} code={code}")
            require(all(value > 0 and value & 1 for value in h),
                    f"positive odd induced H n={n} code={code}")

            tournament_count += 1
            induced_mask_count += 1 << n
            cycle_instance_count += len(cycles)
            parity_checks += len(h)

            if n > 5:
                continue

            full = (1 << n) - 1
            for base in range(1, 1 << n):
                inserted = tuple(
                    vertex for vertex in range(n) if not ((base >> vertex) & 1)
                )
                rank = len(inserted)
                face = tuple(
                    h[base | local_mask_to_global(mask, inserted)]
                    for mask in range(1 << rank)
                )
                mobius = boolean_mobius(face, rank)
                exact = [0] * (1 << rank)
                for union, weight in enumerate(union_atoms):
                    if weight:
                        support = global_mask_to_local(union, inserted)
                        exact[support] += weight
                require(mobius == tuple(exact),
                        f"Möbius exact-support atoms n={n} code={code} base={base}")
                require(mobius[0] > 0 and mobius[0] & 1,
                        "empty-support atom positive odd")
                require(all(value >= 0 for value in mobius),
                        "Möbius atoms nonnegative")
                require(all(value % 2 == 0 for value in mobius[1:]),
                        "nonempty-support atoms even")

                base_splits += 1
                atom_checks += len(mobius)
                max_atom = max(max_atom, *mobius)

                local_full = (1 << rank) - 1
                for active in range(1 << rank):
                    available = local_full ^ active
                    for context in submasks(available):
                        observed = mixed_difference(face, active, context)
                        expected = sum(
                            mobius[active | extra]
                            for extra in submasks(context)
                        )
                        require(observed == expected,
                                "mixed difference equals supported-atom sum")
                        require(observed >= 0, "mixed difference nonnegative")
                        min_mixed = (
                            observed
                            if min_mixed is None
                            else min(min_mixed, observed)
                        )
                        mixed_checks += 1

    require(tournament_count == 33_867, "OCF tournament universe")
    require(induced_mask_count == 2_131_018, "OCF induced-mask universe")
    require(cycle_instance_count == 314_690, "directed-cycle universe")
    require(parity_checks == 2_131_018, "face-coefficient parity universe")
    require(base_splits == 32_767, "Möbius base-split universe")
    require(atom_checks == 220_387, "Möbius coefficient universe")
    require(mixed_checks == 811_255, "mixed-difference universe")
    require(min_mixed == 0, "mixed-difference minimum")

    return {
        "tournaments_n1_to_n6": tournament_count,
        "induced_masks": induced_mask_count,
        "directed_cycle_instances": cycle_instance_count,
        "parity_checks": parity_checks,
        "mobius_base_splits_n1_to_n5": base_splits,
        "mobius_atom_checks": atom_checks,
        "mixed_difference_checks": mixed_checks,
        "minimum_mixed_difference": min_mixed,
        "maximum_atom": max_atom,
    }


def bareiss_determinant(matrix: tuple[tuple[int, ...], ...]) -> int:
    n = len(matrix)
    require(n > 0 and all(len(row) == n for row in matrix),
            "square determinant matrix")
    work = [list(row) for row in matrix]
    sign = 1
    previous = 1
    for pivot_index in range(n - 1):
        if work[pivot_index][pivot_index] == 0:
            swap = next(
                (
                    row
                    for row in range(pivot_index + 1, n)
                    if work[row][pivot_index] != 0
                ),
                None,
            )
            if swap is None:
                return 0
            work[pivot_index], work[swap] = work[swap], work[pivot_index]
            sign *= -1
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, n):
            for column in range(pivot_index + 1, n):
                numerator = (
                    work[row][column] * pivot
                    - work[row][pivot_index] * work[pivot_index][column]
                )
                require(numerator % previous == 0, "Bareiss exact division")
                work[row][column] = numerator // previous
        for row in range(pivot_index + 1, n):
            work[row][pivot_index] = 0
        previous = pivot
    return sign * work[-1][-1]


def rank_mod(matrix: tuple[tuple[int, ...], ...], prime: int = 1_000_003) -> int:
    require(matrix and all(len(row) == len(matrix[0]) for row in matrix),
            "rectangular rank matrix")
    work = [[entry % prime for entry in row] for row in matrix]
    row_count = len(work)
    column_count = len(work[0])
    rank = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(rank, row_count) if work[row][column]), None
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], prime - 2, prime)
        work[rank] = [(entry * inverse) % prime for entry in work[rank]]
        for row in range(row_count):
            if row == rank or not work[row][column]:
                continue
            multiple = work[row][column]
            work[row] = [
                (left - multiple * right) % prime
                for left, right in zip(work[row], work[rank])
            ]
        rank += 1
        if rank == min(row_count, column_count):
            break
    return rank


def c3_tensor_matrix(k: int) -> tuple[tuple[int, ...], ...]:
    size = 1 << k
    return tuple(
        tuple(3 ** ((left & right).bit_count()) for right in range(size))
        for left in range(size)
    )


def transitive_c3_join(k: int) -> tuple[int, ...]:
    n = 3 * k
    adjacency = [0] * n
    for gadget in range(k):
        base = 3 * gadget
        u = base + 1
        v = base + 2
        adjacency[base] |= 1 << u
        adjacency[u] |= 1 << v
        adjacency[v] |= 1 << base
    for first in range(k):
        for second in range(first + 1, k):
            for source in range(3 * first, 3 * first + 3):
                for target in range(3 * second, 3 * second + 3):
                    adjacency[source] |= 1 << target
    out = tuple(adjacency)
    require(tournament_well_formed(out), "transitive C3 join tournament")
    return out


def probe_c3_tensors() -> dict[str, object]:
    rows = []
    actual_face_checks = 0
    for k in range(1, 7):
        matrix = c3_tensor_matrix(k)
        size = 1 << k
        expected_determinant = 2 ** (k * (1 << (k - 1)))
        require(rank_mod(matrix) == size, f"C3 tensor full rank k={k}")
        require(
            {entry for row in matrix for entry in row}
            == {3**power for power in range(k + 1)},
            f"C3 tensor sparse range k={k}",
        )
        if k <= 4:
            require(
                abs(bareiss_determinant(matrix)) == expected_determinant,
                f"C3 tensor determinant k={k}",
            )

        if k <= 3:
            adjacency = transitive_c3_join(k)
            h = all_induced_h(adjacency)
            base_mask = sum(1 << (3 * gadget) for gadget in range(k))
            for left in range(size):
                for right in range(size):
                    selected = base_mask
                    for gadget in range(k):
                        if (left >> gadget) & 1:
                            selected |= 1 << (3 * gadget + 1)
                        if (right >> gadget) & 1:
                            selected |= 1 << (3 * gadget + 2)
                    require(
                        h[selected] == matrix[left][right],
                        f"literal C3 tensor face k={k}",
                    )
                    actual_face_checks += 1

        rows.append(
            {
                "k": k,
                "size": size,
                "rank_over_Q_via_mod_prime": size,
                "determinant_absolute": expected_determinant,
                "distinct_values": k + 1,
            }
        )

    return {"rows": rows, "actual_tournament_face_checks_k1_to_k3": actual_face_checks}


def local_internal_factor(
    adjacency: tuple[int, ...], left: int, right: int, inserted: tuple[int, ...]
) -> tuple[int, ...]:
    coefficients = []
    for mask in range(1 << len(inserted)):
        selected = tuple(
            inserted[index]
            for index in range(len(inserted))
            if (mask >> index) & 1
        )
        count = 0
        for ordering in permutations(selected):
            word = (left,) + ordering + (right,)
            count += all(
                arc(adjacency, word[index], word[index + 1])
                for index in range(len(word) - 1)
            )
        coefficients.append(count)
    return tuple(coefficients)


def probe_strict_degree_delay() -> dict[str, object]:
    adjacency = decode_raw(42, 4)
    factor = local_internal_factor(adjacency, 0, 1, (2, 3))
    defect = int(not arc(adjacency, 0, 1))
    live_degrees = [
        mask.bit_count() for mask, coefficient in enumerate(factor) if coefficient
    ]
    require(factor == (0, 0, 0, 1), "strict degree-delay local factor")
    require(defect == 1 and min(live_degrees) == 2,
            "strict degree delay exceeds defect")
    expected_arcs = (
        "1>0",
        "0>2",
        "3>0",
        "1>2",
        "3>1",
        "2>3",
    )
    observed_arcs = tuple(
        f"{left}>{right}" if arc(adjacency, left, right) else f"{right}>{left}"
        for left in range(4)
        for right in range(left + 1, 4)
    )
    require(observed_arcs == expected_arcs, "strict delay arc portrait")
    return {
        "n": 4,
        "code": 42,
        "base_word": [0, 1],
        "inserted": [2, 3],
        "defect": defect,
        "local_coefficients_by_mask": list(factor),
        "first_live_degree": min(live_degrees),
        "arcs": list(observed_arcs),
    }


def polynomial_add(left: list[int], right: list[int]) -> list[int]:
    out = [0] * max(len(left), len(right))
    for index, value in enumerate(left):
        out[index] += value
    for index, value in enumerate(right):
        out[index] += value
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def polynomial_multiply(left: list[int], right: list[int]) -> list[int]:
    out = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return out


def permutation_sign(word: tuple[int, ...]) -> int:
    inversions = sum(
        word[i] > word[j]
        for i in range(len(word))
        for j in range(i + 1, len(word))
    )
    return -1 if inversions & 1 else 1


def characteristic_polynomial(matrix: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    """Ascending coefficients of det(tI-M), by exact Leibniz expansion."""
    n = len(matrix)
    require(n > 0 and all(len(row) == n for row in matrix),
            "square characteristic matrix")
    result = [0]
    for word in permutations(range(n)):
        term = [permutation_sign(word)]
        for row, column in enumerate(word):
            entry = [-matrix[row][column], 1] if row == column else [-matrix[row][column]]
            term = polynomial_multiply(term, entry)
        result = polynomial_add(result, term)
    result += [0] * (n + 1 - len(result))
    return tuple(result)


def matrix_vector(
    matrix: tuple[tuple[int, ...], ...], vector: tuple[int, ...]
) -> tuple[int, ...]:
    return tuple(sum(a * b for a, b in zip(row, vector)) for row in matrix)


def face_table(
    adjacency: tuple[int, ...], base_mask: int, inserted: tuple[int, ...]
) -> tuple[int, ...]:
    h = all_induced_h(adjacency)
    return tuple(
        h[base_mask | local_mask_to_global(mask, inserted)]
        for mask in range(1 << len(inserted))
    )


def probe_algebraic_hostiles() -> dict[str, object]:
    c3_face = face_table(decode_raw(2, 3), 1 << 0, (1, 2))
    require(c3_face == (1, 1, 1, 3), "C3 face table")
    a, b, c, d = c3_face
    rayleigh = b * c - a * d
    require(rayleigh == -2, "C3 negative Rayleigh slack")

    c3_hessian = ((2, 1, 1), (1, 0, 3), (1, 3, 0))
    require(
        characteristic_polynomial(c3_hessian) == (12, -11, -2, 1),
        "C3 homogenized Hessian characteristic polynomial",
    )
    eigenpairs = (
        (-3, (0, -1, 1)),
        (1, (-2, 1, 1)),
        (4, (1, 1, 1)),
    )
    for eigenvalue, vector in eigenpairs:
        require(
            matrix_vector(c3_hessian, vector)
            == tuple(eigenvalue * entry for entry in vector),
            "C3 exact Hessian eigenpair",
        )
    require(sum(eigenvalue > 0 for eigenvalue, _ in eigenpairs) == 2,
            "C3 Hessian has two positive eigenvalues")

    slice_face = face_table(decode_raw(10, 5), 1 << 0, (1, 2, 3, 4))
    degree_two = tuple(
        tuple(
            0 if left == right else slice_face[(1 << left) | (1 << right)]
            for right in range(4)
        )
        for left in range(4)
    )
    expected_matrix = (
        (0, 3, 1, 3),
        (3, 0, 1, 1),
        (1, 1, 0, 3),
        (3, 1, 3, 0),
    )
    require(degree_two == expected_matrix, "degree-two slice Hessian")
    first_factor = [-13, -4, 1]
    second_factor = [-1, 4, 1]
    expected_charpoly = tuple(polynomial_multiply(first_factor, second_factor))
    require(
        characteristic_polynomial(degree_two) == expected_charpoly,
        "degree-two hostile characteristic factorization",
    )
    # t^2-4t-13 has one positive root; t^2+4t-1 also has one because each
    # has negative constant term.  Thus the symmetric Hessian has two positive
    # eigenvalues, which is the explicit quadratic non-Lorentzian criterion.
    require(first_factor[0] < 0 and second_factor[0] < 0,
            "two quadratic factors each straddle zero")

    log_concave_face = face_table(decode_raw(4, 4), (1 << 1) | (1 << 2), (0, 3))
    log_super_face = face_table(decode_raw(4, 4), (1 << 0) | (1 << 3), (1, 2))
    require(log_concave_face == (1, 1, 1, 5), "rank log-concavity hostile")
    require(2 * 2 < 1 * 5, "ordinary rank sequence not log-concave")
    require(log_super_face == (1, 3, 3, 5), "log-supermodular hostile")
    require(3 * 3 > 1 * 5, "log-supermodularity fails")

    return {
        "c3_face": list(c3_face),
        "c3_rayleigh_at_origin": rayleigh,
        "c3_homogenized_hessian_eigenvalues": [-3, 1, 4],
        "degree_two_slice": {
            "n": 5,
            "code": 10,
            "base": [0],
            "inserted": [1, 2, 3, 4],
            "hessian": [list(row) for row in degree_two],
            "characteristic_factors_ascending": [first_factor, second_factor],
            "positive_eigenvalue_count": 2,
        },
        "rank_log_concavity_face": list(log_concave_face),
        "log_supermodularity_face": list(log_super_face),
    }


def extend(adjacency: tuple[int, ...], signature: int) -> tuple[int, ...]:
    n = len(adjacency)
    require(0 <= signature < (1 << n), "ear signature range")
    child = list(adjacency) + [0]
    for old in range(n):
        if (signature >> old) & 1:
            child[n] |= 1 << old
        else:
            child[old] |= 1 << n
    out = tuple(child)
    require(tournament_well_formed(out), "ear extension tournament")
    return out


def boundary_data_literal(
    adjacency: tuple[int, ...]
) -> tuple[tuple[int, ...], tuple[int, ...], tuple[tuple[int, ...], ...]]:
    """Literal Start, End, Q with the exposed step allowed either direction."""
    n = len(adjacency)
    starts = [0] * n
    ends = [0] * n
    q = [[0] * n for _ in range(n)]
    for word in permutations(range(n)):
        legal_steps = tuple(
            arc(adjacency, word[index], word[index + 1])
            for index in range(n - 1)
        )
        if all(legal_steps):
            starts[word[0]] += 1
            ends[word[-1]] += 1
        for exposed in range(n - 1):
            if all(
                legal_steps[index]
                for index in range(n - 1)
                if index != exposed
            ):
                q[word[exposed]][word[exposed + 1]] += 1
    return tuple(starts), tuple(ends), tuple(tuple(row) for row in q)


def directed_cut_table(
    starts: tuple[int, ...],
    ends: tuple[int, ...],
    q: tuple[tuple[int, ...], ...],
) -> tuple[int, ...]:
    n = len(starts)
    require(len(ends) == n and len(q) == n, "cut boundary dimensions")
    values = []
    for signature in range(1 << n):
        total = sum(
            starts[vertex]
            for vertex in range(n)
            if (signature >> vertex) & 1
        )
        total += sum(
            ends[vertex]
            for vertex in range(n)
            if not ((signature >> vertex) & 1)
        )
        total += sum(
            q[left][right]
            for left in range(n)
            if not ((signature >> left) & 1)
            for right in range(n)
            if (signature >> right) & 1
        )
        values.append(total)
    return tuple(values)


def is_strong(adjacency: tuple[int, ...]) -> bool:
    require(tournament_well_formed(adjacency), "strongness tournament input")
    n = len(adjacency)
    full = (1 << n) - 1
    for start in range(n):
        seen = 1 << start
        frontier = seen
        while frontier:
            next_frontier = 0
            vertices = frontier
            while vertices:
                bit = vertices & -vertices
                vertex = bit.bit_length() - 1
                vertices -= bit
                next_frontier |= adjacency[vertex]
            next_frontier &= ~seen
            seen |= next_frontier
            frontier = next_frontier
        if seen != full:
            return False
    return True


def probe_ear_cut_curvature() -> dict[str, object]:
    tournament_count = 0
    cut_value_checks = 0
    second_difference_checks = 0
    third_difference_checks = 0
    strong_census: dict[int, dict[str, int]] = {}
    c3_second_values: tuple[int, ...] | None = None
    witness_record: dict[str, object] | None = None

    for n in range(1, 6):
        strong_count = 0
        interval_failures = 0
        first_failure_code: int | None = None
        code_count = 1 << (n * (n - 1) // 2)
        for code in range(code_count):
            adjacency = decode_raw(code, n)
            starts, ends, q = boundary_data_literal(adjacency)
            cut_values = directed_cut_table(starts, ends, q)
            direct_values = tuple(
                all_induced_h(extend(adjacency, signature))[-1]
                for signature in range(1 << n)
            )
            require(cut_values == direct_values,
                    f"directed-cut ear table n={n} code={code}")
            require(all(value > 0 and value & 1 for value in cut_values),
                    "ear responses positive odd")
            require(sum(starts) == sum(ends) == all_induced_h(adjacency)[-1],
                    "boundary Hamilton totals")

            tournament_count += 1
            cut_value_checks += len(cut_values)

            full = (1 << n) - 1
            for left in range(n):
                for right in range(left + 1, n):
                    pair = (1 << left) | (1 << right)
                    expected = -(q[left][right] + q[right][left])
                    for context in submasks(full ^ pair):
                        observed = mixed_difference(cut_values, pair, context)
                        require(observed == expected,
                                "cut second difference equals negative Q pair")
                        require(observed <= 0, "ear cut submodularity")
                        second_difference_checks += 1

            for triple_vertices in combinations(range(n), 3):
                triple = sum(1 << vertex for vertex in triple_vertices)
                for context in submasks(full ^ triple):
                    require(
                        mixed_difference(cut_values, triple, context) == 0,
                        "ear cut third difference vanishes",
                    )
                    third_difference_checks += 1

            if n == 3 and code == 2:
                values = []
                for left, right in combinations(range(3), 2):
                    pair = (1 << left) | (1 << right)
                    values.append(mixed_difference(cut_values, pair, 0))
                c3_second_values = tuple(values)

            if n >= 3 and is_strong(adjacency):
                strong_count += 1
                nonconstant = set(cut_values[1:-1])
                solid = set(range(min(nonconstant), max(nonconstant) + 1, 2))
                if nonconstant != solid:
                    interval_failures += 1
                    if first_failure_code is None:
                        first_failure_code = code

            if n == 5 and code == 8:
                require(is_strong(adjacency), "code-8 base strong")
                nonconstant = sorted(set(cut_values[1:-1]))
                interaction = tuple(
                    tuple(
                        0 if left == right else q[left][right] + q[right][left]
                        for right in range(n)
                    )
                    for left in range(n)
                )
                expected_interaction = (
                    (0, 6, 4, 4, 8),
                    (6, 0, 10, 10, 4),
                    (4, 10, 0, 10, 4),
                    (4, 10, 10, 0, 6),
                    (8, 4, 4, 6, 0),
                )
                require(nonconstant == [15, 17, 19, 23, 25, 27, 29, 33, 37, 41],
                        "code-8 nonconstant ear image")
                require(interaction == expected_interaction,
                        "code-8 symmetric Q interaction")
                require(bareiss_determinant(interaction) == -24_320,
                        "code-8 interaction determinant")
                require(rank_mod(interaction) == 5,
                        "code-8 interaction full rank")
                witness_record = {
                    "n": 5,
                    "code": 8,
                    "H": all_induced_h(adjacency)[-1],
                    "nonconstant_ear_image": nonconstant,
                    "first_internal_missing_odd": 21,
                    "interaction": [list(row) for row in interaction],
                    "interaction_determinant": -24_320,
                    "interaction_rank": 5,
                }

        if n >= 3:
            strong_census[n] = {
                "strong_labeled": strong_count,
                "non_solid_individual_ear_images": interval_failures,
                "first_failure_code": first_failure_code
                if first_failure_code is not None
                else -1,
            }

    require(tournament_count == 1_099, "ear tournament universe")
    require(cut_value_checks == 33_866, "ear cut-value universe")
    require(second_difference_checks == 83_506,
            "ear second-difference universe")
    require(third_difference_checks == 41_480,
            "ear third-difference universe")
    require(strong_census[3] == {
        "strong_labeled": 2,
        "non_solid_individual_ear_images": 0,
        "first_failure_code": -1,
    }, "strong order-three census")
    require(strong_census[4] == {
        "strong_labeled": 24,
        "non_solid_individual_ear_images": 0,
        "first_failure_code": -1,
    }, "strong order-four census")
    require(strong_census[5] == {
        "strong_labeled": 544,
        "non_solid_individual_ear_images": 544,
        "first_failure_code": 8,
    }, "strong order-five census")
    require(c3_second_values is not None and all(value < 0 for value in c3_second_values),
            "C3 ear-cut negative curvature")
    require(witness_record is not None, "code-8 witness recorded")

    return {
        "tournaments_n1_to_n5": tournament_count,
        "cut_value_checks": cut_value_checks,
        "second_difference_checks": second_difference_checks,
        "third_difference_checks": third_difference_checks,
        "strong_census": {str(key): value for key, value in strong_census.items()},
        "c3_presence_mixed_difference": 2,
        "c3_ear_cut_second_differences": list(c3_second_values),
        "code8_witness": witness_record,
    }


def main() -> None:
    ocf = probe_ocf_and_mobius()
    tensors = probe_c3_tensors()
    delay = probe_strict_degree_delay()
    hostiles = probe_algebraic_hostiles()
    ear_cut = probe_ear_cut_curvature()

    ledger = {
        "schema": "tournament-gap-polynomial-ocf-curvature-candidate-v1",
        "arithmetic": "exact-integer-no-floating-point",
        "ocf_mobius": ocf,
        "c3_transfer_tensors": tensors,
        "strict_degree_delay": delay,
        "algebraic_hostiles": hostiles,
        "presence_vs_ear_cut": ear_cut,
    }
    semantic = sha256(
        json.dumps(ledger, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, "semantic digest")

    print("THM4099 OCF-CURVATURE CANDIDATE PRIMARY EXACT REFEREE")
    print(
        "ocf_universe=",
        ocf["tournaments_n1_to_n6"],
        "tournaments_n1_to_n6",
        ocf["induced_masks"],
        "induced_masks",
        ocf["directed_cycle_instances"],
        "directed_cycle_instances",
    )
    print(
        "mobius_universe=",
        ocf["mobius_base_splits_n1_to_n5"],
        "base_splits_n1_to_n5",
        ocf["mobius_atom_checks"],
        "atoms",
        ocf["mixed_difference_checks"],
        "mixed_differences",
        "minimum=",
        ocf["minimum_mixed_difference"],
    )
    print(
        "parity_checks=",
        ocf["parity_checks"],
        "all_face_coefficients_positive_odd_nonempty_atoms_even=PASS",
    )
    print("c3_tensor_rows=", tensors["rows"])
    print(
        "c3_tensor_actual_tournament_face_checks=",
        tensors["actual_tournament_face_checks_k1_to_k3"],
    )
    print("strict_degree_delay=", delay)
    print("algebraic_hostiles=", hostiles)
    print(
        "ear_cut_universe=",
        ear_cut["tournaments_n1_to_n5"],
        "tournaments_n1_to_n5",
        ear_cut["cut_value_checks"],
        "cut_values",
        ear_cut["second_difference_checks"],
        "second_differences",
        ear_cut["third_difference_checks"],
        "third_differences",
    )
    print("presence_vs_ear_cut=", {
        "c3_presence_mixed_difference": ear_cut["c3_presence_mixed_difference"],
        "c3_ear_cut_second_differences": ear_cut["c3_ear_cut_second_differences"],
    })
    print("strong_ear_census=", ear_cut["strong_census"])
    print("strong_code8_hostile=", ear_cut["code8_witness"])
    print("semantic_sha256=", semantic)
    print("ALL_CHECKS_PASS")


if __name__ == "__main__":
    main()
