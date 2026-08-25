#!/usr/bin/env python3
"""Independent exact audit for THM-4099.

This verifier does not import the primary implementation.  Hamiltonian-path
counts are obtained by literal permutation enumeration.  The gap polynomial
is assembled independently as a generic subset-incidence convolution: local
blocks are enumerated for every inserted subset, then disjoint subset states
are propagated across the gaps of each base word.

All labelled tournaments and all nonempty prefix-base splittings are checked
through final order five.  Because the universe contains every labelled
tournament, relabelling any chosen base onto the prefix shows that these are
representatives of every labelled split orbit.  The full two-insertion slice
is additionally checked at final order six.  The audit also reconstructs the
two-variable matrix, the sequential L(1)+L'(1)+F partition, the degree/defect
cutoff, and the sharp singleton-base hostile.  Every gate uses ``require`` so
``-O`` removes no verification.
"""

from __future__ import annotations

from functools import lru_cache
from hashlib import sha256
from itertools import permutations
import json
from math import factorial


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError("FAILED: " + label)


def pair_index(n: int, first: int, second: int) -> int:
    require(0 <= first < second < n, "pair-index domain")
    return first * (2 * n - first - 1) // 2 + second - first - 1


def decode(mask: int, n: int) -> tuple[int, ...]:
    """LSB-first lexicographic upper-triangle code; bit one means i -> j."""
    out = [0] * n
    for first in range(n):
        for second in range(first + 1, n):
            if (mask >> pair_index(n, first, second)) & 1:
                out[first] |= 1 << second
            else:
                out[second] |= 1 << first
    return tuple(out)


def has_arc(out: tuple[int, ...], source: int, target: int) -> bool:
    return bool((out[source] >> target) & 1)


def directed(out: tuple[int, ...], word: tuple[int, ...]) -> bool:
    return all(has_arc(out, word[index], word[index + 1])
               for index in range(len(word) - 1))


@lru_cache(maxsize=None)
def subset_vertices(n: int, subset: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(n) if (subset >> vertex) & 1)


@lru_cache(maxsize=None)
def ordered_words(vertices: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(permutations(vertices))


def literal_paths(
    out: tuple[int, ...], subset: int, cache: dict[int, tuple[tuple[int, ...], ...]]
) -> tuple[tuple[int, ...], ...]:
    if subset not in cache:
        vertices = subset_vertices(len(out), subset)
        if not vertices:
            cache[subset] = ((),)
        else:
            cache[subset] = tuple(
                word for word in ordered_words(vertices) if directed(out, word)
            )
    return cache[subset]


def literal_h(
    out: tuple[int, ...], subset: int, cache: dict[int, tuple[tuple[int, ...], ...]]
) -> int:
    return len(literal_paths(out, subset, cache))


def actual_vertices(inserted: tuple[int, ...], local_subset: int) -> tuple[int, ...]:
    return tuple(
        inserted[index]
        for index in range(len(inserted))
        if (local_subset >> index) & 1
    )


def local_gap_vector(
    out: tuple[int, ...],
    inserted: tuple[int, ...],
    left: int | None,
    right: int | None,
) -> tuple[int, ...]:
    """Literal local block counts, indexed by subsets of inserted vertices."""
    require(left is not None or right is not None, "nonempty base gap")
    values = []
    for local_subset in range(1 << len(inserted)):
        vertices = actual_vertices(inserted, local_subset)
        count = 0
        for block in ordered_words(vertices):
            word = (() if left is None else (left,)) + block
            word += () if right is None else (right,)
            count += int(directed(out, word))
        values.append(count)
    return tuple(values)


def squarefree_step(state: tuple[int, ...], local: tuple[int, ...]) -> tuple[int, ...]:
    """Generic disjoint-subset convolution, not a hard-coded two-variable ring."""
    require(len(state) == len(local), "subset-state dimensions")
    size = len(state)
    answer = [0] * size
    for used, old_count in enumerate(state):
        if old_count == 0:
            continue
        for added, local_count in enumerate(local):
            if local_count and not (used & added):
                answer[used | added] += old_count * local_count
    return tuple(answer)


def incidence_matrix(local: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    """Matrix of disjoint union by a local squarefree coefficient vector."""
    size = len(local)
    matrix = []
    for target in range(size):
        row = []
        for source in range(size):
            if source & ~target:
                row.append(0)
            else:
                added = target ^ source
                row.append(local[added] if not (source & added) else 0)
        matrix.append(tuple(row))
    return tuple(matrix)


def apply_matrix(
    matrix: tuple[tuple[int, ...], ...], state: tuple[int, ...]
) -> tuple[int, ...]:
    return tuple(
        sum(entry * coefficient for entry, coefficient in zip(row, state))
        for row in matrix
    )


def explicit_two_matrix(local: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    require(len(local) == 4, "two-insertion local vector")
    a, b, c, d = local
    return (
        (a, 0, 0, 0),
        (b, a, 0, 0),
        (c, 0, a, 0),
        (d, c, b, a),
    )


def explicit_two_gap_formula(
    out: tuple[int, ...],
    left: int | None,
    right: int | None,
    first_inserted: int,
    second_inserted: int,
) -> tuple[int, int, int, int]:
    """Equations (11)--(13), organized as four explicit directed chains."""
    def chain(block: tuple[int, ...]) -> int:
        word = (() if left is None else (left,)) + block
        word += () if right is None else (right,)
        return int(directed(out, word))

    return (
        chain(()),
        chain((first_inserted,)),
        chain((second_inserted,)),
        chain((first_inserted, second_inserted))
        + chain((second_inserted, first_inserted)),
    )


def word_gaps(word: tuple[int, ...]) -> tuple[tuple[int | None, int | None], ...]:
    require(bool(word), "nonempty base word")
    return (
        ((None, word[0]),)
        + tuple((word[index], word[index + 1]) for index in range(len(word) - 1))
        + ((word[-1], None),)
    )


def base_word_polynomial(
    out: tuple[int, ...],
    word: tuple[int, ...],
    inserted: tuple[int, ...],
    gap_cache: dict[tuple[int | None, int | None], tuple[int, ...]],
) -> tuple[tuple[int, ...], int, int]:
    """Expand one base word by generic subset states across all of its gaps."""
    state = (1,) + (0,) * ((1 << len(inserted)) - 1)
    matrix_checks = 0
    for boundary in word_gaps(word):
        if boundary not in gap_cache:
            gap_cache[boundary] = local_gap_vector(out, inserted, *boundary)
        local = gap_cache[boundary]
        generic = squarefree_step(state, local)
        if len(inserted) == 2:
            require(
                local == explicit_two_gap_formula(
                    out, *boundary, inserted[0], inserted[1]
                ),
                "two-insertion local coefficient formulas",
            )
            matrix = incidence_matrix(local)
            require(matrix == explicit_two_matrix(local), "two-insertion matrix formula")
            require(apply_matrix(matrix, state) == generic,
                    "two-insertion matrix contraction")
            matrix_checks += 1
        state = generic

    defect = sum(
        not has_arc(out, word[index], word[index + 1])
        for index in range(len(word) - 1)
    )
    for local_subset, coefficient in enumerate(state):
        require(coefficient == 0 or local_subset.bit_count() >= defect,
                "degree-defect cutoff")
    if defect > len(inserted):
        require(not any(state), "above-degree base word vanishes")
    return state, defect, matrix_checks


def z_coefficients(
    out: tuple[int, ...], base: tuple[int, ...], inserted: tuple[int, ...]
) -> tuple[tuple[int, ...], dict[str, object]]:
    """Sum the squarefree gap products over all permutations of the base."""
    total = [0] * (1 << len(inserted))
    matrix_checks = 0
    cutoff_checks = 0
    dead_words = 0
    defect_histogram: dict[int, int] = {}
    live_layer_examples: dict[int, dict[str, object]] = {}
    gap_cache: dict[tuple[int | None, int | None], tuple[int, ...]] = {}
    for word in ordered_words(base):
        coefficients, defect, checks = base_word_polynomial(
            out, word, inserted, gap_cache
        )
        matrix_checks += checks
        cutoff_checks += len(coefficients)
        defect_histogram[defect] = defect_histogram.get(defect, 0) + 1
        if defect > len(inserted):
            dead_words += 1
        for local_subset, coefficient in enumerate(coefficients):
            total[local_subset] += coefficient
            degree = local_subset.bit_count()
            if coefficient and degree == defect and defect not in live_layer_examples:
                live_layer_examples[defect] = {
                    "word": list(word),
                    "subset": local_subset,
                    "coefficient": coefficient,
                }
    return tuple(total), {
        "matrix_checks": matrix_checks,
        "cutoff_checks": cutoff_checks,
        "dead_words": dead_words,
        "defects": tuple(sorted(defect_histogram.items())),
        "live_layers": live_layer_examples,
    }


def expected_face_counts(
    out: tuple[int, ...],
    base: tuple[int, ...],
    inserted: tuple[int, ...],
    path_cache: dict[int, tuple[tuple[int, ...], ...]],
) -> tuple[int, ...]:
    base_mask = sum(1 << vertex for vertex in base)
    values = []
    for local_subset in range(1 << len(inserted)):
        subset = base_mask
        for index, vertex in enumerate(inserted):
            if (local_subset >> index) & 1:
                subset |= 1 << vertex
        values.append(literal_h(out, subset, path_cache))
    return tuple(values)


def insertion_tau(out: tuple[int, ...], path: tuple[int, ...], vertex: int) -> int:
    signature = tuple(int(has_arc(out, vertex, old)) for old in path)
    return sum(
        signature[index] == 1 and signature[index + 1] == 0
        for index in range(len(signature) - 1)
    )


def insert_vertex(
    path: tuple[int, ...], vertex: int, position: int
) -> tuple[int, ...]:
    return path[:position] + (vertex,) + path[position:]


def sequential_audit(
    out: tuple[int, ...],
    path_cache: dict[int, tuple[tuple[int, ...], ...]],
) -> dict[str, object]:
    """Literal partition proving H(A+v)=L(1)+L'(1)+F."""
    n = len(out)
    require(n >= 3, "two-insertion sequential order")
    vertex = n - 1
    a_mask = (1 << vertex) - 1
    full_mask = (1 << n) - 1
    paths_a = literal_paths(out, a_mask, path_cache)
    paths_final = literal_paths(out, full_mask, path_cache)

    l_histogram: dict[int, int] = {}
    generated_nonorphans: set[tuple[int, ...]] = set()
    for path in paths_a:
        tau = insertion_tau(out, path, vertex)
        l_histogram[tau] = l_histogram.get(tau, 0) + 1
        legal = []
        for position in range(len(path) + 1):
            candidate = insert_vertex(path, vertex, position)
            if directed(out, candidate):
                legal.append(candidate)
                generated_nonorphans.add(candidate)
        require(len(legal) == 1 + tau, "one plus descent insertion law")

    literal_nonorphans = set()
    literal_orphans = set()
    for final_path in paths_final:
        reduced = tuple(old for old in final_path if old != vertex)
        if directed(out, reduced):
            literal_nonorphans.add(final_path)
        else:
            literal_orphans.add(final_path)
    require(generated_nonorphans == literal_nonorphans,
            "literal nonorphan partition")

    repaired_orphans = set()
    a_vertices = subset_vertices(n, a_mask)
    for word in ordered_words(a_vertices):
        bad = tuple(
            index
            for index in range(len(word) - 1)
            if not has_arc(out, word[index], word[index + 1])
        )
        if len(bad) != 1:
            continue
        cut = bad[0]
        if (has_arc(out, word[cut], vertex)
                and has_arc(out, vertex, word[cut + 1])):
            candidate = insert_vertex(word, vertex, cut + 1)
            require(directed(out, candidate), "failed-gap repair is directed")
            repaired_orphans.add(candidate)
    require(repaired_orphans == literal_orphans, "literal orphan repair partition")

    l_at_one = sum(l_histogram.values())
    derivative = sum(exponent * count for exponent, count in l_histogram.items())
    failed = len(literal_orphans)
    require(l_at_one + derivative + failed == len(paths_final),
            "L(1)+L'(1)+F")
    require(l_at_one + derivative == len(literal_nonorphans),
            "derivative counts nonorphans")
    return {
        "L_histogram": tuple(sorted(l_histogram.items())),
        "L1": l_at_one,
        "Lprime1": derivative,
        "F": failed,
        "H": len(paths_final),
        "nonorphans": len(literal_nonorphans),
    }


def merge_histogram(target: dict[int, int], source: tuple[tuple[int, int], ...]) -> None:
    for key, value in source:
        target[key] = target.get(key, 0) + value


def audit_universe(
    n: int, base_sizes: tuple[int, ...], do_sequential: bool
) -> tuple[dict[str, object], dict[str, object] | None, dict[int, dict[str, object]]]:
    tournament_count = 1 << (n * (n - 1) // 2)
    coefficient_checks = 0
    matrix_checks = 0
    cutoff_checks = 0
    split_cases = 0
    base_words = 0
    dead_words = 0
    defect_histogram: dict[int, int] = {}
    layer_examples: dict[int, dict[str, object]] = {}
    sequential_totals = {
        "tournaments": tournament_count,
        "L1": 0,
        "Lprime1": 0,
        "F": 0,
        "H": 0,
        "nonorphans": 0,
        "L_histogram": {},
    } if do_sequential else None

    for mask in range(tournament_count):
        out = decode(mask, n)
        path_cache: dict[int, tuple[tuple[int, ...], ...]] = {}
        for base_size in base_sizes:
            require(1 <= base_size < n, "nonempty proper base")
            base = tuple(range(base_size))
            inserted = tuple(range(base_size, n))
            actual, stats = z_coefficients(out, base, inserted)
            expected = expected_face_counts(out, base, inserted, path_cache)
            require(actual == expected, "[X_S]Z equals induced Hamilton count")
            coefficient_checks += len(actual)
            matrix_checks += int(stats["matrix_checks"])
            cutoff_checks += int(stats["cutoff_checks"])
            dead_words += int(stats["dead_words"])
            merge_histogram(defect_histogram, stats["defects"])
            for defect, example in stats["live_layers"].items():
                if n - base_size == 2 and defect not in layer_examples:
                    layer_examples[defect] = {
                        "n": n,
                        "mask": mask,
                        "base_size": base_size,
                        **example,
                    }
            split_cases += 1
            base_words += factorial(base_size)

        if do_sequential:
            require(sequential_totals is not None, "sequential accumulator")
            row = sequential_audit(out, path_cache)
            for key in ("L1", "Lprime1", "F", "H", "nonorphans"):
                sequential_totals[key] += int(row[key])
            merge_histogram(sequential_totals["L_histogram"], row["L_histogram"])

    row = {
        "final_n": n,
        "tournaments": tournament_count,
        "base_sizes": base_sizes,
        "split_cases": split_cases,
        "coefficient_checks": coefficient_checks,
        "base_words": base_words,
        "matrix_checks": matrix_checks,
        "cutoff_checks": cutoff_checks,
        "dead_words": dead_words,
        "defect_histogram": tuple(sorted(defect_histogram.items())),
        "failures": 0,
    }
    if sequential_totals is not None:
        sequential_totals["L_histogram"] = tuple(
            sorted(sequential_totals["L_histogram"].items())
        )
        sequential_totals["failures"] = 0
    return row, sequential_totals, layer_examples


def first_insertion_profile(out: tuple[int, ...]) -> tuple[int, tuple[int, ...], int, int]:
    require(len(out) == 3, "singleton hostile order")
    base_path = (0,)
    paths_a = tuple(
        word for word in ordered_words((0, 1)) if directed(out, word)
    )
    fibre = tuple(word for word in paths_a if tuple(x for x in word if x != 1) == base_path)
    orphans = tuple(word for word in paths_a if word not in fibre)
    return 1, (len(fibre),), len(orphans), len(paths_a)


def hostile_audit() -> dict[str, object]:
    # Pair bits are (0,1),(0,2),(1,2).  Masks 3 and 2 have the same arcs
    # 0->2 and 2->1; only the already-inserted/base arc is reversed.
    rows = []
    for name, mask in (("transitive", 3), ("C3", 2)):
        out = decode(mask, 3)
        cache: dict[int, tuple[tuple[int, ...], ...]] = {}
        z, _ = z_coefficients(out, (0,), (1, 2))
        expected = expected_face_counts(out, (0,), (1, 2), cache)
        require(z == expected, "hostile gap identity")
        sequential = sequential_audit(out, cache)
        rows.append({
            "name": name,
            "mask": mask,
            "Z": z,
            "first_profile": first_insertion_profile(out),
            "L_histogram": sequential["L_histogram"],
            "Lprime1": sequential["Lprime1"],
            "F": sequential["F"],
            "H": sequential["H"],
        })
    require(rows[0]["Z"][:3] == rows[1]["Z"][:3] == (1, 1, 1),
            "singleton hostile agrees on proper squarefree faces")
    require(rows[0]["first_profile"] == rows[1]["first_profile"]
            == (1, (1,), 0, 1), "singleton hostile first-insertion profile")
    require(rows[0]["Z"][3] == 1 and rows[1]["Z"][3] == 3,
            "singleton hostile differs only at mixed coefficient")
    for mask in (3, 2):
        out = decode(mask, 3)
        require(has_arc(out, 0, 2) and has_arc(out, 2, 1),
                "singleton hostile common future arcs")

    empty_base_rows = []
    for mask in range(2):
        out = decode(mask, 2)
        cache: dict[int, tuple[tuple[int, ...], ...]] = {}
        face = tuple(literal_h(out, subset, cache) for subset in range(4))
        empty_base_rows.append((mask, face))
    require(empty_base_rows == [(0, (1, 1, 1, 1)), (1, (1, 1, 1, 1))],
            "empty-base two-vertex boundary")
    return {"singleton_rows": rows, "empty_base_rows": empty_base_rows}


def main() -> None:
    print("THM-4099 INDEPENDENT SQUAREFREE GAP AUDIT")
    general_rows = []
    sequential_rows = []
    layer_examples: dict[int, dict[str, object]] = {}

    # Every labelled tournament and every nonempty prefix-base splitting.
    for n in range(2, 6):
        row, sequential, examples = audit_universe(
            n, tuple(range(1, n)), n >= 3
        )
        general_rows.append(row)
        if sequential is not None:
            sequential_rows.append(sequential)
        for defect, example in examples.items():
            layer_examples.setdefault(defect, example)
        print(
            "general=", row["final_n"], row["tournaments"], row["base_sizes"],
            row["split_cases"], row["coefficient_checks"], row["matrix_checks"],
            row["cutoff_checks"], row["dead_words"], row["failures"]
        )

    # Stronger two-insertion boundary: all 2^15 labelled order-six tournaments.
    row6, sequential6, examples6 = audit_universe(6, (4,), True)
    general_rows.append(row6)
    require(sequential6 is not None, "order-six sequential row")
    sequential_rows.append(sequential6)
    for defect, example in examples6.items():
        layer_examples.setdefault(defect, example)
    print(
        "two_insert_n6=", row6["tournaments"], row6["split_cases"],
        row6["coefficient_checks"], row6["matrix_checks"],
        row6["cutoff_checks"], row6["dead_words"], row6["failures"]
    )

    require([row["tournaments"] for row in general_rows]
            == [2, 8, 64, 1024, 32768], "labelled universe sizes")
    require(set(layer_examples) >= {0, 1, 2}, "live defect layers through two")
    require(sum(row["dead_words"] for row in general_rows) > 0,
            "strict above-degree cutoff is exercised")

    print("SEQUENTIAL L(1)+L'(1)+F")
    for row in sequential_rows:
        print(
            "sequential=", row["tournaments"], row["L1"], row["Lprime1"],
            row["F"], row["H"], row["nonorphans"],
            row["L_histogram"], row["failures"]
        )

    print("DEGREE-DEFECT LIVE LAYERS")
    for defect in sorted(layer_examples):
        print("layer=", defect, layer_examples[defect])

    hostile = hostile_audit()
    print("MINIMAL SINGLETON HOSTILE")
    for row in hostile["singleton_rows"]:
        print(
            "hostile=", row["name"], row["mask"], row["Z"],
            row["first_profile"], row["L_histogram"],
            row["Lprime1"], row["F"], row["H"]
        )
    print("empty_base_boundary=", hostile["empty_base_rows"])

    ledger = {
        "theorem": "THM-4099",
        "method": "literal permutations plus generic subset-incidence gap convolution",
        "general_rows": general_rows,
        "sequential_rows": sequential_rows,
        "layer_examples": layer_examples,
        "hostile": hostile,
    }
    semantic = sha256(
        json.dumps(ledger, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    print("semantic_sha256=", semantic)
    print("ALL_CHECKS_PASS")


if __name__ == "__main__":
    main()
