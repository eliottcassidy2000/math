#!/usr/bin/env python3
"""Exact direct-realizer audit for the path-relative kappa operation.

The all-tile flip is indexed by the marked Hamiltonian path P.  Relabelling
therefore obeys

    sigma kappa_P = kappa_(sigma P) sigma,

not commutation with a fixed kappa_P.  If sigma T = kappa_P T, then

    sigma^2 T = kappa_(sigma P) kappa_P T
              = Flip_(E(P) symmetric_difference E(sigma P)) T.

This script imports THM-849's exact paired equitable-refinement engine and
adds witness extraction.  It enumerates every black direct realizer through
n=8, verifies the displayed identity arc by arc, and records permutation
cycle type and two-path holonomy size.  In particular it tests, rather than
assumes, the square-in-Aut hypothesis used conditionally in THM-852.

Tournament Analysis boundary: the pairwise observable is the orientation of
an unordered vertex pair in T, the relabelling sigma is the gauge transport,
and P is the declared directed Hamiltonian path.  The new obstruction does
not live on those vertices: its vertices are the two observer paths P and
sigma P, and its carrier is their undirected edge symmetric difference.
Forcing that two-path datum back to a tournament quotient would discard the
quantity being measured.  THM-849 already reports the score, C3, SCC, and H
fingerprints of the same black Klein orbits; this verifier reports the new
edge-flip fingerprint exactly.
"""

from __future__ import annotations

import hashlib
from collections import Counter
from pathlib import Path
from time import perf_counter

import n8_black_self_line_obstruction_codex_S15 as s15


HERE = Path(__file__).resolve().parent
DEPENDENCY = HERE / "n8_black_self_line_obstruction_codex_S15.py"


EXPECTED = {
    5: {
        "black": 8,
        "cycles": Counter({(3, 1, 1): 4, (5,): 4}),
        "holonomy": Counter({6: 8}),
        "joint": Counter({(6, (3, 1, 1)): 4, (6, (5,)): 4}),
    },
    6: {
        "black": 12,
        "cycles": Counter({(4, 2): 8, (5, 1): 4}),
        "holonomy": Counter({6: 8, 8: 4}),
        "joint": Counter({(6, (4, 2)): 8, (8, (5, 1)): 4}),
    },
    7: {
        "black": 88,
        "cycles": Counter({(5, 2): 48, (6, 1): 32, (4, 3): 8}),
        "holonomy": Counter({8: 32, 10: 40, 12: 16}),
        "joint": Counter(
            {
                (8, (5, 2)): 16,
                (8, (6, 1)): 16,
                (10, (5, 2)): 32,
                (10, (6, 1)): 8,
                (12, (4, 3)): 8,
                (12, (6, 1)): 8,
            }
        ),
    },
    8: {
        "black": 404,
        "cycles": Counter(
            {(5, 2, 1): 188, (4, 2, 2): 128, (8,): 56, (4, 3, 1): 32}
        ),
        "holonomy": Counter({6: 136, 8: 8, 10: 84, 12: 128, 14: 48}),
        "joint": Counter(
            {
                (6, (4, 2, 2)): 128,
                (6, (8,)): 8,
                (8, (8,)): 8,
                (10, (5, 2, 1)): 64,
                (10, (8,)): 20,
                (12, (5, 2, 1)): 124,
                (12, (8,)): 4,
                (14, (4, 3, 1)): 32,
                (14, (8,)): 16,
            }
        ),
    },
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_isomorphisms(
    left: tuple[int, ...], right: tuple[int, ...]
) -> tuple[tuple[int, ...], ...]:
    """Witness-returning companion to S15's exact isomorphism counter."""
    colors = s15.paired_equitable_colors(left, right)
    if colors is None:
        return ()
    left_colors, right_colors = colors
    n = len(left)
    candidates = {
        vertex: [
            image
            for image in range(n)
            if right_colors[image] == left_colors[vertex]
        ]
        for vertex in range(n)
    }
    image = [-1] * n
    witnesses: list[tuple[int, ...]] = []

    def search(used: int, assigned: int) -> None:
        if assigned == n:
            witnesses.append(tuple(image))
            return

        best_vertex = -1
        best_options: list[int] | None = None
        for vertex in range(n):
            if image[vertex] >= 0:
                continue
            options = []
            for target in candidates[vertex]:
                if used >> target & 1:
                    continue
                if all(
                    image[other] < 0
                    or ((left[vertex] >> other) & 1)
                    == ((right[target] >> image[other]) & 1)
                    for other in range(n)
                ):
                    options.append(target)
            if not options:
                return
            if best_options is None or len(options) < len(best_options):
                best_vertex, best_options = vertex, options

        assert best_options is not None
        for target in best_options:
            image[best_vertex] = target
            search(used | (1 << target), assigned + 1)
            image[best_vertex] = -1

    search(0, 0)
    return tuple(witnesses)


def relabel(rows: tuple[int, ...], image: tuple[int, ...]) -> tuple[int, ...]:
    """Relabel a tournament by the vertex map image[v]."""
    n = len(rows)
    result = [0] * n
    for vertex in range(n):
        for other in range(n):
            if rows[vertex] >> other & 1:
                result[image[vertex]] |= 1 << image[other]
    return tuple(result)


def compose(
    left: tuple[int, ...], right: tuple[int, ...]
) -> tuple[int, ...]:
    """Return left after right."""
    return tuple(left[right[vertex]] for vertex in range(len(left)))


def path_edges(path: tuple[int, ...]) -> frozenset[frozenset[int]]:
    return frozenset(
        frozenset((path[index], path[index + 1]))
        for index in range(len(path) - 1)
    )


def flip_edges(
    rows: tuple[int, ...], edges: frozenset[frozenset[int]]
) -> tuple[int, ...]:
    result = list(rows)
    for edge in edges:
        left, right = tuple(edge)
        assert ((result[left] >> right) & 1) ^ ((result[right] >> left) & 1)
        result[left] ^= 1 << right
        result[right] ^= 1 << left
    return tuple(result)


def cycle_type(image: tuple[int, ...]) -> tuple[int, ...]:
    seen: set[int] = set()
    lengths = []
    for vertex in range(len(image)):
        if vertex in seen:
            continue
        cursor = vertex
        length = 0
        while cursor not in seen:
            seen.add(cursor)
            length += 1
            cursor = image[cursor]
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def format_counter(counter: Counter[object]) -> str:
    return ", ".join(f"{key}:{counter[key]}" for key in sorted(counter, key=str))


def audit(n: int) -> dict[str, object]:
    cells, reflection = s15.staircase(n)
    dimension = len(cells)
    full = (1 << dimension) - 1
    base_path = tuple(range(n - 1, -1, -1))
    base_edges = path_edges(base_path)

    cycles: Counter[tuple[int, ...]] = Counter()
    holonomy: Counter[int] = Counter()
    overlap: Counter[int] = Counter()
    joint: Counter[tuple[int, tuple[int, ...]]] = Counter()
    black = score_survivors = square_in_aut = involutions = 0

    for word in range(1 << dimension):
        if s15.sigma_word(word, reflection) == word:
            continue
        tournament = s15.tournament_rows(n, cells, word)
        kappa_tournament = s15.tournament_rows(n, cells, word ^ full)
        if sorted(row.bit_count() for row in tournament) != sorted(
            row.bit_count() for row in kappa_tournament
        ):
            continue
        score_survivors += 1

        realizers = exact_isomorphisms(tournament, kappa_tournament)
        if not realizers:
            continue
        automorphisms = s15.isomorphism_count(tournament, tournament)
        assert len(realizers) == automorphisms
        assert automorphisms == 1
        black += 1

        realizer = realizers[0]
        assert relabel(tournament, realizer) == kappa_tournament
        transported_edges = path_edges(tuple(realizer[v] for v in base_path))
        defect = base_edges ^ transported_edges
        square = compose(realizer, realizer)

        # This is the path-relative equivariance identity, checked arc by arc.
        assert relabel(tournament, square) == flip_edges(tournament, defect)
        square_is_aut = relabel(tournament, square) == tournament
        assert square_is_aut == (not defect)
        square_in_aut += int(square_is_aut)
        involutions += int(square == tuple(range(n)))

        kind = cycle_type(realizer)
        cycles[kind] += 1
        holonomy[len(defect)] += 1
        overlap[len(base_edges & transported_edges)] += 1
        joint[(len(defect), kind)] += 1

    expected = EXPECTED[n]
    assert black == expected["black"]
    assert cycles == expected["cycles"]
    assert holonomy == expected["holonomy"]
    assert joint == expected["joint"]
    assert square_in_aut == 0
    assert involutions == 0
    assert sum(holonomy.values()) == black
    assert all(count % 4 == 0 for count in holonomy.values())

    return {
        "n": n,
        "dimension": dimension,
        "black": black,
        "selfK": black // 2,
        "klein_orbits": black // 4,
        "score_survivors": score_survivors,
        "cycles": cycles,
        "holonomy": holonomy,
        "overlap": overlap,
        "joint": joint,
        "square_in_aut": square_in_aut,
        "involutions": involutions,
    }


def main() -> None:
    started = perf_counter()
    rows = [audit(n) for n in range(5, 9)]

    print("THM-854: PATH-RELATIVE KAPPA HOLONOMY EXACT AUDIT")
    print("identity: sigma*kappa_P = kappa_(sigma P)*sigma")
    print("direct realizer: sigma*T=kappa_P*T")
    print("square law: sigma^2*T=Flip_(E(P) symmetric_difference E(sigma P))*T")
    print("method: S15 paired equitable refinement + complete witness backtracking")
    print()
    for row in rows:
        orbit_holonomy = Counter(
            {size: count // 4 for size, count in row["holonomy"].items()}
        )
        print(
            "n={n}: black={black}, selfK={selfK}, Klein_orbits={klein_orbits}, "
            "score_pass={score_survivors}, square_in_Aut={square_in_aut}, "
            "involutive_realizers={involutions}".format(**row)
        )
        print(f"  cycle types: {format_counter(row['cycles'])}")
        print(f"  |path symmetric difference|: {format_counter(row['holonomy'])}")
        print(f"  |path edge intersection|: {format_counter(row['overlap'])}")
        print(f"  Klein orbits by holonomy: {format_counter(orbit_holonomy)}")
        print(f"  holonomy x cycle type: {format_counter(row['joint'])}")

    print()
    print("conditional odd-coset hypothesis: failed on every direct black realizer")
    print("preserved datum: ordered pair of path edge sets (E(P), E(sigma P))")
    print("destroyed by class quotient: their symmetric-difference flip holonomy")
    print(f"source_sha256={sha256(Path(__file__).resolve())}")
    print(f"s15_dependency_sha256={sha256(DEPENDENCY)}")
    print(f"ALL ASSERTIONS PASSED ({perf_counter() - started:.3f}s)")


if __name__ == "__main__":
    main()
