#!/usr/bin/env python3
"""Exact THM-773 prime-seven sheet-monodromy / tournament-atlas bridge.

At c=7 and away from strict-danger endpoints, each of seven unramified
owners occupies one sheet.  This script audits three layers which must not be
conflated:

* the unlabelled exact-cover shape (one seven-sheet partition);
* the 7! owner-labelled assignments and their endpoint monodromy;
* ordinary tournament isomorphism classes and the HYP-6825 staircase fibre.

The pairwise observable for Tournament Analysis is owner precedence.  The
switch changes from a marked linear cut to circular three-step precedence.
Ties use owner order 0<...<6 (there are no pairwise ties in either gauge).
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from collections import Counter, defaultdict
from fractions import Fraction
from math import floor
from pathlib import Path

from merged_metagraph_lines_n7_klein_S161 import canon_tournament
from tournament_tiling_metagraph_address_codex_S4 import tile_schema


P = 7
DELTA = Fraction(1, 14)
ALL_PERMS = tuple(itertools.permutations(range(P)))


def mod1(x: Fraction) -> Fraction:
    return x - floor(x)


def circle_norm(x: Fraction) -> Fraction:
    y = mod1(x)
    return min(y, 1 - y)


def bad_sheets(w: int, x: Fraction) -> tuple[int, ...]:
    return tuple(k for k in range(P) if circle_norm(Fraction(w, P) * (x + k)) < DELTA)


def token(w: int, x: Fraction) -> int:
    """Unique bad sheet off the endpoint set."""
    q = floor(w * x + Fraction(1, 2))
    return (-pow(w, -1, P) * q) % P


def moments(position: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sum(pow(k, m, P) for k in position) % P for m in range(1, P))


FULL_GRID_MOMENTS = moments(tuple(range(P)))


def root_polynomial(position: tuple[int, ...]) -> tuple[int, ...]:
    """Ascending coefficients of product (X-k) over F_7."""
    coeff = [1]
    for k in position:
        nxt = [0] * (len(coeff) + 1)
        for i, value in enumerate(coeff):
            nxt[i] = (nxt[i] - k * value) % P
            nxt[i + 1] = (nxt[i + 1] + value) % P
        coeff = nxt
    return tuple(coeff)


FULL_GRID_POLYNOMIAL = (0, P - 1, 0, 0, 0, 0, 0, 1)


def linear_tournament(position: tuple[int, ...]) -> list[list[int]]:
    return [
        [int(a != b and position[a] < position[b]) for b in range(P)]
        for a in range(P)
    ]


def circular_tournament(position: tuple[int, ...]) -> list[list[int]]:
    return [
        [int(a != b and (position[b] - position[a]) % P in (1, 2, 3)) for b in range(P)]
        for a in range(P)
    ]


def canonical_code(adjacency: list[list[int]]) -> int:
    converse = [[adjacency[j][i] for j in range(P)] for i in range(P)]
    return min(
        int.from_bytes(canon_tournament(adjacency, P), "big"),
        int.from_bytes(canon_tournament(converse, P), "big"),
    )


def hamiltonian_paths(adjacency: list[list[int]]) -> tuple[tuple[int, ...], ...]:
    return tuple(
        path
        for path in ALL_PERMS
        if all(adjacency[path[i]][path[i + 1]] for i in range(P - 1))
    )


def least_hamiltonian_path(adjacency: list[list[int]]) -> tuple[int, ...]:
    for path in ALL_PERMS:
        if all(adjacency[path[i]][path[i + 1]] for i in range(P - 1)):
            return path
    raise AssertionError("every tournament has a Hamiltonian path")


def encode_fixed_path(adjacency: list[list[int]], path: tuple[int, ...]) -> int:
    """Relabel path source..sink to atlas labels 7..1 and return its mask."""
    tiles, _ = tile_schema(P)
    old_at_new = [0] * P
    for i, old in enumerate(path):
        old_at_new[P - 1 - i] = old
    relabelled = [
        [adjacency[old_at_new[i]][old_at_new[j]] for j in range(P)]
        for i in range(P)
    ]
    assert all(relabelled[i][i - 1] for i in range(P - 1, 0, -1))
    mask = 0
    for bit, (x, y) in enumerate(tiles):
        if relabelled[y - 1][x - 1]:
            mask |= 1 << bit
    return mask


def directed_triangles(adjacency: list[list[int]]) -> int:
    total = 0
    for a, b, c in itertools.combinations(range(P), 3):
        total += int(adjacency[a][b] and adjacency[b][c] and adjacency[c][a])
        total += int(adjacency[a][c] and adjacency[c][b] and adjacency[b][a])
    return total


def scc_sizes(adjacency: list[list[int]]) -> tuple[int, ...]:
    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for u in range(P):
            if adjacency[v][u] and u not in seen:
                dfs(u)
        order.append(v)

    for v in range(P):
        if v not in seen:
            dfs(v)
    seen.clear()
    sizes = []

    def rdfs(v: int) -> int:
        seen.add(v)
        size = 1
        for u in range(P):
            if adjacency[u][v] and u not in seen:
                size += rdfs(u)
        return size

    for v in reversed(order):
        if v not in seen:
            sizes.append(rdfs(v))
    return tuple(sorted(sizes, reverse=True))


def tournament_fingerprint(adjacency: list[list[int]]) -> dict:
    paths = hamiltonian_paths(adjacency)
    return {
        "score_histogram": dict(sorted(Counter(sum(row) for row in adjacency).items())),
        "directed_3cycles": directed_triangles(adjacency),
        "scc_sizes": scc_sizes(adjacency),
        "hamiltonian_paths": len(paths),
        "tie_hamiltonian_path": tuple(range(P)),
    }


def event_map(speeds: tuple[int, ...]) -> dict[Fraction, tuple[int, ...]]:
    owners: dict[Fraction, list[int]] = defaultdict(list)
    for a, w in enumerate(speeds):
        for n in range(w):
            owners[Fraction(2 * n + 1, 2 * w)].append(a)
    return {x: tuple(indices) for x, indices in owners.items()}


def chamber_movie(speeds: tuple[int, ...]) -> dict:
    events = event_map(speeds)
    points = sorted(events)
    chambers = []
    for i, left in enumerate(points):
        right_lift = points[i + 1] if i + 1 < len(points) else points[0] + 1
        # Keep the lifted coordinate on the final chamber.  Replacing x=1 by
        # x=0 reindexes every sheet and would erase the deck monodromy.
        midpoint = (left + right_lift) / 2
        position = tuple(token(w, midpoint) for w in speeds)
        for w, k in zip(speeds, position):
            assert bad_sheets(w, midpoint) == (k,)
        chambers.append(
            {
                "left": left,
                "right": mod1(right_lift),
                "right_lift": right_lift,
                "midpoint": midpoint,
                "position": position,
                "free": P - len(set(position)),
                "moments": moments(position),
                "polynomial": root_polynomial(position),
                "left_owners": events[left],
                "right_owners": events[mod1(right_lift)],
            }
        )

    transition_failures = []
    for i, chamber in enumerate(chambers):
        before = chambers[i - 1]["position"]
        # The first event is crossed at x=first+1 when closing the lifted
        # movie.  Its after-state is the first chamber with the global sheet
        # carry retained, not the gauge-reset state stored near x=first.
        after = (
            tuple(token(w, chamber["midpoint"] + 1) for w in speeds)
            if i == 0
            else chamber["position"]
        )
        owners = chamber["left_owners"]
        for a, w in enumerate(speeds):
            expected = (-pow(w, -1, P)) % P if a in owners else 0
            if (after[a] - before[a]) % P != expected:
                transition_failures.append((i, a))
        expected_first = -sum(pow(speeds[a], -1, P) for a in owners) % P
        actual_first = (sum(after) - sum(before)) % P
        if actual_first != expected_first:
            transition_failures.append((i, "first_moment"))

    exact_indices = [i for i, chamber in enumerate(chambers) if chamber["free"] == 0]
    return_holonomy = []
    for z, i in enumerate(exact_indices):
        j = exact_indices[(z + 1) % len(exact_indices)]
        event_indices = []
        cursor = (i + 1) % len(chambers)
        while True:
            event_indices.append(cursor)
            if cursor == j:
                break
            cursor = (cursor + 1) % len(chambers)
        inverse_sum = sum(
            pow(speeds[a], -1, P)
            for q in event_indices
            for a in chambers[q]["left_owners"]
        ) % P
        return_holonomy.append((i, j, inverse_sum, len(event_indices)))
        assert inverse_sum == 0

    exact_records = []
    for i in exact_indices:
        chamber = chambers[i]
        boundary = chamber["right"]
        boundary_lift = chamber["right_lift"]
        at_boundary = [bad_sheets(w, boundary_lift) for w in speeds]
        occupied = {k for row in at_boundary for k in row}
        exact_records.append(
            {
                "left": chamber["left"],
                "right": boundary,
                "width": chamber["right_lift"] - chamber["left"],
                "position": chamber["position"],
                "next_event_owners": tuple(speeds[a] for a in chamber["right_owners"]),
                "next_free_sheets": tuple(k for k in range(P) if k not in occupied),
            }
        )

    assert not transition_failures
    assert all(chambers[i]["moments"] == FULL_GRID_MOMENTS for i in exact_indices)
    assert all(chambers[i]["polynomial"] == FULL_GRID_POLYNOMIAL for i in exact_indices)
    return {
        "speeds": speeds,
        "event_count": len(points),
        "chamber_count": len(chambers),
        "exact_chambers": exact_records,
        "exact_assignments": len({chambers[i]["position"] for i in exact_indices}),
        "return_holonomy": return_holonomy,
        "transition_failures": transition_failures,
    }


def exhaustive_moment_audit() -> dict:
    powers = [[pow(k, m, P) for k in range(P)] for m in range(1, P)]
    exact = moment_matches = false_positive = false_negative = 0
    for position in itertools.product(range(P), repeat=P):
        is_exact = len(set(position)) == P
        signature = tuple(sum(table[k] for k in position) % P for table in powers)
        is_match = signature == FULL_GRID_MOMENTS
        exact += is_exact
        moment_matches += is_match
        false_positive += is_match and not is_exact
        false_negative += is_exact and not is_match
    assert exact == 5040 == moment_matches
    assert false_positive == false_negative == 0
    return {
        "configurations": P**P,
        "exact": exact,
        "moment_matches": moment_matches,
        "false_positive": false_positive,
        "false_negative": false_negative,
    }


def atlas_bridge(atlas_path: Path) -> dict:
    atlas = json.loads(atlas_path.read_text())
    code_to_node = {
        int(code, 16): node
        for node in atlas["nodes"]
        for code in node["canonical_orbit_codes"]
    }
    identity = tuple(range(P))
    linear = linear_tournament(identity)
    circular = circular_tournament(identity)
    linear_node = code_to_node[canonical_code(linear)]
    circular_node = code_to_node[canonical_code(circular)]
    assert linear_node["id"] == "n7-a000"
    assert circular_node["id"] == "n7-a267"

    records = []
    mask_counts: Counter[int] = Counter()
    labelled_tournaments = set()
    for position in ALL_PERMS:
        adjacency = circular_tournament(position)
        labelled_tournaments.add(tuple(tuple(row) for row in adjacency))
        path = least_hamiltonian_path(adjacency)
        mask = encode_fixed_path(adjacency, path)
        mask_counts[mask] += 1
        records.append(
            {
                "owner_to_sheet": position,
                "least_hamiltonian_path": path,
                "mask": mask,
            }
        )

    atlas_masks = tuple(circular_node["tiling_masks"])
    fibre_index = {mask: i for i, mask in enumerate(atlas_masks)}
    for record in records:
        record["fibre_index"] = fibre_index[record["mask"]]
    assert set(mask_counts) == set(atlas_masks)
    assert len(mask_counts) == circular_node["tiling_count"] == 25
    assert len(labelled_tournaments) == 720

    flips = sum(
        linear[i][j] != circular[i][j]
        for i in range(P)
        for j in range(i + 1, P)
    )
    assert flips == 6 == circular_node["local_depth"]
    assignment_digest = hashlib.sha256(
        json.dumps(records, separators=(",", ":")).encode()
    ).hexdigest()
    return {
        "assignments": len(records),
        "rotation_orbits": len(records) // 7,
        "dihedral_orbits": len(records) // 14,
        "labelled_circular_tournaments": len(labelled_tournaments),
        "linear_node": {k: linear_node[k] for k in ("id", "rank", "tiling_count")},
        "circular_node": {
            k: circular_node[k]
            for k in (
                "id",
                "rank",
                "tiling_count",
                "local_depth",
                "local_path_word",
                "blueblack_root_distance",
                "blueblack_root_word",
                "self_converse",
            )
        },
        "gauge_edge_flips": flips,
        "linear_fingerprint": tournament_fingerprint(linear),
        "circular_fingerprint": tournament_fingerprint(circular),
        "mask_counts": dict(sorted(mask_counts.items(), key=lambda item: fibre_index[item[0]])),
        "multiplicity_histogram": dict(sorted(Counter(mask_counts.values()).items())),
        "assignment_digest": assignment_digest,
        "records": records,
    }


def printable(value) -> str:
    if isinstance(value, Fraction):
        return str(value)
    return str(value)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--atlas",
        type=Path,
        default=Path("05-knowledge/results/tournament_tiling_metagraph_address_n7_codex_S4.json"),
    )
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()

    moment_audit = exhaustive_moment_audit()
    bridge = atlas_bridge(args.atlas)
    movies = [
        chamber_movie((1, 4, 5, 6, 8, 9, 10)),
        chamber_movie((12, 38, 72, 96, 151, 169, 188)),
    ]

    print("THM-773 PRIME-SEVEN SHEET MONODROMY / METAGRAPH FIBRE")
    print("=" * 76)
    print("MOMENT CERTIFICATE")
    print(f"  full-grid moments m=1..6: {FULL_GRID_MOMENTS}")
    print(f"  full-grid polynomial: {FULL_GRID_POLYNOMIAL}  [ascending]")
    print(f"  exhaustive token configurations: {moment_audit}")
    print()
    print("EXACT TOURNAMENT-ATLAS BRIDGE")
    print(f"  owner-labelled assignments: {bridge['assignments']}")
    print(f"  rotation / dihedral orbits: {bridge['rotation_orbits']} / {bridge['dihedral_orbits']}")
    print(f"  distinct labelled circular tournaments: {bridge['labelled_circular_tournaments']}")
    print(f"  cut gauge node: {bridge['linear_node']}")
    print(f"  circular gauge node: {bridge['circular_node']}")
    print(f"  gauge edge flips: {bridge['gauge_edge_flips']}")
    print(f"  cut fingerprint: {bridge['linear_fingerprint']}")
    print(f"  circular fingerprint: {bridge['circular_fingerprint']}")
    print(f"  25-mask multiplicity histogram: {bridge['multiplicity_histogram']}")
    print(f"  assignment map sha256: {bridge['assignment_digest']}")
    print("  atlas fibre mask -> assignment count:")
    for mask, count in bridge["mask_counts"].items():
        print(f"    {mask:5d} {mask:015b} weight={mask.bit_count():2d} count={count:3d}")
    print()
    print("CHAMBER MONODROMY")
    for movie in movies:
        print(
            f"  W={movie['speeds']}: events={movie['event_count']}, "
            f"exact_chambers={len(movie['exact_chambers'])}, "
            f"exact_assignments={movie['exact_assignments']}"
        )
        for chamber in movie["exact_chambers"]:
            print(
                "    "
                f"({printable(chamber['left'])},{printable(chamber['right'])}) "
                f"width={printable(chamber['width'])} pos={chamber['position']} "
                f"next_owner={chamber['next_event_owners']} "
                f"boundary_free={chamber['next_free_sheets']}"
            )
        print(f"    consecutive exact-return holonomy: {movie['return_holonomy']}")
        print(f"    transition failures: {movie['transition_failures']}")
    print()
    print("PRESERVATION AUDIT")
    print("  preserves at iso node: unlabelled transitive/heptagon shape only")
    print("  destroyed by iso node: owner->sheet assignment, inverse step, event phase, next free sheet")
    print("  challenged vertex assumption: owners, sheets, events, and fibre masks were all tested")
    print("  minimal theorem carrier here: labelled token assignment + inverse steps + endpoint schedule")
    print("  tie Hamiltonian path: owner order 0<...<6 (unused; both gauges are tournaments)")

    if args.json:
        payload = {
            "schema_version": 1,
            "theorem": "THM-773",
            "moment_audit": moment_audit,
            "full_grid_moments": FULL_GRID_MOMENTS,
            "full_grid_polynomial": FULL_GRID_POLYNOMIAL,
            "atlas_bridge": bridge,
            "chamber_movies": movies,
        }
        args.json.write_text(json.dumps(payload, separators=(",", ":"), default=str) + "\n")


if __name__ == "__main__":
    main()
