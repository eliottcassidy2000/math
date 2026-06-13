#!/usr/bin/env python3
"""S634: apply perspective-flip carriers to LRC, unit distance, and counts.

This scout turns the S629/S630 complement-perspective machinery into practical
checks for three live problems:

1. LRC source/sink target compression.  THM-381 makes loneliness a marked
   source event; converse pairs that with the marked sink obstruction.
2. Rooted tournament class counting modulo converse.  Burnside reduces the
   quotient to all rooted perspectives plus the SC perspectives fixed by the
   S629 flip.
3. Unit-distance n=21 exact cores.  A Hamiltonian unit spine is an actual graph
   property of the five published 57-edge cores, and the remaining bulk is the
   centered-hex C_hex(3)=37 carrier isolated by S630.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import sc_perspective_flip_cyclotomic_s629 as s629
import unit_distance_n22_tournament_lrc_s614 as s614


def score(mask: int, n: int, v: int) -> int:
    return sum(s629.edge(mask, n, v, u) for u in range(n) if u != v)


def orbit_score(mask: int, n: int, orbit: tuple[int, ...]) -> int:
    return score(mask, n, orbit[0])


def rooted_converse_lrc_table(max_n: int = 6) -> list[dict[str, object]]:
    rows = []
    for n in range(1, max_n + 1):
        reps = s629.classes(n)
        total_rooted = 0
        converse_fixed_rooted = 0
        sc_classes = 0
        source_targets = 0
        sink_targets = 0
        source_sink_merged = 0
        perspective_size_hist: Counter[int] = Counter()

        for mask in reps:
            auts = s629.automorphisms(mask, n)
            orbits = s629.vertex_orbits_from_auts(auts, n)
            total_rooted += len(orbits)
            perspective_size_hist.update([len(orbits)])

            has_source = False
            has_sink = False
            for orbit in orbits:
                out = orbit_score(mask, n, orbit)
                if out == n - 1:
                    source_targets += 1
                    has_source = True
                if out == 0:
                    sink_targets += 1
                    has_sink = True

            if n == 1:
                source_sink_merged += int(has_source or has_sink)
            else:
                # Converse pairs the marked source class with the marked sink
                # class.  A nontrivial tournament cannot have the same root be
                # both source and sink.
                source_sink_merged += int(has_source)

            antis = s629.anti_automorphisms(mask, n)
            if antis:
                sc_classes += 1
                flip = s629.perspective_flip(orbits, antis[0])
                if not s629.is_involution(flip):
                    raise AssertionError(("non-involutive perspective flip", n, mask, flip))
                converse_fixed_rooted += sum(1 for i, image in enumerate(flip) if i == image)

        if (total_rooted + converse_fixed_rooted) % 2:
            raise AssertionError(("Burnside parity failure", n))
        rooted_mod_converse = (total_rooted + converse_fixed_rooted) // 2

        expected_source = 1 if n == 1 else s629.KNOWN_UNROOTED[n - 1]
        if source_targets != expected_source or sink_targets != expected_source:
            raise AssertionError(
                (
                    "source/sink target count",
                    n,
                    source_targets,
                    sink_targets,
                    expected_source,
                )
            )

        rows.append(
            {
                "n": n,
                "classes": len(reps),
                "rooted_perspectives": total_rooted,
                "sc_classes": sc_classes,
                "converse_fixed_rooted": converse_fixed_rooted,
                "rooted_mod_converse": rooted_mod_converse,
                "source_targets": source_targets,
                "sink_targets": sink_targets,
                "source_sink_merged": source_sink_merged,
                "perspective_count_hist": dict(sorted(perspective_size_hist.items())),
            }
        )
    return rows


def graph_bits(adj: list[list[int]]) -> list[int]:
    bits = []
    for i, row in enumerate(adj):
        mask = 0
        for j, value in enumerate(row):
            if i != j and value:
                mask |= 1 << j
        bits.append(mask)
    return bits


def low_bit_index(mask: int) -> int:
    return (mask & -mask).bit_length() - 1


def hamiltonian_path(adj: list[list[int]]) -> tuple[int, ...] | None:
    n = len(adj)
    adj_bits = graph_bits(adj)
    full = (1 << n) - 1
    ends = [0] * (1 << n)
    for v in range(n):
        ends[1 << v] = 1 << v

    for mask in range(1 << n):
        endpoint_bits = ends[mask]
        while endpoint_bits:
            v_bit = endpoint_bits & -endpoint_bits
            v = v_bit.bit_length() - 1
            endpoint_bits ^= v_bit
            candidates = adj_bits[v] & ~mask
            while candidates:
                u_bit = candidates & -candidates
                candidates ^= u_bit
                ends[mask | u_bit] |= u_bit

    if not ends[full]:
        return None

    end = low_bit_index(ends[full])
    mask = full
    path = [end]
    while mask != (1 << end):
        prev_mask = mask ^ (1 << end)
        prev_candidates = ends[prev_mask] & adj_bits[end]
        if not prev_candidates:
            raise AssertionError(("lost predecessor", mask, end))
        prev = low_bit_index(prev_candidates)
        path.append(prev)
        mask = prev_mask
        end = prev
    return tuple(reversed(path))


def centered_hex(k: int) -> int:
    return 3 * k * (k + 1) + 1


def n21_core_spine_table() -> list[dict[str, object]]:
    rows = []
    for idx, code in enumerate(s614.N21_GRAPH6, start=1):
        adj = s614.graph6_decode(code)
        n = len(adj)
        edges = s614.edge_count(adj)
        degrees = s614.degrees(adj)
        path = hamiltonian_path(adj)
        bulk = None if path is None else edges - (n - 1)
        shell = "-"
        if bulk is not None:
            for k in range(1, 8):
                if bulk == centered_hex(k):
                    shell = f"C_hex({k})"
        rows.append(
            {
                "core": idx,
                "n": n,
                "edges": edges,
                "traceable": path is not None,
                "spine_edges": None if path is None else n - 1,
                "bulk_edges": bulk,
                "bulk_shell": shell,
                "degree_hist": s614.degree_hist(adj),
                "min_degree": min(degrees),
                "max_degree": max(degrees),
                "delete_edge_deck": dict(sorted(Counter(edges - d for d in degrees).items())),
                "triangles": s614.count_triangles(adj),
                "four_cycles": s614.count_4cycles(adj),
                "path_endpoints": None if path is None else (path[0], path[-1]),
                "path_endpoint_degrees": None
                if path is None
                else (degrees[path[0]], degrees[path[-1]]),
                "path": None if path is None else path,
            }
        )
    return rows


def majority_route_tournament() -> dict[str, object]:
    routes = [
        {
            "name": "lrc_source_sink_pairing",
            "proof_power": 5,
            "compression": 5,
            "computation_ready": 5,
            "cross_problem": 4,
            "risk": 1,
        },
        {
            "name": "rooted_converse_quotient_counts",
            "proof_power": 4,
            "compression": 4,
            "computation_ready": 5,
            "cross_problem": 3,
            "risk": 1,
        },
        {
            "name": "n21_unit_spine_bulk_audit",
            "proof_power": 3,
            "compression": 3,
            "computation_ready": 5,
            "cross_problem": 5,
            "risk": 2,
        },
        {
            "name": "n22_frontier_extension_constraints",
            "proof_power": 5,
            "compression": 3,
            "computation_ready": 3,
            "cross_problem": 5,
            "risk": 4,
        },
        {
            "name": "raw_scalar_phi3_matching",
            "proof_power": 1,
            "compression": 2,
            "computation_ready": 4,
            "cross_problem": 4,
            "risk": 5,
        },
    ]
    criteria = ["proof_power", "compression", "computation_ready", "cross_problem"]
    n = len(routes)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        wins_i = sum(routes[i][c] > routes[j][c] for c in criteria)
        wins_j = sum(routes[j][c] > routes[i][c] for c in criteria)
        if wins_i > wins_j or (wins_i == wins_j and routes[i]["risk"] < routes[j]["risk"]):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    scores = Counter(sum(row) for row in adj)
    cycles = 0
    for a, b, c in combinations(range(n), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            cycles += 1
    return {
        "routes": [route["name"] for route in routes],
        "score_hist": dict(sorted(scores.items())),
        "directed_3_cycles": cycles,
        "leaders": [
            routes[i]["name"]
            for i in range(n)
            if sum(adj[i]) == max(sum(row) for row in adj)
        ],
    }


def format_dict(d: dict[object, object]) -> str:
    return "{" + ", ".join(f"{k}: {v}" for k, v in d.items()) + "}"


def main() -> None:
    print("S634 APPLIED PERSPECTIVE-CARRIER SCOUT")
    print("=" * 78)
    print()

    print("LRC / rooted tournament counting")
    print("-" * 78)
    print(
        "n classes rooted fixed_by_converse rooted_mod_converse source sink "
        "source_sink_merged perspective_hist"
    )
    for row in rooted_converse_lrc_table(6):
        print(
            f"{row['n']:>1} {row['classes']:>7} {row['rooted_perspectives']:>6} "
            f"{row['converse_fixed_rooted']:>17} {row['rooted_mod_converse']:>19} "
            f"{row['source_targets']:>6} {row['sink_targets']:>4} "
            f"{row['source_sink_merged']:>18} "
            f"{format_dict(row['perspective_count_hist'])}"
        )
    print()
    print("Interpretation:")
    print("  * THM-381 source targets are exactly U(n-1); sink targets match by")
    print("    converse.  Merging the source/sink pair therefore keeps U(n-1)")
    print("    marked targets instead of carrying two separate obstruction lists.")
    print("  * Rooted perspectives modulo converse are counted by Burnside as")
    print("    (all rooted perspectives + SC perspective-flip fixed roots)/2.")
    print()

    print("Unit-distance n=21 exact-core Hamiltonian spine audit")
    print("-" * 78)
    print(
        "core edges traceable spine bulk shell min/max-degree degree_hist "
        "delete_edge_deck triangles C4 endpoint_degrees"
    )
    for row in n21_core_spine_table():
        print(
            f"{row['core']:>4} {row['edges']:>5} {str(row['traceable']):>9} "
            f"{row['spine_edges']:>5} {row['bulk_edges']:>4} {row['bulk_shell']:>8} "
            f"{row['min_degree']}/{row['max_degree']} "
            f"{format_dict(row['degree_hist'])} "
            f"{format_dict(row['delete_edge_deck'])} "
            f"{row['triangles']:>9} {row['four_cycles']:>5} "
            f"{row['path_endpoint_degrees']}"
        )
    print()
    print("Interpretation:")
    print("  * Every traceable 57-edge n=21 core has an actual unit Hamiltonian")
    print("    spine, so the S630 decomposition is graph-real: 57 = 20 + 37.")
    print("  * The 37-edge complement to a spine is C_hex(3), a centered-hex")
    print("    carrier rather than an H=21 tournament scalar.")
    print()

    print("Application route tournament")
    print("-" * 78)
    route = majority_route_tournament()
    print(f"routes={route['routes']}")
    print(f"score_hist={route['score_hist']}")
    print(f"directed_3_cycles={route['directed_3_cycles']}")
    print(f"leaders={route['leaders']}")
    print()
    print("Challenged assumptions:")
    print("  * LRC tournament vertices need not be runners; the useful quotient here")
    print("    is marked rooted perspectives/source obligations.")
    print("  * Unit-distance tournament vertices need not be points; the useful")
    print("    witness is a path-plus-bulk carrier split inside the unit graph.")
    print("  * Counting tournament classes should not count a rooted perspective and")
    print("    its converse obstruction separately when the perspective flip already")
    print("    gives the involution.")


if __name__ == "__main__":
    main()
