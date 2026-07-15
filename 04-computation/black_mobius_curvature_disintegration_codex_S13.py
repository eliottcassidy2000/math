#!/usr/bin/env python3
"""Exact black-edge curvature disintegration for THM-811.

The audit joins THM-785's longitudinal complement flux, THM-801's boundary
Mobius curvature q, and the antisymmetric Smith current epsilon.  It checks a
single four-variable line polynomial, node-level curvature classifiers, black
reflection-orbit codecs, and source-normalized boundary flow by q.

Tournament Analysis uses node information carriers as vertices.  Its pairwise
observable is the number of merged-node pairs separated; retention and
retention per log-cell are the switches; increasing curvature content is the
tie Hamiltonian path.

Assumption challenge: a black edge is first a reflection orbit of literal
complement lines.  Signed endpoints, mirror positions, and multiplicity are
retained.  Curvature preserves flow information but deliberately does not
pretend to preserve literal edge identity or an LRC loneliness predicate.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

from mobius_cech_metagraph_codec_codex_S12 import (
    b2_signature,
    b3_signature,
    boundary_c3_curvature,
    compact_partition,
    counter_json,
    endpoint_epsilon,
    tiling_c3,
)
from three_sorted_metagraph_recursion_codex_S9 import (
    N7_ATLAS,
    SMALL_ATLAS,
    apex_zero_endpoint,
    blue_fraction_exponent,
    complement,
    is_grid_symmetric,
    line_index,
    load_atlases,
    m_tiles,
    reflect,
    reflection_orbits,
    tile_index,
    tile_schema,
)
from tournament_tiling_metagraph_address_codex_S4 import carrier_tournament


OUTPUT = Path("05-knowledge/results/black_mobius_curvature_disintegration_codex_S13.out")
JSON_OUTPUT = Path("05-knowledge/results/black_mobius_curvature_disintegration_codex_S13.json")


def predicted_master(n: int) -> Counter[tuple[int, int, int, int]]:
    """Coefficients indexed by (DeltaC3,q0,q1,epsilon)."""
    dp = Counter({(n - 2, 0, 0, 0): 1})
    for epsilon_step in (1, -1):
        nxt = Counter()
        for (flux, q0, q1, epsilon), count in dp.items():
            nxt[(flux, q0, q1, epsilon)] += count
            nxt[(flux - 1, q0, q1 + 1, epsilon + epsilon_step)] += count
        dp = nxt
    for _ in range(n - 4):
        nxt = Counter()
        for (flux, q0, q1, epsilon), count in dp.items():
            nxt[(flux, q0, q1, epsilon)] += count
            nxt[(flux - 1, q0, q1, epsilon + 1)] += count
            nxt[(flux - 1, q0, q1, epsilon - 1)] += count
            nxt[(flux - 2, q0 + 1, q1 + 1, epsilon)] += count
        dp = nxt
    factor = 1 << (m_tiles(n) - 2 * n + 5)
    return Counter({key: factor * value for key, value in dp.items()})


def covariance(rows: list[tuple[int, int]]) -> Fraction:
    total = len(rows)
    sx = sum(x for x, _ in rows)
    sy = sum(y for _, y in rows)
    sxy = sum(x * y for x, y in rows)
    return Fraction(sxy, total) - Fraction(sx, total) * Fraction(sy, total)


def primitive_polynomial(row: Counter[tuple[int, int]]) -> tuple:
    divisor = 0
    from math import gcd
    for value in row.values():
        divisor = gcd(divisor, value)
    return tuple(sorted((key, value // divisor) for key, value in row.items()))


def strict_rate(count: int, mass: int) -> Fraction:
    return Fraction(count, mass) if mass else Fraction()


def size_audit(n: int, node_by_mask: dict[int, list[int]]) -> tuple[dict, dict[str, dict[int, object]]]:
    node_map = node_by_mask[n]
    full = (1 << m_tiles(n)) - 1
    line_count = 1 << (m_tiles(n) - 1)
    _, sigma = tile_schema(n)

    masks_by_node: dict[int, list[int]] = defaultdict(list)
    for mask, node in enumerate(node_map):
        masks_by_node[node].append(mask)
    node_colour = {}
    for node, masks in masks_by_node.items():
        blue = sum(is_grid_symmetric(mask, sigma) for mask in masks)
        node_colour[node] = "pure_blue" if blue == len(masks) else "pure_black" if blue == 0 else "mixed"

    master = Counter()
    black_moment_rows = []
    all_moment_rows = []
    black_lines = []
    b3_current_failures = 0
    flux_current_failures = 0

    node_joint: dict[int, Counter[tuple[int, int]]] = defaultdict(Counter)
    node_q: dict[int, Counter[int]] = defaultdict(Counter)
    node_eps: dict[int, Counter[int]] = defaultdict(Counter)
    node_c3 = {}
    for node in masks_by_node:
        node_joint[node]
        node_q[node]
        node_eps[node]

    source_mass: Counter[tuple[str, int]] = Counter()
    for mask, node in enumerate(node_map):
        if not is_grid_symmetric(mask, sigma):
            q = boundary_c3_curvature(mask, n)
            eps = abs(endpoint_epsilon(mask, n))
            node_joint[node][(q, eps)] += 1
            node_q[node][q] += 1
            node_eps[node][eps] += 1
            source_mass[(node_colour[node], q)] += 1
        c3 = tiling_c3(mask, n)
        if node in node_c3:
            assert node_c3[node] == c3
        node_c3[node] = c3

    boundary_count: Counter[tuple[str, int]] = Counter()
    raw_boundary = Counter()
    loop_codec = {}
    loop_codec_collisions = Counter()

    for line in range(line_count):
        mask = apex_zero_endpoint(line, n)
        other = complement(mask, n)
        q0 = boundary_c3_curvature(mask, n)
        q1 = boundary_c3_curvature(other, n)
        epsilon = endpoint_epsilon(mask, n)
        flux = tiling_c3(other, n) - tiling_c3(mask, n)
        master[(flux, q0, q1, epsilon)] += 1
        all_moment_rows.append((q0, epsilon * epsilon))

        b = sum((mask >> tile_index(n)[(x, 1)]) & 1 for x in range(3, n))
        t = sum((mask >> tile_index(n)[(n, y)]) & 1 for y in range(2, n - 1))
        s3 = b3_signature(mask, n)
        b3_current_failures += b != s3[0] + s3[3]
        b3_current_failures += t != s3[2] + s3[5]
        b3_current_failures += epsilon != b - t
        flux_current_failures += flux != n - 2 - b - t

        blue = is_grid_symmetric(mask, sigma)
        if blue:
            continue
        black_moment_rows.append((q0, epsilon * epsilon))
        u, v = node_map[mask], node_map[other]
        pair = tuple(sorted((u, v)))
        s2 = b2_signature(mask, n)
        s3 = b3_signature(mask, n)
        black_lines.append(
            {
                "line": line,
                "mask": mask,
                "other": other,
                "pair": pair,
                "q0": q0,
                "q1": q1,
                "epsilon": epsilon,
                "s2": s2,
                "s3": s3,
                "reflection_line": line_index(reflect(mask, n), n),
            }
        )

        if u == v:
            key = (pair, s2, s3)
            loop_codec_collisions[key] += 1

        c0, c1 = tiling_c3(mask, n), tiling_c3(other, n)
        if c0 == c1:
            continue
        source_mask, target_mask = (mask, other) if c0 < c1 else (other, mask)
        source_node, target_node = node_map[source_mask], node_map[target_mask]
        source_cat, target_cat = node_colour[source_node], node_colour[target_node]
        if {source_cat, target_cat} != {"mixed", "pure_black"}:
            continue
        direction = "outward" if source_cat == "mixed" else "reverse"
        source_q = boundary_c3_curvature(source_mask, n)
        boundary_count[(direction, source_q)] += 1
        raw_boundary[direction] += 1

    assert master == predicted_master(n)
    assert b3_current_failures == flux_current_failures == 0

    r, p = n - 4, n - 3
    all_cov = covariance(all_moment_rows)
    assert all_cov == Fraction(-r, 8)
    a = Fraction(1 << (reflection_orbits(n) - 1), line_count)
    delta = Fraction(0) if n % 2 == 0 else Fraction(1, 4)
    predicted_black_cov = (
        -Fraction(r, 8) * (1 - a) + a * Fraction(p, 2) * delta
    ) / (1 - a) ** 2
    black_cov = covariance(black_moment_rows)
    assert black_cov == predicted_black_cov

    q_poly = {node: tuple(sorted(row.items())) for node, row in node_q.items()}
    eps_poly = {node: tuple(sorted(row.items())) for node, row in node_eps.items()}
    joint_poly = {node: tuple(sorted(row.items())) for node, row in node_joint.items()}
    joint_orbit_poly = {
        node: tuple(sorted((key, value // 2) for key, value in row.items()))
        for node, row in node_joint.items()
    }
    assert all(value % 2 == 0 for row in node_joint.values() for value in row.values())
    node_partitions = {
        "q": compact_partition(q_poly),
        "q_plus_C3": compact_partition({node: (node_c3[node], q_poly[node]) for node in q_poly}),
        "absolute_epsilon": compact_partition(eps_poly),
        "q_absolute_epsilon": compact_partition(joint_poly),
        "C3_q_absolute_epsilon": compact_partition(
            {node: (node_c3[node], joint_orbit_poly[node]) for node in joint_orbit_poly}
        ),
        "C3_primitive_joint_plus_black_mass": compact_partition(
            {
                node: (
                    node_c3[node],
                    sum(node_joint[node].values()),
                    primitive_polynomial(node_joint[node]),
                )
                for node in node_joint
            }
        ),
    }
    assert node_partitions["C3_q_absolute_epsilon"]["cells"] == len(masks_by_node)
    assert node_partitions["C3_primitive_joint_plus_black_mass"]["cells"] == len(masks_by_node)

    by_line = {record["line"]: record for record in black_lines}
    orbit_reps = sorted({min(line, record["reflection_line"]) for line, record in by_line.items()})
    orbit_carriers = {
        "node_q_pair": {},
        "node_q_pair_absolute_epsilon": {},
        "node_orbit_B2_B3": {},
    }
    orbit_signature_failures = 0
    for rep in orbit_reps:
        record = by_line[rep]
        mate = by_line[record["reflection_line"]]
        signatures = tuple(sorted(((record["s2"], record["s3"]), (mate["s2"], mate["s3"]))))
        orbit_carriers["node_q_pair"][rep] = (record["pair"], record["q0"], record["q1"])
        orbit_carriers["node_q_pair_absolute_epsilon"][rep] = (
            record["pair"], record["q0"], record["q1"], abs(record["epsilon"])
        )
        orbit_carriers["node_orbit_B2_B3"][rep] = (record["pair"], signatures)
        orbit_signature_failures += record["pair"] != mate["pair"]
        orbit_signature_failures += (record["q0"], record["q1"]) != (mate["q0"], mate["q1"])
    orbit_stats = {name: compact_partition(values) for name, values in orbit_carriers.items()}
    assert orbit_signature_failures == 0
    assert orbit_stats["node_orbit_B2_B3"]["cells"] == len(orbit_reps)

    literal_join_counts = Counter(
        (record["pair"], record["s2"], record["s3"]) for record in black_lines
    )
    residual_cells = [key for key, value in literal_join_counts.items() if value > 1]
    reflection_pair_failures = 0
    for key in residual_cells:
        members = [record for record in black_lines if (record["pair"], record["s2"], record["s3"]) == key]
        reflection_pair_failures += len(members) != 2
        if len(members) == 2:
            reflection_pair_failures += members[0]["reflection_line"] != members[1]["line"]
    assert reflection_pair_failures == 0
    assert max(loop_codec_collisions.values(), default=1) == 1

    q_rates = {}
    outward_dominance_failures = 0
    all_q = sorted({q for _, q in source_mass} | {q for _, q in boundary_count})
    for q in all_q:
        outward = boundary_count[("outward", q)]
        reverse = boundary_count[("reverse", q)]
        mixed_mass = source_mass[("mixed", q)]
        black_mass = source_mass[("pure_black", q)]
        out_rate = strict_rate(outward, mixed_mass)
        rev_rate = strict_rate(reverse, black_mass)
        if n >= 6 and (outward or reverse):
            outward_dominance_failures += not (out_rate > rev_rate)
        q_rates[str(q)] = {
            "outward_count": outward,
            "reverse_count": reverse,
            "mixed_source_mass": mixed_mass,
            "pure_black_source_mass": black_mass,
            "outward_rate": str(out_rate),
            "reverse_rate": str(rev_rate),
        }
    if n >= 6:
        assert outward_dominance_failures == 0

    carriers = {
        "C3": node_c3,
        "black_q_polynomial": q_poly,
        "black_absolute_epsilon_polynomial": eps_poly,
        "black_joint_curvature": joint_orbit_poly,
        "C3_black_joint_curvature": {
            node: (node_c3[node], joint_orbit_poly[node]) for node in joint_orbit_poly
        },
    }
    return {
        "n": n,
        "lines": line_count,
        "black_lines": len(black_lines),
        "black_reflection_orbits": len(orbit_reps),
        "master_polynomial": (
            "2^(M-2n+5) u^(n-2)(1+w*v/u)(1+w/(u*v))"
            "(1+v/u+1/(u*v)+z*w/u^2)^(n-4)"
        ),
        "master_polynomial_failures": 0,
        "B3_determines_epsilon_lambda_failures": b3_current_failures,
        "flux_equals_minus_lambda_failures": flux_current_failures,
        "cov_q0_epsilon_squared_all": str(all_cov),
        "cov_q0_epsilon_squared_black": str(black_cov),
        "predicted_black_covariance": str(predicted_black_cov),
        "node_partitions": node_partitions,
        "black_orbit_carriers": orbit_stats,
        "literal_B2_B3_residual_collision_cells": len(residual_cells),
        "literal_B2_B3_residual_collision_excess": sum(
            value - 1 for value in literal_join_counts.values() if value > 1
        ),
        "residual_cells_are_reflection_pairs_failures": reflection_pair_failures,
        "black_loop_B2_B3_collision_excess": sum(value - 1 for value in loop_codec_collisions.values()),
        "black_boundary_raw_counts": dict(raw_boundary),
        "black_boundary_source_q_rates": q_rates,
        "source_q_outward_dominance_failures": outward_dominance_failures,
    }, carriers


def audit(small_atlas: Path, n7_atlas: Path) -> dict:
    _sizes, node_by_mask = load_atlases(small_atlas, n7_atlas)
    result = {"schema_version": 1, "theorem": "THM-811", "sizes": []}
    carriers_n7 = None
    for n in range(4, 8):
        size, carriers = size_audit(n, node_by_mask)
        result["sizes"].append(size)
        if n == 7:
            carriers_n7 = carriers
    assert carriers_n7 is not None
    stats = {name: compact_partition(values) for name, values in carriers_n7.items()}
    retention = carrier_tournament(stats, "retention")
    economy = carrier_tournament(stats, "economy")
    flips = sum(
        retention["adjacency"][i][j] != economy["adjacency"][i][j]
        for i in range(len(stats)) for j in range(i + 1, len(stats))
    )
    result["tournament_analysis_n7"] = {
        "vertices": list(stats),
        "pairwise_observable": "number of unordered merged-node pairs separated by the carrier",
        "switches": ["retention", "retention_per_log2_cells"],
        "tie_hamiltonian_path": list(stats),
        "carrier_stats": stats,
        "retention": retention,
        "economy": economy,
        "edge_flips_between_gauges": flips,
    }
    return result


def render(result: dict) -> str:
    lines = [
        "THM-811 BLACK MOBIUS/SMITH CURVATURE DISINTEGRATION",
        "=" * 82,
        "master: 2^(M-2n+5) u^(n-2)(1+w v/u)(1+w/(uv))",
        "        * (1+v/u+1/(uv)+zw/u^2)^(n-4)",
        "B3 gives leg totals (lambda,epsilon); q gives positional cross-leg overlap.",
        "",
    ]
    for size in result["sizes"]:
        p = size["node_partitions"]
        o = size["black_orbit_carriers"]
        lines.extend(
            [
                f"n={size['n']}: lines={size['lines']} black={size['black_lines']} "
                f"black reflection orbits={size['black_reflection_orbits']}",
                f"  master/B3-current/flux failures={size['master_polynomial_failures']}/"
                f"{size['B3_determines_epsilon_lambda_failures']}/"
                f"{size['flux_equals_minus_lambda_failures']}",
                f"  Cov(q0,epsilon^2) all/black="
                f"{size['cov_q0_epsilon_squared_all']}/"
                f"{size['cov_q0_epsilon_squared_black']}",
                f"  node cells q/q+C3/|eps|/joint/joint+C3="
                f"{p['q']['cells']}/{p['q_plus_C3']['cells']}/"
                f"{p['absolute_epsilon']['cells']}/{p['q_absolute_epsilon']['cells']}/"
                f"{p['C3_q_absolute_epsilon']['cells']}",
                f"  black orbit cells node+qpair/+|eps|/orbit(B2,B3)="
                f"{o['node_q_pair']['cells']}/{o['node_q_pair_absolute_epsilon']['cells']}/"
                f"{o['node_orbit_B2_B3']['cells']}",
                f"  literal B2+B3 residual cells/excess="
                f"{size['literal_B2_B3_residual_collision_cells']}/"
                f"{size['literal_B2_B3_residual_collision_excess']}; "
                f"reflection failures={size['residual_cells_are_reflection_pairs_failures']}",
                f"  loop B2+B3 excess={size['black_loop_B2_B3_collision_excess']}; "
                f"source-q outward-dominance failures={size['source_q_outward_dominance_failures']}",
            ]
        )
        if size["n"] >= 6:
            lines.append("  source-q rates:")
            for q, row in size["black_boundary_source_q_rates"].items():
                if row["outward_count"] or row["reverse_count"]:
                    lines.append(
                        f"    q={q}: outward={row['outward_count']}/{row['mixed_source_mass']}="
                        f"{row['outward_rate']} reverse={row['reverse_count']}/"
                        f"{row['pure_black_source_mass']}={row['reverse_rate']}"
                    )
        lines.append("")
    ta = result["tournament_analysis_n7"]
    lines.extend(
        [
            "TOURNAMENT ANALYSIS (node information carriers)",
            f"  vertices={tuple(ta['vertices'])}",
            f"  observable={ta['pairwise_observable']}",
            f"  switches={tuple(ta['switches'])}; edge flips={ta['edge_flips_between_gauges']}",
            f"  retention scores={ta['retention']['score_hist']} "
            f"C3={ta['retention']['directed_3cycles']} SCC={ta['retention']['scc_sizes']} "
            f"Hpaths={ta['retention']['hamiltonian_paths']}",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--small-atlas", type=Path, default=SMALL_ATLAS)
    parser.add_argument("--n7-atlas", type=Path, default=N7_ATLAS)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--json", type=Path, default=JSON_OUTPUT)
    args = parser.parse_args()
    result = audit(args.small_atlas, args.n7_atlas)
    text = render(result)
    print(text, end="")
    args.output.write_text(text)
    args.json.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
