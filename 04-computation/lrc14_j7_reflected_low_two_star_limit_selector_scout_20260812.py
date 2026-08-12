#!/usr/bin/env python3
"""Exact finite scout for the homogeneous-limit two-star selector.

The reflected low-channel bank left by the upper-median two-star analysis has
649 bodies and 11,856 ordered primitive rays.  For each body/ray assignment,
this program computes the exact homogeneous overlap limit on all nine edges of
the canonical ``K_{2,4}`` plus its centre edge and selects a maximizing edge.
It then checks, exactly, that the selected limit and the selected physical
scale-one overlap both exceed the full six-clause singleton debt.

This is a finite exact compression, not an all-dilation theorem.  Its point is
to replace 7,694,544 unrelated scale-one winners by one canonical selector and
140,711 body/edge/channel configurations to which the fixed-channel
quasipolynomial machinery can subsequently be applied.
"""

from __future__ import annotations

import argparse
import hashlib
from collections import Counter
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "04-computation/lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.py"
LOW = ROOT / "04-computation/lrc14_j7_reflected_low_phase_clique_robust_body_closure_thm2941.py"
MEDIAN = ROOT / "04-computation/lrc14_j7_reflected_median_star_chord_scout_20260801.py"
EXPECTED_DEPENDENCY_HASHES = {
    BASE: "2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31",
    LOW: "416c36f16f7c821feb8d260882711d2717069147b8604a93ba60432785cf1d1c",
    MEDIAN: "b0eba7785cd73ecf76d2887f3d36316eae6980a76652058daa2587ecbe031276",
}
EXPECTED_COUNTS = (649, 11856, 7694544, 140711)
EXPECTED_WEAKEST = (
    (2, 4, 6, 9, 11, 12),
    5544,
    2970,
    0,
    1,
    (12, 6, 4, 18, 3, 9),
    (0, 5),
    F(1, 3168),
    F(20120188378425768804, 99637927188239340409565),
    F(
        35897170405386504838493,
        315652953332342230417501920,
    ),
)
EXPECTED_SEMANTIC_SHA256 = "08b86820b1562f50f317be81981b762dc80c0f438da07048f2b0df6000a2a458"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def load(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


for dependency, expected_hash in EXPECTED_DEPENDENCY_HASHES.items():
    actual_hash = hashlib.sha256(
        dependency.read_bytes().replace(b"\r\n", b"\n")
    ).hexdigest()
    require(
        actual_hash == expected_hash,
        ("dependency hash", dependency, actual_hash, expected_hash),
    )


R = load("limit_selector_base", BASE)
LP = load("limit_selector_low", LOW)
MS = load("limit_selector_median", MEDIAN)


def fracpart(value: F) -> F:
    return value - value.numerator // value.denominator


def bernoulli3(value: F) -> F:
    value = fracpart(value)
    return value**3 - F(3, 2) * value**2 + F(1, 2) * value


def interval_intersection_mass(first, second) -> F:
    i = j = 0
    total = F(0)
    while i < len(first) and j < len(second):
        total += max(
            F(0),
            min(first[i][1], second[j][1])
            - max(first[i][0], second[j][0]),
        )
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return total


def limiting_arcs(multiplier: int, phase: int, ruler: int):
    """Arcs ``||multiplier*x-phase/ruler||<1/14`` on the unit interval."""
    radius = F(1, 14 * multiplier)
    arcs = []
    for k in range(-1, multiplier + 2):
        centre = F(phase, ruler * multiplier) + F(k, multiplier)
        left, right = max(F(0), centre - radius), min(F(1), centre + radius)
        if left < right:
            arcs.append((left, right))
    return tuple(arcs)


def homogeneous_limit(
    ruler: int, cell: int, e: int, p: int, f: int, q: int
) -> F:
    """Exact common-dilation limit of one reflected pair overlap."""
    divisor = gcd(p, q)
    P, Q = p // divisor, q // divisor
    if P > Q:
        P, Q, e, f = Q, P, f, e
    phase_e, phase_f = (e * cell) % ruler, (f * cell) % ruler
    cross = Q * e - P * f
    if cross == 0:
        return interval_intersection_mass(
            limiting_arcs(P, phase_e, ruler),
            limiting_arcs(Q, phase_f, ruler),
        )
    determinant = Q * phase_e - P * phase_f
    a, b = F(P, 14), F(Q, 14)
    u = F(determinant + cross, ruler)
    v = F(-determinant, ruler)
    psi = (
        bernoulli3(u + a - b)
        + bernoulli3(u - a + b)
        + bernoulli3(v + a - b)
        + bernoulli3(v - a + b)
        - bernoulli3(u + a + b)
        - bernoulli3(u - a - b)
        - bernoulli3(v + a + b)
        - bernoulli3(v - a - b)
    )
    return F(1, 49) + F(ruler, 6 * P * Q * cross) * psi


def floor_moments(n: int, modulus: int, a: int, b: int):
    if n == 0:
        return 0, 0, 0
    first = n * (n - 1) // 2
    second = n * (n - 1) * (2 * n - 1) // 6
    qa, a0 = divmod(a, modulus)
    qb, b0 = divmod(b, modulus)
    base0 = qa * first + qb * n
    base1 = qa * second + qb * first
    base2 = qa * qa * second + 2 * qa * qb * first + qb * qb * n
    if a0 == 0:
        return base0, base1, base2
    height = (a0 * (n - 1) + b0) // modulus
    if height == 0:
        return base0, base1, base2
    u0, u1, u2 = floor_moments(height, a0, modulus, modulus - b0 + a0 - 1)
    r0 = n * height - u0
    r1 = height * first - (u2 - u0) // 2
    r2 = n * height * height - 2 * u1 - u0
    return (
        base0 + r0,
        base1 + r1,
        base2 + 2 * qa * r1 + 2 * qb * r0 + r2,
    )


def residue_prefix(n, modulus, a, b, threshold, base):
    shifted = floor_moments(n, modulus, a, b + modulus - threshold)
    d0, d1 = shifted[0] - base[0], shifted[1] - base[1]
    y0d = (shifted[2] - base[2] - d0) // 2
    high_sum = a * d1 + b * d0 - modulus * y0d
    total = a * n * (n - 1) // 2 + b * n - modulus * base[0]
    return n - d0, total - high_sum


def triangle_sum(n, modulus, a, b, peak, ruler, base, total):
    if peak <= 0:
        return 0
    radius = (peak - 1) // ruler
    # When 2*radius >= modulus the low and high residue ranges overlap, but
    # they represent two distinct neighbouring lifts and must both be added.
    # The formula remains exact until radius reaches a full modulus, which
    # cannot occur on the projective low-ratio bank (all ratios are <=6).
    require(radius < modulus, ("three triangle lifts", n, modulus, peak))
    low_count, low_sum = residue_prefix(n, modulus, a, b, radius + 1, base)
    before_count, before_sum = residue_prefix(
        n, modulus, a, b, modulus - radius, base
    )
    high_count, high_sum = n - before_count, total - before_sum
    return (
        peak * low_count
        - ruler * low_sum
        + (peak - ruler * modulus) * high_count
        + ruler * high_sum
    )


def physical_mass(ruler: int, cell: int, e: int, p: int, f: int, q: int) -> F:
    """Exact pair mass from Euclidean floor moments, including two-lift tails."""
    if p > q:
        return physical_mass(ruler, cell, f, q, e, p)
    z, w = ruler * p - e, ruler * q - f
    phase_e, phase_f = (e * cell) % ruler, (f * cell) % ruler
    determinant = phase_e * w - phase_f * z
    require(determinant % ruler == 0, ("nonintegral phase", ruler, cell, e, p, f, q))
    b, a = (determinant // ruler) % z, w % z
    base = floor_moments(p, z, a, b)
    total = a * p * (p - 1) // 2 + b * p - z * base[0]
    unit = ruler // 14
    outer, inner = unit * (z + w), unit * (w - z)
    numerator = triangle_sum(p, z, a, b, outer, ruler, base, total)
    numerator -= triangle_sum(p, z, a, b, inner, ruler, base, total)
    return F(numerator, z * w)


def body_geometry(body, ruler):
    check_ruler, ranges = R.safe_cell_ranges(body)
    require(check_ruler == ruler, (body, ruler, check_ruler))
    cells = tuple(j for left, right in ranges for j in range(left, right))
    cell = cells[len(cells) // 2]
    centres = tuple(
        i for i, label in enumerate(body) if (label * cell) % ruler == ruler // 14
    )
    require(len(centres) == 1, (body, ruler, cell, centres))
    centre = centres[0]
    second = next(i for i in range(6) if i != centre)
    leaves = tuple(i for i in range(6) if i not in (centre, second))
    positions = (centre, second, *leaves)
    edges = tuple(edge for edge in MS.EDGES if centre in edge or second in edge)
    require(len(edges) == 9, (body, centre, second, edges))
    return cell, centre, second, positions, edges


def run_census(progress: bool = False, probe_scale_max: int = 1):
    bodies = MS.body_universe()
    rays, common_hist = MS.projective_low_two_star_rays()
    require((len(bodies), len(rays), len(bodies) * len(rays)) == EXPECTED_COUNTS[:3],
            (len(bodies), len(rays)))
    weakest = None
    weakest_count = 0
    scale_one_weakest = None
    selected_groups = set()
    primitive_selected_groups = set()
    selector_edge_hist = Counter()
    ratio_hist = Counter()
    semantic = hashlib.sha256()
    half_limit_weakest = None

    for body_index, (body, ruler) in enumerate(bodies):
        cell, centre, second, positions, edges = body_geometry(body, ruler)
        words = []
        pair_keys = {edge: set() for edge in edges}
        for ray in rays:
            word_list = [0] * 6
            for position, value in zip(positions, ray):
                word_list[position] = value
            word = tuple(word_list)
            words.append(word)
            for i, j in edges:
                pair_keys[i, j].add((word[i], word[j]))

        limits = {}
        for edge, pairs in pair_keys.items():
            i, j = edge
            for p, q in pairs:
                limits[edge, p, q] = homogeneous_limit(
                    ruler, cell, body[i], p, body[j], q
                )

        masses = {}
        body_selected_groups = set()
        debt_terms = {
            (i, p): F(body[i], 7 * (p * ruler - body[i]))
            for i in range(6)
            for p in {value for word in words for value in word}
        }
        for ray, word in zip(rays, words):
            # The selector is lexicographic: maximize the exact homogeneous
            # limit and break exact ties by the lexicographically largest
            # edge.  The companion all-dilation theorem uses this frozen rule.
            selected_limit, selected_edge = max(
                (limits[edge, word[edge[0]], word[edge[1]]], edge)
                for edge in edges
            )
            debt = sum((debt_terms[i, word[i]] for i in range(6)), F(0))
            margin = selected_limit - debt
            require(margin > 0, ("limit selector failure", body, ray, selected_edge, margin))
            row = (
                margin,
                body,
                ruler,
                cell,
                centre,
                second,
                ray,
                selected_edge,
                selected_limit,
                debt,
            )
            if weakest is None or row < weakest:
                weakest, weakest_count = row, 1
            elif row == weakest:
                weakest_count += 1

            i, j = selected_edge
            mass_key = selected_edge, word[i], word[j]
            if mass_key not in masses:
                masses[mass_key] = physical_mass(
                    ruler, cell, body[i], word[i], body[j], word[j]
                )
            physical = masses[mass_key]
            physical_margin = physical - debt
            require(
                physical_margin > 0,
                ("scale-one selector failure", body, ray, selected_edge, physical, debt),
            )
            scale_row = (
                physical_margin,
                body,
                ruler,
                cell,
                centre,
                second,
                ray,
                selected_edge,
                physical,
                debt,
            )
            if scale_one_weakest is None or scale_row < scale_one_weakest:
                scale_one_weakest = scale_row

            pair_divisor = gcd(word[i], word[j])
            ratio_hist[(word[i] // pair_divisor, word[j] // pair_divisor)] += 1
            selector_edge_hist[(centre, second, selected_edge)] += 1
            selected_groups.add((body, selected_edge, word[i], word[j]))
            pair_divisor = gcd(word[i], word[j])
            primitive_selected_groups.add(
                (body, selected_edge, word[i] // pair_divisor, word[j] // pair_divisor)
            )
            body_selected_groups.add((selected_edge, word[i], word[j]))
            semantic.update(
                repr((body, ray, selected_edge, selected_limit, debt, physical)).encode()
            )

        for selected_edge, p, q in body_selected_groups:
            i, j = selected_edge
            limit = limits[selected_edge, p, q]
            for scale in range(2, probe_scale_max + 1):
                mass = physical_mass(
                    ruler, cell, body[i], scale * p, body[j], scale * q
                )
                row = (
                    mass - limit / 2,
                    body,
                    ruler,
                    cell,
                    selected_edge,
                    p,
                    q,
                    scale,
                    mass,
                    limit,
                )
                if half_limit_weakest is None or row < half_limit_weakest:
                    half_limit_weakest = row

        if progress and (body_index + 1) % 25 == 0:
            print(
                f"progress={body_index + 1}/{len(bodies)};"
                f"selected_groups={len(selected_groups)};weakest_gap={weakest[0]}",
                flush=True,
            )

    require(len(selected_groups) == EXPECTED_COUNTS[3], len(selected_groups))
    require(weakest is not None and scale_one_weakest is not None, "empty census")
    expected_row = EXPECTED_WEAKEST
    actual_row = (
        weakest[1], weakest[2], weakest[3], weakest[4], weakest[5],
        weakest[6], weakest[7], weakest[8], weakest[9], weakest[0],
    )
    require(actual_row == expected_row, ("weakest mismatch", actual_row))
    semantic_sha256 = semantic.hexdigest()
    require(
        semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
        ("semantic digest", semantic_sha256, EXPECTED_SEMANTIC_SHA256),
    )
    return {
        "body_count": len(bodies),
        "ray_count": len(rays),
        "assignment_count": len(bodies) * len(rays),
        "selected_group_count": len(selected_groups),
        "primitive_selected_group_count": len(primitive_selected_groups),
        "common_neighbour_hist": common_hist,
        "selector_edge_hist": tuple(sorted(selector_edge_hist.items())),
        "primitive_ratio_hist": tuple(sorted(ratio_hist.items())),
        "weakest": weakest,
        "weakest_count": weakest_count,
        "scale_one_weakest": scale_one_weakest,
        "selected_groups": tuple(sorted(selected_groups)),
        "half_limit_probe_scale_max": probe_scale_max,
        "half_limit_weakest": half_limit_weakest,
        "semantic_sha256": semantic_sha256,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--progress", action="store_true")
    parser.add_argument("--probe-scale-max", type=int, default=1)
    args = parser.parse_args()
    require(args.probe_scale_max >= 1, args.probe_scale_max)
    row = run_census(args.progress, args.probe_scale_max)
    weakest = row["weakest"]
    physical = row["scale_one_weakest"]
    print("LRC14 REFLECTED LOW TWO-STAR HOMOGENEOUS-LIMIT SELECTOR SCOUT")
    print("status=FINITE-EXACT;not an all-dilation or reflected-branch closure")
    print(
        f"bodies={row['body_count']};rays={row['ray_count']};"
        f"assignments={row['assignment_count']};failures=0"
    )
    print(f"selected_body_edge_channel_groups={row['selected_group_count']}")
    print(
        "selected_body_edge_primitive_channel_groups="
        f"{row['primitive_selected_group_count']}"
    )
    print(
        f"weakest_body={weakest[1]};L={weakest[2]};cell={weakest[3]};"
        f"centre={weakest[4]};second={weakest[5]};ray={weakest[6]};"
        f"edge={weakest[7]}"
    )
    print(
        f"weakest_limit={weakest[8]};debt={weakest[9]};"
        f"limit_gap={weakest[0]};multiplicity={row['weakest_count']}"
    )
    print(
        f"weakest_selected_scale_one_mass={physical[8]};"
        f"scale_one_margin={physical[0]}"
    )
    print(f"primitive_ratio_hist={row['primitive_ratio_hist']}")
    if row["half_limit_weakest"] is not None:
        half = row["half_limit_weakest"]
        print(
            f"half_limit_probe=2..{row['half_limit_probe_scale_max']};"
            f"weakest_margin={half[0]};body={half[1]};L={half[2]};"
            f"cell={half[3]};edge={half[4]};base_pair=({half[5]},{half[6]});"
            f"scale={half[7]};mass={half[8]};limit={half[9]}"
        )
    print(f"semantic_sha256={row['semantic_sha256']}")
    source_hash = hashlib.sha256(
        Path(__file__).read_bytes().replace(b"\r\n", b"\n")
    ).hexdigest()
    print(f"source_sha256={source_hash}")


if __name__ == "__main__":
    main()
