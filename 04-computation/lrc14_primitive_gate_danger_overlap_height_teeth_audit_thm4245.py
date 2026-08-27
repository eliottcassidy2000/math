#!/usr/bin/env python3
"""Independent ordered-safe-teeth audit for the strengthened height obstruction.

This implementation does not construct the common endpoint arrangement and
does not classify midpoints.  It intersects the two ordered safe-tooth lists,
integrates the centered primitive directly at rational-grid endpoints, and
reconstructs the post-THM-4242 graph using an independent text-table parser.
All checks are explicit and remain active under Python ``-O``.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import math
import re
from fractions import Fraction
from pathlib import Path


THRESHOLD = Fraction(1650, 28710227)
POOL = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
}
LITERAL_4231 = {
    (744, 824), (744, 822), (744, 821), (744, 820), (744, 818),
    (744, 817), (744, 815), (744, 814), (744, 813), (744, 812),
    (744, 811), (744, 810), (744, 809), (744, 805), (744, 803),
    (744, 800), (766, 800), (744, 798), (744, 794), (744, 793),
    (744, 791), (744, 790), (744, 789), (744, 787), (744, 780),
    (765, 780), (766, 780), (768, 780), (616, 777), (616, 776),
    (616, 775), (744, 775), (616, 774), (616, 773), (744, 773),
    (616, 772), (616, 771), (721, 771), (616, 770), (721, 770),
    (744, 770), (750, 770), (765, 770), (766, 770), (768, 770),
    (1, 542), (49, 50), (50, 51),
}
LITERAL_4238 = (
    {(6, r) for r in range(590, 614)}
    | {(25, r) for r in range(590, 598)}
)
LITERAL_4242 = {(50, r) for r in range(590, 626)}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def independent_fnv(edges: list[tuple[int, int]]) -> int:
    state = 0xCBF29CE484222325
    for edge in edges:
        for word in edge:
            require(word >= 0, "negative edge label")
            for byte in word.to_bytes(8, byteorder="little", signed=False):
                state ^= byte
                state = (state * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
    return state


def independent_sha256(edges: list[tuple[int, int]]) -> str:
    digest = hashlib.sha256()
    for a, b in edges:
        digest.update(str(a).encode("ascii"))
        digest.update(b",")
        digest.update(str(b).encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def reconstruct_residual(repo: Path) -> list[tuple[int, int]]:
    table_path = repo / (
        "04-computation/"
        "lrc14_all_fixed_outsider_ray_symmetric_postprocess_thm4231.py"
    )
    require(table_path.is_file(), f"missing K(q) table: {table_path}")
    text = table_path.read_text(encoding="utf-8")
    matches = re.findall(r'VALUES\s*=\s*"""(.*?)"""', text, flags=re.DOTALL)
    require(len(matches) == 1, "could not isolate the literal K(q) table")
    values = [int(token) for token in re.findall(r"[0-9]+", matches[0])]
    outsiders = [value for value in range(1, 1307) if value not in POOL]
    require(len(outsiders) == 1276, "outsider universe changed")
    require(len(values) == 1276, "K(q) value count changed")
    kappa = {outsider: values[index] for index, outsider in enumerate(outsiders)}

    original = []
    for a, b in itertools.combinations(outsiders, 2):
        if b < kappa[a] and a < kappa[b]:
            original.append((a, b))
    require(len(original) == 181242, "raw residual count changed")
    require(independent_fnv(original) == 0x8A4E1370FB023907,
            "raw residual FNV changed")
    require(len(LITERAL_4231) == 48, "THM-4231 literal closure count changed")
    require(LITERAL_4231 <= set(original), "THM-4231 literal edge missing")

    after_4231 = sorted(set(original).difference(LITERAL_4231))
    require(len(after_4231) == 181194, "post-THM-4231 count changed")
    require(independent_fnv(after_4231) == 0x3874FECAC4ECBD8A,
            "post-THM-4231 FNV changed")
    require(len(LITERAL_4238) == 32 and LITERAL_4238 <= set(after_4231),
            "THM-4238 closure changed")

    after_4238 = sorted(set(after_4231).difference(LITERAL_4238))
    require(len(after_4238) == 181162, "post-THM-4238 count changed")
    require(independent_fnv(after_4238) == 0x7E5F6AF58A370E3A,
            "post-THM-4238 FNV changed")
    require(len(LITERAL_4242) == 36 and LITERAL_4242 <= set(after_4238),
            "THM-4242 closure changed")

    residual = sorted(set(after_4238).difference(LITERAL_4242))
    require(len(residual) == 181126, "post-THM-4242 residual count changed")
    require(independent_fnv(residual) == 0xBDF59726990A6C92,
            "post-THM-4242 residual FNV changed")
    require(independent_sha256(residual) ==
            "c0e2fe1c69cfe8cfe6e633a1eca0d8d37ca991ecdaa04b98d7c595a99b9be6bf",
            "post-THM-4242 residual SHA-256 changed")

    layers: dict[int, list[tuple[int, int]]] = {}
    for edge in residual:
        layers.setdefault(edge[1], []).append(edge)
    require(max(layers) == 769, "residual maximum endpoint changed")
    require(layers[769] == [(616, 769), (721, 769)],
            "residual top layer changed")
    fixed = [edge for edge in residual if edge[0] == 50 or edge[1] == 50]
    require(len(fixed) == 556, "fixed-50 residual count changed")
    require(max(edge[1] if edge[0] == 50 else edge[0] for edge in fixed) == 589,
            "fixed-50 residual maximum changed")
    return residual


def exact_pair_teeth(u: int, v: int) -> tuple[int, int, int, int]:
    require(u > 0 and v > 0 and u != v, "invalid reduced pair")
    require(math.gcd(u, v) == 1, "pair is not primitive")
    grid = 14 * u * v

    # On the common grid, these are already ordered safe intervals.  No wall
    # sorting or midpoint decision is used.
    u_teeth = [(v * (14 * i + 1), v * (14 * i + 13)) for i in range(u)]
    v_teeth = [(u * (14 * j + 1), u * (14 * j + 13)) for j in range(v)]
    i = 0
    j = 0
    components: list[tuple[int, int]] = []
    while i < len(u_teeth) and j < len(v_teeth):
        low = max(u_teeth[i][0], v_teeth[j][0])
        high = min(u_teeth[i][1], v_teeth[j][1])
        if low < high:
            if components and components[-1][1] >= low:
                if high > components[-1][1]:
                    components[-1] = (components[-1][0], high)
            else:
                components.append((low, high))
        if u_teeth[i][1] < v_teeth[j][1]:
            i += 1
        elif v_teeth[j][1] < u_teeth[i][1]:
            j += 1
        else:
            i += 1
            j += 1

    require(0 < len(components) <= u + v - 1,
            "shared-danger-component census failed")
    for index, (low, high) in enumerate(components):
        require(0 <= low < high <= grid, "invalid joint-safe component")
        if index:
            require(components[index - 1][1] < low,
                    "joint-safe components were not maximally merged")

    safe_ticks = sum(high - low for low, high in components)
    require(0 < safe_ticks < grid, "degenerate joint-safe density")
    require(14 * safe_ticks <= 11 * grid,
            "primitive density upper bound failed")
    covered = 0
    minimum = 0
    maximum = 0
    for low, high in components:
        at_low = grid * covered - safe_ticks * low
        covered += high - low
        at_high = grid * covered - safe_ticks * high
        minimum = min(minimum, at_low, at_high)
        maximum = max(maximum, at_low, at_high)
    require(covered == safe_ticks, "covered length changed")
    require(grid * covered - safe_ticks * grid == 0,
            "centered primitive did not return to zero")
    return grid, safe_ticks, maximum - minimum, len(components)


def scan(
    residual: list[tuple[int, int]],
) -> tuple[int, int, int, Fraction, tuple[int, ...]]:
    cache: dict[tuple[int, int], tuple[int, int, int, int]] = {}
    density_failures = 0
    gate_passes = 0
    closest: tuple[int, ...] | None = None
    for a, b in residual:
        dilation = math.gcd(a, b)
        u = a // dilation
        v = b // dilation
        key = (u, v)
        if key not in cache:
            cache[key] = exact_pair_teeth(u, v)
        grid, safe_ticks, oscillation_raw, component_count = cache[key]
        density_margin = 91 * safe_ticks - 66 * grid
        if density_margin < 0:
            density_failures += 1
            continue
        left = oscillation_raw * THRESHOLD.denominator
        right = THRESHOLD.numerator * dilation * grid * grid
        if left <= right:
            gate_passes += 1
        record = (
            left, right, a, b, dilation, u, v, grid, safe_ticks,
            oscillation_raw, component_count, density_margin,
        )
        if closest is None or record[0] * closest[1] < closest[0] * record[1]:
            closest = record
    require(closest is not None, "no closest density-admissible edge")
    for (u, v), (grid, safe_ticks, _, component_count) in cache.items():
        require(14 * safe_ticks <= 11 * grid,
                f"beta upper bound failed for primitive pair {(u, v)}")
        require(component_count <= u + v - 1,
                f"component bound failed for primitive pair {(u, v)}")
    maximum_beta = max(Fraction(row[1], row[0]) for row in cache.values())
    return len(cache), density_failures, gate_passes, maximum_beta, closest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    options = parser.parse_args()
    repo = options.repo.resolve()

    n_113, s_113, o_113, _ = exact_pair_teeth(1, 13)
    n_12, s_12, o_12, c_12 = exact_pair_teeth(1, 2)
    n_cert, s_cert, o_cert, _ = exact_pair_teeth(3713, 5149)
    require(Fraction(s_113, n_113) == Fraction(66, 91),
            "(1,13) density hostile failed")
    require(Fraction(o_113, n_113 * n_113) == Fraction(990, 8281),
            "(1,13) oscillation hostile failed")
    require(Fraction(s_12, n_12) == Fraction(11, 14),
            "(1,2) density-equality control failed")
    require(c_12 == 2, "(1,2) shared-component control failed")
    require(Fraction(o_12, n_12 * n_12) == Fraction(11, 98),
            "(1,2) oscillation control failed")
    cert_beta = Fraction(s_cert, n_cert)
    cert_omega = Fraction(o_cert, n_cert * n_cert)
    require(cert_beta == Fraction(98322360, 133827659),
            "(3713,5149) density certificate changed")
    require(cert_omega == Fraction(277071798, 4823550313337),
            "(3713,5149) oscillation certificate changed")
    require(cert_beta >= Fraction(66, 91),
            "(3713,5149) failed the density gate")
    require(cert_omega <= THRESHOLD,
            "(3713,5149) failed the oscillation gate")

    residual = reconstruct_residual(repo)
    (ratio_count, density_failures, gate_passes,
     maximum_beta, closest) = scan(residual)
    require(ratio_count == 115429, "primitive-ratio count changed")
    require(maximum_beta == Fraction(11, 14),
            "full-universe beta maximum changed")
    require(density_failures == 0, "density failure count changed")
    require(gate_passes == 0, "a post-THM-4242 edge passed")

    (left, right, a, b, dilation, u, v, grid, safe_ticks,
     oscillation_raw, component_count, density_margin) = closest
    beta = Fraction(safe_ticks, grid)
    omega = Fraction(oscillation_raw, grid * grid)
    omega_over_g = omega / dilation
    factor = Fraction(left, right)
    require((a, b, dilation, u, v) == (466, 699, 233, 2, 3),
            "closest edge changed")
    require(component_count == 4, "closest component count changed")
    require(beta == Fraction(16, 21), "closest beta changed")
    require(omega == Fraction(67, 882), "closest omega changed")
    require(omega_over_g == Fraction(67, 205506),
            "closest omega/g changed")
    require(factor == Fraction(39256841, 6920100),
            "closest gate factor changed")

    minimum_height = Fraction(33, 196) / THRESHOLD
    cap_lower = Fraction(33, 196 * 1536)
    integer_height = (
        minimum_height.numerator + minimum_height.denominator - 1
    ) // minimum_height.denominator
    integer_sum = integer_height + 1
    integer_max = (integer_sum + 1) // 2
    require(minimum_height == Fraction(585923, 200),
            "analytic a+b-g threshold changed")
    require((integer_height, integer_sum, integer_max) == (2930, 2931, 1466),
            "integer height consequences changed")
    require(cap_lower - THRESHOLD == Fraction(3065953, 58798544896),
            "height-cap margin changed")

    print("LRC14_PRIMITIVE_GATE_HEIGHT_OBSTRUCTION_AUDIT")
    print("PATH ORDERED_SAFE_TEETH_INTERSECTION EXPLICIT_REQUIRE")
    print(
        "CONTROLS PAIR_1_13 PASS PAIR_1_2_BETA_EQUALITY PASS "
        "PAIR_3713_5149_GATE PASS"
    )
    print(
        "RESIDUAL COUNT 181126 FNV bdf59726990a6c92 "
        "SHA256 c0e2fe1c69cfe8cfe6e633a1eca0d8d37ca991ecdaa04b98d7c595a99b9be6bf"
    )
    print(
        f"PRIMITIVE_RATIOS {ratio_count} BETA_UPPER_CHECKS {ratio_count} "
        f"BETA_UPPER_VIOLATIONS 0 BETA_MAX {maximum_beta} "
        f"COMPONENT_BOUND_CHECKS {ratio_count} DENSITY_FAILURES {density_failures} "
        f"OSCILLATION_GATE_PASSES {gate_passes}"
    )
    print(
        f"CLOSEST EDGE {a},{b} G {dilation} PRIMITIVE {u},{v} "
        f"COMPONENTS {component_count} BETA {beta} DENSITY_MARGIN_NUM {density_margin} "
        f"OMEGA {omega} OMEGA_OVER_G {omega_over_g} "
        f"FACTOR_OVER_GATE {factor}"
    )
    print(
        f"ANALYTIC MIN_A_PLUS_B_MINUS_G {minimum_height} "
        f"INTEGER_MIN_A_PLUS_B_MINUS_G {integer_height} "
        f"INTEGER_MIN_SUM {integer_sum} INTEGER_MIN_MAX {integer_max} "
        f"RESIDUAL_MAX_A_PLUS_B_MINUS_G 1536 "
        f"HEIGHT_MARGIN {cap_lower - THRESHOLD}"
    )
    print("VERDICT NO_POST_THM4242_RESIDUAL_EDGE_CAN_PASS_THM4233_GATE")


if __name__ == "__main__":
    main()
