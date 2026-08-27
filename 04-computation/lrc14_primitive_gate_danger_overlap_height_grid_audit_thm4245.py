#!/usr/bin/env python3
"""Exact common-grid audit for the strengthened THM-4233 height obstruction.

Run from any directory with

    python3 lrc14_primitive_gate_height_grid_audit.py --repo /path/to/math

The script reconstructs the post-THM-4242 residual from the frozen THM-4231
K(q) table, applies every proved literal closure through THM-4242, and then
computes the THM-4233 pair observable by a common-grid endpoint/midpoint
method.  Every verification uses explicit ``require`` calls, so Python's
``-O`` flag cannot disable the audit.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import math
from fractions import Fraction
from pathlib import Path


THRESHOLD = Fraction(1650, 28710227)
POOL_EXPECTED = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
}
BOUNDARY_LITERAL = {
    (744, 824), (744, 822), (744, 821), (744, 820), (744, 818),
    (744, 817), (744, 815), (744, 814), (744, 813), (744, 812),
    (744, 811), (744, 810), (744, 809), (744, 805), (744, 803),
    (744, 800), (766, 800), (744, 798), (744, 794), (744, 793),
    (744, 791), (744, 790), (744, 789), (744, 787), (744, 780),
    (765, 780), (766, 780), (768, 780), (616, 777), (616, 776),
    (616, 775), (744, 775), (616, 774), (616, 773), (744, 773),
    (616, 772), (616, 771), (721, 771), (616, 770), (721, 770),
    (744, 770), (750, 770), (765, 770), (766, 770), (768, 770),
}
INHERITED_LITERAL = {(1, 542), (49, 50), (50, 51)}
THM4238_LITERAL = (
    {(6, r) for r in range(590, 614)}
    | {(25, r) for r in range(590, 598)}
)
THM4242_LITERAL = {(50, r) for r in range(590, 626)}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fnv_add(state: int, word: int) -> int:
    for shift in range(0, 64, 8):
        state ^= (word >> shift) & 0xFF
        state = (state * 0x100000001B3) & ((1 << 64) - 1)
    return state


def edge_fnv(edges: list[tuple[int, int]]) -> int:
    state = 0xCBF29CE484222325
    for a, b in edges:
        state = fnv_add(state, a)
        state = fnv_add(state, b)
    return state


def edge_sha256(edges: list[tuple[int, int]]) -> str:
    encoded = b"".join(f"{a},{b}\n".encode("ascii") for a, b in edges)
    return hashlib.sha256(encoded).hexdigest()


def ast_assignment(module: ast.Module, name: str):
    matches = []
    for node in module.body:
        if isinstance(node, ast.Assign):
            if any(isinstance(target, ast.Name) and target.id == name
                   for target in node.targets):
                matches.append(ast.literal_eval(node.value))
    require(len(matches) == 1, f"expected one literal assignment for {name}")
    return matches[0]


def reconstruct_residual(repo: Path) -> list[tuple[int, int]]:
    source_path = repo / (
        "04-computation/"
        "lrc14_all_fixed_outsider_ray_symmetric_postprocess_thm4231.py"
    )
    require(source_path.is_file(), f"missing canonical K(q) source: {source_path}")
    source = source_path.read_text(encoding="utf-8")
    module = ast.parse(source, filename=str(source_path))
    pool = set(ast_assignment(module, "POOL"))
    values_text = ast_assignment(module, "VALUES")
    require(pool == POOL_EXPECTED, "THM-4231 pool changed")
    require(isinstance(values_text, str), "VALUES is not a literal string")

    qs = [q for q in range(1, 1307) if q not in pool]
    values = [int(word) for word in values_text.split()]
    require(len(qs) == 1276, "outsider ray count changed")
    require(len(values) == len(qs), "K(q) table length changed")
    kappa = dict(zip(qs, values, strict=True))

    original = [
        (a, b)
        for a in qs
        for b in qs
        if a < b and b < kappa[a] and a < kappa[b]
    ]
    require(len(original) == 181242, "THM-4231 certificate residual count changed")
    require(edge_fnv(original) == 0x8A4E1370FB023907,
            "THM-4231 certificate residual FNV changed")

    original_set = set(original)
    require(len(BOUNDARY_LITERAL) == 45, "boundary-literal set size changed")
    require(BOUNDARY_LITERAL == {e for e in original if e[1] >= 770},
            "boundary-literal layer changed")
    require(INHERITED_LITERAL <= original_set, "inherited literal missing")
    require(BOUNDARY_LITERAL.isdisjoint(INHERITED_LITERAL),
            "literal closure classes overlap")

    thm4231 = [
        edge for edge in original
        if edge not in BOUNDARY_LITERAL and edge not in INHERITED_LITERAL
    ]
    require(len(thm4231) == 181194, "patched THM-4231 count changed")
    require(edge_fnv(thm4231) == 0x3874FECAC4ECBD8A,
            "patched THM-4231 FNV changed")
    require(THM4238_LITERAL <= set(thm4231), "THM-4238 closure missing")

    thm4238 = [edge for edge in thm4231 if edge not in THM4238_LITERAL]
    require(len(thm4238) == 181162, "post-THM-4238 count changed")
    require(edge_fnv(thm4238) == 0x7E5F6AF58A370E3A,
            "post-THM-4238 FNV changed")
    require(THM4242_LITERAL <= set(thm4238), "THM-4242 closure missing")

    residual = [edge for edge in thm4238 if edge not in THM4242_LITERAL]
    require(len(residual) == 181126, "post-THM-4242 count changed")
    require(edge_fnv(residual) == 0xBDF59726990A6C92,
            "post-THM-4242 FNV changed")
    require(edge_sha256(residual) ==
            "c0e2fe1c69cfe8cfe6e633a1eca0d8d37ca991ecdaa04b98d7c595a99b9be6bf",
            "post-THM-4242 SHA-256 changed")
    require(max(b for _, b in residual) == 769, "residual height changed")
    require([e for e in residual if e[1] == 769] == [(616, 769), (721, 769)],
            "top residual layer changed")
    fixed50 = [e for e in residual if 50 in e]
    require(len(fixed50) == 556, "fixed-50 remainder count changed")
    require(max(b if a == 50 else a for a, b in fixed50) == 589,
            "fixed-50 remainder height changed")
    return residual


def midpoint_safe(speed: int, grid: int, left: int, right: int) -> bool:
    residue = speed * (left + right) % (2 * grid)
    return grid <= 7 * residue <= 13 * grid


def exact_pair_grid(u: int, v: int) -> tuple[int, int, int, int]:
    require(u > 0 and v > 0 and u != v, "invalid primitive pair")
    require(math.gcd(u, v) == 1, "pair was not reduced")
    grid = 14 * u * v
    points = {0, grid}
    points.update((v * (14 * i + sign)) % grid
                  for i in range(u) for sign in (-1, 1))
    points.update((u * (14 * j + sign)) % grid
                  for j in range(v) for sign in (-1, 1))
    ordered = sorted(points)
    require(ordered[0] == 0 and ordered[-1] == grid,
            "common grid lost a sentinel")

    cells: list[tuple[int, bool]] = []
    safe_ticks = 0
    components = 0
    previous_safe = False
    first_safe = False
    for index, (left, right) in enumerate(zip(ordered, ordered[1:])):
        safe = (midpoint_safe(u, grid, left, right)
                and midpoint_safe(v, grid, left, right))
        length = right - left
        require(length > 0, "nonpositive common-grid cell")
        cells.append((length, safe))
        safe_ticks += length * int(safe)
        if index == 0:
            first_safe = safe
        if safe and not previous_safe:
            components += 1
        previous_safe = safe
    if first_safe and previous_safe:
        components -= 1
    require(0 < safe_ticks < grid, "degenerate joint safe set")
    require(0 < components <= u + v - 1,
            "shared-danger-component bound failed")
    require(14 * safe_ticks <= 11 * grid,
            "primitive density upper bound failed")

    primitive = 0
    minimum = 0
    maximum = 0
    for length, safe in cells:
        primitive += length * (grid * int(safe) - safe_ticks)
        minimum = min(minimum, primitive)
        maximum = max(maximum, primitive)
    require(primitive == 0, "centered primitive did not close")
    return grid, safe_ticks, maximum - minimum, components


def scan(
    residual: list[tuple[int, int]],
) -> tuple[int, int, int, Fraction, tuple[int, ...]]:
    cache: dict[tuple[int, int], tuple[int, int, int, int]] = {}
    density_failures = 0
    oscillation_passes = 0
    best: tuple[int, ...] | None = None
    for a, b in residual:
        dilation = math.gcd(a, b)
        key = (a // dilation, b // dilation)
        if key not in cache:
            cache[key] = exact_pair_grid(*key)
        grid, safe_ticks, oscillation_raw, components = cache[key]
        density_margin = 91 * safe_ticks - 66 * grid
        if density_margin < 0:
            density_failures += 1
            continue
        lhs = oscillation_raw * THRESHOLD.denominator
        rhs = THRESHOLD.numerator * dilation * grid * grid
        if lhs <= rhs:
            oscillation_passes += 1
        row = (
            lhs, rhs, a, b, dilation, key[0], key[1], grid, safe_ticks,
            oscillation_raw, components, density_margin,
        )
        if best is None or row[0] * best[1] < best[0] * row[1]:
            best = row
    require(best is not None, "no density-admissible residual edge")
    for (u, v), (grid, safe_ticks, _, components) in cache.items():
        require(14 * safe_ticks <= 11 * grid,
                f"beta upper bound failed for primitive pair {(u, v)}")
        require(components <= u + v - 1,
                f"component bound failed for primitive pair {(u, v)}")
    maximum_beta = max(Fraction(row[1], row[0]) for row in cache.values())
    return (len(cache), density_failures, oscillation_passes,
            maximum_beta, best)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.resolve()

    # Canonical THM-4233 controls, checked by this implementation.
    n_113, s_113, o_113, _ = exact_pair_grid(1, 13)
    n_12, s_12, o_12, c_12 = exact_pair_grid(1, 2)
    n_cert, s_cert, o_cert, _ = exact_pair_grid(3713, 5149)
    require(Fraction(s_113, n_113) == Fraction(66, 91),
            "(1,13) density control failed")
    require(Fraction(o_113, n_113 * n_113) == Fraction(990, 8281),
            "(1,13) oscillation control failed")
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
    (ratio_count, density_failures, oscillation_passes,
     maximum_beta, best) = scan(residual)
    require(ratio_count == 115429, "primitive-ratio count changed")
    require(maximum_beta == Fraction(11, 14),
            "full-universe beta maximum changed")
    require(density_failures == 0, "unexpected density-gate failure")
    require(oscillation_passes == 0, "a residual edge passed the gate")

    (lhs, rhs, a, b, dilation, u, v, grid, safe_ticks,
     oscillation_raw, components, density_margin) = best
    beta = Fraction(safe_ticks, grid)
    omega = Fraction(oscillation_raw, grid * grid)
    omega_over_g = omega / dilation
    factor = Fraction(lhs, rhs)
    require((a, b, dilation, u, v) == (466, 699, 233, 2, 3),
            "closest residual edge changed")
    require(beta == Fraction(16, 21), "closest density changed")
    require(omega == Fraction(67, 882), "closest oscillation changed")
    require(omega_over_g == Fraction(67, 205506),
            "closest scaled oscillation changed")
    require(factor == Fraction(39256841, 6920100),
            "closest gate factor changed")

    lower = Fraction(33, 196 * 1536)
    threshold_height = Fraction(33, 196) / THRESHOLD
    integer_height = (
        threshold_height.numerator + threshold_height.denominator - 1
    ) // threshold_height.denominator
    integer_sum = integer_height + 1
    integer_max = (integer_sum + 1) // 2
    require(threshold_height == Fraction(585923, 200),
            "analytic a+b-g threshold changed")
    require((integer_height, integer_sum, integer_max) == (2930, 2931, 1466),
            "integer height consequences changed")
    require(lower - THRESHOLD == Fraction(3065953, 58798544896),
            "height-cap strict margin changed")

    print("LRC14_PRIMITIVE_GATE_HEIGHT_OBSTRUCTION_AUDIT")
    print("PATH COMMON_GRID_ENDPOINT_MIDPOINT EXPLICIT_REQUIRE")
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
        f"OSCILLATION_GATE_PASSES {oscillation_passes}"
    )
    print(
        f"CLOSEST EDGE {a},{b} G {dilation} PRIMITIVE {u},{v} "
        f"COMPONENTS {components} BETA {beta} DENSITY_MARGIN_NUM {density_margin} "
        f"OMEGA {omega} OMEGA_OVER_G {omega_over_g} "
        f"FACTOR_OVER_GATE {factor}"
    )
    print(
        f"ANALYTIC MIN_A_PLUS_B_MINUS_G {threshold_height} "
        f"INTEGER_MIN_A_PLUS_B_MINUS_G {integer_height} "
        f"INTEGER_MIN_SUM {integer_sum} INTEGER_MIN_MAX {integer_max} "
        f"RESIDUAL_MAX_A_PLUS_B_MINUS_G 1536 "
        f"HEIGHT_MARGIN {lower - THRESHOLD}"
    )
    print("VERDICT NO_POST_THM4242_RESIDUAL_EDGE_CAN_PASS_THM4233_GATE")


if __name__ == "__main__":
    main()
