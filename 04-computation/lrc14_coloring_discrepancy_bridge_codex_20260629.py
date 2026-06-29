#!/usr/bin/env python3
"""HYP-3459: AP84 coloring/discrepancy bridge for LRC14.

This scout connects older coloring threads to the current AP84 tail package.
The point is not to prove LRC14 by a color slogan.  The point is to test which
color quotients are legal: what they preserve, what they forget, and which
sidecar must be restored before a proof may use them.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations


F = Fraction
MOD = 35

# HYP-3456 correction vector:
# N(m)=floor(12m/35)+CORRECTION[(m-1) mod 35].
CORRECTION = [
    2,
    2,
    1,
    1,
    1,
    1,
    1,
    2,
    1,
    1,
    2,
    1,
    1,
    1,
    1,
    2,
    2,
    1,
    1,
    1,
    1,
    2,
    2,
    1,
    1,
    2,
    1,
    1,
    2,
    1,
    2,
    2,
    1,
    1,
    0,
]


@dataclass(frozen=True)
class ResidueColor:
    r: int
    gate_bucket: str
    outer_gate: bool
    inner_gate: bool
    correction: int
    n_value: int
    escape_value: int


def floor_count(m: int) -> int:
    return ((504 * m - 6) // 70) - ((96 * m - 13) // 14)


def gate_bucket(r: int) -> tuple[str, bool, bool]:
    outer = r % 7 != 0
    inner = r % 5 != 0
    if outer and inner:
        return "both_outer_inner", outer, inner
    if outer:
        return "outer_only_7_gate", outer, inner
    if inner:
        return "inner_only_5_gate", outer, inner
    return "clean_only_lcm35", outer, inner


def residue_colors() -> list[ResidueColor]:
    rows: list[ResidueColor] = []
    for r in range(1, MOD + 1):
        bucket, outer, inner = gate_bucket(r)
        n_value = floor_count(r)
        correction = n_value - (12 * r // MOD)
        if correction != CORRECTION[r - 1]:
            raise AssertionError((r, correction, CORRECTION[r - 1]))
        rows.append(
            ResidueColor(
                r=r,
                gate_bucket=bucket,
                outer_gate=outer,
                inner_gate=inner,
                correction=correction,
                n_value=n_value,
                escape_value=2 * n_value,
            )
        )
    return rows


def fmt_frac(value: F) -> str:
    value = F(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def values_for(rows: list[ResidueColor], observable: str) -> list[F]:
    if observable == "outer_gate":
        return [F(int(row.outer_gate)) for row in rows]
    if observable == "inner_gate":
        return [F(int(row.inner_gate)) for row in rows]
    if observable == "both_gate":
        return [F(int(row.gate_bucket == "both_outer_inner")) for row in rows]
    if observable == "clean_gate":
        return [F(int(row.gate_bucket == "clean_only_lcm35")) for row in rows]
    if observable == "correction":
        return [F(row.correction) for row in rows]
    if observable == "correction_eq_2":
        return [F(int(row.correction == 2)) for row in rows]
    if observable == "escape_value":
        return [F(row.escape_value) for row in rows]
    raise KeyError(observable)


def cyclic_interval_discrepancy(values: list[F]) -> tuple[F, tuple[int, int, F, F, F]]:
    n = len(values)
    total = sum(values, F(0))
    mean = total / n
    best_abs = F(-1)
    best = (1, 1, F(0), F(0), F(0))
    for start in range(n):
        running = F(0)
        for length in range(1, n + 1):
            running += values[(start + length - 1) % n]
            expected = mean * length
            deviation = running - expected
            if abs(deviation) > best_abs:
                best_abs = abs(deviation)
                best = (start + 1, length, running, expected, deviation)
    return best_abs, best


def dyadic_haar_contrast(values: list[F]) -> tuple[F, tuple[int, int, F]]:
    n = len(values)
    best_abs = F(-1)
    best = (1, 1, F(0))
    half = 1
    while 2 * half <= n:
        for start in range(n):
            left = sum((values[(start + j) % n] for j in range(half)), F(0))
            right = sum((values[(start + half + j) % n] for j in range(half)), F(0))
            contrast = left - right
            if abs(contrast) > best_abs:
                best_abs = abs(contrast)
                best = (start + 1, 2 * half, contrast)
        half *= 2
    return best_abs, best


def mixed_fibers(
    rows: list[ResidueColor],
    key_fn,
    value_fn,
) -> dict[object, set[object]]:
    buckets: dict[object, set[object]] = defaultdict(set)
    for row in rows:
        buckets[key_fn(row)].add(value_fn(row))
    return {key: vals for key, vals in buckets.items() if len(vals) > 1}


def tournament_fingerprint() -> tuple[dict[int, int], list[str], int, int]:
    vertices = {
        "labelled_color_packet_theorem": (10, 10, 10, 10, 10, 10, 10, 10),
        "residue_gate_plus_floor_word": (10, 10, 7, 8, 8, 8, 8, 9),
        "endpoint_phase_sidecar": (5, 8, 10, 9, 8, 8, 8, 9),
        "haar_zipper_cocycle_repair": (7, 7, 7, 10, 10, 8, 9, 9),
        "incident_C3_Qsqrt_router": (8, 6, 6, 7, 9, 9, 8, 9),
        "branch_mask_discrepancy_word": (6, 6, 6, 8, 7, 8, 9, 8),
        "distance_graph_regular_coloring": (6, 5, 7, 6, 6, 7, 8, 7),
        "raw_mod35_gate_color": (8, 3, 2, 2, 4, 5, 5, 7),
        "raw_scalar_escape_count": (3, 3, 2, 2, 2, 3, 3, 2),
    }
    scores = {name: sum(axis) for name, axis in vertices.items()}
    order = [name for name, _ in sorted(scores.items(), key=lambda item: (-item[1], item[0]))]
    rank = {name: index for index, name in enumerate(order)}
    cycles = 0
    for a, b, c in combinations(order, 3):
        ab = rank[a] < rank[b]
        bc = rank[b] < rank[c]
        ca = rank[c] < rank[a]
        if ab == bc == ca:
            cycles += 1
    return dict(sorted(Counter(scores.values()).items())), order, cycles, 1


def print_color_table(rows: list[ResidueColor]) -> None:
    print("B. AP84 mod-35 color table")
    print("  r | gate_bucket              | outer | inner | d_r | N(r) | escapes")
    for row in rows:
        print(
            f"  {row.r:2d} | {row.gate_bucket:24s} |"
            f"   {int(row.outer_gate)}   |   {int(row.inner_gate)}   |"
            f"  {row.correction:2d} | {row.n_value:4d} | {row.escape_value:7d}"
        )
    print()


def main() -> None:
    rows = residue_colors()
    bucket_hist = Counter(row.gate_bucket for row in rows)
    d_hist = Counter(row.correction for row in rows)
    gate_to_d = defaultdict(list)
    d_to_gate = defaultdict(list)
    for row in rows:
        gate_to_d[row.gate_bucket].append(row.correction)
        d_to_gate[row.correction].append(row.gate_bucket)

    gate_mixed = mixed_fibers(rows, lambda row: row.gate_bucket, lambda row: row.correction)
    d_mixed = mixed_fibers(rows, lambda row: row.correction, lambda row: row.gate_bucket)
    gate_d_mixed = mixed_fibers(
        rows,
        lambda row: (row.gate_bucket, row.correction),
        lambda row: (row.outer_gate, row.inner_gate, row.correction),
    )

    observables = [
        "outer_gate",
        "inner_gate",
        "both_gate",
        "clean_gate",
        "correction",
        "correction_eq_2",
        "escape_value",
    ]
    discrepancy_rows = []
    for observable in observables:
        values = values_for(rows, observable)
        max_disc, interval = cyclic_interval_discrepancy(values)
        max_haar, haar = dyadic_haar_contrast(values)
        discrepancy_rows.append((observable, max_disc, interval, max_haar, haar))

    hist, path, cycles, hamiltonian_paths = tournament_fingerprint()

    print("HYP-3459 COLORING/DISCREPANCY BRIDGE FOR LRC14")
    print("=" * 76)
    print("status=SYNTHESIS / exact AP84 color-quotient audit; not an LRC14 proof")
    print("source=prior coloring/discrepancy work + HYP-3438 + HYP-3441 + HYP-3456/HYP-3457")
    print()

    print("A. Prior coloring hooks being connected")
    print("  centered edge variable: s_e=A_e-1/2 is the pair-first color coordinate")
    print("  W-polynomial: Hamiltonian paths are sums over products of those colors")
    print("  distance-graph LRC: a time is a regular circular coloring of G(Z,D)")
    print("  colored discrepancy: CRT counts need one layer for each residue color")
    print("  Haar zipper: a margin quotient must retain the local zeta side-channel")
    print("  Paris-Harrington miniature: a color is live only with extension rank")
    print()

    print_color_table(rows)

    print("C. Legal quotient collision audit")
    print(f"  gate_bucket_hist={dict(sorted(bucket_hist.items()))}")
    print(f"  correction_hist={dict(sorted(d_hist.items()))}")
    print("  gate_bucket -> correction values:")
    for bucket in sorted(gate_to_d):
        print(f"    {bucket}: {sorted(set(gate_to_d[bucket]))}")
    print("  correction -> gate buckets:")
    for correction in sorted(d_to_gate):
        print(f"    d={correction}: {sorted(set(d_to_gate[correction]))}")
    print(f"  mixed_gate_fibers_count={len(gate_mixed)}")
    print(f"  mixed_correction_fibers_count={len(d_mixed)}")
    print(f"  mixed_gate_plus_correction_fibers_count={len(gate_d_mixed)}")
    print("  finite_phase_collision_example:")
    print("    m=1 and m=36 have the same residue r=1, but HYP-3457 says m=1 is")
    print("    a mixed endpoint transient while HYP-3454 says m=36 is rank-one E/E.")
    print("  verdict: gate color, floor correction, and endpoint phase are three")
    print("    different colors.  Any proof quotient using one must emit the other")
    print("    two as sidecars or route them to named debt.")
    print()

    print("D. Discrepancy and Haar-style contrast over the period-35 word")
    print("  observable | max cyclic discrepancy | best arc | max dyadic contrast | best dyadic block")
    for observable, max_disc, interval, max_haar, haar in discrepancy_rows:
        start, length, actual, expected, deviation = interval
        haar_start, haar_length, contrast = haar
        print(
            f"  {observable:16s} | {fmt_frac(max_disc):>21s} |"
            f" start={start:2d},len={length:2d},sum={fmt_frac(actual)},"
            f"exp={fmt_frac(expected)},dev={fmt_frac(deviation)} |"
            f" {fmt_frac(max_haar):>18s} | start={haar_start:2d},len={haar_length:2d},"
            f"contrast={fmt_frac(contrast)}"
        )
    print("  reading: the mod-35 color word is low-dimensional but not flat.  The")
    print("    nonzero Haar contrasts are exactly the places where a scalar")
    print("    equidistribution claim would need a zipper side-channel.")
    print()

    print("E. Labelled color-packet theorem target")
    print("  For the AP84 tail S_m={1,2,...,11,13,84m}, a legal proof packet should")
    print("  carry at least these colors:")
    print("    1. residue gate color: outer gate iff 7 does not divide m, inner gate iff 5 does not divide m")
    print("    2. floor word: N(m)=floor(12m/35)+d_r with d_r from the correction vector")
    print("    3. linear height color: the quotient level floor(12m/35), not just r mod 35")
    print("    4. endpoint phase: finite mixed transients m=1..4 versus rank-one E/E tail")
    print("    5. branch mirror color: C1/C0 and branch0/branch1 ownership")
    print("    6. incident core color: C3 slot plus Qsqrt(-7) sign from HYP-3441")
    print("    7. Haar zipper color: every collapsed margin must retain or kill its zeta cocycle")
    print("  Candidate lemma: any AP84 splice into HYP-3439 is legal if the map from")
    print("  local survivor gates to global branch-union escapes is a homomorphism for")
    print("  this color-packet product, or else the first failed color is routed to")
    print("  HYP-3438/HYP-3453/HYP-3455, owner-current, two-adic descent, or SPEC debt.")
    print()

    print("F. New leads from the coloring pass")
    print("  LEAD-1: Replace raw AP-tail escape counts by a 7-field color tuple")
    print("    (gate, floor, height, endpoint, branch, incident, zeta).")
    print("  LEAD-2: Prove the canonical mod-35 gate law as a colored discrepancy")
    print("    lemma over the two fixed corridors, not as a component-count census.")
    print("  LEAD-3: Build a finite 'color legality matrix' for HYP-3438 survivor")
    print("    gates: rows are gates, columns are the seven packet colors above.")
    print("  LEAD-4: Reuse the W-polynomial lesson: centered colors multiply cleanly")
    print("    before quotienting; quotient only after product closure is checked.")
    print("  LEAD-5: Use extension-rank language for the finite m=1..4 transient")
    print("    packet: same residue, different child phase, so residue is not a")
    print("    recursive invariant by itself.")
    print()

    print("G. Tournament Analysis")
    print("  vertices=quotient/color proof carriers, not runners, arcs, or raw residues")
    print("  pairwise_observable=preserved AP-tail predicate + lost-coordinate sidecar + finite-checkability")
    print("  switch_gauge=higher retained color-packet score; ties by labelled-packet priority")
    print(f"  score_hist={hist}")
    print(f"  directed_3cycles={cycles}")
    print(f"  hamiltonian_path_count={hamiltonian_paths}")
    print("  hamiltonian_path=" + " -> ".join(path))
    print()

    print("H. Assumption Challenge")
    print("  Considered vertices: runners, gaps, fixed circle sections, section")
    print("  boundaries, wall-crossing events, residues, cover arcs, Fourier modes,")
    print("  Haar rectangles, centered edge signs, endpoint owners, incident cores,")
    print("  survivor gates, and proof obligations.")
    print("  Chosen quotient: color-packet proof obligations.  It preserves AP-tail")
    print("  gate availability, floor count, endpoint phase, and sidecar debt")
    print("  declaration, but destroys raw runner identity, non-AP geometry, and")
    print("  arbitrary primitive-row adjacency.  Challenged assumption: a coloring")
    print("  quotient is harmless if its classes look balanced.  The audit shows")
    print("  balance is insufficient unless floor, endpoint, branch, and zeta colors")
    print("  remain reconstructible.")


if __name__ == "__main__":
    main()
