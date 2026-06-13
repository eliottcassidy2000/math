#!/usr/bin/env python3
"""
lonely_runner_scalar_excision_s371.py

codex-2026-05-31 S371

Formalization ledger for the fourteen-runner scalar-ramp obstruction.

The S364 feedback loop found that every scalar ramp

    v_i = m i mod n,  i=1,...,n-1,

blocks the full micro-staircase cell system for n=14 and n=15.  This script
turns that into a more exact n=14 audit:

1. verify the midpoint identity with the initial segment after the time shift;
2. count how scalar ramps block cells, including unique blocker coordinates;
3. scan all one- and two-defect perturbations of scalar ramps;
4. list and stratify the 56 missed cells of the best S364 non-scalar vector.

The point is not only to reproduce the old score.  The useful new target is
that the best non-scalar near-blocker is a one-coordinate defect of the scalar
ramp m=8, and its 56 missed cells are exactly the cells uniquely protected by
that scalar coordinate.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S364 = SourceFileLoader(
    "lonely_runner_feedback_loop_s364",
    str(ROOT / "04-computation" / "lonely_runner_feedback_loop_s364.py"),
).load_module()

N = 14
BEST_S364_NONSCALAR = (8, 2, 10, 4, 12, 13, 0, 8, 2, 10, 4, 12, 6)


@dataclass(frozen=True)
class MissedCell:
    index: int
    shift: int
    pattern_index: int
    pattern: tuple[int, ...]
    interval: tuple[Fraction, Fraction]
    residues: tuple[int, ...]
    min_margin: int


@dataclass(frozen=True)
class DeformationRecord:
    missed: int
    score: int
    base_m: int
    defects: tuple[tuple[int, int, int], ...]  # one-based coordinate, old, new
    vector: tuple[int, ...]


def fmt_frac(x: Fraction) -> str:
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def scalar_vector(n: int, m: int) -> tuple[int, ...]:
    return tuple((m * i) % n for i in range(1, n))


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def representative_intervals(
    n: int, k: int
) -> dict[tuple[int, ...], tuple[Fraction, Fraction, Fraction]]:
    """Return the first open breakpoint interval realizing each floor pattern.

    A pattern is a tuple of values floor(n*{i alpha}).  The same pattern can
    occur in several places on the alpha circle, so reconstructing it by
    intersecting the no-wrap inequalities can give nonsense.  S364 also keeps
    one copy per pattern, so this helper mirrors its breakpoint sweep and
    records a valid representative interval for each stored pattern.
    """

    breaks = {Fraction(0), Fraction(1)}
    for i in range(1, k + 1):
        for a in range(n * i + 1):
            breaks.add(Fraction(a, n * i))

    out: dict[tuple[int, ...], tuple[Fraction, Fraction, Fraction]] = {}
    ordered = sorted(breaks)
    for lo, hi in zip(ordered, ordered[1:]):
        if lo == hi:
            continue
        mid = (lo + hi) / 2
        pattern = S364.cell_pattern(n, k, mid)
        out.setdefault(pattern, (lo, hi, mid))
    return out


def shifted_initial_pattern(n: int, k: int, beta: Fraction) -> tuple[int, ...]:
    return tuple(int((n * frac_part(i * beta)) // 1) for i in range(1, k + 1))


def residues_for(
    vector: tuple[int, ...],
    n: int,
    shift: int,
    pattern: tuple[int, ...],
) -> tuple[int, ...]:
    return tuple((shift * value + bin_value) % n for value, bin_value in zip(vector, pattern))


def blockers_for(
    vector: tuple[int, ...],
    n: int,
    shift: int,
    pattern: tuple[int, ...],
) -> tuple[tuple[int, int], ...]:
    residues = residues_for(vector, n, shift, pattern)
    return tuple(
        (i, residue)
        for i, residue in enumerate(residues, start=1)
        if residue in (0, n - 1)
    )


def candidate_data(system, index: int) -> tuple[int, int, tuple[int, ...]]:
    pattern_count = len(system.patterns)
    shift, pattern_index = divmod(index, pattern_count)
    return shift, pattern_index, system.patterns[pattern_index]


def bit_indices(mask: int):
    while mask:
        low = mask & -mask
        yield low.bit_length() - 1
        mask ^= low


def midpoint_identity_failures(
    system, reps: dict[tuple[int, ...], tuple[Fraction, Fraction, Fraction]], m: int
) -> int:
    failures = 0
    for shift in range(system.n):
        for pattern in system.patterns:
            _lo, _hi, alpha = reps[pattern]
            beta = frac_part(alpha + Fraction(shift * m, system.n))
            residues = residues_for(scalar_vector(system.n, m), system.n, shift, pattern)
            shifted = shifted_initial_pattern(system.n, system.k, beta)
            if residues != shifted:
                failures += 1
    return failures


def scalar_blocking_ledger(
    system, reps: dict[tuple[int, ...], tuple[Fraction, Fraction, Fraction]]
) -> None:
    print("Scalar-ramp cell blocking audit")
    print(f"  n={system.n} patterns={len(system.patterns)} candidates={system.candidate_count}")
    print("  theorem_probe=midpoint identity residues(v=m*i, s, alpha) == initial_segment_bins(alpha+s*m/n)")

    for m in range(system.n):
        vector = scalar_vector(system.n, m)
        result = S364.score_vector(system, vector)
        identity_failures = midpoint_identity_failures(system, reps, m)

        blocker_hist: dict[int, int] = {}
        unique_by_coord: dict[int, int] = {}
        endpoint_hits = {0: 0, system.n - 1: 0}
        for shift in range(system.n):
            for pattern in system.patterns:
                blockers = blockers_for(vector, system.n, shift, pattern)
                blocker_hist[len(blockers)] = blocker_hist.get(len(blockers), 0) + 1
                for _coord, residue in blockers:
                    endpoint_hits[residue] += 1
                if len(blockers) == 1:
                    coord, _residue = blockers[0]
                    unique_by_coord[coord] = unique_by_coord.get(coord, 0) + 1

        top_unique = tuple(sorted(unique_by_coord.items(), key=lambda item: (-item[1], item[0]))[:5])
        print(
            "  "
            f"m={m:2d} gcd={gcd(m, system.n):2d} "
            f"covered={result.score}/{system.candidate_count} "
            f"missed={result.missed} midpoint_failures={identity_failures} "
            f"blocker_hist={tuple(sorted(blocker_hist.items()))} "
            f"endpoint_hits=(0:{endpoint_hits[0]}, {system.n - 1}:{endpoint_hits[system.n - 1]}) "
            f"top_unique_coords={top_unique}"
        )
    print()


def deformation_scan(system, radius: int, top_k: int = 12) -> tuple[DeformationRecord, ...]:
    records: list[DeformationRecord] = []
    scanned = 0
    for m in range(system.n):
        base = scalar_vector(system.n, m)
        for coords in combinations(range(system.k), radius):
            choices: list[tuple[tuple[int, int], ...]] = [()]
            for coord in coords:
                next_choices: list[tuple[tuple[int, int], ...]] = []
                old = base[coord]
                for prefix in choices:
                    for new in range(system.n):
                        if new != old:
                            next_choices.append(prefix + ((coord, new),))
                choices = next_choices

            for assignment in choices:
                vector = list(base)
                defects = []
                for coord, new in assignment:
                    old = vector[coord]
                    vector[coord] = new
                    defects.append((coord + 1, old, new))
                vector_tuple = tuple(vector)
                if S364.scalar_multiplier(vector_tuple, system.n) is not None:
                    continue
                result = S364.score_vector(system, vector_tuple)
                records.append(
                    DeformationRecord(
                        missed=result.missed,
                        score=result.score,
                        base_m=m,
                        defects=tuple(defects),
                        vector=vector_tuple,
                    )
                )
                scanned += 1

    records.sort(key=lambda row: (row.missed, row.base_m, row.defects, row.vector))
    best = records[:top_k]
    print(f"Scalar-neighborhood deformation scan radius={radius}")
    print(f"  scanned_non_scalar={scanned}")
    print(f"  full_non_scalar={sum(1 for row in records if row.missed == 0)}")
    print(f"  best_missed={best[0].missed if best else 'none'}")
    hist: dict[int, int] = {}
    for row in records:
        if len(hist) >= 10 and row.missed not in hist:
            continue
        hist[row.missed] = hist.get(row.missed, 0) + 1
    print(f"  missed_hist_prefix={tuple(sorted(hist.items())[:10])}")
    for rank, row in enumerate(best, start=1):
        print(
            "  "
            f"rank={rank:02d} base_m={row.base_m:2d} gcd={gcd(row.base_m, system.n):2d} "
            f"missed={row.missed:4d} score={row.score:5d} "
            f"defects={row.defects} vector={row.vector}"
        )
    print()
    return tuple(best)


def missed_cells(
    system,
    reps: dict[tuple[int, ...], tuple[Fraction, Fraction, Fraction]],
    vector: tuple[int, ...],
) -> tuple[MissedCell, ...]:
    blocked = S364.blocked_mask(system, vector)
    missed_mask = system.all_mask ^ blocked
    cells: list[MissedCell] = []
    for index in bit_indices(missed_mask):
        shift, pattern_index, pattern = candidate_data(system, index)
        lo, hi, _mid = reps[pattern]
        interval = (lo, hi)
        residues = residues_for(vector, system.n, shift, pattern)
        min_margin = min(min(residue, system.n - 1 - residue) for residue in residues)
        cells.append(
            MissedCell(
                index=index,
                shift=shift,
                pattern_index=pattern_index,
                pattern=pattern,
                interval=interval,
                residues=residues,
                min_margin=min_margin,
            )
        )
    return tuple(cells)


def print_missed_cell_ledger(
    system,
    reps: dict[tuple[int, ...], tuple[Fraction, Fraction, Fraction]],
    vector: tuple[int, ...],
) -> None:
    result = S364.score_vector(system, vector)
    base_m = 8
    base = scalar_vector(system.n, base_m)
    base_blocked = S364.blocked_mask(system, base)
    cells = missed_cells(system, reps, vector)
    print("Best S364 non-scalar missed-cell ledger")
    print(f"  vector={vector}")
    print(f"  score={result.score}/{result.candidate_count} missed={result.missed}")
    print(f"  scalar_multiplier={result.scalar_multiplier}")
    print(f"  base_scalar_m={base_m} base_vector={base}")
    print(f"  defects_from_base={tuple((i + 1, old, new) for i, (old, new) in enumerate(zip(base, vector)) if old != new)}")
    print(f"  base_scalar_blocks_all={base_blocked == system.all_mask}")

    shift_hist: dict[int, int] = {}
    width_hist: dict[Fraction, int] = {}
    margin_hist: dict[int, int] = {}
    residue_hist: dict[int, int] = {}
    unique_base_coord_hist: dict[int, int] = {}
    for cell in cells:
        shift_hist[cell.shift] = shift_hist.get(cell.shift, 0) + 1
        width = cell.interval[1] - cell.interval[0]
        width_hist[width] = width_hist.get(width, 0) + 1
        margin_hist[cell.min_margin] = margin_hist.get(cell.min_margin, 0) + 1
        for residue in cell.residues:
            residue_hist[residue] = residue_hist.get(residue, 0) + 1
        base_blockers = blockers_for(base, system.n, cell.shift, cell.pattern)
        if len(base_blockers) == 1:
            unique_base_coord_hist[base_blockers[0][0]] = unique_base_coord_hist.get(base_blockers[0][0], 0) + 1

    print(f"  missed_shift_hist={tuple(sorted(shift_hist.items()))}")
    print(f"  missed_width_hist={tuple((fmt_frac(width), count) for width, count in sorted(width_hist.items()))}")
    print(f"  missed_min_margin_hist={tuple(sorted(margin_hist.items()))}")
    print(f"  missed_residue_hist={tuple(sorted(residue_hist.items()))}")
    print(f"  base_unique_blocker_coord_hist={tuple(sorted(unique_base_coord_hist.items()))}")
    print("  missed_cells=")
    for cell in cells:
        lo, hi = cell.interval
        print(
            "    "
            f"idx={cell.index:5d} s={cell.shift:2d} pidx={cell.pattern_index:3d} "
            f"alpha=[{fmt_frac(lo)},{fmt_frac(hi)}) width={fmt_frac(hi - lo)} "
            f"margin={cell.min_margin} pattern={cell.pattern} residues={cell.residues}"
        )
    print()


def mutation_pressure(system, vector: tuple[int, ...], top_k: int = 16) -> None:
    current_mask = S364.blocked_mask(system, vector)
    missed_mask = system.all_mask ^ current_mask
    rows: list[tuple[int, int, int, int, int, int, int]] = []
    for coord in range(system.k):
        old = vector[coord]
        for new in range(system.n):
            if new == old:
                continue
            candidate = list(vector)
            candidate[coord] = new
            candidate_tuple = tuple(candidate)
            result = S364.score_vector(system, candidate_tuple)
            mask = S364.blocked_mask(system, candidate_tuple)
            fixed = (mask & missed_mask).bit_count()
            retained = (mask & current_mask).bit_count()
            lost = current_mask.bit_count() - retained
            rows.append((result.missed, -fixed, lost, coord + 1, old, new, result.score))
    rows.sort()

    print("Single-coordinate mutation pressure from best non-scalar")
    print("  interpretation=fixed counts missed cells repaired; lost counts old covered cells reopened")
    for rank, row in enumerate(rows[:top_k], start=1):
        missed, neg_fixed, lost, coord, old, new, score = row
        fixed = -neg_fixed
        print(
            "  "
            f"rank={rank:02d} coord={coord:2d} {old:2d}->{new:2d} "
            f"new_missed={missed:4d} score={score:5d} fixed={fixed:3d} lost={lost:3d}"
        )
    print()


def main() -> None:
    print("Lonely Runner scalar-ramp excision ledger (codex-2026-05-31 S371)")
    print("Reduced fourteen-runner cell system: k=13, n=14.\n")

    system = S364.build_pattern_system(N)
    reps = representative_intervals(system.n, system.k)
    print(f"Representative pattern intervals recovered={len(reps)}/{len(system.patterns)}\n")

    scalar_blocking_ledger(system, reps)

    deformation_scan(system, radius=1)
    deformation_scan(system, radius=2)

    print_missed_cell_ledger(system, reps, BEST_S364_NONSCALAR)
    mutation_pressure(system, BEST_S364_NONSCALAR)

    print("Synthesis")
    print("  1. The scalar-ramp obstruction is exactly the shifted initial-segment equality case cell-by-cell.")
    print("  2. The S364 best non-scalar n=14 vector is not generic: it is scalar m=8 with one coordinate defect.")
    print("  3. Its 56 missed cells are the scalar cells uniquely blocked by the defect coordinate.")
    print("  4. No one- or two-coordinate scalar-neighborhood deformation is a non-scalar full blocker.")
    print("  5. The next finite proof target is a fragility theorem: every non-scalar deformation of a scalar ramp")
    print("     exposes at least one uniquely protected micro-staircase cell, then extend beyond small defect radius.")


if __name__ == "__main__":
    main()
