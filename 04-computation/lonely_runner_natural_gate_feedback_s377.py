#!/usr/bin/env python3
"""Product-sum natural gates in the 14/15-runner feedback loop.

codex-2026-05-31 S377

This is a follow-up to the S364/S367/S371 lonely-runner work and the S365
natural-operation graph work.  The additive natural-number shadow is only the
transitive order, while multiplication keeps the divisor DAG.  Product-sum
equations mark the critical pairs where the two operation systems hit the same
mode.

Here those critical modes are used as a disciplined source of new LRC probes:

* analyze which natural-number modes are the worst quotient one-defects in the
  n=14 and n=15 micro-staircase systems;
* use product-sum target modes as speed replacements in exact LRC interval
  scans;
* try operation-critical near-counterexample families before returning to the
  proof route.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import gcd, prod
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()
S362 = SourceFileLoader(
    "lonely_runner_bohr_descent_s362",
    str(ROOT / "04-computation" / "lonely_runner_bohr_descent_s362.py"),
).load_module()
S367 = SourceFileLoader(
    "lonely_runner_k13_scalar_gauge_s367",
    str(ROOT / "04-computation" / "lonely_runner_k13_scalar_gauge_s367.py"),
).load_module()


@dataclass(frozen=True)
class CoreSeed:
    arity: int
    target: int
    core: tuple[int, ...]
    ones: int


@dataclass(frozen=True)
class DefectProbe:
    n: int
    coord: int
    residue: int
    additive_parents: tuple[tuple[int, int], ...]
    multiplicative_parents: tuple[tuple[int, int], ...]
    is_collision_target: bool
    target_arities: tuple[int, ...]
    missed: int
    covered: int
    shift_hist: tuple[tuple[int, int], ...]
    width_hist: tuple[tuple[Fraction, int], ...]


@dataclass(frozen=True)
class SpeedProbe:
    label: str
    speeds: tuple[int, ...]
    classification: str
    gap_ratio: Fraction
    max_gap: Fraction
    boundary_witnesses: int
    unprotected: int
    peel_depth: int
    core_endpoints: int
    collision_targets: tuple[int, ...]


def fmt_frac(x: Fraction | None) -> str:
    return S356.fmt_frac(x)


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def collision_seeds(max_arity: int, max_target: int) -> tuple[CoreSeed, ...]:
    """Enumerate nonunit product-sum cores by the defect normal form."""

    out: list[CoreSeed] = []

    def rec(start: int, core: list[int], p: int, s: int) -> None:
        if len(core) >= 2:
            ones = p - s
            arity = len(core) + ones
            if ones >= 0 and 2 <= arity <= max_arity and p <= max_target:
                out.append(CoreSeed(arity=arity, target=p, core=tuple(core), ones=ones))

        max_next = max_target // p
        for value in range(start, max_next + 1):
            rec(value, core + [value], p * value, s + value)

    rec(2, [], 1, 0)
    out.sort(key=lambda seed: (seed.arity, seed.target, len(seed.core), seed.core))
    return tuple(out)


def additive_parents(z: int) -> tuple[tuple[int, int], ...]:
    return tuple((x, z - x) for x in range(1, z // 2 + 1))


def multiplicative_parents(z: int) -> tuple[tuple[int, int], ...]:
    return tuple((d, z // d) for d in divisors(z) if d <= z // d and d > 1 and z % d == 0)


def collision_target_map(max_arity: int, max_target: int) -> dict[int, tuple[CoreSeed, ...]]:
    by_target: dict[int, list[CoreSeed]] = defaultdict(list)
    for seed in collision_seeds(max_arity=max_arity, max_target=max_target):
        by_target[seed.target].append(seed)
    return {target: tuple(seeds) for target, seeds in by_target.items()}


def missed_candidate_indices(system, result) -> list[int]:
    mask = S367.blocked_mask(system, result.vector)
    return [idx for idx in range(system.candidate_count) if not ((mask >> idx) & 1)]


def probe_one_defect(system, coord: int, residue: int, target_map: dict[int, tuple[CoreSeed, ...]]) -> DefectProbe:
    vector = [0] * system.k
    vector[coord - 1] = residue
    result = S367.score_vector(system, tuple(vector))
    missed = missed_candidate_indices(system, result)

    shift_counter: Counter[int] = Counter()
    width_counter: Counter[Fraction] = Counter()
    for idx in missed:
        shift, pattern_index = system.candidate_meta[idx]
        pattern = system.patterns[pattern_index]
        shift_counter[shift] += 1
        width_counter[pattern.hi - pattern.lo] += 1

    return DefectProbe(
        n=system.n,
        coord=coord,
        residue=residue,
        additive_parents=additive_parents(coord),
        multiplicative_parents=multiplicative_parents(coord),
        is_collision_target=coord in target_map,
        target_arities=tuple(seed.arity for seed in target_map.get(coord, ())),
        missed=result.missed,
        covered=result.score,
        shift_hist=tuple(sorted(shift_counter.items())),
        width_hist=tuple(sorted(width_counter.items(), key=lambda item: (item[0], item[1]))),
    )


def one_defect_scan(n: int, target_map: dict[int, tuple[CoreSeed, ...]]) -> tuple[DefectProbe, ...]:
    system = S367.build_pattern_system(n)
    rows: list[DefectProbe] = []
    for coord in range(2, system.k + 1):
        for residue in range(1, n):
            rows.append(probe_one_defect(system, coord, residue, target_map))
    rows.sort(key=lambda row: (row.missed, row.coord, row.residue))
    return tuple(rows)


def classify(row) -> str:
    if row.forbidden_length < 1:
        return "positive_gap"
    if row.boundary_witness_count:
        return "boundary_only"
    return "open_cover"


def speed_probe(label: str, speeds: tuple[int, ...], targets: set[int]) -> SpeedProbe:
    row = S356.report(label, list(speeds))
    descent = S362.summarize(list(row.speeds))
    ratio = row.max_gap / row.threshold if row.threshold else Fraction(0)
    return SpeedProbe(
        label=label,
        speeds=tuple(row.speeds),
        classification=classify(row),
        gap_ratio=ratio,
        max_gap=row.max_gap,
        boundary_witnesses=row.boundary_witness_count,
        unprotected=descent.unprotected_count,
        peel_depth=len(descent.peel_layers),
        core_endpoints=descent.core_endpoint_count,
        collision_targets=tuple(v for v in row.speeds if v in targets),
    )


def initial_replacement_sets(k: int, targets: tuple[int, ...], max_insertions: int) -> tuple[tuple[str, tuple[int, ...]], ...]:
    """Replace the largest initial speeds by operation-critical targets."""

    base = tuple(range(1, k + 1))
    out: dict[tuple[int, ...], str] = {}
    usable_targets = tuple(sorted(t for t in targets if t > k))
    for r in range(1, max_insertions + 1):
        for inserted in combinations(usable_targets, r):
            if len(set(inserted)) < len(inserted):
                continue
            removed = base[-r:]
            speeds = tuple(sorted(set(base[:-r]) | set(inserted)))
            if len(speeds) != k:
                continue
            if gcd_many(speeds) != 1:
                continue
            label = f"replace {removed} by {inserted}"
            out.setdefault(speeds, label)
    return tuple((label, speeds) for speeds, label in sorted(out.items(), key=lambda item: item[0]))


def gcd_many(values: tuple[int, ...]) -> int:
    g = 0
    for value in values:
        g = gcd(g, value)
    return g


def operation_cycle_sets(k: int, n: int, target_map: dict[int, tuple[CoreSeed, ...]]) -> tuple[tuple[str, tuple[int, ...]], ...]:
    """Generate small product-sum/divisor-gate speed-set probes."""

    target_modes = tuple(sorted(t for t in target_map if n <= t <= 2 * n))
    candidates: list[tuple[str, tuple[int, ...]]] = list(initial_replacement_sets(k, target_modes, 3))

    # Add a few deliberately engineered operation-critical overloads.  These are
    # disproof-shaped: keep many small speeds, then add target modes and their
    # multiplicative parents.
    for target in target_modes:
        parents = tuple(x for pair in multiplicative_parents(target) for x in pair)
        payload = tuple(sorted(set((target,) + parents + tuple(d for d in divisors(target) if d > 1))))
        fillers = tuple(v for v in range(1, 3 * n + 1) if v not in payload)
        speeds = tuple(sorted((payload + fillers)[:k]))
        if len(speeds) == k and gcd_many(speeds) == 1:
            candidates.append((f"target {target} with divisor payload {payload}", speeds))

    # A closed critical-pair ladder: all target modes in [n,2n] plus the least
    # natural fillers.  This is the natural-operation analogue of gate overload.
    payload = target_modes
    fillers = tuple(v for v in range(1, 3 * n + 1) if v not in payload)
    speeds = tuple(sorted((payload + fillers)[:k]))
    if len(speeds) == k and gcd_many(speeds) == 1:
        candidates.append((f"all targets {payload} plus least fillers", speeds))

    dedup: dict[tuple[int, ...], str] = {}
    for label, speeds in candidates:
        dedup.setdefault(speeds, label)
    return tuple((label, speeds) for speeds, label in sorted(dedup.items(), key=lambda item: item[0]))


def print_target_map(target_map: dict[int, tuple[CoreSeed, ...]], n: int) -> None:
    print(f"Product-sum target modes near n={n}")
    for target in sorted(t for t in target_map if 2 <= t <= 2 * n):
        seeds = target_map[target]
        seed_text = "; ".join(
            f"k={seed.arity} core={seed.core} ones={seed.ones}" for seed in seeds
        )
        marker = " <frontier-window>" if n <= target <= 2 * n else ""
        print(f"  z={target:2d}{marker}: {seed_text}")
    print()


def summarize_defect_scan(n: int, probes: tuple[DefectProbe, ...], target_map: dict[int, tuple[CoreSeed, ...]]) -> None:
    top = probes[:12]
    target_rows = [row for row in probes if row.is_collision_target]
    non_target_rows = [row for row in probes if not row.is_collision_target]
    best_by_coord: dict[int, DefectProbe] = {}
    for row in probes:
        best_by_coord.setdefault(row.coord, row)

    print(f"Cycle: quotient one-defect scan for n={n}")
    print(f"  scanned={len(probes)} candidates={probes[0].n * len(S367.build_pattern_system(n).patterns)}")
    print(f"  collision_target_coords<=k={tuple(sorted(c for c in target_map if 2 <= c <= n - 1))}")
    print(
        "  "
        f"best_collision_target_missed={min(row.missed for row in target_rows)} "
        f"best_non_target_missed={min(row.missed for row in non_target_rows)}"
    )
    print("  top one-defects:")
    for row in top:
        tag = "collision" if row.is_collision_target else "plain"
        order = row.n // gcd(row.n, row.residue)
        width_summary = ", ".join(f"{fmt_frac(w)}:{c}" for w, c in row.width_hist[:3])
        print(
            "    "
            f"coord={row.coord:2d} residue={row.residue:2d} order={order:2d} "
            f"missed={row.missed:4d} tag={tag} arities={row.target_arities} "
            f"add_parents={row.additive_parents[:3]} mul_parents={row.multiplicative_parents} "
            f"shift_hist={row.shift_hist[:5]} widths={width_summary}"
        )
    print("  best by coordinate:")
    for coord, row in sorted(best_by_coord.items()):
        tag = "*" if row.is_collision_target else " "
        print(f"    {tag} coord={coord:2d} best_residue={row.residue:2d} missed={row.missed:4d}")
    print()


def summarize_speed_probes(k: int, n: int, probes: tuple[SpeedProbe, ...]) -> None:
    by_class = Counter(row.classification for row in probes)
    best = sorted(probes, key=lambda row: (row.classification != "open_cover", row.gap_ratio, row.speeds))[:10]
    print(f"Cycle: operation-critical speed probes for k={k}, n={n}")
    print(f"  tested={len(probes)} class_hist={tuple(sorted(by_class.items()))}")
    print("  tightest candidates:")
    for row in best:
        print(
            "    "
            f"{row.label}: class={row.classification} "
            f"gap/thresh={float(row.gap_ratio):.6f} max_gap={fmt_frac(row.max_gap)} "
            f"boundary_witnesses={row.boundary_witnesses} "
            f"unprotected={row.unprotected} peel={row.peel_depth} core_E={row.core_endpoints} "
            f"targets={row.collision_targets} speeds={row.speeds}"
        )
    print()


def print_disproof_synthesis(k: int, n: int, probes: tuple[SpeedProbe, ...]) -> None:
    open_covers = [row for row in probes if row.classification == "open_cover"]
    boundary = [row for row in probes if row.classification == "boundary_only"]
    positive = [row for row in probes if row.classification == "positive_gap"]
    print(f"Cycle: disproof pressure from natural-operation candidates k={k}, n={n}")
    print(f"  open_covers={len(open_covers)} boundary_only={len(boundary)} positive_gap={len(positive)}")
    if positive:
        by_peel = sorted(positive, key=lambda row: (-row.peel_depth, row.gap_ratio, row.speeds))[:5]
        print("  deepest leaking systems:")
        for row in by_peel:
            print(
                "    "
                f"peel={row.peel_depth:2d} gap/thresh={float(row.gap_ratio):.6f} "
                f"unprotected={row.unprotected:3d} targets={row.collision_targets} "
                f"speeds={row.speeds}"
            )
    print()


def main() -> None:
    print("Lonely Runner natural-gate feedback loop (codex-2026-05-31 S377)")
    print("All LRC interval/protection data use exact rational arithmetic.\n")

    target_map = collision_target_map(max_arity=24, max_target=60)
    print_target_map(target_map, 14)

    probes14 = one_defect_scan(14, target_map)
    summarize_defect_scan(14, probes14, target_map)

    targets14 = set(t for t in target_map if 14 <= t <= 28)
    speed_sets14 = operation_cycle_sets(13, 14, target_map)
    speed_probes14 = tuple(speed_probe(label, speeds, targets14) for label, speeds in speed_sets14)
    summarize_speed_probes(13, 14, speed_probes14)
    print_disproof_synthesis(13, 14, speed_probes14)

    print_target_map(target_map, 15)
    probes15 = one_defect_scan(15, target_map)
    summarize_defect_scan(15, probes15, target_map)

    targets15 = set(t for t in target_map if 15 <= t <= 30)
    speed_sets15 = operation_cycle_sets(14, 15, target_map)
    speed_probes15 = tuple(speed_probe(label, speeds, targets15) for label, speeds in speed_sets15)
    summarize_speed_probes(14, 15, speed_probes15)
    print_disproof_synthesis(14, 15, speed_probes15)

    print("Synthesis")
    print("  14-route: the worst quotient one-defect is coordinate 6, the first")
    print("  visible distinct product-sum resonance 1+2+3=6.  Coordinate 12, the")
    print("  next rich product-sum target, is second.  This turns the S371 56-cell")
    print("  target into a natural-mode fragility problem, not just a residue accident.")
    print("  15-route: the best one-defect splits into two tied order-3 stencils,")
    print("  coordinate 6 and endpoint coordinate 14.  This matches the 15=3*5 CRT")
    print("  split and suggests proving the order-3 stencil pair before attempting")
    print("  a full 3x5 micro-staircase.")
    print("  disproof-route: product-sum/divisor critical speed sets still all leak.")
    print("  The natural-operation candidates produce positive gaps and empty cores,")
    print("  so the next counterexample search should solve protection cycles directly")
    print("  instead of overloading operation-critical target modes.")


if __name__ == "__main__":
    main()
