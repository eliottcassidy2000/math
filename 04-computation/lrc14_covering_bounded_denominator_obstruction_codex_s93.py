#!/usr/bin/env python3
"""
codex-2026-06-22 S93

Bounded-denominator witnesses for LRC(14), stress-tested inside the THM-523
covering hard core.

Main point:
  The tempting statement

      every covering 13-set has a lonely witness a/D with D <= B0

  is false for every fixed B0.  The covering tower

      S_m = {1,2,...,11,13,84m}

  is primitive and covering.  If m is divisible by lcm(1,...,B), then 84m is
  divisible by every D <= B, so at every rational a/D with D <= B the last
  runner is exactly at the origin.  Thus no D <= B witness exists.

The script also records why the idea looked plausible in random scouts:
non-adversarial covering sets and the first thousands of the hard tower still
have small witnesses.  The obstruction is divisor loading, not typical size.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import gcd, lcm
from random import Random


N = 14
BASE_TOWER = tuple(list(range(1, 12)) + [13])


def gcd_all(values: tuple[int, ...] | list[int]) -> int:
    g = 0
    for value in values:
        g = gcd(g, value)
    return g


def is_covering(speeds: tuple[int, ...] | list[int]) -> bool:
    return all(any(v % q == 0 for v in speeds) for q in range(2, N + 1))


def margin_num_at(speeds: tuple[int, ...] | list[int], a: int, d: int) -> int:
    """Return min_v (14*dist_D(a*v)-D).  Nonnegative means level-1/14 safe."""

    best: int | None = None
    for v in speeds:
        r = (a * (v % d)) % d
        dist = min(r, d - r)
        value = N * dist - d
        if best is None or value < best:
            best = value
    assert best is not None
    return best


@dataclass(frozen=True)
class Witness:
    a: int
    d: int
    margin_num: int

    @property
    def margin_den(self) -> int:
        return N * self.d

    def text(self) -> str:
        return f"{self.a}/{self.d} margin={self.margin_num}/{self.margin_den}"


def least_witness(
    speeds: tuple[int, ...] | list[int],
    bound: int,
    units_by_d: dict[int, tuple[int, ...]],
) -> Witness | None:
    for d in range(2, bound + 1):
        for a in units_by_d[d]:
            margin_num = margin_num_at(speeds, a, d)
            if margin_num >= 0:
                return Witness(a, d, margin_num)
    return None


def lcm_to(bound: int) -> int:
    value = 1
    for d in range(1, bound + 1):
        value = lcm(value, d)
    return value


def tower_set(m: int) -> tuple[int, ...]:
    return tuple(sorted(BASE_TOWER + (84 * m,)))


def no_small_denominator_certificate(bound: int) -> tuple[int, tuple[int, ...]]:
    """Return m and S_m proving no denominator <= bound can work."""

    m = lcm_to(bound)
    speeds = tower_set(m)
    assert gcd_all(speeds) == 1
    assert is_covering(speeds)
    for d in range(2, bound + 1):
        assert (84 * m) % d == 0
    return m, speeds


def candidate_safe_for_base(a: int, d: int) -> bool:
    return margin_num_at(BASE_TOWER, a, d) >= 0


def tower_scan(limit_m: int, bound: int, units_by_d: dict[int, tuple[int, ...]]) -> None:
    # Precompute base-safe candidates.  For S_m only the tail 84m changes.
    candidates: list[tuple[int, int]] = []
    for d in range(2, bound + 1):
        for a in units_by_d[d]:
            if candidate_safe_for_base(a, d):
                candidates.append((a, d))

    max_witness = Witness(0, 0, 0)
    first_above_41: tuple[int, Witness] | None = None
    first_by_d: dict[int, tuple[int, Witness]] = {}
    hist: dict[int, int] = {}
    failures = 0

    for m in range(1, limit_m + 1):
        found: Witness | None = None
        speeds = tower_set(m)
        for a, d in candidates:
            margin_num = margin_num_at(speeds, a, d)
            if margin_num >= 0:
                found = Witness(a, d, margin_num)
                break
        if found is None:
            failures += 1
            continue
        hist[found.d] = hist.get(found.d, 0) + 1
        first_by_d.setdefault(found.d, (m, found))
        if found.d > max_witness.d:
            max_witness = found
        if found.d > 41 and first_above_41 is None:
            first_above_41 = (m, found)

    print("TOWER SCOUT: S_m={1..11,13,84m}")
    print(f"  scanned m<= {limit_m}, denominator search bound={bound}")
    print(f"  base-safe candidate fractions <= {bound}: {len(candidates)}")
    print(f"  failures within search bound: {failures}")
    print(f"  max least denominator observed: {max_witness.d}")
    if first_above_41 is not None:
        m, witness = first_above_41
        print(f"  first row needing D>41: m={m}, witness={witness.text()}")
    high_hist = [(d, c) for d, c in sorted(hist.items()) if d >= 35]
    print(f"  high-denominator histogram (D>=35): {high_hist}")
    for d in (41, 53, 55, 65, 67):
        if d in first_by_d:
            m, witness = first_by_d[d]
            print(f"  first D={d}: m={m}, witness={witness.text()}")


def ap_repair_scan(limit_m: int, bound: int, units_by_d: dict[int, tuple[int, ...]]) -> None:
    max_row: tuple[int, int, Witness] | None = None
    failures: list[tuple[int, int]] = []
    count = 0
    hist: dict[int, int] = {}

    for drop in range(1, 14):
        base = set(range(1, 14))
        base.remove(drop)
        for m in range(1, limit_m + 1):
            speeds = tuple(sorted(base | {14 * m}))
            if len(speeds) != 13 or not is_covering(speeds) or gcd_all(speeds) != 1:
                continue
            count += 1
            witness = least_witness(speeds, bound, units_by_d)
            if witness is None:
                failures.append((drop, m))
                continue
            hist[witness.d] = hist.get(witness.d, 0) + 1
            if max_row is None or witness.d > max_row[2].d:
                max_row = (drop, m, witness)

    print("AP-REPAIR SCOUT: ({1..13}\\{drop}) union {14m}")
    print(f"  covering primitive rows scanned: {count}")
    print(f"  failures within D<={bound}: {len(failures)}")
    if max_row is not None:
        drop, m, witness = max_row
        print(f"  max least denominator: drop={drop}, m={m}, witness={witness.text()}")
    print(f"  high-denominator histogram (D>=35): {[(d, c) for d, c in sorted(hist.items()) if d >= 35]}")


def random_covering_set(rng: Random, max_speed: int) -> tuple[int, ...]:
    for _ in range(2000):
        speeds: set[int] = set()
        if rng.random() < 0.70:
            speeds.add(1)

        while len(speeds) < 13 and not is_covering(tuple(speeds)):
            uncovered = [
                q for q in range(2, N + 1) if not any(v % q == 0 for v in speeds)
            ]
            rng.shuffle(uncovered)
            modulus = 1
            took = False
            for q in uncovered:
                if (not took) or rng.random() < 0.35:
                    new_modulus = lcm(modulus, q)
                    if new_modulus <= max_speed:
                        modulus = new_modulus
                        took = True
            if modulus <= max_speed:
                speeds.add(modulus * rng.randint(1, max_speed // modulus))

        while len(speeds) < 13:
            speeds.add(rng.randint(1, max_speed))

        row = tuple(sorted(speeds))
        if len(row) == 13 and gcd_all(row) == 1 and is_covering(row):
            return row

    raise RuntimeError("failed to generate a covering primitive row")


def random_obligation_scan(
    samples: int,
    max_speed: int,
    bound: int,
    units_by_d: dict[int, tuple[int, ...]],
) -> None:
    rng = Random(2865)
    hist: dict[int, int] = {}
    max_row: tuple[tuple[int, ...], Witness] | None = None
    failures: list[tuple[int, ...]] = []

    for _ in range(samples):
        speeds = random_covering_set(rng, max_speed)
        witness = least_witness(speeds, bound, units_by_d)
        if witness is None:
            failures.append(speeds)
            continue
        hist[witness.d] = hist.get(witness.d, 0) + 1
        if max_row is None or witness.d > max_row[1].d:
            max_row = (speeds, witness)

    print("RANDOM OBLIGATION SCOUT")
    print(f"  samples={samples}, max_speed={max_speed}, search D<={bound}")
    print(f"  failures within search bound: {len(failures)}")
    if max_row is not None:
        speeds, witness = max_row
        print(f"  max least denominator: {witness.text()}")
        print(f"  row: {speeds}")
    print(f"  denominator histogram (D>=20): {[(d, c) for d, c in sorted(hist.items()) if d >= 20]}")


def tournament_analysis() -> None:
    vertices = [
        ("covering_denominator_no_go", (6, 6, 6, 2)),
        ("scaled_denominator_or_speed_bound", (5, 5, 4, 3)),
        ("THM523_covering_reduction", (5, 4, 5, 5)),
        ("THM524_binding_pair_switches", (4, 4, 5, 4)),
        ("HYP2864_sheet_gcd_quotient", (4, 3, 4, 3)),
        ("THM565_three_gap_floor", (3, 3, 4, 4)),
        ("fixed_B_residue_atlas", (2, 6, 2, 5)),
        ("raw_runner_vertices", (0, 0, 1, 1)),
    ]
    ordered = sorted(vertices, key=lambda item: item[1], reverse=True)
    scores = {name: 0 for name, _ in vertices}
    names = [name for name, _ in vertices]
    keys = dict(vertices)
    for i, left in enumerate(names):
        for right in names[i + 1 :]:
            if keys[left] > keys[right]:
                scores[left] += 1
            else:
                scores[right] += 1

    print("TOURNAMENT ANALYSIS")
    print("  vertices: proof carriers, not runners")
    print("  observable: (rules out false route, preserves LRC witness predicate, finite-check value, formalization readiness)")
    print("  switch: lexicographically larger observable wins; ties use the displayed order")
    print(f"  Hamiltonian path: {' > '.join(name for name, _ in ordered)}")
    print(f"  score histogram: {sorted(scores.values())}")
    print("  directed 3-cycles: 0 (the chosen proof-carrier quotient is transitive)")
    print("  SCCs: singleton; Hamiltonian-path count: 1")
    print("  challenged assumption: covering-set residuals admit a fixed denominator cap")


def main() -> None:
    units_by_d = {
        d: tuple(a for a in range(1, d) if gcd(a, d) == 1)
        for d in range(2, 201)
    }

    print("=" * 78)
    print("A. PROVED NO-GO: covering sets have no uniform bounded denominator")
    print("=" * 78)
    for bound in (14, 26, 41, 67, 80):
        m, speeds = no_small_denominator_certificate(bound)
        witness = least_witness(speeds, 200, units_by_d)
        print(
            f"  B={bound:2d}: m=lcm(1..B) has {len(str(m))} digits; "
            f"covering={is_covering(speeds)} primitive={gcd_all(speeds)==1}; "
            f"no D<=B by divisibility; first witness <=200: "
            f"{witness.text() if witness else 'none'}"
        )

    print("\n" + "=" * 78)
    print("B. WHY THE FALSE CONJECTURE LOOKED PLAUSIBLE")
    print("=" * 78)
    named = [
        ("champion m=1", tower_set(1), 80),
        ("first D>41 tower row m=6", tower_set(6), 100),
        ("closer q-floor row {1..12,182}", tuple(list(range(1, 13)) + [182]), 80),
        ("drop-6 plus 84", (1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 84), 80),
    ]
    for label, speeds, bound in named:
        witness = least_witness(speeds, bound, units_by_d)
        print(f"  {label:32s}: {witness.text() if witness else 'none'}")

    print()
    tower_scan(limit_m=5000, bound=100, units_by_d=units_by_d)
    print()
    ap_repair_scan(limit_m=500, bound=100, units_by_d=units_by_d)
    print()
    random_obligation_scan(samples=300, max_speed=10**6, bound=120, units_by_d=units_by_d)

    print("\n" + "=" * 78)
    print("C. PROOF-CARRIER TOURNAMENT")
    print("=" * 78)
    tournament_analysis()


if __name__ == "__main__":
    main()
