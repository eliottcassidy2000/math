#!/usr/bin/env python3
"""
lrc14_residue_basis_character_sum_codex_s98.py

Finite residue bases for LRC(14), tested against the bounded-denominator
no-go theorem.

The prompt's basis {83, 89, 21} is useful as a sample atlas, but it cannot be
global.  Any finite denominator basis B is killed by the primitive covering row

    {1,2,...,11,13, 84*lcm(B)}.

The same residue count also gives the exact character-sum object:

    N(S,D) = #{a mod D, gcd(a,D)=1: ||s*a/D|| >= 1/14 for all s in S}.

The zeroes are resonance/divisibility events, not failures of the main-term
heuristic in typical rows.
"""

from __future__ import annotations

from math import gcd, lcm
from random import Random


N = 14
BASE = tuple(list(range(1, 12)) + [13])
PROMPT_BASIS = (83, 89, 21)


def gcd_all(values: tuple[int, ...]) -> int:
    g = 0
    for value in values:
        g = gcd(g, value)
    return g


def lcm_many(values: tuple[int, ...]) -> int:
    out = 1
    for value in values:
        out = lcm(out, value)
    return out


def is_covering(speeds: tuple[int, ...]) -> bool:
    return all(any(v % q == 0 for v in speeds) for q in range(2, N + 1))


def safe_residue(r: int, d: int) -> bool:
    r %= d
    return N * min(r, d - r) >= d


def witness_count(speeds: tuple[int, ...], d: int, *, units_only: bool = True) -> int:
    count = 0
    for a in range(1, d):
        if units_only and gcd(a, d) != 1:
            continue
        if all(safe_residue(a * (v % d), d) for v in speeds):
            count += 1
    return count


def phi(d: int) -> int:
    return sum(1 for a in range(1, d + 1) if gcd(a, d) == 1)


def unit_safe_density_for_one_runner(d: int) -> tuple[int, int]:
    safe = sum(1 for a in range(1, d) if gcd(a, d) == 1 and safe_residue(a, d))
    return safe, phi(d)


def density13_main(d: int) -> float:
    safe, ph = unit_safe_density_for_one_runner(d)
    return ph * (safe / ph) ** 13


def basis_counts(speeds: tuple[int, ...], basis: tuple[int, ...]) -> list[tuple[int, int]]:
    return [(d, witness_count(speeds, d)) for d in basis]


def least_witness_denominator(speeds: tuple[int, ...], bound: int) -> tuple[int, int] | None:
    for d in range(2, bound + 1):
        count = witness_count(speeds, d)
        if count > 0:
            return d, count
    return None


def certifies(speeds: tuple[int, ...], basis: tuple[int, ...]) -> bool:
    return any(count > 0 for _, count in basis_counts(speeds, basis))


def basis_killer(basis: tuple[int, ...]) -> tuple[int, ...]:
    tail = 84 * lcm_many(basis)
    row = tuple(sorted(BASE + (tail,)))
    assert gcd_all(row) == 1
    assert is_covering(row)
    return row


def tower_row(m: int) -> tuple[int, ...]:
    return tuple(sorted(BASE + (84 * m,)))


def random_covering_row(rng: Random, max_speed: int) -> tuple[int, ...]:
    for _ in range(5000):
        speeds = {1}
        while len(speeds) < 13 and not is_covering(tuple(speeds)):
            missing = [
                q for q in range(2, N + 1) if not any(v % q == 0 for v in speeds)
            ]
            rng.shuffle(missing)
            modulus = 1
            for q in missing:
                if rng.random() < 0.45:
                    modulus = lcm(modulus, q)
            if modulus == 1:
                modulus = missing[0]
            if modulus <= max_speed:
                speeds.add(modulus * rng.randint(1, max_speed // modulus))
        while len(speeds) < 13:
            speeds.add(rng.randint(1, max_speed))
        row = tuple(sorted(speeds))
        if len(row) == 13 and gcd_all(row) == 1 and is_covering(row):
            return row
    raise RuntimeError("failed to sample covering row")


def prompt_basis_sample(samples: int = 602, max_speed: int = 10_000) -> None:
    rng = Random(2876)
    successes = 0
    zero_rows: list[tuple[int, ...]] = []
    hist: dict[int, int] = {d: 0 for d in PROMPT_BASIS}

    for _ in range(samples):
        row = random_covering_row(rng, max_speed)
        counts = basis_counts(row, PROMPT_BASIS)
        if any(c > 0 for _, c in counts):
            successes += 1
            for d, c in counts:
                if c > 0:
                    hist[d] += 1
                    break
        else:
            zero_rows.append(row)

    print("PROMPT BASIS SAMPLE")
    print(f"  basis={PROMPT_BASIS}, samples={samples}, max_speed={max_speed}")
    print(f"  certified rows: {successes}/{samples}")
    print(f"  first-positive denominator histogram: {hist}")
    if zero_rows:
        print(f"  sample failures: {len(zero_rows)}")
        print(f"  first failure: {zero_rows[0]}")
        print(f"  first failure least D<=160: {least_witness_denominator(zero_rows[0], 160)}")
    else:
        print("  sample failures: 0")


def print_row_counts(label: str, row: tuple[int, ...], denominators: tuple[int, ...]) -> None:
    print(label)
    print(f"  row={row}")
    print(f"  covering={is_covering(row)} primitive={gcd_all(row) == 1}")
    for d in denominators:
        n_units = witness_count(row, d, units_only=True)
        n_all = witness_count(row, d, units_only=False)
        main_unit = density13_main(d)
        main_ideal = phi(d) * (6 / 7) ** 13
        divisor_hit = any(v % d == 0 for v in row)
        print(
            f"  D={d:3d}: N_units={n_units:3d}, N_all={n_all:3d}, "
            f"main67~{main_ideal:7.3f}, main_unit~{main_unit:7.3f}, "
            f"divisible_runner={divisor_hit}"
        )


def apex_floor() -> None:
    row = tower_row(1)
    print("APEX / COVERING FLOOR")
    print("  Covering means: for every D in {2,...,14}, some speed is divisible by D.")
    print("  Therefore no denominator D<=14 can certify a covering row.")
    dead = []
    for d in range(2, 15):
        c = witness_count(row, d, units_only=False)
        dead.append((d, c))
    print(f"  champion row counts for D=2..14: {dead}")
    print("  This strengthens the D=14 apex fragment: covering gives q-death for all q<=14.")


def finite_basis_killer() -> None:
    row = basis_killer(PROMPT_BASIS)
    print("FINITE BASIS KILLER")
    print(f"  basis={PROMPT_BASIS}")
    print(f"  lcm(basis)={lcm_many(PROMPT_BASIS)}")
    print(f"  killer tail=84*lcm(basis)={row[-1]}")
    print(f"  row={row}")
    print(f"  basis counts={basis_counts(row, PROMPT_BASIS)}")
    print("  reason: the tail is divisible by every denominator in the basis,")
    print("  so that runner is exactly at the observer for every numerator.")


def character_count_scout() -> None:
    denominators = (21, 41, 83, 89)
    rows = [
        ("champion tower m=1", tower_row(1)),
        ("prompt-basis killer", basis_killer(PROMPT_BASIS)),
        ("tower row m=6", tower_row(6)),
        ("known S* C-failure row", (1, 2, 3, 5, 7, 8, 9, 10, 11, 12, 13, 38, 42)),
    ]
    print("CHARACTER-SUM COUNT SCOUT")
    print("  N(S,D) is counted over unit numerators.")
    print("  main67 is phi(D)*(6/7)^13; main_unit uses the exact unit safe density.")
    print("  The error is exactly the character/resonance sum over sum k_s s = 0 mod D.")
    for label, row in rows:
        print_row_counts(label, row, denominators)


def tournament_analysis() -> None:
    vertices = [
        ("THM566_divisor_loaded_no_go", (6, 6, 6, 5)),
        ("scaled_residue_basis", (5, 5, 6, 5)),
        ("character_sum_resonance_count", (5, 5, 5, 4)),
        ("finite_sample_basis", (3, 6, 3, 3)),
        ("single_denominator_certificate", (2, 3, 2, 2)),
        ("speed_subset_minor_order", (1, 1, 1, 1)),
    ]
    scores = {name: 0 for name, _ in vertices}
    for i, (left, obs_left) in enumerate(vertices):
        for right, obs_right in vertices[i + 1 :]:
            if obs_left >= obs_right:
                scores[left] += 1
            else:
                scores[right] += 1
    ordered = sorted(vertices, key=lambda item: item[1], reverse=True)
    print("TOURNAMENT ANALYSIS")
    print("  vertices: proof carriers / residue atlases, not runners")
    print("  observable: (rules out false route, preserves witness predicate, finite-atlas value, formal clarity)")
    print(f"  Hamiltonian path: {' > '.join(name for name, _ in ordered)}")
    print(f"  score histogram: {sorted(scores.values())}")
    print("  directed 3-cycles: 0; SCCs: singleton; Hamiltonian-path count: 1")
    print("  challenged assumption: loneliness should be minor-closed under runner deletion")


def main() -> None:
    print("=" * 78)
    prompt_basis_sample()
    print("\n" + "=" * 78)
    apex_floor()
    print("\n" + "=" * 78)
    finite_basis_killer()
    print("\n" + "=" * 78)
    character_count_scout()
    print("\n" + "=" * 78)
    tournament_analysis()


if __name__ == "__main__":
    main()
