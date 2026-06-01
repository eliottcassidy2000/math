#!/usr/bin/env python3
"""Restricted tournament mappings for the Lonely Runner Conjecture.

codex-2026-06-01-S537

The S535 mappings showed a compression ladder, but the strongest compression
can throw away the arithmetic that constrains the LRC walk.  This script tests
new mappings by two criteria:

  1. How many pointed isomorphism classes are exhibited?
  2. Are the good classes pure, or mixed with unsafe states?

Every mapping here keeps LRC as a class-exhibition problem: a speed set gives a
finite wall/cell walk, the walk exhibits classes in a restricted tournament
space, and LRC asks for the walk to exhibit a class in a target subfamily.

Tournament Analysis declaration:
  pairwise observables:
    - half-turn phase difference,
    - safety danger shell,
    - safe-arc sentinel separation,
    - CRT channel majority,
    - safe-only subconfiguration.
  switches/gauges:
    - source observer threshold,
    - lexicographic danger override,
    - fixed boundary sentinels,
    - majority over residue channels.
  tie Hamiltonian path:
    - vertex index order, with observer/sentinel vertices fixed first.
  fingerprints:
    - total pointed classes,
    - good-only / bad-only / mixed fiber counts,
    - target-class counts.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from functools import reduce
from itertools import combinations, permutations
from math import gcd


ZERO = Fraction(0)
ONE = Fraction(1)
HALF = Fraction(1, 2)


def frac(x: Fraction) -> Fraction:
    return x - Fraction(x.numerator // x.denominator)


def dist0(x: Fraction) -> Fraction:
    f = frac(x)
    return min(f, ONE - f)


def clockwise_less_half(a: Fraction, b: Fraction) -> bool:
    d = frac(b - a)
    return ZERO < d < HALF


def canonicalize_fixed(adj: tuple[tuple[int, ...], ...], fixed: int = 1) -> tuple[tuple[int, ...], ...]:
    """Canonicalize a small tournament with vertices [0, fixed) marked.

    For this probe we keep fixed marks in exact positions and permute the rest.
    This gives pointed or multi-pointed isomorphism classes.
    """
    n = len(adj)
    if n - fixed <= 1:
        return adj
    best = None
    tail = range(fixed, n)
    for perm_tail in permutations(tail):
        perm = tuple(range(fixed)) + perm_tail
        candidate = tuple(tuple(adj[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if best is None or candidate < best:
            best = candidate
    return best


def compute_walls(speeds: tuple[int, ...], n: int) -> list[Fraction]:
    """LRC observer walls plus runner half-turn/collision walls."""
    walls = {ZERO}
    for v in speeds:
        for a in (1, n - 1):
            for k in range(v):
                walls.add(Fraction(k * n + a, v * n))
        for k in range(v):
            walls.add(Fraction(k, v))

    for vi, vj in combinations(speeds, 2):
        diff = abs(vi - vj)
        for k in range(diff):
            walls.add(Fraction(k, diff))
            t = Fraction(2 * k + 1, 2 * diff)
            if ZERO <= t < ONE:
                walls.add(t)
    return sorted(walls)


def phase_edge(vi: int, vj: int, t: Fraction) -> int:
    return 1 if clockwise_less_half(frac(Fraction(vi) * t), frac(Fraction(vj) * t)) else 0


def lrc_good(speeds: tuple[int, ...], n: int, t: Fraction) -> bool:
    thr = Fraction(1, n)
    return all(dist0(Fraction(v) * t) >= thr for v in speeds)


def standard_source_tournament(speeds: tuple[int, ...], n: int, t: Fraction):
    """S511/S512 baseline: observer is source exactly at lonely times."""
    thr = Fraction(1, n)
    m = len(speeds) + 1
    adj = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i == j:
                continue
            if i == 0:
                adj[i][j] = 1 if dist0(Fraction(speeds[j - 1]) * t) >= thr else 0
            elif j == 0:
                adj[i][j] = 1 if dist0(Fraction(speeds[i - 1]) * t) < thr else 0
            else:
                adj[i][j] = phase_edge(speeds[i - 1], speeds[j - 1], t)
    return tuple(tuple(r) for r in adj), 1


def binary_danger_lex_tournament(speeds: tuple[int, ...], n: int, t: Fraction):
    """Compress standard by forcing all blockers above all safe runners.

    Danger levels:
      blocker = 2, observer = 1, safe = 0.
    Larger danger beats smaller danger.  Ties among runners use half-turn phase.
    The observer is a source exactly when every runner is safe.
    """
    thr = Fraction(1, n)
    levels = [1] + [2 if dist0(Fraction(v) * t) < thr else 0 for v in speeds]
    m = len(levels)
    adj = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i == j:
                continue
            if levels[i] > levels[j]:
                adj[i][j] = 1
            elif levels[i] < levels[j]:
                adj[i][j] = 0
            elif i == 0 or j == 0:
                adj[i][j] = 1 if i < j else 0
            else:
                adj[i][j] = phase_edge(speeds[i - 1], speeds[j - 1], t)
    return tuple(tuple(r) for r in adj), 1


def ternary_danger_lex_tournament(speeds: tuple[int, ...], n: int, t: Fraction):
    """A shell version: blockers > observer > near-safe > deep-safe.

    This gives a stricter restricted family than the baseline while preserving
    more geometry than the binary collapse.
    """
    thr = Fraction(1, n)
    levels = [2]
    for v in speeds:
        d = dist0(Fraction(v) * t)
        if d < thr:
            levels.append(3)
        elif d < 2 * thr:
            levels.append(1)
        else:
            levels.append(0)
    m = len(levels)
    adj = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i == j:
                continue
            if levels[i] > levels[j]:
                adj[i][j] = 1
            elif levels[i] < levels[j]:
                adj[i][j] = 0
            elif i == 0 or j == 0:
                adj[i][j] = 1 if i < j else 0
            else:
                adj[i][j] = phase_edge(speeds[i - 1], speeds[j - 1], t)
    return tuple(tuple(r) for r in adj), 1


def two_sentinel_safe_arc_tournament(speeds: tuple[int, ...], n: int, t: Fraction):
    """Half-turn tournament with two fixed safe-arc boundary sentinels.

    Vertices:
      0 = left safe boundary 1/n
      1 = right safe boundary 1 - 1/n
      2.. = runners

    LRC is target membership: all runner positions lie in the marked safe arc
    from vertex 0 to vertex 1.  The tournament has two fixed marks, so the
    isomorphism universe is far smaller than unmarked tournaments on n+1
    vertices in any proof that respects the sentinels.
    """
    positions = [Fraction(1, n), ONE - Fraction(1, n)] + [frac(Fraction(v) * t) for v in speeds]
    m = len(positions)
    adj = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i == j:
                continue
            if clockwise_less_half(positions[i], positions[j]):
                adj[i][j] = 1
            elif frac(positions[j] - positions[i]) == HALF:
                adj[i][j] = 1 if i < j else 0
    return tuple(tuple(r) for r in adj), 2


def safe_only_subtournament(speeds: tuple[int, ...], n: int, t: Fraction):
    """Variable-size tournament on currently safe runners only.

    The class key is (safe_count, safe half-turn iso class).  LRC target is the
    top layer safe_count = n-1.
    """
    thr = Fraction(1, n)
    safe = [v for v in speeds if dist0(Fraction(v) * t) >= thr]
    m = len(safe)
    adj = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i != j:
                adj[i][j] = phase_edge(safe[i], safe[j], t)
    return tuple(tuple(r) for r in adj), 0


def crt_channel_tournament(speeds: tuple[int, ...], n: int, t: Fraction):
    """CRT/channel majority tournament, using the largest proper divisor channel.

    For n=14 this is the mod-7 channel picture; for n=6 it is mod-3.
    The observer is fixed and a channel is safe only if all its runners are safe.
    """
    divisors = [d for d in range(2, n) if n % d == 0]
    if not divisors:
        return None, 0
    modulus = max(divisors)
    thr = Fraction(1, n)
    classes: dict[int, list[int]] = defaultdict(list)
    for v in speeds:
        classes[v % modulus].append(v)
    residues = sorted(classes)
    m = len(residues) + 1
    adj = [[0] * m for _ in range(m)]
    for idx, r in enumerate(residues, start=1):
        all_safe = all(dist0(Fraction(v) * t) >= thr for v in classes[r])
        if all_safe:
            adj[0][idx] = 1
        else:
            adj[idx][0] = 1
    for ai, ra in enumerate(residues, start=1):
        for bi, rb in enumerate(residues, start=1):
            if ai == bi:
                continue
            wins = 0
            total = 0
            for va in classes[ra]:
                for vb in classes[rb]:
                    total += 1
                    wins += phase_edge(va, vb, t)
            if 2 * wins > total:
                adj[ai][bi] = 1
            elif 2 * wins == total:
                adj[ai][bi] = 1 if ai < bi else 0
    return tuple(tuple(r) for r in adj), 1


MAPPINGS = [
    ("standard_source", standard_source_tournament),
    ("binary_danger_lex", binary_danger_lex_tournament),
    ("ternary_danger_lex", ternary_danger_lex_tournament),
    ("two_sentinel_safe_arc", two_sentinel_safe_arc_tournament),
    ("safe_only_subtournament", safe_only_subtournament),
    ("crt_channel", crt_channel_tournament),
]


def class_key(name: str, adj: tuple[tuple[int, ...], ...], fixed: int, speeds: tuple[int, ...], n: int, t: Fraction):
    if name == "safe_only_subtournament":
        return (len(adj), canonicalize_fixed(adj, 0))
    return canonicalize_fixed(adj, fixed)


def summarize_mapping(name: str, seen: dict[object, set[bool]]) -> str:
    good_only = sum(flags == {True} for flags in seen.values())
    bad_only = sum(flags == {False} for flags in seen.values())
    mixed = sum(flags == {True, False} for flags in seen.values())
    return (
        f"{name:24s} classes={len(seen):4d} "
        f"good_only={good_only:4d} bad_only={bad_only:4d} mixed={mixed:4d}"
    )


def run_probe() -> None:
    print("Restricted LRC tournament mappings — codex S537")
    print("=" * 72)
    print("Counts are over open wall cells; boundary compactification is not included.")
    print("A good-only class is a target class never seen in an unsafe state.")
    print()

    configs = {
        4: (12, None),
        5: (9, None),
        6: (8, 80),
        7: (9, 40),
    }
    for n, (max_speed, max_sets) in configs.items():
        speed_sets = [
            combo
            for combo in combinations(range(1, max_speed + 1), n - 1)
            if reduce(gcd, combo) == 1
        ]
        if max_sets is not None:
            speed_sets = speed_sets[:max_sets]
        print(f"n={n}: max_speed={max_speed}, speed_sets={len(speed_sets)}")
        mapping_seen: dict[str, dict[object, set[bool]]] = {name: defaultdict(set) for name, _ in MAPPINGS}
        state_count = 0
        good_state_count = 0
        for speeds in speed_sets:
            walls = compute_walls(speeds, n)
            walls_ext = walls + [ONE]
            for a, b in zip(walls_ext, walls_ext[1:]):
                if b <= a:
                    continue
                t = (a + b) / 2
                good = lrc_good(speeds, n, t)
                state_count += 1
                good_state_count += int(good)
                for name, maker in MAPPINGS:
                    adj, fixed = maker(speeds, n, t)
                    if adj is None:
                        continue
                    key = class_key(name, adj, fixed, speeds, n, t)
                    mapping_seen[name][key].add(good)

        print(f"  open states={state_count}, good states={good_state_count}")
        for name, _ in MAPPINGS:
            print("  " + summarize_mapping(name, mapping_seen[name]))
        print()

    print("Interpretation")
    print("-" * 72)
    print("standard_source is exact but still too close to the raw observer-marked")
    print("  quotient.  Its good classes are pure because observer-source is visible.")
    print("binary_danger_lex forces all blockers above all safe runners.  It usually")
    print("  collapses many unsafe states while preserving pure good target classes.")
    print("ternary_danger_lex keeps near-safe versus deep-safe information; it trades")
    print("  compression for a more geometric source neighborhood.")
    print("two_sentinel_safe_arc makes the target a two-marked safe-arc class problem.")
    print("safe_only_subtournament is the harshest quotient: LRC is exhibition of")
    print("  the top safe-count layer.  It is pure by construction but loses blockers.")
    print("crt_channel is the arithmetic quotient.  It can be tiny, but its coarse")
    print("  channel memory needs coupled debt labels to become a proof tool.")


if __name__ == "__main__":
    run_probe()
