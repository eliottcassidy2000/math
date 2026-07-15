#!/usr/bin/env python3
"""Independent exact endpoint and Tournament Analysis replay for THM-845."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations
from pathlib import Path


DELTA = F(1, 13)

# (missing labels, sorted (owner, speed) word, expected exact maximin)
TERMINALS = (
    ((1, 3, 4, 6, 10), ((1, 14), (4, 17), (6, 19), (10, 23), (3, 29)), F(3, 28)),
    ((1, 3, 5, 8, 9), ((1, 14), (3, 16), (5, 18), (8, 21), (9, 22)), F(2, 17)),
    ((1, 3, 5, 9, 11), ((1, 14), (3, 16), (5, 18), (9, 22), (11, 24)), F(1, 10)),
    ((2, 3, 4, 10, 12), ((3, 16), (4, 17), (10, 23), (12, 25), (2, 28)), F(1, 10)),
    ((2, 3, 6, 8, 10), ((10, 23), (2, 28), (3, 29), (6, 32), (8, 34)), F(2, 15)),
    ((3, 4, 5, 6, 8), ((3, 16), (4, 17), (5, 18), (6, 19), (8, 21)), F(3, 26)),
    ((3, 4, 5, 6, 8), ((4, 17), (6, 19), (3, 29), (5, 31), (8, 34)), F(1, 8)),
    ((3, 4, 6, 8, 10), ((4, 17), (6, 19), (8, 21), (10, 23), (3, 29)), F(1, 8)),
    ((4, 5, 6, 9, 12), ((4, 17), (5, 18), (6, 19), (12, 25), (9, 35)), F(4, 37)),
)


def safe_bands(speed: int) -> tuple[tuple[F, F], ...]:
    return tuple(
        (F(13 * k + 1, 13 * speed), F(13 * (k + 1) - 1, 13 * speed))
        for k in range(speed)
    )


def meet(
    left: tuple[tuple[F, F], ...], right: tuple[tuple[F, F], ...]
) -> tuple[tuple[F, F], ...]:
    out: list[tuple[F, F]] = []
    i = j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo < hi:
            out.append((lo, hi))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return tuple(out)


def safe_components(speeds: tuple[int, ...]) -> tuple[tuple[F, F], ...]:
    current = ((F(0), F(1)),)
    for speed in sorted(speeds):
        current = meet(current, safe_bands(speed))
    return current


def measure(intervals: tuple[tuple[F, F], ...]) -> F:
    return sum((hi - lo for lo, hi in intervals), F(0))


def circle_norm(numerator: int, denominator: int) -> F:
    residue = numerator % denominator
    return F(min(residue, denominator - residue), denominator)


def exact_maximin(speeds: tuple[int, ...]) -> tuple[F, tuple[F, ...]]:
    denominators = {2 * u for u in speeds}
    denominators |= {
        q
        for u, v in combinations(speeds, 2)
        for q in (u + v, abs(u - v))
        if q
    }
    best = F(0)
    witnesses: set[F] = set()
    for q in sorted(denominators):
        for a in range(q):
            value = min(circle_norm(u * a, q) for u in speeds)
            t = F(a, q)
            if value > best:
                best = value
                witnesses = {t}
            elif value == best:
                witnesses.add(t)
    return best, tuple(sorted(witnesses))


def tournament(
    values: tuple[F, ...], tie_keys: tuple[tuple[int, int], ...]
) -> tuple[tuple[bool, ...], ...]:
    n = len(values)
    edge = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        i_wins = values[i] > values[j] or (
            values[i] == values[j] and tie_keys[i] < tie_keys[j]
        )
        edge[i][j] = i_wins
        edge[j][i] = not i_wins
    return tuple(tuple(row) for row in edge)


def fingerprint(edge: tuple[tuple[bool, ...], ...]) -> tuple:
    n = len(edge)
    scores = tuple(sorted(Counter(sum(row) for row in edge).items()))
    cycles = sum(
        (edge[i][j] and edge[j][k] and edge[k][i])
        or (edge[j][i] and edge[i][k] and edge[k][j])
        for i, j, k in combinations(range(n), 3)
    )
    reach = [list(row) for row in edge]
    for i in range(n):
        reach[i][i] = True
    for k in range(n):
        for i in range(n):
            for j in range(n):
                reach[i][j] |= reach[i][k] and reach[k][j]
    unused = set(range(n))
    scc: list[int] = []
    while unused:
        i = min(unused)
        component = {j for j in unused if reach[i][j] and reach[j][i]}
        unused -= component
        scc.append(len(component))
    hp = sum(
        all(edge[path[i]][path[i + 1]] for i in range(n - 1))
        for path in permutations(range(n))
    )
    return scores, cycles, tuple(sorted(scc, reverse=True)), hp


def flips(
    first: tuple[tuple[bool, ...], ...], second: tuple[tuple[bool, ...], ...]
) -> int:
    return sum(first[i][j] != second[i][j] for i, j in combinations(range(5), 2))


def main() -> None:
    print("THM845_HAMMING_FIVE_TERMINAL_ENDPOINT_CROSSCHECK")
    print(f"terminal_rows={len(TERMINALS)} delta={DELTA}")
    certificate: list[str] = []
    total_flips = 0
    raw_cycle_rows = conditional_cycle_rows = 0

    for index, (missing, word, expected) in enumerate(TERMINALS):
        assert tuple(sorted(missing)) == missing
        owners = tuple(owner for owner, _ in word)
        replacements = tuple(speed for _, speed in word)
        assert set(owners) == set(missing)
        assert replacements == tuple(sorted(replacements))
        assert all(speed % 13 == owner for owner, speed in word)
        core = tuple(q for q in range(1, 13) if q not in missing)
        packet = tuple(sorted(core + replacements))

        components = safe_components(packet)
        assert components
        maximum, witnesses = exact_maximin(packet)
        assert maximum == expected and maximum > DELTA

        core_components = safe_components(core)
        longest = max(core_components, key=lambda interval: (interval[1] - interval[0], -interval[0]))
        longest_word = (longest,)
        raw_losses = tuple(
            (longest[1] - longest[0]) - measure(meet(longest_word, safe_bands(speed)))
            for speed in replacements
        )
        full_measure = measure(components)
        conditional_losses = tuple(
            measure(safe_components(tuple(sorted(core + replacements[:i] + replacements[i + 1 :]))))
            - full_measure
            for i in range(5)
        )

        tie_keys = tuple((speed, owner) for owner, speed in word)
        raw_tournament = tournament(raw_losses, tie_keys)
        conditional_tournament = tournament(conditional_losses, tie_keys)
        raw_fp = fingerprint(raw_tournament)
        conditional_fp = fingerprint(conditional_tournament)
        edge_flips = flips(raw_tournament, conditional_tournament)
        total_flips += edge_flips
        raw_cycle_rows += raw_fp[1] > 0
        conditional_cycle_rows += conditional_fp[1] > 0

        row = (
            f"row[{index}] missing={missing} word={word} M={maximum} "
            f"witnesses={witnesses} components={len(components)} "
            f"longest={max(hi-lo for lo,hi in components)} "
            f"raw_losses={raw_losses} conditional_losses={conditional_losses} "
            f"raw_fp={raw_fp} conditional_fp={conditional_fp} flips={edge_flips}"
        )
        certificate.append(row)
        print(row)

    digest = sha256("\n".join(certificate).encode()).hexdigest()
    source_digest = sha256(Path(__file__).read_bytes()).hexdigest()
    print("\nTOURNAMENT_ANALYSIS")
    print("vertices=five_remaining_comb_obligations_not_runners")
    print("pair_observable=raw_erosion_of_current_longest_core_component")
    print("switch=conditional_marginal_erosion_after_the_other_four_combs")
    print("tie_path=increasing_speed_then_label")
    print(f"edge_flips_total={total_flips}")
    print(f"rows_with_raw_cycles={raw_cycle_rows} rows_with_conditional_cycles={conditional_cycle_rows}")
    print("preserves=planning_priority destroys=higher_order_overlap_component_splitting_cover_truth")
    print("challenged_vertex_set=proof_obligations_are_more_faithful_than_runners_or_residues")
    print(f"certificate_sha256={digest}")
    print(f"source_sha256={source_digest}")
    print("PASS: nine terminal rows are independently loose; tournaments remain telemetry")


if __name__ == "__main__":
    main()
