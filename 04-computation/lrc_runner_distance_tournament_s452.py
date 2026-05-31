#!/usr/bin/env python3
"""
lrc_runner_distance_tournament_s452.py

codex-2026-05-31 S452

Explore how tournament structure appears when LRC keeps more than the usual
stationary-runner distance data.

Data lifts:
  star     = distances from runner 0 to every other runner (usual LRC frame)
  cycle    = adjacent circular gaps, i.e. the two nearest neighbors
  complete = all pairwise distances, one coordinate per tournament edge

The two-nearest-neighbor lift is exact for a static configuration: a runner is
lonely iff the two adjacent circular gaps around it are both at least 1/n.
No-lonely safe-gap masks are independent sets in the cycle graph C_n.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations


def mod1(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def norm_dist(x: Fraction) -> Fraction:
    y = mod1(x)
    return min(y, 1 - y)


def fmt(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def hist(values: tuple[int, ...]) -> str:
    pieces = []
    for value in sorted(set(values)):
        count = sum(1 for x in values if x == value)
        pieces.append(f"{value}^{count}" if count > 1 else str(value))
    return " ".join(pieces)


def initial(n: int) -> tuple[int, ...]:
    return tuple(range(1, n))


def forbidden_intervals(speeds: tuple[int, ...]) -> list[tuple[Fraction, Fraction]]:
    n = len(speeds) + 1
    pieces: list[tuple[Fraction, Fraction]] = []
    for v in speeds:
        half = Fraction(1, n * v)
        for m in range(v):
            center = Fraction(m, v)
            lo = center - half
            hi = center + half
            if lo < 0:
                pieces.append((lo + 1, Fraction(1)))
                pieces.append((Fraction(0), hi))
            elif hi > 1:
                pieces.append((lo, Fraction(1)))
                pieces.append((Fraction(0), hi - 1))
            else:
                pieces.append((lo, hi))
    return pieces


def merge_intervals(pieces: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    pieces = sorted(pieces)
    if not pieces:
        return []
    merged: list[list[Fraction]] = [[pieces[0][0], pieces[0][1]]]
    for lo, hi in pieces[1:]:
        cur = merged[-1]
        if lo <= cur[1]:
            cur[1] = max(cur[1], hi)
        else:
            merged.append([lo, hi])
    return [(lo, hi) for lo, hi in merged]


def complement_gaps(merged: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    if not merged:
        return [(Fraction(0), Fraction(1))]
    gaps: list[tuple[Fraction, Fraction]] = []
    for (_, hi), (lo2, _) in zip(merged, merged[1:]):
        if lo2 > hi:
            gaps.append((hi, lo2))
    first_lo = merged[0][0]
    last_hi = merged[-1][1]
    if first_lo > 0 or last_hi < 1:
        gaps.append((last_hi, first_lo + 1))
    return gaps


def widest_lonely_time(speeds: tuple[int, ...]) -> tuple[Fraction, Fraction]:
    merged = merge_intervals(forbidden_intervals(speeds))
    gaps = complement_gaps(merged)
    lo, hi = max(gaps, key=lambda ab: ab[1] - ab[0])
    return mod1((lo + hi) / 2), hi - lo


def cycle_independence(n: int, x: int) -> int:
    total = 0
    for mask in range(1 << n):
        ok = True
        for i in range(n):
            if (mask >> i) & 1 and (mask >> ((i + 1) % n)) & 1:
                ok = False
                break
        if ok:
            total += x ** mask.bit_count()
    return total


@dataclass(frozen=True)
class Snapshot:
    label: str
    n: int
    time: Fraction
    threshold: Fraction
    stationary_left_gap: Fraction
    stationary_right_gap: Fraction
    stationary_two_nearest: tuple[Fraction, Fraction]
    lonely_vertices: tuple[int, ...]
    close_degrees: tuple[tuple[int, int], ...]
    circular_score_sequence: tuple[int, ...]
    circular_three_cycles: int
    safe_gap_mask: str
    max_gap: Fraction
    gap_ratio: Fraction


def snapshot(label: str, speeds: tuple[int, ...], time: Fraction) -> Snapshot:
    velocities = (0,) + tuple(speeds)
    n = len(velocities)
    threshold = Fraction(1, n)
    positions = [(i, mod1(Fraction(v) * time)) for i, v in enumerate(velocities)]
    ordered = sorted(positions, key=lambda item: (item[1], item[0]))
    order_labels = [i for i, _ in ordered]
    pos_by_label = {i: p for i, p in positions}
    idx_by_label = {label_i: idx for idx, label_i in enumerate(order_labels)}

    gaps: list[Fraction] = []
    for idx in range(n):
        here = ordered[idx][1]
        nxt = ordered[(idx + 1) % n][1]
        gaps.append(mod1(nxt - here))

    lonely: list[int] = []
    for label_i in order_labels:
        idx = idx_by_label[label_i]
        left = gaps[(idx - 1) % n]
        right = gaps[idx]
        if left >= threshold and right >= threshold:
            lonely.append(label_i)

    stationary_idx = idx_by_label[0]
    stationary_left = gaps[(stationary_idx - 1) % n]
    stationary_right = gaps[stationary_idx]
    distances_from_zero = sorted(norm_dist(pos_by_label[j] - pos_by_label[0]) for j in range(1, n))
    stationary_two = (distances_from_zero[0], distances_from_zero[1])

    close_degrees = []
    for i in range(n):
        deg = sum(
            1
            for j in range(n)
            if i != j and norm_dist(pos_by_label[i] - pos_by_label[j]) < threshold
        )
        close_degrees.append((i, deg))

    # Circular distance tournament: i beats j if j lies within the clockwise
    # half-circle from i.  Antipodal ties, if any, are omitted from the score.
    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        cw = mod1(pos_by_label[j] - pos_by_label[i])
        if cw == 0 or cw == Fraction(1, 2):
            continue
        if cw < Fraction(1, 2):
            adj[i][j] = True
        else:
            adj[j][i] = True
    scores = tuple(sum(1 for j in range(n) if adj[i][j]) for i in range(n))

    tri = 0
    for a, b, c in combinations(range(n), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[a][c] and adj[c][b] and adj[b][a]):
            tri += 1

    safe_mask = "".join("1" if g >= threshold else "0" for g in gaps)
    max_gap = max(gaps)
    return Snapshot(
        label=label,
        n=n,
        time=time,
        threshold=threshold,
        stationary_left_gap=stationary_left,
        stationary_right_gap=stationary_right,
        stationary_two_nearest=stationary_two,
        lonely_vertices=tuple(sorted(lonely)),
        close_degrees=tuple(close_degrees),
        circular_score_sequence=tuple(sorted(scores)),
        circular_three_cycles=tri,
        safe_gap_mask=safe_mask,
        max_gap=max_gap,
        gap_ratio=max_gap / threshold,
    )


def print_data_lift_table() -> None:
    print("Runner-data lifts")
    print("=" * 88)
    print(f"{'n':>3} {'star distances':>14} {'cycle gaps':>12} {'complete distances':>20} polygonal name")
    for n in (8, 14, 16):
        star = n - 1
        cycle = n
        complete = n * (n - 1) // 2
        print(f"{n:>3} {star:>14} {cycle:>12} {complete:>20} T_{n-1}={complete}")
    print()
    print("Interpretation")
    print("=" * 88)
    print("- star: the usual stationary-runner frame")
    print("- cycle: the two nearest neighbors; loneliness is two adjacent safe gaps")
    print("- complete: all pairwise distances, indexed by tournament edges")


def print_static_mask_table() -> None:
    print()
    print("Static two-nearest obstruction masks")
    print("=" * 88)
    print(f"{'n':>3} {'I(C_n,1)':>12} {'feasible no-lonely':>20} {'I(C_n,2)':>12} {'formula':>12}")
    for n in (8, 14, 16):
        c1 = cycle_independence(n, 1)
        c2 = cycle_independence(n, 2)
        print(f"{n:>3} {c1:>12} {c1 - 1:>20} {c2:>12} {2 ** n + (-1) ** n:>12}")
    print()
    print(
        "A safe gap is one with length at least 1/n.  A runner is lonely "
        "exactly when the two adjacent gaps around it are safe.  Thus a "
        "configuration with no lonely runner has safe gaps forming an "
        "independent set in C_n; the empty mask is not geometrically feasible "
        "because all gaps would be shorter than average."
    )


def print_snapshot_table() -> None:
    examples: list[tuple[str, tuple[int, ...], Fraction]] = []
    examples.append(("initial n=14 boundary", initial(14), Fraction(1, 14)))
    seven = (1, 7, 14, 21, 28, 35, 49, 56, 63, 70, 77, 84, 91)
    s380 = (1, 14, 28, 42, 56, 70, 98, 112, 126, 140, 154, 168, 182)
    for name, speeds in (("n14 seven-ladder max gap", seven), ("n14 S380 max gap", s380)):
        t, _ = widest_lonely_time(speeds)
        examples.append((name, speeds, t))
    examples.append(("initial n=16 boundary", initial(16), Fraction(1, 16)))

    print()
    print("Exact runner-distance snapshots")
    print("=" * 120)
    print(
        f"{'example':<28} {'t':>12} {'left0':>10} {'right0':>10} "
        f"{'two-nearest0':>20} {'lonely':>18} {'gap/th':>8} {'3cyc':>6} safe-gap-mask"
    )
    for label, speeds, time in examples:
        row = snapshot(label, speeds, time)
        two = f"{fmt(row.stationary_two_nearest[0])},{fmt(row.stationary_two_nearest[1])}"
        lonely = ",".join(str(v) for v in row.lonely_vertices[:8])
        if len(row.lonely_vertices) > 8:
            lonely += ",..."
        print(
            f"{row.label:<28} {fmt(row.time):>12} {fmt(row.stationary_left_gap):>10} "
            f"{fmt(row.stationary_right_gap):>10} {two:>20} {lonely:>18} "
            f"{fmt(row.gap_ratio):>8} {row.circular_three_cycles:>6} {row.safe_gap_mask}"
        )

    print()
    print("Complete pairwise-distance / tournament summaries")
    print("=" * 120)
    print(f"{'example':<28} {'close-degree hist':<28} {'score sequence hist':<28} {'3cyc':>6} {'comment'}")
    for label, speeds, time in examples:
        row = snapshot(label, speeds, time)
        close_values = tuple(deg for _, deg in row.close_degrees)
        print(
            f"{row.label:<28} {hist(close_values):<28} "
            f"{hist(row.circular_score_sequence):<28} {row.circular_three_cycles:>6} "
            "complete lift"
        )


def print_methodology_table() -> None:
    print()
    print("Methodology reframes")
    print("=" * 88)
    rows = [
        (
            "endpoint cover",
            "single stationary star",
            "THM-357 boundary peel",
        ),
        (
            "two-neighbor cycle",
            "gap polygon C_n",
            "safe masks are independent sets; import Zeckendorf/OCF",
        ),
        (
            "pairwise lift",
            "complete graph K_n",
            "distances indexed by tournament edges T_{n-1}",
        ),
        (
            "order chambers",
            "circular orders between collisions",
            "adjacent swaps act like tournament arc flips",
        ),
        (
            "distance tournament",
            "i -> j when j is in i's clockwise half-circle",
            "score/cycle data measures global crowding",
        ),
    ]
    print(f"{'route':<22} {'object':<36} payoff")
    for route, obj, payoff in rows:
        print(f"{route:<22} {obj:<36} {payoff}")


def print_synthesis() -> None:
    print()
    print("S452 synthesis")
    print("=" * 88)
    print(
        "The two-nearest-neighbor lift is exact: a runner is lonely iff the "
        "two adjacent gaps around it are both at least 1/n."
    )
    print(
        "Therefore no-lonely static configurations are controlled by "
        "independent sets in the gap cycle C_n.  This is the clean place "
        "where Zeckendorf/OCF enters LRC."
    )
    print(
        "The full pairwise-distance lift has one coordinate per edge of K_n, "
        "the same triangular-number count as tournament data.  This suggests "
        "a proof route through chamber graphs, adjacent swaps, good cuts, and "
        "cycle/score invariants rather than scalar distance alone."
    )


def main() -> None:
    print_data_lift_table()
    print_static_mask_table()
    print_snapshot_table()
    print_methodology_table()
    print_synthesis()


if __name__ == "__main__":
    main()
