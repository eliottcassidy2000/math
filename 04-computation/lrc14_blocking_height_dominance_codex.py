#!/usr/bin/env python3
"""Does dominance grow with LRC14 blocking height?

The prompt asks whether dominance grows with blocking height, remembering that
tournaments model dominance relations.  Here "blocking height" is the
dilated-band ladder height from the latest HYP-2471/KPS work:

    h(S) = least q with an uncovered unit numerator a mod q.

For q < h(S), the row S blocks every unit numerator.  For two speeds v,w, define
a dominance comparison over the pre-height band shells by:

    v -> w iff sum_q |U_v(q) \ U_w(q)| >= sum_q |U_w(q) \ U_v(q)|,

where U_v(q) is the dilated danger band of v on (Z/q)^*.  The fixed tie
Hamiltonian path is the row order.

This is a deliberately local, proof-facing tournament.  Vertices are speeds,
not runners in continuous time.  It preserves the cover/load relation that makes
a shell blocked, and it discards continuous owner geometry; Bprime/owner checks
remain a separate side channel.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from math import gcd, sqrt
from pathlib import Path
import random
import sys


HERE = Path(__file__).resolve().parent
sys.path.append(str(HERE))

from lrc14_pisano_band_ladder_codex import CORE, MAX_R  # noqa: E402


N = 14
Q_START = 14
HMAX = 80
CORE_T = tuple(CORE)


def band_set(q: int, n: int = N) -> set[int]:
    h = q // n
    return {r for r in range(q) if min(r, q - r) <= h}


def unit_list(q: int) -> tuple[int, ...]:
    return tuple(a for a in range(1, q) if gcd(a, q) == 1)


def cover_mask_for_speed(v: int, q: int, units: tuple[int, ...], band: set[int]) -> int:
    mask = 0
    for i, a in enumerate(units):
        if (v * a) % q in band:
            mask |= 1 << i
    return mask


def q_cover_masks(row: tuple[int, ...], q: int) -> tuple[int, tuple[int, ...], int]:
    units = unit_list(q)
    band = band_set(q)
    masks = tuple(cover_mask_for_speed(v, q, units, band) for v in row)
    full = (1 << len(units)) - 1
    return full, masks, len(units)


def witness_deficit(row: tuple[int, ...], q: int) -> int:
    full, masks, _ = q_cover_masks(row, q)
    covered = 0
    for mask in masks:
        covered |= mask
    return (full & ~covered).bit_count()


def ladder_height(row: tuple[int, ...], qmax: int = HMAX) -> int | None:
    for q in range(2, qmax + 1):
        if witness_deficit(row, q) > 0:
            return q
    return None


def hamiltonian_paths(adj: tuple[tuple[int, ...], ...]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if not count:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += count
    return sum(dp[-1])


def scc_sizes(adj: tuple[tuple[int, ...], ...]) -> list[int]:
    n = len(adj)
    reach = [[bool(adj[i][j]) for j in range(n)] for i in range(n)]
    for i in range(n):
        reach[i][i] = True
    for k in range(n):
        for i in range(n):
            if reach[i][k]:
                for j in range(n):
                    reach[i][j] = reach[i][j] or reach[k][j]
    seen: set[int] = set()
    out = []
    for i in range(n):
        if i in seen:
            continue
        comp = {j for j in range(n) if reach[i][j] and reach[j][i]}
        seen |= comp
        out.append(len(comp))
    return sorted(out, reverse=True)


@dataclass(frozen=True)
class DominanceProfile:
    label: str
    row: tuple[int, ...]
    height: int | None
    pre_shells: int
    blocked_unit_mass: int
    top_speed: int | None
    top_cover_share: float
    top_unique_share: float
    mean_pair_margin: float
    mean_pair_margin_per_shell: float
    mean_pair_margin_norm: float
    score_spread: int
    top_outdegree: int
    directed_3cycles: int
    score_hist: dict[int, int]
    scc_sizes: list[int]
    hamiltonian_paths: int | None
    cover_loads: tuple[tuple[int, int], ...]
    unique_loads: tuple[tuple[int, int], ...]


def dominance_profile(
    label: str,
    row: tuple[int, ...],
    q_start: int = Q_START,
    hmax: int = HMAX,
    full_fingerprint: bool = False,
) -> DominanceProfile:
    h = ladder_height(row, hmax)
    if h is None or h <= q_start:
        n = len(row)
        return DominanceProfile(
            label,
            row,
            h,
            0,
            0,
            None,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0,
            0,
            0,
            {},
            [],
            None,
            tuple(),
            tuple(),
        )

    n = len(row)
    pair_margin = [[0] * n for _ in range(n)]
    cover_load = [0] * n
    unique_load = [0] * n
    blocked_unit_mass = 0
    pre_shells = 0

    for q in range(q_start, h):
        full, masks, unit_count = q_cover_masks(row, q)
        covered = 0
        for mask in masks:
            covered |= mask
        if covered != full:
            # h should be first leak; if a smaller q leaks, stop trusting this row.
            continue
        pre_shells += 1
        blocked_unit_mass += unit_count
        for i, mi in enumerate(masks):
            cover_load[i] += mi.bit_count()
            others = 0
            for j, mj in enumerate(masks):
                if i != j:
                    others |= mj
            unique_load[i] += (mi & ~others).bit_count()
        for i, j in combinations(range(n), 2):
            mi = masks[i]
            mj = masks[j]
            pair_margin[i][j] += (mi & ~mj).bit_count() - (mj & ~mi).bit_count()
            pair_margin[j][i] = -pair_margin[i][j]

    adj = [[0] * n for _ in range(n)]
    outdeg = [0] * n
    margins = []
    for i, j in combinations(range(n), 2):
        margin = pair_margin[i][j]
        margins.append(abs(margin))
        if margin >= 0:
            adj[i][j] = 1
            outdeg[i] += 1
        else:
            adj[j][i] = 1
            outdeg[j] += 1

    directed_3cycles = 0
    for i, j, k in combinations(range(n), 3):
        degs = [0, 0, 0]
        verts = (i, j, k)
        for a, b in ((0, 1), (0, 2), (1, 2)):
            x = verts[a]
            y = verts[b]
            if adj[x][y]:
                degs[a] += 1
            else:
                degs[b] += 1
        if sorted(degs) == [1, 1, 1]:
            directed_3cycles += 1

    max_cover = max(cover_load) if cover_load else 0
    total_cover = sum(cover_load)
    max_unique = max(unique_load) if unique_load else 0
    total_unique = sum(unique_load)
    top_index = cover_load.index(max_cover) if max_cover else None
    top_speed = row[top_index] if top_index is not None else None
    mean_margin = sum(margins) / len(margins) if margins else 0.0
    norm = blocked_unit_mass if blocked_unit_mass else 1
    adj_tuple = tuple(tuple(x for x in row_adj) for row_adj in adj)

    return DominanceProfile(
        label=label,
        row=row,
        height=h,
        pre_shells=pre_shells,
        blocked_unit_mass=blocked_unit_mass,
        top_speed=top_speed,
        top_cover_share=max_cover / total_cover if total_cover else 0.0,
        top_unique_share=max_unique / total_unique if total_unique else 0.0,
        mean_pair_margin=mean_margin,
        mean_pair_margin_per_shell=mean_margin / pre_shells if pre_shells else 0.0,
        mean_pair_margin_norm=mean_margin / norm,
        score_spread=max(outdeg) - min(outdeg) if outdeg else 0,
        top_outdegree=max(outdeg) if outdeg else 0,
        directed_3cycles=directed_3cycles,
        score_hist=dict(sorted(Counter(outdeg).items())),
        scc_sizes=scc_sizes(adj_tuple) if full_fingerprint else [],
        hamiltonian_paths=hamiltonian_paths(adj_tuple) if full_fingerprint else None,
        cover_loads=tuple(sorted(zip(row, cover_load), key=lambda x: (-x[1], x[0]))),
        unique_loads=tuple(sorted(zip(row, unique_load), key=lambda x: (-x[1], x[0]))),
    )


def pearson(xs: list[float], ys: list[float]) -> float:
    if len(xs) < 2:
        return 0.0
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0:
        return 0.0
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sqrt(vx * vy)


def ranks(vals: list[float]) -> list[float]:
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    out = [0.0] * len(vals)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and vals[order[j]] == vals[order[i]]:
            j += 1
        rank = (i + j - 1) / 2 + 1
        for k in range(i, j):
            out[order[k]] = rank
        i = j
    return out


def spearman(xs: list[float], ys: list[float]) -> float:
    return pearson(ranks(xs), ranks(ys))


def one_stranger_rows() -> list[tuple[str, tuple[int, ...]]]:
    rows = []
    for r in range(1, MAX_R + 1):
        if r % 7 == 0:
            continue
        row = tuple(sorted(CORE_T + (r,)))
        if len(set(row)) != 13:
            continue
        g = 0
        for v in row:
            g = gcd(g, v)
        if g != 1:
            continue
        rows.append((f"S({r})", row))
    return rows


def random_primitive_mult14_rows(count: int = 120, seed: int = 2481) -> list[tuple[str, tuple[int, ...]]]:
    rng = random.Random(seed)
    rows = []
    attempts = 0
    while len(rows) < count and attempts < 20000:
        attempts += 1
        row = tuple(sorted(rng.sample(range(1, 181), 13)))
        g = 0
        for v in row:
            g = gcd(g, v)
        if g == 1 and any(v % 14 == 0 for v in row):
            rows.append((f"R{len(rows):03d}", row))
    return rows


def named_rows() -> list[tuple[str, tuple[int, ...]]]:
    return [
        ("one-stranger-611", tuple(sorted(CORE_T + (611,)))),
        ("one-stranger-702", tuple(sorted(CORE_T + (702,)))),
        ("one-stranger-793", tuple(sorted(CORE_T + (793,)))),
        ("one-stranger-962", tuple(sorted(CORE_T + (962,)))),
        ("one-stranger-1053", tuple(sorted(CORE_T + (1053,)))),
        ("HYP2470-E1-q33", (7, 14, 21, 35, 49, 63, 70, 77, 91, 322, 350, 504, 936)),
        ("HYP2470-E2-q31", (7, 14, 21, 28, 35, 49, 63, 77, 91, 119, 700, 1008, 1066)),
        ("HYP2471-Q31-E1", (7, 14, 21, 35, 49, 63, 70, 77, 91, 322, 350, 504, 832)),
        ("HYP2471-Q31-E2", (7, 14, 21, 28, 35, 49, 63, 77, 91, 119, 156, 700, 1008)),
    ]


def print_correlation_block(title: str, profiles: list[DominanceProfile]) -> None:
    active = [p for p in profiles if p.height is not None and p.pre_shells > 0]
    print(title)
    print(f"  rows={len(profiles)} active_preheight_rows={len(active)}")
    metrics = [
        ("top_cover_share", [p.top_cover_share for p in active]),
        ("top_unique_share", [p.top_unique_share for p in active]),
        ("mean_pair_margin", [p.mean_pair_margin for p in active]),
        ("mean_pair_margin_per_shell", [p.mean_pair_margin_per_shell for p in active]),
        ("mean_pair_margin_norm", [p.mean_pair_margin_norm for p in active]),
        ("score_spread", [float(p.score_spread) for p in active]),
        ("top_outdegree", [float(p.top_outdegree) for p in active]),
        ("directed_3cycles", [float(p.directed_3cycles) for p in active]),
    ]
    heights = [float(p.height or 0) for p in active]
    for name, vals in metrics:
        print(
            f"    corr(height,{name}) pearson={pearson(heights, vals): .3f} "
            f"spearman={spearman(heights, vals): .3f}"
        )

    by_h: dict[int, list[DominanceProfile]] = defaultdict(list)
    for p in active:
        assert p.height is not None
        by_h[p.height].append(p)
    print("  height buckets with mean dominance metrics:")
    for h in sorted(by_h):
        bucket = by_h[h]
        if len(bucket) < 3 and h < 30:
            continue
        print(
            f"    h={h:2d} n={len(bucket):3d} "
            f"top_cover={sum(p.top_cover_share for p in bucket)/len(bucket):.3f} "
            f"pair_margin={sum(p.mean_pair_margin for p in bucket)/len(bucket):.2f} "
            f"spread={sum(p.score_spread for p in bucket)/len(bucket):.2f} "
            f"c3={sum(p.directed_3cycles for p in bucket)/len(bucket):.1f}"
        )
    print()


def print_named_profiles() -> None:
    print("C. Named hard packets: dominance tournaments")
    for label, row in named_rows():
        p = dominance_profile(label, row, full_fingerprint=True)
        print(f"  {label}: h={p.height} pre_shells={p.pre_shells} blocked_unit_mass={p.blocked_unit_mass}")
        print(
            f"    top_speed={p.top_speed} top_cover_share={p.top_cover_share:.3f} "
            f"top_unique_share={p.top_unique_share:.3f}"
        )
        print(
            f"    score_hist={p.score_hist} spread={p.score_spread} top_outdegree={p.top_outdegree} "
            f"directed_3cycles={p.directed_3cycles} scc={p.scc_sizes} Hpaths={p.hamiltonian_paths}"
        )
        print(f"    top cover loads={p.cover_loads[:6]}")
        print(f"    top unique loads={p.unique_loads[:6]}")
    print()


def print_tournament_analysis_method() -> None:
    print("D. Tournament Analysis setup and proof use")
    print("  vertices: speeds inside a fixed row; alternate vertices considered were")
    print("    denominator shells, unit numerators, dilated bands, deleted-core addresses,")
    print("    cover obligations, owner/Bprime exits, and proof tactics.")
    print("  selected quotient: speed dominance over pre-height dilated-band cover.")
    print("  pairwise observable: sum_q |U_v(q)\\U_w(q)| - |U_w(q)\\U_v(q)| for q<h(S).")
    print("  switch/gauge: q starts at the band-1 threshold 14; q<h means the row blocks q.")
    print("  tie Hamiltonian path: row order.")
    print("  preserves: who carries the blocked-unit cover before the first witness.")
    print("  destroys: continuous-time owner geometry and Bprime interval widths.")
    print("  proof interpretation: high blocking height is a long perfect-cover streak;")
    print("    if raw dominance grows but normalized dominance thins, the proof target")
    print("    becomes a dichotomy between peelable carriers and balanced-cover rigidity.")
    print()


def main() -> None:
    print("=" * 78)
    print("Codex LRC14 blocking-height dominance atlas")
    print("=" * 78)
    print(f"n={N}; q_start={Q_START}; hmax={HMAX}; core={CORE_T}")
    print()

    one_profiles = [
        dominance_profile(label, row)
        for label, row in one_stranger_rows()
    ]
    random_profiles = [
        dominance_profile(label, row)
        for label, row in random_primitive_mult14_rows()
    ]

    print_correlation_block("A. One-stranger family S(r)=7*{1..12} union {r}", one_profiles)
    print_correlation_block("B. Random primitive multiple-of-14 rows", random_profiles)
    print_named_profiles()
    print_tournament_analysis_method()
    print("E. Takeaway")
    print("  Dominance grows with blocking height only in the raw cumulative sense:")
    print("  mean pair margins rise strongly because a high-h row blocks many shells.")
    print("  Height-normalized dominance moves the other way, and the binary tournament")
    print("  saturates to a transitive order in these samples (score spread 12, no")
    print("  directed 3-cycles).  The hard rows therefore look less like arbitrary")
    print("  hierarchies and more like long balanced-cover chains with accumulated")
    print("  but diluted dominance.")
    print("  Proof route: show that every long blocking-height row either has a")
    print("  peelable carrier whose cumulative excess forces a witness after deletion,")
    print("  or else satisfies balanced-cover congruences tight enough to enter the")
    print("  ramified Q31/band-2 portal already isolated by HYP-2471/HYP-2480.")


if __name__ == "__main__":
    main()
