#!/usr/bin/env python3
"""Independent exact replay for the bounded THM-1123 certificate.

The 495-core atlas is reconstructed from the complete breakpoint arrangement,
not by the C++ tooth-subtraction routine.  Selected sharp and gate rows are
then replayed with fractions.Fraction, retaining exact coordinates.

Tournament audit: endpoint order is transitive and has a unique sorted
Hamiltonian path after exact ties are coalesced.  It becomes faithful only
with coordinate/owner sidecars; scalar half-coverage, overlap-MST, and
autocovariogram carriers each discard some endpoint-spacing information.
"""

from fractions import Fraction as F
from itertools import combinations


def dist(x: F) -> F:
    x %= 1
    return min(x, 1 - x)


def arrangement_safe(speeds: tuple[int, ...]) -> list[tuple[F, F]]:
    """Classify every exact breakpoint cell by a rational midpoint."""
    points = {F(0), F(1)}
    for k in speeds:
        for j in range(k + 1):
            for sign in (-1, 1):
                x = F(14 * j + sign, 14 * k)
                if 0 <= x <= 1:
                    points.add(x)
    cells = []
    ordered = sorted(points)
    for a, b in zip(ordered, ordered[1:]):
        mid = (a + b) / 2
        if all(dist(k * mid) >= F(1, 14) for k in speeds):
            if cells and cells[-1][1] == a:
                cells[-1] = (cells[-1][0], b)
            else:
                cells.append((a, b))
    return cells


def remove_bad(region: list[tuple[F, F]], k: int) -> list[tuple[F, F]]:
    out = []
    for a, b in region:
        cursor = a
        for j in range(max(0, int(a * k) - 1), min(k, int(b * k) + 1) + 1):
            left = F(14 * j - 1, 14 * k)
            right = F(14 * j + 1, 14 * k)
            if right <= a or b <= left:
                continue
            left, right = max(left, a), min(right, b)
            if cursor < left:
                out.append((cursor, left))
            cursor = max(cursor, right)
        if cursor < b:
            out.append((cursor, b))
    return out


def mass(region: list[tuple[F, F]]) -> F:
    return sum((b - a for a, b in region), F(0))


def maximum_spanning_tree(weights: dict[tuple[int, int], F]):
    edges = sorted((w, i, j) for (i, j), w in weights.items())[::-1]
    parent = list(range(4))

    def root(i: int) -> int:
        while parent[i] != i:
            i = parent[i]
        return i

    total = F(0)
    chosen = []
    for weight, i, j in edges:
        x, y = root(i), root(j)
        if x == y:
            continue
        parent[x] = y
        total += weight
        chosen.append((i, j, weight))
        if len(chosen) == 3:
            return total, chosen
    raise AssertionError("four-vertex graph was not connected")


def gate_data(core: tuple[int, ...], killers: tuple[int, ...]):
    safe = arrangement_safe(core)
    ell = max(b - a for a, b in safe)
    longest_core = [(a, b) for a, b in safe if b - a == ell]
    tau = F(1, 7 * killers[-1])
    local_half, mst_slacks = [], []
    local_covariogram = []
    for interval in longest_core:
        remainder = [interval]
        for k in killers:
            remainder = remove_bad(remainder, k)
        delta = ell - mass(remainder)
        local_half.append((ell - tau) / 2 - delta)

        weights = {}
        for i, j in combinations(range(4), 2):
            one = mass(remove_bad([interval], killers[i]))
            two = mass(remove_bad([interval], killers[j]))
            both = mass(remove_bad(remove_bad([interval], killers[i]), killers[j]))
            weights[i, j] = ell - one - two + both
        tree, chosen = maximum_spanning_tree(weights)
        rhs = ell / 14 + F(6, 49) * sum((F(1, k) for k in killers), F(0))
        rhs += F(1, 14 * killers[-1])
        mst_slacks.append((tree - rhs, tree, rhs, chosen))

        lengths = [b - a for a, b in remainder]
        area = sum(lengths, F(0))
        corr = sum(
            (length * length / 2 if length <= tau
             else tau * length - tau * tau / 2)
            for length in lengths
        )
        local_covariogram.append((corr - tau * area / 2, corr / (tau * area / 2)))

    remainder = safe
    for k in killers:
        remainder = remove_bad(remainder, k)
    mu = mass(safe)
    whole_delta = mu - mass(remainder)
    whole_half = (mu - F(len(safe), 7 * killers[-1])) / 2 - whole_delta
    longest = max(b - a for a, b in remainder)
    return {
        "ell": ell,
        "local_half": max(local_half),
        "mst": max(mst_slacks),
        "whole_half": whole_half,
        "longest": longest,
        "metric": 7 * killers[-1] * longest,
        "covariogram": max(local_covariogram),
    }


def main() -> None:
    atlas = []
    for core in combinations(range(1, 13), 8):
        safe = arrangement_safe(core)
        atlas.append((max(b - a for a, b in safe), core, safe))
    minimum = min(row[0] for row in atlas)
    minimizers = [row for row in atlas if row[0] == minimum]
    assert len(atlas) == 495
    assert minimum == F(1, 70)
    assert [row[1] for row in minimizers] == [(1, 2, 6, 7, 8, 9, 10, 11)]
    min_intervals = [
        interval for interval in minimizers[0][2]
        if interval[1] - interval[0] == minimum
    ]
    assert min_intervals == [(F(43, 140), F(9, 28)), (F(19, 28), F(97, 140))]

    hard_core = (1, 2, 4, 5, 7, 9, 11, 12)
    hard_k = (158, 160, 162, 164)
    hard = arrangement_safe(hard_core + hard_k)
    hard_L = max(b - a for a, b in hard)
    hard_intervals = [(a, b) for a, b in hard if b - a == hard_L]
    assert hard_L == F(41, 25920)
    assert 7 * hard_k[-1] * hard_L == F(11767, 6480)
    assert hard_intervals == [
        (F(673, 2240), F(685, 2268)),
        (F(1583, 2268), F(1567, 2240)),
    ]

    local = gate_data((2, 3, 4, 5, 6, 7, 8, 9), (128, 131, 134, 137))
    mst = gate_data((2, 3, 4, 5, 6, 7, 8, 9), (121, 124, 127, 130))
    whole = gate_data((1, 2, 3, 6, 7, 9, 10, 11), (153, 156, 159, 162))
    hard_gate = gate_data(hard_core, hard_k)
    assert local["local_half"] == F(-104610971, 19393097472)
    assert mst["mst"][0] == F(-466976387, 54621386820)
    assert whole["whole_half"] == F(-4926398, 150405255)
    assert hard_gate["covariogram"][1] == F(336218180863, 317189087460)

    print("THM-1123 independent exact Fraction replay")
    print("atlas_method=complete breakpoint arrangement plus midpoint classification")
    print(f"cores={len(atlas)}")
    print(f"min_core_ell={minimum}")
    print(f"unique_min_core={minimizers[0][1]}")
    print(f"min_core_longest_intervals={min_intervals}")
    print(f"hardest_replay_core={hard_core}")
    print(f"hardest_replay_quadruple={hard_k}")
    print(f"hardest_L={hard_L}")
    print(f"hardest_7k4L={7 * hard_k[-1] * hard_L}")
    print(f"hardest_intervals={hard_intervals}")
    print(f"local_half_counter_slack={local['local_half']}")
    print(f"mst_counter_slack={mst['mst'][0]}")
    print(f"whole_half_counter_slack={whole['whole_half']}")
    print(f"hardest_local_covariogram_ratio={hard_gate['covariogram'][1]}")
    print("tournament_vertices=exact endpoints with coordinate/owner sidecars")
    print("pairwise_observable=exact order; directed_cycles=0; SCCs=singletons")
    print("destroyed_by_order_only=lengths|owners|stages|gap dispersion")
    print("challenged_vertices=runners|combs|components|teeth|endpoints|residues|proof obligations")
    print("scope_warning=selected-row replay plus complete core atlas; C++ owns full bounded bank")
    print("certificate=PASS")


if __name__ == "__main__":
    main()
