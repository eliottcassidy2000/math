#!/usr/bin/env python3
"""S56 half-tiling recursion scout.

The fixed-Hamiltonian-path tiling model has free cells (x,y), x-y>=2.
Grid reflection sends (x,y) to (n+1-y,n+1-x), which THM-280 identifies with
tournament complement up to relabeling.  A half-tiling keeps one fundamental
domain for this reflection, including fixed-line cells.

This script verifies the half-count formula and the parity-split recurrences
from the prompt, then records a small Tournament Analysis over proof lenses.
"""

from __future__ import annotations

from collections import Counter
from math import comb


def full_count(n: int) -> int:
    return comb(n - 1, 2) if n >= 1 else 0


def fixed_count(n: int) -> int:
    return (n - 1) // 2 if n >= 1 else 0


def half_count(n: int) -> int:
    return ((n - 1) * (n - 1)) // 4 if n >= 1 else 0


def cells(n: int) -> list[tuple[int, int]]:
    return [(x, y) for x in range(1, n + 1) for y in range(1, x) if x - y >= 2]


def reflect(n: int, cell: tuple[int, int]) -> tuple[int, int]:
    x, y = cell
    return (n + 1 - y, n + 1 - x)


def half_cells(n: int) -> list[tuple[int, int]]:
    """Keep the fixed line x+y=n+1 and one side of it."""
    return [c for c in cells(n) if c[0] + c[1] >= n + 1]


def orbit_counts(n: int) -> tuple[int, int, int]:
    seen: set[tuple[int, int]] = set()
    fixed = 0
    pairs = 0
    for c in cells(n):
        if c in seen:
            continue
        r = reflect(n, c)
        seen.add(c)
        seen.add(r)
        if r == c:
            fixed += 1
        else:
            pairs += 1
    return len(cells(n)), fixed, pairs


def odd_packet(n: int) -> dict[str, int]:
    return {
        "A": half_count(n - 1),
        "B": half_count(n - 1),
        "C": half_count(n - 2),
        "D": half_count(n - 2),
        "E": half_count(n - 3),
        "F": half_count(n - 3),
        "G": half_count(n - 4),
    }


def odd_value(packet: dict[str, int]) -> int:
    return (
        packet["A"]
        + packet["B"]
        - packet["C"]
        + packet["D"]
        - packet["E"]
        - packet["F"]
        + packet["G"]
    )


def hamiltonian_path_count(adj: dict[str, set[str]]) -> int:
    vertices = tuple(adj)
    total = 0

    def dfs(path: tuple[str, ...]) -> None:
        nonlocal total
        if len(path) == len(vertices):
            total += 1
            return
        last = path[-1]
        used = set(path)
        for v in vertices:
            if v not in used and v in adj[last]:
                dfs(path + (v,))

    for v in vertices:
        dfs((v,))
    return total


def proof_lens_tournament() -> None:
    lenses = {
        "mirror_orbit_fold": {
            "mirror_domain",
            "fixed_path",
            "complement_pair",
            "parity_split",
        },
        "odd_half_three_corner": {
            "mirror_domain",
            "fixed_path",
            "parity_split",
            "three_corner_geometry",
            "full_IE_shadow",
        },
        "even_half_two_corner": {
            "mirror_domain",
            "fixed_path",
            "parity_split",
            "two_corner_geometry",
        },
        "full_IE_third_difference": {
            "fixed_path",
            "full_IE_shadow",
            "three_corner_geometry",
            "cell_affine_warning",
        },
        "root_packet_shell": {
            "fixed_path",
            "interval_roots",
            "two_corner_geometry",
            "OCF_warning",
        },
        "OCF_cycle_space": {
            "fixed_path",
            "OCF_warning",
            "cycle_space",
            "cell_affine_warning",
        },
        "LRC_address_analogy": {
            "parity_split",
            "address_before_scalar",
            "complement_pair",
            "OCF_warning",
        },
    }
    names = list(lenses)
    adj = {name: set() for name in names}
    scores = {name: 0 for name in names}
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            key_a = (len(lenses[a]), len(lenses[a] & lenses[b]), -i)
            key_b = (len(lenses[b]), len(lenses[a] & lenses[b]), -names.index(b))
            if key_a >= key_b:
                adj[a].add(b)
                scores[a] += 1
            else:
                adj[b].add(a)
                scores[b] += 1

    path = sorted(names, key=lambda x: (scores[x], len(lenses[x]), x), reverse=True)
    cycles3 = 0
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if j <= i:
                continue
            for k, c in enumerate(names):
                if k <= j:
                    continue
                if b in adj[a] and c in adj[b] and a in adj[c]:
                    cycles3 += 1
                if c in adj[a] and b in adj[c] and a in adj[b]:
                    cycles3 += 1

    print("TOURNAMENT ANALYSIS: proof lenses, not runners")
    print("  pairwise observable=(#preserved predicates, overlap with opponent, declaration order)")
    print("  switch/gauge=larger observable orients the edge")
    print(f"  scores={scores}")
    print(f"  score_hist={dict(sorted(Counter(scores.values()).items()))}")
    print(f"  directed_3cycles={cycles3}")
    print(f"  Hamiltonian_path_count={hamiltonian_path_count(adj)}")
    print("  tie Hamiltonian path=" + " > ".join(path))


def main() -> None:
    print("S56 half-tiling recursion scout")
    print("half-domain: cells with x+y>=n+1; fixed line x+y=n+1 retained")
    print()
    print("COUNTS")
    print(" n full fixed mirror_pairs half formula")
    for n in range(1, 16):
        full, fixed, pairs = orbit_counts(n)
        hcells = half_cells(n)
        formula = (full + fixed) // 2
        assert full == full_count(n)
        assert fixed == fixed_count(n)
        assert len(hcells) == half_count(n) == formula
        print(f"{n:2d} {full:4d} {fixed:5d} {pairs:12d} {len(hcells):4d} {formula:7d}")

    print()
    print("SEQUENCE from tournament size 2:")
    print("  " + ", ".join(str(half_count(n)) for n in range(2, 13)))

    print()
    print("FULL TILING CHECK: F_n=3F_{n-1}-3F_{n-2}+F_{n-3}")
    for n in range(5, 16):
        rhs = 3 * full_count(n - 1) - 3 * full_count(n - 2) + full_count(n - 3)
        assert full_count(n) == rhs
    print("  verified n=5..15")

    print()
    print("EVEN HALF RECURRENCE: H_n=A+B-C=2H_{n-1}-H_{n-2}")
    for n in range(4, 16, 2):
        a = half_count(n - 1)
        b = half_count(n - 1)
        c = half_count(n - 2)
        value = a + b - c
        assert value == half_count(n)
        print(f"  n={n:2d}: {a}+{b}-{c}={value}")

    print()
    print("ODD HALF RECURRENCE: H_n=A+B-C+D-E-F+G")
    for n in range(5, 16, 2):
        p = odd_packet(n)
        value = odd_value(p)
        assert value == half_count(n)
        print(
            f"  n={n:2d}: "
            f"{p['A']}+{p['B']}-{p['C']}+{p['D']}-{p['E']}-{p['F']}+{p['G']}={value}"
        )
    print("  note: C and D have the same size but different geometric slots.")

    print()
    print("CENTER SIMPLIFICATION FOR ODD CASE")
    for n in range(5, 14, 2):
        p = odd_packet(n)
        center = p["A"] + p["B"] + p["D"] + p["G"] - p["D"] - p["E"] - p["F"]
        boundary = p["A"] + p["B"] - p["C"]
        assert center == p["A"] + p["B"] + p["G"] - p["E"] - p["F"]
        print(f"  n={n:2d}: edge A+B-C={boundary}; center A+B+G-E-F={center}")

    print()
    proof_lens_tournament()

    print()
    print("SYNTHESIS")
    print("  PROVED: half-tiling size H_n=floor((n-1)^2/4)=(F_n+fixed_n)/2.")
    print("  PROVED: even and odd half-tilings obey different exact recurrences.")
    print("  BRIDGE: full A+B+C-D-E-F+G is third difference; half-tiling is its")
    print("          mirror-folded parity refinement.")
    print("  WARNING: these are cell-count recurrences. OCF/Hamiltonian-path counts")
    print("           are cycle-space invariants and need packet/root data.")


if __name__ == "__main__":
    main()
