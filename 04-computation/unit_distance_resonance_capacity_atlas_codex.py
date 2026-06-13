#!/usr/bin/env python3
"""Exact small-factor resonance capacity atlas for the unit-distance N* gate.

codex-2026-06-13 -- HYP-2467 / OPEN-Q-085 / T807

This script takes THM-493's resonant-product formula seriously as a finite
capacity ledger.  It enumerates every connected triangular-lattice patch through
size 9 up to translation and D6 symmetry, then asks:

    For a product N=a*b with a,b<=9, how many unit distances can the best
    resonant product A (+)_t B carry?

The count is exact inside the carrier class:

    U(A (+)_t B) = e(A)|B| + |A|e(B) + Delta_t(A,B)

where Delta_t is the norm-t displacement-spectrum correlation from THM-493.
Unlike the curated search, this is exhaustive for connected factor patches up
to size 9 and tries every relative D6 orientation of the second factor.  It is
not a proof of u(27)<=81, because arbitrary 27-point planar configurations need
not be two connected triangular factors.  It is a proof-facing obstruction:
the 3x9 resonant product lane cannot even get close to 3N, while the 4x7 lane
crosses immediately at t=3.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations


Point = tuple[int, int]
Patch = tuple[Point, ...]


UNITS: tuple[Point, ...] = (
    (1, 0),
    (0, 1),
    (-1, 1),
    (-1, 0),
    (0, -1),
    (1, -1),
)


def add(a: Point, b: Point) -> Point:
    return (a[0] + b[0], a[1] + b[1])


def sub(a: Point, b: Point) -> Point:
    return (a[0] - b[0], a[1] - b[1])


def eisen_norm(p: Point) -> int:
    m, n = p
    return m * m + m * n + n * n


def rotate60(p: Point) -> Point:
    m, n = p
    return (-n, m + n)


def reflect(p: Point) -> Point:
    # Complex conjugation in coordinates m+n*zeta6, where zeta6_bar=1-zeta6.
    m, n = p
    return (m + n, -n)


def transform_point(p: Point, transform_index: int) -> Point:
    q = reflect(p) if transform_index >= 6 else p
    turns = transform_index % 6
    for _ in range(turns):
        q = rotate60(q)
    return q


def translate_to_anchor(points: list[Point]) -> Patch:
    anchor = min(points)
    return tuple(sorted((p[0] - anchor[0], p[1] - anchor[1]) for p in points))


def canonical(points: set[Point] | list[Point] | tuple[Point, ...]) -> Patch:
    variants: list[Patch] = []
    for idx in range(12):
        variants.append(translate_to_anchor([transform_point(p, idx) for p in points]))
    return min(variants)


def enumerate_connected_patches(max_size: int) -> dict[int, list[Patch]]:
    """All connected site animals on the triangular lattice through max_size."""
    levels: dict[int, set[Patch]] = {1: {((0, 0),)}}
    for size in range(2, max_size + 1):
        current: set[Patch] = set()
        for patch in levels[size - 1]:
            pts = set(patch)
            frontier = {add(p, d) for p in pts for d in UNITS}
            for q in frontier - pts:
                current.add(canonical(pts | {q}))
        levels[size] = current
    return {size: sorted(patches) for size, patches in levels.items()}


def edge_count(patch: Patch) -> int:
    return sum(1 for a, b in combinations(patch, 2) if eisen_norm(sub(a, b)) == 1)


def spectra_by_norm(patch: Patch) -> dict[int, Counter[Point]]:
    """Ordered displacement counts grouped by Eisenstein norm."""
    spectra: dict[int, Counter[Point]] = defaultdict(Counter)
    for a in patch:
        for b in patch:
            if a == b:
                continue
            v = sub(a, b)
            spectra[eisen_norm(v)][v] += 1
    return dict(spectra)


def transform_counter(counter: Counter[Point], transform_index: int) -> Counter[Point]:
    return Counter({transform_point(v, transform_index): c for v, c in counter.items()})


@dataclass(frozen=True)
class PatchData:
    patch: Patch
    edges: int
    spectra: dict[int, Counter[Point]]


@dataclass(frozen=True)
class ProductBest:
    n: int
    a: int
    b: int
    total: int
    generic: int
    bonus: int
    rung: int | None
    orientation: int | None
    patch_a: Patch
    patch_b: Patch
    edges_a: int
    edges_b: int


def best_bonus(
    left: PatchData, right: PatchData
) -> tuple[int, int | None, int | None]:
    """Maximum Delta_t over all shared t and all D6 orientations of right."""
    best = 0
    best_t: int | None = None
    best_orientation: int | None = None
    for t in sorted((set(left.spectra) & set(right.spectra)) - {1}):
        left_counter = left.spectra[t]
        for orientation in range(12):
            right_counter = transform_counter(right.spectra[t], orientation)
            raw = sum(
                multiplicity * right_counter.get(vector, 0)
                for vector, multiplicity in left_counter.items()
            )
            bonus = raw // 2
            if bonus > best:
                best = bonus
                best_t = t
                best_orientation = orientation
    return best, best_t, best_orientation


def capacity_for_factor_pair(
    n: int, a: int, b: int, data: dict[int, list[PatchData]]
) -> ProductBest:
    best: ProductBest | None = None
    for left in data[a]:
        for right in data[b]:
            generic = left.edges * b + a * right.edges
            bonus, rung, orientation = best_bonus(left, right)
            total = generic + bonus
            candidate = ProductBest(
                n=n,
                a=a,
                b=b,
                total=total,
                generic=generic,
                bonus=bonus,
                rung=rung,
                orientation=orientation,
                patch_a=left.patch,
                patch_b=right.patch,
                edges_a=left.edges,
                edges_b=right.edges,
            )
            key = (
                candidate.total,
                candidate.bonus,
                candidate.generic,
                candidate.edges_a + candidate.edges_b,
            )
            if best is None:
                best = candidate
            else:
                best_key = (best.total, best.bonus, best.generic, best.edges_a + best.edges_b)
                if key > best_key:
                    best = candidate
    assert best is not None
    return best


def factorizations(n: int, max_size: int) -> list[tuple[int, int]]:
    out = []
    for a in range(2, int(n**0.5) + 1):
        if n % a == 0:
            b = n // a
            if a <= max_size and b <= max_size:
                out.append((a, b))
    return out


def patch_norm_summary(data: PatchData) -> str:
    entries = []
    for t in sorted(data.spectra):
        unordered_pairs = sum(data.spectra[t].values()) // 2
        if unordered_pairs:
            entries.append(f"{t}:{unordered_pairs}")
    return "{" + ", ".join(entries[:8]) + ("..." if len(entries) > 8 else "") + "}"


def directed_cycles(names: list[str], beats: dict[tuple[str, str], bool]) -> list[tuple[str, str, str]]:
    cycles = []
    for a, b, c in combinations(names, 3):
        triples = ((a, b, c), (a, c, b))
        for x, y, z in triples:
            if beats[(x, y)] and beats[(y, z)] and beats[(z, x)]:
                cycles.append((x, y, z))
    return cycles


def scc_sizes(names: list[str], edges: dict[str, set[str]]) -> list[int]:
    reverse = {name: set() for name in names}
    for u, outs in edges.items():
        for v in outs:
            reverse[v].add(u)

    seen: set[str] = set()
    order: list[str] = []

    def dfs(u: str) -> None:
        seen.add(u)
        for v in edges[u]:
            if v not in seen:
                dfs(v)
        order.append(u)

    for name in names:
        if name not in seen:
            dfs(name)

    seen.clear()
    sizes: list[int] = []

    def rdfs(u: str) -> int:
        seen.add(u)
        total = 1
        for v in reverse[u]:
            if v not in seen:
                total += rdfs(v)
        return total

    for name in reversed(order):
        if name not in seen:
            sizes.append(rdfs(name))
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(names: list[str], wins: dict[tuple[str, str], bool]) -> int:
    # Small dynamic program over subsets; tie Hamiltonian path is declaration order.
    index = {name: i for i, name in enumerate(names)}
    dp: dict[tuple[int, int], int] = {}
    for name in names:
        dp[(1 << index[name], index[name])] = 1
    full = (1 << len(names)) - 1
    for mask in range(1 << len(names)):
        for last in range(len(names)):
            val = dp.get((mask, last), 0)
            if not val:
                continue
            last_name = names[last]
            for nxt, nxt_name in enumerate(names):
                if mask & (1 << nxt):
                    continue
                if wins[(last_name, nxt_name)]:
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + val
    return sum(dp.get((full, last), 0) for last in range(len(names)))


def tournament_analysis() -> None:
    print()
    print("Tournament Analysis: proof-obligation vertices")
    print("  challenged vertex choices considered:")
    print("    points, unit edges, connected patches, factor sizes, displacement vectors,")
    print("    Loeschian rungs, carrier lattices, and proof obligations.")
    print("  chosen quotient: proof obligations; preserves route leverage toward u(27)<=81,")
    print("    destroys individual coordinates, but keeps which obstruction can own a proof step.")
    print("  challenged assumption: the right vertices need not be points or directions;")
    print("    the capacity resource itself can be the vertex.")
    obligations = {
        "small_factor_capacity": {
            "exact": 5,
            "upper_bound": 4,
            "transfer": 5,
            "computation": 5,
            "scope": 3,
            "proof_ready": 4,
        },
        "generic_edge_budget": {
            "exact": 5,
            "upper_bound": 3,
            "transfer": 4,
            "computation": 5,
            "scope": 4,
            "proof_ready": 4,
        },
        "resonance_bonus_spectrum": {
            "exact": 5,
            "upper_bound": 5,
            "transfer": 5,
            "computation": 4,
            "scope": 4,
            "proof_ready": 3,
        },
        "free_patch_annealing": {
            "exact": 2,
            "upper_bound": 2,
            "transfer": 4,
            "computation": 5,
            "scope": 5,
            "proof_ready": 2,
        },
        "carrier_ladder_gate": {
            "exact": 4,
            "upper_bound": 3,
            "transfer": 4,
            "computation": 4,
            "scope": 5,
            "proof_ready": 3,
        },
        "global_upper_bound": {
            "exact": 1,
            "upper_bound": 5,
            "transfer": 3,
            "computation": 1,
            "scope": 5,
            "proof_ready": 1,
        },
    }
    names = list(obligations)
    criteria = list(next(iter(obligations.values())).keys())
    wins: dict[tuple[str, str], bool] = {}
    edges = {name: set() for name in names}
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if i == j:
                continue
            votes_a = 0
            votes_b = 0
            for criterion in criteria:
                av = obligations[a][criterion]
                bv = obligations[b][criterion]
                votes_a += av > bv
                votes_b += bv > av
            if votes_a == votes_b:
                a_wins = i < j
            else:
                a_wins = votes_a > votes_b
            wins[(a, b)] = a_wins
            if a_wins:
                edges[a].add(b)
    scores = Counter({name: len(edges[name]) for name in names})
    hist = Counter(scores.values())
    cycles = directed_cycles(names, wins)
    print(f"  score sequence: {sorted(scores.items(), key=lambda kv: (-kv[1], kv[0]))}")
    print(f"  score histogram: {dict(sorted(hist.items()))}")
    print(f"  directed 3-cycles: {len(cycles)}")
    if cycles:
        print(f"  sample cycles: {cycles[:4]}")
    print(f"  SCC sizes: {scc_sizes(names, edges)}")
    print(f"  Hamiltonian path count: {hamiltonian_path_count(names, wins)}")


def main() -> None:
    max_size = 9
    patches = enumerate_connected_patches(max_size)
    data: dict[int, list[PatchData]] = {}
    for size, shape_list in patches.items():
        data[size] = [
            PatchData(patch=patch, edges=edge_count(patch), spectra=spectra_by_norm(patch))
            for patch in shape_list
        ]

    print("Unit-distance resonance capacity atlas (connected triangular factors)")
    print(f"  max factor size: {max_size}")
    print("  connected patch counts up to translation and D6:")
    for size in range(1, max_size + 1):
        max_edges = max(d.edges for d in data[size])
        print(f"    size {size}: {len(data[size])} patches, max internal edges {max_edges}")

    print()
    print("Factor-product capacity scan")
    best_by_n: dict[int, ProductBest] = {}
    for n in range(24, 33):
        pairs = factorizations(n, max_size)
        if not pairs:
            print(f"  N={n:2d}  3N={3*n:3d}  no factorization inside <= {max_size} atlas")
            continue
        best_for_n: ProductBest | None = None
        for a, b in pairs:
            candidate = capacity_for_factor_pair(n, a, b, data)
            if best_for_n is None or candidate.total > best_for_n.total:
                best_for_n = candidate
        assert best_for_n is not None
        best_by_n[n] = best_for_n
        status = "BEATS" if best_for_n.total > 3 * n else ("TIES" if best_for_n.total == 3 * n else "below")
        rung = "-" if best_for_n.rung is None else str(best_for_n.rung)
        print(
            f"  N={n:2d}  3N={3*n:3d}  best={best_for_n.total:3d} ({status})"
            f"  [{best_for_n.a}x{best_for_n.b}, generic={best_for_n.generic},"
            f" bonus={best_for_n.bonus}, t={rung}, e=({best_for_n.edges_a},{best_for_n.edges_b})]"
        )

    print()
    print("27-vs-28 gate witness")
    n27 = best_by_n[27]
    n28 = best_by_n[28]
    print(
        f"  N=27 connected two-factor lane maxes at {n27.total},"
        f" deficit to 3N is {3*n27.n - n27.total}."
    )
    print(
        f"  N=28 connected two-factor lane reaches {n28.total},"
        f" excess over 3N is {n28.total - 3*n28.n}."
    )
    print(f"  N=28 champion A={n28.patch_a}")
    print(f"  N=28 champion B={n28.patch_b}")

    print()
    print("Size-3 stress test against all connected 9-patches")
    for left in sorted(data[3], key=lambda d: (-d.edges, d.patch)):
        best: ProductBest | None = None
        for right in data[9]:
            generic = left.edges * 9 + 3 * right.edges
            bonus, rung, orientation = best_bonus(left, right)
            candidate = ProductBest(
                n=27,
                a=3,
                b=9,
                total=generic + bonus,
                generic=generic,
                bonus=bonus,
                rung=rung,
                orientation=orientation,
                patch_a=left.patch,
                patch_b=right.patch,
                edges_a=left.edges,
                edges_b=right.edges,
            )
            if best is None or candidate.total > best.total:
                best = candidate
        assert best is not None
        rung = "-" if best.rung is None else str(best.rung)
        print(
            f"  shape {left.patch}: e={left.edges}, norms={patch_norm_summary(left)}"
            f" -> best total {best.total} [generic={best.generic}, bonus={best.bonus}, t={rung}]"
        )

    print()
    print("Interpretation")
    print("  The densest 3-point factor is K3, but K3 has no norm-t displacement with t>1,")
    print("  so it is resonance-free. The 3-point factors that do carry a transverse")
    print("  displacement lose an internal edge; even their best bonus against all 9-patches")
    print("  does not recover the K3 generic budget. The first small factor that can be both")
    print("  edge-dense and norm-3-bearing is size 4, and 4x7 is exactly the Moser crossing.")

    tournament_analysis()


if __name__ == "__main__":
    main()
